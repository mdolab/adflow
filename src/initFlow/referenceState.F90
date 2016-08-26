!
!       File:          referenceState.f90                              
!       Author:        Edwin van der Weide, Seonghyeon Hahn            
!       Starting date: 05-29-2003                                      
!       Last modified: 04-22-2006                                      
!
subroutine referenceState
  !
  !       The original version has been nuked since the computations are 
  !       no longer necessary when calling from python                   
  !       This is the most compliclated routine in all of SUMb. It is    
  !       stupidly complicated. This is most likely the reason your      
  !       derivatives are wrong. You don't understand this routine       
  !       and its effects.                                               
  !       This routine *requries* the following as input:                
  !       Mach, pInfDim, TInfDim, rhoInfDim, rGasDim (machCoef non-SA    
  !        turbulence only)                                              
  !       Optionally, pRef, rhoRef and Tref are used if they are         
  !       are non-negative. This only happens when you want the equations
  !       normalized by values other than the freestream                 
  !      * This routine computes as output:  
  !      *   muInfDim, (unused anywhere in code)
  !         pRef, rhoRef, Tref, muRef, timeRef ('dimensional' reference) 
  !         pInf, pInfCorr, rhoInf, uInf, rGas, muInf, gammaInf and wInf 
  !         (Non-dimensionalized values used in actual computations)     
  !
  use constants
  use paramTurb
  use inputPhysics, only : equations, Mach, machCoef, &
       muSuthDim, TSuthDim, velDirFreeStream, &
       rGasDim, SSuthDim, eddyVisInfRatio, turbModel, turbIntensityInf
  use flowVarRefState, only : pInfDim, TinfDim, rhoInfDim,  &
       muInfDim, &
       pRef, rhoRef, Tref, muRef, timeRef, &
       pInf, pInfCorr, rhoInf, uInf, rGas, muInf, gammaInf, wInf, &
       nw, nwf, kPresent, wInf
  use flowUtils, only : computeGamma, eTot

  implicit none
 
  integer(kind=intType) :: sps, nn, mm, ierr
  real(kind=realType) :: gm1, ratio
  real(kind=realType) :: nuInf, ktmp, uInf2
  real(kind=realType) :: saNuKnownEddyRatio
  real(kind=realType) :: vinf, zinf, tmp1(1), tmp2(1)

  ! Compute the dimensional viscosity from Sutherland's law
  muInfDim = muSuthDim &
       * ((TSuthDim + SSuthDim)/(TInfDim + SSuthDim)) &
       * ((TInfDim/TSuthDim)**1.5_realType)

  ! Set the reference values. They *COULD* be different from the
  ! free-stream values for an internal flow simulation. For now,
  ! we just use the actual free stream values. 
  pref = PInfDim
  tref = TInfDim
  rhoref = rhoInfDim

  ! Compute the value of muRef, such that the nonDimensional
  ! equations are identical to the dimensional ones.
  ! Note that in the non-dimensionalization of muRef there is
  ! a reference length. However this reference length is 1.0
  ! in this code, because the coordinates are converted to
  ! meters.

  muRef = sqrt(pRef*rhoRef)

  ! Compute timeRef for a correct nonDimensionalization of the
  ! unsteady equations. Some story as for the reference viscosity
  ! concerning the reference length.

  timeRef = sqrt(rhoRef/pRef)

  ! Compute the nonDimensional pressure, density, velocity,
  ! viscosity and gas constant.

  pInf   = pInfDim/pRef
  rhoInf = rhoInfDim/rhoRef
  uInf   = Mach*sqrt(gammaInf*pInf/rhoInf)
  RGas   = RGasDim*rhoRef*TRef/pRef
  muInf  = muInfDim/muRef
  tmp1(1) = TinfDim
  call computeGamma(tmp1, tmp2, 1)
  gammaInf = tmp2(1)

  ! ----------------------------------------
  !      Compute the final wInf
  ! ----------------------------------------

  ! Allocate the memory for wInf if necessary
#ifndef USE_TAPENADE
  if( allocated(wInf)) deallocate(wInf)
  allocate(wInf(nw), stat=ierr)
#endif

  ! zero out the winf first
  wInf(:) = zero

  ! Set the reference value of the flow variables, except the total
  ! energy. This will be computed at the end of this routine.
  
  wInf(irho) = rhoInf
  wInf(ivx)  = uInf*velDirFreestream(1)
  wInf(ivy)  = uInf*velDirFreestream(2)
  wInf(ivz)  = uInf*velDirFreestream(3)
       
  ! Compute the velocity squared based on MachCoef. This gives a
  ! better indication of the 'speed' of the flow so the turubulence
  ! intensity ration is more meaningful especially for moving
  ! geometries. (Not used in SA model)

  uInf2 = MachCoef*MachCoef*gammaInf*pInf/rhoInf
  
  ! Set the turbulent variables if transport variables are to be
  ! solved. We should be checking for RANS equations here,
  ! however, this code is included in block res. The issue is
  ! that for frozen turbulence (or ANK jacobian) we call the
  ! block_res with equationType set to Laminar even though we are
  ! actually solving the rans equations. The issue is that, the
  ! freestream turb variables will be changed to zero, thus
  ! changing the solution. Insteady we check if nw > nwf which
  ! will accomplish the same thing. 
       
  if(nw > nwf) then 

     nuInf  = muInf/rhoInf

     select case(turbModel)

     case (spalartAllmaras, spalartAllmarasEdwards)
        
        wInf(itu1) = saNuKnownEddyRatio(eddyVisInfRatio, nuInf)
        
        !=============================================================
        
     case (komegaWilcox, komegaModified, menterSST)
        
        wInf(itu1) = 1.5_realType*uInf2*turbIntensityInf**2
        wInf(itu2) = wInf(itu1)/(eddyVisInfRatio*nuInf)

        !=============================================================
        
     case (ktau)
        
        wInf(itu1) = 1.5_realType*uInf2*turbIntensityInf**2
        wInf(itu2) = eddyVisInfRatio*nuInf/wInf(itu1)
        
        !=============================================================
        
     case (v2f)
        
        wInf(itu1) = 1.5_realType*uInf2*turbIntensityInf**2
        wInf(itu2) = 0.09_realType*wInf(itu1)**2 &
             / (eddyVisInfRatio*nuInf)
        wInf(itu3) = 0.666666_realType*wInf(itu1)
        wInf(itu4) = 0.0_realType

     end select
     
  endif

  ! Set the value of pInfCorr. In case a k-equation is present
  ! add 2/3 times rho*k.
       
  pInfCorr = pInf
  if( kPresent ) pInfCorr = pInf + two*third*rhoInf*wInf(itu1)

  ! Compute the free stream total energy.
  
  ktmp = zero
  if( kPresent ) ktmp = wInf(itu1)
  vInf = zero
  zInf = zero
  call etot(rhoInf, uInf, vInf, zInf, pInfCorr, ktmp, &
       wInf(irhoE), kPresent)

end subroutine referenceState

