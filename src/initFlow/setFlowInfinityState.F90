!
!      ******************************************************************
!      *                                                                *
!      * File:          setFlowInfinityState.f90                        *
!      * Author:        Edwin van der Weide, Georgi Kalitzin            *
!      * Starting date: 02-21-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setFlowInfinityState
!
!      ******************************************************************
!      *                                                                *
!      * setFlowInfinityState sets the free stream state vector of      *
!      * the flow variables. If nothing is specified for each of the    *
!      * farfield boundaries, these values will be taken to define the  *
!      * free stream.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use flowVarRefState
       use inputPhysics
       use paramTurb
       implicit none
!
!      Local variables
!
       integer(kind=intType) :: ierr

       real(kind=realType) :: nuInf, ktmp, uInf2
!
!      Function definition
!
       real(kind=realType) :: saNuKnownEddyRatio

       ! Dummy parameters
       real(kind=realType) :: vinf, zinf
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Compute the velocity squared based on MachCoef;
       ! needed for the initialization of the turbulent energy,
       ! especially for moving geometries.

       uInf2 = MachCoef*MachCoef*gammaInf*pInf/rhoInf

       ! Allocate the memory for wInf.
#ifndef USE_TAPENADE
       if( allocated(wInf)) deallocate(wInf)
      
       allocate(wInf(nw), stat=ierr)
       if(ierr /= 0)                             &
         call returnFail("setFlowReferenceState", &
                        "Memory allocation failure for wInf")
#endif

       ! zero out the winf first
       wInf(:) = zero

       ! Set the reference value of the flow variables, except the total
       ! energy. This will be computed at the end of this routine.

       wInf(irho) = rhoInf
       wInf(ivx)  = uInf*velDirFreestream(1)
       wInf(ivy)  = uInf*velDirFreestream(2)
       wInf(ivz)  = uInf*velDirFreestream(3)
       
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
       call etotArray(rhoInf, uInf, vInf, zInf, pInfCorr, ktmp, &
                       wInf(irhoE), kPresent, 1)
       
       end subroutine setFlowInfinityState
