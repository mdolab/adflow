module initializeFlow

  use constants, only : intType, realType, maxStringLen

  implicit none
  save

contains

  subroutine referenceState
    !
    !       The original version has been nuked since the computations are
    !       no longer necessary when calling from python
    !       This is the most compliclated routine in all of ADflow. It is
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
         pRef, rhoRef, Tref, muRef, timeRef, uRef, hRef, &
         pInf, pInfCorr, rhoInf, uInf, rGas, muInf, gammaInf, wInf, &
         nw, nwf, kPresent, wInf
    use flowUtils, only : computeGamma, eTot
    use turbUtils, only : saNuKnownEddyRatio
    implicit none

    integer(kind=intType) :: sps, nn, mm, ierr
    real(kind=realType) :: gm1, ratio
    real(kind=realType) :: nuInf, ktmp, uInf2
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
    hRef = pRef/rhoRef
    uRef = sqrt(hRef)

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




  ! ----------------------------------------------------------------------
  !                                                                      |
  !                    No Tapenade Routine below this line               |
  !                                                                      |
  ! ----------------------------------------------------------------------

#ifndef  USE_TAPENADE
  subroutine infChangeCorrection(oldWinf)
    ! Adjust the flow states to a change in wInf
    use constants
    use blockPointers, only : il, jl, kl, w, nDom, d2wall
    use flowVarRefState, only : wInf, nwf, nw
    use inputPhysics , only : equations
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use haloExchange, only : whalo2
    use flowUtils, only : adjustInflowAngle
    use oversetData, only : oversetPresent
    use iteration, only : currentLevel
    use turbbcRoutines, only : applyallTurbBCthisblock, bcTurbTreatment
    use BCRoutines, only : applyallBC_block
    use utils, only : setPointers, mynorm2
    implicit none

    real(kind=realType), intent(in), dimension(nwf) :: oldWinf
    integer(kind=intType) :: sps, nn, i, j, k, l
    real(kind=realType) :: deltaWinf(nwf)

    ! Make sure we have the updated wInf
    call adjustInflowAngle()
    call referenceState

    deltaWinf = Winf(1:nwf) - oldWinf(1:nwf)

    if (mynorm2(deltaWinf) < 1e-12) then
       ! The change deltaWinf is so small, (or zero) don't do the
       ! update and just return. This will save some time when the
       ! solver is called with the same AP conditions multiple times,
       ! such as during a GS AS solution

       return
    end if

    ! Loop over all the blocks, adding the subtracting off the oldWinf
    do sps=1, nTimeIntervalsSpectral
       do nn=1, nDom
          call setPointers(nn, 1_intType, sps)
          do k=2, kl
             do j=2, jl
                do i=2, il
                   do l=1, nwf
                      w(i, j, k, l) = w(i, j, k, l) + deltaWinf(l)
                   end do
                end do
             end do
          end do
          call applyAllBC_block(.True.)
       end do
    end do

    ! Exchange values
    call whalo2(currentLevel, 1_intType, nw, .True., .True., .True.)

    ! Need to re-apply the BCs. The reason is that BC halos behind
    ! interpolated cells need to be recomputed with their new
    ! interpolated values from actual compute cells. Only needed for
    ! overset.
    if (oversetPresent) then
       do sps=1,nTimeIntervalsSpectral
          do nn=1,nDom
             call setPointers(nn, 1, sps)
             if (equations == RANSequations) then
                call BCTurbTreatment
                call applyAllTurbBCthisblock(.True.)
             end if
             call applyAllBC_block(.True.)
          end do
       end do
    end if

end subroutine infChangeCorrection


  ! Section out the BCdata setup so that it can by called from python when needed
  subroutine updateBCDataAllLevels()
    ! sets the prescribed boundary data from the CGNS arrays

    use constants
    use iteration, only : groundLevel
    use bcdata, only : setbcdataFineGrid, setBCDataCoarseGrid
    implicit none

    ! Allocate the memory for the prescribed boundary data at the
    ! boundary faces and determine the data for both the fine grid.

    groundLevel = 1

    ! Determine the reference state.
    call referenceState

    call setBCDataFineGrid(.true.)

    ! Determine the prescribed data on the coarse grid levels
    ! by interpolation.
#ifndef USE_TAPENADE
    call setBCDataCoarseGrid
#endif

  end subroutine updateBCDataAllLevels

  subroutine initFlow
    !
    !       initFlow allocates the
    !       memory for and initializes the flow variables. In case a
    !       restart is performed the owned variables are read from the
    !       previous solution file(s).
    !
    use constants
    use block, only : flowDoms
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use variableReading, only : halosRead

    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: sps, level, nLevels

    ! Determine the number of multigrid levels.

    nLevels = ubound(flowDoms,2)

    ! As some boundary conditions can be treated in multiple ways,
    ! some memory allocated must be released again.

     call releaseExtraMemBcs

    ! Determine for the time spectral mode the matrices for the
    ! time derivatives.
    call timeSpectralMatrices

    ! Loop over the number of spectral solutions to allocate
    ! the memory for the w-variables and p on the fine grid.
    do sps=1,nTimeIntervalsSpectral
       call allocMemFlovarPart1(sps, 1_intType)
    enddo

    ! Allocate the memory for the solution variables on the coarse
    ! grid levels and the memory for the dependent flow variables,
    ! residuals, etc, on all multigrid levels.
    do sps=1,nTimeIntervalsSpectral
       call allocMemFlovarPart2(sps, 1_intType)

       do level=2,nLevels
          call allocMemFlovarPart1(sps, level)
          call allocMemFlovarPart2(sps, level)
       enddo
    enddo

    ! Initialize free stream field
    call initFlowfield

    ! Initialize the dependent flow variables and the halo values.

    call initDepvarAndHalos(halosRead)

  end subroutine initFlow


  subroutine allocMemFlovarPart1(sps,level)
    !
    !       allocMemFlovarPart1 allocates the memory for the flow
    !       variables w and p for all the blocks on the given multigrid
    !       level and spectral solution sps.
    !
    use constants
    use block, only : flowDoms, nDOm
    use flowVarRefState, only : nw, nwf, nt1, nt2
    use inputPhysics, only : equationMode, gammaConstant
    use inputUnsteady, only :timeIntegrationScheme
    use inputIteration, only : mgStartLevel, turbTreatment
    use iteration, only : nOldLevels
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps, level
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn
    integer(kind=intType) :: il, jl, kl, ie, je, ke, ib, jb, kb

    ! Loop over the domains.

    domains: do nn=1,nDom

       ! Store some dimensions a bit easier.

       il = flowDoms(nn,level,sps)%il
       jl = flowDoms(nn,level,sps)%jl
       kl = flowDoms(nn,level,sps)%kl

       ie = flowDoms(nn,level,sps)%ie
       je = flowDoms(nn,level,sps)%je
       ke = flowDoms(nn,level,sps)%ke

       ib = flowDoms(nn,level,sps)%ib
       jb = flowDoms(nn,level,sps)%jb
       kb = flowDoms(nn,level,sps)%kb

       ! Allocate the memory for the independent variables.
       ! Memory is allocated for the turbulent variables (if any) if
       ! the current level is smaller or equal to the multigrid start
       ! level or if the turbulent transport equations are solved in
       ! a coupled manner.

       if(level <= mgStartlevel .or. turbTreatment == coupled) then
          allocate(flowDoms(nn,level,sps)%w(0:ib,0:jb,0:kb,1:nw), &
               stat=ierr)
       else
          allocate(flowDoms(nn,level,sps)%w(0:ib,0:jb,0:kb,1:nwf), &
               stat=ierr)
       endif
       if(ierr /= 0)                           &
            call terminate("allocMemFlovarPart1", &
            "Memory allocation failure for w")

       ! Alloc mem for nodal gradients
       allocate(flowDoms(nn,level,sps)%ux(il,jl,kl), stat=ierr)
       allocate(flowDoms(nn,level,sps)%uy(il,jl,kl), stat=ierr)
       allocate(flowDoms(nn,level,sps)%uz(il,jl,kl), stat=ierr)

       allocate(flowDoms(nn,level,sps)%vx(il,jl,kl), stat=ierr)
       allocate(flowDoms(nn,level,sps)%vy(il,jl,kl), stat=ierr)
       allocate(flowDoms(nn,level,sps)%vz(il,jl,kl), stat=ierr)

       allocate(flowDoms(nn,level,sps)%wx(il,jl,kl), stat=ierr)
       allocate(flowDoms(nn,level,sps)%wy(il,jl,kl), stat=ierr)
       allocate(flowDoms(nn,level,sps)%wz(il,jl,kl), stat=ierr)

       allocate(flowDoms(nn,level,sps)%qx(il,jl,kl), stat=ierr)
       allocate(flowDoms(nn,level,sps)%qy(il,jl,kl), stat=ierr)
       allocate(flowDoms(nn,level,sps)%qz(il,jl,kl), stat=ierr)

       ! Allocate memory for the pressure.
       allocate(flowDoms(nn,level,sps)%p(0:ib,0:jb,0:kb), stat=ierr)
       if(ierr /= 0)                           &
            call terminate("allocMemFlovarPart1", &
            "Memory allocation failure for p")

       ! Allocate memory for the speed of sound squared
       allocate(flowDoms(nn,level,sps)%aa(0:ib,0:jb,0:kb), stat=ierr)
       if(ierr /= 0)                           &
            call terminate("allocMemFlovarPart1", &
            "Memory allocation failure for p")

       ! The eddy viscosity for eddy viscosity models.
       ! Although a dependent variable, it is allocated on all grid
       ! levels, because the eddy viscosity might be frozen in the
       ! multigrid.

       ! Always allocate rev due to reverse mode AD- Peter Lyu. Also
       ! zero so that it doesn't affect laminar cases.
       allocate(flowDoms(nn,level,sps)%rev(0:ib,0:jb,0:kb), &
            stat=ierr)
       flowDoms(nn, level, sps)%rev = zero
       if(ierr /= 0)                           &
            call terminate("allocMemFlovarPart1", &
            "Memory allocation failure for rev")
       !endif

       ! If this is the finest grid some more memory must be allocated.

       fineLevelTest: if(level == 1) then

          ! Allocate the memory for gamma and initialize it to
          ! the constant gamma value.

          allocate(flowDoms(nn,level,sps)%gamma(0:ib,0:jb,0:kb),    &
               stat=ierr)
          if(ierr /= 0)                           &
               call terminate("allocMemFlovarPart1", &
               "Memory allocation failure for gamma.")

          flowDoms(nn,level,sps)%gamma = gammaConstant

          ! The laminar viscosity for viscous computations.
          ! Always allocate rlv due to reverse mode - Peter Lyu
          !if( viscous ) then
          allocate(flowDoms(nn,level,sps)%rlv(0:ib,0:jb,0:kb), &
               stat=ierr)
          if(ierr /= 0)                           &
               call terminate("allocMemFlovarPart1", &
               "Memory allocation failure for rlv")
          !endif

          ! The state vectors in the past for unsteady computations.

          if(equationMode          == unsteady .and. &
               timeIntegrationScheme == BDF) then
             allocate( &
                  flowDoms(nn,level,sps)%wOld(nOldLevels,2:il,2:jl,2:kl,nw), &
                  stat=ierr)
             if(ierr /= 0)                           &
                  call terminate("allocMemFlovarPart1", &
                  "Memory allocation failure for wOld")

             ! Initialize wOld to zero, such that it is initialized.
             ! The actual values do not matter.

             flowDoms(nn,level,sps)%wOld = zero

             ! Added by HDN
          else if ( equationMode == unsteady .and. &
               timeIntegrationScheme == MD) then
             allocate( &
                  flowDoms(nn,level,sps)%wOld(nOldLevels,2:il,2:jl,2:kl,nw), &
                  stat=ierr)
             if(ierr /= 0)                           &
                  call terminate("allocMemFlovarPart1", &
                  "Memory allocation failure for wOld")

             flowDoms(nn,level,sps)%wOld = zero
          endif

          ! If this is the 1st spectral solution (note that we are
          ! already on the finest grid) and the rans equations are
          ! solved, allocate the memory for the arrays used for the
          ! implicit boundary condition treatment. Normally this should
          ! only be allocated for RANS but the derivative calcs require
          ! these be allocated.

          sps1RansTest: if(sps == 1) then
             allocate(flowDoms(nn,level,sps)%bmti1(je,ke,nt1:nt2,nt1:nt2), &
                  flowDoms(nn,level,sps)%bmti2(je,ke,nt1:nt2,nt1:nt2), &
                  flowDoms(nn,level,sps)%bmtj1(ie,ke,nt1:nt2,nt1:nt2), &
                  flowDoms(nn,level,sps)%bmtj2(ie,ke,nt1:nt2,nt1:nt2), &
                  flowDoms(nn,level,sps)%bmtk1(ie,je,nt1:nt2,nt1:nt2), &
                  flowDoms(nn,level,sps)%bmtk2(ie,je,nt1:nt2,nt1:nt2), &
                  flowDoms(nn,level,sps)%bvti1(je,ke,nt1:nt2), &
                  flowDoms(nn,level,sps)%bvti2(je,ke,nt1:nt2), &
                  flowDoms(nn,level,sps)%bvtj1(ie,ke,nt1:nt2), &
                  flowDoms(nn,level,sps)%bvtj2(ie,ke,nt1:nt2), &
                  flowDoms(nn,level,sps)%bvtk1(ie,je,nt1:nt2), &
                  flowDoms(nn,level,sps)%bvtk2(ie,je,nt1:nt2), &
                  stat=ierr)
             if(ierr /= 0)                           &
                  call terminate("allocMemFlovarPart1", &
                  "Memory allocation failure for bmti1, etc")

          endif sps1RansTest
       endif fineLevelTest

    enddo domains

  end subroutine allocMemFlovarPart1

  !      ==================================================================

  subroutine allocMemFlovarPart2(sps, level)
    !
    !       AllocMemFlovarPart2 allocates the memory for the dependent
    !       flow variables and iteration variables for all the blocks on
    !       the given multigrid level and spectral solution sps. Some
    !       variables are only allocated on the coarser grids, e.g. the
    !       multigrid forcing terms and the state vector upon entrance on
    !       the mg level. Other variables are only allocated on the finest
    !       mesh. These are typically dependent variables like laminar
    !       viscosity, or residuals, time step, etc. Exceptions are
    !       pressure and eddy viscosity. Although these are dependent
    !       variables, they are allocated on all grid levels.
    !
    use block
    use constants
    use flowVarRefState
    use inputPhysics
    use inputDiscretization
    use inputIteration
    use inputUnsteady
    use iteration
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps, level
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, mm
    integer(kind=intType) :: il, jl, kl, ie, je, ke, ib, jb, kb

    ! Loop over the domains.

    domains: do nn=1,nDom

       ! Store some dimensions a bit easier.

       il = flowDoms(nn,level,sps)%il
       jl = flowDoms(nn,level,sps)%jl
       kl = flowDoms(nn,level,sps)%kl

       ie = flowDoms(nn,level,sps)%ie
       je = flowDoms(nn,level,sps)%je
       ke = flowDoms(nn,level,sps)%ke

       ib = flowDoms(nn,level,sps)%ib
       jb = flowDoms(nn,level,sps)%jb
       kb = flowDoms(nn,level,sps)%kb

       ! Block is moving. Allocate the memory for s, sFaceI,
       ! sFaceJ and sFaceK.

       allocate(flowDoms(nn,level,sps)%s(ie,je,ke,3),  &
            flowDoms(nn,level,sps)%sFaceI(0:ie,je,ke), &
            flowDoms(nn,level,sps)%sFaceJ(ie,0:je,ke), &
            flowDoms(nn,level,sps)%sFaceK(ie,je,0:ke), stat=ierr)
       if(ierr /= 0)                              &
            call terminate("allocMemFlovarPart2", &
            "Memory allocation failure for s, &
            &sFaceI, sFaceJ and sFaceK.")

       ! Extra face velocities for ALE
       if (equationMode == unSteady .and. useALE) then
          allocate( &
               flowDoms(nn,level,sps)%sVeloIALE(0:ie,je,ke,3), &
               flowDoms(nn,level,sps)%sVeloJALE(ie,0:je,ke,3), &
               flowDoms(nn,level,sps)%sVeloKALE(ie,je,0:ke,3), &
               flowDoms(nn,level,sps)%sFaceIALE(0:nALEsteps,0:ie,je,ke), &
               flowDoms(nn,level,sps)%sFaceJALE(0:nALEsteps,ie,0:je,ke), &
               flowDoms(nn,level,sps)%sFaceKALE(0:nALEsteps,ie,je,0:ke), stat=ierr)
          if(ierr /= 0)                              &
               call terminate("allocMemFlovarPart2", &
               "Memory allocation failure for &
               sVeloIALE, sVeloJALE and sVeloKALE; &
               sFaceIALE, sFaceJALE and sFaceKALE.")
       end if


       ! Test if we are on the finest mesh.

       fineLevelTest: if(level == 1) then

          ! Allocate the memory that must always be allocated.

          allocate( &
               flowDoms(nn,level,sps)%dw(0:ib,0:jb,0:kb,1:nw),  &
               flowDoms(nn,level,sps)%fw(0:ib,0:jb,0:kb,1:nwf), &
               flowDoms(nn,level,sps)%dtl(1:ie,1:je,1:ke),      &
               flowDoms(nn,level,sps)%radI(1:ie,1:je,1:ke),     &
               flowDoms(nn,level,sps)%radJ(1:ie,1:je,1:ke),     &
               flowDoms(nn,level,sps)%radK(1:ie,1:je,1:ke),     &
               flowDoms(nn,level,sps)%scratch(0:ib,0:jb,0:kb,10), &
               flowDoms(nn,level,sps)%shockSensor(0:ib, 0:jb, 0:kb), &
               stat=ierr)
          if(ierr /= 0)                              &
               call terminate("allocMemFlovarPart2", &
               "Memory allocation failure for dw, fw, dwOld, fwOld, &
               &gamma, dtl and the spectral radii.")

          ! Initialize dw and fw to zero.

          flowDoms(nn,level,sps)%dw = zero
          flowDoms(nn,level,sps)%fw = zero


          ! Extra variables for ALE
          if (equationMode == unSteady .and. useALE) then
             allocate( &
                  flowDoms(nn,level,sps)%dwALE(0:nALEsteps,0:ib,0:jb,0:kb,1:nw),  &
                  flowDoms(nn,level,sps)%fwALE(0:nALEsteps,0:ib,0:jb,0:kb,1:nwf), &
                  stat=ierr)
             if(ierr /= 0)                              &
                  call terminate("allocMemFlovarPart2", &
                  "Memory allocation failure for dwALE, fwALE.")

             flowDoms(nn,level,sps)%dwALE = zero
             flowDoms(nn,level,sps)%fwALE = zero
          end if

          ! Allocate the memory for the zeroth runge kutta stage
          allocate(flowDoms(nn,level,sps)%wn(2:il,2:jl,2:kl,1:nwf), &
               flowDoms(nn,level,sps)%pn(2:il,2:jl,2:kl), stat=ierr)
          if(ierr /= 0)                              &
               call terminate("allocMemFlovarPart2", &
               "Memory allocation failure for wn and pn")

          ! For unsteady mode using Runge-Kutta schemes allocate the
          ! memory for dwOldRK.

          if(equationMode          == unsteady .and. &
               timeIntegrationScheme == explicitRK) then

             mm = nRKStagesUnsteady - 1
             allocate(flowDoms(nn,level,sps)%dwOldRK(mm,il,jl,kl,nw), &
                  stat=ierr)
             if(ierr /= 0) &
                  call terminate("allocMemFlovarPart2", &
                  "Memory allocation failure for dwOldRK.")
          endif

       else fineLevelTest

          ! Coarser level. Allocate the memory for the multigrid
          ! forcing term and the state variables upon entry.

          allocate(flowDoms(nn,level,sps)%p1(1:ie,1:je,1:ke),          &
               flowDoms(nn,level,sps)%w1(1:ie,1:je,1:ke,1:nwf), &
               flowDoms(nn,level,sps)%wr(2:il,2:jl,2:kl,1:nwf), &
               stat=ierr)
          if(ierr /= 0)                              &
               call terminate("allocMemFlovarPart2", &
               "Memory allocation failure for p1, w1 &
               &and wr")

          ! Initialize w1 and p1 to zero, just that they
          ! are initialized.

          flowDoms(nn,level,sps)%p1 = zero
          flowDoms(nn,level,sps)%w1 = zero

       endif fineLevelTest

    enddo domains

  end subroutine allocMemFlovarPart2

  subroutine allocRestartFiles(nFiles)
    !
    !          Allocate memory for the restartfles                         *
    !          The array is populated from Python using setRestartFiles
    !          If memeory has been allocated for the array there exist at
    !          least one element in the array.
    !
    use constants
    use inputIO, only : restartFiles
    use utils, only : terminate
    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType) :: nFiles
    !
    !      Local variables.
    !
    integer :: ierr

    if (allocated(restartFiles)) then
       deallocate(restartFiles)
    end if

    allocate(restartFiles(nFiles), stat=ierr)
    if(ierr /= 0)                          &
         call terminate("allocRestartFiles", &
         "Memory allocation failure for restartFiles")

    ! Zero Array with empty strings
    restartFiles = ""

  end subroutine allocRestartFiles

  subroutine copySpectralSolution
    !
    !       copySpectralSolution copies the solution of the 1st spectral
    !       solution to all spectral solutions. This typically occurs when
    !       a for the spectral mode a restart is made from a steady or an
    !       unsteady solution. Possible rotation effects are taken into
    !       account for the velocity components.
    !
    use constants
    use block, only : flowDoms, nDom
    use flowVarRefState, only : nw
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use IOModule, only : IOVar
    use monitor, only: timeUnsteadyRestart
    use section, only : sections, nSections
    use utils, only : rotMatrixRigidBody
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: sps, spsm1, mm, nn, i, j, k, l

    real(kind=realType) :: dt, tnew, told, tmp
    real(kind=realType) :: theta, cosTheta, sinTheta

    real(kind=realType), dimension(3)   :: rotPoint, uu
    real(kind=realType), dimension(3)   :: xt, yt, zt
    real(kind=realType), dimension(3,3) :: rotMat

    real(kind=realType), dimension(nSections,3,3) :: rotMatSec

    ! Determine the rotation matrix from one spectral solution to the
    ! other for every section.

    sectionLoop: do nn=1,nSections

       ! Test if the section is rotating.

       testRotating: if( sections(nn)%rotating ) then

          ! Section is rotating. Determine the angle between the
          ! spectral solutions and its sine and cosine.

          i = sections(nn)%nSlices*nTimeIntervalsSpectral
          theta     = two*pi/real(i,realType)
          cosTheta = cos(theta)
          sinTheta = sin(theta)

          ! Transform to a frame where the xt-axis points in the
          ! direction of the rotation vector.

          xt(1) = sections(nn)%rotAxis(1)
          xt(2) = sections(nn)%rotAxis(2)
          xt(3) = sections(nn)%rotAxis(3)

          ! Construct the yt axis. It does not matter exactly as long
          ! as it is normal to xt.

          if(abs(xt(2)) < 0.707107_realType) then
             yt(1) = zero
             yt(2) = one
             yt(3) = zero
          else
             yt(1) = zero
             yt(2) = zero
             yt(3) = one
          endif

          ! Make sure that yt is normal to xt.

          tmp   = xt(1)*yt(1) + xt(2)*yt(2) + xt(3)*yt(3)
          yt(1) = yt(1) - tmp*xt(1)
          yt(2) = yt(2) - tmp*xt(2)
          yt(3) = yt(3) - tmp*xt(3)

          ! And create a unit vector.

          tmp   = one/sqrt(yt(1)**2 + yt(2)**2 + yt(3)**2)
          yt(1) = tmp*yt(1)
          yt(2) = tmp*yt(2)
          yt(3) = tmp*yt(3)

          ! Create the vector zt by taking the cross product xt*yt.

          zt(1) = xt(2)*yt(3) - xt(3)*yt(2)
          zt(2) = xt(3)*yt(1) - xt(1)*yt(3)
          zt(3) = xt(1)*yt(2) - xt(2)*yt(1)

          ! The rotation matrix in the xt,yt,zt frame is given by
          !
          ! R = | 1      0           0      |
          !     | 0  cos(theta) -sin(theta) |
          !     | 0  sin(theta)  cos(theta) |
          !
          ! The rotation matrix in the standard cartesian frame is then
          ! given by t * r * t^t, where the colums of the transformation
          ! matrix t are the unit vectors xt,yt,zt. One can easily check
          ! this by checking rotation around the y- and z-axis. The
          ! result of this is the expression below.

          rotMatSec(nn,1,1) = xt(1)*xt(1)                          &
               + cosTheta*(yt(1)*yt(1) + zt(1)*zt(1))
          rotMatSec(nn,1,2) = xt(1)*xt(2)                          &
               + cosTheta*(yt(1)*yt(2) + zt(1)*zt(2)) &
               - sinTheta*(yt(1)*zt(2) - yt(2)*zt(1))
          rotMatSec(nn,1,3) = xt(1)*xt(3)                          &
               + cosTheta*(yt(1)*yt(3) + zt(1)*zt(3)) &
               - sinTheta*(yt(1)*zt(3) - yt(3)*zt(1))

          rotMatSec(nn,2,1) = xt(1)*xt(2)                          &
               + cosTheta*(yt(1)*yt(2) + zt(1)*zt(2)) &
               + sinTheta*(yt(1)*zt(2) - yt(2)*zt(1))
          rotMatSec(nn,2,2) = xt(2)*xt(2)                          &
               + cosTheta*(yt(2)*yt(2) + zt(2)*zt(2))
          rotMatSec(nn,2,3) = xt(2)*xt(3)                          &
               + cosTheta*(yt(2)*yt(3) + zt(2)*zt(3)) &
               - sinTheta*(yt(2)*zt(3) - yt(3)*zt(2))

          rotMatSec(nn,3,1) = xt(1)*xt(3)                          &
               + cosTheta*(yt(1)*yt(3) + zt(1)*zt(3)) &
               + sinTheta*(yt(1)*zt(3) - yt(3)*zt(1))
          rotMatSec(nn,3,2) = xt(2)*xt(3)                          &
               + cosTheta*(yt(2)*yt(3) + zt(2)*zt(3)) &
               + sinTheta*(yt(2)*zt(3) - yt(3)*zt(2))
          rotMatSec(nn,3,3) = xt(3)*xt(3)                          &
               + cosTheta*(yt(3)*yt(3) + zt(3)*zt(3))

       else testRotating

          ! Section is not rotating. Set the rotation matrix to the
          ! identity matrix.

          rotMatSec(nn,1,1) = one
          rotMatSec(nn,1,2) = zero
          rotMatSec(nn,1,3) = zero

          rotMatSec(nn,2,1) = zero
          rotMatSec(nn,2,2) = one
          rotMatSec(nn,2,3) = zero

          rotMatSec(nn,3,1) = zero
          rotMatSec(nn,3,2) = zero
          rotMatSec(nn,3,3) = one

       endif testRotating

    enddo sectionLoop

    ! Initialize told to timeUnsteadyRestart. This takes the
    ! possibility into account that the spectral mode is restarted
    ! from an unsteady computation. Although not likely this
    ! possibility is allowed and should therefore be taken into
    ! account. Anyway told corresponds to the time of the 1st
    ! spectral solution. Also determine the time step between the
    ! spectral solutions. Both told and dt are only used to determine
    ! the rigid body motion and if these are specified it is assumed
    ! that there is only one section present in the grid.

    told = timeUnsteadyRestart
    dt   = sections(1)%timePeriod &
         / real(nTimeIntervalsSpectral,realType)

    ! Loop over the number of spectral modes, starting at 2.

    spectralLoop: do sps=2,nTimeIntervalsSpectral

       ! Determine the corresponding time for this spectral solution
       ! and store sps - 1 a bit easier.

       tnew  = told + dt
       spsm1 = sps - 1

       ! Determine the rotation matrix and rotation point between the
       ! told and tnew for the rigid body rotation of the entire mesh.
       ! The rotation point is not needed for the transformation of the
       ! velocities, but rotMatrixRigidBody happens to compute it.

       call rotMatrixRigidBody(tnew, told, rotMat, rotPoint)

       ! Loop over the local number of blocks.

       domains: do nn=1,nDom

          ! Store the section ID of this block a bit easier in mm.

          mm = flowDoms(nn,1,1)%sectionId

          ! Loop over the owned cells of this block. As the number of
          ! cells is identical for all spectral solutions, it does not
          ! matter which mode is taken for the upper dimensions.

          do k=2,flowDoms(nn,1,1)%kl
             do j=2,flowDoms(nn,1,1)%jl
                do i=2,flowDoms(nn,1,1)%il

                   ! Step 1. Copy the solution variables w from
                   ! the previous spectral solution.

                   do l=1,nw
                      IOVar(nn,sps)%w(i,j,k,l) = IOVar(nn,spsm1)%w(i,j,k,l)
                   enddo

                   ! Step 2. Apply the rigid body motion rotation matrix
                   ! to the velocity. Use uu as a temporary storage.

                   uu(1) = rotMat(1,1)*IOVar(nn,sps)%w(i,j,k,ivx) &
                        + rotMat(1,2)*IOVar(nn,sps)%w(i,j,k,ivy) &
                        + rotMat(1,3)*IOVar(nn,sps)%w(i,j,k,ivz)

                   uu(2) = rotMat(2,1)*IOVar(nn,sps)%w(i,j,k,ivx) &
                        + rotMat(2,2)*IOVar(nn,sps)%w(i,j,k,ivy) &
                        + rotMat(2,3)*IOVar(nn,sps)%w(i,j,k,ivz)

                   uu(3) = rotMat(3,1)*IOVar(nn,sps)%w(i,j,k,ivx) &
                        + rotMat(3,2)*IOVar(nn,sps)%w(i,j,k,ivy) &
                        + rotMat(3,3)*IOVar(nn,sps)%w(i,j,k,ivz)

                   ! Step 3. Apply the rotation matrix of the section to
                   ! the velocity.

                   IOVar(nn,sps)%w(i,j,k,ivx) = rotMatSec(mm,1,1)*uu(1) &
                        + rotMatSec(mm,1,2)*uu(2) &
                        + rotMatSec(mm,1,3)*uu(3)

                   IOVar(nn,sps)%w(i,j,k,ivy) = rotMatSec(mm,2,1)*uu(1) &
                        + rotMatSec(mm,2,2)*uu(2) &
                        + rotMatSec(mm,2,3)*uu(3)

                   IOVar(nn,sps)%w(i,j,k,ivz) = rotMatSec(mm,3,1)*uu(1) &
                        + rotMatSec(mm,3,2)*uu(2) &
                        + rotMatSec(mm,3,3)*uu(3)
                enddo
             enddo
          enddo

       enddo domains

       ! Set told to tnew for the next spectral solution.

       told = tnew

    enddo spectralLoop

  end subroutine copySpectralSolution

  subroutine determineSolFileNames
    !
    !       determineSolFileNames determines the number and names of the
    !       files that contain the solutions. For steady computations only
    !       one file must be present. For unsteady the situation is a
    !       little more complicated. It is attempted to read as many
    !       solutions as needed for a consistent restart. If not possible
    !       as many as possible solutions are read. For an unsteady
    !       computation the order will be reduced; for time spectral mode
    !       the solution will be interpolated.
    !
    use constants
    use communication, only : myID
    use inputIO, only : restartFiles
    use inputPhysics, only : equationMode
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use iteration, only : nOldSolAvail, oldSolWritten, nOldLevels
    use variableReading, only : solFiles, nSolsRead, interpolSpectral, copySpectral
    use utils, only : terminate
    implicit none
    !
    !      Local variables
    !
    integer :: ierr

    integer(kind=intType) :: ii, nn

    character(len=7)            :: integerString
    character(len=maxStringLen) :: tmpName

    ! Initialize copySpectral and interpolSpectral to .false.

    copySpectral     = .false.
    interpolSpectral = .false.

    ! Determine the desired number of files to be read. This depends
    ! on the equation mode we have to solve for. Also set the
    ! corresponding file names.

    select case(equationMode)

    case (steady)

       ! Steady computation. Only one solution needs to be read.
       ! In case a list of restart files were provided in the python
       ! script, we force a read from only the first solution file.

       nSolsRead = 1
       allocate(solFiles(nSolsRead), stat=ierr)
       if(ierr /= 0)                              &
            call terminate("determineSolFileNames", &
            "Memory allocation failure for solFiles")

       solFiles(1) = restartFiles(1)

       ! Check if the files can be opened, exit if that fails.
       call checkSolFileNames()

       !===============================================================

    case (unsteady)

       ! Unsteady computation. For a consistent restart nOldLevels
       ! solutions must be read. All restart files are provided explicitly
       ! from python script.
       call setSolFileNames()

       ! Check if the files can be opened, exit if that fails.
       call checkSolFileNames()

       ! Set nOldSolAvail to nSolsRead and check if a consistent
       ! restart can be made. If not, processor 0 prints a warning.

       nOldSolAvail = nSolsRead
       if(myID == 0 .and. nOldSolAvail < nOldLevels) then

          print "(a)", "#"
          print "(a)", "#                 Warning"
          print "(a)", "# Not enough data found for a consistent &
               &time accurate restart."
          print "(a)", "# Order is reduced in the first time steps &
               &until enough data is available again."
          print "(a)", "#"

       endif

       ! Set the logicals oldSolWritten, such that nothing is
       ! overwritten when solution files are dumped for this
       ! computation.

       ii = min(nSolsRead,nOldLevels-1)
       do nn=1,ii
          oldSolWritten(nn) = .true.
       enddo

       !===============================================================

    case (timeSpectral)

       ! Time spectral computation. For a consistent restart
       ! nTimeIntervalsSpectral solutions must be read. First
       ! determine the the restart files.
       call setSolFileNames()

       ! Check if the files can be opened, exit if that fails.
       call checkSolFileNames()

       ! Check whether or not the spectral solution must be copied
       ! or interpolated.

       if(nSolsRead == 1) then
          copySpectral = .true.
       else if(nSolsRead /= nTimeIntervalsSpectral) then
          interpolSpectral = .true.
       endif

    end select

  end subroutine determineSolFileNames

  subroutine setSolFileNames
    !
    !       setSolFileNames allocates and set the solution files that
    !       will be read and loaded in the restart
    !
    use communication
    use inputIO
    use utils, only : terminate
    use variableReading, only : solFiles, nSolsRead
    implicit none
    !
    !      Local variables
    !
    integer :: ierr
    integer(kind=intType) :: nn


    ! The length of the array provided gives the number of nSolsRead.
    nSolsRead = SIZE(restartFiles,1)

    ! Allocate the memory for the file names and set them.
    allocate(solFiles(nSolsRead), stat=ierr)
    if(ierr /= 0)                             &
         call terminate("determineSolFileNames", &
         "Memory allocation failure for solFiles")

    do nn=1,nSolsRead
       solFiles(nn) = restartFiles(nn)
    enddo

  end subroutine setSolFileNames



  subroutine checkSolFileNames
    !
    !       checkSolFileNames will check if the provided restart files
    !       are readable on disk. If not readable return fail, if readable
    !       message will be printed to let the user know that restart file
    !       will be tried to read.
    !
    use communication
    use utils, only : terminate
    use variableReading, only : solFiles, nSolsRead
    implicit none
    !
    !      Local variables
    !
    integer :: ierr
    character(len=maxStringLen) :: errorMessage
    integer(kind=intType) :: nn

    do nn=1,nSolsRead
       open(unit=21,file=solFiles(nn),status="old",iostat=ierr)
       if(ierr /= 0) then
          write(errorMessage,*) "Restart file ", trim(solFiles(nn)), &
               " could not be opened for reading"
          call terminate("checkSolFileNames", errorMessage)
          exit
       end if
       close(unit=21)

       if(myID == 0) then
          print "(a)", "#"
          write (*,100) trim(solFiles(nn))
          print "(a)", "#"
100       format("# Found restart file: ", A, 1X)
       end if
    enddo

  end subroutine checkSolFileNames

  subroutine initDepvarAndHalos(halosRead)
    !
    !       InitDepvarAndHalos computes the dependent flow variables,
    !       like viscosities, and initializes the halo cells by applying
    !       the boundary conditions and exchanging the internal halo's.
    !       This is all done on the start level grid.
    !
    use blockPointers
    use flowVarRefState
    use inputIO
    use inputIteration
    use inputPhysics
    use inputTimeSpectral
    use iteration
    use monitor
    use section
    use utils, only : setPointers
    use haloExchange, only : whalo2
    use turbUtils, only : computeEddyViscosity
    use turbBCRoutines, only :applyAllTurbBC
    use solverUtils, only : gridVelocitiesFineLevel, gridVelocitiesCoarseLevels, &
         normalVelocitiesAllLevels, slipVelocitiesFineLevel, slipVelocitiesCoarseLevels
    use flowUtils, only : computeLamViscosity
    use BCRoutines, only : applyAllBC
    use residuals, only : residual
    implicit none
    !
    !      Subroutine arguments.
    !
    logical, intent(in) :: halosRead
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm
    real(kind=realType)   :: relaxBleedsOr

    real(kind=realType), dimension(nSections) :: t

    logical :: initBleeds

    ! Set the logical whether or not to initialize the prescribed
    ! data for the bleed regions. If the halos were read the bleeds
    ! have been initialized already and nothing needs to be done.

    if( halosRead ) then
       initBleeds = .false.
    else
       initBleeds = .true.
    endif

    ! Compute the face velocities and for viscous walls the slip
    ! velocities. This is done for all the mesh levels.

    currentLevel = 1
    groundLevel  = 1

    do mm=1,nTimeIntervalsSpectral

       ! Compute the time, which corresponds to this spectral solution.
       ! For steady and unsteady mode this is simply the restart time;
       ! for the spectral mode the periodic time must be taken into
       ! account, which can be different for every section.

       t = timeUnsteadyRestart

       if(equationMode == timeSpectral) then
          do nn=1,nSections
             t(nn) = t(nn) + (mm-1)*sections(nn)%timePeriod &
                  /         real(nTimeIntervalsSpectral,realType)
          enddo
       endif
       call gridVelocitiesFineLevel(.false., t, mm)
       call gridVelocitiesCoarseLevels(mm)
       call normalVelocitiesAllLevels(mm)

       call slipVelocitiesFineLevel(.false., t, mm)
       call slipVelocitiesCoarseLevels(mm)

    enddo

    ! Loop over the number of spectral solutions and blocks
    ! to compute the laminar viscosity.

    do mm=1,nTimeIntervalsSpectral
       do nn=1,nDom
          call setPointers(nn,mgStartlevel,mm)
          call computeLamViscosity(.False.)
       enddo
    enddo

    ! Exchange the solution on the multigrid start level.
    ! It is possible that the halo values are needed for the boundary
    ! conditions. Viscosities are not exchanged.

    call whalo2(mgStartlevel, 1_intType, nw, .true., .true., &
         .false.)

    ! Apply all flow boundary conditions to be sure that the halo's
    ! contain the correct values. These might be needed to compute
    ! the eddy-viscosity. Also the data for the outflow bleeds
    ! is determined.

    currentLevel = mgStartlevel
    groundLevel  = mgStartlevel

    call applyAllBC(.true.)

    ! Loop over the number of spectral solutions and blocks
    ! to compute the eddy viscosities.

    do mm=1,nTimeIntervalsSpectral
       do nn=1,nDom

          ! Set the pointers for this block.

          call setPointers(nn,mgStartlevel,mm)

          ! Compute the eddy viscosity for rans computations using
          ! an eddy viscosity model.

          call computeEddyViscosity(.False.)

          ! In case of a rans computation and no restart, initialize
          ! the turbulent variables a bit better for some turbulence
          ! models.

          ! if(equations == RANSEquations .and. (.not. restart)) then
          !
          !   select case (turbModel)
          !
          !     case (komegaWilcox, komegaModified, menterSST)
          !       call initKOmega(0_intType)
          !       call computeEddyViscosity
          !
          !   end select
          !
          ! endif

       enddo
    enddo

    ! Exchange the laminar and eddy viscosities.

    call whalo2(mgStartlevel, 1_intType, 0_intType, .false., &
         .false., .true.)

    if (equations == RANSEquations) then
       call applyAllTurbBC(.true.)
    end if
    call applyAllBC(.true.)

    ! Exchange the solution for the second time to be sure that all
    ! halo's are initialized correctly. As this is the initialization
    ! phase, this is not critical.

    call whalo2(mgStartlevel, 1_intType, nw, .true., .true., &
         .true.)

  end subroutine initDepvarAndHalos

  subroutine initFlowRestart
    !
    !       initFlowRestart loads restart information from the restart
    !       file into the state variables.
    !
    use constants
    use IOModule, only : IOVar
    use variableReading, only : solFiles, interpolSpectral, copySpectral, &
         halosRead
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: ierr

    ! Initialize halosRead to .false. This will be overwritten
    ! if halo values are read during the restart.

    halosRead = .false.

    ! Determine the number and names of the solution files and
    ! determine the contents of IOVar, the data structure to
    ! generalize the IO.

    call determineSolFileNames
    call setIOVar

    ! Determine the format of the files and read them.
    ! Note that halosRead is possibly overwritten in the
    ! folloing select case statement below

    call readRestartFile()

    ! Copy or interpolate the spectral solution, if needed.

    if( copySpectral )     call copySpectralSolution
    if( interpolSpectral ) call interpolateSpectralSolution

    ! Release the memory of the file names and IOVar.

    deallocate(solFiles, IOVar, stat=ierr)
    if(ierr /= 0)                &
         call terminate("initFlow", &
         "Deallocation failure for solFiles and IOVar")

    ! At the moment the pressure is stored at the location of the
    ! total energy. Copy the pressure to its own arrays and
    ! compute the total energy.

    call setPressureAndComputeEnergy(halosRead)

    ! Initialize the halo cells if a restart is performed and the
    ! entire flow field if this is not the case.

    call initializeHalos(halosRead)

    ! Initialize the dependent flow variables and the halo values.
    call initDepvarAndHalos(halosRead)

  end subroutine initFlowRestart

  subroutine initFlowfield
    !
    !       initFlowfield initializes the flow field to a uniform flow on
    !       the start level grid. Exception may be some turbulence
    !       variables, which are initialized a bit smarter.
    !
    use constants
    use communication, only : myID
    use inputIteration, only : nCycles, nSGStartup
    use inputPhysics, only : equationMode
    use inputUnsteady, only: nTimeStepsFine
    use iteration, only : nOldSolAvail
    use monitor, only : nTimeStepsRestart, timeUnsteadyRestart
    use utils, only : allocConvArrays, allocTimeArrays
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn

    ! Initialize nTimeStepsRestart to 0 (no restart
    ! is performed) and allocate the memory for the arrays to store
    ! convergence history. This allocation is only to be done by
    ! processor 0. For an unsteady computation the entire convergence
    ! history is stored and the values after every time step is
    ! stored in a separate array.

    nTimeStepsRestart   = 0
    timeUnsteadyRestart = zero

    nn = nsgStartup + nCycles
    if(equationMode == unsteady) nn = nTimeStepsFine*nn
    if(myID == 0) call allocConvArrays(nn)
    if(equationMode == unsteady .and. myId == 0) &
         call allocTimeArrays(nTimeStepsFine)

    ! Set nOldSolAvail to 1, to indicate that an unsteady
    ! computation should be started with a lower order
    ! time integration scheme, because not enough states
    ! in the past are available.

    nOldSolAvail = 1

    ! Initialize the flow field to uniform flow.

    call setUniformFlow

  end subroutine initFlowfield

  subroutine initializeHalos(halosRead)
    !
    !       initializeHalos sets the flow variables in the halo cells
    !       using a constant extrapolation. If the halos are read only the
    !       second halos are initialized, otherwise both.
    !
    use constants
    use blockPointers, only : nDom, w, p, rlv, ib, jb, kb, &
         il, ie, jl, je, kl, ke, rev
    use flowVarRefState, only : viscous, nw, eddyModel, muInf
    use inputIteration, only : mgStartLevel
    use inputPhysics, only : eddyVisInfRatio
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    logical, intent(in) :: halosRead
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm, i, j, k, l
    integer(kind=intType) :: jj, kk

    ! Loop over the number of spectral solutions and blocks.

    spectralLoop: do mm=1,nTimeIntervalsSpectral
       domains: do nn=1,nDom

          ! Set the pointers for this block.

          call setPointers(nn,mgStartlevel,mm)

          ! Determine the situation we are dealing with.

          testHalosRead: if( halosRead ) then

             ! The first layer of halo cells have been read. Initialize
             ! the second layer from the first layer.

             ! Halo cells in the i-direction.

             do k=0,kb
                kk = max(1_intType,min(k,ke))
                do j=0,jb
                   jj = max(1_intType,min(j,je))

                   do l=1,nw
                      w(0, j,k,l) = w(1, jj,kk,l)
                      w(ib,j,k,l) = w(ie,jj,kk,l)
                   enddo

                   p(0, j,k) = p(1, jj,kk)
                   p(ib,j,k) = p(ie,jj,kk)

                enddo
             enddo

             ! Halo cells in j-direction. Note that the i-halo's have
             ! already been set.

             do k=0,kb
                kk = max(1_intType,min(k,ke))
                do i=1,ie

                   do l=1,nw
                      w(i,0, k,l) = w(i,1, kk,l)
                      w(i,jb,k,l) = w(i,je,kk,l)
                   enddo

                   p(i,0 ,k) = p(i,1, kk)
                   p(i,jb,k) = p(i,je,kk)

                enddo
             enddo

             ! Halo cells in k-direction. Note that the halo's in both
             ! i and j direction have already been set.

             do j=1,je
                do i=1,ie

                   do l=1,nw
                      w(i,j,0, l) = w(i,j,1, l)
                      w(i,j,kb,l) = w(i,j,ke,l)
                   enddo

                   p(i,j,0)  = p(i,j,1)
                   p(i,j,kb) = p(i,j,ke)

                enddo
             enddo

          else testHalosRead

             ! No halo cells have been read. Initialize both layers
             ! using the internal value.

             ! Halo cells in the i-direction.

             do k=0,kb
                kk = max(2_intType,min(k,kl))
                do j=0,jb
                   jj = max(2_intType,min(j,jl))

                   do l=1,nw
                      w(0, j,k,l) = w(2, jj,kk,l)
                      w(1, j,k,l) = w(2, jj,kk,l)
                      w(ie,j,k,l) = w(il,jj,kk,l)
                      w(ib,j,k,l) = w(il,jj,kk,l)
                   enddo

                   p(0, j,k) = p(2, jj,kk)
                   p(1, j,k) = p(2, jj,kk)
                   p(ie,j,k) = p(il,jj,kk)
                   p(ib,j,k) = p(il,jj,kk)

                enddo
             enddo

             ! Halo cells in j-direction. Note that the i-halo's have
             ! already been set.

             do k=0,kb
                kk = max(2_intType,min(k,kl))
                do i=2,il

                   do l=1,nw
                      w(i,0, k,l) = w(i,2, kk,l)
                      w(i,1, k,l) = w(i,2, kk,l)
                      w(i,je,k,l) = w(i,jl,kk,l)
                      w(i,jb,k,l) = w(i,jl,kk,l)
                   enddo

                   p(i,0 ,k) = p(i,2, kk)
                   p(i,1 ,k) = p(i,2, kk)
                   p(i,je,k) = p(i,jl,kk)
                   p(i,jb,k) = p(i,jl,kk)

                enddo
             enddo

             ! Halo cells in k-direction. Note that the halo's in both
             ! i and j direction have already been set.

             do j=2,jl
                do i=2,il

                   do l=1,nw
                      w(i,j,0, l) = w(i,j,2, l)
                      w(i,j,1, l) = w(i,j,2, l)
                      w(i,j,ke,l) = w(i,j,kl,l)
                      w(i,j,kb,l) = w(i,j,kl,l)
                   enddo

                   p(i,j,0)  = p(i,j,2)
                   p(i,j,1)  = p(i,j,2)
                   p(i,j,ke) = p(i,j,kl)
                   p(i,j,kb) = p(i,j,kl)

                enddo
             enddo

          endif testHalosRead

          ! Initialize the laminar and eddy viscosity, if appropriate,
          ! such that no uninitialized memory is present.
          ! As the viscosities are dependent variables their values
          ! are not read during a restart.

          if( viscous )   rlv = muInf
          if( eddyModel ) rev = eddyVisInfRatio*muInf

       enddo domains
    enddo spectralLoop

  end subroutine initializeHalos
  subroutine interpolateSpectralSolution
    !
    !       interpolateSpectralSolution uses a spectral interpolation to
    !       determine the initialization of the flow solution.
    !       The solution is interpolated from the solution read, which
    !       contains a different number of time instances and is stored in
    !       IOVar()%w. This variable can be found in IOModule.
    !
    use constants
    use blockPointers, only : w, il, jl, kl, nDom, sectionID
    use flowVarRefState, only : nw
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use IOModule, only : IOVar
    use section, only: sections, nSections
    use utils, only : setPointers, terminate, spectralInterpolCoef
    use variableReading, only : nSolsRead
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: jj, nn, ll, sps, i, j, k, l

    real(kind=realType) :: t

    real(kind=realType), dimension(nSolsRead) :: alpScal
    real(kind=realType), dimension(nSections,nSolsRead,3,3) :: alpMat

    ! Loop over the number of spectral solutions to be interpolated.

    spectralLoop: do sps=1,nTimeIntervalsSpectral

       ! Determine the ratio of the time of this solution and the
       ! periodic time.

       nn = sps - 1
       t  = real(nn,realType)/real(nTimeIntervalsSpectral,realType)

       ! Determine the interpolation coefficients for both the scalar
       ! and the vector quantities.

       call spectralInterpolCoef(nSolsRead, t, alpScal, alpMat)

       ! Loop over the local number of blocks.

       domains: do nn=1,nDom

          ! Set the pointers for this block to the finest grid level.

          call setPointers(nn, 1_intType, sps)

          ! Loop over the number of variables to be interpolated.

          varLoop: do l=1,nw

             ! Check if this is a velocity variable.

             velTest: if(l == ivx .or. l == ivy .or. l == ivz) then

                ! Velocity variable. Set ll, which is the row in the
                ! matrix coefficients.

                select case (l)
                case (ivx)
                   ll = 1
                case (ivy)
                   ll = 2
                case (ivz)
                   ll = 3
                end select

                ! Loop over the owned cells to interpolate the variable.

                do k=2,kl
                   do j=2,jl
                      do i=2,il

                         ! Initialization to zero and loop over the number of
                         ! spectral solutions used in the interpolation and
                         ! update the variable accordingly. Note that for the
                         ! vector variables the matrix coefficients must be
                         ! used; these matrices can be different for the
                         ! sections present in the grid.

                         w(i,j,k,l) = zero
                         do jj=1,nSolsRead
                            w(i,j,k,l) = w(i,j,k,l)                &
                                 + alpMat(sectionID,jj,ll,1) &
                                 * IOVar(nn,jj)%w(i,j,k,ivx) &
                                 + alpMat(sectionID,jj,ll,2) &
                                 * IOVar(nn,jj)%w(i,j,k,ivy) &
                                 + alpMat(sectionID,jj,ll,3) &
                                 * IOVar(nn,jj)%w(i,j,k,ivz)
                         enddo

                      enddo
                   enddo
                enddo

             else velTest

                ! Scalar variable.
                ! Loop over the owned cells to interpolate the variable.

                do k=2,kl
                   do j=2,jl
                      do i=2,il

                         ! Initialization to zero and loop over the number of
                         ! spectral solutions used in the interpolation and
                         ! update the variable accordingly.

                         w(i,j,k,l) = zero
                         do jj=1,nSolsRead
                            w(i,j,k,l) = w(i,j,k,l)  &
                                 + alpScal(jj)*IOVar(nn,jj)%w(i,j,k,l)
                         enddo

                      enddo
                   enddo
                enddo

             endif velTest

          enddo varLoop

       enddo domains

    enddo spectralLoop

    ! Release the memory of w of IOVar.

    do sps=1,nSolsRead
       do nn=1,nDom
          deallocate(IOVar(nn,sps)%w, stat=ierr)
          if(ierr /= 0) &
               call terminate("interpolateSpectralSolution", &
               "Deallocation failure for w.")
       enddo
    enddo

  end subroutine interpolateSpectralSolution

  subroutine releaseExtraMemBCs
    !
    !       releaseExtraMemBCs releases the extra memory allocated in
    !       allocMemBcdata. This additional memory was allocated, such
    !       that alternative boundary condition treatments can be handled
    !       in setBCDataFineGrid.
    !
    use constants
    use blockPointers, only : flowDoms, nDom, BCType, BCData, nBocos
    use inputTimeSpectral, only: nTimeIntervalsSpectral
    use utils, only : setPointers, terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: mm, nn, sps, level, nLevels, ii

    ! Determine the number of multigrid levels.

    nLevels = ubound(flowDoms,2)

    ! Loop over the number of multigrid level, spectral solutions
    ! and local blocks.

    levelLoop: do level=1,nLevels
       spectralLoop: do sps=1,nTimeIntervalsSpectral
          domainsLoop: do nn=1,nDom

             ! Have the pointers in blockPointers point to the
             ! current block to make everything more readable.

             call setPointers(nn, level, sps)

             ! Loop over the number of boundary subfaces for this block.

             bocoLoop: do mm=1,nBocos

                ! Determine the boundary condition we are having here.

                inflowType: select case (BCType(mm))

                case (SubsonicInflow)

                   ! Subsonic inflow. Determine the boundary condition
                   ! treatment and release the accordingly. Note that
                   ! the boundary condition treatment of the finest mesh
                   ! must be used, because this data is not available yet
                   ! on the coarse grid.

                   ii = &
                        flowDoms(nn,1,sps)%BCData(mm)%subsonicInletTreatment

                   select case (ii)

                   case (totalConditions)

                      ! Total conditions used. Release the memory for
                      ! the density and velocities and nullify their
                      ! pointers.

                      deallocate(BCData(mm)%rho,  BCData(mm)%velx, &
                           BCData(mm)%vely, BCData(mm)%velz, &
                           stat=ierr)
                      if(ierr /= 0) &
                           call terminate("releaseExtraMemBCs", &
                           "Deallocation failure for rho, &
                           &velx, vely and velz")

                      nullify(BCData(mm)%rho)
                      nullify(BCData(mm)%velx)
                      nullify(BCData(mm)%vely)
                      nullify(BCData(mm)%velz)
                      !===================================================

                   case (massFlow)

                      ! Full velocity vector and density prescribed at
                      ! inlet boundaries. Release the memory for the
                      ! total conditions and flow directions and nullify
                      ! their pointers.

                      deallocate(BCData(mm)%ptInlet,        &
                           BCData(mm)%ttInlet,        &
                           BCData(mm)%htInlet,        &
                           BCData(mm)%flowXdirInlet, &
                           BCData(mm)%flowYdirInlet, &
                           BCData(mm)%flowZdirInlet, stat=ierr)
                      if(ierr /= 0) &
                           call terminate("releaseExtraMemBCs", &
                           "Deallocation failure for the &
                           &total conditions.")

                      nullify(BCData(mm)%ptInlet)
                      nullify(BCData(mm)%ttInlet)
                      nullify(BCData(mm)%htInlet)
                      nullify(BCData(mm)%flowXdirInlet)
                      nullify(BCData(mm)%flowYdirInlet)
                      nullify(BCData(mm)%flowZdirInlet)

                   end select

                end select inflowType

             enddo bocoLoop
          enddo domainsLoop
       enddo spectralLoop
    enddo levelLoop

  end subroutine releaseExtraMemBCs

  subroutine setIOVar
    !
    !       setIOVar allocates the memory for the derived data type IOVar,
    !       which is simplifies the reading. If an interpolation must be
    !       performed for the time spectral method also the solution of
    !       this IO type is allocated. For all other cases the pointers of
    !       IOVar are set the the appropriate entries of flowDoms, with
    !       possible offset due to the usage of pointers.
    !
    use constants
    use block, only : flowDoms, nDom
    use flowVarRefState, only : nw
    use inputPhysics, only : equationMode
    use IOModule, only : IOVar
    use utils, only : terminate
    use variableReading, only : interpolSpectral, halosRead, nSolsRead
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, mm, il, jl, kl

    ! Allocate the memory for IOVar.

    allocate(IOVar(nDom,nSolsRead), stat=ierr)
    if(ierr /= 0)                &
         call terminate("setIOVar", &
         "Memory allocation failure for solRead")

    ! Determine the equation mode we are solving and set the pointers
    ! of coorRead accordingly, or even allocate the memory, if needed.

    select case(equationMode)

    case (steady)

       ! Steady computation. Only one solution needs to be read.
       ! Loop over the number of blocks and set the pointers.
       ! No pointer offset is needed.

       do nn=1,nDom
          IOVar(nn,1)%pointerOffset = 0
          IOVar(nn,1)%w => flowDoms(nn,1,1)%w(1:,1:,1:,:)
       enddo

       !===============================================================

    case (unsteady)

       ! Unsteady computation. The first solution should be stored in
       ! w. For the others the pointers point to wOld. As the
       ! starting indices of wOld are 2, a pointer shift takes place
       ! here. I know this is a pain in the butt, but that's what
       ! we have to live with.

       do nn=1,nDom
          IOVar(nn,1)%pointerOffset = 0
          IOVar(nn,1)%w => flowDoms(nn,1,1)%w(1:,1:,1:,:)

          do mm=2,nSolsRead
             IOVar(nn,mm)%pointerOffset = -1
             IOVar(nn,mm)%w => flowDoms(nn,1,1)%wOld(mm-1,2:,2:,2:,:)
          enddo
       enddo

       !===============================================================

    case (timeSpectral)

       ! Time spectral mode. A further check is required.

       testAllocSolRead: if( interpolSpectral ) then

          ! A restart is performed using a different number of time
          ! instances than the previous computation. Consequently the
          ! solutions will be interpolated later on. Hence some
          ! additional storage is required for the solutions and thus
          ! the w variables of IOVar are allocated. No halo data is
          ! needed here. No pointer offset either due to the explicit
          ! allocation.

          do nn=1,nDom
             il = flowDoms(nn,1,1)%il
             jl = flowDoms(nn,1,1)%jl
             kl = flowDoms(nn,1,1)%kl

             do mm=1,nSolsRead
                IOVar(nn,mm)%pointerOffset = 0

                allocate(IOVar(nn,mm)%w(2:il,2:jl,2:kl,nw), stat=ierr)
                if(ierr /= 0)                &
                     call terminate("setIOVar", &
                     "Memory allocation failure for w")
             enddo

          enddo

       else testAllocSolRead

          ! A restart is made using either 1 solution or the correct
          ! number of solution instances. In both cases simply set
          ! the pointers of IOVar. No pointer offset is needed.

          do nn=1,nDom
             do mm=1,nSolsRead
                IOVar(nn,mm)%pointerOffset = 0
                IOVar(nn,mm)%w => flowDoms(nn,1,mm)%w(1:,1:,1:,:)
             enddo
          enddo

       endif testAllocSolRead

    end select

  end subroutine setIOVar

  subroutine setPressureAndComputeEnergy(halosRead)
    !
    !       Due to the usage of the variable IOVar, which generalizes the
    !       IO and leads to reuse of code, currently the pressure is
    !       stored at the position of rhoE. In this routine that data is
    !       copied to the pressure array and the total energy is computed.
    !       Note that this routine is only called when a restart is done.
    !
    use constants
    use blockPointers, only : nDom, p, w, il, jl, kl
    use flowVarRefState, only : kPresent
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use utils, only : setPointers
    use flowUtils, only : computeEtotBlock
    implicit none
    !
    !      Subroutine arguments.
    !
    logical, intent(in) :: halosRead
    !
    !      Local variables.
    !
    integer(kind=intType) :: sps, nn, nHalo
    integer(kind=intType) :: i, j, k
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    ! Set the value of nHalo, depending whether or not the halo cells
    ! have been read from the restart file.

    nHalo = 0
    if( halosRead ) nHalo = 1

    ! Loop over the number of time instances and the local blocks.
    ! As this routine is only called when a restart is performed,
    ! the MG start level is the finest level.

    do sps=1,nTimeIntervalsSpectral
       do nn=1,nDom

          ! Set the pointers to the correct block. As this routine is
          ! only called when a restart is performed, the MG start level
          ! is the finest level.

          call setPointers(nn,1_intType,sps)

          ! Determine the range for which the pressure must be computed.

          iBeg = 2 -nHalo; jBeg =  2-nHalo; kBeg =  2-nHalo
          iEnd = il+nHalo; jEnd = jl+nHalo; kEnd = kl+nHalo

          ! Copy the pressure for the required cells.

          do k=kBeg,kEnd
             do j=jBeg,jEnd
                do i=iBeg,iEnd
                   p(i,j,k) = w(i,j,k,irhoE)
                enddo
             enddo
          enddo

          ! Compute the total energy as well.

          call computeEtotBlock(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd, kPresent)
       enddo
    enddo

  end subroutine setPressureAndComputeEnergy
  subroutine setRestartFiles(fileName, i)
    !
    !          Populates the restartfiles
    !          The array is populated from Python using setRestartFiles
    !
    use constants
    use inputIO, only : restartFiles
    implicit none
    !
    !      Subroutine argument.
    !
    character(len=*), intent(inout) :: fileName
    integer(kind=intType) :: i



    restartFiles(i) = fileName

  end subroutine setRestartFiles

  subroutine setUniformFlow
    !
    !       setUniformFlow set the flow variables of all local blocks on
    !       the start level to the uniform flow field.
    !
    use constants
    use blockPointers, only : w, dw, fw, flowDoms, ib, jb, kb, &
         rev, rlv, nDom, BCData, p
    use communication
    use flowVarRefState, only : eddyModel, viscous, muInf, nw, nwf, &
         pInfCorr, wInf
    use inputIteration, only : mgStartLevel
    use inputPhysics, only : equationMode, flowType, eddyVisInfRatio
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use utils, only : setPointers
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, mm, i, j, k, l

    real(kind=realType) :: tmp

    real(kind=realType), dimension(3) :: dirLoc, dirGlob

    ! Loop over the number of spectral solutions and blocks.
    spectralLoop: do mm=1,nTimeIntervalsSpectral
       domains: do nn=1,nDom

          ! Set the pointers for this block.

          call setPointers(nn,mgStartlevel,mm)

          ! Set the w-variables to the ones of the uniform flow field.
          do l=1,nw
             do k=0,kb
                do j=0,jb
                   do i=0,ib
                      w(i,j,k,l) = wInf(l)
                      dw(i,j,k,l) = zero
                   enddo
                enddo
             enddo
          enddo
          !set this here for a reinitialize flow to eliminate possible NAN's
          do l=1,nwf
             do k=0,kb
                do j=0,jb
                   do i=0,ib
                      fw(i,j,k,l) = zero
                   enddo
                enddo
             enddo
          enddo

          ! Set the pressure.

          p = pInfCorr

          ! Initialize the laminar and eddy viscosity, if appropriate,
          ! such that no uninitialized memory is present.

          if( viscous )   rlv = muInf
          if( eddyModel ) rev = eddyVisInfRatio*muInf

       enddo domains
    enddo spectralLoop

    ! Correct for the time spectral method in combination with an
    ! internal flow computation the velocity direction.
    ! It is possible that the prescribed direction is different
    ! for every time instance.

    testCorrection: if(equationMode == timeSpectral .and. &
         flowType     == internalFlow) then

       ! Loop over the number of spectral solutions.

       spectralLoopCorr: do mm=1,nTimeIntervalsSpectral

          ! Initialize the local direction to zero. In the loop
          ! below this direction will be accumulated.

          dirLoc = zero

          ! Loop over the number of blocks and determine the average
          ! velocity direction prescribed at the inlets.

          domainLoop1: do nn=1,nDom

             ! Set the pointer for the BCData to make the code more
             ! readable.

             BCData => flowDoms(nn,mgStartlevel,mm)%BCData

             ! Loop over the number of boundary subfaces and update
             ! the local prescribed velocity direction.

             do l=1,flowDoms(nn,mgStartlevel,mm)%nBocos
                call velMagnAndDirectionSubface(tmp, dirLoc, BCData, l)
             enddo
          enddo domainLoop1

          ! Determine the sum of dirLoc and create a unit vector
          ! for the global direction.

          call mpi_allreduce(dirLoc, dirGlob, 3, adflow_real, mpi_sum, &
               ADflow_comm_world, ierr)

          tmp         = one/max(eps,sqrt(dirGlob(1)**2 &
               +                  dirGlob(2)**2 &
               +                  dirGlob(3)**2))
          dirGlob(1) = tmp*dirGlob(1)
          dirGlob(2) = tmp*dirGlob(2)
          dirGlob(3) = tmp*dirGlob(3)

          ! Loop again over the local domains and correct the
          ! velocity direction.

          domainsLoop2: do nn=1,nDom

             ! Set the pointers for this block.

             call setPointers(nn,mgStartlevel,mm)

             do k=0,kb
                do j=0,jb
                   do i=0,ib
                      tmp          = sqrt(w(i,j,k,ivx)**2 &
                           +      w(i,j,k,ivy)**2 &
                           +      w(i,j,k,ivz)**2)
                      w(i,j,k,ivx) = tmp*dirGlob(1)
                      w(i,j,k,ivy) = tmp*dirGlob(2)
                      w(i,j,k,ivz) = tmp*dirGlob(3)
                   enddo
                enddo
             enddo

          enddo domainsLoop2
       enddo spectralLoopCorr

    endif testCorrection

  end subroutine setUniformFlow

  !=================================================================

  subroutine velMagnAndDirectionSubface(vmag, dir, BCData, mm)
    !
    !       VelMagnAndDirectionSubface determines the maximum value
    !       of the magnitude of the velocity as well as the sum of the
    !       flow directions for the currently active subface.
    !
    use constants
    use block
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: mm

    real(kind=realType), intent(out) :: vmag
    real(kind=realType), dimension(3), intent(inout) :: dir

    type(BCDataType), dimension(:), pointer :: BCData
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j
    real(kind=realType)   :: vel

    ! Initialize vmag to -1.0.

    vmag = -one

    ! Check if the velocity is prescribed.

    if( associated(BCData(mm)%velx) .and. &
         associated(BCData(mm)%vely) .and. &
         associated(BCData(mm)%velz) ) then

       ! Loop over the owned faces of the subface. As the cell range
       ! may contain halo values, the nodal range is used.

       do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
          do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd

             ! Compute the magnitude of the velocity and compare it
             ! with the current maximum. Store the maximum of the two.

             vel  = sqrt(BCData(mm)%velx(i,j)**2 &
                  +      BCData(mm)%vely(i,j)**2 &
                  +      BCData(mm)%velz(i,j)**2)
             vmag = max(vmag, vel)

             ! Compute the unit vector of the velocity and add it to dir.

             vel    = one/max(eps,vel)
             dir(1) = dir(1) + vel*BCData(mm)%velx(i,j)
             dir(2) = dir(2) + vel*BCData(mm)%vely(i,j)
             dir(3) = dir(3) + vel*BCData(mm)%velz(i,j)

          enddo
       enddo
    endif

    ! Check if the velocity direction is prescribed.

    if( associated(BCData(mm)%flowXdirInlet) .and. &
         associated(BCData(mm)%flowYdirInlet) .and. &
         associated(BCData(mm)%flowZdirInlet) ) then

       ! Add the unit vectors to dir by looping over the owned
       ! faces of the subfaces. Again the nodal range must be
       ! used for this.

       do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
          do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd

             dir(1) = dir(1) + BCData(mm)%flowXdirInlet(i,j)
             dir(2) = dir(2) + BCData(mm)%flowYdirInlet(i,j)
             dir(3) = dir(3) + BCData(mm)%flowZdirInlet(i,j)

          enddo
       enddo

    endif

  end subroutine velMagnAndDirectionSubface
  subroutine timeSpectralCoef(coefSpectral, matrixCoefSpectral, &
       diagMatCoefSpectral)
    !
    !       timeSpectralCoef computes the time integration coefficients
    !       for the time spectral method. As it is possible that sections
    !       have different periodic times these coefficients are
    !       determined for all the sections. For vector quantities, such
    !       as momentum, these coefficients can also be different due to
    !       rotation and the fact that only a part of the wheel is
    !       simulated.
    !
    use constants
    use flowVarRefState, only : timeRef
    use inputTimeSpectral, only : nTimeIntervalsSpectral, rotMatrixSpectral
    use section, only : nSections, sections
    implicit none
    !
    !      Subroutine arguments.
    !
    real(kind=realType),                               &
         dimension(nSections,nTimeIntervalsSpectral-1), &
         intent(out) :: coefSpectral
    real(kind=realType),                                    &
         dimension(nSections,nTimeIntervalsSpectral-1,3,3), &
         intent(out) :: matrixCoefSpectral
    real(kind=realType), dimension(nSections,3,3), &
         intent(out) :: diagMatCoefSpectral
    !
    !      Local variables.
    !
    integer(kind=intType) :: pp, nn, mm, ii, i, j, ntot
    real(kind=realType)   :: coef, dAngle, angle, fact, slicesFact

    real(kind=realType), dimension(3,3) :: rotMat, tmp

    ! Loop over the number of sections.

    sectionLoop: do mm=1,nSections

       ! Initialize dAngle (smallest angle in the cotangent function)
       ! and coef, which is the multiplication factor in front of the
       ! cotangent/cosecant function. Coef is a combination of the 1/2
       ! and the 2*pi/timePeriod

       dAngle = pi/real(nTimeIntervalsSpectral,realType)
       coef   = pi*timeRef/sections(mm)%timePeriod

       ! Computation of the scalar coefficients.

       scalarLoop: do nn=1,(nTimeIntervalsSpectral-1)

          angle = nn*dAngle

          ! The coefficient for an odd and even number of time
          ! instances are different; the former is 1/sin, the
          ! latter cos/sin or 1/tan.

          coefSpectral(mm,nn) = coef/sin(angle)

          if (mod(nTimeIntervalsSpectral,2_intType) == 0) &
               coefSpectral(mm,nn) = coefSpectral(mm,nn)*cos(angle)

          ! Negate coef for the next spectral coefficient.

          coef = -coef

       enddo scalarLoop

       ! Initialize dAngle to the smallest angle in the cotangent
       ! or cosecant function. Now this angle is for the entire wheel,
       ! i.e. the number of slices must be taken into account.

       ntot   = nTimeIntervalsSpectral*sections(mm)%nSlices
       dAngle = pi/real(ntot,realType)

       ! Initialize the rotation matrix to the unity matrix.

       rotMat(1,1) = one;  rotMat(1,2) = zero; rotMat(1,3) = zero
       rotMat(2,1) = zero; rotMat(2,2) = one;  rotMat(2,3) = zero
       rotMat(3,1) = zero; rotMat(3,2) = zero; rotMat(3,3) = one

       ! Loop over the number of spectral coefficient to initialize the
       ! matrix coefficients; this is basically pp == 0 in the loop
       ! over the number slices. Use is made of the fact that the
       ! rotation matrix is the identity for pp == 0.
       ! coef changes sign at every time instance

       slicesFact = one/real(sections(mm)%nSlices,realType)
       fact = one

       do nn=1,(nTimeIntervalsSpectral-1)

          ! Determine the scalar coefficient. This value depends now
          ! whether the total number of time instances in the wheel is
          ! odd or even.

          angle = nn*dAngle
          coef  = one/sin(angle)

          if (mod(ntot,2_intType) == 0) &
               coef = coef*cos(angle)

          coef = coef*fact*slicesFact

          ! The first part of matrixCoefSpectral is a diagonal matrix,
          ! because this indicates the contribution of the current
          ! slice to the time derivative.

          matrixCoefSpectral(mm,nn,1,1) = coef
          matrixCoefSpectral(mm,nn,1,2) = zero
          matrixCoefSpectral(mm,nn,1,3) = zero

          matrixCoefSpectral(mm,nn,2,1) = zero
          matrixCoefSpectral(mm,nn,2,2) = coef
          matrixCoefSpectral(mm,nn,2,3) = zero

          matrixCoefSpectral(mm,nn,3,1) = zero
          matrixCoefSpectral(mm,nn,3,2) = zero
          matrixCoefSpectral(mm,nn,3,3) = coef

          fact = -fact

       enddo

       ! Initialize diagMatCoefSpectral to zero, because the
       ! starting index in the loop over the number of slices -1 is
       ! 1, i.e. the slice where the actual computation takes places
       ! does not contribute to diagMatCoefSpectral.

       do j=1,3
          do i=1,3
             diagMatCoefSpectral(mm,i,j) = zero
          enddo
       enddo

       ! Loop over the additional slices which complete an entire
       ! revolution. To be able to compute the coefficients a bit
       ! easier the loop runs from 1 to nSlices-1 and not from
       ! 2 to nSlices.

       slicesLoop: do pp=1,(sections(mm)%nSlices-1)

          ! Compute the rotation matrix for this slice. This is the
          ! old one multiplied by the transformation matrix going from
          ! one slices to the next. Use tmp as temporary storage.

          do j=1,3
             do i=1,3
                tmp(i,j) = rotMatrixSpectral(mm,i,1)*rotMat(1,j) &
                     + rotMatrixSpectral(mm,i,2)*rotMat(2,j) &
                     + rotMatrixSpectral(mm,i,3)*rotMat(3,j)
             enddo
          enddo

          rotMat = tmp

          slicesFact = one/real(sections(mm)%nSlices,realType)

          ! Loop over the number of spectral coefficients and update
          ! matrixCoefSpectral. The multiplication with (-1)**nn
          ! takes place here too.

          ! Multiply also by the term (-1)**(pN+1)

          fact = one
          if (mod(pp*nTimeIntervalsSpectral,2_intType) /= 0) &
               fact = -one
          slicesFact = fact*slicesFact

          fact = one
          ii = pp*nTimeIntervalsSpectral
          do nn=1,(nTimeIntervalsSpectral-1)

             ! Compute the coefficient multiplying the rotation matrix.
             ! Again make a distinction between an odd and an even
             ! number of time instances for the entire wheel.

             angle = (nn+ii)*dAngle
             coef  = one/sin(angle)

             if (mod(ntot,2_intType) == 0) &
                  coef = coef*cos(angle)

             coef = coef*fact*slicesFact

             ! Update matrixCoefSpectral.

             do j=1,3
                do i=1,3
                   matrixCoefSpectral(mm,nn,i,j) = &
                        matrixCoefSpectral(mm,nn,i,j) + coef*rotMat(i,j)
                enddo
             enddo

             fact = -fact

          enddo

          ! Update diagMatCoefSpectral. Also here the distinction
          ! between odd and even number of time instances.

          angle = ii*dAngle
          coef  = one/sin(angle)

          if (mod(ntot,2_intType) == 0) &
               coef = coef*cos(angle)

          coef = coef*slicesFact

          do j=1,3
             do i=1,3
                diagMatCoefSpectral(mm,i,j) = &
                     diagMatCoefSpectral(mm,i,j) - coef*rotMat(i,j)
             enddo
          enddo

       enddo slicesLoop

       ! The matrix coefficients must be multiplied by the leading
       ! coefficient, which depends on the actual periodic time.

       coef = pi*timeRef/sections(mm)%timePeriod

       do j=1,3
          do i=1,3
             diagMatCoefSpectral(mm,i,j) = &
                  coef*diagMatCoefSpectral(mm,i,j)
          enddo
       enddo

       do nn=1,(nTimeIntervalsSpectral-1)

          do j=1,3
             do i=1,3
                matrixCoefSpectral(mm,nn,i,j) = &
                     coef*matrixCoefSpectral(mm,nn,i,j)
             enddo
          enddo

       end do

    enddo sectionLoop

  end subroutine timeSpectralCoef

  subroutine timeSpectralMatrices
    !
    !       timeSpectralMatrices computes the matrices for the time
    !       derivative of the time spectral method for all sections. For
    !       scalar quantities these matrices only differ if sections have
    !       different periodic times. For vector quantities, such as
    !       momentum, these matrices can be different depending on whether
    !       the section is rotating or not and the number of slices
    !       present.
    !
    use constants
    use inputPhysics, only : equationMode
    use inputTimeSpectral, only : nTimeIntervalsSpectral, dscalar, dvector, &
         rotMatrixSpectral
    use section, only: nSections
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, mm, ll, kk, ii
    integer(kind=intType) :: i, j

    real(kind=realType), dimension(3,3) :: tmpMat

    real(kind=realType), dimension(:,:), allocatable :: coefSpectral
    real(kind=realType), dimension(:,:,:,:), allocatable :: &
         matrixCoefSpectral
    real(kind=realType), dimension(:,:,:), allocatable :: &
         diagMatCoefSpectral
    !
    ! This routine is only used for the spectral solutions. Return
    ! immediately if a different mode is solved.

    if(equationMode /= timeSpectral) return

    ! Allocate the memory for the matrices as well as the help
    ! variables needed to construct these matrices.

    !added to allow second call in mdUpdateRoutines
    if( allocated(dscalar))deallocate(dscalar)
    if( allocated(dvector))deallocate(dvector)

    nn = nTimeIntervalsSpectral
    mm = 3*nn
    kk = nn - 1



    allocate(dscalar(nSections,nn,nn),             &
         dvector(nSections,mm,mm),             &
         coefSpectral(nSections,kk),           &
         matrixCoefSpectral(nSections,kk,3,3), &
         diagMatCoefSpectral(nSections,3,3), stat=ierr)
    if(ierr /= 0)                              &
         call terminate("timeSpectralMatrices", &
         "Memory allocation failure for the matrices of &
         &the spectral time derivatives.")

    ! Determine the help variables needed to construct the
    ! actual matrices.

    call timeSpectralCoef(coefSpectral, matrixCoefSpectral, &
         diagMatCoefSpectral)
    !
    !       Determine the time derivative matrices for the sections.
    !
    ! Loop over the number of sections.

    sectionLoop: do ii=1,nSections
       !
       !         Matrix for scalar quantities.
       !
       ! Loop over the number of rows.

       do nn=1,nTimeIntervalsSpectral

          ! Set the diagonal element to zero, i.e. there is no
          ! contribution to the own time derivative.

          dscalar(ii,nn,nn) = zero

          ! Loop over the rest of the columns.

          do mm=1,(nTimeIntervalsSpectral - 1)

             ! Determine the corresponding column index.

             ll = nn + mm
             if(ll > nTimeIntervalsSpectral) &
                  ll = ll - nTimeIntervalsSpectral

             ! Store the corresponding coefficient in dscalar.

             dscalar(ii,nn,ll) = coefSpectral(ii,mm)

          enddo
       enddo
       !
       !         Matrices for vector quantities.
       !
       ! Loop over the number of time intervals; the number of rows
       ! is 3 times this number.

       rowLoop: do nn=1,nTimeIntervalsSpectral

          ! Initialize the diagonal block to diagMatCoefSpectral,
          ! the additional diagonal entry needed for the rotational
          ! periodicity.

          kk = 3*(nn-1)
          do j=1,3
             do i=1,3
                dvector(ii,kk+i,kk+j) = diagMatCoefSpectral(ii,i,j)
             enddo
          enddo

          ! Loop over the other time intervals, which contribute to
          ! the time derivative.

          columnLoop: do mm=1,(nTimeIntervalsSpectral - 1)

             ! Determine the corresponding column index and check the
             ! situation we are having here.

             ll = nn + mm
             if(ll > nTimeIntervalsSpectral) then

                ! Index is outside the range and a shift must be applied.

                ll = ll - nTimeIntervalsSpectral

                ! The vector must be rotated. This effect is incorporated
                ! directly in the matrix of time derivatives.

                do j=1,3
                   do i=1,3
                      tmpMat(i,j) = matrixCoefSpectral(ii,mm,i,1) &
                           * rotMatrixSpectral(ii,1,j)     &
                           + matrixCoefSpectral(ii,mm,i,2) &
                           * rotMatrixSpectral(ii,2,j)     &
                           + matrixCoefSpectral(ii,mm,i,3) &
                           * rotMatrixSpectral(ii,3,j)
                   enddo
                enddo

             else

                ! Index is in the range. Copy the matrix coefficient
                ! into tmpMat.

                do j=1,3
                   do i=1,3
                      tmpMat(i,j) = matrixCoefSpectral(ii,mm,i,j)
                   enddo
                enddo

             endif

             ! Determine the offset for the column index and store
             ! this submatrix in the correct place of dvector.

             ll = 3*(ll-1)
             do j=1,3
                do i=1,3
                   dvector(ii,kk+i,ll+j) = tmpMat(i,j)
                enddo
             enddo

          enddo columnLoop
       enddo rowLoop
    enddo sectionLoop

    ! Release the memory of the help variables needed to construct
    ! the matrices of the time derivatives.

    deallocate(coefSpectral, matrixCoefSpectral, &
         diagMatCoefSpectral, stat=ierr)
    if(ierr /= 0)                            &
         call terminate("timeSpectralMatrices", &
         "Deallocation failure for the help variables.")

  end subroutine timeSpectralMatrices


  subroutine readRestartFile()
    !
    !       readRestartFile reads the fine grid solution(s) from the
    !       restart file(s). If the restart file(s) do not correspond to
    !       the current mesh, the solution(s) are interpolated onto this
    !       mesh. It is also allowed to change boundary conditions, e.g.
    !       an alpha and/or Mach sweep is possible. Furthermore there is
    !       some support when starting from a different turbulence model,
    !       although this should be used with care.
    !
    use constants
    use cgnsGrid
    use su_cgns
    use variableReading ! Full import since we need basically everything
    use blockPointers, only: iBegOr, jBegOr, kBegOr, il, jl, kl, nDom, &
         nBKGlobal, nx, ny, nz
    use communication, only : adflow_comm_world, myid
    use inputPhysics, only : equationMode
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use monitor, only : nTimeStepsRestart, timeUnsteadyRestart
    use utils, only : terminate, setPointers
    use sorting, only : bsearchStrings
    implicit none
    !
    !      Local variables.
    !
    integer :: nZones, cellDim, physDim, ierr, nSols

    integer(cgsize_t), dimension(9) :: sizes
    integer, dimension(9) :: rindSizes
    integer, dimension(nSolsRead) :: fileIDs

    integer(kind=intType) :: ii, jj, nn
    integer(kind=intType) :: nTypeMismatch
    integer(kind=intType) :: nHiMin, nHjMin, nHkMin
    integer(kind=intType) :: nHiMax, nHjMax, nHkMax

    character(len=7)              :: integerString
    character(len=maxCGNSNameLen) :: cgnsName
    character(len=2*maxStringLen) :: errorMessage

    ! Initialize halosRead to .true. This will be overwritten if
    ! there is at least one block present for which the halo data
    ! cannot be read.

    halosRead = .true.

    ! Initialize nTypeMismatch to 0.

    nTypeMismatch = 0

    ! Loop over the number of files to be read and open them.

    fileOpenLoop: do solID=1,nSolsRead

       ! Open the restart file for reading.

       call cg_open_f(solFiles(solID), mode_read, cgnsInd, ierr)
       if(ierr /= all_ok) then
          write(errorMessage,*) "File ", trim(solFiles(solID)), &
               " could not be opened for reading"
          call terminate("readRestartFile", errorMessage)
       endif

       fileIDs(solID) = cgnsInd

       ! Determine the number of bases in the cgns file.
       ! This must be at least 1.

       call cg_nbases_f(cgnsInd, cgnsBase, ierr)
       if(ierr /= all_ok)                  &
            call terminate("readRestartFile", &
            "Something wrong when calling cg_nbases_f")

       if(CGNSBase < 1) then
          write(errorMessage,*) "CGNS file ", trim(solFiles(solID)), &
               " does not contain a base"
          call terminate("readRestartFile", errorMessage)
       endif

       ! Only data from the first base is read. Information from
       ! higher bases is ignored.

       cgnsBase = 1

       ! Read the cell and physical dimensions as well as the name for
       ! this base.

       call cg_base_read_f(cgnsInd, cgnsBase, cgnsName, cellDim, &
            physDim, ierr)
       if(ierr /= all_ok)                  &
            call terminate("readRestartFile", &
            "Something wrong when calling cg_base_read_f")

       ! Check the cell and physical dimensions. Both must be 3 for
       ! this code to work.

       if(cellDim /= 3 .or. physDim /= 3) then
          write(errorMessage,100) cellDim, physDim
100       format("Both the number of cell and physical dimensions &
               &should be 3, not",1X,I1,1X,"and",1X,I1)
          call terminate("readRestartFile", errorMessage)
       endif

    enddo fileOpenLoop


    ! Broadcast nTimeStepsRestart and timeUnsteadyRestart to all
    ! processors. These values are needed to perform a consistent
    ! unsteady restart.

    call mpi_bcast(nTimeStepsRestart, 1, adflow_integer, 0, &
         ADflow_comm_world, ierr)
    call mpi_bcast(timeUnsteadyRestart, 1, adflow_real, 0, &
         ADflow_comm_world, ierr)

    ! Get the scaling factors for density, pressure and velocity
    ! by reading the reference state.

    call scaleFactors(fileIDs)

    ! Loop over the number of files to be read and read the solution.

    solLoop: do solID=1,nSolsRead

       ! Store the file index a bit easier and set the base to 1.

       cgnsInd  = fileIDs(solID)
       cgnsBase = 1

       ! Determine the number of zones (blocks) in the restart file
       ! and check if this is identical to the number in the grid file.

       call cg_nzones_f(cgnsInd, cgnsBase, nZones, ierr)
       if(ierr /= all_ok)                  &
            call terminate("readRestartFile", &
            "Something wrong when calling cg_nzones_f")

       if(nZones /= cgnsNdom)              &
            call terminate("readRestartFile", &
            "Number of blocks in grid file and restart &
            &file differ")

       ! Create a sorted version of the zone names of the restart file
       ! and store its corresponding zone numbers in zoneNumbers.

       call getSortedZoneNumbers

       ! Loop over the number of blocks stored on this processor.

       domains: do nn=1,nDom

          ! Set the pointers for this block. Make sure that the
          ! correct data is set.

          ii = min(solID,nTimeIntervalsSpectral)
          call setPointers(nn, 1_intType, ii)

          ! Store the zone name of the original grid a bit easier.

          cgnsName = cgnsDoms(nbkGlobal)%zoneName

          ! Search in the sorted zone names of the restart file for
          ! cgnsName. The name must be found; otherwise the restart
          ! is pointless. If found, the zone number is set accordingly.

          jj = bsearchStrings(cgnsname, zoneNames)
          if(jj == 0) then
             write(errorMessage,*) "Zone name ", trim(cgnsName),  &
                  " not found in restart file ", &
                  trim(solFiles(solID))
             call terminate("readRestartFile", errorMessage)
          else
             jj = zoneNumbers(jj)
          endif

          cgnsZone = jj

          ! Determine the dimensions of the zone and check if these are
          ! identical to the dimensions of the block in the grid file.

          call cg_zone_read_f(cgnsInd, cgnsBase, cgnsZone, &
               cgnsname, sizes, ierr)
          if(ierr /= all_ok)                  &
               call terminate("readRestartFile", &
               "Something wrong when calling &
               &cg_zone_read_f")

          if(cgnsDoms(nbkGlobal)%il /= sizes(1) .or. &
               cgnsDoms(nbkGlobal)%jl /= sizes(2) .or. &
               cgnsDoms(nbkGlobal)%kl /= sizes(3))     &
               call terminate("readRestartFile", &
               "Corresponding zones in restart file and &
               &grid file have different dimensions")

          ! Determine the number of flow solutions in this zone and
          ! check if there is a solution stored.

          call cg_nsols_f(cgnsInd, cgnsBase, cgnsZone, nSols, ierr)
          if(ierr /= all_ok)                  &
               call terminate("readRestartFile", &
               "Something wrong when calling cg_nsols_f")

          if(nSols == 0)                      &
               call terminate("readRestartFile", &
               "No solution present in restart file")

          ! Check for multiple solutions. A distinction is needed for
          ! overset cases because there will be an extra solution node
          ! for the nodal iblanks.

          if((nSols > 1) .or. nSols > 2) &
               call terminate("readRestartFile", &
               "Multiple solutions present in restart file")

          ! Determine the location of the solution variables. A loop is
          ! done over the solution nodes which is either 1 or 2. In the
          ! latter case, pick the node not named "Nodal Blanks".

          do cgnsSol=1,nSols
             call cg_sol_info_f(cgnsInd, cgnsBase, cgnsZone, cgnsSol, &
                  cgnsName, location, ierr)
             if(ierr /= all_ok)                  &
                  call terminate("readRestartFile", &
                  "Something wrong when calling &
                  &cg_sol_info_f")

             if (trim(cgnsName) /= "Nodal Blanks") exit
          end do

          ! Determine the rind info.

          call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", &
               cgnsZone, "FlowSolution_t", cgnsSol, "end")
          if(ierr /= all_ok)                  &
               call terminate("readRestartFile", &
               "Something wrong when calling cg_goto_f")

          call cg_rind_read_f(rindSizes, ierr)
          if(ierr /= all_ok)                  &
               call terminate("readRestartFile", &
               "Something wrong when calling &
               &cg_rind_read_f")

          ! Check if halo's are present. If not, set halosRead to .false.
          ! This only needs to be done if this is not an older state for
          ! an unsteady computation.

          if(solID == 1 .or. equationMode == timeSpectral) then
             if(rindSizes(1) == 0 .or. rindSizes(2) == 0 .or. rindSizes(3) == 0 .or. &
                  rindSizes(4) == 0 .or. rindSizes(5) == 0 .or. rindSizes(6) == 0)     &
                  halosRead = .false.
          endif

          ! Initialize the number of halo cells to read to 0.

          nHiMin = 0; nHjMin = 0; nHkMin = 0
          nHiMax = 0; nHjMax = 0; nHkMax = 0

          ! Determine the range which must be read. A few things must be
          ! taken into account: - in iBegor, iEndor, etc. The nodal
          !                       range is stored. As in CGNS the cell
          !                       range start at 1, 1 must be subtracted
          !                       from the upper bound.
          !                     - the rind info must be taken into
          !                       account, because only the upper bound
          !                       is changed in cgns; the lower bound
          !                       remains 1.
          !                     - in case the solution is stored in the
          !                       vertices one extra variable in each
          !                       direction is read. An averaging will
          !                       take place to obtain cell centered
          !                       values.
          ! Also when vertex data is present, set halosRead to .false.,
          ! because it is not possible to determine the halos.

          if(location == CellCenter) then

             ! Correct the number of halo cells to be read.
             ! Only if this is not an older state in time for an
             ! unsteady computation.

             if(solID == 1 .or. equationMode == timeSpectral) then
                if(rindSizes(1) > 0) nHiMin = 1; if(rindSizes(2) > 0) nHiMax = 1
                if(rindSizes(3) > 0) nHjMin = 1; if(rindSizes(4) > 0) nHjMax = 1
                if(rindSizes(5) > 0) nHkMin = 1; if(rindSizes(6) > 0) nHkMax = 1
             endif

             ! Set the cell range to be read from the CGNS file.

             rangeMin(1) = iBegOr + rindSizes(1) - nHiMin
             rangeMin(2) = jBegOr + rindSizes(3) - nHjMin
             rangeMin(3) = kBegOr + rindSizes(5) - nHkMin

             rangeMax(1) = rangeMin(1) + nx-1 + nHiMin + nHiMax
             rangeMax(2) = rangeMin(2) + ny-1 + nHjMin + nHjMax
             rangeMax(3) = rangeMin(3) + nz-1 + nHkMin + nHkMax

          else if(location == Vertex) then

             ! Set the nodal range such that enough info is present
             ! to average the nodal data to create the cell centered
             ! data in the owned cells. No halo cells will be
             ! initialized.

             halosRead   = .false.

             rangeMin(1) = iBegor + rindSizes(1)
             rangeMin(2) = jBegor + rindSizes(3)
             rangeMin(3) = kBegor + rindSizes(5)

             rangeMax(1) = rangeMin(1) + nx
             rangeMax(2) = rangeMin(2) + ny
             rangeMax(3) = rangeMin(3) + nz
          else
             call terminate("readRestartFile", &
                  "Only CellCenter or Vertex data allowed in &
                  &restart file")
          endif

          ! Allocate the memory for buffer, needed to store the variable
          ! to be read, and bufferVertex in case the solution is stored
          ! in the vertices.

          allocate(buffer(2-nHiMin:il+nHiMax, &
               2-nHjMin:jl+nHjMax, &
               2-nHkMin:kl+nHkMax), stat=ierr)
          if(ierr /= 0)                       &
               call terminate("readRestartFile", &
               "Memory allocation failure for buffer")

          if(location == Vertex) then
             allocate(bufferVertex(1:il,1:jl,1:kl), stat=ierr)
             if(ierr /= 0)                       &
                  call terminate("readRestartFile", &
                  "Memory allocation failure for bufferVertex")
          endif

          ! Create a sorted version of the variable names and store the
          ! corresponding type in varTypes.

          call getSortedVarNumbers

          ! Read the density and the turbulence variables.

          call readDensity(nTypeMismatch)
          call readTurbvar(nTypeMismatch)

          ! Read the other variables, depending on the situation.

          testPrim: if(solID == 1 .or. equationMode == timeSpectral) then

             ! Either the first solution or time spectral mode. Read
             ! the primitive variables from the restart file.

             call readXvelocity(nTypeMismatch)
             call readYvelocity(nTypeMismatch)
             call readZvelocity(nTypeMismatch)
             call readPressure(nTypeMismatch)

          else testPrim

             ! Old solution in unsteady mode. Read the conservative
             ! variables.

             call readXmomentum(nTypeMismatch)
             call readYmomentum(nTypeMismatch)
             call readZmomentum(nTypeMismatch)
             call readEnergy(nTypeMismatch)

          endif testPrim

          ! Release the memory of buffer, varNames and varTypes.

          deallocate(buffer, varNames, varTypes, stat=ierr)
          if(ierr /= 0)                       &
               call terminate("readRestartFile", &
               "Deallocation error for buffer, varNames &
               &and varTypes.")

          ! In case bufferVertex is allocated, release it.

          if(location == Vertex) then
             deallocate(bufferVertex, stat=ierr)
             if(ierr /= 0)                       &
                  call terminate("readRestartFile", &
                  "Deallocation error for bufferVertex")
          endif

       enddo domains

       ! Release the memory of zoneNames and zoneNumbers.

       deallocate(zoneNames, zoneNumbers, stat=ierr)
       if(ierr /= 0)                       &
            call terminate("readRestartFile", &
            "Deallocation failure for zoneNames &
            &and zoneNumbers.")

       ! Close the cgns solution file.

       call cg_close_f(cgnsInd, ierr)
       if(ierr /= all_ok)                  &
            call terminate("readRestartFile", &
            "Something wrong when calling cg_close_f")

    enddo solLoop


    ! Determine the global sum of nTypeMismatch; the result only
    ! needs to be known on processor 0. Use ii as the global buffer
    ! to store the result. If a type mismatch occured,
    ! print a warning.

    call mpi_reduce(nTypeMismatch, ii, 1, adflow_integer, &
         mpi_sum, 0, ADflow_comm_world, ierr)
    if(myID == 0 .and. ii > 0) then

       write(integerString,"(i6)") ii
       integerString = adjustl(integerString)

       print "(a)", "#"
       print "(a)", "#                      Warning"
       print 120, trim(integerString)
       print "(a)", "#"
120    format("# ",a," type mismatches occured when reading the &
            &solution of the blocks")
    endif

  end subroutine readRestartFile

  subroutine getSortedZoneNumbers
    !
    !       getSortedZoneNumbers reads the names of the zones of the
    !       cgns file given by cgnsInd and cgnsBase. Afterwards the
    !       zonenames are sorted in increasing order, such that a binary
    !       search algorithm can be employed. The original zone numbers
    !       are stored in zoneNumbers.
    !       If the zone contains a link to a zone containing the
    !       coordinates the name of the linked zone is taken.
    !

    use constants
    use cgnsGrid, only : cgnsNDom
    use su_cgns
    use variableReading, only : zoneNames, zoneNumbers, cgnsInd, cgnsBase
    use sorting, only : qsortStrings, bsearchStrings
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr
    integer :: zone, zonetype, ncoords, pathLength
    integer :: pos

    integer(kind=cgsize_t), dimension(9) :: sizesBlock

    integer(kind=intType) :: nn, ii

    character(len=maxStringLen) :: errorMessage, linkPath
    character(len=maxCGNSNameLen), dimension(cgnsNdom) :: tmpNames

    logical :: nameFound

    character(len=7) :: int1String, int2String

    ! Allocate the memory for zoneNames and zoneNumbers.

    allocate(zoneNames(cgnsNdom), zoneNumbers(cgnsNdom), stat=ierr)
    if(ierr /= 0)                            &
         call terminate("getSortedZoneNumbers", &
         "Memory allocation failure for zoneNames &
         &and zoneNumbers")

    ! Loop over the number of zones in the file.

    cgnsDomains: do nn=1,cgnsNdom

       ! Initialize nameFound to .false.

       nameFound = .false.

       ! Check if the zone is structured.

       zone = nn
       call cg_zone_type_f(cgnsInd, cgnsBase, zone, zonetype, ierr)
       if(ierr /= all_ok)                       &
            call terminate("getSortedZoneNumbers", &
            "Something wrong when calling cg_zone_type_f")

       if(zonetype /= structured) then

          write(int1String,"(i7)") cgnsBase
          int1String = adjustl(int1String)
          write(int2String,"(i7)") zone
          int2String = adjustl(int2String)

          write(errorMessage,100) trim(int1String), trim(int2String)
100       format("Base",1X,A,": Zone",1X,A, " of the cgns restart &
               &file is not structured")
          call terminate("getSortedZoneNumbers", errorMessage)

       endif

       ! Determine the number of grid coordinates of this zone.

       call cg_ncoords_f(cgnsInd, cgnsBase, zone, ncoords, ierr)
       if(ierr /= all_ok)                       &
            call terminate("getSortedZoneNumbers", &
            "Something wrong when calling cg_ncoords_f")

       ! If ncoords == 3, there are coordinates present. Then check if
       ! it is a link. If so, take the zone name of the link.

       if(ncoords == 3) then

          ! Go to the coordinates node.

          call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", zone, &
               "GridCoordinates_t", 1, "end")
          if(ierr /= all_ok)                       &
               call terminate("getSortedZoneNumbers", &
               "Something wrong when calling cg_goto_f")

          ! Check if this node is a link.

          call cg_is_link_f(pathLength, ierr)
          if(ierr /= all_ok)                       &
               call terminate("getSortedZoneNumbers", &
               "Something wrong when calling cg_is_link_f")

          if(pathLength > 0) then

             ! Determine the name of the linkPath.

             call cg_link_read_f(errorMessage, linkPath, ierr)
             if(ierr /= all_ok)                       &
                  call terminate("getSortedZoneNumbers", &
                  "Something wrong when calling &
                  &cg_link_read_f")

             ! Find the zone name.
             ! Find, starting from the back, the forward slash.

             pos = index(linkPath, "/", .true.)
             if(pos > 0) then
                linkPath = linkPath(:pos-1)

                ! Find the next forward slash from the back and
                ! remove the leading part from the path name.

                pos = index(linkPath, "/", .true.)
                if(pos > 0) linkPath = linkPath(pos+1:)
             endif

             ! Create the zone name and set nameFound to .true..

             linkPath      = adjustl(linkPath)
             zoneNames(nn) = trim(linkPath)
             nameFound     = .true.

          endif
       endif

       ! If no name was yet found set it to the name of the
       ! current zone.

       if(.not. nameFound) then
          call cg_zone_read_f(cgnsInd, cgnsBase, zone, &
               zoneNames(nn), sizesBlock, ierr)
          if(ierr /= all_ok)                       &
               call terminate("getSortedZoneNumbers", &
               "Something wrong when calling &
               &cg_zone_read_f")
       endif

    enddo cgnsDomains

    ! Set tmpNames to zoneNames and sort the latter
    ! in increasing order.

    do nn=1,cgnsNdom
       tmpNames(nn) = zoneNames(nn)
    enddo

    ! Sort zoneNames in increasing order.

    call qsortStrings(zoneNames, cgnsNdom)

    ! Initialize zoneNumbers to -1. This serves as a check during
    ! the search.

    do nn=1,cgnsNdom
       zoneNumbers(nn) = -1
    enddo

    ! Find the original zone numbers for the sorted zone names.

    do nn=1,cgnsNdom
       ii = bsearchStrings(tmpNames(nn), zoneNames)

       ! Check if the zone number is not already taken. If this is the
       ! case, this means that two identical zone names are present.

       if(zoneNumbers(ii) /= -1)                &
            call terminate("getSortedZoneNumbers", &
            "Error occurs only when two identical zone &
            &names are present")

       ! And set the zone number.

       zoneNumbers(ii) = nn
    enddo

  end subroutine getSortedZoneNumbers

  subroutine getSortedVarNumbers
    !
    !       getSortedVarNumbers reads the names of variables stored in
    !       the given solution node of the cgns file, indicated by
    !       cgnsInd, cgnsBase and cgnsZone. Afterwards the variable
    !       names are sorted in increasing order, such that they can be
    !       used in a binary search. Their original variable number and
    !       type is stored.
    !
    use constants
    use su_cgns
    use variableReading, only : varTypes, varNames, cgnsBase, cgnsInd, &
         cgnsZone, cgnsSol, nVar
    use sorting, only : qsortStrings, bsearchStrings
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: i, ierr
    integer, dimension(:), allocatable :: tmpTypes

    integer(kind=intType) :: nn, ii

    integer(kind=intType), dimension(:), allocatable :: varNumbers

    character(len=maxCGNSNameLen), allocatable, dimension(:) :: &
         tmpNames

    ! Determine the number of solution variables stored.

    call cg_nfields_f(cgnsInd, cgnsBase, cgnsZone, cgnsSol, &
         nVar, ierr)
    if(ierr /= all_ok)                      &
         call terminate("getSortedVarNumbers", &
         "Something wrong when calling cg_nfield_f")

    ! Allocate the memory for varnames, vartypes and varnumber

    allocate(varNames(nVar), varTypes(nVar), varNumbers(nVar), &
         stat=ierr)
    if(ierr /= 0)                           &
         call terminate("getSortedVarNumbers", &
         "Memory allocation failure for varNames, etc.")

    ! Loop over the number of variables and store their names and
    ! types.

    do i=1,nVar
       call cg_field_info_f(cgnsInd, cgnsBase, cgnsZone, cgnsSol, &
            i, varTypes(i), varNames(i), ierr)
       if(ierr /= 0)                           &
            call terminate("getSortedVarNumbers", &
            "Something wrong when calling cg_field_info_f")
    enddo

    ! Allocate the memory for tmpTypes and tmpNames and initialize
    ! their values.

    allocate(tmpTypes(nVar), tmpNames(nVar), stat=ierr)
    if(ierr /= 0)                           &
         call terminate("getSortedVarNumbers", &
         "Memory allocation failure for tmp variables")

    do i=1,nVar
       tmpTypes(i) = varTypes(i)
       tmpNames(i) = varNames(i)
    enddo

    ! Sort varNames in increasing order.

    nn = nVar
    call qsortStrings(varNames, nn)

    ! Initialize varNumbers to -1. This serves as a check during
    ! the search.

    do i=1,nVar
       varNumbers(i) = -1
    enddo

    ! Find the original types and numbers for the just sorted
    ! variable names.

    do i=1,nVar
       ii = bsearchStrings(tmpNames(i), varNames)

       ! Check if the variable number is not already taken. If this is
       ! the case, this means that two identical var names are present.

       if(varNumbers(ii) /= -1)                &
            call terminate("getSortedVarNumbers", &
            "Error occurs only when two identical &
            &variable names are present")

       ! And set the variable number and type.

       varNumbers(ii) = i
       varTypes(ii)   = tmpTypes(i)
    enddo

    ! Release the memory of varNumbers, tmpNames and tmpTypes.

    deallocate(varNumbers, tmpTypes, tmpNames, stat=ierr)
    if(ierr /= 0)                           &
         call terminate("getSortedVarNumbers", &
         "Deallocation error for tmp variables")

  end subroutine getSortedVarNumbers
#endif
end module initializeFlow
