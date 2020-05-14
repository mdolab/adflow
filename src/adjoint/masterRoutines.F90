module masterRoutines
contains
  subroutine master(useSpatial, famLists, funcValues, forces, &
       bcDataNames, bcDataValues, bcDataFamLists)

    use constants
    use communication, only : adflow_comm_world
    use BCRoutines, only : applyallBC_block
    use bcdata, only : setBCData, setBCDataFineGrid
    use turbbcRoutines, only : applyallTurbBCthisblock, bcTurbTreatment
    use iteration, only : currentLevel
    use inputAdjoint,  only : viscPC
    use flowVarRefState, only : nwf, nw
    use blockPointers, only : nDom, il, jl, kl
    use flowVarRefState, only : viscous
    use inputPhysics , only : turbProd, equationMode, equations, turbModel
    use inputDiscretization, only : lowSpeedPreconditioner, lumpedDiss, spaceDiscr, useAPproxWallDistance
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use initializeFlow, only : referenceState
    use section, only: sections, nSections
    use monitor, only : timeUnsteadyRestart
    use sa, only : saSource, saViscous, saResScale, qq
    use haloExchange, only : exchangeCoor, whalo2
    use wallDistance, only : updateWallDistancesQuickly
    use solverUtils, only : timeStep_block
    use flowUtils, only : allNodalGradients, computeLamViscosity, computePressureSimple, &
         computeSpeedOfSoundSquared, adjustInflowAngle
    use fluxes, only : inviscidDissFluxScalarApprox, inviscidDissFluxMatrixApprox, &
         inviscidUpwindFlux, inviscidDissFluxScalar, inviscidDissFluxMatrix, &
         viscousFlux, viscousFluxApprox, inviscidCentralFlux
    use utils, only : setPointers, EChk
    use turbUtils, only : turbAdvection, computeEddyViscosity
    use residuals, only : initRes_block, sourceTerms_block
    use surfaceIntegrations, only : getSolution
    use adjointExtra, only : volume_block, metric_block, boundaryNormals,&
         xhalo_block, sumdwandfw, resScale
    use oversetData, only : oversetPresent
    use inputOverset, only : oversetUpdateMode
    use oversetCommUtilities, only : updateOversetConnectivity
    use actuatorRegionData, only : nActuatorRegions
    implicit none

    ! Input Arguments:
    logical, intent(in) :: useSpatial
    integer(kind=intType), optional, dimension(:, :), intent(in) :: famLists
    real(kind=realType), optional, dimension(:, :), intent(out) :: funcValues
    character, optional, dimension(:, :), intent(in) :: bcDataNames
    real(kind=realType), optional, dimension(:), intent(in) :: bcDataValues
    integer(kind=intType), optional, dimension(:, :) :: bcDataFamLists

    ! Output Variables
    real(kind=realType), intent(out), optional, dimension(:, :, :) :: forces

    ! Working Variables
    integer(kind=intType) :: ierr, nn, sps, fSize, iRegion
    real(kind=realType), dimension(nSections) :: t
    real(kind=realType) :: dummyReal

    if (useSpatial) then
       call adjustInflowAngle()

       ! Update all the BCData
       call referenceState
       if (present(bcDataNames)) then
          do sps=1,nTimeIntervalsSpectral
             call setBCData(bcDataNames, bcDataValues, bcDataFamLists, sps, &
                  size(bcDataValues), size(bcDataFamLIsts, 2))
          end do
          call setBCDataFineGrid(.true.)
       end if

        do sps=1,nTimeIntervalsSpectral
          do nn=1,nDom
             call setPointers(nn, 1, sps)
             call xhalo_block()
          end do
       end do

       ! Now exchange the coordinates (fine level only)
       call exchangecoor(1)

       do sps=1, nTimeIntervalsSpectral
          ! Update overset connectivity if necessary
          if (oversetPresent .and. &
               (oversetUpdateMode == updateFast .or. &
               oversetUpdateMode == updateFull)) then
             call updateOversetConnectivity(1_intType, sps)
          end if
       end do
    end if

    do sps=1,nTimeIntervalsSpectral
       do nn=1,nDom
          call setPointers(nn, 1, sps)

          if (useSpatial) then
             call volume_block
             call metric_block
             call boundaryNormals

             if (equations == RANSEquations .and. useApproxWallDistance) then
                call updateWallDistancesQuickly(nn, 1, sps)
             end if
          end if

          ! Compute the pressures/viscositites
          call computePressureSimple(.False.)

          ! Compute Laminar/eddy viscosity if required
          call computeLamViscosity(.False.)
          call computeEddyViscosity(.False.)

          ! Make sure to call the turb BC's first incase we need to
          ! correct for K
          if (equations == RANSequations) then
             call BCTurbTreatment
             call applyAllTurbBCthisblock(.True.)
          end if
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

    ! Main loop for the residual
    do sps=1, nTimeIntervalsSpectral
       do nn=1, nDom
          call setPointers(nn, 1, sps)
          call initRes_block(1, nw, nn, sps)
          do iRegion=1, nActuatorRegions
             call sourceTerms_block(nn, .True., iRegion, dummyReal)
          end do

          ! Compute turbulence residual for RANS equations
          if( equations == RANSEquations) then

             ! Initialize only the Turblent Variables
             !call unsteadyTurbSpectral_block(itu1, itu1, nn, sps)

             select case (turbModel)

             case (spalartAllmaras)
                allocate(qq(2:il, 2:jl, 2:kl))
                call saSource
                call turbAdvection(1_intType, 1_intType, itu1-1, qq)
                !call unsteadyTurbTerm(1_intType, 1_intType, itu1-1, qq)
                call saViscous
                call saResScale
                deallocate(qq)
             end select
          endif

          ! Compute the mean flow residuals
          call timeStep_block(.false.)
          call inviscidCentralFlux
          if (lumpedDiss) then
             select case (spaceDiscr)
             case (dissScalar)
                call inviscidDissFluxScalarApprox
             case (dissMatrix)
                call inviscidDissFluxMatrixApprox
             case (upwind)
                call inviscidUpwindFlux(.True.)
             end select
          else
             select case (spaceDiscr)
             case (dissScalar)
                call inviscidDissFluxScalar
             case (dissMatrix)
                call inviscidDissFluxMatrix
             case (upwind)
                call inviscidUpwindFlux(.True.)
             end select
          end if

          if (viscous) then
             call computeSpeedOfSoundSquared
             if (.not. lumpedDiss .or. viscPC) then
                call allNodalGradients
                call viscousFlux
             else
                call viscousFluxApprox
             end if
          end if
          call sumDwAndFw
          ! if (lowSpeedPreconditioner) then
          !    call applyLowSpeedPreconditioner
          ! end if
          call resScale
       end do
    end do

    ! Compute the final solution values
    if (present(famLists)) then
       call getSolution(famLists, funcValues)
    end if

    do sps=1, nTimeIntervalsSpectral
       if (present(forces)) then
          ! Now we can retrieve the forces/tractions for this spectral instance
          fSize = size(forces, 2)
          call getForces(forces(:, :, sps), fSize, sps)
       end if
    end do

  end subroutine master
#ifndef USE_COMPLEX
  subroutine master_d(wdot, xdot, forcesDot, dwDot, famLists, funcValues, funcValuesd, &
       bcDataNames, bcDataValues, bcDataValuesd, bcDataFamLists)

    use constants
    use diffsizes, only :  ISIZE1OFDrfbcdata, ISIZE1OFDrfviscsubface
    use communication, only : adflow_comm_world
    use iteration, only : currentLevel
    use BCExtra_d, only : applyAllBC_Block_d
    use inputAdjoint,  only : viscPC
    use flowVarRefState, only : nw, nwf
    use blockPointers, only : nDom, il, jl, kl, wd, xd, dw, dwd, nBocos, nViscBocos
    use flowVarRefState, only : viscous, timerefd
    use inputPhysics , only : turbProd, equationMode, equations, turbModel, wallDistanceNeeded
    use inputDiscretization, only : lowSpeedPreconditioner, lumpedDiss, spaceDiscr, useAPproxWallDistance
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use section, only: sections, nSections
    use monitor, only : timeUnsteadyRestart
    use utils, only : isWallType, setPointers, setPointers_d, EChk
    use sa_d, only : saSource_d, saViscous_d, saResScale_d, qq
    use turbutils_d, only : turbAdvection_d, computeEddyViscosity_d
    use fluxes_d, only :inviscidDissFluxScalarApprox_d, inviscidDissFluxMatrixApprox_d, &
         inviscidUpwindFlux_d, inviscidDissFluxScalar_d, inviscidDissFluxMatrix_d, &
         inviscidUpwindFlux_d, viscousFlux_d, viscousFluxApprox_d, inviscidCentralFlux_d
    use residuals_d, only : sourceTerms_block_d
    use adjointPETSc, only : x_like
    use haloExchange, only : whalo2_d, exchangeCoor_d, exchangeCoor, whalo2
    use wallDistance_d, only : updateWallDistancesQuickly_d
    use wallDistanceData, only : xSurfVec, xSurfVecd, xSurf, xSurfd, wallScatter
    use flowutils_d, only : computePressureSimple_d, computeLamViscosity_d, &
         computeSpeedOfSoundSquared_d, allNodalGradients_d, adjustInflowAngle_d
    use solverutils_d, only : timeStep_Block_d
    use turbbcroutines_d, only : applyAllTurbBCthisblock_d,  bcTurbTreatment_d
    use initializeflow_d, only : referenceState_d
    use surfaceIntegrations, only : getSolution_d
    use adjointExtra_d, only : xhalo_block_d, volume_block_d, metric_BLock_d, boundarynormals_d
    use adjointextra_d, only : resscale_D, sumdwandfw_d
    use bcdata, only : setBCData_d, setBCDataFineGrid_d
    use oversetData, only : oversetPresent
    use inputOverset, only : oversetUpdateMode
    use oversetCommUtilities, only : updateOversetConnectivity_d
    use actuatorRegionData, only : nActuatorRegions
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input Arguments:
    real(kind=realType), intent(in), dimension(:) :: wDot, xDot
    integer(kind=intType), optional, dimension(:, :), intent(in) :: famLists
    real(kind=realType), optional, dimension(:, :), intent(out) :: funcValues, funcValuesd
    character, optional, dimension(:, :), intent(in) :: bcDataNames
    real(kind=realType), optional, dimension(:), intent(in) :: bcDataValues, bcDataValuesd
    integer(kind=intType), optional, dimension(:, :) :: bcDataFamLists

    ! Output variables:
    real(kind=realType), intent(out), dimension(:) :: dwDot
    real(kind=realType), intent(out), dimension(:, :, :) :: forcesDot

    ! Working Variables
    real(kind=realType), dimension(:, :, :), allocatable :: forces
    integer(kind=intType) :: ierr, nn, sps, mm,i,j,k, l, fSize, ii, jj, iRegion
    real(kind=realType), dimension(nSections) :: t
    real(kind=realType) :: dummyReal, dummyReald

    fSize = size(forcesDot, 2)
    allocate(forces(3, fSize, nTimeIntervalsSpectral))

    call VecPlaceArray(x_like, xdot, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set the provided w and x seeds:
    ii = 0
    jj = 0
    domainLoop1: do nn=1,nDom
       spectalLoop1: do sps=1,nTimeIntervalsSpectral
          call setPointers_d(nn, 1, sps)
          do k=1, kl
             do j=1,jl
                do i=1,il
                   do l=1,3
                      ii = ii + 1
                      xd(i, j, k, l) = xdot(ii)
                   end do
                end do
             end do
          end do
          do k=2, kl
             do j=2,jl
                do i=2,il
                   do l = 1, nw
                      jj = jj + 1
                      wd(i, j, k, l) = wDot(jj)
                   end do
                end do
             end do
          end do
       end do spectalLoop1
    end do domainLoop1

    do sps=1,nTimeIntervalsSpectral
       do nn=1,nDom
          call setPointers_d(nn, 1, sps)
          call xhalo_block_d()
       end do
    end do

    ! Now exchange the coordinates. Note that we *must* exhchange the
    ! actual coordinates as well becuase xhao_block overwrites all
    ! halo nodes and exchange coor corrects them.
    call exchangecoor_d(1)
    call exchangecoor(1)

    do sps=1, nTimeIntervalsSpectral
       ! Update overset connectivity if necessary
       if (oversetPresent .and. &
            (oversetUpdateMode == updateFast .or. &
            oversetUpdateMode == updateFull)) then
          print *,'Full overset update derivatives not implemented'
          !call updateOversetConnectivity_d(1_intType, sps)
       end if
    end do

    ! Now set the xsurfd contribution from the full x perturbation.
    ! scatter from the global seed (in x_like) to xSurfVecd...but only
    ! if wallDistances were used
    if (wallDistanceNeeded .and. useApproxWallDistance) then
       do sps=1, nTimeIntervalsSpectral
          call VecScatterBegin(wallScatter(1, sps), x_like, xSurfVecd(sps), INSERT_VALUES, SCATTER_FORWARD, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          call VecScatterEnd(wallScatter(1, sps), x_like, xSurfVecd(sps),  INSERT_VALUES, SCATTER_FORWARD, ierr)
          call EChk(ierr,__FILE__,__LINE__)
       end do
    end if

    call adjustInflowAngle_d
    call referenceState_d
    if (present(bcDataNames)) then
       do sps=1,nTimeIntervalsSpectral
          call setBCData_d(bcDataNames, bcDataValues, bcDataValuesd, &
               bcDataFamLists, sps, size(bcDataValues), size(bcDataFamLists, 2))
       end do
       call setBCDataFineGrid_d(.true.)
    end if

    do sps=1,nTimeIntervalsSpectral
       do nn=1,nDom

          call setPointers_d(nn, 1, sps)
          ISIZE1OFDrfbcdata = nBocos
          ISIZE1OFDrfviscsubface = nViscBocos

          ! Get the pointers from the petsc vector for the surface
          ! perturbation for wall distance.
          call VecGetArrayF90(xSurfVec(1, sps), xSurf, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! And it's derivative
          call VecGetArrayF90(xSurfVecd(sps), xSurfd, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          call volume_block_d()
          call metric_block_d()
          call boundaryNormals_d()
          if (equations == RANSEquations .and. useApproxWallDistance) then
             call updateWallDistancesQuickly_d(nn, 1, sps)
          end if

          call computePressureSimple_d(.False.)
          call computeLamViscosity_d(.False.)
          call computeEddyViscosity_d(.False.)

          ! Make sure to call the turb BC's first incase we need to
          ! correct for K
          if (equations == RANSequations) then
             call BCTurbTreatment_d
             call applyAllTurbBCthisblock_d(.True.)
          end if

          call applyAllBC_block_d(.True.)

          ! These arrays need to be restored before we can move to the next spectral instance.
          call VecRestoreArrayF90(xSurfVec(1, sps), xSurf, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! And it's derivative
          call VecRestoreArrayF90(xSurfVecd(sps), xSurfd, ierr)
          call EChk(ierr,__FILE__,__LINE__)
       end do
    end do

    ! Just exchange the derivative values.
    call whalo2_d(1, 1, nw, .True., .True., .True.)

    ! Need to re-apply the BCs. The reason is that BC halos behind
    ! interpolated cells need to be recomputed with their new
    ! interpolated values from actual compute cells. Only needed for
    ! overset.
    if (oversetPresent) then
       do sps=1,nTimeIntervalsSpectral
          do nn=1,nDom
             call setPointers_d(nn, 1, sps)
             if (equations == RANSequations) then
                call BCTurbTreatment_d
                call applyAllTurbBCthisblock_d(.True.)
             end if
             call applyAllBC_block_d(.True.)
          end do
       end do
    end if

    do sps=1, nTimeIntervalsSpectral
       do nn=1, nDom
          call setPointers_d(nn, 1, sps)
          ISIZE1OFDrfbcdata = nBocos
          ISIZE1OFDrfviscsubface = nViscBocos

          call timeStep_block_d(.false.)
          dw = zero
          dwd = zero
          ! Compute any source terms
          do iRegion=1, nActuatorRegions
             call sourceTerms_block_d(nn, .True. , iRegion, dummyReal, dummyReald)
          end do

          !Compute turbulence residual for RANS equations
          if( equations == RANSEquations) then
             !call unsteadyTurbSpectral_block(itu1, itu1, nn, sps)

             select case (turbModel)
             case (spalartAllmaras)
                call saSource_d
                call turbAdvection_d(1_intType, 1_intType, itu1-1, qq)
                !!call unsteadyTurbTerm_d(1_intType, 1_intType, itu1-1, qq)
                call saViscous_d
                call saResScale_d
             end select
          end if

          ! compute the mean flow residual
          call inviscidCentralFlux_d

          if (lumpedDiss) then
             select case (spaceDiscr)
             case (dissScalar)
                call inviscidDissFluxScalarApprox_d
             case (dissMatrix)
                call inviscidDissFluxMatrixApprox_d
             case (upwind)
                call inviscidUpwindFlux_d(.True.)
             end select
          else
             select case (spaceDiscr)
             case (dissScalar)
                call inviscidDissFluxScalar_d
             case (dissMatrix)
                call inviscidDissFluxMatrix_d
             case (upwind)
                call inviscidUpwindFlux_d(.True.)
             end select
          end if

          if (viscous) then
             call computeSpeedOfSoundSquared_d
             if (.not. lumpedDiss .or. viscPC) then
                call allNodalGradients_d
                call viscousFlux_d
             else
                call viscousFluxApprox_d
             end if
          end if

          ! if (lowSpeedPreconditioner) then
          !    call applyLowSpeedPreconditioner_d
          ! end if
          call sumDwAndFw_d
          call resscale_d
       end do
    end do

    ! Compute final solution values
    if (present(famLists)) then
       call getSolution_d(famLists, funcValues, funcValuesd)
    end if

    do sps=1, nTimeIntervalsSpectral
       call getForces_d(forces(:, :, sps), forcesDot(:, :, sps), fSize, sps)
    end do

    ! Copy out the residual derivative into the provided dwDot
    ii =0
    do nn=1, nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers_d(nn, 1, sps)
          do k=2, kl
             do j=2, jl
                do i=2, il
                   do l=1, nw
                      ii = ii + 1
                      dwDot(ii) = dwd(i,j,k,l)
                   end do
                end do
             end do
          end do
       end do
    end do
    call VecResetArray(x_like, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    deallocate(forces)
  end subroutine master_d

  subroutine master_b(wbar, xbar, extraBar, forcesBar, dwBar, nState, famLists, &
       funcValues, funcValuesd, bcDataNames, bcDataValues, bcDataValuesd, bcDataFamLists)

    ! This is the main reverse mode differentiaion of master. It
    ! compute the reverse mode sensitivity of *all* outputs with
    ! respect to *all* inputs. Anything that needs to be
    ! differentiated for the adjoint method should be included in this
    ! function. This routine is written by had, assembling the various
    ! individually differentiated tapenade routines.

    use constants
    use adjointVars, only : iAlpha, iBeta, iMach, iMachGrid, iTemperature, iDensity, &
         iPointrefX, iPointRefY, iPointRefZ, iPressure
    use communication, only : adflow_comm_world, myid
    use iteration, only : currentLevel
    use inputAdjoint,  only : viscPC
    use fluxes, only : viscousFlux
    use flowVarRefState, only : nw, nwf, viscous,pInfDimd, rhoInfDimd, TinfDimd
    use blockPointers, only : nDom, il, jl, kl, wd, xd, dw, dwd
    use inputPhysics, only :pointRefd, alphad, betad, equations, machCoefd, &
         machd, machGridd, rgasdimd, equationMode, turbModel, wallDistanceNeeded
    use inputDiscretization, only : lowSpeedPreconditioner, lumpedDiss, spaceDiscr, useAPproxWallDistance
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputAdjoint, only : frozenTurbulence
    use utils, only : isWallType, setPointers_b, EChk
    use adjointPETSc, only : x_like
    use haloExchange, only : whalo2_b, exchangeCoor_b, exchangeCoor, whalo2
    use wallDistanceData, only : xSurfVec, xSurfVecd, xSurf, xSurfd, wallScatter
    use surfaceIntegrations, only : getSolution_b
    use flowUtils, only : fixAllNodalGradientsFromAD
    use adjointextra_b, only : resscale_B, sumdwandfw_b
    use adjointExtra_b, only : xhalo_block_b, volume_block_b, metric_block_b, boundarynormals_b
    use flowutils_b, only : computePressureSimple_b, computeLamViscosity_b, &
         computeSpeedOfSoundSquared_b, allNodalGradients_b, adjustInflowAngle_b
    use solverutils_b, only : timeStep_Block_b
    use turbbcroutines_b, only : applyAllTurbBCthisblock_b,  bcTurbTreatment_b
    use initializeflow_b, only : referenceState_b
    use wallDistance_b, only : updateWallDistancesQuickly_b
    use sa_b, only : saSource_b, saViscous_b, saResScale_b, qq
    use turbutils_b, only : turbAdvection_b, computeEddyViscosity_b
    use residuals_b, only : sourceTerms_block_b
    use fluxes_b, only :inviscidUpwindFlux_b, inviscidDissFluxScalar_b, &
         inviscidDissFluxMatrix_b, viscousFlux_b, inviscidCentralFlux_b
    use BCExtra_b, only : applyAllBC_Block_b
    use bcdata, only : setBCData_b, setBCDataFineGrid_b
    use oversetData, only : oversetPresent
    use inputOverset, only : oversetUpdateMode
    use oversetCommUtilities, only : updateOversetConnectivity_b
    use BCRoutines, only : applyAllBC_block
    use actuatorRegionData, only : nActuatorRegions
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input variables:
    real(kind=realType), intent(in), dimension(:) :: dwBar
    real(kind=realType), intent(in), dimension(:, :, :) :: forcesBar
    integer(kind=intType), intent(in) :: nState
    integer(kind=intType), optional, dimension(:, :), intent(in) :: famLists
    real(kind=realType), optional, dimension(:, :) :: funcValues, funcValuesd
    character, optional, dimension(:, :), intent(in) :: bcDataNames
    real(kind=realType), optional, dimension(:), intent(in) :: bcDataValues
    integer(kind=intType), optional, dimension(:, :) :: bcDataFamLists

    ! Output Arguments:
    real(kind=realType), optional, intent(out), dimension(:) :: wBar, xBar, extraBar
    real(kind=realType), optional, dimension(:), intent(out) :: bcDataValuesd

    ! Working Variables
    integer(kind=intType) :: ierr, nn, sps, mm,i,j,k, l, fSize, ii, jj,  level, iRegion
    real(kind=realType), dimension(:), allocatable :: extraLocalBar, bcDataValuesdLocal
    real(kind=realType) :: dummyReal, dummyReald
    logical ::resetToRans

    ! extraLocalBar accumulates the seeds onto the extra variables
    allocate(extraLocalBar(size(extrabar)))
    extraLocalBar = zero

    ! Place the output spatial seed into the temporary petsc x-like
    ! vector.
    xBar = zero
    call VecPlaceArray(x_like, xBar, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Set the residual seeds.
    ii = 0
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral

          ! Set pointers and derivative pointers
          call setPointers_b(nn, 1, sps)

          ! Set the dw seeds
          do k=2, kl
             do j=2, jl
                do i=2, il
                   do l=1, nState
                      ii = ii + 1
                      dwd(i, j, k, l) = dwbar(ii)
                   end do
                end do
             end do
          end do
       end do
    end do

    forceSpsLoop: do sps=1, nTimeIntervalsSpectral
       fSize = size(forcesBar, 2)
       call getForces_b(forcesBar(:, :, sps), fSize, sps)
    end do forceSpsLoop

    ! Call the final getSolution_b routine
    if (present(famLists)) then
       call getSolution_b(famLists, funcValues, funcValuesd)
    end if

    spsLoop1: do sps=1, nTimeIntervalsSpectral

       domainLoop1: do nn=1, nDom
          call setPointers_b(nn, 1, sps)

          ! Now we start running back through the main residual code:
          call resScale_b
          call sumDwAndFw_b

          ! if (lowSpeedPreconditioner) then
          !    call applyLowSpeedPreconditioner_b
          ! end if

          ! Note that master_b does not include the first order flux
          ! approxation codes as those are never needed in reverse.
          if (viscous) then
             call viscousFlux_b
             call allNodalGradients_b
             call computeSpeedOfSoundSquared_b
          end if

          ! So the all nodal gradients doesnt' perform the final
          ! scaling by the volume since it isn't necessary for the
          ! derivative. We have a special routine to fix that.
          if (viscous) then
             call fixAllNodalGradientsFromAD()
             call viscousFlux
          end if

          select case (spaceDiscr)
          case (dissScalar)
             call inviscidDissFluxScalar_b
          case (dissMatrix)
             call inviscidDissFluxMatrix_b
          case (upwind)
             call inviscidUpwindFlux_b(.True.)
          end select

          call inviscidCentralFlux_b
          call timeStep_block_b(.false.)
          ! Compute turbulence residual for RANS equations
          if( equations == RANSEquations) then
             select case (turbModel)
             case (spalartAllmaras)
                call saResScale_b
                call saViscous_b
                !call unsteadyTurbTerm_b(1_intType, 1_intType, itu1-1, qq)
                call turbAdvection_b(1_intType, 1_intType, itu1-1, qq)
                call saSource_b
             end select

             !call unsteadyTurbSpectral_block_b(itu1, itu1, nn, sps)
          end if

          ! Just to be safe, zero the pLocald value...should not matter
          dummyReald = zero
          do iRegion=1, nActuatorRegions
             call sourceTerms_block_b(nn, .True. , iRegion, dummyReal, dummyReald)
          end do

          ! Initres_b should be called here. For steady just zero:
          dwd = zero
       end do domainLoop1
    end do spsLoop1

    ! Need to re-apply the BCs. The reason is that BC halos behind
    ! interpolated cells need to be recomputed with their new
    ! interpolated values from actual compute cells. Only needed for
    ! overset.
    if (oversetPresent) then
       do sps=1, nTimeIntervalsSpectral
          do nn=1,nDom
             call setPointers_b(nn, 1, sps)
             call applyAllBC_block_b(.True.)
             call applyAllBC_block(.true.)

             if (equations == RANSequations) then
                call applyAllTurbBCThisBlock_b(.True.)
                call bcTurbTreatment_b
             end if
          end do
       end do
    end if

    ! Exchange the adjoint values.
    call whalo2_b(currentLevel, 1_intType, nw, .True., .True., .True.)

    spsLoop2: do sps=1,nTimeIntervalsSpectral

       ! Get the pointers from the petsc vector for the wall
       ! surface and it's accumulation. Only necessary for wall
       ! distance.
       call VecGetArrayF90(xSurfVec(1, sps), xSurf, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! And it's derivative
       call VecGetArrayF90(xSurfVecd(sps), xSurfd, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       !Zero the accumulation vector on a per-time-spectral instance basis
       xSurfd = zero

       domainLoop2: do nn=1,nDom
          call setPointers_b(nn, 1, sps)
          call applyAllBC_block_b(.True.)

          ! Run the forward application of the BCs. The reason is that
          ! the reverse application of the BCs can result in
          ! inconsisent values in the halos. This has only been
          ! observed to caluse problems with hot subsonic flow in an
          ! engine core. This is the same reason we need the applyBCs
          ! after the _b version above.
          call applyAllBC_block(.true.)

          if (equations == RANSequations) then
             call applyAllTurbBCThisBlock_b(.True.)
             call bcTurbTreatment_b
          end if
          call computeEddyViscosity_b(.false.)
          call computeLamViscosity_b(.false.)
          call computePressureSimple_b(.false.)

          if (equations == RANSEquations .and. useApproxWallDistance) then
             call updateWallDistancesQuickly_b(nn, 1, sps)
          end if

          call boundaryNormals_b
          call metric_block_b
          call volume_block_b

       end do domainLoop2

       ! Restore the petsc pointers.
       call VecGetArrayF90(xSurfVec(1, sps), xSurf, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! And it's derivative
       call VecGetArrayF90(xSurfVecd(sps), xSurfd, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Now accumulate the xsurfd accumulation by using the wall scatter
       ! in reverse.
       if (wallDistanceNeeded .and. useApproxWallDistance) then

          call VecScatterBegin(wallScatter(1, sps), xSurfVecd(sps), x_like, ADD_VALUES, SCATTER_REVERSE, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          call VecScatterEnd(wallScatter(1, sps), xSurfVecd(sps), x_like, ADD_VALUES, SCATTER_REVERSE, ierr)
          call EChk(ierr,__FILE__,__LINE__)
       end if
    end do spsLoop2


    if (present(bcDataNames)) then
       allocate(bcDataValuesdLocal(size(bcDataValuesd)))
       bcDataValuesdLocal = zero
       call setBCDataFineGrid_b(.true.)
       do sps=1, nTimeIntervalsSpectral
          call setBCData_b(bcDataNames, bcDataValues, bcDataValuesdLocal, bcDataFamLists, &
               sps, size(bcDataValues), size(bcDataFamLIsts, 2))
       end do
       ! Reverse seeds need to accumulated across all processors:
       call mpi_allreduce(bcDataValuesdLocal, bcDataValuesd, size(bcDataValuesd), adflow_real, &
            mpi_sum, ADflow_comm_world, ierr)
       deallocate(bcDataValuesdLocal)
    end if
    call referenceState_b
    call adjustInflowAngle_b

    do sps=1, nTimeIntervalsSpectral
       ! Update overset connectivity if necessary
       if (oversetPresent .and. &
            (oversetUpdateMode == updateFast .or. &
            oversetUpdateMode == updateFull)) then
          print *,'Full overset update derivatives not implemented'
          !call updateOversetConnectivity_b(1_intType, sps)
       end if
    end do
    ! Now the adjoint of the coordinate exhcange
    call exchangecoor_b(1)
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers_b(nn, 1, sps)
          call xhalo_block_b()
       end do
    end do

    call VecResetArray(x_like, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! =========================================
    !            End of reverse pass
    ! =========================================

    ! Store the extra derivatives
    extraLocalBar(iAlpha) = alphad
    extraLocalBar(iBeta) =  betad
    extraLocalBar(iMach) =  machd + machcoefd
    extraLocalBar(iMachGrid) = machgridd
    extraLocalBar(iPressure) =  pinfdimd
    extraLocalBar(iTemperature) = tinfdimd
    extraLocalBar(iDensity) =  rhoinfdimd
    extraLocalBar(iPointRefX) = pointrefd(1)
    extraLocalBar(iPointRefY) = pointrefd(2)
    extraLocalBar(iPointRefZ) = pointrefd(3)

    ! Finally put the output seeds into the provided vectors.
    ii = 0
    jj = 0
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral

          ! Set pointers and derivative pointers
          call setPointers_b(nn, 1, sps)
          ! Set the wbar accumulation
          do k=2, kl
             do j=2, jl
                do i=2, il
                   do l=1, nState
                      ii = ii + 1
                      wbar(ii) = wd(i, j, k, l)
                   end do
                end do
             end do
          end do

          ! Set the xvbar accumulation. Note that this must be a sum,
          ! becuase we may already have wall distance accumulation
          ! from the wallScatter directly into xbar (through x_like).
          do k=1, kl
             do j=1, jl
                do i=1, il
                   do l=1, 3
                      jj = jj + 1
                      xbar(jj) = xbar(jj) + xd(i, j, k, l)
                   end do
                end do
             end do
          end do
       end do
    end do

    ! Finally get the full contribution of the extra variables by
    ! summing all local contributions.
    extraBar = zero
    call mpi_allreduce(extraLocalBar, extraBar, size(extraBar), adflow_real, &
         mpi_sum, ADflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

  end subroutine master_b

  subroutine master_state_b(wbar, dwBar, nState)

    ! This is specialized form of master that *ONLY* computes drdw
    ! products. It uses a few specialzed routines that are
    ! differentiated without including spatial dependencies. This
    ! results in slightly faster code. This specialization is
    ! justififed since this routine is needed for the transpose
    ! matrix-vector products in solving the adjoint system and thus
    ! this routine is called several orders of magniutde more than
    ! master_b. This routine has to be fast!

    use constants
    use iteration, only : currentLevel
    use flowVarRefState, only : nw, viscous
    use blockPointers, only : nDom, il, jl, kl, wd, dwd, iblank
    use inputPhysics, only : equationMode, turbModel, equations
    use inputDiscretization, only : lowSpeedPreconditioner, spaceDiscr
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use utils, only : setPointers_d
    use haloExchange, only : whalo2_b
    use flowUtils, only : fixAllNodalGradientsFromAD
    use adjointextra_b, only : resscale_B, sumdwandfw_b
    use flowutils_b, only : computePressureSimple_b, computeLamViscosity_b, &
         computeSpeedOfSoundSquared_b
    use turbbcroutines_b, only : applyAllTurbBCthisblock_b,  bcTurbTreatment_b
    use turbUtils_b, only : computeEddyViscosity_b
    use BCExtra_b, only : applyAllBC_Block_b

    use sa_fast_b, only : saresscale_fast_b, saviscous_fast_b, &
         sasource_fast_b, qq
    use turbutils_fast_b, only : turbAdvection_fast_b
    use fluxes_fast_b, only :inviscidUpwindFlux_fast_b, inviscidDissFluxScalar_fast_b, &
         inviscidDissFluxMatrix_fast_b, viscousFlux_fast_b, inviscidCentralFlux_fast_b
    use solverutils_fast_b, only : timeStep_block_fast_b
    use flowutils_fast_b, only : allnodalgradients_fast_b
    use residuals_fast_b, only : sourceTerms_block_fast_b
    use oversetData, only : oversetPresent
    use bcroutines, only : applyallbc_block
    use actuatorRegionData, only : nActuatorRegions
    implicit none

    ! Input variables:
    real(kind=realType), intent(in), dimension(:) :: dwBar
    integer(kind=intType), intent(in) :: nState

    ! Input Arguments:
    real(kind=realType), intent(out), dimension(:) :: wBar

    ! Working Variables
    integer(kind=intType) :: ierr, nn, sps, mm,i,j,k, l, fSize, ii, jj,  level, iRegion
    real(kind=realType) :: dummyReal

    ! Set the residual seeds.
    ii = 0
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers_d(nn, 1, sps)
          do k=2, kl
             do j=2, jl
                do i=2, il
                   do l=1, nState
                      ii = ii + 1
                      dwd(i, j, k, l) = dwbar(ii)
                   end do
                end do
             end do
          end do
       end do
    end do

    ! ============================================
    !  reverse the order of calls from master
    ! ============================================

    spsLoop1: do sps=1, nTimeIntervalsSpectral
       domainLoop1: do nn=1, nDom
          call setPointers_d(nn, 1, sps)

          ! Now we start running back through the main residual code:
          call resScale_b
          call sumDwAndFw_b

          ! if (lowSpeedPreconditioner) then
          !    call applyLowSpeedPreconditioner_b
          ! end if

          ! Note that master_b does not include the approximation codes
          ! as those are never needed in reverse.
          if (viscous) then
             call viscousFlux_fast_b
             call allNodalGradients_fast_b
             call computeSpeedOfSoundSquared_b
          end if

          select case (spaceDiscr)
          case (dissScalar)
             call inviscidDissFluxScalar_fast_b
          case (dissMatrix)
             call inviscidDissFluxMatrix_fast_b
          case (upwind)
             call inviscidUpwindFlux_fast_b(.True.)
          end select

          call inviscidCentralFlux_fast_b

          ! Compute turbulence residual for RANS equations
          if( equations == RANSEquations) then
             select case (turbModel)
             case (spalartAllmaras)
                call saResScale_fast_b
                call saViscous_fast_b
                !call unsteadyTurbTerm_b(1_intType, 1_intType, itu1-1, qq)
                call turbAdvection_fast_b(1_intType, 1_intType, itu1-1, qq)
                call saSource_fast_b
             end select

             !call unsteadyTurbSpectral_block_b(itu1, itu1, nn, sps)
          end if

          call timeStep_block_fast_b(.false.)
          do iRegion=1, nActuatorRegions
             call sourceTerms_block_fast_b(nn, .True. , iRegion, dummyReal)
          end do

          ! Initres_b should be called here. For steady just zero:
          dwd = zero
       end do domainLoop1
    end do spsLoop1

    ! Need to re-apply the BCs. The reason is that BC halos behind
    ! interpolated cells need to be recomputed with their new
    ! interpolated values from actual compute cells. Only needed for
    ! overset.
    if (oversetPresent) then
       do sps=1, nTimeIntervalsSpectral
          do nn=1,nDom
             call setPointers_d(nn, 1, sps)
             call applyAllBC_block_b(.True.)
             call applyAllBC_block(.true.)

             if (equations == RANSequations) then
                call applyAllTurbBCThisBlock_b(.True.)
                call bcTurbTreatment_b
             end if
          end do
       end do
    end if

    ! Exchange the adjoint values.
    call whalo2_b(currentLevel, 1_intType, nw, .True., .True., .True.)

    spsLoop2: do sps=1,nTimeIntervalsSpectral
       domainLoop2: do nn=1,nDom
          call setPointers_d(nn, 1, sps)

          call applyAllBC_block_b(.True.)
          call applyAllBC_block(.true.)

          if (equations == RANSequations) then
             call applyAllTurbBCThisBlock_b(.True.)
             call bcTurbTreatment_b
          end if

          call computeEddyViscosity_b(.false.)
          call computeLamViscosity_b(.false.)
          call computePressureSimple_b(.false.)

       end do domainLoop2
    end do spsLoop2

    ! Finally put the output seeds into wbar
    ii = 0
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers_d(nn, 1, sps)
          do k=2, kl
             do j=2, jl
                do i=2, il
                   do l=1, nState
                      ii = ii + 1
                      wbar(ii) = wd(i, j, k, l)!*max(real(iblank(i,j,k)), zero)
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine master_state_b
#endif
  subroutine block_res_state(nn, sps, useFlowRes, useTurbRes)

    ! This is a special state-only routine used only for finite
    ! differce computations of the jacobian
    use constants
    use BCRoutines, only : applyAllBC_Block
    use inputAdjoint,  only : viscPC
    use blockPointers, only : nDom, wd, xd, dw, il, jl, kl
    use flowVarRefState, only : viscous
    use inputPhysics, only : equations, turbModel
    use inputDiscretization, only : lowSpeedPreconditioner, lumpedDiss, spaceDiscr
    use utils, only :  setPointers, EChk
    use residuals, only : sourceTerms_block
    use turbutils, only : computeEddyViscosity
    use turbbcroutines, only : applyAllTurbBCthisblock,  bcTurbTreatment
    use adjointExtra, only :  sumdwandfw, resScale
    use inputDiscretization, only : useBlockettes
    use blockette, only : blocketteResCore, blockResCore
    use flowUtils, only : computeLamViscosity, computePressureSimple
    use actuatorRegionData, only : nActuatorRegions
    implicit none

    ! Input Arguments:
    integer(kind=intType), intent(in) :: nn, sps
    logical, optional ,intent(in) :: useFlowRes, useTurbRes

    ! Working Variables
    integer(kind=intType) :: ierr, mm,i,j,k, l, fSize, ii, jj, iRegion
    real(kind=realType) ::  pLocal
    logical :: dissApprox, viscApprox, updateIntermed, flowRes, turbRes, storeWall

    flowRes = .True.
    if (present(useFlowRes)) then
       flowRes = useFlowRes
    end if

    turbRes = .True.
    if (present(useTurbRes)) then
       turbRes = useTurbRes
    end if

    call computePressureSimple(.True.)
    call computeLamViscosity(.True.)
    call computeEddyViscosity(.True.)

    ! Make sure to call the turb BC's first incase we need to
    ! correct for K
    if (equations == RANSequations) then
       call BCTurbTreatment
       call applyAllTurbBCthisblock(.True.)
    end if

    call applyAllBC_block(.True.)

    ! Set the flags:
    dissApprox = lumpedDiss
    viscApprox = lumpedDiss
    updateIntermed = .False.
    storeWall = .True.
    blockettes: if (useBlockettes) then
       call blocketteResCore(dissApprox, viscApprox, updateIntermed, flowRes, turbRes, storeWall)
    else
       call blockResCore(dissApprox, viscApprox, updateIntermed, flowRes, turbRes, storeWall, nn, sps)
    end if blockettes
    do iRegion=1, nActuatorRegions
       call sourceTerms_block(nn, .True., iRegion, pLocal)
    end do
    call resscale

  end subroutine block_res_state
#ifndef USE_COMPLEX
  subroutine block_res_state_d(nn, sps)

    ! This is a special state-only forward mode linearization
    ! computation used to assemble the jacobian.
    use constants
    use BCExtra_d, only : applyAllBC_Block_d
    use inputAdjoint,  only : viscPC
    use blockPointers, only : nDom, wd, xd, dw, dwd
    use flowVarRefState, only : viscous
    use inputPhysics, only : equations, turbModel
    use inputDiscretization, only : lowSpeedPreconditioner, lumpedDiss, spaceDiscr
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use utils, only : setPointers_d, EChk
    use sa_d, only : saSource_d, saViscous_d, saResScale_d, qq
    use turbutils_d, only : turbAdvection_d, computeEddyViscosity_d
    use fluxes_d, only :inviscidDissFluxScalarApprox_d, inviscidDissFluxMatrixApprox_d, &
         inviscidUpwindFlux_d, inviscidDissFluxScalar_d, inviscidDissFluxMatrix_d, &
         inviscidUpwindFlux_d, viscousFlux_d, viscousFluxApprox_d, inviscidCentralFlux_d
    use flowutils_d, only : computePressureSimple_d, computeLamViscosity_d, &
         computeSpeedOfSoundSquared_d, allNodalGradients_d
    use solverutils_d, only : timeStep_Block_d
    use turbbcroutines_d, only : applyAllTurbBCthisblock_d,  bcTurbTreatment_d
    use adjointextra_d, only : resscale_D, sumdwandfw_d
    use residuals_d, only: sourceterms_block_d
    use actuatorRegionData, only : nActuatorRegions
    implicit none

    ! Input Arguments:
    integer(kind=intType), intent(in) :: nn, sps

    ! Working Variables
    integer(kind=intType) :: ierr, mm,i,j,k, l, fSize, ii, jj, iRegion
    real(kind=realType) :: dummyReal, dummyReald

    call computePressureSimple_d(.True.)
    call computeLamViscosity_d(.True.)
    call computeEddyViscosity_d(.True.)

    ! Make sure to call the turb BC's first incase we need to
    ! correct for K
    if (equations == RANSequations) then
       call BCTurbTreatment_d
       call applyAllTurbBCthisblock_d(.True.)
    end if

    call applyAllBC_block_d(.True.)
    call timeStep_block_d(.false.)

    dw = zero ! These two lines are init_res
    dwd = zero

    ! Compute any source terms
    do iRegion=1, nActuatorRegions
       call sourceTerms_block_d(nn, .True. , iRegion, dummyReal, dummyReald)
    end do

    !Compute turbulence residual for RANS equations
    if( equations == RANSEquations) then
       !call unsteadyTurbSpectral_block(itu1, itu1, nn, sps)

       select case (turbModel)
       case (spalartAllmaras)
          call saSource_d
          call turbAdvection_d(1_intType, 1_intType, itu1-1, qq)
          !!call unsteadyTurbTerm_d(1_intType, 1_intType, itu1-1, qq)
          call saViscous_d
          call saResScale_d
       end select
    end if

    ! compute the mean flow residual
    call inviscidCentralFlux_d

    if (lumpedDiss) then
       select case (spaceDiscr)
       case (dissScalar)
          call inviscidDissFluxScalarApprox_d
       case (dissMatrix)
          call inviscidDissFluxMatrixApprox_d
       case (upwind)
          call inviscidUpwindFlux_d(.True.)
       end select
    else
       select case (spaceDiscr)
       case (dissScalar)
          call inviscidDissFluxScalar_d
       case (dissMatrix)
          call inviscidDissFluxMatrix_d
       case (upwind)
          call inviscidUpwindFlux_d(.True.)
       end select
    end if

    if (viscous) then
       call computeSpeedOfSoundSquared_d
       if (.not. lumpedDiss .or. viscPC) then
          call allNodalGradients_d
          call viscousFlux_d
       else
          call viscousFluxApprox_d
       end if
    end if

    ! if (lowSpeedPreconditioner) then
    !    call applyLowSpeedPreconditioner_d
    ! end if
    call sumDwAndFw_d
    call resscale_d
  end subroutine block_res_state_d
#endif

end module masterRoutines
