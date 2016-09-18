module masterRoutines
contains
  subroutine master(forces, fSize)

    use constants
    use communication, only : sumb_comm_world
    use BCRoutines, only : applyallBC_block
    use turbbcRoutines, only : applyallTurbBCthisblock, bcTurbTreatment
    use iteration, only : currentLevel
    use inputAdjoint,  only : viscPC
    use flowVarRefState, only : nwf, nw
    use costFunctions, only : nLocalValues
    use blockPointers, only : nDom
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
    use residuals, only : initRes_block
    use surfaceIntegrations, only : integrateSurfaces
    use adjointExtra, only : volume_block, metric_block, boundaryNormals,&
         xhalo_block, sumdwandfw, resScale, getCostFunctions
    implicit none

    ! Input Arguments:
    integer(kind=intType), intent(in) ::  fSize

    ! Output Variables
    real(kind=realType), intent(out),  dimension(3, fSize, nTimeIntervalsSpectral) :: forces

    ! Working Variables
    integer(kind=intType) :: ierr, nn, sps
    real(kind=realType), dimension(nSections) :: t
    real(kind=realType), dimension(nLocalValues, nTimeIntervalsSpectral) :: localVal, globalVal

    call adjustInflowAngle()
    call referenceState

    do sps=1,nTimeIntervalsSpectral
       do nn=1,nDom
          call setPointers(nn, 1, sps)
          call xhalo_block()
       end do
    end do

    ! Now exchange the coordinates (fine level only)
    call exchangecoor(1)
    do sps=1,nTimeIntervalsSpectral
       do nn=1,nDom
          call setPointers(nn, 1, sps)

          call volume_block
          call metric_block
          call boundaryNormals

          if (equations == RANSEquations .and. useApproxWallDistance) then 
             call updateWallDistancesQuickly(nn, 1, sps)
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

    ! Main loop for the residual
    do sps=1, nTimeIntervalsSpectral
       do nn=1, nDom
          call setPointers(nn, 1, sps)
          call initRes_block(1, nw, nn, sps)

          ! Compute turbulence residual for RANS equations
          if( equations == RANSEquations) then

             ! Initialize only the Turblent Variables
             !call unsteadyTurbSpectral_block(itu1, itu1, nn, sps)

             select case (turbModel)

             case (spalartAllmaras)
                call saSource
                call turbAdvection(1_intType, 1_intType, itu1-1, qq)
                !call unsteadyTurbTerm(1_intType, 1_intType, itu1-1, qq)
                call saViscous
                call saResScale
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

          ! Now compute the forces and moments for this block. 
          call integrateSurfaces(localval(:, sps))

          ! Also comptue mass flows
          !call massFlux(localVal)
       end do

       ! Now we can retrieve the forces/tractions for this spectral instance
       call getForces(forces(:, :, sps), fSize, sps)
    end do

    ! Now we need to reduce all the cost functions
    call mpi_allreduce(localval, globalVal, nLocalValues*nTimeIntervalsSpectral, sumb_real, &
         MPI_SUM, sumb_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Call the final routine that will comptue all of our functions of
    ! interest.
    call getCostFunctions(globalVal)

  end subroutine master

  subroutine master_d(wdot, xdot, forcesDot, dwDot)
    use constants
    use costFunctions
    use diffsizes, only :  ISIZE1OFDrfbcdata, ISIZE1OFDrfviscsubface
    use communication, only : sumb_comm_world
    use iteration, only : currentLevel
    use BCExtra_d, only : applyAllBC_Block_d
    use inputAdjoint,  only : viscPC
    use flowVarRefState, only : nw, nwf
    use blockPointers, only : nDom, il, jl, kl, wd, xd, dw, dwd, nBocos, BCType, nViscBocos, BCData
    use flowVarRefState, only : viscous, timerefd
    use inputPhysics , only : turbProd, equationMode, equations, turbModel, wallDistanceNeeded
    use inputDiscretization, only : lowSpeedPreconditioner, lumpedDiss, spaceDiscr, useAPproxWallDistance
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use section, only: sections, nSections
    use monitor, only : timeUnsteadyRestart
    use utils, only : isWallType, setPointers, setPointers_d, EChk
    use utils, only : setBCPointers_d
    use sa_d, only : saSource_d, saViscous_d, saResScale_d, qq
    use turbutils_d, only : turbAdvection_d, computeEddyViscosity_d
    use fluxes_d, only :inviscidDissFluxScalarApprox_d, inviscidDissFluxMatrixApprox_d, &
         inviscidUpwindFlux_d, inviscidDissFluxScalar_d, inviscidDissFluxMatrix_d, &
         inviscidUpwindFlux_d, viscousFlux_d, viscousFluxApprox_d, inviscidCentralFlux_d
    use adjointPETSc, only : x_like
    use haloExchange, only : whalo2_d, exchangeCoor_d, exchangeCoor, whalo2
    use wallDistance_d, only : updateWallDistancesQuickly_d
    use wallDistanceData, only : xSurfVec, xSurfVecd, xSurf, xSurfd, wallScatter
    use flowutils_d, only : computePressureSimple_d, computeLamViscosity_d, &
         computeSpeedOfSoundSquared_d, allNodalGradients_d, adjustInflowAngle_d
    use solverutils_d, only : timeStep_Block_d
    use turbbcroutines_d, only : applyAllTurbBCthisblock_d,  bcTurbTreatment_d
    use initializeflow_d, only : referenceState_d
    use surfaceIntegrations_d, only : forcesAndMomentsFace_d
    use surfaceFamilies, only : famGroups
    use sorting, only : bsearchIntegers
    use adjointExtra_d, only : xhalo_block_d, volume_block_d, metric_BLock_d, boundarynormals_d
    use adjointextra_d, only : getcostfunctions_D, resscale_D, sumdwandfw_d

    implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"

    ! Input Arguments:
    real(kind=realType), intent(in), dimension(:) :: wDot, xDot

    ! Output variables:
    real(kind=realType), intent(out), dimension(:) :: dwDot
    real(kind=realType), intent(out), dimension(:, :, :) :: forcesDot

    ! Working Variables
    real(kind=realType), dimension(:, :, :), allocatable :: forces
    integer(kind=intType) :: ierr, nn, sps, mm,i,j,k, l, fSize, ii, jj
    real(kind=realType), dimension(nSections) :: t
    real(kind=realType), dimension(nLocalValues, nTimeIntervalsSpectral) :: localVal, globalVal
    real(kind=realType), dimension(nLocalValues, nTimeIntervalsSpectral) :: localVald, globalVald

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
    
    localVal = zero
    localVald = zero

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

    do sps=1, nTimeIntervalsSpectral
       do nn=1, nDom
          call setPointers_d(nn, 1, sps)
          ISIZE1OFDrfbcdata = nBocos
          ISIZE1OFDrfviscsubface = nViscBocos

          call timeStep_block_d(.false.)
          dw = zero
          dwd = zero

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
                          
          ! Now compute the forces and moments for this block. 
          do mm=1, nBocos
             ! Determine if this boundary condition is to be incldued in the
             ! currently active group
             famInclude: if (bsearchIntegers(BCdata(mm)%famID, &
                  famGroups, size(famGroups)) > 0) then
                
                ! Set a bunch of pointers depending on the face id to make
                ! a generic treatment possible. 
                call setBCPointers_d(mm, .True.)
                
                isWall: if( isWallType(BCType(mm))) then 
                   call forcesAndMomentsFace_d(localVal, localVald, mm)
                end if isWall
             end if famInclude
          end do
       end do
       ! Now we can retrieve the forces/tractions for this spectral instance
       call getForces_d(forces(:, :, sps), forcesDot(:, :, sps), fSize, sps)
    end do

    ! Now we need to reduce all the cost functions
    call mpi_allreduce(localval, globalVal, nLocalValues*nTimeIntervalsSpectral, sumb_real, &
         MPI_SUM, sumb_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call mpi_allreduce(localvald, globalVald, nLocalValues*nTimeIntervalsSpectral, sumb_real, &
         MPI_SUM, sumb_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Call the final routine that will comptue all of our functions of
    ! interest.
    call getCostFunctions_d(globalVal, globalVald)

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

  subroutine block_res_state(nn, sps)
    
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
    use sa, only : sa_block, saSource, saViscous, saResScale, qq
    use turbutils, only : turbAdvection, computeEddyViscosity
    use fluxes, only :inviscidDissFluxScalarApprox, inviscidDissFluxMatrixApprox, &
         inviscidUpwindFlux, inviscidDissFluxScalar, inviscidDissFluxMatrix, &
         inviscidUpwindFlux, viscousFlux, viscousFluxApprox, inviscidCentralFlux
    use flowutils, only : computePressureSimple, computeLamViscosity, &
         computeSpeedOfSoundSquared, allNodalGradients
    use solverutils, only : timeStep_Block
    use turbbcroutines, only : applyAllTurbBCthisblock,  bcTurbTreatment
    use adjointExtra, only :  sumdwandfw, resScale
    use iteration
    implicit none

    ! Input Arguments:
    integer(kind=intType), intent(in) :: nn, sps

    ! Working Variables
    integer(kind=intType) :: ierr, mm,i,j,k, l, fSize, ii, jj
  
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
    call timeStep_block(.false.)
    dw = zero

    !Compute turbulence residual for RANS equations
    if( equations == RANSEquations) then
       !call unsteadyTurbSpectral_block(itu1, itu1, nn, sps)
       
       select case (turbModel)
       case (spalartAllmaras)
          call sa_block(.True.)
          allocate(qq(2:il,2:jl,2:kl))
          call saSource
          call turbAdvection(1_intType, 1_intType, itu1-1, qq)
          !!call unsteadyTurbTerm_d(1_intType, 1_intType, itu1-1, qq)
          call saViscous
          call saResScale
          deallocate(qq)
       end select
    end if

    rFil = one
    ! compute the mean flow residual
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
    
    ! if (lowSpeedPreconditioner) then 
    !    call applyLowSpeedPreconditioner_d
    ! end if
    call sumDwAndFw
    call resscale

  end subroutine block_res_state

  subroutine block_res_state_d(nn, sps)

    ! This is a special state-only forward mode linearization
    ! computation used to assemble the jacobian. 
    use constants
    use costFunctions
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
    implicit none

    ! Input Arguments:
    integer(kind=intType), intent(in) :: nn, sps

    ! Working Variables
    integer(kind=intType) :: ierr, mm,i,j,k, l, fSize, ii, jj

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
    dw = zero
    dwd = zero

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


end module masterRoutines
