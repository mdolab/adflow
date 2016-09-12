module adjointAPI

contains

#ifndef USE_COMPLEX
  subroutine computeMatrixFreeProductFwd(xvdot, extradot, wdot, useSpatial, useState, dwdot, funcsDot, &
       fDot, spatialSize, extraSize, stateSize, costSize, fSize, nTime)

    ! This is the main matrix-free forward mode computation
    use constants
    use block, only : flowDomsd
    use communication, only : sumb_comm_world, myid
    use costfunctions
    use blockPointers
    use inputDiscretization 
    use inputTimeSpectral 
    use inputPhysics
    use iteration         
    use flowVarRefState     
    use inputAdjoint 
    use adjointvars
    use stencils
    use diffSizes
    use wallDistanceData, only : xSurfVec, xSurfVecd, xSurf, xSurfd, wallScatter
    use adjointPETSc, only : x_like
    use surfaceFamilies, only: wallFamilies, totalWallFamilies
    use utils, only : setPointers, EChk, getDirAngle, setPointers_d
    use haloExchange, only : whalo2_d, exchangeCoor_d
    use adjointextra_d, only : xhalo_block_d, block_res_d
    use adjointUtils, only : allocDerivativeValues, zeroADSeeds
    implicit none

#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"

    ! Input Variables
    integer(kind=intType), intent(in) :: spatialSize, extraSize, stateSize, costSize, fSize, nTime
    real(kind=realType), dimension(spatialSize), intent(in) :: xvdot
    real(kind=realType), dimension(extraSize), intent(in) :: extradot
    real(kind=realType), dimension(stateSize), intent(in) :: wdot
    logical, intent(in) :: useSpatial, useState

    ! Ouput Variables
    real(kind=realType), dimension(stateSize), intent(out) :: dwDot
    real(kind=realType), dimension(costSize), intent(out) :: funcsDot
    real(kind=realType), dimension(3, fSize, nTime), intent(out) :: fDot

    ! Working Variables
    real(kind=realType), dimension(3, fSize) :: forces
    integer(kind=intType) :: ierr,nn,mm,sps,i,j,k,l,ii,jj,idim,sps2
    integer(kind=intType) ::  level, irow
    real(kind=realType), dimension(costSize) :: funcsLocalDot
    logical :: resetToRans

    ! Determine if we want to use frozenTurbulent Adjoint
    resetToRANS = .False. 
    if (frozenTurbulence .and. equations == RANSEquations) then
       equations = NSEquations 
       resetToRANS = .True.
    end if

    call VecPlaceArray(x_like, xvdot, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    
    ! Need to trick the residual evalution to use coupled (mean flow and
    ! turbulent) together.
    level = 1
    currentLevel = level
    groundLevel = level

    ! Allocate the memory we need for this block to do the forward
    ! mode derivatives and copy reference values
    ! Allocate the memory for reverse
    if (.not. derivVarsAllocated) then 
       call allocDerivativeValues(level)
    end if
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn, level, sps)
          call zeroADSeeds(nn,level, sps)
       end do
    end do

    ! All arrays are zeroed in alloc_deriv_values, but we still need to
    ! zero the extra variables here
    alphad = zero
    betad  = zero
    machd  = zero
    machGridd = zero
    machcoefd = zero
    pointRefd  = zero
    lengthRefd = zero
    pinfdimd = zero
    tinfdimd = zero
    rhoinfdimd = zero
    rgasdimd = zero

    if (useSpatial) then 
       ! Here we set the spatial and extra seeds if necessary.
       ii = 0
       domainLoop1: do nn=1,nDom

          ! Just to get sizes
          call setPointers_d(nn, level, 1)

          spectalLoop1: do sps=1,nTimeIntervalsSpectral
             do k=1, kl
                do j=1,jl
                   do i=1,il
                      do l=1,3
                         ii = ii + 1
                         flowdomsd(nn,level,sps)%x(i, j, k, l) = xvdot(ii)
                      end do
                   end do
                end do
             end do

             ! Do the xhalo stuff that has to be done before the main call to block_res
             call xhalo_block_d()

          end do spectalLoop1
       end do domainLoop1

       ! Now run the halo exchange for the nodes
       call exchangecoor_d(level)

       ! And now do the extra ones
       ! Set the seeds we have, this is zero-based:
       if (nDesignAoA >= 0) &
            alphad = extraDot(nDesignAoA+1)
       if (nDesignSSA >= 0) &
            betad = extraDot(nDesignSSA+1)
       if (nDesignMach  >= 0) then
          machd = extraDot(nDesignMach+1)
          machCoefd = extraDot(nDesignMach+1)
       end if
       if (nDesignMachGrid >= 0) then
          machGridd = extraDot(nDesignMachGrid+1)
          machCoefd = extraDot(nDesignMachGrid+1)
       end if
       if (nDesignPressure >= 0) &
            PinfDimd = extraDot(nDesignPressure+1)
       if (nDesignDensity >= 0) &
            rhoinfDimd = extraDot(nDesignDensity+1)
       if (nDesignTemperature >= 0) &
            tinfdimd = extraDot(nDesignTemperature+1)
       if (nDesignPointRefX >= 0) &
            pointrefd(1) = extraDot(nDesignPointRefX+1)
       if (nDesignPointRefY >= 0) &
            pointrefd(2) = extraDot(nDesignPointRefY+1)
       if (nDesignPointRefZ >= 0) &
            pointrefd(3) = extraDot(nDesignPointRefZ+1)
    end if

    ! Now set any STATE seeds
    if (useState) then 

       ! Now we have to set all the seeds
       ii = 0
       domainLoop2: do nn=1,nDom

          ! Just to get sizes
          call setPointers_d(nn, level, 1)

          spectalLoop2: do sps=1,nTimeIntervalsSpectral
             do k=2, kl
                do j=2,jl
                   do i=2,il
                      do l = 1, nw
                         ii = ii + 1
                         flowdomsd(nn,level,sps)%w(i, j, k, l) = wDot(ii)
                      end do
                   end do
                end do
             end do
          end do spectalLoop2
       end do domainLoop2

       ! Now run the derivatie halo exchange:
       call whalo2_d(level, 1, nw, .False., .False., .False.)
    end if

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

    funcsLocalDot = zero

    ! Now we are ready to call block_res_d with all the correct seeds
    ii = 0
    jj = 0
    domainLoopAD: do nn=1,nDom

       ! Set pointers to the first timeInstance...just to getSizes
       call setPointers(nn, level, 1)
       ! Set unknown sizes in diffSizes for AD routine
       ISIZE1OFDrfbcdata = nBocos
       ISIZE1OFDrfviscsubface = nViscBocos

       spectalLoopAD: do sps=1,nTimeIntervalsSpectral
          ! Set pointers and derivative pointers
          call setPointers_d(nn, level, sps)

          ! Get the pointers from the petsc vector for the surface
          ! perturbation for wall distance. 
          call VecGetArrayF90(xSurfVec(level, sps), xSurf, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! And it's derivative
          call VecGetArrayF90(xSurfVecd(sps), xSurfd, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          call BLOCK_RES_D(nn, level, useSpatial, frozenTurbulence)

          ! Now extract dw
          do sps2=1,nTimeIntervalsSpectral
             do k=2, kl
                do j=2, jl
                   do i=2, il
                      do l=1, nw
                         ii = ii + 1
                         dwDot(ii) = flowdomsd(nn, level, sps2)%dw(i, j, k, l)
                      end do
                   end do
                end do
             end do
          end do

          ! We need to SUM the funcs into the local array
          funcsLocalDot = funcsLocalDot + funcValuesd

          ! These arrays need to be restored before we can move to the next spectral instance. 
          call VecRestoreArrayF90(xSurfVec(level, sps), xSurf, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! And it's derivative
          call VecRestoreArrayF90(xSurfVecd(sps), xSurfd, ierr)
          call EChk(ierr,__FILE__,__LINE__)

       end do spectalLoopAD
    end do domainLoopAD

    ! This only works currently with 1 spectral instance
    sps = 1
    call getForces_d(forces, fDot(:, :, sps), fSize, sps)

    ! We need to allreduce the function values to all processors
    call mpi_allreduce(funcsLocalDot, funcsDot, costSize, sumb_real, mpi_sum, SUmb_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Reset the correct equation parameters if we were useing the frozen
    ! Turbulent 
    if (resetToRANS) then
       equations = RANSEquations
    end if

    call VecResetArray(x_like, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Just to be sure, we'll zero everything when we're done.
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn, level, sps)
          call zeroADSeeds(nn,level, sps)
       end do
    end do
  end subroutine computeMatrixFreeProductFwd

  subroutine computeMatrixFreeProductBwd(dwbar, funcsbar, fbar, useSpatial, useState, xvbar, &
       extrabar, wbar, spatialSize, extraSize, stateSize, costSize, fSize, nTime)
    use constants
    use block, only : flowDomsd
    use communication, only : sumb_comm_world
    use blockPointers
    use inputDiscretization 
    use inputTimeSpectral 
    use inputPhysics
    use iteration         
    use flowVarRefState     
    use inputAdjoint       
    use diffSizes
    use ADjointPETSc, only : x_like, psi_like3
    use adjointvars
    use costfunctions
    use wallDistanceData, only : xSurfVec, xSurfVecd, xSurf, xSurfd, wallScatter
    use utils, only : setPointers, EChk, getDirAngle, setPointers_d
    use solverUtils, only : timeStep_block
    use flowUtils, only : allNodalGradients
    use fluxes, only : viscousFlux
    use haloExchange, only : exchangeCoor_b
    use adjointextra_b, only : block_res_b, xhalo_block_b
    use adjointUtils, only : allocDerivativeValues, zeroADSeeds, setDiffSizes
    implicit none

#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"

    ! Input Variables
    integer(kind=intType), intent(in) :: spatialSize, extraSize, stateSize, costSize, fSize, nTime
    real(kind=realType), dimension(stateSize), intent(in) :: dwbar
    real(kind=realType), dimension(costSize), intent(in) :: funcsbar
    real(kind=realType), dimension(3, fSize, nTime), intent(in) :: fBar
    logical, intent(in) :: useSpatial, useState

    ! Ouput Variables
    real(kind=realType), dimension(stateSize), intent(out) :: wbar
    real(kind=realType), dimension(extraSize), intent(out) :: extrabar
    real(kind=realType), dimension(spatialSize), intent(out) :: xvbar

    ! Working variables
    integer(kind=intType) :: ierr,nn,mm,sps,i,j,k,l,ii,jj,sps2, idim
    integer(kind=intType) ::  level, irow, nState
    logical :: resetToRans
    real(kind=realType), dimension(extraSize) :: extraLocalBar
    real(kind=realType), dimension(:, :), allocatable :: xSurfbSum

    ! Setup number of state variable based on turbulence assumption
    if ( frozenTurbulence ) then
       nState = nwf
    else
       nState = nw
    endif

    ! Place output arrays in psi_like and x_like vectors if necessary
    wbar = zero
    if (useState) then 
       call VecPlaceArray(psi_like3, wbar, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end if

    ! Place adjoint in Vector
    xvbar = zero
    if (useSpatial) then 
       call VecPlaceArray(x_like, xvbar, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end if

    ! Need to trick the residual evalution to use coupled (mean flow and
    ! turbulent) together.
    level = 1
    currentLevel = level
    groundLevel = level

    ! Determine if we want to use frozenTurbulent Adjoint
    resetToRANS = .False. 
    if (frozenTurbulence .and. equations == RANSEquations) then
       equations = NSEquations 
       resetToRANS = .True.
    end if

    ! Allocate the memory for reverse
    if (.not. derivVarsAllocated) then 
       call allocDerivativeValues(level)
    end if
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn, level, sps)
          call zeroADSeeds(nn,level, sps)
       end do
    end do

    ! Zero the function seeds
    funcValuesd= zero

    call VecGetLocalSize(xSurfVec(1, 1), i, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    allocate(xSurfbSum(nTimeIntervalsSpectral, i))
    xSurfbSum = zero

    ! Zero out extraLocal
    extraLocalBar = zero

    ii = 0 ! Residual bar counter
    jj = 0 ! Force bar counter

    ! Before we start, perform the reverse of getForces() for each
    ! spectral instance. This set the bcDatad%F and bcData%area
    ! seeds.
    do sps=1, nTimeIntervalsSpectral
       call getForces_b(fBar(:, :, sps), fSize, sps)
    end do

    domainLoopAD: do nn=1,nDom

       ! Just to get sizes
       call setPointers(nn,1_intType,1)
       call setDiffSizes

       do sps=1,nTimeIntervalsSpectral
          ! Set pointers and derivative pointers
          call setPointers_d(nn, level, sps)

          ! Get the pointers from the petsc vectors
          call VecGetArrayF90(xSurfVec(level, sps), xSurf, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! And it's derivative
          call VecGetArrayF90(xSurfVecd(sps), xSurfd, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Set the dw seeds
          do k=2, kl
             do j=2,jl
                do i=2,il
                   do l=1,nState
                      ii = ii + 1
                      flowdomsd(nn, level, sps)%dw(i, j, k, l) = dwbar(ii)
                   end do
                end do
             end do
          end do

          ! And the function value seeds
          funcValuesd = funcsBar

          ! For some reason tapenade forget to zero rgasd when it
          ! should, therefore we must manually zero here before calling the code. 
          rgasd = zero
          call BLOCK_RES_B(nn, sps, useSpatial, frozenTurbulence)

          ! Currently these tapenade has decided the output values from
          ! these routines do not matter, these need to be recomputed to
          ! consistent.

          call timestep_block(.false.)
          if (viscous) then 
             call allNodalGradients
             call viscousflux
          end if

          ! Assmeble the vectors requested:
          do sps2=1,nTimeIntervalsSpectral

             if (useSpatial) then 
                ! We need to acculumate the contribution from this block into xSurfbSum
                xSurfbSum(sps, :) = xSurfbSum(sps, :) + xSurfd

                ! Also need the extra variables, those are zero-based:
                if (nDesignAoA >= 0) &
                     extraLocalBar(nDesignAoA+1) = extraLocalBar(nDesignAoA+1) + alphad
                if (nDesignSSA >= 0) &
                     extraLocalBar(nDesignSSA+1) = extraLocalBar(nDesignSSA+1) + alphad
                if (nDesignMach  >= 0) &
                     extraLocalBar(nDesignMach+1) = extraLocalBar(nDesignMach+1) + machd + machcoefd
                if (nDesignMachGrid >= 0) &
                     extraLocalBar(nDesignMachGrid+1) = extraLocalBar(nDesignMachGrid+1) + machgridd + machcoefd
                if (nDesignPressure >= 0) &
                     extraLocalBar(nDesignPressure+1) = extraLocalBar(nDesignPressure+1) + pinfdimd
                if (nDesignTemperature >= 0) &
                     extraLocalBar(nDesignTemperature+1) = extraLocalBar(nDesignTemperature+1) + tinfdimd
                if (nDesignDensity >= 0) &
                     extraLocalBar(nDesignDensity+1) = extraLocalBar(nDesignDensity+1) + rhoinfdimd
                if (nDesignPointRefX >= 0) &
                     extraLocalBar(nDesignPointRefX+1) = extraLocalBar(nDesignPointRefX+1) + pointrefd(1)
                if (nDesignPointRefY >= 0) &
                     extraLocalBar(nDesignPointRefY+1) = extraLocalBar(nDesignPointRefY+1) + pointrefd(2)
                if (nDesignPointRefZ >= 0) &
                     extraLocalBar(nDesignPointRefZ+1) = extraLocalBar(nDesignPointRefZ+1) + pointrefd(3)
             end if

             if (useState) then 
                do k=0, kb
                   do j=0,jb
                      do i=0,ib
                         do l=1,nState
                            irow = flowDoms(nn, 1, sps2)%globalCell(i,j,k)*nState + l -1
                            if (irow >= 0) then 
                               call VecSetValues(psi_like3, 1, (/irow/), &
                                    (/flowdomsd(nn, level, sps)%w(i, j, k, l)/), ADD_VALUES, ierr)
                               call EChk(ierr,__FILE__,__LINE__)
                            end if
                         end do
                      end do
                   end do
                end do
             end if

          end do

          ! These arrays need to be restored before we can move to the next spectral instance. 
          call VecRestoreArrayF90(xSurfVec(level, sps), xSurf, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! And it's derivative
          call VecRestoreArrayF90(xSurfVecd(sps), xSurfd, ierr)
          call EChk(ierr,__FILE__,__LINE__)

       end do
    end do domainLoopAD

    if (useSpatial) then

       ! We are not done with the spatial code just yet, we need to run
       ! the backwards exchange coor and then run xhalo_block_b

       call exchangeCoor_b(level)

       do nn=1,nDom
          do sps=1,nTimeIntervalsSpectral
             call setPointers_d(nn, level, sps)

             ! Complete the xhalo dependance
             call xhalo_block_b

             ! We can now just need to loop over the owned nodes...no halos now
             do k=1,kl
                do j=1,jl
                   do i=1,il

                      do l=1,3
                         irow = flowDoms(nn, 1, sps)%globalNode(i,j,k)*3 + l -1
                          call VecSetValues(x_like, 1, (/irow/), &
                               (/flowdomsd(nn, level, sps)%x(i, j, k, l)/), INSERT_VALUES, ierr)
                         call EChk(ierr,__FILE__,__LINE__)
                      end do
                   end do
                end do
             end do
          end do
       end do

       ! Finally we have to do an mpi all reduce on the local parts:

       extraBar = zero
       call mpi_allreduce(extraLocalBar, extraBar, extraSize, sumb_real, mpi_sum, SUmb_comm_world, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end if

    ! Copy the Sum value back to xSurfb which is actually the PETSc
    ! array. And now we are done with the xSurfBSum value

    do sps=1,nTimeIntervalsSpectral

       ! And it's derivative
       call VecGetArrayF90(xSurfVecd(sps), xSurfd, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       xSurfd = xSurfbSum(sps, :)
       call VecRestoreArrayF90(xSurfVecd(sps), xSurfd, ierr)
       call EChk(ierr,__FILE__,__LINE__)

    end do

    deallocate(xsurfbSum)

    ! And perform assembly on the w vectors if 
    if (useState) then 
       call VecAssemblyBegin(psi_like3, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecAssemblyEnd(psi_like3, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecResetArray(psi_like3, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end if

    if (useSpatial) then 
       ! And perform assembly on the x vectors
       call VecAssemblyBegin(x_like, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecAssemblyEnd(x_like, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Now scatter the xsurb contribution back into x_like Perform the
       ! scatter from the global x vector to xSurf...but only if
       ! wallDistances were used
       if (wallDistanceNeeded .and. useApproxWallDistance) then 
          do sps=1, nTimeIntervalsSpectral
             call VecScatterBegin(wallScatter(1, sps), xSurfVecd(sps), x_like, ADD_VALUES, SCATTER_REVERSE, ierr)
             call EChk(ierr,__FILE__,__LINE__)

             call VecScatterEnd(wallScatter(1, sps), xSurfVecd(sps), x_like, ADD_VALUES, SCATTER_REVERSE, ierr)
             call EChk(ierr,__FILE__,__LINE__)
          end do
       end if
       call VecResetArray(x_like, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end if

    ! Reset the correct equation parameters if we were useing the frozen
    ! Turbulent 
    if (resetToRANS) then
       equations = RANSEquations
    end if

    ! Just to be sure, we'll zero everything when we're done.
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn, level, sps)
          call zeroADSeeds(nn,level, sps)
       end do
    end do
  end subroutine computeMatrixFreeProductBwd

  subroutine computeMatrixFreeProductBwdFast(dwbar, wbar, stateSize)
    ! This is the "Fast" ie. State variable only version of the reverse
    ! mode computation. It is intended to compute dRdw^T product
    ! ONLY. The main purpose is for fast matrix-vector products for the
    ! actual adjoint solve. 
    use block, only : flowDomsd
    use blockPointers
    use inputDiscretization
    use inputTimeSpectral 
    use inputPhysics
    use flowVarRefState
    use inputAdjoint       
    use iteration
    use inputIteration
    use sa_fast_b, only : saresscale_fast_b, saviscous_fast_b, &
         sasource_fast_b, cb3Inv, cv13, cw36, kar2inv, qq
    use adjointvars
    use communication, only : sumb_comm_world
    use paramTurb
    use utils, only : terminate, setPointers_d
    use haloExchange, only : whalo2_b
    use flowutils_fast_b
    use turbutils_fast_b
    use turbbcroutines_b
    use turbutils_b, only : saeddyviscosity_b
    use fluxes_fast_b
    use solverutils_fast_b
    use flowutils_fast_b, only : computelamviscosity_fast_b, allnodalgradients_fast_b
    implicit none

    ! Input Variables
    integer(kind=intType), intent(in) :: stateSize
    real(kind=realType), dimension(stateSize), intent(in) :: dwbar

    ! Ouput Variables
    real(kind=realType), dimension(stateSize), intent(out) :: wbar

    ! Working variables
    integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
    integer(kind=intType) :: nState, level
    logical :: resetToRans
    real(kind=realType) :: ovol, timea, timeb, totaltime

    ! Setup number of state variable based on turbulence assumption
    if ( frozenTurbulence ) then
       nState = nwf
    else
       nState = nw
    endif
    ! Assembling matrix on coarser levels is not entirely implemented yet. 
    level = 1
    currentLevel = level
    groundLevel = level

    ! Determine if we want to use frozenTurbulent Adjoint
    resetToRANS = .False. 
    if (frozenTurbulence .and. equations == RANSEquations) then
       equations = NSEquations 
       resetToRANS = .True.
    end if

    ! Note: The calling routine is responsible for ensuring that the
    ! derivative values are allocated AND ZEROED! This routine makes use
    ! of the fact that only wbar needs to be zeroed since all other
    ! required seeds are zeroed in the individual fast routines. This is
    ! slightly unsafe, but it necessary for speed. 
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          flowDomsd(nn, level, sps)%w = zero
       end do
    end do

    ii = 0

    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          ! Set pointers and derivative pointers
          call setPointers_d(nn, level, sps)
          do k=2,kl
             do j=2,jl
                do i=2,il
                   ovol = one/vol(i,j,k)
                   do l=1,nwf
                      dwd(i,j,k,l) = dwbar(ii+ l)*ovol
                      fwd(i,j,k,l) = dwd(i,j,k,l)
                   end do
                   do l=nt1,nState
                      dwd(i,j,k,l) = dwbar(ii+ l)*ovol*turbresscale(l-nt1+1)
                      fwd(i,j,k,l) = dwd(i,j,k,l)
                   end do
                   ii = ii + nState
                end do
             end do
          end do

          call pushreal8array(radk, size(radk, 1)*size(radk, 2)*size(radk, 3))
          call pushreal8array(radj, size(radj, 1)*size(radj, 2)*size(radj, 3))
          call pushreal8array(radi, size(radi, 1)*size(radi, 2)*size(radi, 3))

          if (viscous) then 
             call viscousFlux_fast_b
             call allnodalgradients_fast_b
             call computespeedofsoundsquared_fast_b
          end if

          select case (spaceDiscr)
          case(dissScalar) 
             call inviscidDissFluxScalar_fast_b
          case(dissMatrix)
             call inviscidDissFluxMatrix_fast_b
          end select

          call inviscidcentralflux_fast_b

          if (equations == RANSEquations) then 
             select case(turbModel)
             case (spalartAllmaras)
                cv13    = rsaCv1**3
                kar2Inv = one/(rsaK**2)
                cw36    = rsaCw3**6
                cb3Inv  = one/rsaCb3
                call saresscale_fast_b()
                call saviscous_fast_b()
                call turbadvection_fast_b(1_inttype, 1_inttype, itu1-1, qq)
                call sasource_fast_b()
             case default
                call terminate("matrixFreeRoutines", &
                     "Only SA turbulence adjoint implemented")
             end select

             ! Do the production term
             select case  (turbprod) 
             case (strain) 
                call prodsmag2_fast_b()
             case (vorticity) 
                call prodwmag2_fast_b()
             case (katolaunder) 
                call prodkatolaunder_fast_b()
             end select

             ! And the turbulence BCs
             call applyallturbbcthisblock_b(.true.)
             call bcturbtreatment_b()
          end if

          call timestep_block_fast_b(.False.)
          call applyAllBC_block_fast_b(.True.)

          call popreal8array(radi, size(radi, 1)*size(radi, 2)*size(radi, 3))
          call popreal8array(radj, size(radj, 1)*size(radj, 2)*size(radj, 3))
          call popreal8array(radk, size(radk, 1)*size(radk, 2)*size(radk, 3))

       end do
    end do

    ! Communicate all the derivative values in reverse
    call whalo2_b(1, 1, nw, .True., .True., .True.)

    ii = 0
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers_d(nn, level, sps)

          if (equations == RANSEquations) then 
             select case(turbModel)
             case (spalartAllmaras)
                call saeddyviscosity_b(0, ib, 0, jb, 0, kb)
             end select
          end if

          if (viscous) then 
             call computelamviscosity_fast_b(.True.)
          end if

          call computepressuresimple_fast_b(.true.)

          ! We can put stuff directly into wbar with no assembly; the
          ! whalo_b already takes care of it. 
          do k=2, kl
             do j=2,jl
                do i=2,il
                   do l=1,nState
                      ii =ii + 1
                      wbar(ii) = flowdomsd(nn, level, sps)%w(i,j,k,l)
                   end do
                end do
             end do
          end do
       end do
    end do

    ! Reset the correct equation parameters if we are using the frozen
    ! Turbulent
    if (resetToRANS) then
       equations = RANSEquations
    end if
  end subroutine computeMatrixFreeProductBwdFast
#endif 
  ! if def for complex

  subroutine solveAdjointForRHS(inVec, outVec, nDOF, relativeTolerance)

    use ADJointPETSc
    use inputADjoint
    use adjointvars
    use killsignals
    use constants
    use blockPointers
    use inputTimeSpectral
    use utils, only : EChk
    use adjointUtils, only : allocDerivativeValues, zeroADSeeds
    implicit none

    ! Input Variables
    real(kind=realType), dimension(ndof), intent(in) :: inVec
    real(kind=realType), dimension(ndof), intent(out) :: outVec
    real(kind=realType), intent(in) :: relativeTolerance
    integer(kind=intType), intent(in) :: nDOF
    integer(kind=intTYpe) :: adjointConvergedReason
    ! Working variables
    integer(kind=intType) :: ierr, nn, sps
    real(kind=realType) :: val

#ifndef USE_COMPLEX

    ! Place the arrays
    call VecPlaceArray(psi_like1, inVec, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecPlaceArray(psi_like2, outVec, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Zero out initial solution
    call VecSet(psi_like2, zero, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Set desired realtive tolerance
    call KSPSetTolerances(adjointKSP, relativeTolerance, adjAbsTol, adjDivTol, &
         adjMaxIter, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Make sure the derivative memory is allocated and zeroed. 
    if (.not. derivVarsAllocated) then 
       call allocDerivativeValues(1_intType)
    end if

    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call zeroADSeeds(nn, 1_intType, sps)
       end do
    end do

    ! Solve (remember this is actually a transpose solve)
    call KSPSolve(adjointKSP, psi_like1, psi_like2, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call KSPGetConvergedReason(adjointKSP, adjointConvergedReason, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (adjointConvergedReason ==  KSP_CONVERGED_RTOL .or. &
         adjointConvergedReason ==  KSP_CONVERGED_ATOL .or. &
         adjointConvergedReason ==  KSP_CONVERGED_HAPPY_BREAKDOWN) then
       adjointFailed = .False.
    else
       adjointFailed = .True.
    end if

    ! Rest arrays
    call VecResetArray(psi_like1,  ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecResetArray(psi_like2,  ierr)
    call EChk(ierr,__FILE__,__LINE__)

#endif

  end subroutine solveAdjointForRHS

  subroutine solveDirectForRHS(inVec, outVec, nDOF, relativeTolerance)

    use ADJointPETSc
    use inputADjoint
    use adjointVars
    use constants
    use killsignals
    use blockPointers
    use inputTimeSpectral
    use utils, only : EChk
    use adjointUtils, only : allocDerivativeValues, zeroADSeeds
    implicit none

    ! Input Variables
    real(kind=realType), dimension(ndof), intent(in) :: inVec
    real(kind=realType), dimension(ndof), intent(out) :: outVec
    real(kind=realType), intent(in) :: relativeTolerance
    integer(kind=intType), intent(in) :: nDOF
    integer(kind=intTYpe) :: adjointConvergedReason
    ! Working variables
    integer(kind=intType) :: ierr, nn, sps

#ifndef USE_COMPLEX

    ! Place the arrays
    call VecPlaceArray(psi_like1, inVec, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecPlaceArray(psi_like2, outVec, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Zero out initial solution
    call VecSet(psi_like2, zero, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Set desired realtive tolerance
    call KSPSetTolerances(adjointKSP, relativeTolerance, adjAbsTol, adjDivTol, &
         adjMaxIter, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Make sure the derivative memory is allocated and zeroed. 
    if (.not. derivVarsAllocated) then 
       call allocDerivativeValues(1_intType)
    end if

    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call zeroADSeeds(nn, 1_intType, sps)
       end do
    end do

    ! Solve (this is the transpose solve of a transpose matrix, so it's direct)
    call KSPSolveTranspose(adjointKSP, psi_like1, psi_like2, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call KSPGetConvergedReason(adjointKSP, adjointConvergedReason, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    if (adjointConvergedReason ==  KSP_CONVERGED_RTOL .or. &
         adjointConvergedReason ==  KSP_CONVERGED_ATOL .or. &
         adjointConvergedReason ==  KSP_CONVERGED_HAPPY_BREAKDOWN) then
       adjointFailed = .False.
    else
       adjointFailed = .True.
    end if

    ! Rest arrays
    call VecResetArray(psi_like1, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecResetArray(psi_like2, ierr)
    call EChk(ierr,__FILE__,__LINE__)

#endif

  end subroutine solveDirectForRHS

  subroutine saveADjointMatrix(fileName)

    use constants
    use ADjointPETSc, only: drdwt
    use communication, only : sumb_comm_world
    use utils, only : EChk
    implicit none

#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"

    ! Input params
    character*(*), intent(in) :: fileName

    ! Working parameters
    PetscViewer binViewer
    integer(kind=intType) :: ierr

    call PetscViewerBinaryOpen(sumb_comm_world, fileName, FILE_MODE_WRITE, binViewer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call MatView(dRdwT, binViewer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call PetscViewerDestroy(binViewer,ierr)
    call EChk(ierr, __FILE__, __LINE__)

  end subroutine saveADjointMatrix

  subroutine saveAdjointPC(fileName)

    use constants
    use ADjointPETSc, only: drdwpret
    use communication, only : sumb_comm_world
    use utils, only : EChk
    implicit none

    PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"

    ! Input params
    character*(*), intent(in) :: fileName

    ! Working parameters
    PetscViewer binViewer
    integer(kind=intType) :: ierr

    call PetscViewerBinaryOpen(sumb_comm_world, fileName, FILE_MODE_WRITE, binViewer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call MatView(dRdwPreT, binViewer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call PetscViewerDestroy(binViewer,ierr)
    call EChk(ierr, __FILE__, __LINE__)

  end subroutine saveAdjointPC

  

  subroutine spectralPrecscribedMotion(input, nin, dXv, nout)

    use constants
    use blockPointers, only : il, jl, kl, nDom
    use section, only : sections, nSections
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use monitor , only : timeUnsteadyRestart, timeUnsteady
    use utils, only : setPointers, rotMatrixRigidBody
    implicit none
    ! Input/Output Variables
    integer(kind=intType), intent(in) :: nin, nout
    real(kind=realType), intent(out)  :: dXv(nout)
    real(kind=realType), intent(in)   :: input(nin)

    ! Local Variables
    integer(kind=intType) :: ierr, sps, i, nn, mm, counter0, counter1
    integer(kind=intType) :: nodes_on_block, cum_nodes_on_block
    real(kind=realType), dimension(3)   :: rotationPoint, r
    real(kind=realType), dimension(3, 3) :: rotationMatrix  
    real(kind=realType) :: t(nSections), dt(nSections)
    real(kind=realType) :: tOld, tNew, pt(3)
    real(kind=realType), pointer :: xvec_pointer(:)
    real(kind=realType) :: time(3)

    !       For the TimeSpectral case, we need to include    *
    !      the operation that rotates the base grid to each time instance 
    !      This is basically the reverse of the operation that is done in 
    !      setGrid.f90                                                    
    !      The operation in setGrid.f90 is the following                  
    !      X_sps = M(X - rotPoint) + rotPoint                             
    !      where                                                          
    !      X_sps is the set of coordinates at each time instance          
    !      M is the rotation matrix calculated by rotMatrixRigidBody      
    !      rotPoint is the point about which the motion takes place       
    !      It is easy to see dX_sps/dX = M                                
    !      What we are actually computing is the following:               
    !                 T          T                                        
    !        /dX_sps \ /   dR   \                                         
    !        |-------| |------- |  psi                                    
    !        \  dX   / \ dX_sps /                                         

    ! Zero dXv for time spectral case since we add to array.
    dXv = zero

    do nn=1, nSections
       dt(nn) = sections(nn)%timePeriod &
            / real(nTimeIntervalsSpectral, realType)
    enddo

    timeUnsteady = zero
    counter0 = 0
    cum_nodes_on_block = 0
    ! The nDom loop followed by the sps loop is required to follow
    ! the globalNode ordering such that we can use the pointer from
    ! vecGetArrayF90

    do nn=1, nDom
       do sps = 1, nTimeIntervalsSpectral

          call setPointers(nn, 1, sps)
          nodes_on_block = il*jl*kl

          do mm=1, nSections
             t(mm) = (sps-1)*dt(mm)
          enddo

          ! Compute the displacements due to the rigid motion of the mesh.

          tNew = timeUnsteady + timeUnsteadyRestart
          tOld = tNew - t(1)

          call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

          ! Take rotation Matrix Transpose
          rotationMatrix = transpose(rotationMatrix)

          counter1 = cum_nodes_on_block        

          ! Loop over the localally owned nodes:
          do i=1, nodes_on_block
             pt = (/input(3*counter0+1), &
                  input(3*counter0+2), &
                  input(3*counter0+3)/)

             dXv(3*counter1+1:3*counter1+3) = &
                  dXv(3*counter1+1:3*counter1+3) + &
                  matmul(rotationMatrix, pt)

             counter0 = counter0 + 1
             counter1 = counter1 + 1
          end do

       end do
       ! Increment the cumulative number of nodes by the nodes on the
       ! block we just did
       cum_nodes_on_block = cum_nodes_on_block + nodes_on_block
    end do

  end subroutine spectralPrecscribedMotion

  subroutine setupAllResidualMatricesfwd

    use constants
    use ADjointPETSc, only : dRdwT
    use communication, only : sumb_comm_world, myid
    use inputADjoint, only : frozenTurbulence, useMatrixFreedRdw
    use adjointUtils, only : setupStateResidualMatrix
    use utils, only : EChk
    implicit none

    logical :: useAD, useTranspose, usePC, useObjective
    real(kind=realType) :: timeAdjLocal, timeAdj, time(2)
    integer(kind=intType) :: ierr

    ! If we are assembling matrices...we ned to assemble the
    ! 'transpose', with 'AD', we want the exact matrix not the 'PC',
    ! and will compute objective RHS
    useAD = .True.
    usePC = .False.
    useTranspose = .True.
    useObjective = .True.

    if (.not. useMatrixFreedRdw) then 
       if( myid ==0 ) then
          write(*, 10) "Assembling State Residual Matrix in Forward mode..."
       end if
       time(1) = mpi_wtime()
       call setupStateResidualMatrix(drdwT, useAD, usePC, useTranspose, &
            useObjective, frozenTurbulence, 1_intType)
       time(2) = mpi_wtime()
       timeAdjLocal = time(2)-time(1)

       call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
            mpi_max, 0, SUMB_COMM_WORLD, ierr)
       call EChk(ierr,  __FILE__, __LINE__)

       if(myid ==0)  then 
          write(*, 20) "Assembling State Residaul Matrices Fwd time (s) = ", timeAdj
       end if
    end if

    ! Output formats.
10  format(a)
20  format(a, 1x, f8.2)

  end subroutine setupAllResidualMatricesfwd

  subroutine solveAdjoint(RHS, psi, checkSolution, nState)
    !
    !      Solve the linear discrete ADjoint system of equations          
    !          [dR/dW]T . psi = {RHS}                                     
    !      using preconditioned GMRES provided by PETSc. The values in psi
    !      are significant as they are used as the inital guess.          
    !
    use constants
    use ADjointPETSc, only : dRdwT, psi_like1, psi_like2, adjointKSP, &
         adjResInit, adjResStart, adjResFinal

    use killsignals, only : adjointFailed
    use inputADjoint, only : adjAbsTol, adjDivTol, adjMaxIter, adjRelTol, &
         adjRelTolRel, printTiming
    use adjointVars, only: derivVarsAllocated
    use communication, only : myid, sumb_comm_world
    use blockPointers, only : nDom
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use adjointUtils, only : allocDerivativeValues, zeroADSeeds
    use utils, only : EChk
    implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"

    ! Input Parameters
    real(kind=realType), dimension(nState) :: RHS, psi
    integer(kind=intType) :: nState
    logical :: checkSolution 
    !
    !     Local variables.
    real(kind=alwaysRealType)   :: norm
    real(kind=alwaysRealType), dimension(2) :: time
    real(kind=alwaysRealType)               :: timeAdjLocal, timeAdj
    real(kind=realType) :: l2abs, l2rel
    integer(kind=intType) :: ierr, nn, sps
    integer(kind=intType) :: adjConvIts
    KSPConvergedReason adjointConvergedReason
    Vec adjointRes, RHSVec

    ! Send some feedback to screen.

    if(myid ==0 .and. printTiming)  &
         write(*,10) "Solving ADjoint Transpose with PETSc..."

    call cpu_time(time(1))

    ! Make sure the derivative memory is allocated and zeroed. 
    if (.not. derivVarsAllocated) then 
       call allocDerivativeValues(1_intType)
    end if

    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call zeroADSeeds(nn, 1_intType, sps)
       end do
    end do

    ! Dump psi into psi_like1 and RHS into psi_like2
    call VecPlaceArray(psi_like1, psi, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecPlaceArray(psi_like2, RHS, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecDuplicate(psi_like1, adjointRes, ierr)
    call EChk(ierr,__FILE__,__LINE__)
    if (checkSolution) then 
       call VecDuplicate(psi_like1, RHSVec, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call vecCopy(psi_like2, RHSVec, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end if

    ! Get the RHS norm....this is the 'init' norm:
    call VecNorm(psi_like2, NORM_2, adjResInit, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Get Current Residual -- we always solve for the delta
    call MatMult(dRdWT, psi_like1, adjointRes, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! AdjointRes = AdjointRes - adjointRHS
    call VecAXPY(adjointRes, -one, psi_like2, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Norm of adjoint Residual
    call VecNorm(adjointRes, NORM_2, adjResStart,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! The way we use tolerances are as follows: The residual must
    ! statify:
    ! res < adjRelTol * adjResInit OR 
    ! res < adjRelTolRel * adjResStart OR
    ! res < adjAbsTol

    ! L2Abs is used to stipulate an exit criteria for adjreltolrel
    L2abs = adjResStart * adjreltolrel

    ! If L2Abs is less that what we actually want as the absolute
    ! tolerance, clip it
    if (L2Abs < adjAbsTol) then
       L2abs = adjabstol
    end if

    ! L2Rel is a little tricky since if the start residual is *larger*
    ! than the init residual, it won't converge enough. While this seems
    ! strange this is *always* the case for restarted RANS-based
    ! adjoints.
    L2Rel = (adjReltol * adjResInit) / adjResStart

    ! We need to clip L2Rel such that it can never be greater than one. 
    L2Rel = min(L2Rel, 0.9)

    ! Set the tolerances
    call KSPSetTolerances(adjointKSP, L2Rel, L2Abs, adjDivTol, &
         adjMaxIter, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Solve the update (psi_like2)
    call KSPSolve(adjointKSP, adjointRes, psi_like2, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Now compute the update to psi_like1 (psi)
    call VecAXPY(psi_like1, -one, psi_like2, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    if (checkSolution) then 

       ! Get new time and compute the elapsed time.
       call cpu_time(time(2))
       timeAdjLocal = time(2)-time(1)

       ! Determine the maximum time using MPI reduce
       ! with operation mpi_max.

       call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
            mpi_max, 0, SUMB_COMM_WORLD, ierr)

       call MatMult(dRdWT, psi_like1, adjointRes, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecAXPY(adjointRes, -one, RHSVec, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecNorm(adjointRes, NORM_2, norm,ierr)
       call EChk(ierr,__FILE__,__LINE__)
       adjResFinal = norm

       call KSPGetIterationNumber(adjointKSP,adjConvIts,ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Use the root processor to display the output summary, such as
       ! the norm of error and the number of iterations

       if( myid ==0 .and. printTiming) then
          write(*,20) "Solving ADjoint Transpose with PETSc time (s) =", timeAdj
          write(*,30) "Norm of error =",norm,"Iterations =",adjConvIts
          write(*,*) "------------------------------------------------"
          if( adjConvIts.lt.0 ) then
             write(*,40) "PETSc solver diverged after", -adjConvIts, &
                  "iterations..."
          else
             write(*,40) "PETSc solver converged after", adjConvIts, &
                  "iterations."
          endif
          write(*,*) "------------------------------------------------"
       endif

       call VecDestroy(RHSVec, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end if

    ! Destroy the temporary vector and reset the arrays
    call VecDestroy(adjointRes, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecResetArray(psi_like1, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecResetArray(psi_like2, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Get the petsc converged reason and set the fail flag
    call KSPGetConvergedReason(adjointKSP, adjointConvergedReason,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    if (adjointConvergedReason ==  KSP_CONVERGED_RTOL .or. &
         adjointConvergedReason ==  KSP_CONVERGED_ATOL .or. &
         adjointConvergedReason ==  KSP_CONVERGED_HAPPY_BREAKDOWN) then
       adjointFailed = .False.
    else
       adjointFailed = .True.
    end if

    ! Output formats.

10  format(a)
20  format(a,1x,f8.2)
30  format(1x,a,1x,e10.4,4x,a,1x,i4)
40  format(1x,a,1x,i5,1x,a)

  end subroutine solveAdjoint

  subroutine setupPETScKsp

    use ADjointPETSc, only: drdwpret, drdwt, adjointKSP
    use inputADjoint
    use utils, only : ECHk, terminate
    use adjointUtils, only : mykspmonitor
    use adjointUtils, only : setupStateResidualMatrix, setupStandardKSP

    implicit none

#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"

    !     Local variables.
    logical :: useAD, usePC, useTranspose, useObjective
    integer(kind=intType) :: ierr

    if (ApproxPC)then
       !setup the approximate PC Matrix
       useAD = ADPC
       useTranspose = .True.
       usePC = .True.
       useObjective = .False.
       call setupStateResidualMatrix(drdwpret, useAD, usePC, useTranspose, &
            useObjective, frozenTurbulence, 1_intType)
       call KSPSetOperators(adjointKSP, dRdwT, dRdWPreT, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    else
       ! Use the exact jacobian.  Here the matrix that defines the
       ! linear system also serves as the preconditioning matrix. This
       ! is only valid if useMatrixFree is flase. 
       if (useMatrixfreedRdw) then 
          call terminate("setupPETScKSP", "useMatrixFreedRdW option cannot be true when the approxPC option is False")
       end if
       call KSPSetOperators(adjointKSP, dRdWt, dRdWT, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end if

    if (PreCondType == 'asm') then
       ! Run the super-dee-duper function to setup the ksp object:
       call setupStandardKSP(adjointKSP, ADjointSolverType, adjRestart, adjointpcside, &
            PreCondType, overlap, outerPreConIts, localPCType, &
            matrixOrdering, FillLevel, innerPreConIts)
    else if (PreCondType == 'mg') then
       print *,'Only ASM precondtype is usable'
       stop
    end if

    ! Setup monitor if necessary:
    if (setMonitor) then
       call KSPMonitorSet(adjointKSP, MyKSPMonitor, PETSC_NULL_OBJECT, &
            PETSC_NULL_FUNCTION, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    endif

  end subroutine setupPETScKsp

  subroutine saveCellCenters(fileName)

    use blockPointers
    use iteration
    use inputTimeSpectral
    use adjointVars, only: nCellsLocal
    use communication
    use utils, only : setPointers, EChk
    implicit none

#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"

    ! Input params
    character*(*), intent(in) :: fileName

    ! Working parameters
    PetscViewer binViewer
    integer(kind=intType) :: nn,sps,i,j,k,n
    integer(kind=intType) :: ierr, iRow, level
    real(kind=realType),dimension(3)::cellCenter

    Vec cellCenters

    call VecCreateMPI(SUMB_COMM_WORLD, nCellsLocal(1)*3, &
         PETSC_DETERMINE, cellCenters, ierr)
    call EChk(ierr, __FILE__, __LINE__)
    call VecSetBlockSize(cellCenters, 3, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! compute and store the cell centers
    level = 1
    domainLoop: do nn=1, nDom
       spectralLoop: do sps=1, nTimeIntervalsSpectral
          call setPointers(nn, level, sps)
          do k=2,kl
             do j=2,jl
                do i=2,il
                   iRow = flowDoms(nn, level, sps)%globalCell(i, j, k)
                   ! The location of the cell center is determined
                   ! by averaging the cell coordinates.
                   do n=1,3
                      cellCenter(n) = (x(i-1,j-1,k-1,n) + x(i,j-1,k-1,n)  &
                           +  x(i-1,j,  k-1,n) + x(i,j,  k-1,n)  &
                           +  x(i-1,j-1,k,  n) + x(i,j-1,k,  n)  &
                           +  x(i-1,j,  k,  n) + x(i,j,  k,  n))/8
                   end do
                   call VecSetValues(cellCenters, 1, iRow, cellCenter, INSERT_VALUES, ierr)
                   call EChk(ierr, __FILE__, __LINE__)
                end do
             end do
          end do
       end do spectralLoop
    end do domainLoop

    ! PETSc Matrix Assembly begin
    call VecAssemblyBegin(cellCenters, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecAssemblyEnd  (cellCenters, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call PetscViewerBinaryOpen(sumb_comm_world, fileName, FILE_MODE_WRITE, binViewer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call VecView(cellCenters, binViewer, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call PetscViewerDestroy(binViewer,ierr)
    call EChk(ierr, __FILE__, __LINE__)

  end subroutine saveCellCenters

 subroutine dRdwTMatMult(A, vecX,  vecY, ierr)

    ! PETSc user-defied call back function for computing the product of
    ! dRdwT with a vector. Here we just call the much more broadly
    ! useful routine computeMatrixFreeProductBwdFast()

    use constants
    use communication
    use blockPointers
    use iteration         
    use flowVarRefState     
    use inputAdjoint       
    use ADjointVars
    use inputTimeSpectral
    use utils, only : EChk
    implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"

    ! PETSc Arguments
    Mat   A
    Vec   vecX, vecY
    integer(kind=intType) ::ierr

    real(kind=realType), pointer :: dwb_pointer(:)
    real(kind=realType), pointer :: wb_pointer(:)

#ifndef USE_COMPLEX
    call VecGetArrayReadF90(vecX, dwb_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecGetArrayF90(VecY, wb_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call computeMatrixFreeProductBwdFast(dwb_pointer, wb_pointer, size(wb_pointer))

    call VecRestoreArrayF90(vecX, dwb_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ierr = 0
#endif

  end subroutine dRdwTMatMult

  subroutine dRdwMatMult(A, vecX,  vecY, ierr)

    ! PETSc user-defied call back function for computing the product of
    ! dRdw with a vector. Here we just call the much more broadly
    ! useful routine computeMatrixFreeProductFwd()

    use constants
    use communication
    use blockPointers
    use iteration         
    use flowVarRefState     
    use inputAdjoint       
    use ADjointVars
    use inputTimeSpectral  
    use surfaceFamilies, only: wallFamilies, totalWallFamilies
    use utils, only : EChk
    use surfaceUtils, only : getSurfaceSize, setFulLFamilyList
    use costFunctions
    implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"

    ! PETSc Arguments
    Mat   A
    Vec   vecX, vecY
    integer(kind=intType) ::ierr

    real(kind=realType), pointer :: wd_pointer(:)
    real(kind=realType), pointer :: dwd_pointer(:)
    real(kind=realType) :: funcsBar(nCostFunction)
    logical :: useState, useSpatial
    real(kind=realType) :: extraBarDummy
    integer(kind=intType) :: spatialSize, extraSize
    integer(kind=intType) :: stateSize, costSize, fSize, fSIzeCell
    real(kind=realType), dimension(:), allocatable :: Xvdot
    real(kind=realType), dimension(:, :, :), allocatable :: fDot
    real(kind=realType) :: extraDot(nDesignExtra)
    real(kind=realType) ::funcsDot(nCostFunction)
#ifndef USE_COMPLEX

    call VecGetArrayReadF90(vecX, wd_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecGetArrayF90(VecY, dwd_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    funcsBar = zero
    useSpatial  = .False.
    useState    = .True.
    spatialSize =  3 * nNodesLocal(1_intType)*nTimeIntervalsSpectral
    extraSize   = size(extraDot)
    stateSize   = size(wd_pointer)
    costSize    = nCostFunction

    call getSurfaceSize(fSize, fSizeCell, wallFamilies, totalWallFamilies)
    call setFullFamilyList()
    allocate(xvdot(spatialSize))
    allocate(fdot(3, fSize, nTimeIntervalsSpectral))

    xvdot = zero
    extradot = zero
    fdot = zero
    call computeMatrixFreeProductFwd(xvdot, extradot, wd_pointer, &
         useSpatial, useState, dwd_pointer, funcsDot, fDot, &
         spatialSize, extraSize, stateSize, costSize, fSize, nTimeIntervalsSpectral)
    deallocate(xvdot)

    call VecRestoreArrayReadF90(vecX, wd_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecRestoreArrayF90(VecY, dwd_pointer, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ierr = 0
#endif
  end subroutine dRdwMatMult

  subroutine createPETScVars
    !
    !      Create the matrices/vectors that are required for the adjoint  
    !
    use constants
    use ADjointPETSc, only: dRdwT, dRdwPreT, &
         adjointKSP, matfreectx, x_like, psi_like1, adjointPETScVarsAllocated
    use ADjointVars   
    use communication, only : sumb_comm_world
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use flowVarRefState, only : nwf, nw, viscous
    use inputADjoint, only : approxPC, frozenTurbulence, useMatrixFreedRdw, viscPC
    use stencils, only : N_visc_drdw, n_euler_drdw, visc_drdw_stencil,  euler_drdw_stencil, &
         visc_drdw_stencil, visc_pc_stencil, N_visc_PC, N_euler_PC, euler_PC_stencil
    use utils, only : EChk, setPointers
    use adjointUtils, only : myMatCreate, destroyPETScVars, statePreAllocation
    implicit none

#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petsc.h"

    !     Local variables.
    integer(kind=intType)  :: nDimW, nDimX
    integer(kind=intType) :: i, n_stencil, nState
    integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
    integer(kind=intType), dimension(:), allocatable :: nnzDiagonal2, nnzOffDiag2
    integer(kind=intType), dimension(:, :), pointer :: stencil
    integer(kind=intType) :: level, ierr, nlevels
    integer(kind=intType) :: rows(4), iCol, nn, sps, ii
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, iDim, iStride, j, mm
    integer(kind=intType) :: npts, ncells, nTS

    ! Destroy variables if they already exist
    call destroyPETScVars()

    ! DETERMINE ALL SIZES HERE!
    if ( frozenTurbulence ) then
       nState = nwf
    else
       nState = nw
    endif

    nDimW = nState * nCellsLocal(1_intType)*nTimeIntervalsSpectral
    nDimX = 3 * nNodesLocal(1_intType)*nTimeIntervalsSpectral

    if (.not. useMatrixFreedRdw) then 
       ! Setup matrix-based dRdwT
       allocate(nnzDiagonal(nCellsLocal(1_intType)*nTimeIntervalsSpectral), &
            nnzOffDiag(nCellsLocal(1_intType)*nTimeIntervalsSpectral) )

       if (viscous) then
          n_stencil = N_visc_drdw
          stencil => visc_drdw_stencil
       else
          n_stencil = N_euler_drdw
          stencil => euler_drdw_stencil 
       end if

       level = 1

       call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nState, stencil, n_stencil, &
            level)
       call myMatCreate(dRdwT, nState, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
            __FILE__, __LINE__)

       call matSetOption(dRdwT, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       deallocate(nnzDiagonal, nnzOffDiag)
    else
       ! Setup matrix-free dRdwT
       call MatCreateShell(SUMB_COMM_WORLD, nDimW, nDimW, PETSC_DETERMINE, &
            PETSC_DETERMINE, matfreectx, dRdwT, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set the shell operation for doing matrix vector multiplies
       call MatShellSetOperation(dRdwT, MATOP_MULT, dRdwTMatMult, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Set the shell operation for doing TRNASPOSE matrix vector
       ! multiplies
       call MatShellSetOperation(dRdwT, MATOP_MULT_TRANSPOSE, dRdwMatMult, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       call MatSetup(dRdwT, ierr)
       call EChk(ierr, __FILE__, __LINE__)
    end if

    ! Create the approxPC if required
    if (ApproxPC) then
       ! ------------------- Determine Preallocation for dRdwPre -------------
       allocate(nnzDiagonal(nCellsLocal(1_intType)*nTimeIntervalsSpectral), &
            nnzOffDiag(nCellsLocal(1_intType)*nTimeIntervalsSpectral) )

       if (viscous .and. viscPC) then
          stencil => visc_pc_stencil
          n_stencil = N_visc_pc
       else
          stencil => euler_pc_stencil
          n_stencil = N_euler_pc
       end if

       level = 1
       call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nState, stencil, n_stencil, &
            level)
       call myMatCreate(dRdwPreT, nState, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
            __FILE__, __LINE__)

       call matSetOption(dRdwPreT, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       deallocate(nnzDiagonal, nnzOffDiag)
    end if

    ! Create the KSP Object
    call KSPCreate(SUMB_COMM_WORLD, adjointKSP, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    adjointPETScVarsAllocated = .True.
  end subroutine createPETScVars

end module adjointAPI
