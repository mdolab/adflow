#ifndef USE_COMPLEX
subroutine computeMatrixFreeProductFwd(xvdot, extradot, wdot, useSpatial, useState, dwdot, funcsDot, &
     fDot, spatialSize, extraSize, stateSize, costSize, fSize)

  ! This is the main matrix-free forward mode computation
  use constants
  use bcTypes
  use communication
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

  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

  ! Input Variables
  integer(kind=intType), intent(in) :: spatialSize, extraSize, stateSize, costSize, fSize
  real(kind=realType), dimension(spatialSize), intent(in) :: xvdot
  real(kind=realType), dimension(extraSize), intent(in) :: extradot
  real(kind=realType), dimension(stateSize), intent(in) :: wdot
  logical, intent(in) :: useSpatial, useState

  ! Ouput Variables
  real(kind=realType), dimension(stateSize), intent(out) :: dwDot
  real(kind=realType), dimension(costSize), intent(out) :: funcsDot
  real(kind=realType), dimension(3, fSize), intent(out) :: fDot

  ! Working Variables
  integer(kind=intType) :: ierr,nn,mm,sps,i,j,k,l,ii,jj,idim,sps2
  real(kind=realType) :: alpha, beta, alphad, betad
  integer(kind=intType) ::  level, irow, liftIndex
  real(kind=realType), dimension(costSize) :: funcsLocalDot
  logical :: resetToRans

  ! Determine if we want to use frozenTurbulent Adjoint
  resetToRANS = .False. 
  if (frozenTurbulence .and. equations == RANSEquations) then
     equations = NSEquations 
     resetToRANS = .True.
  end if

  ! Need to trick the residual evalution to use coupled (mean flow and
  ! turbulent) together.
  level = 1
  currentLevel = level
  groundLevel = level

  call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)

  ! Allocate the memory we need for this block to do the forward
  ! mode derivatives and copy reference values
  ! Allocate the memory for reverse
  if (.not. derivVarsAllocated) then 
     call alloc_derivative_values(level)
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
  prefd = zero
  tempfreestreamd = zero
  reynoldsd = zero

  if (useSpatial) then 
     ! Here we set the spatial and extra seeds if necessary.
     ii = 0
     domainLoop1: do nn=1,nDom

        ! Just to get sizes
        call setPointers(nn, level, 1)

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
          Prefd = extraDot(nDesignPressure+1)
     if (nDesignTemperature >= 0) &
          Tempfreestreamd = extraDot(nDesignTemperature+1)
     if (nDesignReynolds >= 0) &
          reynoldsd= extraDot(nDesignReynolds+1)
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
        call setPointers(nn, level, 1)

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

        call BLOCK_RES_D(nn, level, useSpatial, alpha, alphad, beta, betad, &
             & liftindex, frozenTurbulence)

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

        ! And extract fDot
        bocos: do mm=1,nBocos
           if (bctype(mm) .eq. eulerwall .or. &
               bctype(mm) .eq. nswalladiabatic  .or. &
               bctype(mm) .eq. nswallisothermal) then
                 
              ! Loop over the nodes since that's where the forces get
              ! defined.
              do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
                 do i=BCData(mm)%inBeg,BCData(mm)%inEnd
                    jj = jj + 1
                    do iDim=1,3
                       fDot(idim, jj) = bcDatad(mm)%F(i, j, iDim)
                    end do
                 end do
              end do
           end if
        end do bocos
     end do spectalLoopAD
  end do domainLoopAD

  ! We need to allreduce the function values to all processors
  call mpi_allreduce(funcsLocalDot, funcsDot, costSize, sumb_real, mpi_sum, SUmb_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Reset the correct equation parameters if we were useing the frozen
  ! Turbulent 
  if (resetToRANS) then
     equations = RANSEquations
  end if
end subroutine computeMatrixFreeProductFwd

subroutine computeMatrixFreeProductBwd(dwbar, funcsbar, fbar, useSpatial, useState, xvbar, &
     extrabar, wbar, spatialSize, extraSize, stateSize, costSize, fSize)
  use constants
  use bcTypes
  use communication
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
  use walldistancedata, only : xSurfVec, xSurfVecd, xSurf, xSurfd, wallScatter
  implicit none

#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif


  ! Input Variables
  integer(kind=intType), intent(in) :: spatialSize, extraSize, stateSize, costSize, fSize
  real(kind=realType), dimension(stateSize), intent(in) :: dwbar
  real(kind=realType), dimension(costSize), intent(in) :: funcsbar
  real(kind=realType), dimension(3, fSize), intent(in) :: fBar
  logical, intent(in) :: useSpatial, useState

  ! Ouput Variables
  real(kind=realType), dimension(stateSize), intent(out) :: wbar
  real(kind=realType), dimension(extraSize), intent(out) :: extrabar
  real(kind=realType), dimension(spatialSize), intent(out) :: xvbar

  ! Working variables
  integer(kind=intType) :: ierr,nn,mm,sps,i,j,k,l,ii,jj,sps2, idim
  real(kind=realType) :: alpha, beta, alphad, betad
  integer(kind=intType) ::  level, irow, liftIndex, nState
  logical :: resetToRans
  real(kind=realType), dimension(extraSize) :: extraLocalBar
  real(kind=realType), dimension(:), allocatable :: xSurfbSum

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
     call alloc_derivative_values(level)
  end if
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn, level, sps)
        call zeroADSeeds(nn,level, sps)
     end do
  end do

  ! Zero the function seeds
  funcValuesd= zero
  call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)

  ! Now extract the vector of the surface data we need
  call VecGetArrayF90(xSurfVec(level), xSurf, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! And it's derivative
  call VecGetArrayF90(xSurfVecd, xSurfd, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  allocate(xSurfbSum(size(xSurfd)))
  xSurfbSum = zero

  ! Zero out extraLocal
  extraLocalBar = zero

  ii = 0 ! Residual bar counter
  jj = 0 ! Force bar counter
  domainLoopAD: do nn=1,nDom

     ! Just to get sizes
     call setPointers(nn,1_intType,1)
     call setDiffSizes

     do sps=1,nTimeIntervalsSpectral
        ! Set pointers and derivative pointers
        call setPointers_d(nn, level, sps)

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
        
        ! And the Force seeds is spatial is true
        bocos: do mm=1,nBocos
           if (bctype(mm) .eq. eulerwall .or. &
                bctype(mm) .eq. nswalladiabatic  .or. &
                bctype(mm) .eq. nswallisothermal) then
              
              ! Loop over the nodes since that's where the forces get
              ! defined.
              do j=(BCData(mm)%jnBeg),BCData(mm)%jnEnd
                 do i=(BCData(mm)%inBeg),BCData(mm)%inEnd
                    jj = jj + 1
                    do iDim=1,3
                       bcDatad(mm)%F(i, j, iDim) = fBar(idim, jj)
                    end do
                 end do
              end do
           end if
        end do bocos

        call BLOCK_RES_B(nn, sps, useSpatial, alpha, alphad, beta, betad, &
             & liftindex, frozenTurbulence)

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
              xSurfbSum = xSurfbSum + xSurfd
              
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
                   extraLocalBar(nDesignPressure+1) = extraLocalBar(nDesignPressure+1) + prefd
              if (nDesignTemperature >= 0) &
                   extraLocalBar(nDesignTemperature+1) = extraLocalBar(nDesignTemperature+1) + tempfreestreamd
              if (nDesignReynolds >= 0) &
                   extraLocalBar(nDesignReynolds+1) = extraLocalBar(nDesignReynolds+1) + reynoldsd
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
  xSurfd = xSurfbSum
  deallocate(xsurfbSum)

  ! These arrays need to be restored before we can do the spatial scatter below:
  call VecRestoreArrayF90(xSurfVec(level), xSurf, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! And it's derivative
  call VecRestoreArrayF90(xSurfVecd, xSurfd, ierr)
  call EChk(ierr,__FILE__,__LINE__)

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
        call VecScatterBegin(wallScatter(1), xSurfVecd, x_like, ADD_VALUES, SCATTER_REVERSE, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        call VecScatterEnd(wallScatter(1), xSurfVecd, x_like, ADD_VALUES, SCATTER_REVERSE, ierr)
        call EChk(ierr,__FILE__,__LINE__)
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

  use blockPointers
  use inputDiscretization
  use inputTimeSpectral 
  use inputPhysics
  use flowVarRefState
  use inputAdjoint       
  use iteration
  use inputIteration
  use saModule_fast_b
  use bcroutines_b
  use adjointvars
  use communication
  use paramTurb
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

  ! Allocate the memory we need for this block to do the forward
  ! mode derivatives and copy reference values
  !call alloc_derivative_values( level)
  if (.not. derivVarsAllocated) then 
     call alloc_derivative_values(level)
  end if

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
              call returnFail("matrixFreeRoutines", &
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
              call saeddyviscosity_b
           end select
        end if

        if (viscous) then 
           call computelamviscosity_fast_b
        end if

        call computepressuresimple_fast_b

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
  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif


  ! PETSc Arguments
  Mat   A
  Vec   vecX, vecY
  integer(kind=intType) ::ierr

  real(kind=realType), pointer :: dwb_pointer(:)
  real(kind=realType), pointer :: wb_pointer(:)

#ifndef USE_COMPLEX

#if PETSC_VERSION_MINOR > 5
  call VecGetArrayReadF90(vecX, dwb_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#else
  call VecGetArrayF90(vecX, dwb_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif

  call VecGetArrayF90(VecY, wb_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call computeMatrixFreeProductBwdFast(dwb_pointer, wb_pointer, size(wb_pointer))

#if PETSC_VERSION_MINOR > 5
  call VecRestoreArrayReadF90(vecX, dwb_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#else
  call VecRestoreArrayF90(vecX, dwb_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
  call VecRestoreArrayF90(VecY, wb_pointer, ierr)
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

  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h90"
#endif

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
  real(kind=realType), dimension(:, :), allocatable :: fDot
  real(kind=realType) :: extraDot(nDesignExtra)
  real(kind=realType) ::funcsDot(nCostFunction)
#ifndef USE_COMPLEX

#if PETSC_VERSION_MINOR > 5
  call VecGetArrayReadF90(vecX, wd_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#else
  call VecGetArrayF90(vecX, wd_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif

  call VecGetArrayF90(VecY, dwd_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  funcsBar = zero
  useSpatial  = .False.
  useState    = .True.
  spatialSize =  3 * nNodesLocal(1_intType)*nTimeIntervalsSpectral
  extraSize   = size(extraDot)
  stateSize   = size(wd_pointer)
  costSize    = nCostFunction
  call getForceSize(fSize, fSizeCell)
  allocate(xvdot(spatialSize))
  allocate(fdot(3, fSize))
  xvdot = zero
  extradot = zero
  fdot = zero
  call computeMatrixFreeProductFwd(xvdot, extradot, wd_pointer, &
       useSpatial, useState, dwd_pointer, funcsDot, fDot, &
       spatialSize, extraSize, stateSize, costSize, fSize)
  deallocate(xvdot)
#if PETSC_VERSION_MINOR > 5
  call VecRestoreArrayReadF90(vecX, wd_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#else
  call VecRestoreArrayF90(vecX, wd_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
  call VecRestoreArrayF90(VecY, dwd_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ierr = 0
#endif
end subroutine dRdwMatMult

subroutine solveAdjointForRHS(inVec, outVec, nDOF, relativeTolerance)

  !   use ADJointPETSc
  !   use inputADjoint
  !   use adjointvars
  !   use killsignals
  !   use constants
  !   implicit none

  !   ! Input Variables
  !   real(kind=realType), dimension(ndof), intent(in) :: inVec
  !   real(kind=realType), dimension(ndof), intent(out) :: outVec
  !   real(kind=realType), intent(in) :: relativeTolerance
  !   integer(kind=intType), intent(in) :: nDOF

  !   ! Working variables
  !   integer(kind=intType) :: ierr

  ! #ifndef USE_COMPLEX

  !   ! Place the arrays
  !   call VecPlaceArray(psi_like1, inVec, ierr)
  !   call EChk(ierr,__FILE__,__LINE__)

  !   call VecPlaceArray(psi_like2, outVec, ierr)
  !   call EChk(ierr,__FILE__,__LINE__)

  !   ! Zero out initial solution
  !   call VecSet(psi_like2, zero, ierr)
  !   call EChk(ierr,__FILE__,__LINE__)

  !   ! Set desired realtive tolerance
  !   call KSPSetTolerances(adjointKSP, relativeTolerance, adjAbsTol, adjDivTol, &
  !        adjMaxIter, ierr)
  !   call EChk(ierr, __FILE__, __LINE__)

  !   ! Solve (remember this is actually a transpose solve)
  !   call KSPSolve(adjointKSP, psi_like1, psi_like2, ierr)
  !   call EChk(ierr, __FILE__, __LINE__)

  !   call KSPGetConvergedReason(adjointKSP, adjointConvergedReason, ierr)
  !   call EChk(ierr, __FILE__, __LINE__)

  !   if (adjointConvergedReason ==  KSP_CONVERGED_RTOL .or. &
  !        adjointConvergedReason ==  KSP_CONVERGED_ATOL .or. &
  !        adjointConvergedReason ==  KSP_CONVERGED_HAPPY_BREAKDOWN) then
  !      adjointFailed = .False.
  !   else
  !      adjointFailed = .True.
  !   end if

  !   ! Rest arrays
  !   call VecResetArray(psi_like1,  ierr)
  !   call EChk(ierr,__FILE__,__LINE__)

  !   call VecResetArray(psi_like2,  ierr)
  !   call EChk(ierr,__FILE__,__LINE__)

  ! #endif

end subroutine solveAdjointForRHS

subroutine solveDirectForRHS(inVec, outVec, nDOF, relativeTolerance)

  !   use ADJointPETSc
  !   use inputADjoint
  !   use adjointVars
  !   use constants
  !   use killsignals
  !   implicit none

  !   ! Input Variables
  !   real(kind=realType), dimension(ndof), intent(in) :: inVec
  !   real(kind=realType), dimension(ndof), intent(out) :: outVec
  !   real(kind=realType), intent(in) :: relativeTolerance
  !   integer(kind=intType), intent(in) :: nDOF

  !   ! Working variables
  !   integer(kind=intType) :: ierr

  ! #ifndef USE_COMPLEX

  !   ! Place the arrays
  !   call VecPlaceArray(psi_like1, inVec, ierr)
  !   call EChk(ierr,__FILE__,__LINE__)

  !   call VecPlaceArray(psi_like2, outVec, ierr)
  !   call EChk(ierr,__FILE__,__LINE__)

  !   ! Zero out initial solution
  !   call VecSet(psi_like2, zero, ierr)
  !   call EChk(ierr,__FILE__,__LINE__)

  !   ! Set desired realtive tolerance
  !   call KSPSetTolerances(adjointKSP, relativeTolerance, adjAbsTol, adjDivTol, &
  !        adjMaxIter, ierr)
  !   call EChk(ierr, __FILE__, __LINE__)

  !   ! Solve (this is the transpose solve of a transpose matrix, so it's direct)
  !   call KSPSolveTranspose(adjointKSP, psi_like1, psi_like2, ierr)
  !   call EChk(ierr, __FILE__, __LINE__)

  !   call KSPGetConvergedReason(adjointKSP, adjointConvergedReason, ierr)
  !   call EChk(ierr, __FILE__, __LINE__)

  !   if (adjointConvergedReason ==  KSP_CONVERGED_RTOL .or. &
  !        adjointConvergedReason ==  KSP_CONVERGED_ATOL .or. &
  !        adjointConvergedReason ==  KSP_CONVERGED_HAPPY_BREAKDOWN) then
  !      adjointFailed = .False.
  !   else
  !      adjointFailed = .True.
  !   end if

  !   ! Rest arrays
  !   call VecResetArray(psi_like1, ierr)
  !   call EChk(ierr,__FILE__,__LINE__)

  !   call VecResetArray(psi_like2, ierr)
  !   call EChk(ierr,__FILE__,__LINE__)

  ! #endif

end subroutine solveDirectForRHS


subroutine testdRdWTSpeed(X, Y, n)


  use ADjointPETSc, only : dRdwT, psi_like1, psi_like2
  use constants
  implicit none

  real(kind=realType), dimension(n) :: X, Y
  integer(kind=intType) :: ierr, n

  ! Simple debugging routine to just do a matrix-vector product with
  ! dRdwT to see how long it takes. 

  call VecPlaceArray(psi_like1, X, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecPlaceArray(psi_like2, Y, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Get Current Residual -- we always solve for the delta
  call MatMult(dRdWT, psi_like1, psi_like2, ierr) 
  call EChk(ierr,__FILE__,__LINE__)


  call VecResetArray(psi_like1, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  call VecResetArray(psi_like2, ierr)
  call EChk(ierr,__FILE__,__LINE__)


end subroutine testdRdWTSpeed
