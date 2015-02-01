subroutine computeMatrixFreeProductFwd(xvdot, extradot, wdot, useSpatial, useState, dwdot, funcsDot, &
     spatialSize, extraSize, stateSize, costSize)

  ! This is the main matrix-free forward mode computation
  use constants
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
#include "finclude/petsc.h"

  ! Input Variables
  integer(kind=intType), intent(in) :: spatialSize, extraSize, stateSize, costSize
  real(kind=realType), dimension(spatialSize), intent(in) :: xvdot
  real(kind=realType), dimension(extraSize), intent(in) :: extradot
  real(kind=realType), dimension(stateSize), intent(in) :: wdot
  logical, intent(in) :: useSpatial, useState

  ! Ouput Variables
  real(kind=realType), dimension(stateSize), intent(out) :: dwDot
  real(kind=realType), dimension(costSize), intent(out) :: funcsDot

  ! Working Variables
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii, sps2
  real(kind=realType) :: alpha, beta, force(3), moment(3), sepSensor, cavitation, cavitationd
  real(kind=realType) :: alphad, betad, forced(3), momentd(3), sepSensord
  integer(kind=intType) ::  level, irow, liftIndex
  real(kind=realType), dimension(costSize) :: funcsLocalDot

  logical :: resetToRans

#ifndef USE_COMPLEX

  if (equations == RANSEquations) then
     nMGVar = nw
     nt1MG = nt1
     nt2MG = nt2

     turbSegregated = .False.
     turbCoupled = .True.
  end if

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
  call alloc_derivative_values(level)

  ! All arrays are zeroed in alloc_deriv_values, but we still need to
  ! zero the extra variables here
  alphad = zero
  betad = zero
  machd = zero
  machGridd = zero
  machcoefd = zero
  pointrefd = zero
  lengthrefd = zero
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
        end do spectalLoop1
     end do domainLoop1

     ! Now run the halo exchange for the nodes
     call exchangecoord(level)

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

     ! Now run the halo exchange:
     call whalo1to1d(level, 1, nw, .False., .False., .False., .False., &
          commPatternCell_2nd, internalCell_2nd)
     ! Technically we need to call all the *other* whalo calls (sliding,
     ! overset) etc, however we have not implemneted those. 
  end if

  funcsLocalDot = zero

  ! Now we are ready to call block_res_d with all the correct seeds
  ii = 0
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
             & liftindex, force, forced, moment, momentd, sepsensor, sepsensord, &
             & cavitation, cavitationd, frozenTurbulence)

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

     end do spectalLoopAD
  end do domainLoopAD

  ! We need to allreduce the function values to all processors
  call mpi_allreduce(funcsLocalDot, funcsDot, costSize, sumb_real, mpi_sum, SUmb_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Finish up the tear-down of this routine
  call dealloc_derivative_values(level)

  ! Reset the correct equation parameters if we were useing the frozen
  ! Turbulent 
  if (resetToRANS) then
     equations = RANSEquations
  end if

  ! Reset the paraters to use segrated turbulence solve. 
  if (equations == RANSEquations) then
     nMGVar = nwf
     nt1MG = nwf + 1
     nt2MG = nwf

     turbSegregated = .True.
     turbCoupled = .False.
     restrictEddyVis = .false.
     if( eddyModel ) restrictEddyVis = .true.
  end if

#endif

end subroutine computeMatrixFreeProductFwd

subroutine computeMatrixFreeProductBwd(dwbar, funcsbar, useSpatial, useState, xvbar, extrabar, wbar,&
     spatialSize, extraSize, stateSize, costSize)
  use constants
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
#include "finclude/petsc.h"
#include "finclude/petscvec.h90"

  ! Input Variables
  integer(kind=intType), intent(in) :: spatialSize, extraSize, stateSize, costSize
  real(kind=realType), dimension(stateSize), intent(in) :: dwbar
  real(kind=realType), dimension(costSize), intent(in) :: funcsbar
  logical, intent(in) :: useSpatial, useState

  ! Ouput Variables
  real(kind=realType), dimension(stateSize), intent(out) :: wbar
  real(kind=realType), dimension(extraSize), intent(out) :: extrabar
  real(kind=realType), dimension(spatialSize), intent(out) :: xvbar

  ! Working variables
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii, sps2
  real(kind=realType) :: alpha, beta, force(3), moment(3), sepSensor, cavitation
  real(kind=realType) :: alphad, betad, forced(3), momentd(3), sepSensord, cavitationd
  integer(kind=intType) ::  level, irow, liftIndex, nState
  logical :: resetToRans
  real(kind=realType), dimension(extraSize) :: extraLocalBar
  real(kind=realType), dimension(:), allocatable :: xSurfbSum
#ifndef USE_COMPLEX
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

  ! If we are computing the jacobian for the RANS equations, we need
  ! to make block_res think that we are evauluating the residual in a
  ! fully coupled sense.  This is reset after this routine is
  ! finished.
  if (equations == RANSEquations) then
     nMGVar = nw
     nt1MG = nt1
     nt2MG = nt2

     turbSegregated = .False.
     turbCoupled = .True.
  end if

  ! Determine if we want to use frozenTurbulent Adjoint
  resetToRANS = .False. 
  if (frozenTurbulence .and. equations == RANSEquations) then
     equations = NSEquations 
     resetToRANS = .True.
  end if

  ! Allocate the memory for reverse
  call alloc_derivative_values(level)

  ! Zero the function seeds
  forced= zero
  momentd= zero
  sepSensord= zero
  cavitationd= zero
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
  
  ii = 0
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

        call BLOCK_RES_B(nn, 1, useSpatial, alpha, alphad, beta, betad, &
             & liftindex, force, forced, moment, momentd, sepsensor, sepsensord, &
             & cavitation, cavitationd, frozenTurbulence)

        ! Assmeble the vectors requested:
        
        do sps2=1,nTimeIntervalsSpectral
           if (useSpatial) then 
              do k=0,ke
                 do j=0,je
                    do i=0,ie
                       do l=1,3
                          irow = flowDoms(nn, 1, sps2)%globalNode(i,j,k)*3 + l -1
                          if (irow >= 0) then 
                             call VecSetValues(x_like, 1, (/irow/), &
                                  (/flowdomsd(nn, level, sps)%x(i, j, k, l)/), ADD_VALUES, ierr)
                             call EChk(ierr,__FILE__,__LINE__)
                          end if
                       end do
                    end do
                 end do
              end do

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
 

  ! Finally we have to do an mpi all reduce on the local parts:
  if (useSpatial) then 
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

  call dealloc_derivative_values(level)

  ! Reset the correct equation parameters if we were useing the frozen
  ! Turbulent 
  if (resetToRANS) then
     equations = RANSEquations
  end if

  ! Reset the paraters to use segrated turbulence solve. 
  if (equations == RANSEquations) then
     nMGVar = nwf
     nt1MG = nwf + 1
     nt2MG = nwf

     turbSegregated = .True.
     turbCoupled = .False.
     restrictEddyVis = .false.
     if( eddyModel ) restrictEddyVis = .true.
  end if
#endif

end subroutine computeMatrixFreeProductBwd

subroutine whalo1to1d(level, start, end, commPressure,       &
     commVarGamma, commLamVis, commEddyVis, &
     commPattern, internal)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Modification to whalo1to1 that communicates derivative values  *
  !      *                                                                *
  !      ******************************************************************
  !
  use block
  use communication
  use inputTimeSpectral
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level, start, end
  logical, intent(in) :: commPressure, commVarGamma
  logical, intent(in) :: commLamVis, commEddyVis

  type(commType), dimension(*), intent(in)         :: commPattern
  type(internalCommType), dimension(*), intent(in) :: internal
  !
  !      Local variables.
  !
  integer :: size, procID, ierr, index
  integer, dimension(mpi_status_size) :: status

  integer(kind=intType) :: nVar, mm
  integer(kind=intType) :: i, j, k, ii, jj
  integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

  logical :: correctPeriodic
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Set the logical correctPeriodic. Only if a momentum variable
  ! is communicated it is needed to apply the periodic
  ! transformations.

  correctPeriodic = .false.
  if(start <= ivx .and. end >= ivz) correctPeriodic = .true.

  ! Determine the number of variables per cell to be sent.

  nVar = max(0_intType,(end - start + 1))
  if( commPressure )  nVar = nVar + 1
  if( commVarGamma ) nVar = nVar + 1
  if( commLamVis )   nVar = nVar + 1
  if( commEddyVis )  nVar = nVar + 1

  if(nVar == 0) return

  ! Loop over the number of spectral solutions.

  spectralModes: do mm=1,nTimeIntervalsSpectral

     ! Send the variables. The data is first copied into
     ! the send buffer after which the buffer is sent asap.

     ii = 1
     sends: do i=1,commPattern(level)%nProcSend

        ! Store the processor id and the size of the message
        ! a bit easier.

        procID = commPattern(level)%sendProc(i)
        size    = nVar*commPattern(level)%nsend(i)

        ! Copy the data in the correct part of the send buffer.

        jj = ii
        do j=1,commPattern(level)%nsend(i)

           ! Store the block id and the indices of the donor
           ! a bit easier.

           d1 = commPattern(level)%sendList(i)%block(j)
           i1 = commPattern(level)%sendList(i)%indices(j,1)
           j1 = commPattern(level)%sendList(i)%indices(j,2)
           k1 = commPattern(level)%sendList(i)%indices(j,3)

           ! Copy the given range of the working variables for
           ! this cell in the buffer. Update the counter jj.

           do k=start,end
              sendBuffer(jj) = flowDomsd(d1,level,mm)%w(i1,j1,k1,k)
              jj = jj + 1
           enddo

           ! The pressure, if needed.

           if( commPressure ) then
              sendBuffer(jj) = flowDomsd(d1,level,mm)%p(i1,j1,k1)
              jj = jj + 1
           endif

           ! The specific heat ratio, if needed. Note that level == 1.

           if( commVarGamma ) then
              sendBuffer(jj) = flowDomsd(d1,1,mm)%gamma(i1,j1,k1)
              jj = jj + 1
           endif

           ! The laminar viscosity for a viscous computation.
           ! Again level == 1.

           if( commLamVis ) then
              sendBuffer(jj) = flowDomsd(d1,1,mm)%rlv(i1,j1,k1)
              jj = jj + 1
           endif

           ! The eddy viscosity for eddy viscosity models.
           ! Level is the true multigrid level, because the eddy
           ! viscosity is allocated on all grid levels.

           if( commEddyVis ) then
              sendBuffer(jj) = flowDomsd(d1,level,mm)%rev(i1,j1,k1)
              jj = jj + 1
           endif

        enddo

        ! Send the data.

        call mpi_isend(sendBuffer(ii), size, sumb_real, procID,  &
             procID, SUmb_comm_world, sendRequests(i), &
             ierr)

        ! Set ii to jj for the next processor.

        ii = jj

     enddo sends

     ! Post the nonblocking receives.

     ii = 1
     receives: do i=1,commPattern(level)%nProcRecv

        ! Store the processor id and the size of the message
        ! a bit easier.

        procID = commPattern(level)%recvProc(i)
        size    = nVar*commPattern(level)%nrecv(i)

        ! Post the receive.

        call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
             myID, SUmb_comm_world, recvRequests(i), ierr)

        ! And update ii.

        ii = ii + size

     enddo receives

     ! Copy the local data.

     localCopy: do i=1,internal(level)%ncopy

        ! Store the block and the indices of the donor a bit easier.

        d1 = internal(level)%donorBlock(i)
        i1 = internal(level)%donorIndices(i,1)
        j1 = internal(level)%donorIndices(i,2)
        k1 = internal(level)%donorIndices(i,3)

        ! Idem for the halo's.

        d2 = internal(level)%haloBlock(i)
        i2 = internal(level)%haloIndices(i,1)
        j2 = internal(level)%haloIndices(i,2)
        k2 = internal(level)%haloIndices(i,3)

        ! Copy the given range of working variables.

        do k=start,end
           flowDomsd(d2,level,mm)%w(i2,j2,k2,k) = &
                flowDomsd(d1,level,mm)%w(i1,j1,k1,k)
        enddo

        ! The pressure, if needed.

        if( commPressure )                   &
             flowDomsd(d2,level,mm)%p(i2,j2,k2) = &
             flowDomsd(d1,level,mm)%p(i1,j1,k1)

        ! The specific heat ratio, if needed. Note that level == 1.

        if( commVarGamma )                   &
             flowDomsd(d2,1,mm)%gamma(i2,j2,k2) = &
             flowDomsd(d1,1,mm)%gamma(i1,j1,k1)

        ! The laminar viscosity for viscous computations.
        ! Again level == 1.

        if( commLamVis )                   &
             flowDomsd(d2,1,mm)%rlv(i2,j2,k2) = &
             flowDomsd(d1,1,mm)%rlv(i1,j1,k1)

        ! The eddy viscosity for eddy viscosity models.
        ! Level is the true multigrid level, because the eddy
        ! viscosity is allocated on all grid levels.

        if( commEddyVis )                      &
             flowDomsd(d2,level,mm)%rev(i2,j2,k2) = &
             flowDomsd(d1,level,mm)%rev(i1,j1,k1)

     enddo localCopy

     ! Correct the periodic halo's of the internal communication
     ! pattern, if needed.

     !if(correctPeriodic .and. internal(level)%nPeriodic > 0)   &
                                !NOT IMPLEMENTED YET
                                ! call correctPeriodicVelocity(level, mm,                 &
                                !                              internal(level)%nPeriodic, &
                                !                              internal(level)%periodicData)
          
                                ! Complete the nonblocking receives in an arbitrary sequence and
                                ! copy the variables from the buffer into the halo's.
          
          size = commPattern(level)%nProcRecv
     completeRecvs: do i=1,commPattern(level)%nProcRecv

        ! Complete any of the requests.

        call mpi_waitany(size, recvRequests, index, status, ierr)

        ! Copy the data just arrived in the halo's.

        ii = index
        jj = nVar*commPattern(level)%nrecvCum(ii-1)
        do j=1,commPattern(level)%nrecv(ii)

           ! Store the block and the indices of the halo a bit easier.

           d2 = commPattern(level)%recvList(ii)%block(j)
           i2 = commPattern(level)%recvList(ii)%indices(j,1)
           j2 = commPattern(level)%recvList(ii)%indices(j,2)
           k2 = commPattern(level)%recvList(ii)%indices(j,3)

           ! Copy the conservative variables.

           do k=start,end
              jj = jj + 1
              flowDomsd(d2,level,mm)%w(i2,j2,k2,k) = recvBuffer(jj)
           enddo

           ! The pressure, if needed.

           if( commPressure ) then
              jj = jj + 1
              flowDomsd(d2,level,mm)%p(i2,j2,k2) = recvBuffer(jj)
           endif

           ! The specific heat ratio, if needed. Note that level == 1.

           if( commVarGamma ) then
              jj = jj + 1
              flowDomsd(d2,1,mm)%gamma(i2,j2,k2) = recvBuffer(jj)
           endif

           ! The laminar viscosity for viscous computations.
           ! Again level == 1.

           if( commLamVis ) then
              jj = jj + 1
              flowDomsd(d2,1,mm)%rlv(i2,j2,k2) = recvBuffer(jj)
           endif

           ! The eddy viscosity ratio for eddy viscosity models.
           ! Level is the true multigrid level, because the eddy
           ! viscosity is allocated on all grid levels.

           if( commEddyVis ) then
              jj = jj + 1
              flowDomsd(d2,level,mm)%rev(i2,j2,k2) = recvBuffer(jj)
           endif

        enddo

     enddo completeRecvs

     ! Correct the periodic halo's of the external communication
     ! pattern, if needed.

     !if(correctPeriodic .and. commPattern(level)%nPeriodic > 0)   &
                                !NOT IMLEMENTED
                                ! call correctPeriodicVelocity(level, mm,                    &
                                !                              commPattern(level)%nPeriodic, &
                                !                              commPattern(level)%periodicData)
          
                                ! Complete the nonblocking sends.
          
          size = commPattern(level)%nProcSend
     do i=1,commPattern(level)%nProcSend
        call mpi_waitany(size, sendRequests, index, status, ierr)
     enddo

  enddo spectralModes

end subroutine whalo1to1d

subroutine exchangeCoord(level)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * ExchangeCoor exchanges the coordinates of the given grid       *
  !      * level.                                                         *
  !      *                                                                *
  !      ******************************************************************
  !
  use block
  use communication
  use inputTimeSpectral
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level
  !
  !      Local variables.
  !
  integer :: size, procID, ierr, index
  integer, dimension(mpi_status_size) :: status

  integer(kind=intType) :: i, j, ii, jj, mm
  integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Loop over the number of spectral solutions.

  spectralLoop: do mm=1,nTimeIntervalsSpectral

     ! Send the coordinates i have to send. The data is first copied
     ! into the send buffer and this buffer is sent.

     ii = 1
     sends: do i=1,commPatternNode_1st(level)%nProcSend

        ! Store the processor id and the size of the message
        ! a bit easier.

        procID = commPatternNode_1st(level)%sendProc(i)
        size   = 3*commPatternNode_1st(level)%nSend(i)

        ! Copy the data in the correct part of the send buffer.

        jj = ii
        do j=1,commPatternNode_1st(level)%nSend(i)

           ! Store the block id and the indices of the donor
           ! a bit easier.

           d1 = commPatternNode_1st(level)%sendList(i)%block(j)
           i1 = commPatternNode_1st(level)%sendList(i)%indices(j,1)
           j1 = commPatternNode_1st(level)%sendList(i)%indices(j,2)
           k1 = commPatternNode_1st(level)%sendList(i)%indices(j,3)

           ! Copy the coordinates of this point in the buffer.
           ! Update the counter jj accordingly.

           sendBuffer(jj)   = flowDomsd(d1,level,mm)%x(i1,j1,k1,1)
           sendBuffer(jj+1) = flowDomsd(d1,level,mm)%x(i1,j1,k1,2)
           sendBuffer(jj+2) = flowDomsd(d1,level,mm)%x(i1,j1,k1,3)
           jj = jj + 3

        enddo

        ! Send the data.

        call mpi_isend(sendBuffer(ii), size, sumb_real, procID,    &
             procID, SUmb_comm_world, sendRequests(i), &
             ierr)

        ! Set ii to jj for the next processor.

        ii = jj

     enddo sends

     ! Post the nonblocking receives.

     ii = 1
     receives: do i=1,commPatternNode_1st(level)%nProcRecv

        ! Store the processor id and the size of the message
        ! a bit easier.

        procID = commPatternNode_1st(level)%recvProc(i)
        size   = 3*commPatternNode_1st(level)%nRecv(i)

        ! Post the receive.

        call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
             myID, SUmb_comm_world, recvRequests(i), ierr)

        ! And update ii.

        ii = ii + size

     enddo receives

     ! Copy the local data.

     localCopy: do i=1,internalNode_1st(level)%nCopy

        ! Store the block and the indices of the donor a bit easier.

        d1 = internalNode_1st(level)%donorBlock(i)
        i1 = internalNode_1st(level)%donorIndices(i,1)
        j1 = internalNode_1st(level)%donorIndices(i,2)
        k1 = internalNode_1st(level)%donorIndices(i,3)
        ! Idem for the halo's.

        d2 = internalNode_1st(level)%haloBlock(i)
        i2 = internalNode_1st(level)%haloIndices(i,1)
        j2 = internalNode_1st(level)%haloIndices(i,2)
        k2 = internalNode_1st(level)%haloIndices(i,3)
        ! Copy the coordinates.
        flowDomsd(d2,level,mm)%x(i2,j2,k2,1) = &
             flowDomsd(d1,level,mm)%x(i1,j1,k1,1)
        flowDomsd(d2,level,mm)%x(i2,j2,k2,2) = &
             flowDomsd(d1,level,mm)%x(i1,j1,k1,2)
        flowDomsd(d2,level,mm)%x(i2,j2,k2,3) = &
             flowDomsd(d1,level,mm)%x(i1,j1,k1,3)

     enddo localCopy

     ! Correct the periodic halos of the internal communication
     ! pattern

     ! NOT IMPLEMENTED
     ! call correctPeriodicCoor(level, mm,                          &
     !      internalNode_1st(level)%nPeriodic,  &
     !      internalNode_1st(level)%periodicData)

     ! Complete the nonblocking receives in an arbitrary sequence and
     ! copy the coordinates from the buffer into the halo's.

     size = commPatternNode_1st(level)%nProcRecv
     completeRecvs: do i=1,commPatternNode_1st(level)%nProcRecv

        ! Complete any of the requests.

        call mpi_waitany(size, recvRequests, index, status, ierr)

        ! Copy the data just arrived in the halo's.

        ii = index
        jj = 3*commPatternNode_1st(level)%nRecvCum(ii-1) +1
        do j=1,commPatternNode_1st(level)%nRecv(ii)

           ! Store the block and the indices of the halo a bit easier.

           d2 = commPatternNode_1st(level)%recvList(ii)%block(j)
           i2 = commPatternNode_1st(level)%recvList(ii)%indices(j,1)
           j2 = commPatternNode_1st(level)%recvList(ii)%indices(j,2)
           k2 = commPatternNode_1st(level)%recvList(ii)%indices(j,3)

           ! Copy the data.

           flowDomsd(d2,level,mm)%x(i2,j2,k2,1) = recvBuffer(jj)
           flowDomsd(d2,level,mm)%x(i2,j2,k2,2) = recvBuffer(jj+1)
           flowDomsd(d2,level,mm)%x(i2,j2,k2,3) = recvBuffer(jj+2)
           jj = jj + 3

        enddo

     enddo completeRecvs

     ! Correct the periodic halos of the external communication
     ! pattern.
     ! NOT IMLEMENTED
     ! call correctPeriodicCoor(level, mm,                            &
     !      commPatternNode_1st(level)%nPeriodic, &
     !      commPatternNode_1st(level)%periodicData)

     ! Complete the nonblocking sends.

     size = commPatternNode_1st(level)%nProcSend
     do i=1,commPatternNode_1st(level)%nProcSend
        call mpi_waitany(size, sendRequests, index, status, ierr)
     enddo

  enddo spectralLoop

end subroutine exchangeCoord

subroutine dRdwTMatMult(A, vecX,  vecY, ierr)

  ! PETSc user-defied call back function for computing the product of
  ! dRdwT with a vector. Here we just call the much more broadly
  ! useful routine computeMatrixFreeProductBwd()

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
#include "finclude/petsc.h"
#include "finclude/petscvec.h90"
  ! PETSc Arguments
  Mat   A
  Vec   vecX, vecY
  integer(kind=intType) ::ierr

  real(kind=realType), pointer :: dwb_pointer(:)
  real(kind=realType), pointer :: wb_pointer(:)
  real(kind=realType) :: funcsBar(nCostFunction)
  logical :: useState, useSpatial
  real(kind=realType) :: extraBar(0)
  integer(kind=intType) :: spatialSize, extraSize
  integer(kind=intType) :: stateSize, costSize
  real(kind=realType), dimension(:), allocatable :: Xvbar
#ifndef USE_COMPLEX

  call VecGetArrayF90(vecX, dwb_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecGetArrayF90(VecY, wb_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  funcsBar = zero
  useSpatial  = .False.
  useState    = .True.
  spatialSize =  3 * nNodesLocal(1_intType)*nTimeIntervalsSpectral
  extraSize   = 0
  stateSize   = size(wb_pointer)
  costSize    = nCostFunction
  allocate(xvbar(spatialSize))
  call computeMatrixFreeProductBwd(dwb_pointer, funcsbar, &
       useSpatial, useState, xvbar, extrabar, wb_pointer, &
       spatialSize, extraSize, stateSize, costSize)
  deallocate(xvbar)
  call VecRestoreArrayF90(vecX, dwb_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)

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
#include "finclude/petsc.h"
#include "finclude/petscvec.h90"

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
  integer(kind=intType) :: stateSize, costSize
  real(kind=realType), dimension(:), allocatable :: Xvdot
  real(kind=realType) :: extraDot(nDesignExtra)
  real(kind=realType) ::funcsDot(nCostFunction)
#ifndef USE_COMPLEX

  call VecGetArrayF90(vecX, wd_pointer, ierr)
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
  
  allocate(xvdot(spatialSize))
  xvdot = zero
  extradot = zero
  call computeMatrixFreeProductFwd(xvdot, extradot, wd_pointer, &
       useSpatial, useState, dwd_pointer, funcsDot, &
       spatialSize, extraSize, stateSize, costSize)
  deallocate(xvdot)

  call VecRestoreArrayF90(vecX, wd_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)

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
