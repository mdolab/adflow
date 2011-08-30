! This files contains several utilitiy functions that are used with
! the NK solver.

! 1. setWVec: set the PETSc vector wVec from the current w. Only used
!             for initialization
! 2. setRVec: sets the PETSc residual vector rVec from the computed dw
! 3. setW : sets the SUmb state vector w from the petsc vector wVec

! 4. getStates: (Python Wrapped) Returns all the states, w, to Python
! 5. getRes : (Python Wrapped) Computes the residual and returns the
!             result to Python
! 6. setStates: (Python Wrapped) Sets an input vector into the states,
! w. Only used for the coupled NK Aerostructural solver. 

subroutine setWVec(wVec)

  ! Set the petsc vector wVec from the current SUmb soltuion in w
  use communication
  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate 
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  Vec     wVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l
  real(kind=realType) :: states(nw)

  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointersAdj(nn,1_intType,sps)
        ! Copy off w to wVec
        do k=2,kl
           do j=2,jl
              do i=2,il
                 do l=1,nw
                    states(l) = w(i,j,k,l)
                 end do
                 call VecSetValuesBlocked(wVec,1,globalCell(i,j,k),states,&
                      INSERT_VALUES,ierr)
                 call EChk(ierr,__FILE__,__LINE__)
              end do
           end do
        end do
     end do
  end do
  call VecAssemblyBegin(wVec,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecAssemblyEnd(wVec,ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine setWVec

subroutine setRVec(rVec)
  ! Set the current residual in dw into the PETSc Vector

  use communication
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  Vec     rVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l
  real(kind=realType) :: ovv,temp(nw)

  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointersAdj(nn,1_intType,sps)
        ! Copy off dw/vol to rVec
        do k=2,kl
           do j=2,jl
              do i=2,il
                 ovv = 1/vol(i,j,k)
                 do l=1,nwf
                    temp(l) = dw(i,j,k,l)*ovv
                 end do

                 do l=nt1,nt2
                    temp(l) = dw(i,j,k,l)*ovv!*1e-3
                 end do

                 call VecSetValuesBlocked(rVec,1,globalCell(i,j,k),&
                      dw(i,j,k,:)*ovv, INSERT_VALUES,ierr)
                 call EChk(ierr,__FILE__,__LINE__)
              end do
           end do
        end do
     end do
  end do

  call VecAssemblybegin(rVec,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecAssemblyEnd(rVec,ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine setRVec

subroutine setW(wVec)

  ! Set the SUmb state vector, w, from the petsc vec wVec
  use communication
  use blockPointers
  use inputTimeSpectral
  use flowVarRefState
  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  Vec     wVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l
  real(kind=realType) :: temp,diff

  ! Note this is not ideal memory access but the values are stored by
  ! block in PETSc (grouped in nw) but are stored separately in SUmb
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1_intType,sps)
        ! Copy off w to wVec
        do k=2,kl
           do j=2,jl
              do i=2,il
                 do l=1,nw
                    call VecGetValues(wVec,1,globalCell(i,j,k)*nw+l-1,&
                         w(i,j,k,l),ierr)
                 end do
              end do
           end do
        end do
     end do
  end do
end subroutine setW


subroutine getStates(states,ndimw)

  ! Return the state vector, w to Python

  use ADjointPETSc
  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate
  implicit none

  integer(kind=intType),intent(in):: ndimw
  real(kind=realType),dimension(ndimw),intent(out) :: states(ndimw)

  ! Local Variables
  integer(kind=intType) :: nn,i,j,k,l,counter,sps

  counter = 0 
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1,sps)
        do k=2,kl
           do j=2,jl
              do i=2,il
                 do l=1,nw
                    counter = counter + 1
                    states(counter) = w(i,j,k,l)
                 end do
              end do
           end do
        end do
     end do
  end do
end subroutine getStates

subroutine getRes(res,ndimw)

  ! Compute the residual and return result to Python

  use ADjointPETSc
  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate
  implicit none

  integer(kind=intType),intent(in):: ndimw
  real(kind=realType),dimension(ndimw),intent(out) :: res(ndimw)

  ! Local Variables
  integer(kind=intType) :: nn,i,j,k,l,counter,sps

  call computeResidualNK()
  counter = 0 
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1,sps)
        do k=2,kl
           do j=2,jl
              do i=2,il
                 do l=1,nw
                    counter = counter + 1
                    res(counter) = dw(i,j,k,l)/vol(i,j,k)
                 end do
              end do
           end do
        end do
     end do
  end do
end subroutine getRes

subroutine setStates(states,ndimw)

  ! Take in externallly generated states and set them in SUmb

  use ADjointPETSc
  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate
  implicit none

  integer(kind=intType),intent(in):: ndimw
  real(kind=realType),dimension(ndimw),intent(in) :: states(ndimw)

  ! Local Variables
  integer(kind=intType) :: nn,i,j,k,l,counter,sps

  counter = 0 
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1,sps)
        do k=2,kl
           do j=2,jl
              do i=2,il
                 do l=1,nw
                    counter = counter + 1
                    w(i,j,k,l) = states(counter)
                 end do
              end do
           end do
        end do
     end do
  end do
end subroutine setStates


subroutine weak_scaling_test(useComm,setPETSCVecs,niterations)


  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate
  use iteration
  use inputPhysics 
  use nksolvervars
  use communication
  implicit none

  ! do a weak scaling test by running the computeResidual NK function a number of times:

  logical, intent(in) :: useComm, setPETScVecs
  integer(kind=intType) :: niterations

  ! Local Variables
  integer(kind=intType) :: ierr,i,j,k,l,sps,nn,iter
  logical secondHalo ,correctForK
  real(kind=realType) :: gm1,v2,val
  real(kind=realType) :: comm_time,time(4),total_time_local,total_time
  real(kind=realType),dimension(:),allocatable :: all_times
  secondHalo = .True. 
  currentLevel = 1_intType
  groundLevel = 1_intTYpe
  ! Next we need to compute the pressures
  gm1 = gammaConstant - one
  correctForK = .False.

  ! --------------- Master Loop -------------------
  comm_time  = 0.0
  if (myid == 0) then
     print *,'Running ',niterations, ' of the residual routine'
  end if

  call mpi_barrier(sumb_comm_world,ierr)
  
  time(3) = mpi_wtime()

  do iter=1,niterations

     if (setPETScVecs) then
        call setW(wVec)
     end if

     spectralLoop: do sps=1,nTimeIntervalsSpectral
        domainsState: do nn=1,nDom
           ! Set the pointers to this block.
           call setPointers(nn, currentLevel, sps)

           do k=2,kl
              do j=2,jl
                 do i=2,il

                    v2 = w(i,j,k,ivx)**2 + w(i,j,k,ivy)**2 &
                         + w(i,j,k,ivz)**2

                    p(i,j,k) = gm1*(w(i,j,k,irhoE) &
                         - half*w(i,j,k,irho)*v2)
                    p(i,j,k) = max(p(i,j,k), 1.e-4_realType*pInfCorr)
                 enddo
              enddo
           enddo

           call computeEtot(2_intType,il, 2_intType,jl, &
                2_intType,kl, correctForK)

        end do domainsState
     end do spectralLoop


     call computeLamViscosity
     call computeEddyViscosity

     !   Apply BCs
     call applyAllBC(secondHalo)

     ! Exchange solution -- always the fine level
     if (useComm) then
        time(1) = mpi_wtime()
        call whalo1(1_intType, 1_intType, nMGVar, .true., &
             .true., .true.)

        if (equations == RANSEquations) then
           call whalo2(1_intType, nt1, nt2, .false., .false., .true.)  
        end if
        time(2) = mpi_wtime()
        comm_time = comm_time + time(2)-time(1)
     end if

     ! Why does this need to be set?
     rkStage = 0

     ! Compute the skin-friction velocity
     call computeUtau

     ! Compute time step
     call timestep(.false.)

     ! Possible Turblent Equations
     if( equations == RANSEquations ) then
        call initres(nt1MG, nMGVar) ! Initialize only the Turblent Variables
        call turbResidual
     endif

     ! Initialize Flow residuals
     call initres(1_intType, nwf)

     ! Actual Residual Calc
     call residual 

     if (setPETScVecs) then
        call setRVec(rVec)
     end if

  end do ! Iteration Loop

  time(4) = mpi_wtime()

  ! Do a nice reduction on the times:

  allocate(all_times(nProc))
  call MPI_Gather (comm_time,1,sumb_real,all_times,1,sumb_real,0,sumb_comm_world,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  

  total_time_local = time(4)-time(3)
  
  call  MPI_Reduce(total_time_local,total_time,1,sumb_real,MPI_MAX,sumb_comm_world,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (myid == 0) then
     do i=1,nProc
        print *,'Proc',i,' comm time:',all_times(i)
     end do
     print *,  ' '
     print *,'Total Time:',total_time
  end if




  deallocate(all_times)
end subroutine weak_scaling_test
