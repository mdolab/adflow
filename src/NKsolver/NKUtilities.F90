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
  ! Set the current residual in dw into the PETSc Vector
#ifndef USE_NO_PETSC
#define PETSC_AVOID_MPIF_
#include "finclude/petscdef.h"

  use communication
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration

  implicit none

  PetscInt     wVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType),pointer :: wvec_pointer(:)

  call VecGetArrayF90(wVec,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ii = 1
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1_intType,sps)
        ! Copy off w to wVec
        do k=2,kl
           do j=2,jl
              do i=2,il
                 do l=1,nwf
                    wvec_pointer(ii) = w(i,j,k,l)
                    ii = ii + 1
                 end do
              end do
           end do
        end do
     end do
  end do

  call VecRestoreArrayF90(wVec,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
end subroutine setWVec

subroutine setRVec(rVec)
#ifndef USE_NO_PETSC
  ! Set the current residual in dw into the PETSc Vector
#define PETSC_AVOID_MPIF_H
#include "finclude/petscdef.h"

  use communication
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration
  use NKsolvervars, only: scalevec
  implicit none

  PetscInt     rVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType) :: ovv
  real(kind=realType),pointer :: rvec_pointer(:)

  call VecGetArrayF90(rVec,rvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ii = 1
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1_intType,sps)
        ! Copy off dw/vol to rVec
        do k=2,kl
           do j=2,jl
              do i=2,il
                 !ovv = 1/vol(i,j,k)
                 do l=1,nw
                    rvec_pointer(ii) = dw(i,j,k,l)!*ovv
                    ii = ii + 1
                 end do
              end do
           end do
        end do
     end do
  end do
  
  call VecRestoreArrayF90(rVec,rvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ! VecPointwiseMult(Vec w, Vec x,Vec y) w = x *y
  call VecPointwiseMult(rVec, rVec, scaleVec, ierr)
  call EChk(ierr,__FILE__,__LINE__)

#endif
end subroutine setRVec

subroutine setW(wVec)
#ifndef USE_NO_PETSC
#define PETSC_AVOID_MPIF_H
#include "finclude/petscdef.h"

  ! Set the SUmb state vector, w, from the petsc vec wVec
  use communication
  use blockPointers
  use inputTimeSpectral
  use flowVarRefState

  implicit none

  PetscInt     wVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType) :: temp,diff
  real(kind=realType),pointer :: wvec_pointer(:)
  logical :: commPressure, commGamma, commViscous
  ! Note this is not ideal memory access but the values are stored by
  ! block in PETSc (grouped in nw) but are stored separately in SUmb

  call VecGetArrayF90(wVec,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ii = 1
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1_intType,sps)

        do k=2,kl
           do j=2,jl
              do i=2,il
                 do l=1,nw
                    w(i,j,k,l) = wvec_pointer(ii) 
                    ii = ii + 1
                 end do
              end do
           end do
        end do
     end do
  end do

  call VecRestoreArrayF90(wVec,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
 
  ! Run the double halo exchange:
!   commPressure = .False.
!   commGamma = .False.
!   commViscous = .True. 
!   call whalo2(1_intType, 1_intType, nw, commPressure, commGamma, commViscous)
#endif
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
  real(kind=realType),dimension(ndimw),intent(inout) :: res(ndimw)

  ! Local Variables
  integer(kind=intType) :: nn,i,j,k,l,counter,sps

  call whalo2(1_intType, 1_intType, nw, .False., &
       .False.,.False.)

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


subroutine calcScaling(scaleVec)

#ifndef USE_NO_PETSC


  use communication 
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration
  use NKSolverVars, only : resSum
  implicit none
#define PETSC_AVOID_MPIF_H
#include "finclude/petscdef.h"
#include "include/finclude/petsc.h"

  PetscInt     scaleVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType) :: ovv
  real(kind=realType),pointer :: scale_pointer(:)
  real(kind=realType) :: resSum_l(nw)
  real(kind=realTYpe) :: norm


  ! Loop over current residual and determine the scaling for each of
  ! the nw equations:
  
  resSum(:) = zero
  resSum_l(:) = zero
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1_intType,sps)
        do l=1,nw
           do k=2,kl
              do j=2,jl
                 do i=2,il
                    resSum_l(l) = resSum_l(l) + (dw(i,j,k,l)/vol(i,j,k))**2
                 end do
              end do
           end do
        end do
     end do
  end do

  do l=1,nw
     call mpi_allreduce(resSum_l(l),resSum(l),1,sumb_real,mpi_sum,&
          SUmb_comm_world, ierr)
  end do

!   if (myid == 0) then
!      print *,'resSum after reduce:',sqrt(resSum(1:5)/nCellGlobal(1))
!   end if

  ! determine scaline wrt resSum(1) (density residual)

!   do l=2,nw
!      resSum(l) = sqrt(resSum(l)/resSum(1))
!   end do
!   resSum(1) = one

!   do l=1,nw
!      resSum(l) = min(resSum(l),10.0)
!      resSum(l) = max(resSum(l),0.1)
!   end do

!   if (myid == 0) then
!      print *,'resSUm:',resSum(1:3)
!   end if
!   do l=1,nw
!      resSum(l) = one / resSum(l)
!   end do

  resSum(:) = one
  call VecGetArrayF90(scaleVec,scale_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ii = 1
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1_intType,sps)
        ! Copy off dw/vol to rVec
        do k=2,kl
           do j=2,jl
              do i=2,il
                 ovv = 1/vol(i,j,k)
                 do l=1,nw
                    scale_pointer(ii) = ovv * resSum(l)
                    ii = ii + 1
                 end do
              end do
           end do
        end do
     end do
  end do
  
  call VecRestoreArrayF90(scaleVec,scale_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecNorm(scaleVec, NORM_2,norm,ierr)


#endif
end subroutine calcScaling
