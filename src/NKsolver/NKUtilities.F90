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
  use communication
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration

  implicit none
#define PETSC_AVOID_MPIF_H
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  Vec   wVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType),pointer :: wvec_pointer(:)

  call VecGetArrayF90(wVec,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ii = 1
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1_intType,sps)
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

subroutine setWVec2(wVec)
  ! Set the current residual in dw into the PETSc Vector
#ifndef USE_NO_PETSC
  use communication
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration

  implicit none
#define PETSC_AVOID_MPIF_H
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  Vec   wVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType),pointer :: wvec_pointer(:)

  call VecGetArrayF90(wVec,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ii = 1
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1_intType,sps)
        ! Copy off w to wVec
        do l=1,nwf
           do k=2,kl
              do j=2,jl
                 do i=2,il
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

end subroutine setWVec2

subroutine setRVec(rVec)
#ifndef USE_NO_PETSC
  ! Set the current residual in dw into the PETSc Vector
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration
  use NKsolvervars, only: scalevec
  implicit none
#define PETSC_AVOID_MPIF_H
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  Vec    rVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType) :: ovv
  real(kind=realType),pointer :: rvec_pointer(:)

  call VecGetArrayF90(rVec,rvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ii = 1
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1_intType,sps)
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

subroutine setRVec2(rVec)
#ifndef USE_NO_PETSC
  ! Set the current residual in dw into the PETSc Vector
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration
  use NKsolvervars, only: scalevec
  implicit none
#define PETSC_AVOID_MPIF_H
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  Vec    rVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType) :: ovv
  real(kind=realType),pointer :: rvec_pointer(:)

  call VecGetArrayF90(rVec,rvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ii = 1
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1_intType,sps)
        ! Copy off dw/vol to rVec
        do l=1,nw
           do k=2,kl
              do j=2,jl
                 do i=2,il
                    rvec_pointer(ii) = dw(i,j,k,l)
                    ii = ii + 1
                 end do
              end do
           end do
        end do
     end do
  end do
  
  call VecRestoreArrayF90(rVec,rvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
end subroutine setRVec2

subroutine setW(wVec)
#ifndef USE_NO_PETSC

  use communication
  use blockPointers
  use inputTimeSpectral
  use flowVarRefState

  implicit none
#define PETSC_AVOID_MPIF_H
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  Vec  wVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType) :: temp,diff
  real(kind=realType),pointer :: wvec_pointer(:)

  call VecGetArrayF90(wVec,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ii = 1
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1_intType,sps)

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
#endif
end subroutine setW

subroutine setW2(wVec)
#ifndef USE_NO_PETSC

  use communication
  use blockPointers
  use inputTimeSpectral
  use flowVarRefState

  implicit none
#define PETSC_AVOID_MPIF_H
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  Vec  wVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType) :: temp,diff
  real(kind=realType),pointer :: wvec_pointer(:)

  call VecGetArrayF90(wVec,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ii = 1
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1_intType,sps)
        do l=1,nw
           do k=2,kl
              do j=2,jl
                 do i=2,il
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
#endif
end subroutine setW2


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
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  Vec   scaleVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType) :: ovv
  real(kind=realType) :: resSum_l(nw)
  real(kind=realTYpe) :: norm
  real(kind=realType),pointer :: scale_pointer(:)

  ! Loop over current residual and determine the scaling for each of
  ! the nw equations:
  
  resSum(:) = zero
  resSum_l(:) = zero
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1_intType,sps)
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
        call setPointers(nn,1_intType,sps)
        ! Copy off dw/vol to rVec
        do k=2,kl
           do j=2,jl
              do i=2,il
                 ovv = 1/vol(i,j,k)
                 do l=1,nw
                    scale_pointer(ii) = ovv !* resSum(l)
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

subroutine  MyMult(matrix, X, F, ierr)

#ifndef USE_NO_PETSC

  use constants
  use communication
  use NKSolverVars, only: wBase, rBase, wVec, diag

  implicit none
#define PETSC_AVOID_MPIF_H
#include "finclude/petsc.h"
  
  Mat   matrix
  Vec   X, F
  
  real(kind=realType) :: h, sum, nrm, dot, umin, err_rel, value
  integer(kind=intType) :: ierr
  umin = 1e-6
  err_rel = 1e-8
  ! Get the step size we should use:
  call VecDotBegin(wBase,X,dot,ierr)
  call VecNormBegin(X,NORM_1,sum,ierr)
  call VecNormBegin(X,NORM_2,nrm,ierr)
  call VecDotEnd(wBase,X,dot,ierr)
  call VecNormEnd(X,NORM_1,sum,ierr)
  call VecNormEnd(X,NORM_2,nrm,ierr)

  !   Safeguard for step sizes that are "too small"
  if (dot < umin*sum .and. dot >= zero) then
     dot = umin*sum
  else if (dot < 0.0 .and. dot > -umin*sum) then
     dot = -umin*sum;
  end if

  ! Step is
   h = err_rel*dot/(nrm*nrm)
!   if (myid == 0) then
!      print *,'h is:',h
!   end if
   !h = 1e-5
  ! Compute the peturbed wVec =  h*X + wBase
  call VecWAXPY(wVec, h, X, wBase, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ! Evalute Residual
  call setW(wVec)
  call computeResidualNK()
  call setRVec(F)

  ! Compute the difference: F = F -one*rVec
  call VecAXPY(F, -one, rBase, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Scale by 1/h
  call VecScale(F, one/h, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif

!   call setdtl(diag)
!  call VecPointwiseMult(diag, diag, X, ierr)
!  call VecAXPY(F, one, diag, ierr) !

end subroutine MyMult

subroutine setdtl(D)
#ifndef USE_NO_PETSC
  use communication
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration
  
  implicit none
#define PETSC_AVOID_MPIF_H
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  Vec  D
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType),pointer :: d_pointer(:)

  call VecGetArrayF90(D,d_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ii = 1
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1_intType,sps)
        ! Copy off w to wVec
        do k=2,kl
           do j=2,jl
              do i=2,il
                 do l=1,nwf
                    D_pointer(ii) = 1/(5000000*dtl(i,j,k))
                    ii = ii + 1
                 end do
              end do
           end do
        end do
     end do
  end do

  call VecRestoreArrayF90(D,D_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
end subroutine setdtl

subroutine setBase(w)

#ifndef USE_NO_PETSC

  use constants
  use NKSolverVars, only: wBase, rBase

  implicit none
#define PETSC_AVOID_MPIF_H
#include "finclude/petsc.h"
  
  Vec w
  integer(kind=intType) :: ierr

  ! Copy to wBase 
  call VecCopy(w, wBase, ierr)

  ! Evalue residual at wBase to get rBase:
   call setW(wBase)
  call computeResidualNK()
  call setRVec(rBase)
#endif

end subroutine setBase

