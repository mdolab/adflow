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
