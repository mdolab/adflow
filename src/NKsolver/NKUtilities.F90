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
                 do l=1,nw
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
  real(kind=realType),pointer :: rvec_pointer(:)
  real(Kind=realType) :: ovv
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
                 ovv = 1/vol(i,j,k)
                 do l=1,nwf
                    rvec_pointer(ii) = dw(i,j,k,l)*ovv
                    ii = ii + 1
                 end do
                 do l=nt1,nt2
                    rvec_pointer(ii) = dw(i,j,k,l)*ovv*turbResScale
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

subroutine getInfoSize(iSize)
  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate
  use inputphysics

  implicit none
  integer(kind=intType), intent(out) :: iSize
  integer(kind=intType) :: nn, sps, nc
  ! Determine the size of a flat array needed to store w, P, ( and
  ! rlv, rev if necessary) with full double halos. 
  iSize = 0
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1_intType,sps)
        nc = (kb+1)*(jb+1)*(ib+1)
        iSize = iSize + nc*(nw + 1) ! plus 1 for the P
        if (viscous) then
           iSize = iSize + nc
        end if
        if (eddyModel) then
           iSize = iSize + nc
        end if
     end do
  end do
end subroutine getInfoSize

subroutine setInfo(info, iSize)

  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate
  use inputphysics
  implicit none

  real(kind=realType), intent(in), dimension(iSize) :: info
  integer(kind=intType), intent(in) :: iSize
  integer(kind=intType) :: nn, counter, i, j, k, l, sps
  ! Determine the size of a flat array needed to store w, P, ( and
  ! rlv, rev if necessary) with full double halos. 
  counter = 0
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1,sps)
        do k=0,kb
           do j=0,jb
              do i=0,ib
                 do l=1,nw
                    counter = counter + 1
                    w(i,j,k,l) = info(counter)
                 end do
                 
                 counter = counter + 1
                 P(i,j,k) = info(counter)
                 
                 if (viscous) then
                    counter = counter + 1
                    rlv(i,j,k) = info(counter)
                 end if

                 if (eddyModel) then
                    counter = counter + 1
                    rev(i,j,k) = info(counter)
                 end if
              end do
           end do
        end do
     end do
  end do
end subroutine setInfo

subroutine getInfo(info, iSize)
  use blockPointers
  use inputTimeSpectral
  use flowvarrefstate
  use inputphysics
  implicit none

  real(kind=realType), intent(out), dimension(iSize) :: info
  integer(kind=intType), intent(in) :: iSize
  integer(kind=intType) ::  nn, counter, i, j, k, l, sps
  ! Determine the size of a flat array needed to store w, P, ( and
  ! rlv, rev if necessary) with full double halos. 
  counter = 0
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1,sps)
        do k=0,kb
           do j=0,jb
              do i=0,ib
                 do l=1,nw
                    counter = counter + 1
                    info(counter) = w(i,j,k,l)
                 end do
                 
                 counter = counter + 1
                 info(counter) = P(i,j,k) 
                 
                 if (viscous) then
                    counter = counter + 1
                    info(counter) = rlv(i,j,k) 
                 end if

                 if (eddyModel) then
                    counter = counter + 1
                    info(counter) = rev(i,j,k)
                 end if
              end do
           end do
        end do
     end do
  end do
end subroutine getInfo

