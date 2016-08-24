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
  use constants
  use blockPointers, only : nDom, il, jl, kl, w
  use inputtimespectral, only : ntimeIntervalsSpectral
  use flowvarrefstate, only : nw
  use utils, only : setPointers, EChk
  implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

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

end subroutine setWVec

subroutine setRVec(rVec)

  ! Set the current residual in dw into the PETSc Vector
  use constants
  use blockPointers, only : nDom, volRef, il, jl, kl, dw
  use inputtimespectral, only : nTimeIntervalsSpectral
  use flowvarrefstate, only : nw, nwf, nt1, nt2
  use inputIteration, only : turbResScale
  use utils, only : setPointers, EChk
  implicit none

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

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
        do k=2, kl
           do j=2, jl
              do i=2, il
                 ovv = 1/volRef(i, j, k)
                 do l=1,nwf
                    rvec_pointer(ii) = dw(i, j, k, l)*ovv
                    ii = ii + 1
                 end do
                 do l=nt1,nt2
                    rvec_pointer(ii) = dw(i, j, k, l)*ovv*turbResScale(l-nt1+1)
                    ii = ii + 1
                 end do
              end do
           end do
        end do
     end do
  end do
  
  call VecRestoreArrayF90(rVec,rvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine setRVec

subroutine setW(wVec)

  use constants
  use blockPointers, only : nDom, il, jl, kl, w
  use inputTimeSpectral, only : nTimeIntervalsSpectral
  use flowVarRefState, only : nw
  use utils, only : setPointers, EChk

  implicit none
#define PETSC_AVOID_MPIF_H
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  Vec  wVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType),pointer :: wvec_pointer(:)

#if PETSC_VERSION_MINOR > 5
  call VecGetArrayReadF90(wVec,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#else
  call VecGetArrayF90(wVec,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
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
#if PETSC_VERSION_MINOR > 5
  call VecRestoreArrayReadF90(wVec,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#else
  call VecRestoreArrayF90(wVec,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif

end subroutine setW

subroutine getStates(states,ndimw)
  ! Return the state vector, w to Python

  use constants
  use blockPointers, only : il, jl, kl, nDom, w
  use inputTimeSpectral, only : nTimeIntervalsSpectral
  use flowvarrefstate, only : nw
  use utils, only : setPointers

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
  use constants
  use blockPointers, only : il, jl, kl, nDom, dw, volRef
  use inputTimeSpectral, only : nTimeIntervalsSpectral
  use flowvarrefstate, only : nw
  use utils, only : setPointers

  implicit none
 
  integer(kind=intType),intent(in):: ndimw
  real(kind=realType),dimension(ndimw),intent(inout) :: res(ndimw)

  ! Local Variables
  integer(kind=intType) :: nn,i,j,k,l,counter,sps
  real(kind=realType) :: ovv

  call computeResidualNK()
  counter = 0 
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,1,sps)
        do k=2,kl
           do j=2,jl
              do i=2,il
                 ovv = one/volRef(i,j,k)
                 do l=1,nw
                    counter = counter + 1
                    res(counter) = dw(i,j,k,l)*ovv
                 end do
              end do
           end do
        end do
     end do
  end do

end subroutine getRes

subroutine setStates(states,ndimw)

  ! Take in externallly generated states and set them in SUmb
  use constants
  use blockPointers, only : il, jl, kl, nDom, w
  use inputTimeSpectral, only : nTimeIntervalsSpectral
  use flowvarrefstate, only : nw
  use utils, only : setPointers

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
  use constants
  use blockPointers, only : ib, jb, kb, nDom
  use inputTimeSpectral, only : nTimeIntervalsSpectral
  use flowvarrefstate, only : nw, viscous, eddymodel
  use utils, only : setPointers

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

  use constants
  use blockPointers, only : w, p, ib, jb, kb, rlv, rev, nDom
  use inputTimeSpectral, only : nTimeIntervalsSpectral
  use flowvarrefstate, only : nw, viscous, eddymodel
  use utils, only : setPointers
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

  use constants
  use blockPointers, only : w, p, ib, jb, kb, rlv, rev, nDom
  use inputTimeSpectral, only : nTimeIntervalsSpectral
  use flowvarrefstate, only : nw, viscous, eddymodel
  use utils, only : setPointers

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

subroutine getEWTol(norm, old_norm, rtol_last, rtol)

  use constants
  implicit none

  ! There are the default EW Parameters from PETSc. They seem to work well
  !version:           2
  !rtol_0:  0.300000000000000     
  !rtol_max:  0.900000000000000     
  !gamma:   1.00000000000000     
  !alpha:   1.61803398874989     
  !alpha2:   1.61803398874989     
  !threshold:  0.100000000000000     

  real(kind=realType), intent(in) :: norm, old_norm, rtol_last
  real(kind=realType), intent(out) :: rtol
  real(kind=realType) :: rtol_max, gamma, alpha, alpha2, threshold, stol

  rtol_max  = 0.5_realType
  gamma     = 1.0_realType
  alpha     = (1.0_realType+sqrt(five))/2.0_realType
  alpha2    = (1.0_realType+sqrt(five))/2.0_realType
  threshold = 0.10_realType
  ! We use version 2:
  rtol = gamma*(norm/old_norm)**alpha
  stol = gamma*rtol_last**alpha
  
  if (stol > threshold) then
     rtol = max(rtol, stol)
  end if
  
  ! Safeguard: avoid rtol greater than one
  rtol = min(rtol, rtol_max)
  
end subroutine getEWTol

