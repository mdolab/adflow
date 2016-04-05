! This files contains several utilitiy functions that are used with
! the ANK solver.

subroutine setWVecANK(wVec)
  ! Set the current FLOW variables in the PETSc Vector
#ifndef USE_NO_PETSC
  use communication
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration

  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#endif

  Vec   wVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType),pointer :: wvec_pointer(:)

  call VecGetArrayF90(wVec,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ii = 0
  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, 1_intType, sps)
        ! Copy off w to wVec
        do k=2, kl
           do j=2, jl
              do i=2, il
                 do l=1, nwf
                    ii = ii + 1
                    wvec_pointer(ii) = w(i, j, k, l)
                 end do
              end do
           end do
        end do
     end do
  end do

  call VecRestoreArrayF90(wVec, wvec_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif

end subroutine setWVecANK

subroutine setRVecANK(rVec)
#ifndef USE_NO_PETSC

  ! Set the current FLOW residual in dw into the PETSc Vector
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration
  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#endif

  Vec    rVec
  integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
  real(kind=realType),pointer :: rvec_pointer(:)
  real(Kind=realType) :: ovv
  call VecGetArrayF90(rVec,rvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ii = 0
  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn,1_intType,sps)
        ! Copy off dw/vol to rVec
        do k=2, kl
           do j=2, jl
              do i=2, il
                 ovv = one/vol(i,j,k)
                 do l=1, nwf
                    ii = ii + 1        
                    rvec_pointer(ii) = dw(i, j, k, l)*ovv
                 end do
              end do
           end do
        end do
     end do
  end do
  
  call VecRestoreArrayF90(rVec, rvec_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
end subroutine setRVecANK

subroutine setWANK(wVec)
#ifndef USE_NO_PETSC

  use communication
  use blockPointers
  use inputTimeSpectral
  use flowVarRefState

  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#endif

  Vec  wVec
  integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
  real(kind=realType), pointer :: wvec_pointer(:)

#if PETSC_VERSION_MINOR > 5
  call VecGetArrayReadF90(wVec, wvec_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#else
  call VecGetArrayF90(wVec, wvec_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
  ii = 0
  do nn=1, nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn, 1_intType, sps)

        do k=2, kl
           do j=2, jl
              do i=2, il
                 do l=1, nwf
                    ii = ii + 1
                    w(i, j, k, l) = wvec_pointer(ii) 
                 end do
              end do
           end do
        end do
     end do
  end do
#if PETSC_VERSION_MINOR > 5
  call VecRestoreArrayReadF90(wVec, wvec_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#else
  call VecRestoreArrayF90(wVec, wvec_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
#endif
end subroutine setWANK


subroutine setWVecANKTurb(wVec)
  ! Set the current Turbulence variables in the PETSc Vector
#ifndef USE_NO_PETSC
  use communication
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration

  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#endif

  Vec   wVec
  integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
  real(kind=realType),pointer :: wvec_pointer(:)

  call VecGetArrayF90(wVec, wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ii = 0
  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, 1_intType, sps)
        ! Copy off w to wVec
        do k=2, kl
           do j=2, jl
              do i=2, il
                 do l=nt1, nt2
                    ii = ii + 1
                    wvec_pointer(ii) = w(i, j, k, l)
                 end do
              end do
           end do
        end do
     end do
  end do

  call VecRestoreArrayF90(wVec, wvec_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif

end subroutine setWVecANKTurb

subroutine setRVecANKTurb(rVec)
#ifndef USE_NO_PETSC

  ! Set the current FLOW residual in dw into the PETSc Vector
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration
  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#endif

  Vec    rVec
  integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
  real(kind=realType),pointer :: rvec_pointer(:)
  real(Kind=realType) :: ovv
  call VecGetArrayF90(rVec,rvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  ii = 0
  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, 1_intType, sps)
        ! Copy off dw/vol to rVec
        do k=2, kl
           do j=2, jl
              do i=2, il
                 do l=nt1, nt2
                    ii = ii + 1
                    rvec_pointer(ii) = dw(i,j,k,l)*turbResScale(l-nt1+1)/vol(i, j, k)
                 end do
              end do
           end do
        end do
     end do
  end do
  
  call VecRestoreArrayF90(rVec, rvec_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
end subroutine setRVecANKTurb

subroutine setWANKTurb(wVec)
#ifndef USE_NO_PETSC

  use communication
  use blockPointers
  use inputTimeSpectral
  use flowVarRefState

  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#endif

  Vec  wVec
  integer(kind=intType) :: ierr, nn, sps, i, j, k, l, ii
  real(kind=realType), pointer :: wvec_pointer(:)

#if PETSC_VERSION_MINOR > 5
  call VecGetArrayReadF90(wVec, wvec_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#else
  call VecGetArrayF90(wVec, wvec_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
  ii = 0
  do nn=1, nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn, 1_intType, sps)

        do k=2, kl
           do j=2, jl
              do i=2, il
                 do l=nt1, nt2
                    ii = ii + 1
                    w(i, j, k, l) = wvec_pointer(ii) 
                 end do
              end do
           end do
        end do
     end do
  end do
#if PETSC_VERSION_MINOR > 5
  call VecRestoreArrayReadF90(wVec, wvec_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#else
  call VecRestoreArrayF90(wVec, wvec_pointer, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
#endif
end subroutine setWANKTurb


