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
#define PETSC_AVOID_MPIF_H
#include "finclude/petscdef.h"

  use communication
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration
  use NKsolvervars, only : times
  use petscvec

  implicit none

  Vec     wVec
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

end subroutine setWVec

subroutine setRVec(rVec)
  ! Set the current residual in dw into the PETSc Vector
#define PETSC_AVOID_MPIF_H
#include "finclude/petscdef.h"

  use communication
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration
  use NKsolvervars, only : times
  use petscvec

  implicit none

  Vec     rVec
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
                 ovv = 1/vol(i,j,k)
                 do l=1,nwf
                    rvec_pointer(ii) = dw(i,j,k,l)*ovv
                    ii = ii + 1
                 end do

                 do l=nt1,nt2
                    rvec_pointer(ii) = dw(i,j,k,l)*ovv
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
#define PETSC_AVOID_MPIF_H
#include "finclude/petscdef.h"

  ! Set the SUmb state vector, w, from the petsc vec wVec
  use communication
  use blockPointers
  use inputTimeSpectral
  use flowVarRefState
  use NKsolvervars, only : times
  use petscvec

  implicit none

  Vec     wVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType) :: temp,diff
  real(kind=realType),pointer :: wvec_pointer(:)

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
 
end subroutine setW

subroutine setW_ghost(wVec)
#define PETSC_AVOID_MPIF_H
#include "finclude/petscdef.h"

  ! Set the SUmb state vector, w, from the petsc vec wVec
  use communication
  use blockPointers
  use inputTimeSpectral
  use flowVarRefState
  use NKsolvervars, only : times
  use petscvec

  implicit none

  Vec     wVec
  Vec     wVec_l
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType) :: temp,diff
  real(kind=realType),pointer :: wvec_pointer(:)

  ! Note this is not ideal memory access but the values are stored by
  ! block in PETSc (grouped in nw) but are stored separately in SUmb

  call cpu_time(times(1))
  call VecGhostUpdateBegin(wVec,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecGhostUpdateEnd(wVec,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call cpu_time(times(2))

  times(10) = times(10) + times(2)-times(1)
  call VecGhostGetLocalForm(wVec,wVec_l,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecGetArrayF90(wVec_l,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ii = 1
  ! First do the "Owned" cells
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

  ! Next do the halos:
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,1_intType,sps)
        
        ! Loop over all 6 faces doubly extruded faces and add to
        ! list if necessary:
        
        ! I-Low Face
        do k=0,kb
           do j=0,jb
              do i=0,1
                 if (globalCell(i,j,k) .ge. 0) then
                    do l=1,nw
                       w(i,j,k,l) = wvec_pointer(ii)
                       ii = ii + 1
                    end do
                 end if
              end do
           end do
        end do
           
        ! I-High Face
        do k=0,kb
           do j=0,jb
              do i=ib-1,ib
                 if (globalCell(i,j,k) .ge. 0) then
                    do l=1,nw
                       w(i,j,k,l) = wvec_pointer(ii)
                       ii = ii + 1
                    end do
                 end if
              end do
           end do
        end do

        ! J-Low Face
        do k=0,kb
           do j=0,1
              do i=0,ib
                 if (globalCell(i,j,k) .ge. 0) then
                    do l=1,nw
                       w(i,j,k,l) = wvec_pointer(ii)
                       ii = ii + 1
                    end do
                 end if
              end do
           end do
        end do

        ! J-High Face
        do k=0,kb
           do j=jb-1,jb
              do i=0,ib
                 if (globalCell(i,j,k) .ge. 0) then
                    do l=1,nw
                       w(i,j,k,l) = wvec_pointer(ii)
                       ii = ii + 1
                    end do
                 end if
              end do
           end do
        end do

        ! K-Low Face
        do k=0,1
           do j=0,jb
              do i=0,ib
                 if (globalCell(i,j,k) .ge. 0) then
                    do l=1,nw
                       w(i,j,k,l) = wvec_pointer(ii)
                       ii = ii + 1
                    end do
                 end if
              end do
           end do
        end do

        ! K-High Face
        do k=kb-1,kb
           do j=0,jb
              do i=0,ib
                 if (globalCell(i,j,k) .ge. 0) then
                    do l=1,nw
                       w(i,j,k,l) = wvec_pointer(ii)
                       ii = ii + 1
                    end do
                 end if
              end do
           end do
        end do
     end do
  end do
  
  ! Restore pointer to local array
  call VecRestoreArrayF90(wVec_l,wvec_pointer,ierr)
  call EChk(ierr,__FILE__,__LINE__)
 
  ! Restore the ghosted form of the array
  call VecGhostRestoreLocalForm(wVec,wVec_l,ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine setW_ghost

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

subroutine setRVec_old(rVec)
  ! Set the current residual in dw into the PETSc Vector

  use communication
  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputiteration
  use NKsolvervars, only : times
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  Vec     rVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l
  real(kind=realType) :: ovv,temp(nw)
  
  times(1) = mpi_wtime()

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
                    temp(l) = dw(i,j,k,l)*ovv!*1e4
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

  times(2) = mpi_wtime()
  times(10) = times(10) + times(2)-times(1)

end subroutine setRVec_old

subroutine setW_old(wVec)

  ! Set the SUmb state vector, w, from the petsc vec wVec
  use communication
  use blockPointers
  use inputTimeSpectral
  use flowVarRefState
  use NKsolvervars, only : times
  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  Vec     wVec
  integer(kind=intType) :: ierr,nn,sps,i,j,k,l,ii
  real(kind=realType) :: temp,diff

  times(1) = mpi_wtime()
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

  times(2) = mpi_wtime()
  times(20) = times(20) + times(2)-times(1)

end subroutine setW_old
