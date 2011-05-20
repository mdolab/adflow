subroutine LSCheck(snes,VecX,VecY,lsctx,changed_Y,ierr)
  ! ---------------------------------------------------------------------
  !
  !  FormPreCheck - checks the validity of a new direction given by the linear solve 
  !                 before the line search is called
  !
  !  Input Parameters:
  !  snes  - the SNES context
  !  X     - previous iterate
  !  Y     - New search direction 
  !  lsctx - line search context (empty)
  !  changed_y - Flag to indicated if the direction was changed

  use precision
  use flowVarRefState
  use inputTimeSpectral
  use blockPointers
  use communication
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! PETSc Variables

  SNES           snes
  Vec            VecX,VecY
  Vec            VecXTemp
  PetscFortranAddr   lsctx
  logical :: changed_Y
  !PetscErrorCode ierr
  integer(kind=intType) :: ierr
  integer(kind=intType) :: foundNan
  integer(kind=intType) :: iter
  integer(kind=intType) :: i,j,k,l,nn,sps
  real(kind=realType) :: value
  print *,'Inside formPreCheck,ierr is:',ierr
  ierr = 0
  ! Basically what we need to know is does Residual(X+Y) give us a nan?

  ! First check to see if petsc is fucking with us, and run it on the just VecX
!   call setW(vecX)
     
!   print *,'Computing Residual on x vec'
!   call computeResidualNK2()

  ! First generate the actual new iterate:
  call VecDuplicate(vecX,VecXTemp,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecCopy(vecX,VecXTemp,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecAXPY(vecXTemp,1.0,vecY,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  changed_y = .False.
  backTrack: do iter=1,10
     ! Set the vector and compute the residual
     print *,'---------------- Iter ---------------'
     print *,iter
     
     print *,'Setting w'
     call setW(vecXTemp)
     
     print *,'Computing Residual'
     call computeResidualNK2()

     ! Check for Nan's:
     print *,'Checking for nan'
     call checkdwForNan(foundNan)
     print *,'Found Nan is:',foundNan
     if (foundNan .ne. 0) then
        ! We have a found a nan, so we have to keep looping by halving
        ! the y value until we find a point that works
        ! We will scale the value of Y by 0.5, then SUBTRACT it from
        ! X+Y. The series of iterates is therefore
        ! X+Y, X+Y-0.5Y, X+Y-0.5Y-0.25Y, X+Y-0.5Y-0.25Y-.125Y... = 
        ! X+Y, X+0.5Y, X+0.25Y, X+0.125Y ...
        ! Which is precisely our 0.5 backtrack we want
        !print *,'scaling y'
        
        call VecScale(vecY,0.5,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call VecNorm(vecY,NORM_2,value,ierr)

        print *,'Y VecNorm:',value

        call EChk(ierr,__FILE__,__LINE__)
        call VecAXPY(vecXTemp,-1.0,vecY,ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! Try to get rid of fucking nan:
        call setUniformFlow()
        ! Set dw to zero 
        do sps=1,nTimeIntervalsSpectral
           do nn=1,nDom
              call setPointers(nn,1,sps)
              dw = 0.0
              fw = 0.0
           end do
        end do
        call EChk(ierr,__FILE__,__LINE__)
        changed_y = .True.
     else
        ! Its ok, so just break loop
        exit backTrack
     end if
  end do backTrack
  print *,'ierr,ichanged_y:',ierr,changed_y

  ! Reset the correct 'x' vector back, since weve been screwing with it
  call setUniformFlow()
  ! Set dw to zero 
  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointers(nn,1,sps)
        dw = 0.0
        fw = 0.0
     end do
  end do
  call setW(VecX)

  ! Destroy the temporary vector
  call VecDestroy(VecXTemp,ierr)
  call EChk(ierr,__FILE__,__LINE__)


end subroutine LSCheck

subroutine checkdwForNan(foundNan)
  use precision
  use flowVarRefState
  use inputTimeSpectral
  use blockPointers
  use communication
  implicit none

  integer(kind=intType) :: foundNan
  integer(kind=intType) :: i,j,k,l,nn,sps
  real(kind=realType) :: sum
  logical :: local_nan_logical
  integer(kind=intType) :: local_nan
  integer(kind=intType) :: ierr
  logical :: myIsNAN
  sum = 0.0
  foundNan = 0
  ! We want a 'fast' way to check for nan...blindly sum all residuals in dw:
  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointers(nn,1,sps)
        do l=1,nw
           do k=2,kl
              do j=2,jl
                 do i=2,il
                    sum = sum + dw(i,j,k,l)
                 end do ! iloop
              end do ! jloop
           end do ! kloop
        end do ! state loop
     end do ! domain loop
  end do ! sps loop

  ! Check the nan on this processor
  local_nan_logical = myIsNAN(sum)
  if (local_nan_logical) then
     local_nan = 1
  else
     local_nan= 0
  end if
  print *,'dw local Nan:',myid,local_nan
  ! All reduce with a logical or:
  call mpi_allreduce(local_nan, foundNan, 1, sumb_integer, &
       mpi_sum, SUmb_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)  
  print *,'dw global Nan:',foundNan
end subroutine checkdwForNan

subroutine checkwForNan(foundNan)
  use precision
  use flowVarRefState
  use inputTimeSpectral
  use blockPointers
  use communication
  implicit none

  integer(kind=intType) :: foundNan
  integer(kind=intType) :: i,j,k,l,nn,sps
  real(kind=realType) :: sum
  logical :: local_nan_logical
  integer(kind=intType) :: local_nan
  integer(kind=intType) :: ierr
  logical :: myIsNAN
  sum = 0.0
  foundNan = 0
  ! We want a 'fast' way to check for nan...blindly sum all residuals in dw:
  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointers(nn,1,sps)
        do l=1,nw
           do k=2,kl
              do j=2,jl
                 do i=2,il
                    sum = sum + w(i,j,k,l)
                 end do ! iloop
              end do ! jloop
           end do ! kloop
        end do ! state loop
     end do ! domain loop
  end do ! sps loop

  ! Check the nan on this processor
  local_nan_logical = myIsNAN(sum)
  if (local_nan_logical) then
     local_nan = 1
  else
     local_nan= 0
  end if
  print *,'w local Nan:',myid,local_nan
  ! All reduce with a logical or:
  call mpi_allreduce(local_nan, foundNan, 1, sumb_integer, &
       mpi_sum, SUmb_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)      
  print *,'w global Nan:',foundNan
end subroutine checkwForNan

subroutine checkPForNan(foundNan)
  use precision
  use flowVarRefState
  use inputTimeSpectral
  use blockPointers
  use communication
  implicit none

  integer(kind=intType) :: foundNan
  integer(kind=intType) :: i,j,k,l,nn,sps
  real(kind=realType) :: sum
  logical :: local_nan_logical
  integer(kind=intType) :: local_nan
  integer(kind=intType) :: ierr
  logical :: myIsNAN
  sum = 0.0
  foundNan = 0
  ! We want a 'fast' way to check for nan...blindly sum all residuals in dw:
  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointers(nn,1,sps)
        do k=2,kl
           do j=2,jl
              do i=2,il
                 sum = sum + p(i,j,k)
              end do ! iloop
           end do ! jloop
        end do ! kloop
     end do ! domain loop
  end do ! sps loop
  
  ! Check the nan on this processor
  local_nan_logical = myIsNAN(sum)
  if (local_nan_logical) then
     local_nan = 1
  else
     local_nan= 0
  end if

  print *,'P Local Nan:',local_nan
  ! All reduce with a logical or:
  call mpi_allreduce(local_nan, foundNan, 1, sumb_integer, &
       mpi_sum, SUmb_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)  
  print *,'P GLobal Nan:',foundNan
end subroutine checkPForNan
