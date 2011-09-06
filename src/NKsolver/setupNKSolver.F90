subroutine setupNKsolver

  ! Setup the PETSc objects for the Newton-Krylov
  ! solver. destroyNKsolver can be used to destroy the objects created
  ! in this function
  use blockPointers
  use communication
  use inputTimeSpectral
  use flowVarRefState
  use iteration
  use inputPhysics
  use stencils
  use ADjointVars , only: nCellsLocal
  use NKSolverVars, only: dRdw,dRdwPre,dRdwPseudo, ctx, wVec,rVec,deltaW,&
       NKsolvedOnce,nksolversetup,ksp_subspace,ksp_solver_type,global_ksp

  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Working Variables
  integer(kind=intType) :: ierr,nDimw,totalCells
  integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  integer(kind=intType) :: n_stencil
  integer(kind=intType), dimension(:,:), allocatable :: stencil
  integer(kind=intType),dimension(:),allocatable :: ghostInd
  integer(kind=intType) :: nghost

  integer(kind=intType) :: i,j,k,nn,i2,j2,k2,d2,l,sps
  external FormFunction_mf

  nGhost = 0
  if (not(NKSolverSetup)) then
     nDimW = nw * nCellsLocal * nTimeIntervalsSpectral

!      allocate(ghostInd(size(recvBuffer)),stat=ierr)
 
!      ! W is going to be a ghosted vector. This means that it will take
!      ! care of the communication of the halo data instead of using
!      ! whalo2.  We already have all the required information for the
!      ! transfers in the communication object. 
     
!      ! Loop over all procs second-level halos will receive from:
!      do nn=1,nDom
!         do sps=1,nTimeIntervalsSpectral
!            call setPointersAdj(nn,1_intType,sps)
           
!            ! Loop over all 6 faces doubly extruded faces and add to
!            ! list if necessary:
           
!            ! I-Low Face
!            do k=0,kb
!               do j=0,jb
!                  do i=0,1
!                     if (globalCell(i,j,k) .ge. 0) then
!                        nGhost = nGhost + 1
!                        ghostInd(nGhost) = globalCell(i,j,k)*nw
!                     end if
!                  end do
!               end do
!            end do
           
!            ! I-High Face
!            do k=0,kb
!               do j=0,jb
!                  do i=ib-1,ib
!                     if (globalCell(i,j,k) .ge. 0) then
!                        nGhost = nGhost + 1
!                        ghostInd(nGhost) = globalCell(i,j,k)*nw
!                     end if
!                  end do
!               end do
!            end do

!            ! J-Low Face
!            do k=0,kb
!               do j=0,1
!                  do i=0,ib
!                     if (globalCell(i,j,k) .ge. 0) then
!                        nGhost = nGhost + 1
!                        ghostInd(nGhost) = globalCell(i,j,k)*nw 
!                     end if
!                  end do
!               end do
!            end do

!            ! J-High Face
!            do k=0,kb
!               do j=jb-1,jb
!                  do i=0,ib
!                     if (globalCell(i,j,k) .ge. 0) then
!                        nGhost = nGhost + 1
!                        ghostInd(nGhost) = globalCell(i,j,k)*nw
!                     end if
!                  end do
!               end do
!            end do

!            ! K-Low Face
!            do k=0,1
!               do j=0,jb
!                  do i=0,ib
!                     if (globalCell(i,j,k) .ge. 0) then
!                        nGhost = nGhost + 1
!                        ghostInd(nGhost) = globalCell(i,j,k)*nw                    
!                     end if
!                  end do
!               end do
!            end do

!            ! K-High Face
!            do k=kb-1,kb
!               do j=0,jb
!                  do i=0,ib
!                     if (globalCell(i,j,k) .ge. 0) then
!                        nGhost = nGhost + 1
!                        ghostInd(nGhost) = globalCell(i,j,k)*nw
!                    end if
!                  end do
!               end do
!            end do
!         end do
!      end do

! !      do sps=1,nTimeIntervalsSpectral
! !         completeRecvs: do i=1,commPatternCell_2nd(1_intType)%nProcRecv
           
! !            do j=1,commPatternCell_2nd(1_intType)%nrecv(i)
              
! !               d2 = commPatternCell_2nd(1_intType)%recvList(i)%block(j)
! !               i2 = commPatternCell_2nd(1_intType)%recvList(i)%indices(j,1)
! !               j2 = commPatternCell_2nd(1_intType)%recvList(i)%indices(j,2)
! !               k2 = commPatternCell_2nd(1_intType)%recvList(i)%indices(j,3)

! !               nGhost = nGhost + 1
! !               ghostInd(nGhost) = flowDoms(d2,1_intType,sps)%globalCell(i2,j2,k2)
             

! !               if (myid == 0) then
! !                  !print *,d2,flowDoms(d2,1,sps)%il,flowDoms(d2,1,sps)%jl,flowDoms(d2,1,sps)%kl
! !                  print *,i2,j2,k2,ghostInd(nGhost)
! !               end if
! !            end do
           
! !         end do completeRecvs
! !      end do

! !      if (myid == 0) then
! !         do i=1,nGhost
! !            print *,ghostInd(i)
! !         end do
! !      end if
     

!      call VecCreateGhostBlock(SUMB_PETSC_COMM_WORLD,nw,nDimw,PETSC_DECIDE,&
!           nghost,ghostInd,wVec,ierr)
!      call EChk(ierr,__FILE__,__LINE__)

     
!      !call VecDestroy(wVec,ierr)
!      !call EChk(ierr,__FILE__,__LINE__)
!      deallocate(ghostInd)

     call VecCreateMPI(SUMB_PETSC_COMM_WORLD,nDimw,PETSC_DETERMINE,wVec,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call VecSetBlockSize(wVec,nw,ierr)
     call EChk(ierr,__FILE__,__LINE__)


     !  Create residual and state vectors
     call VecCreateMPI(SUMB_PETSC_COMM_WORLD,nDimw,PETSC_DETERMINE,rVec,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call VecSetBlockSize(rVec,nw,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Use the rVec Template to create deltaW 
     call VecDuplicate(rVec, deltaW, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Create Pre-Conditioning Matrix
     totalCells = nCellsLocal*nTimeIntervalsSpectral
     allocate( nnzDiagonal(totalCells),nnzOffDiag(totalCells))

     call initialize_stencils
     if (not(viscous)) then
        n_stencil = N_euler_drdw
        allocate(stencil(n_stencil,3))
        stencil = euler_drdw_stencil
     else
        n_stencil = N_visc_pc
        allocate(stencil(n_stencil,3))
        stencil = visc_pc_stencil
     end if

     call statePreAllocation(nnzDiagonal,nnzOffDiag,nDimW/nw,stencil,n_stencil)
  
     call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
          nDimW, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dRdWPre, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     deallocate(nnzDiagonal,nnzOffDiag,stencil)

     ! Setup Matrix-Free dRdw matrix and its function
     call MatCreateMFFD(sumb_comm_world,nDimW,nDimW,&
          PETSC_DETERMINE,PETSC_DETERMINE,dRdw,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call MatMFFDSetFunction(dRdw,FormFunction_mf,ctx,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set the mat_row_oriented option to false so that dense
     ! subblocks can be passed in in fortran column-oriented format
     call MatSetOption(dRdWPre, MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call MatSetOption(dRdW   , MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     !  Create the linear solver context
     call KSPCreate(SUMB_PETSC_COMM_WORLD,global_ksp,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set operators for the solver
     call KSPSetOperators(global_ksp,dRdw,dRdWPre, DIFFERENT_NONZERO_PATTERN,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     NKSolverSetup = .True.
     NKSolvedOnce = .False.
     if(equations == RANSEquations)  then
        turbCoupled = .True.
        turbSegregated = .False.
     end if
  end if
end subroutine setupNKsolver


subroutine MyMult(matrix,X,F,ierr)

  !   Input Parameters:
  !.  X - input vector
  !
  !   Output Parameter:
  !.  F - function vector
  !
  use precision
  use NKSolverVars, only: dRdw,diagV
  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
  
  ! Input/Output Vars
  Mat matrix
  Vec X,F
  integer(kind=intType) :: ierr

  ! Do a matmult followed by an addition

  call MatMult(dRdw,X,F,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! VecAXPY : Computes y = alpha x + y. 
  ! VecAXPY(Vec y,PetscScalar alpha,Vec x)

  !call VecAXPY(F,1.0,diagV,ierr)
  !call EChk(ierr,__FILE__,__LINE__)

end subroutine MyMult

