subroutine createStatePETScVars

#ifndef USE_NO_PETSC
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create the matrices/vectors that are required for the adjoint  *
  !     * solution: Matrices: dRdwT (dRdwPreT)                           *
  !     *           Vectors: dRdw, psi, adjointRes, adjointRHS,             *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc, only: dRdwT, drdwPreT, dJdw, psi, adjointRes, adjointRHS, &
       PETScIerr, PETScBlockMatrix, coarsedRdwPreT, restrictionOperator, &
       prolongationOperator
  use ADjointVars   
  use communication  
  use inputTimeSpectral 
  use flowVarRefState 
  use inputADjoint    
  use stencils
  use blockPointers
  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"


  !     Local variables.
  integer(kind=intType)  :: nDimW, nDimX, nDimw_fine, nDimw_coarse, l
  integer(kind=intType) :: i, n_stencil
  integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  integer(kind=intType), dimension(:), allocatable :: nnzDiagonal2, nnzOffDiag2
  integer(kind=intType), dimension(:, :), allocatable :: stencil
  integer(kind=intType) :: level, ierr, nlevels

  ! Define matrix dRdW local size, taking into account the total
  ! number of Cells owned by the processor and the number of 
  ! equations.

  nDimW = nw * nCellsLocal(1_intType)*nTimeIntervalsSpectral
  nDimX = 3 * nNodesLocal(1_intType)*nTimeIntervalsSpectral
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create matrix dRdW that define the adjoint linear system of    *
  !     * equations, dRdW^T psi = dJdW. Matrix dRdW has size [nDimW, nDimW]*
  !     * but is very sparse because of the computational stencil R=R(W).*
  !     *                                                                *
  !     ******************************************************************
  !

  ! ------------------- Determine Preallocation for dRdw --------------

  allocate(nnzDiagonal(nCellsLocal(1_intType)*nTimeIntervalsSpectral), &
            nnzOffDiag(nCellsLocal(1_intType)*nTimeIntervalsSpectral) )

  call initialize_stencils
  if (.not. viscous) then
     n_stencil = N_euler_drdw
     allocate(stencil(n_stencil, 3))
     stencil = euler_drdw_stencil
  else
     n_stencil = N_visc_pc
     allocate(stencil(n_stencil, 3))
     stencil = visc_pc_stencil
  end if

  level = 1
  call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nw, stencil, n_stencil, &
       level)

  if( nw <= 7 ) then

     PETScBlockMatrix = .true.

     ! Create a Block AIJ Matrix with block size nw, (number of states)
     if (PETSC_VERSION_MINOR <  3) then
        call MatCreateMPIBAIJ(SUMB_COMM_WORLD, nw,             &
             nDimW, nDimW,                     &
             PETSC_DETERMINE, PETSC_DETERMINE, &
             0, nnzDiagonal,         &
             0, nnzOffDiag,            &
             dRdWT, PETScIerr)
     else
        call MatCreateBAIJ(SUMB_COMM_WORLD, nw,             &
             nDimW, nDimW,                     &
             PETSC_DETERMINE, PETSC_DETERMINE, &
             0, nnzDiagonal,         &
             0, nnzOffDiag,            &
             dRdWT, PETScIerr)
     end if
     call EChk(PETScIerr, __FILE__, __LINE__)
  else

     PETScBlockMatrix = .false.

     allocate(nnzDiagonal2(nDimw), nnzOffDiag2(nDimw))
     ! The drdw prealloc function is done per block, which we can use
     ! to compute the correct preallocation for the non-block matrix:

     do i=1, nCellsLocal(1_intType)*nTimeIntervalsSpectral
        nnzDiagonal2((i-1)*nw+1:(i-1)*nw+nw) = nnzDiagonal(i)
        nnzOffDiag((i-1)*nw+1:(i-1)*nw+nw) = nnzOffDiag(i)
     end do
     if (PETSC_VERSION_MINOR <  3) then
        call MatCreateMPIAIJ(SUMB_COMM_WORLD,                 &
             nDimW, nDimW,                     &
             PETSC_DETERMINE, PETSC_DETERMINE, &
             8, nnzDiagonal2,         &
             8, nnzOffDiag2,            &
             dRdWT, PETScIerr)
     else
        call MatCreateAIJ(SUMB_COMM_WORLD,                 &
             nDimW, nDimW,                     &
             PETSC_DETERMINE, PETSC_DETERMINE, &
             8, nnzDiagonal2,         &
             8, nnzOffDiag2,            &
             dRdWT, PETScIerr)
     end if
     call EChk(PETScIerr, __FILE__, __LINE__)

     deallocate(nnzDiagonal2, nnzOffDiag2)
  endif

  deallocate(nnzDiagonal, nnzOffDiag, stencil)

  ! Set the matrix dRdW options.

  ! Warning: The array values is logically two-dimensional, 
  ! containing the values that are to be inserted. By default the
  ! values are given in row major order, which is the opposite of
  ! the Fortran convention, meaning that the value to be put in row
  ! idxm[i] and column idxn[j] is located in values[i*n+j]. To allow
  ! the insertion of values in column major order, one can call the
  ! command MatSetOption(Mat A, MAT COLUMN ORIENTED);

  call MatSetOption(dRdWt, MAT_ROW_ORIENTED, PETSC_FALSE, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)

  if (ApproxPC) then
     ! ------------------- Determine Preallocation for dRdwPre -------------
     allocate(nnzDiagonal(nCellsLocal(1_intType)*nTimeIntervalsSpectral), &
              nnzOffDiag(nCellsLocal(1_intType)*nTimeIntervalsSpectral) )

     n_stencil = N_euler_PC
     allocate(stencil(n_stencil, 3))
     stencil = euler_PC_stencil

     level = 1
     call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nw, stencil, &
          n_stencil, level)
     ! --------------------------------------------------------------------

     if( nw <= 7 ) then

        PETScBlockMatrix = .true.
        if (PETSC_VERSION_MINOR <  3) then
           call MatCreateMPIBAIJ(SUMB_COMM_WORLD, nw,             &
                nDimW, nDimW,                     &
                PETSC_DETERMINE, PETSC_DETERMINE, &
                0, nnzDiagonal,         &
                0, nnzOffDiag,            &
                dRdWPreT, PETScIerr)
        else
           call MatCreateBAIJ(SUMB_COMM_WORLD, nw,             &
                nDimW, nDimW,                     &
                PETSC_DETERMINE, PETSC_DETERMINE, &
                0, nnzDiagonal,         &
                0, nnzOffDiag,            &
                dRdWPreT, PETScIerr)
        end if
        call EChk(PETScIerr, __FILE__, __LINE__)
     else

        PETScBlockMatrix = .false.
        allocate(nnzDiagonal2(nDimw), nnzOffDiag2(nDimw))
        ! The drdwPC prealloc function is done per block, which we can use
        ! to compute the correct preallocation for the non-block matrix:

        do i=1, nCellsLocal(1_intType)*nTimeIntervalsSpectral
           nnzDiagonal2((i-1)*nw+1:(i-1)*nw+nw) = nnzDiagonal(i)
           nnzOffDiag((i-1)*nw+1:(i-1)*nw+nw) = nnzOffDiag(i)
        end do
        if (PETSC_VERSION_MINOR <  3) then
           call MatCreateMPIAIJ(SUMB_COMM_WORLD,                 &
                nDimW, nDimW,                     &
                PETSC_DETERMINE, PETSC_DETERMINE, &
                0, nnzDiagonal2,         &
                0, nnzOffDiag2,            &
                dRdWPret, PETScIerr)
        else
           call MatCreateAIJ(SUMB_COMM_WORLD,                 &
                nDimW, nDimW,                     &
                PETSC_DETERMINE, PETSC_DETERMINE, &
                0, nnzDiagonal2,         &
                0, nnzOffDiag2,            &
                dRdWPret, PETScIerr)
        end if
        call EChk(PETScIerr, __FILE__, __LINE__)

        deallocate(nnzDiagonal2, nnzOffDiag2)
     endif

     deallocate(nnzDiagonal, nnzOffDiag, stencil)

     ! Set the matrix dRdWPre options.
     call MatSetOption(dRdWPret, MAT_ROW_ORIENTED, PETSC_FALSE, PETScIerr)
     call EChk(PETScIerr, __FILE__, __LINE__)
  end if ! Approx PC

  ! Vectors:
  ! Get dJdw and psi from one MatGetVecs Call
  call MatGetVecs(dRdwT, dJdW, psi, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)

  ! adjointRes is the same size as dJdw, psi
  call VecDuplicate(dJdW, adjointRes, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)

  call VecDuplicate(dJdW, adjointRHS, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)

  ! If we are using the multigrind PC we alos have to assemble coarse
  ! grid approximations to the Jacobian. 
  nLevels = ubound(flowDoms, 2)

  if (preCondType == 'mg') then

     ! Allocate coarse grid 
     allocate(coarsedRdwPreT(nlevels), stat=ierr)
     call EChk(PETScIerr, __FILE__, __LINE__)

     allocate(restrictionOperator(nlevels), stat=ierr)
     call EChk(PETScIerr, __FILE__, __LINE__)

     allocate(prolongationOperator(nlevels), stat=ierr)
     call EChk(PETScIerr, __FILE__, __LINE__)

     ! Loop over coarse levels:
     do l=2, nLevels
        
        ! Compute the nDimw for this level
        nDimW_coarse = nw * nCellsLocal(l  )*nTimeIntervalsSpectral
        nDimW_fine   = nw * nCellsLocal(l-1)*nTimeIntervalsSpectral

         ! Allocate sizes
         allocate(nnzDiagonal(nCellsLocal(l)*nTimeIntervalsSpectral), &
                   nnzOffDiag(nCellsLocal(l)*nTimeIntervalsSpectral), &
                   nnzDiagonal2(nCellsLocal(l-1)*nTimeIntervalsSpectral), &
                   nnzOffDiag2(nCellsLocal(l-1)*nTimeIntervalsSpectral))
           nnzDiagonal2(:) = 8
           nnzOffDiag2(:) = 0

         call EChk(ierr, __FILE__, __LINE__)

        n_stencil = N_euler_PC
        allocate(stencil(n_stencil, 3))
        stencil = euler_PC_stencil

        call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW_coarse/nw, stencil, &
          n_stencil, l)

        PETScBlockMatrix = .true.
        if (PETSC_VERSION_MINOR <  3) then
           call MatCreateMPIBAIJ(SUMB_COMM_WORLD,  nw, &
                nDimW_coarse, nDimW_coarse, &
                PETSC_DETERMINE, PETSC_DETERMINE, &
                0, nnzDiagonal,         &
                0, nnzOffDiag,            &
                coarsedRdwPreT(l), PETScIerr)
           call EChk(PETScIerr, __FILE__, __LINE__)
        
           nnzDiagonal = 8
           nnzOffDiag = 0
           call MatCreateMPIBAIJ(SUMB_COMM_WORLD, nw, &
                nDimW_coarse, nDimW_fine, &
                PETSC_DETERMINE, PETSC_DETERMINE, &
                0, nnzDiagonal, &
                0, nnzOffDiag, &
                restrictionOperator(l), PETScIerr)
           call EChk(PETScIerr, __FILE__, __LINE__)

           nnzDiagonal2 = 8
           nnzOffDiag = 0
           call MatCreateMPIBAIJ(SUMB_COMM_WORLD, nw, &
                nDimW_fine, nDimW_coarse, &
                PETSC_DETERMINE, PETSC_DETERMINE, &
                0, nnzDiagonal2, &
                0, nnzOffDiag2, &
                prolongationOperator(l), PETScIerr)
           call EChk(PETScIerr, __FILE__, __LINE__)
        else
!            call MatCreateBAIJ(SUMB_COMM_WORLD, nw, &
!                 nDimW_coarse, nDimW_coarse,       &
!                 PETSC_DETERMINE, PETSC_DETERMINE, &
!                 0, nnzDiagonal,         &
!                 0, nnzOffDiag,            &
!                 coarsedRdwPreT(l), PETscIerr)
!            call EChk(PETScIerr, __FILE__, __LINE__)

           nnzDiagonal = 8
           nnzOffDiag = 0
           call MatCreateBAIJ(SUMB_COMM_WORLD, nw, &
                nDimW_coarse, nDimW_fine, &
                PETSC_DETERMINE, PETSC_DETERMINE, &
                0, nnzDiagonal, &
                0, nnzOffDiag, &
                restrictionOperator(l), PETScIerr)
           call EChk(PETScIerr, __FILE__, __LINE__)

           nnzDiagonal2 = 8
           nnzOffDiag2 = 0
           call MatCreateBAIJ(SUMB_COMM_WORLD, nw, &
                nDimW_fine, nDimW_coarse, &
                PETSC_DETERMINE, PETSC_DETERMINE, &
                0, nnzDiagonal2, &
                0, nnzOffDiag2, &
                prolongationOperator(l), PETScIerr)
           call EChk(PETScIerr, __FILE__, __LINE__)
         end if
         deallocate(nnzDiagonal, nnzOffDiag, nnzDiagonal2, nnzOffDiag2, stencil)
      end do
  end if
#endif
end subroutine createStatePETScVars

subroutine createSpatialPETScVars
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create the matrices/vectors that are required for the total    *
  !     * derivative: Matrices: dRdx                                     *
  !     *              Vectors: dJdx                                     *
  !     *                                                                *
  !     ******************************************************************

  use ADjointPETSc, only : dRdx, PETScIerr, dJdx, Xvec
  use ADjointVars     
  use communication   
  use inputTimeSpectral
  use flowVarRefState 
  use inputADjoint  

  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
#include "include/petscversion.h"

  ! Local variables.
  integer(kind=intType) :: nDimW, nDimX
  integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  integer(kind=intType) :: level
#ifndef USE_NO_PETSC

  nDimW = nw * nCellsLocal(1_intType)*nTimeIntervalsSpectral
  nDimX = 3 * nNodesLocal(1_intType)*nTimeIntervalsSpectral

  !     ******************************************************************
  !     *                                                                *
  !     * Create matrix dRdx that is used to compute the total cost /    *
  !     * constraint function sensitivity with respect to the spatial    *
  !     * design variables 'x' as dIdx = dJdx - psi^T dRdx.              *
  !     *                                                                *
  !     * Matrix dRdx has size [nDimW, nDimX] and is generally            *
  !     * sparse for the coordinate design variables.                    *
  !     *                                                                *
  !     * The local dimensions are specified so that the spatial         *
  !     * coordinates x (a) are placed in the local processor. This has  *
  !     * to be consistent with the vectors dIdx and dJdx.               *
  !     *                                                                *
  !     ******************************************************************

  allocate( nnzDiagonal(nDimX), nnzOffDiag(nDimX) )
  ! Create the matrix dRdx.
  level = 1_intType
  call drdxPreAllocation(nnzDiagonal, nnzOffDiag, nDimX, level)
  nnzDiagonal = int(nnzDiagonal * 1.2)
  nnzOffDiag = int(nnzOffDIag * 1.2)
  ! Note we are creating the TRANPOSE of dRdx. It is size dDimX by nDimW
  if (PETSC_VERSION_MINOR <  3) then
     call MatCreateMPIAIJ(SUMB_COMM_WORLD, &
          nDimX, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          8, nnzDiagonal,     &
          8, nnzOffDiag,            &
          dRdx, PETScIerr)
     call EChk(PETScIerr, __FILE__, __LINE__)
  else
     call MatCreateAIJ(SUMB_COMM_WORLD, &
          nDimX, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          8, nnzDiagonal,     &
          8, nnzOffDiag,            &
          dRdx, PETScIerr)
     call EChk(PETScIerr, __FILE__, __LINE__)
  end if

  deallocate( nnzDiagonal, nnzOffDiag )

  ! Set column major order for the matrix dRdx.
  call MatSetOption(dRdx, MAT_ROW_ORIENTED, PETSC_FALSE, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)

  call MatSetOption(dRdx, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)


  ! Vectors
  call VecCreate(SUMB_COMM_WORLD, dJdx, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)
  
  call VecSetSizes(dJdx, nDimX, PETSC_DECIDE, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)
  
  call VecSetBlockSize(dJdx, 3, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)
  
  call VecSetType(dJdx, VECMPI, PETScIerr) 
  call EChk(PETScIerr, __FILE__, __LINE__)

  ! xVec
  call VecDuplicate(dJdx, xVec, PETScIerr)
  call EChk(PETScierr, __FILE__, __LINE__)

#endif
end subroutine createSpatialPETScVars

subroutine createPETScKsp

  use ADjointPETSc
  use communication
  implicit none
  
#ifndef USE_NO_PETSC
  call KSPCreate(SUMB_COMM_WORLD, adjointKSP, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)
#endif

end subroutine createPETScKsp

subroutine createExtraPETScVars

  use ADjointPETSc, only : dRda, dRda_data, PETScIerr
  use ADjointVars     
  use communication   
  use inputTimeSpectral
  use flowvarrefstate
  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
#include "include/petscversion.h"

  !     Local variables.
  integer(kind=intType) :: nDimW
  !
#ifndef USE_NO_PETSC

  nDimW = nw * nCellsLocal(1_intType)*nTimeIntervalsSpectral
 
  ! dRda

  ! Once again, PETSC is royally screwed up. You CANNOT use PETSC_NULL
  ! arguments. They simply do NOT work in Fortran. The PETSc
  ! documentation lies to you. We have to allocate our own data. 
  if (allocated(dRda_data)) then
     deallocate(dRda_data)
  end if
  allocate(dRda_data(nDimw, nDesignExtra))

  if (PETSC_VERSION_MINOR < 3 ) then
     call MatCreateMPIDense(SUMB_COMM_WORLD, nDimW, PETSC_DECIDE, &
          PETSC_DETERMINE, nDesignExtra, dRda_data, dRda, PETScIerr)
  else
     call MatCreateDense(SUMB_COMM_WORLD, nDimW, PETSC_DECIDE, &
          PETSC_DETERMINE, nDesignExtra, dRda_data, dRda, PETScIerr)
  end if
  call EChk(PETScIerr, __FILE__, __LINE__)
#endif
end subroutine createExtraPETScVars
 
subroutine createCouplingPETScVars
#ifndef USE_NO_PETSC

  use ADjointPETSc, only: dFdx, dFdw, PETScIerr
  use ADjointVars   
  use communication
  use inputTimeSpectral
  use flowVarRefState
  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  !     Local variables.
  integer(kind=intType)  :: nDimW, nDimS, nTS
  integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag

  nDimW = nw * nCellsLocal(1_intType)*nTimeIntervalsSpectral
  call getForceSize(nDimS, nTS)
  nDimS = nDimS * 3 *nTimeIntervalsSpectral! Multiply by 3 for each
                                           ! dof on each point

  ! Create dFdx and dFdw

  ! Each nodal force is contribed by (nominally) 4 quadrilateral cells
  ! surrounding it. The pressure on each of these cells depend on the
  ! average of the pressure of the cell above and the halo below. The
  ! 1st halo is computed by linear pressure extrapolation from the
  ! first two cells ABOVE the surface. The results in each node being
  ! affected by 8 cells. All of these cells on on-processors so there
  ! should be zero offdiag entries. For dFdx, the coordinates of each
  ! of the 4 quadrilateral affects the force, so this results in 9
  ! points spatial points affecting the force. Each coordiante has 3
  ! dimension which results in 3x3x3=27 nonzeros per row
  ! don't know where the non-zeros will end up

  ! Create the matrix dFdw

  allocate( nnzDiagonal(nDimS), nnzOffDiag(nDimS) )
  nnzDiagonal = 8*nw
  nnzOffDiag  = 8*nw! Make the off diagonal the same, since we
  if (PETSC_VERSION_MINOR <  3) then
     call MatCreateMPIAIJ(SUMB_COMM_WORLD, &
          nDimS, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dFdw, PETScIerr)
  else
      call MatCreateAIJ(SUMB_COMM_WORLD, &
          nDimS, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dFdw, PETScIerr)
   end if
   call EChk(PETScIerr, __FILE__, __LINE__)

  ! Create the matrix dFdx
  nnzDiagonal = 27
  nnzOffDiag = 27
  if (PETSC_VERSION_MINOR <  3) then
     call MatCreateMPIAIJ(SUMB_COMM_WORLD, &
          nDimS, nDimS,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dFdx, PETScIerr)
  else
     call MatCreateAIJ(SUMB_COMM_WORLD, &
          nDimS, nDimS,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dFdx, PETScIerr)
  end if

  call EChk(PETScIerr, __FILE__, __LINE__)
  deallocate( nnzDiagonal, nnzOffDiag )

  ! Set column major order for the matrix dFdw.
  call MatSetOption(dFdw, MAT_ROW_ORIENTED, PETSC_TRUE, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)
  call MatSetOption(dFdx, MAT_ROW_ORIENTED, PETSC_TRUE, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)

#endif

end subroutine createCouplingPETScVars
