subroutine setupProlongationMatrix(matrix, level)
#ifndef USE_NO_PETSC
  !     ******************************************************************
  !     *                                                                *
  !     * Compute the prolongation matrix from level-1 to level           *
  !                                                                      *
  !     ******************************************************************
  !
  use blockPointers      
  use inputTimeSpectral
  use flowVarRefState

  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! PETSc Matrix Variable
  Mat matrix 

  ! Input Variable
  integer(kind=intType) :: level

  ! Local Variables
  integer(kind=intType) :: nn, sps, i, j, k, l, irow, icol, ierr
  integer(kind=intType) :: ii, ii1, jj, jj1, kk, kk1
  real(kind=realType) :: value

  blockLoop: do nn=1,nDom
     spectralLoop: do sps=1,nTimeIntervalsSpectral

        ! Set the pointer for the fine grid
        call setPointers(nn, level-1, sps)

        ! Loop over row of the fine grid
        do k=2, kl

           kk  = mgKCoarse(k,1)
           kk1 = mgKCoarse(k,2)

           do j=2, jl

              jj  = mgJCoarse(j,1)
              jj1 = mgJCoarse(j,2)

              do i=2, il

                 ii  = mgICoarse(i,1)
                 ii1 = mgICoarse(i,2)

                 do l=1, nw

                    ! -1 is for zero-based ordering
                    irow = globalCell(i,j,k)*nw + l - 1
                    
                    ! Set the 8 values 
                    icol = flowDoms(nn, level, sps)%globalCell(ii, jj, kk)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 0.421875_realType, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    icol = flowDoms(nn, level, sps)%globalCell(ii1, jj, kk)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 0.140625_realType, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    icol = flowDoms(nn, level, sps)%globalCell(ii, jj1, kk)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 0.140625_realType, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    icol = flowDoms(nn, level, sps)%globalCell(ii, jj, kk1)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 0.140625_realType, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    icol = flowDoms(nn, level, sps)%globalCell(ii1, jj1, kk)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 0.046875_realType, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    icol = flowDoms(nn, level, sps)%globalCell(ii1, jj, kk1)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 0.046875_realType, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    icol = flowDoms(nn, level, sps)%globalCell(ii, jj1, kk1)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 0.046875_realType, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    icol = flowDoms(nn, level, sps)%globalCell(ii1, jj1, kk1)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 0.015625_realType, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)                 
                 end do
              end do
           end do
        end do
     end do spectralLoop
  end do blockLoop

  ! Finally assemble the matrix
  call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  call MatAssemblyEnd  (matrix, MAT_FINAL_ASSEMBLY, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call MatSetOption(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif

end subroutine setupProlongationMatrix
