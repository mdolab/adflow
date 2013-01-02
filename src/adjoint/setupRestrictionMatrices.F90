subroutine setupRestrictionMatrix(matrix, level)
#ifndef USE_NO_PETSC
  !     ******************************************************************
  !     *                                                                *
  !     * Compute the restriction matrix from level-1 to level           *
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
  real(kind=realType), dimension(:,:,:), pointer :: vvol
  real(kind=realType) :: vola

  blockLoop: do nn=1,nDom
     spectralLoop: do sps=1,nTimeIntervalsSpectral

        ! Set the pointer for the coarse grid
        call setPointers(nn, level, sps)
        
        ! Set pointer for fine grid volume
        vvol    => flowDoms(nn, level-1, sps)%vol

        ! Loop over rows of the coarse grid
        do k=2, kl

           kk  = mgKFine(k,1)
           kk1 = mgKFine(k,2)

           do j=2, jl

              jj  = mgJFine(j,1)
              jj1 = mgJFine(j,2)

              do i=2, il

                 ii  = mgIFine(i,1)
                 ii1 = mgIFine(i,2)

                 do l=1, nw
                    
                    ! Compute the sum of the fine grid volumes.

                    vola = vvol(ii,jj,kk)   + vvol(ii1,jj,kk)   &
                         + vvol(ii,jj1,kk)  + vvol(ii1,jj1,kk)  &
                         + vvol(ii,jj,kk1)  + vvol(ii1,jj,kk1)  &
                         + vvol(ii,jj1,kk1) + vvol(ii1,jj1,kk1)

                    vola = one/vola

                    ! -1 is for zero-based ordering
                    irow = globalCell(i,j,k)*nw + l - 1
                   
                    ! Set the 8 values 
                    icol =  flowDoms(nn, level-1, sps)%globalCell(ii, jj, kk)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 1.0/8.0, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    icol =  flowDoms(nn, level-1, sps)%globalCell(ii1, jj, kk)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 1.0/8.0, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    icol =  flowDoms(nn, level-1, sps)%globalCell(ii, jj1, kk)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 1.0/8.0, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    icol =  flowDoms(nn, level-1, sps)%globalCell(ii1, jj1, kk)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 1.0/8.0, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    !-------------------------------------------!
                    icol =  flowDoms(nn, level-1, sps)%globalCell(ii, jj, kk1)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 1.0/8.0, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    icol =  flowDoms(nn, level-1, sps)%globalCell(ii1, jj, kk1)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 1.0/8.0, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    icol =  flowDoms(nn, level-1, sps)%globalCell(ii, jj1, kk1)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 1.0/8.0, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    icol =  flowDoms(nn, level-1, sps)%globalCell(ii1, jj1, kk1)*nw + l -1
                    call MatSetValues(matrix, 1, irow, 1, icol, 1.0/8.0, &
                         INSERT_VALUES, ierr)
                    call EChk(ierr, __FILE__, __LINE__)

                    
!                     ! Set the 8 values 
!                     icol =  flowDoms(nn, level-1, sps)%globalCell(ii, jj, kk)*nw + l -1
!                     call MatSetValues(matrix, 1, irow, 1, icol, vvol(ii,jj,kk)*vola, &
!                          INSERT_VALUES, ierr)
!                     call EChk(ierr, __FILE__, __LINE__)

!                     icol =  flowDoms(nn, level-1, sps)%globalCell(ii1, jj, kk)*nw + l -1
!                     call MatSetValues(matrix, 1, irow, 1, icol, vvol(ii1,jj,kk)*vola, &
!                          INSERT_VALUES, ierr)
!                     call EChk(ierr, __FILE__, __LINE__)

!                     icol =  flowDoms(nn, level-1, sps)%globalCell(ii, jj1, kk)*nw + l -1
!                     call MatSetValues(matrix, 1, irow, 1, icol, vvol(ii,jj1,kk)*vola, &
!                          INSERT_VALUES, ierr)
!                     call EChk(ierr, __FILE__, __LINE__)

!                     icol =  flowDoms(nn, level-1, sps)%globalCell(ii1, jj1, kk)*nw + l -1
!                     call MatSetValues(matrix, 1, irow, 1, icol, vvol(ii,jj1,kk)*vola, &
!                          INSERT_VALUES, ierr)
!                     call EChk(ierr, __FILE__, __LINE__)

!                     !-------------------------------------------!
!                     icol =  flowDoms(nn, level-1, sps)%globalCell(ii, jj, kk1)*nw + l -1
!                     call MatSetValues(matrix, 1, irow, 1, icol, vvol(ii,jj,kk1)*vola, &
!                          INSERT_VALUES, ierr)
!                     call EChk(ierr, __FILE__, __LINE__)

!                     icol =  flowDoms(nn, level-1, sps)%globalCell(ii1, jj, kk1)*nw + l -1
!                     call MatSetValues(matrix, 1, irow, 1, icol, vvol(ii1,jj,kk1)*vola, &
!                          INSERT_VALUES, ierr)
!                     call EChk(ierr, __FILE__, __LINE__)

!                     icol =  flowDoms(nn, level-1, sps)%globalCell(ii, jj1, kk1)*nw + l -1
!                     call MatSetValues(matrix, 1, irow, 1, icol, vvol(ii,jj1,kk1)*vola, &
!                          INSERT_VALUES, ierr)
!                     call EChk(ierr, __FILE__, __LINE__)

!                     icol =  flowDoms(nn, level-1, sps)%globalCell(ii1, jj1, kk1)*nw + l -1
!                     call MatSetValues(matrix, 1, irow, 1, icol, vvol(ii1,jj1,kk1)*vola, &
!                          INSERT_VALUES, ierr)
!                     call EChk(ierr, __FILE__, __LINE__)
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
end subroutine setupRestrictionMatrix
