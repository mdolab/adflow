
! @File    :   test_f2py.cuf
! @Time    :   2023/07/19
! @Desc    :   Simple tests for GPU development

module constants
    integer, parameter :: ConstN = 32
end module constants

module test_f2py

    ! private :: my_matmult_gpu
    use matrixType, only: matType
    implicit none
    type(matType), dimension(:),allocatable,device :: mat(:)
    contains

    attributes(global) subroutine my_matmult_gpu
        use constants
        implicit none
        real(8), device,pointer :: a(:,:), b(:,:), c(:,:)
        ! --- Working vars ---
        real(8) :: tmp
        integer :: n
        integer :: ii, jj ! thread indices
        integer :: kk !loop index
        integer :: blockSize
        integer :: tx, ty, bx, by
        ! real(8), shared :: aTile(32, 32), bTile(32, 32)
        integer :: ll, mm ! loop indices for tiling
        integer :: rowGlob, colGlob
        
        n = size(mat(1)%a, dim=1)
        

        ! At a bare minimum, any thread that can do work has at minimum two identifying numbers:
        ! block index
        ! thread index
        ! There is also blockDim which is how many threads are in a block
        ! Convenience
        tx = threadIdx%x
        ty = threadIdx%y
        bx = blockIdx%x
        by = blockIdx%y
        ! ************************************************
        ! !     Naive implementation of matrix multiplication (coalleced too because its fortran :)
        ! ! ************************************************
        ! ! Ok so now...I've made it so there are 2 id numbers within the thread index so there's %x and %y
        ii = tx + (bx - 1) * blockDim%x
        jj = ty + (by - 1) * blockDim%y
        ! ! The blockIdx%x - 1 is because FORTRAN is 1-based
        if (ii <= n .and. jj <= n) then
            tmp = 0.0
            do kk = 1, n
                ! tmp=0.0
                tmp = tmp + mat(1)%a(ii, kk) * mat(1)%b(kk, jj)
            end do
            mat(1)%c(ii, jj) = tmp
        end if

        ! ************************************************
        !     Tiling implementation of matrix multiplication
        ! ************************************************
        ! if (ii <= n .and. jj <= n) then
        !     tmp = 0.0
        !     rowGlob = (by - 1) * 32
        !     colGlob = (bx - 1) * 32

        !     ! ---------------------------
        !     !   Loop over cache blocks
        !     ! ---------------------------
        !     do kk = 0, n - 1, 32
        !         !copy data into shared memory
        !         aTile(ty, tx) = a(rowGlob + ty, kk + tx)
        !         bTile(ty, tx) = b(kk + ty, colGlob + tx)

        !         !print *, "aTile: ", aTile(ty, tx), "bTile: ", bTile(ty, tx), "tx: ", tx, "   ty: ", ty
        !         print *, aTile(ty,tx), bTile(ty,tx),tx,ty

        !         call syncthreads()
        !         ! Add contribution of cache block back to matrix
        !         ! --- Loop within the tile ---
        !         do mm = 1, 32
        !             ! row of aTile times column of bTile
        !             tmp = tmp + aTile(ty, mm) * bTile(mm, tx)
        !         end do

        !         call syncthreads()
        !     end do

        !     c(ii, jj) = tmp

        ! end if

    end subroutine my_matmult_gpu

end module test_f2py

