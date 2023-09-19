
! @File    :   f2py_api.f90
! @Time    :   2023/07/19
! @Desc    :   Simple tests for GPU development

module f2py_api
    use matrixType, only: matType

    public :: my_matmult_api

    contains

    subroutine my_matmult_api(a, b, n, c)
        ! This is the public interface to the GPU matmult that
        ! we call from python
        use cudafor ! this has to be hear in order to make the module call
        use test_f2py
        implicit none
        real(8), intent(in) :: a(:, :)
        real(8), intent(in) :: b(:, :)
        integer, intent(in) :: n
        real(8), intent(out) :: c(n, n)
        integer, i,j,k
        integer, nThreads 
        ! --- GPU DEVICE PARAMETERS ---
        real(8), allocatable,device :: d_a(:, :), d_b(:, :), d_c(:, :)
        type(dim3) :: grid, tBlock
        type(matType) :: h_mat

        ! Variables for timing
        integer :: start, finish, count_rate
        real(8) :: elapsed_time


        allocate(mat(1))

        allocate(h_mat%a(n,n),h_mat%b(n,n),h_mat%c(n,n))
        h_mat%a = a
        h_mat%b = b
        h_mat%c = 0

        mat(1) = h_mat
        ! mat(1)%a = h_mat%a
        ! mat(1)%b = h_mat%b
        ! mat(1)%c = h_mat%c




        ! --- GPU device grid and block size ---
        nThreads = 32
        grid = dim3(ceiling(real(n) / nThreads), ceiling(real(n) / nThreads),1) ! we ask for more from the GPU using the 'ceiling()'
        tBlock = dim3(nThreads, nThreads,1) ! 32 by 32 is the maximum you can ask for
        print*,grid%x,grid%y
        print *,tBlock%x, tBlock%y
        c = 0.0
        ! do i =1,n
        !     do j =1,N
        !         do k = 1,n
        !             c(i,j) = c(i,j)+a(i,k)*b(k,j)
        !         end do 
        !     end do 
        ! end do
        ! --- Allocation ---
        ! allocate
        call system_clock(start, count_rate) ! get start time
        !call cuda
        call my_matmult_gpu<<<grid, tBlock>>>
        ! copy memory back
        h_mat = mat(1)
        c = h_mat%c
        ! c = mat(1)%c
        ! print *, "c(1,1) = ", c(1, 1)
        ! convert time to seconds and print
        call system_clock(finish) ! get finish time

        elapsed_time = real(finish - start, 8) / real(count_rate, 8)
        write (*, '(a,f9.4,a)') "gpu mat-mat product", elapsed_time, " seconds!"

    end subroutine my_matmult_api

end module f2py_api

