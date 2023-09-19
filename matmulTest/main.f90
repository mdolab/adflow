program main
    ! This is the 'host' code

    use cudafor ! enable cuda fortran
    use f2py_api
    implicit none

    real(8), allocatable :: a(:, :), b(:, :), c(:, :)
    integer :: DIM = 32
    ! GPU device variables
    real, allocatable, device :: d_a(:, :), d_b(:, :), d_c(:, :)
    ! real, device :: x_d(NN), y_d(NN)
    type(dim3) :: grid, tBlock
    ! Variables for timing
    integer :: start, finish, count_rate
    real(8) :: elapsed_time

    print *, "Hello World!"

    grid = dim3(ceiling(real(DIM) / 32), ceiling(real(DIM) / 32), 1) ! we ask for more from the GPU using the 'ceiling()'
    tBlock = dim3(32, 32, 1) ! 32 by 32 is the maximum you can ask for

    allocate (a(DIM, DIM), b(DIM, DIM), c(DIM, DIM))
    allocate (d_a(DIM, DIM), d_b(DIM, DIM), d_c(DIM, DIM))

    a = 3.D0
    b = 2.D0
    c = 0.D0

    ! --- Copy input to the GPU ---
    d_a = a
    d_b = b
    d_c = c

    call system_clock(start, count_rate) ! get start time
    call my_matmult_api(a, b, DIM, c) ! so the stuff in the <<<>>> are the launch parameters for the GPU
    ! grid is how many blocks
    ! tBlock is how many threads per block
    call system_clock(finish) ! get finish time

    ! --- Copy results back to the host (CPU) ---
    ! c = d_c

    ! convert time to seconds and print
    elapsed_time = real(finish - start, 8) / real(count_rate, 8)
    write (*, '(a,f9.3,a)') "gpu elementwise mat product", elapsed_time, " seconds!"
    print *, "c(1,1) = ", c(1, 1)
    deallocate(d_a,d_b,d_c)
    deallocate(a,b,c)
    ! use cudafor
    ! implicit none
    ! integer, parameter :: N = 40000
    ! real :: x(N), y(N), a
    ! real, device :: x_d(N), y_d(N)
    ! type(dim3) :: grid, tBlock

    ! tBlock = dim3(256,1,1)
    ! grid = dim3(ceiling(real(N)/tBlock%x),1,1)

    ! x = 1.0; y = 2.0; a = 2.0
    ! x_d = x
    ! y_d = y
    ! call saxpy<<<grid, tBlock>>>(x_d, y_d, a)
    ! y = y_d
    ! write(*,*) 'Max error: ', maxval(abs(y-4.0))

end program main
