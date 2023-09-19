module matrixType
contains
type matType
        real(8), allocatable,device,dimension(:,:) :: a(:,:), b(:,:), c(:,:)
end type matType
end module matrixType