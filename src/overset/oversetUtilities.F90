
subroutine emptyFringe(fringe)

  use block
  implicit none
  ! Input/Output
  type(fringeType), intent(inout) :: fringe

  ! Initialize data in empty fringe
  fringe%x = (/large, large, large/)
  fringe%quality = large
  fringe%forceRecv = .False.
  fringe%donorProc = -1
  fringe%donorBlock = -1
  fringe%dI = -1
  fringe%dJ = -1
  fringe%dK = -1
  fringe%myBlock = -1
  fringe%myI = -1
  fringe%myJ = -1
  fringe%myK = -1
  fringe%donorFrac = -one
  fringe%gInd = -1
  fringe%isDonor = .False.
  fringe%isHole = .False.
  fringe%isCompute = .True.
  fringe%isWall = .False.
  fringe%isWallDonor = .False.
end subroutine emptyFringe

subroutine printOverlapMatrix(overlap)

  ! This is a debugging routine to print out the overlap matrix.

  use communication
  use overset
  implicit none

  ! Input/output
  type(CSRMatrix), intent(in) :: overlap

  ! Working
  integer(kind=intType) :: i, jj

  if (myid == 0) then 
     ! Now dump out who owns what:
     do i=1, overlap%nrow
        write(*, "(a,I4, a)", advance='no'), 'Row:', i, "   "
        do jj=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
           write(*, "(a,I2, a, I6)", advance='no'), "(", overlap%colInd(jj), ")", int(overlap%data(jj))
        end do
        write(*, *) " "
     end do

     print *, '--------------------------------------'
     ! Now dump out who owns what:
     do i=1, overlap%nRow
        write(*, "(a,I4, a)", advance='no'), 'Row:', i, "   "
        do jj=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
           write(*, "(a,I2, a, I6)", advance='no'), "(", overlap%colInd(jj), ")", int(overlap%assignedProc(jj))
        end do
        write(*, *) " "
     end do
  end if
end subroutine printOverlapMatrix

subroutine getCumulativeForm(sizeArray, n, cumArray)
  use precision
  implicit none


  ! Input/Output
  integer(kind=intType), dimension(n), intent(in) :: sizeArray
  integer(kind=intType), intent(in) :: n
  integer(kind=intType), dimension(0:n), intent(out) :: cumArray

  ! Working
  integer(kind=intType) :: i

  cumArray(0) = 0
  do i=1, n
     cumArray(i) = cumArray(i-1) + sizeArray(i)
  end do

end subroutine getCumulativeForm


subroutine transposeOverlap(A, B)

  ! Create the matrix Create the matrix transpose. 
  ! Inspired by: https://people.sc.fsu.edu/~jburkardt/f_src/sparsekit/sparsekit.f90
  use overset
  implicit none

  ! Input/Output
  type(CSRMatrix), intent(in) :: A
  type(CSRMatrix), intent(inout) :: B

  ! Working
  integer(kind=intType) :: col, colp1, i, k, next

  ! A CSR matrix is the same as a CSC matrix of the transpose. So
  ! essentially the algorithm is convert A as a CSC matrix to B
  ! (a CSR matrix)

  ! Allocate space for everything in B
  B%nnz = A%nnz
  B%nRow = A%nCol
  B%nCol = A%nRow
  allocate(B%data(B%nnz), B%colInd(B%nnz), &
       B%assignedProc(B%nnz), B%rowPtr(B%nRow + 1))
  B%allocated = .True.
  !  Compute lengths of rows of B (ie the columns of A)

  B%rowPtr = 0

  do i = 1, A%nRow
     do k = A%rowPTr(i), A%rowPtr(i+1)-1
        colp1 = A%colInd(k) +1
        B%rowPtr(colp1) = B%rowPtr(colp1) + 1

     end do
  end do
  !
  !  Compute pointers from lengths.
  !
  B%rowPtr(1) = 1
  do i = 1, A%nRow
     B%rowPtr(i+1) = B%rowPtr(i) + B%rowPtr(i+1)
  end do
  !
  !  Do the actual copying.
  !
  do i = 1, A%nRow
     do k = A%rowPtr(i), A%rowPtr(i+1)-1
        col = A%colInd(k)
        next = B%rowPtr(col)

        B%data(next) = A%data(k)
        B%assignedProc(next) = A%assignedProc(k) 

        B%colInd(next) = i 
        B%rowPtr(col) = next + 1
     end do
  end do

  !  Reshift B%rowPtr
  do i = A%nRow, 1, -1
     B%rowPtr(i+1) = B%rowPtr(i)
  end do
  B%rowPtr(1) = 1 

end subroutine transposeOverlap

subroutine deallocateCSRMatrix(mat1)

  use overset
  implicit none

  type(CSRMatrix), intent(inout) :: mat1

  if (mat1%allocated) then 
     deallocate(&
          mat1%data, &
          mat1%colInd, &
          mat1%rowPtr, &
          mat1%assignedProc)
  end if

end subroutine deallocateCSRMatrix

subroutine computeFringeProcArray(fringes, n, fringeProc, cumFringeProc, nFringeProc)
  ! Compute the breakpoints "cumFringeProc" for a list of sorted n
  ! fringes "fringes". nFringeProc is the total number of unique
  ! processors. fringeProc is the processor number for each section.
  use block
  use communication
  implicit none

  ! Input/Output
  type(fringeType), intent(in), dimension(n) :: fringes
  integer(kind=intType), intent(in) :: n
  integer(kind=intType), intent(out) :: nFringeProc
  integer(kind=intType), intent(out) :: fringeProc(nProc), cumFringeProc(1:nProc+1)

  ! Working
  integer(kind=intType) :: i, currentProc

  fringeProc = -1
  nFringeProc = 0
  cumFringeProc(1) = 1
  currentProc = -1

  do i=1, n
     if (currentProc /= fringes(i)%donorProc) then 
        nFringeProc = nFringeProc + 1
        cumFringeProc(nFringeProc) = i
        fringeProc(nFringeProc) = fringes(i)%donorProc
        currentProc = fringes(i)%donorProc
     end if
  end do

  ! Finally, the nFringeProc+1 entry is always n+1
  cumFringeProc(nFringeProc+1) = n + 1

end subroutine computeFringeProcArray


subroutine receiveFringes(status, n)
  ! This is a little helper routnine that receives the
  ! oversetFringeType and puts them in tmpFringes. This exact code
  ! shows up in a variety of places. 
  use communication
  use overset
  implicit none

  integer, intent(in) ::  status(MPI_STATUS_SIZE) 
  integer(kind=inttype), intent(out) ::  n
  integer(kind=intType) :: ierr

  ! Get size
  call MPI_Get_count(status, oversetMPIFringe, n, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Allocate space for temporary fringes if not big enough
  if (n > size(tmpFringes)) then 
     deallocate(tmpFringes)
     allocate(tmpFringes(n))
  end if
           
  ! Now actually receive the fringes
  call mpi_recv(tmpFringes, n, oversetMPIFringe, status(MPI_SOURCE), &
       status(MPI_TAG), SUmb_comm_world, status, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

end subroutine receiveFringes


  subroutine fracToWeights(frac, weights)
    use constants
    implicit none
    real(kind=realType), intent(in), dimension(3) :: frac
    real(kind=realType), intent(out), dimension(8) :: weights

    weights(1) = (one-frac(1))*(one-frac(2))*(one-frac(3))
    weights(2) = (    frac(1))*(one-frac(2))*(one-frac(3))
    weights(3) = (one-frac(1))*(    frac(2))*(one-frac(3))
    weights(4) = (    frac(1))*(    frac(2))*(one-frac(3))
    weights(5) = (one-frac(1))*(one-frac(2))*(    frac(3))
    weights(6) = (    frac(1))*(one-frac(2))*(    frac(3))
    weights(7) = (one-frac(1))*(    frac(2))*(    frac(3))
    weights(8) = (    frac(1))*(    frac(2))*(    frac(3))
  end subroutine fracToWeights
