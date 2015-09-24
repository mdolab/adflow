!
!      ******************************************************************
!      *                                                                *
!      * pairSearch searches block 'A' for potential donors in block 'B'*
!      *                                                                *
!      ******************************************************************

subroutine pairSearch(A, B)

  use constants
  use overset
  use adtAPI

  implicit none

  ! Input/Ouput
  type(oversetBlock), intent(inout) :: A, B

  ! Working Varaibles
  integer(kind=intType) :: i, j, k, ii, jj, kk, mm, tmp
  integer(kind=intType) :: nCoor, nHexa, nInterpol, elemID

  integer(kind=intType), dimension(:), allocatable :: procID, elementID
  real(kind=realType), dimension(:, :), allocatable ::  uvw, arrInterpol
  integer(kind=adtElementType), dimension(:), allocatable :: elementType
  real(kind=realType) :: donorQual, xINterp(3)
  real(kind=realType) :: timeA, timeB
  real(kind=realType) :: f1, f2, f3, f4, f5, f6, f7, f8
  real(kind=realType) :: di0, di1, dj0, dj1, dk0, dk1

  ! Now search the coordinates from *A* using the
  ! containment search
  nCoor = A%nx * A%ny * A%nz

  ! Allocate a bunch of stuff we need:
  allocate(procID(nCoor), elementType(nCoor), elementID(nCoor), uvw(3, nCoor), &
       arrInterpol(1, nCoor))
  nInterpol = 1 ! we get the ADT to compute the interpolated volume for us. 

  call adtContainmentSearch(nCoor, A%xSearch, B%adtName, procID, elementType, &
       elementID, uvw, nInterpol, B%qualDonor, arrInterpol)

  ! Check all the cells to see if we found anything:
  mm = 0
  do k=2, A%kl
     do j=2, A%jl
        do i=2, A%il
           mm = mm + 1

           elemFound: if (.not. procID(mm) < 0) then 

              donorQual = arrInterpol(1, mm)

              if (donorQual < .8*A%qualRecv(1, mm)) then 
                 elemID = elementID(mm)
                 
                 ! Increment the number of fringes on this block
                 A%nFringe = A%nFringe + 1

                 ! Save the indices of this fringe:
                 A%fringeIndices(1, A%nFringe) = i
                 A%fringeIndices(2, A%nFringe) = j
                 A%fringeIndices(3, A%nFringe) = k

                 ! Donor is better than me: Flag myself as a fringe.
                 A%iBlank(i, j, k) = -1

                 ! Set my quality as a receiver as my donor
                 A%qualRecv(1, mm) = donorQual
                 B%qualDonor(1, elemID) = 1e-30

                 ! Save the fractions. 
                 A%donorFrac(:, A%nFringe) = uvw(:, mm)

                 ! Getting the global set of indices is slightly
                 ! tricker than it sounds. We frist have to
                 ! reconstruct the i,j,k indices of mm (which is in
                 ! block B ordering). From this "cell" we can get the
                 ! i,j,k coordinates of the nodes, which are actually
                 ! the cells of the primal mesh. We can then index
                 ! into the globalCell array

                 ! Remember we have (il, jl, kl) elements in the dual
                 ! mesh. 
                 tmp = elemID - 1
                 ii = mod(tmp, B%il) + 1
                 jj = mod(tmp/B%il, B%jl) + 1
                 kk = tmp/(B%il*B%jl) + 1

                 ! that unwinds our index. Now do the 2x2x2 loop and reassemble the index
                 A%donorIndices(1, A%nFringe) = B%globalCell(ii  , jj  , kk  )
                 A%donorIndices(2, A%nFringe) = B%globalCell(ii+1, jj  , kk  )
                 A%donorIndices(3, A%nFringe) = B%globalCell(ii  , jj+1, kk  )
                 A%donorIndices(4, A%nFringe) = B%globalCell(ii+1, jj+1, kk  )
                 A%donorIndices(5, A%nFringe) = B%globalCell(ii  , jj  , kk+1)
                 A%donorIndices(6, A%nFringe) = B%globalCell(ii+1, jj  , kk+1)
                 A%donorIndices(7, A%nFringe) = B%globalCell(ii  , jj+1, kk+1)
                 A%donorIndices(8, A%nFringe) = B%globalCell(ii+1, jj+1, kk+1)
           
              end if
           end if elemFound
        end do
     end do
  end do

  ! Ditch temporary memory
  deallocate(procID, elementType, elementID, uvw, arrInterpol)
end subroutine pairSearch
