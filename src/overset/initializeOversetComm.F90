subroutine initializeOversetComm

  ! This routine setups the data structures needed for overset
  ! communication. Specificially it creates the required petsc scatter
  ! context.
  use communication
  use overset
  use blockPointers
  implicit none
  ! Working Variables
  integer(kind=intType) :: ii, i, nn, ierr
  integer(kind=intType) :: nDonorProc, nFringeProc
  integer(kind=intType), dimension(:), allocatable :: fringesProc

  ! First we determine the total number of donors and fringes on this block
  nDonorProc = 0
  nFringeProc = 0
  do nn=1,nDom
     call setPointers(nn, 1, 1)
     nDonorProc = nDonorProc + nDonor
     nFringeProc = nFringeProc + nFringe
  end do
  
  ! Allocate spaces for the fringes indices
  allocate(fringesProc(nFringeProc))
  ii = 0

  ! Copy in the indices each block needs. 
  do nn=1, nDom
     call setPointers(nn, 1, 1)
     do i = 1, nFringe
        ! We convert iBC which is in 1-based indexing to zero based
        ! indexing  or petscc
        ii = ii + 1
        fringesProc(ii) = iBC(i) - 1
     end do
  end do

  ! Create the two arrays needed for the overset comm
  call VecCreateMPI(sumb_comm_world, nDonorProc, PETSC_DETERMINE, oversetDonors, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  call VecCreateMPI(sumb_comm_world, nFringeProc, PETSC_DETERMINE, oversetFringes, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! The only other thing we need to do is create the actual scatter context. 
  call ISCreateGeneral(sumb_comm_world, nFringeProc, fringesProc, PETSC_COPY_VALUES, IS1, ierr)
  call ECHK(ierr, __FILE__, __LINE__)  
  
  ! We are putting all values into oversetFringes in order, hence the
  ! strided index set for IS2 
  call ISCreateStride(sumb_comm_world, nFringeProc, 0, 1, IS2, ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  
  call VecScatterCreate(oversetDonors, IS1, oversetFringes, IS2, oversetScatter, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Clean up the temporary IS1
  deallocate(fringesProc)

end subroutine initializeOversetComm
