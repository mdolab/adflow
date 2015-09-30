subroutine initializeOversetComm

  ! This routine setups the data structures needed for overset
  ! communication. Specificially it creates the required petsc scatter
  ! context.
  use communication
  use overset
  use blockPointers
  use adjointVars
  use inputOverset
  implicit none
  ! Working Variables
  integer(kind=intType) :: ii, i, j, nn, ierr
  integer(kind=intType) :: nFringeProc, nDonorsPerCell
  integer(kind=intType), dimension(:), allocatable :: fringesProc

  ! First we determine the total number of donors and fringes on this block
  nFringeProc = 0
  if (oversetInterpolation == linear) then 
     nDonorsPerCell = 8
  else
     nDonorsPerCell = 27
  end if

  do nn=1, nDom
     ! Count up the total number of fringes (times 8) I need. This is
     ! quite inefficient but we'll fix later
     nFringeProc = nFringeProc + nDonorsPerCell*oBlocks(nn)%nFringe
  end do
  
  ! Allocate spaces for the fringes indices
  allocate(fringesProc(nFringeProc))
  ii = 0

  ! Copy in the indices each block needs. 
  do nn=1, nDom
     do i=1, oBlocks(nn)%nFringe
        do j=1, nDonorsPerCell
           ii = ii + 1
           ! Note that donorInidices are already zero-based and in
           ! SUmb's global petsc ordering.
           fringesProc(ii) = oBlocks(nn)%donorIndices(j, i)
          
        end do
     end do
  end do

  ! Create the two arrays needed for the overset comm. The first one
  ! is just a scalar version of w. That is just the number of local cells. 
  call VecCreateMPI(sumb_comm_world, nCellsLocal(1_intType), PETSC_DETERMINE, oversetDonors, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! The second is equal to the 8*the number of fringes we have. 
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
