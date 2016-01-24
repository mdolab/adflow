subroutine registerOversetFringeTypes

  ! This routine setups the MPI data type for the overset fringe.

  use communication
  use overset
  use blockPointers
  implicit none

  ! Working Variables
  type(fringeType) :: fringe
  integer(kind=intType), parameter :: nfields=19
  integer(kind=intType) :: ierr
  integer :: iblock(nfields)
  integer(kind=MPI_ADDRESS_KIND) :: start, idisp(nfields)
  integer(kind=intType) :: itype(nfields), i

  ! ---------------------------------------------------
  !              Full data Type 
  ! ---------------------------------------------------
  call MPI_Get_address(fringe, start, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  i = 1

  iType(i) = sumb_real   ;  iBlock(i) = 3; call MPI_Get_address(fringe%x, idisp(i), ierr); i = i + 1
  iType(i) = sumb_real   ;  iBlock(i) = 1; call MPI_Get_address(fringe%quality, idisp(i), ierr); i = i + 1
  iType(i) = sumb_integer; iBlock(i) = 1;  call MPI_Get_address(fringe%donorProc, idisp(i), ierr); i = i + 1
  iType(i) = sumb_integer; iBlock(i) = 1;  call MPI_Get_address(fringe%donorBlock, idisp(i), ierr); i = i + 1
  iType(i) = sumb_integer; iBlock(i) = 1;  call MPI_Get_address(fringe%dI, idisp(i), ierr); i = i + 1
  iType(i) = sumb_integer; iBlock(i) = 1;  call MPI_Get_address(fringe%dJ, idisp(i), ierr); i = i + 1
  iType(i) = sumb_integer; iBlock(i) = 1;  call MPI_Get_address(fringe%dK, idisp(i), ierr); i = i + 1
  iType(i) = sumb_real   ; iBlock(i) = 3;  call MPI_Get_address(fringe%donorFrac, idisp(i), ierr); i = i + 1
  iType(i) = sumb_integer ; iBlock(i) = 8; call MPI_Get_address(fringe%gInd, idisp(i), ierr); i = i + 1
  iType(i) = MPI_LOGICAL ; iBlock(i) = 1;  call MPI_Get_address(fringe%isDonor, iDisp(i), ierr); i = i + 1
  iType(i) = MPI_LOGICAL ; iBlock(i) = 1;  call MPI_Get_address(fringe%isHole, iDisp(i), ierr); i = i + 1
  iType(i) = MPI_LOGICAL ; iBlock(i) = 1;  call MPI_Get_address(fringe%isCompute, iDisp(i), ierr); i = i + 1
  iType(i) = sumb_integer; iBlock(i) = 1;  call MPI_Get_address(fringe%myBlock, idisp(i), ierr); i = i + 1
  iType(i) = sumb_integer; iBlock(i) = 1;  call MPI_Get_address(fringe%myI, idisp(i), ierr); i = i + 1
  iType(i) = sumb_integer; iBlock(i) = 1;  call MPI_Get_address(fringe%myJ, idisp(i), ierr); i = i + 1
  iType(i) = sumb_integer; iBlock(i) = 1;  call MPI_Get_address(fringe%myK, idisp(i), ierr); i = i + 1
  iType(i) = MPI_LOGICAL ; iBlock(i) = 1;  call MPI_Get_address(fringe%isWall, iDisp(i), ierr); i = i + 1
  iType(i) = MPI_LOGICAL ; iBlock(i) = 1;  call MPI_Get_address(fringe%isWallDonor, iDisp(i), ierr); 
  
  ! Compute the dispalcements by subtracting off the starting point
  idisp = idisp - start

  call MPI_Type_create_struct(i, iblock, idisp, itype, oversetMPIFringe, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Finally register the new type
  call MPI_Type_commit(oversetMPIFringe, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! ---------------------------------------------------
  !             Search Coords data Type 
  ! ---------------------------------------------------

  i = 1

  iType(i) = sumb_real   ; iBlock(i) = 3; call MPI_Get_address(fringe%x, idisp(i), ierr); i = i + 1
  iType(i) = sumb_real   ; iBlock(i) = 1; call MPI_Get_address(fringe%origQuality, idisp(i), ierr); i = i + 1
  iType(i) = MPI_LOGICAL ; iBlock(i) = 1; call MPI_Get_address(fringe%isWall, iDisp(i), ierr); 

  ! Compute the dispalcements by subtracting off the starting point
  idisp = idisp - start

  call MPI_Type_create_struct(i, iblock, idisp, itype, oversetMPISearchCoord, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Finally register the new type
  call MPI_Type_commit(oversetMPISearchCoord, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

end subroutine registerOversetFringeTypes
