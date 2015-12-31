subroutine registerOversetFringeType

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
  integer(kind=intType) :: itype(nfields)

  call MPI_Get_address(fringe, start, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  iType(1) = sumb_real   ;  iBlock(1) = 3;  call MPI_Get_address(fringe%x, idisp(1), ierr)

  iType(2) = sumb_real   ;  iBlock(2) = 1;  call MPI_Get_address(fringe%quality, idisp(2), ierr)

  iType(3) = MPI_LOGICAL ;  iBlock(3) = 1;  call MPI_Get_address(fringe%forceRecv, idisp(3), ierr)

  iType(4) = sumb_integer; iBlock(4) = 1;  call MPI_Get_address(fringe%donorProc, idisp(4), ierr)

  iType(5) = sumb_integer; iBlock(5) = 1;  call MPI_Get_address(fringe%donorBlock, idisp(5), ierr)

  iType(6) = sumb_integer; iBlock(6) = 1;  call MPI_Get_address(fringe%dI, idisp(6), ierr)

  iType(7) = sumb_integer; iBlock(7) = 1;  call MPI_Get_address(fringe%dJ, idisp(7), ierr)

  iType(8) = sumb_integer; iBlock(8) = 1;  call MPI_Get_address(fringe%dK, idisp(8), ierr)

  iType(9) = sumb_real   ; iBlock(9) = 3;  call MPI_Get_address(fringe%donorFrac, idisp(9), ierr)

  iType(10) = sumb_integer ; iBlock(10) = 8;  call MPI_Get_address(fringe%gInd, idisp(10), ierr)

  iType(11) = MPI_LOGICAL ; iBlock(11) = 1;  call MPI_Get_address(fringe%isDonor, iDisp(11), ierr)

  iType(12) = MPI_LOGICAL ; iBlock(12) = 1;  call MPI_Get_address(fringe%isHole, iDisp(12), ierr)

  iType(13) = MPI_LOGICAL ; iBlock(13) = 1;  call MPI_Get_address(fringe%isCompute, iDisp(13), ierr)

  iType(14) = sumb_integer; iBlock(14) = 1;  call MPI_Get_address(fringe%myBlock, idisp(14), ierr)

  iType(15) = sumb_integer; iBlock(15) = 1;  call MPI_Get_address(fringe%myI, idisp(15), ierr)

  iType(16) = sumb_integer; iBlock(16) = 1;  call MPI_Get_address(fringe%myJ, idisp(16), ierr)

  iType(17) = sumb_integer; iBlock(17) = 1;  call MPI_Get_address(fringe%myK, idisp(17), ierr)

  iType(18) = MPI_LOGICAL ; iBlock(18) = 1;  call MPI_Get_address(fringe%isWall, iDisp(18), ierr)

  iType(19) = MPI_LOGICAL ; iBlock(19) = 1;  call MPI_Get_address(fringe%isWallDonor, iDisp(19), ierr)


  ! Compute the dispalcements by subtracting off the starting point
  idisp = idisp - start

  call MPI_Type_create_struct(nfields, iblock, idisp, itype, oversetMPIFringe, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Finally register the new type
  call MPI_Type_commit(oversetMPIFringe, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

end subroutine registerOversetFringeType
