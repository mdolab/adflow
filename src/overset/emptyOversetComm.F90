!      ******************************************************************
!      *                                                                *
!      * Short cut function to make empty overset comm structure for    *
!      * problems that do not use overset meshes.
!      *                                                                *
!      ******************************************************************

subroutine emptyOversetComm(level, sps)

  use blockPointers
  use communication
  implicit none

  ! Function
  integer(kind=intType), intent(in) :: level, sps

  ! Working
  integer(Kind=intType) :: nn, mm, ierr

  commPatternOverset(level, sps)%nProcRecv = 0
  allocate(commPatternOverset(level, sps)%recvProc(0))
  allocate(commPatternOverset(level, sps)%nRecv(0))
  allocate(commPatternOverset(level, sps)%recvList(0))

  commPatternOverset(level, sps)%nProcSend = 0
  allocate(commPatternOverset(level, sps)%sendProc(0))
  allocate(commPatternOverset(level, sps)%nSend(0))
  allocate(commPatternOverset(level, sps)%sendList(0))

  internalOverset(level, sps)%nCopy = 0
  allocate(internalOverset(level, sps)%donorBlock(0))
  allocate(internalOverset(level, sps)%donorIndices(0, 3))
  allocate(internalOverset(level, sps)%donorInterp(0, 8))
  allocate(internalOverset(level, sps)%haloBlock(0))
  allocate(internalOverset(level, sps)%haloIndices(0, 3))

end subroutine emptyOversetComm
