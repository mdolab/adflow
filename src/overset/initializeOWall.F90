subroutine initializeOWall(oWall)

  ! This routine allocates the data for the supplied oBlock using the
  !  data currently in blockPointers
  use constants
  use overset
  use blockPointers
  use adtAPI
  use BCTypes
  use cgnsGrid
  use communication
  implicit none 

  ! Input Params
  type(oversetBlock), intent(inout) :: oWall

  ! Working paramters
  integer(kind=intType) :: i, j, k, mm, nADT, nHexa, planeOffset
  integer(kind=intType) :: iStart, iEnd, jStart, jEnd, kStart, kEnd

  ! Set all the sizes for this block.
  oWall%il = il
  oWall%jl = jl
  oWall%kl = kl

  ! Call the custom build routine -- Serial only, only Hexa volumes,
  ! we supply our own ADT Type

  !call buildSerialHex(nHexa, nADT, oBlock%xADT, oBlock%hexaConn, oBlock%ADT)

  ! Flag this block as being allocated
  oWall%allocated = .True.

end subroutine initializeOWall
