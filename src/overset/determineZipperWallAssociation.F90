subroutine determineZipperWallAssociation(master, pocketMaster, fullWall)
           
  
  ! This routine determines global cell indices of the primal cells that
  ! contain the cell center of all zipper triangles. Saves them in
  ! master and pocketMaster strings.

  use adtAPI
  use blockPointers
  use communication
  use overset
  use constants
  implicit none

  ! Input Variables
  type(oversetWall), intent(inout) :: fullWall
  type(oversetString), intent(inout) :: master, pocketMaster

  ! Local Variables
  ! ------------------------------------------------------------
  integer(kind=intType) :: i, j, k, ii, kk, n(3), cellID

  ! Data for the ADT
  integer(kind=intType) :: intInfo(3), intInfo2(3)
  real(kind=realType) :: coor(4), uvw(5)
  real(kind=realType), dimension(3, 2) :: dummy
  real(kind=realType), parameter :: tol=1e-12
  integer(kind=intType), dimension(:), pointer :: frontLeaves, frontLeavesNew, BBint
  type(adtBBoxTargetType), dimension(:), pointer :: BB
  real(kind=realType), dimension(3) :: xp

  ! Debug related
  character(80) :: filename
  integer(kind=intType) :: nQuads, iQuad
  integer(kind=intType), dimension(:,:), allocatable :: quadConn
  ! ------------------------------------------------------------

  nQuads = master%nTris + pocketMaster%nTris
  allocate(quadConn(4, nQuads))


  ! Allocate the (pointer) memory that may be resized as necessary for
  ! the singlePoint search routine. 
  allocate(stack(100), BB(20), BBint(20), frontLeaves(25), frontLeavesNew(25))

  ! We need to store the global cell index defining the primal quad that
  ! each point has the closest point wrt.

  ! Allocate surfCellID of master and pocketMaster
  if (.not. associated(master%surfCellID)) then
     allocate(master%surfCellID(master%nTris))
  end if
  if (.not. associated(pocketMaster%surfCellID)) then
     allocate(pocketMaster%surfCellID(pocketMaster%nTris))
  end if

  ! Loop through all master triangles
  iQuad = 0
  do i=1, master%nTris
     n(1:3) = master%tris(:, i) ! global strings nodes
     coor(1) = third*(master%x(1, n(1)) + master%x(1, n(2)) + &
                      master%x(1, n(3)))
     coor(2) = third*(master%x(2, n(1)) + master%x(2, n(2)) + &
                      master%x(2, n(3)))
     coor(3) = third*(master%x(3, n(1)) + master%x(3, n(2)) + &
                      master%x(3, n(3)))

     coor(4) = large
     intInfo(3) = 0 
     call minDistancetreeSearchSinglePoint(fullWall%ADT, coor, &
          intInfo, uvw, dummy, 0, BB, frontLeaves, frontLeavesNew)
     cellID = intInfo(3)

     if (cellID > 0) then
        master%surfCellID(i) = fullWall%indCell(cellID)
     else
        ! Set dummy info. Should not be here mostly.
        master%surfCellID(i) = 0
     end if
 
     iQuad = iQuad + 1
     quadConn(:, iQuad) = fullWall%conn(:, cellID)
  end do

  ! Loop through all pocket triangles
  do i=1, pocketMaster%nTris
     n(1:3) = pocketMaster%tris(:, i) ! global strings nodes
     coor(1) = third*(pocketMaster%x(1, n(1)) + pocketMaster%x(1, n(2)) + &
                      pocketMaster%x(1, n(3)))
     coor(2) = third*(pocketMaster%x(2, n(1)) + pocketMaster%x(2, n(2)) + &
                      pocketMaster%x(2, n(3)))
     coor(3) = third*(pocketMaster%x(3, n(1)) + pocketMaster%x(3, n(2)) + &
                      pocketMaster%x(3, n(3)))

     coor(4) = large
     intInfo(3) = 0 
     call minDistancetreeSearchSinglePoint(fullWall%ADT, coor, &
          intInfo, uvw, dummy, 0, BB, frontLeaves, frontLeavesNew)
     cellID = intInfo(3)

     if (cellID > 0) then
        pocketMaster%surfCellID(i) = fullWall%indCell(cellID)
     else
        ! Set dummy info. Should not be here mostly.
        pocketMaster%surfCellID(i) = 0
     end if

     iQuad = iQuad + 1
     quadConn(:, iQuad) = fullWall%conn(:, cellID)
  end do

  !! Debug containment
  !! ----------------------
  !write (fileName,"(a,I2.2,a)") "fullwallSearch_", myid, ".dat"
  !open(unit=101,file=trim(fileName),form='formatted')
  !write(101,*) 'TITLE = "mywalls"'
  !write(101,*) 'Variables = "X", "Y", "Z"'
  !write(101,*) "Zone T=fullwall"
  !write (101,*) "Nodes = ", fullWall%nNodes, " Elements= ", master%nTris+pocketMaster%nTris, " ZONETYPE=FEQUADRILATERAL"
  !write(101, *) "DATAPACKING=POINT"
  !do i=1, fullWall%nNodes
  !   write(101, '(3(E20.12,x))')fullWall%x(1:3,i)
  !end do
  !do i=1, nQuads
  !   write(101, '(4(I5,x))') quadConn(1:4, i)
  !end do
  !close(101)
  !! ----------------------

  deallocate(quadConn)

end subroutine determineZipperWallAssociation
