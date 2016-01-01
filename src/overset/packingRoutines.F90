subroutine packOBlock(oBlock)

  use overset
  use constants
  implicit none

  ! Pack up everything we need for this block into its own buffer
  ! inlucding the data required for the ADTree

  ! Input/Output Parameters
  type(oversetBlock), intent(inout) :: oBlock

  ! Working paramters
  integer(kind=intType) :: rSize, iSize, i, j, k, nHexa, nADT

  iSize = 0
  rSize = 0

  ! Count up the integers we want to send:

  iSize = iSize + 14 ! All block indices

  iSize = iSize + size(oBlock%hexaConn)

  iSize = iSize + size(oBlock%globalCell)
  
  iSize = iSize + size(oBlock%nearWall)

  iSize = iSize + size(oBlock%invalidDonor)

  iSize = iSize + oBlock%ADT%nLeaves*2 ! The two itegers for the
  ! children in each leaf

  ! Count up the reals we ned to send:
  rSize = rSize + size(oBlock%qualDonor)

  rSize = rSize + size(oBlock%xADT)

  rSize = rSize + oBlock%ADT%nBBoxes*6 ! Cell bounding boxes

  ! Bounding boxes for leaves
  rSize = rSize + oBlock%ADT%nLeaves*12

  rSize = rSize + 1 ! Min block volume

  ! Allocate the buffers
  allocate(oBlock%rBuffer(rSize), oBlock%iBuffer(iSize))

  ! Reset the integer counter and add all the integers on this pass
  iSize = 0

  oBlock%iBuffer(1) = oBlock%ib
  oBlock%iBuffer(2) = oBlock%jb
  oBlock%iBuffer(3) = oBlock%kb
  oBlock%iBuffer(4) = oBlock%ie
  oBlock%iBuffer(5) = oBlock%je
  oBlock%iBuffer(6) = oBlock%ke
  oBlock%iBuffer(7) = oBlock%il
  oBlock%iBuffer(8) = oBlock%jl
  oBlock%iBuffer(9) = oBlock%kl
  oBlock%iBuffer(10)= oBlock%nx
  oBlock%iBuffer(11) = oBlock%ny
  oBlock%iBuffer(12) = oBlock%nz

  oBlock%iBuffer(13) = oBlock%proc
  oBlock%iBuffer(14) = oBlock%block

  iSize = iSize + 14

  nHexa = oBlock%il * oBlock%jl * oBlock%kl
  nADT = oBlock%ie * oBlock%je * oBlock%ke

  do j=1, nHexa
     do i=1, 8
        iSize = iSize + 1
        oBlock%iBuffer(iSize) = oBlock%hexaConn(i, j)
     end do
  end do

  do k=0, oBlock%kb
     do j=0, oBlock%jb
        do i=0, oBlock%ib
           iSize = iSize + 1
           oBlock%iBuffer(iSize) = oBlock%globalCell(i, j, k)
        end do
     end do
  end do

  do k=1, oBlock%ke
     do j=1, oBlock%je
        do i=1, oBlock%ie
           iSize = iSize + 1
           oBlock%iBuffer(iSize) = oBlock%nearWall(i, j, k)
        end do
     end do
  end do

  do k=1, oBlock%ke
     do j=1, oBlock%je
        do i=1, oBlock%ie
           iSize = iSize + 1
           oBlock%iBuffer(iSize) = oBlock%invalidDonor(i, j, k)
        end do
     end do
  end do

  do i=1, oBlock%ADT%nLeaves
     iSize = iSize + 1
     oBlock%iBuffer(iSize) = oBlock%ADT%ADTree(i)%children(1)
     iSize = iSize + 1
     oBlock%iBuffer(iSize) = oBlock%ADT%ADTree(i)%children(2)  
  end do

  ! Reset the real counter and add all the real values on this pass.
  rSize = 0

  do i=1, oBlock%ie * oBlock%je * oBlock%ke
     rSize = rSize + 1
     oBlock%rBuffer(rSize) = oBlock%qualDonor(1, i)
  end do

  do j=1, oBlock%ie * oBlock%je * oBlock%ke
     do i=1, 3
        rSize = rSize + 1
        oBlock%rBuffer(rSize) = oBlock%xADT(i, j)
     end do
  end do

  do i=1, oBlock%ADT%nBboxes
     oBlock%rBuffer(rSize+1:rSize+6) = oBlock%ADT%xBBox(:, i)
     rSize = rSize + 6
  end do

  do i=1, oBlock%ADT%nLeaves
     oBlock%rBuffer(rSize+1:rSize+6) = oBlock%ADT%ADTree(i)%xMin(:)
     rSize = rSize + 6

     oBlock%rBuffer(rSize+1:rSize+6) = oBlock%ADT%ADTree(i)%xMax(:)
     rSize = rSize + 6
  end do

  rSize = rSize + 1
  oBlock%rBuffer(rSize) = oBlock%minVol

end subroutine packOBlock

subroutine unpackOBlock(oBlock)

  use adtData
  use overset
  implicit none

  ! unPack everything we need for this block from its own buffer
  ! and reconstitute the data required for the ADTree. It is assumed
  ! the buffers are already allocated and the data is available. This
  ! does the exact *OPPOSITE* operation as the packBlock() routine

  ! Input/Output Parameters
  type(oversetBlock), intent(inout) :: oBlock

  ! Working paramters
  integer(kind=intType) :: rSize, iSize, i, j, k, nHexa, nADT

  ! Reset the integer counter and add all the integers on this pass
  iSize = 0

  oBlock%ib = oBlock%iBuffer(1) 
  oBlock%jb = oBlock%iBuffer(2)
  oBlock%kb = oBlock%iBuffer(3)
  oBlock%ie = oBlock%iBuffer(4)
  oBlock%je = oBlock%iBuffer(5)
  oBlock%ke = oBlock%iBuffer(6)
  oBlock%il = oBlock%iBuffer(7)
  oBlock%jl = oBlock%iBuffer(8)
  oBlock%kl = oBlock%iBuffer(9)
  oBlock%nx = oBlock%iBuffer(10)
  oBlock%ny = oBlock%iBuffer(11)
  oBlock%nz = oBlock%iBuffer(12)
  oBlock%proc = oBlock%iBuffer(13)
  oBlock%block = oBlock%iBuffer(14)
  iSize = iSize + 14

  nHexa = oBlock%il * oBlock%jl * oBlock%kl
  nADT = oBlock%ie * oBlock%je * oBlock%ke

  ! Allocate the remainder of the arrays in oBlock.
  allocate(oBlock%hexaConn(8, nHexa))
  allocate(oBlock%globalCell(0:oBlock%ib, 0:oBlock%jb, 0:oBlock%kb))
  allocate(oBlock%nearWall(1:oBlock%ie, 1:oBlock%je, 1:oBlock%ke))
  allocate(oBlock%invalidDonor(1:oBlock%ie, 1:oBlock%je, 1:oBlock%ke))
  allocate(oBlock%qualDonor(1, oBlock%ie * oBlock%je * oBlock%ke))
  allocate(oBlock%xADT(3, nADT))

  ! -------------------------------------------------------------------
  ! Once we know the sizes, allocate all the arrays in the
  ! ADTree. Since we are not going to call the *actual* build routine
  ! for the ADT, we need to set all the information ourselves. This
  ! essentially does the same thing as buildSerialHex.
  oBlock%ADT%adtType = adtVolumeADT
  oBlock%ADT%nNodes = nADT
  oBlock%ADT%nTetra = 0
  oBlock%ADT%nPyra = 0
  oBlock%ADT%nPrisms = 0
  oBlock%ADT%nTria = 0
  oBlock%ADT%nQuads = 0
  oBlock%ADT%coor => oBlock%xADT
  oBlock%ADT%hexaConn => oBlock%hexaConn
  nullify(oBlock%ADT%tetraConn, oBlock%ADT%pyraConn, oBlock%ADT%prismsConn)
  oBlock%ADT%nBBoxes = nHexa
  allocate(oBlock%ADT%xBBOX(6, nHexa))
  allocate(oBlock%ADT%elementType(nHexa))
  allocate(oBlock%ADT%elementID(nHexa))
  oBlock%ADT%comm = MPI_COMM_SELF
  oBlock%ADT%nProcs = 1
  oBlock%ADT%myID = 0

  ! All hexas
  oBlock%ADT%elementType = adtHexahedron

  do i=1,nHexa
     oBlock%ADT%elementID(i) = i
  end do

  oBlock%ADT%nLeaves = oBlock%ADT%nBBoxes - 1
  if(oBlock%ADT%nBBoxes <= 1) oBlock%ADT%nLeaves = oBlock%ADT%nLeaves + 1
  allocate(oBlock%ADT%ADTree(oBlock%ADT%nLeaves))

  ! -------------------------------------------------------------------

  ! Now continue copying out the integer values
  do i=1, nHexa
     do j=1, 8
        iSize = iSize + 1
        oBlock%hexaConn(j, i) = oBlock%iBuffer(iSize)
     end do
  end do

  do k=0, oBlock%kb
     do j=0, oBlock%jb
        do i=0, oBlock%ib
           iSize = iSize + 1
           oBlock%globalCell(i, j, k) = oBlock%iBuffer(iSize)
        end do
     end do
  end do

  do k=1, oBlock%ke
     do j=1, oBlock%je
        do i=1, oBlock%ie
           iSize = iSize + 1
           oBlock%nearWall(i, j, k) = oBlock%iBuffer(iSize)
        end do
     end do
  end do


  do k=1, oBlock%ke
     do j=1, oBlock%je
        do i=1, oBlock%ie
           iSize = iSize + 1
           oBlock%invalidDonor(i, j, k) = oBlock%iBuffer(iSize)
        end do
     end do
  end do

  do i=1, oBlock%ADT%nLeaves
     iSize = iSize + 1
     oBlock%ADT%ADTree(i)%children(1) = oBlock%iBuffer(iSize)
     iSize = iSize + 1
     oBlock%ADT%ADTree(i)%children(2) = oBlock%iBuffer(iSize)
  end do

  ! Now copy out the real values
  rSize = 0

  do i=1, oBlock%ie * oBlock%je * oBlock%ke
     rSize = rSize + 1
     oBlock%qualDonor(1, i) =  oBlock%rBuffer(rSize)
  end do

  do j=1, oBlock%ie * oBlock%je * oBlock%ke
     do i=1, 3
        rSize = rSize + 1
        oBlock%xADT(i, j) = oBlock%rBuffer(rSize)
     end do
  end do

  do i=1, oBlock%ADT%nBboxes
     oBlock%ADT%xBBox(:, i) = oBlock%rBuffer(rSize+1:rSize+6)
     rSize = rSize + 6
  end do

  do i=1, oBlock%ADT%nLeaves
     oBlock%ADT%ADTree(i)%xMin(:) = oBlock%rBuffer(rSize+1:rSize+6)
     rSize = rSize + 6

     oBlock%ADT%ADTree(i)%xMax(:) = oBlock%rBuffer(rSize+1:rSize+6)
     rSize = rSize + 6
  end do

  rSize = rSize + 1
  oBlock%minVol = oBlock%rBuffer(rSize)

  ! Flag this oBlock as being allocated:
  oBlock%allocated = .True.

end subroutine unpackOBlock
