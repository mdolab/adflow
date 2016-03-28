subroutine getOBlockBufferSizes(il, jl, kl, iSize, rSize)

  ! Subroutine to get the required buffer sizes. This only uses the
  ! block dimensions and technically doesn't have anything to do with
  ! oBlock. This allows us to figure out the sizes and perform a
  ! global communication before actually building the ADTrees, which
  ! may not be that well load-balanced. 

  use constants
  implicit none

  ! Input/OUtput
  integer(kind=intType) :: il, jl ,kl
  integer(kind=intType) :: rSize, iSize

  ! Working paramters
  integer(kind=intType) :: ie, je, ke, nBBox, nLeaves
  
  ! Initializeation
  iSize = 0
  rSize = 0

  ! Create ie sizes as well
  ie = il + 1; je = jl + 1; ke = kl+1

  ! Count up the integers we want to send:

  iSize = iSize + 5 ! Blocks sizes + proc + nn

  iSize = iSize + 8*il*jl*kl ! hexa conn

  iSize = iSize + (ie+2)*(je+2)*(ke+2) ! global cell
  
  iSize = iSize + ie*je*ke ! nearWall

  iSize = iSize + ie*je*ke ! invalidDonor
  
  ! Number of boxes in the ADT is the same as the number of elements
  ! (like hexa conn (without the *8 obviously)
  nBBox = il * jl * kl
  nLeaves = nBBox-1 ! See ADT/adtBuild.f90

  iSize = iSize + nLeaves*2 ! Two ints for the children in each leaf

  ! Count up the reals we ned to send:
  rSize = rSize + ie*je*ke ! qualDonor

  rSize = rSize + 3*ie*je*ke ! xADT

  rSize = rSize + nBBox*6 ! Cell bounding boxes

  rSize = rSize + nLeaves*12 ! Bounding boxes for leaves

  rSize = rSize + 1 ! Min block volume
end subroutine getOBlockBufferSizes

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
  integer(kind=intType) :: ie, je, ke

  ! If the buffer is already allocated, the block is packed and there
  ! is nothing to do
  if (allocated(oBlock%rBuffer)) then 
     return
  end if

  call getOBlockBufferSizes(oBlock%il, oBlock%jl, oBlock%kl, iSize, rSize)
 
  ! Allocate the buffers
  allocate(oBlock%rBuffer(rSize), oBlock%iBuffer(iSize))

  ! Reset the integer counter and add all the integers on this pass
  iSize = 0

  oBlock%iBuffer(1) = oBlock%il
  oBlock%iBuffer(2) = oBlock%jl
  oBlock%iBuffer(3) = oBlock%kl
  oBlock%iBuffer(4) = oBlock%proc
  oBlock%iBuffer(5) = oBlock%block

  iSize = iSize + 5
  ie = oBlock%il + 1
  je = oBlock%jl + 1
  ke = oBlock%kl + 1

  nHexa = oBlock%il * oBlock%jl * oBlock%kl
  nADT = ie*je*ke

  do j=1, nHexa
     do i=1, 8
        iSize = iSize + 1
        oBlock%iBuffer(iSize) = oBlock%hexaConn(i, j)
     end do
  end do

  do k=0, ke+1
     do j=0, je+1
        do i=0, ie+1
           iSize = iSize + 1
           oBlock%iBuffer(iSize) = oBlock%globalCell(i, j, k)
        end do
     end do
  end do

  do k=1, ke
     do j=1, je
        do i=1, ie
           iSize = iSize + 1
           oBlock%iBuffer(iSize) = oBlock%nearWall(i, j, k)
        end do
     end do
  end do

  do k=1, ke
     do j=1, je
        do i=1, ie
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

  do i=1, ie*je*ke
     rSize = rSize + 1
     oBlock%rBuffer(rSize) = oBlock%qualDonor(1, i)
  end do

  do j=1, ie*je*ke
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
  integer(kind=intType) :: ie, je, ke

  ! Reset the integer counter and add all the integers on this pass
  iSize = 0

  oBlock%il = oBlock%iBuffer(1) 
  oBlock%jl = oBlock%iBuffer(2)
  oBlock%kl = oBlock%iBuffer(3)
  oBlock%proc = oBlock%iBuffer(4)
  oBlock%block = oBlock%iBuffer(5)
  iSize = iSize + 5

  ie = oBlock%il + 1
  je = oBlock%jl + 1
  ke = oBlock%kl + 1
  nHexa = oBlock%il * oBlock%jl * oBlock%kl
  nADT = ie*je*ke

  ! Allocate the remainder of the arrays in oBlock.
  allocate(oBlock%hexaConn(8, nHexa))
  allocate(oBlock%globalCell(0:ie+1, 0:je+1, 0:ke+1))
  allocate(oBlock%nearWall(1:ie, 1:je, 1:ke))
  allocate(oBlock%invalidDonor(1:ie, 1:je, 1:ke))
  allocate(oBlock%qualDonor(1, ie * je * ke))
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

  do k=0, ke+1
     do j=0, je+1
        do i=0, ie+1
           iSize = iSize + 1
           oBlock%globalCell(i, j, k) = oBlock%iBuffer(iSize)
        end do
     end do
  end do

  do k=1, ke
     do j=1, je
        do i=1, ie
           iSize = iSize + 1
           oBlock%nearWall(i, j, k) = oBlock%iBuffer(iSize)
        end do
     end do
  end do

  do k=1, ke
     do j=1, je
        do i=1, ie
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

  do i=1, ie*je*ke
     rSize = rSize + 1
     oBlock%qualDonor(1, i) =  oBlock%rBuffer(rSize)
  end do

  do j=1, ie*je*ke
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
  deallocate(oBlock%iBuffer, oBlock%rBuffer)

end subroutine unpackOBlock

subroutine getOFringeBufferSizes(il, jl, kl, iSize, rSize)

  ! Subroutine to get the required buffer sizes. This one is pretty
  ! easy, but we use a routine to make it look the same as for hte
  ! oBlock. Note that these bufer sizes are over-estimates. They are
  ! the maximum possible amount of data to send. 

  use constants
  implicit none

  ! Input/OUtput
  integer(kind=intType), intent(in) :: il, jl ,kl
  integer(kind=intType), intent(out) :: rSize, iSize

  ! Working
  integer(kind=intType) :: mm
  
  ! All arrays have the same size
  mm = (il-1)*(jl-1)*(kl-1) ! nx*ny*nz

  ! Initializeation
  iSize = mm * 2 + 3! We need isWall and myIndex plus 3 for the sizes
  rSize = mm * 4 ! Need to send x and quality

end subroutine getOFringeBufferSizes

subroutine packOFringe(oFringe)

  use overset
  use constants
  implicit none

  ! Pack up the search coordines in this oFringe into its own buffer
  ! so we are ready to send it.

  ! Input/Output Parameters
  type(oversetFringe), intent(inout) :: oFringe

  ! Working paramters
  integer(kind=intType) :: rSize, iSize, mm, i, ii

  ! If the buffer is already allocated, the block is packed and there
  ! is nothing to do
  if (allocated(oFringe%rBuffer)) then 
     return
  end if

  call getOFringeBufferSizes(oFringe%il, oFringe%jl, oFringe%kl, &
       iSize, rSize)

  ! Allocate the buffers
  allocate(oFringe%rBuffer(rSize), oFringe%iBuffer(iSize))

  mm = (oFringe%il-1)*(oFringe%jl-1)*(oFringe%kl-1)

  oFringe%iBuffer(1) = oFringe%il
  oFringe%iBuffer(2) = oFringe%jl
  oFringe%iBuffer(3) = oFringe%kl
  ii = 3

  ! Copy the integers. Just isWall and myIndex for sending.
  do i=1, mm
     ii = ii +1
     oFringe%iBuffer(ii) = oFringe%isWall(i)

     ii = ii +1
     oFringe%iBuffer(ii) = oFringe%myIndex(i)
  end do

  ! Copy the reals
  ii = 0
  do i=1, mm
     oFringe%rBuffer(ii+1) = oFringe%x(1, i)
     oFringe%rBuffer(ii+2) = oFringe%x(2, i)
     oFringe%rBuffer(ii+3) = oFringe%x(3, i)
     oFringe%rBuffer(ii+4) = oFringe%quality(i)
     ii = ii + 4
  end do

end subroutine packOFringe

subroutine unpackOFringe(oFringe)
  use communication
  use overset
  use constants
  implicit none

  ! Pack up the search coordines in this oFringe into its own buffer
  ! so we are ready to send it.

  ! Input/Output Parameters
  type(oversetFringe), intent(inout) :: oFringe

  ! Working paramters
  integer(kind=intType) :: rSize, iSize, idom, i, ii, mm

  ! Set the sizes of this oFringe
  oFringe%il = oFringe%iBuffer(1)
  oFringe%jl = oFringe%iBuffer(2)
  oFringe%kl = oFringe%iBuffer(3)

  mm = (oFringe%il-1)*(oFringe%jl-1)*(oFringe%kl-1)

  allocate(& 
       oFringe%x(3, mm), &
       oFringe%quality(mm), &
       oFringe%myBlock(mm), &
       oFringe%myIndex(mm), &
       oFringe%donorProc(mm), &
       oFringe%donorBlock(mm), &
       oFringe%dI(mm), &
       oFringe%dJ(mm), &
       oFringe%dK(mm), &
       oFringe%donorFrac(3, mm), &
       oFringe%gInd(8, mm), &    
       oFringe%isWall(mm))

  ! Initialzize default values
  oFringe%donorProc = -1
  oFringe%dI = -1
  oFringe%dJ = -1
  oFringe%dK = -1
  oFringe%donorFrac = -one
  oFringe%gInd = -1
  oFringe%isWall = 0

  ii = 3 ! Already copied out the sizes
  do i=1, mm
     ii = ii + 1
     oFringe%isWall(i) = oFringe%iBuffer(ii)

     ii = ii + 1
     oFringe%myIndex(i) = oFringe%iBuffer(ii)
  end do

  ! Copy the reals
  ii = 0
  do i=1, mm
     oFringe%x(1, i) = oFringe%rBuffer(ii+1)
     oFringe%x(2, i) = oFringe%rBuffer(ii+2)
     oFringe%x(3, i) = oFringe%rBuffer(ii+3)
     oFringe%quality(i) = oFringe%rBuffer(ii+4)
     ii = ii + 4
  end do

  ! Flag this oFringe as being allocated:
  oFringe%allocated = .True.
  deallocate(oFringe%rBuffer, oFringe%iBuffer)

end subroutine unpackOFringe

subroutine getWallSize(nNodes, nCells, dualMesh)
  ! Simple helper routine to return the number of wall nodes and cells
  ! for the block pointed to by blockPointers. 

  use blockPointers
  use BCTypes
  implicit none
  
  ! Input
  logical :: dualMesh

  ! Output
  integer(kind=intType), intent(out) :: nNodes, nCells

  ! Working:
  integer(kind=intType) :: mm, iBeg, iEnd, jBeg, jEnd

  ! Figure out the size the wall is going to be. See
  ! initializeOWall.F90 for why the sizes are the way they are. 
  nNodes = 0
  nCells = 0
  do mm=1, nBocos
      if (BCType(mm) == EulerWall .or. BCType(mm) == NSWallAdiabatic .or. &
          BCType(mm) == NSWallIsoThermal) then 
         if (dualMesh) then 
            jBeg = BCData(mm)%jnBeg-1 ; jEnd = BCData(mm)%jnEnd
            iBeg = BCData(mm)%inBeg-1 ; iEnd = BCData(mm)%inEnd
         else
            jBeg = BCData(mm)%jnBeg; jEnd = BCData(mm)%jnEnd
            iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd
         end if

         nNodes = nNodes + (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
         nCells = nCells + (iEnd - iBeg )*(jEnd - jBeg)
     end if
  end do

end subroutine getWallSize

subroutine getOWallBufferSizes(il, jl, kl, iSize, rSize, dualMesh)

  ! Subroutine to get the required buffer sizes. This one is pretty
  ! easy, but we use a routine to make it look the same as for hte
  ! oBlock. Note that these bufer sizes are over-estimates. They are
  ! the maximum possible amount of data to send. 

  use constants
  implicit none

  ! Input/OUtput
  integer(kind=intType), intent(in) :: il, jl ,kl
  logical, intent(in) :: dualMesh
  integer(kind=intType), intent(out) :: rSize, iSize

  ! Working
  integer(kind=intType) :: mm, nNodes, nCells, nBBox, nLeaves

  ! Initalization
  iSize = 3 ! For the block sizes
  iSize = iSize + 4 ! For the maxCells/nCells/nNodes variables
  rSize = 0
  
  call getWallSize(nNodes, nCells, dualMesh)

  ! Note that nCells here is the maximum number size. This will result
  ! in a slight overestimate of the buffer size. This is ok.

  if (nNodes > 0) then 
     ! Count up the integers we want to send:

     iSize = iSize + nCells*4 ! This is for the connectivity

     iSize = iSize + nCells   ! This is for the iblank array

     iSize = iSize + nCells   ! This is for the cellPtr array
     
     ! Number of boxes in the ADT is the same as the number of elements
     nBBox = nCells
     nLeaves = nBBox-1 ! See ADT/adtBuild.f90
     
     iSize = iSize + nLeaves*2 ! Two ints for the children in each leaf
     
     ! Count up the reals we ned to send:
     rSize = rSize + 3*nNodes ! surface coordinates

     rSize = rSize + nBBox*6 ! Cell bounding boxes
     
     rSize = rSize + nLeaves*12 ! Bounding boxes for leaves
  end if
end subroutine getOWallBufferSizes

subroutine packOWall(oWall, dualMesh)

  use overset
  use constants
  implicit none

  ! Pack up the search coordines in this oWall into its own buffer
  ! so we are ready to send it.

  ! Input/Output Parameters
  type(oversetWall), intent(inout) :: oWall
  logical, intent(in) :: dualMesh
  ! Working paramters
  integer(kind=intType) :: rSize, iSize, mm, i, j, nNodes, nCells

  call getOWallBufferSizes(oWall%il, oWall%kl, oWall%kl, isize, rSize, dualMesh)

  ! Allocate the buffers
  allocate(oWall%rBuffer(rSize), oWall%iBuffer(iSize))

  oWall%iBuffer(1) = oWalL%il
  oWall%iBuffer(2) = oWalL%jl
  oWall%iBuffer(3) = oWalL%kl
  oWall%iBuffer(4) = oWall%nNodes
  oWall%iBuffer(5) = oWall%nCells
  oWall%iBuffer(6) = oWall%maxCells
  oWall%iBuffer(7) = oWall%cluster
  
  if (oWall%nNodes > 0) then 
     iSize = 7
     do j=1, oWall%nCells
        do i=1, 4
           iSize = iSize + 1
           oWall%iBuffer(iSize) = oWall%conn(i, j)
        end do
     end do

     do i=1, oWall%maxCells
        iSize = iSize + 1
        oWall%iBuffer(iSize) = oWall%iBlank(i)
     end do

     do i=1, oWall%nCells
        iSize = iSize + 1
        oWall%iBuffer(iSize) = oWall%cellPtr(i)
     end do
     
     do i=1, oWall%ADT%nLeaves
        iSize = iSize + 1
        oWall%iBuffer(iSize) = oWall%ADT%ADTree(i)%children(1)
        iSize = iSize + 1
        oWall%iBuffer(iSize) = oWall%ADT%ADTree(i)%children(2)
     end do

     ! Done with the integer values, do the real ones
     rSize = 0
     
     do i=1, oWall%nNodes
        do j=1, 3
           rSize = rSize + 1
           oWall%rBuffer(rSize) = oWall%x(j, i)
        end do
     end do
     
     do i=1, oWall%ADT%nBboxes
        oWall%rBuffer(rSize+1:rSize+6) = oWall%ADT%xBBox(:, i)
        rSize = rSize + 6
     end do
     
     do i=1, oWall%ADT%nLeaves
        oWall%rBuffer(rSize+1:rSize+6) = oWall%ADT%ADTree(i)%xMin(:)
        rSize = rSize + 6
        
        oWall%rBuffer(rSize+1:rSize+6) = oWall%ADT%ADTree(i)%xMax(:)
        rSize = rSize + 6
     end do
  end if
end subroutine packOWall

subroutine unpackOWall(oWall)

  use overset
  use constants
  use kdtree2_module
  implicit none

  ! Input/Output Parameters
  type(oversetWall), intent(inout) :: oWall

  ! Working paramters
  integer(kind=intType) :: rSize, iSize, idom, i,  j, k, n, iNode

  ! Set the sizes of this oWall
  oWall%il = oWall%iBuffer(1)
  oWall%jl = oWall%iBuffer(2)
  oWall%kl = oWall%iBuffer(3)

  oWall%nNodes = oWall%iBuffer(4)
  oWall%nCells = oWall%iBuffer(5)
  oWall%maxCells = oWall%iBuffer(6)
  oWall%cluster = oWall%iBuffer(7)

  iSize = 7
  rSize = 0

  ! Allocate the arrays now that we know the sizes
  allocate(oWall%x(3, oWall%nNodes))
  allocate(oWall%conn(4, oWall%nCells))
  allocate(oWall%iBlank(oWall%maxCells))
  allocate(oWall%cellPtr(oWall%nCells))
  allocate(oWall%nte(4, oWall%nNodes))
  ! Once we know the sizes, allocate all the arrays in the
  ! ADTree. Since we are not going to call the *actual* build routine
  ! for the ADT, we need to set all the information ourselves. This
  ! essentially does the same thing as buildSerialHex.
  oWall%ADT%adtType = adtSurfaceADT
  oWall%ADT%nNodes = oWall%nNodes
  oWall%ADT%nTetra = 0
  oWall%ADT%nPyra = 0
  oWall%ADT%nPrisms = 0
  oWall%ADT%nTria = 0
  oWall%ADT%nQuads = oWall%nCells
  oWall%ADT%coor => oWall%x
  oWall%ADT%quadsConn => oWall%conn
  nullify(oWall%ADT%triaConn)
  oWall%ADT%nBBoxes = oWall%nCells
  allocate(oWall%ADT%xBBOX(6, oWall%nCells))
  allocate(oWall%ADT%elementType(oWall%nCells))
  allocate(oWall%ADT%elementID(oWall%nCells))
  oWall%ADT%comm = MPI_COMM_SELF
  oWall%ADT%nProcs = 1
  oWall%ADT%myID = 0

  ! All hexas
  oWall%ADT%elementType = adtQuadrilateral

  do i=1, oWall%nCells
     oWall%ADT%elementID(i) = i
  end do

  oWall%ADT%nLeaves = oWall%ADT%nBBoxes - 1
  if(oWall%ADT%nBBoxes <= 1) oWall%ADT%nLeaves = oWall%ADT%nLeaves + 1

  allocate(oWall%ADT%ADTree(oWall%ADT%nLeaves))

  ! Now continue copying out the values if necessary:
  if (oWall%nNodes > 0) then 
     do j=1, oWall%nCells
        do i=1, 4
           iSize = iSize + 1
           oWall%conn(i, j) = oWall%iBuffer(iSize)
        end do
     end do

     do i=1, oWall%maxCells
        iSize = iSize + 1
        oWall%iBlank(i) = oWall%iBuffer(iSize)
     end do

     do i=1, oWall%nCells
        iSize = iSize + 1
        oWall%cellPtr(i) = oWall%iBuffer(iSize)
     end do

     do i=1, oWall%ADT%nLeaves
        iSize = iSize + 1
        oWall%ADT%ADTree(i)%children(1) = oWall%iBuffer(iSize)
        iSize = iSize + 1
        oWall%ADT%ADTree(i)%children(2) = oWall%iBuffer(iSize)
     end do

     ! Done with the integer values, do the real ones
     rSize = 0
     
     do i=1, oWall%nNodes
        do j=1, 3
           rSize = rSize + 1
           oWall%x(j, i) = oWall%rBuffer(rSize)
        end do
     end do

     do i=1, oWall%ADT%nBboxes
        oWall%ADT%xBBox(:, i) = oWall%rBuffer(rSize+1:rSize+6)
        rSize = rSize + 6
     end do
     
     do i=1, oWall%ADT%nLeaves
        oWall%ADT%ADTree(i)%xMin(:) = oWall%rBuffer(rSize+1:rSize+6)
        rSize = rSize + 6
        
        oWall%ADT%ADTree(i)%xMax(:) = oWall%rBuffer(rSize+1:rSize+6)
        rSize = rSize + 6
     end do
  end if

  ! Build the KDTree
  if (oWall%nNodes > 0) then 
     oWall%tree => kdtree2_create(oWall%x)
  end if

  ! Build the inverse of the connectivity, the nodeToElem array. 
  oWall%nte = 0
  do i=1, oWall%nCells
     do j=1, 4
        n = oWall%conn(j, i)
        inner:do k=1,4
           if (oWall%nte(k, n) == 0) then 
              oWall%nte(k, n) = i
              exit inner
           end if
        end do inner
     end do
  end do

  ! Flag this oWall as being allocated:
  oWall%allocated = .True.
  deallocate(oWall%rBuffer, oWall%iBuffer)
end subroutine unpackOWall
