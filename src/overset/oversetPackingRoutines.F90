module oversetPackingRoutines
contains

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

    iSize = iSize + 6 ! Blocks sizes + proc + nn + cluster

    iSize = iSize + 8*il*jl*kl ! hexa conn

    iSize = iSize + (ie+2)*(je+2)*(ke+2) ! global cell

    iSize = iSize + il*jl*kl ! nearWall

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

    use constants
    use overset, only : oversetBlock
    implicit none

    ! Pack up everything we need for this block into its own buffer
    ! inlucding the data required for the ADTree

    ! Input/Output Parameters
    type(oversetBlock), intent(inout) :: oBlock

    ! Working paramters
    integer(kind=intType) :: rSize, iSize, i, j, k, nHexa, nADT
    integer(kind=intType) :: ie, je, ke, il, jl ,kl

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
    oBlock%iBuffer(6) = oBlock%cluster

    iSize = iSize + 6

    il = oBlock%il
    jl = oBlock%jl
    kl = oBlock%kl

    ie = il + 1
    je = jl + 1
    ke = kl + 1

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

    do k=1, kl
       do j=1, jl
          do i=1, il
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

    use constants
    use overset, only : oversetBlock
    implicit none

    ! unPack everything we need for this block from its own buffer
    ! and reconstitute the data required for the ADTree. It is assumed
    ! the buffers are already allocated and the data is available. This
    ! does the exact *OPPOSITE* operation as the packBlock() routine

    ! Input/Output Parameters
    type(oversetBlock), intent(inout) :: oBlock

    ! Working paramters
    integer(kind=intType) :: rSize, iSize, i, j, k, nHexa, nADT
    integer(kind=intType) :: ie, je, ke, il, jl, kl

    ! Reset the integer counter and add all the integers on this pass
    iSize = 0

    oBlock%il = oBlock%iBuffer(1) 
    oBlock%jl = oBlock%iBuffer(2)
    oBlock%kl = oBlock%iBuffer(3)
    oBlock%proc = oBlock%iBuffer(4)
    oBlock%block = oBlock%iBuffer(5)
    oBlock%cluster = oBlock%iBuffer(6)
    iSize = iSize + 6

    il = oBlock%il
    jl = oBlock%jl
    kl = oBlock%kl
    ie = il + 1
    je = jl + 1
    ke = kl + 1

    nHexa = oBlock%il * oBlock%jl * oBlock%kl
    nADT = ie*je*ke

    ! Allocate the remainder of the arrays in oBlock.
    allocate(oBlock%hexaConn(8, nHexa))
    allocate(oBlock%globalCell(0:ie+1, 0:je+1, 0:ke+1))
    allocate(oBlock%nearWall(1:il, 1:jl, 1:kl))
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

    do k=1, kl
       do j=1, jl
          do i=1, il
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
    ! easy, but we use a routine to make it look the same as for the 
    ! oBlock. 

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
    iSize = mm * 3 + 5 ! We need wallInd, isWall, myIndex plus 5 for the sizes
    rSize = mm * 6 ! Need to send x and xSeed (3 each)

  end subroutine getOFringeBufferSizes

  subroutine packOFringe(oFringe)

    use constants
    use overset, only : oversetFringe

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

    mm = (oFringe%nx)*(oFringe%ny)*(oFringe%nz)

    oFringe%iBuffer(1) = oFringe%il
    oFringe%iBuffer(2) = oFringe%jl
    oFringe%iBuffer(3) = oFringe%kl
    oFringe%iBuffer(4) = oFringe%cluster
    oFringe%iBuffer(5) = oFringe%block
    ii = 5

    ! Copy the integers. Just wallInd and isWall
    do i=1, mm

       ii = ii +1
       oFringe%iBuffer(ii) = oFringe%wallInd(i)

       ii = ii +1
       oFringe%iBuffer(ii) = oFringe%isWall(i)

       ii = ii +1
       oFringe%iBuffer(ii) = oFringe%fringeIntBuffer(5, i) ! myIndex

    end do

    ! Copy the reals. Reset the buffer here.
    ii = 0
    do i=1, mm
       oFringe%rBuffer(ii+1) = oFringe%x(1, i)
       oFringe%rBuffer(ii+2) = oFringe%x(2, i)
       oFringe%rBuffer(ii+3) = oFringe%x(3, i)
       oFringe%rBuffer(ii+4) = oFringe%xSeed(1, i)
       oFringe%rBuffer(ii+5) = oFringe%xSeed(2, i)
       oFringe%rBuffer(ii+6) = oFringe%xSeed(3, i)
       ii = ii + 6
    end do

  end subroutine packOFringe

  subroutine unpackOFringe(oFringe)

    use constants
    use overset, only : oversetFringe
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
    oFringe%nx = oFringe%il-1
    oFringe%ny = oFringe%jl-1
    oFringe%nz = oFringe%kl-1
    oFringe%cluster = oFringe%iBuffer(4)
    oFringe%block = oFringe%ibuffer(5)

    mm = (oFringe%nx)*(oFringe%ny)*(oFringe%nz)

    allocate(& 
         oFringe%x(3, mm), &
         oFringe%xSeed(3, mm), &
         oFringe%wallInd(mm), &
         oFringe%isWall(mm))

    ! Assume each cell will get just one donor. It's just a guess, it
    ! will be expanded if necessary so the exact value doesn't matter.
    allocate(oFringe%fringeIntBuffer(5, mm), oFringe%fringeRealBuffer(4, mm))
    oFringe%nDonor = 0

    ii = 5 ! Already copied out the sizes

    ! Copy out integers
    do i=1, mm
       ii = ii + 1
       oFringe%wallInd(i) = oFringe%iBuffer(ii)

       ii = ii + 1
       oFringe%isWall(i) = oFringe%iBuffer(ii)

       ii = ii + 1
       oFringe%fringeIntBuffer(5, i) = oFringe%iBuffer(ii)

       oFringe%fringeIntBuffer(4, i) = oFringe%block
    end do

    ! Copy the reals. Reset the counter ii counter here.
    ii = 0
    do i=1, mm
       oFringe%x(1, i) = oFringe%rBuffer(ii+1)
       oFringe%x(2, i) = oFringe%rBuffer(ii+2)
       oFringe%x(3, i) = oFringe%rBuffer(ii+3)

       oFringe%xSeed(1, i) = oFringe%rBuffer(ii+4)
       oFringe%xSeed(2, i) = oFringe%rBuffer(ii+5)
       oFringe%xSeed(3, i) = oFringe%rBuffer(ii+6)

       ii = ii + 6
    end do

    ! Flag this oFringe as being allocated:
    oFringe%allocated = .True.
    deallocate(oFringe%rBuffer, oFringe%iBuffer)

  end subroutine unpackOFringe

  subroutine getWallSize(famList, nNodes, nCells, dualMesh)
    ! Simple helper routine to return the number of wall nodes and cells
    ! for the block pointed to by blockPointers. 

    use constants
    use blockPointers, only :BCType, nBocos, BCData
    use sorting, only : bsearchIntegers
    implicit none

    ! Input
    integer(kind=intType), intent(in), dimension(:) :: famList
    logical :: dualMesh

    ! Output
    integer(kind=intType), intent(out) :: nNodes, nCells

    ! Working:
    integer(kind=intType) :: mm, iBeg, iEnd, jBeg, jEnd

    ! Figure out the size the wall is going to be.
    nNodes = 0
    nCells = 0
    do mm=1, nBocos
       famInclude: if (bsearchIntegers(BCData(mm)%famID, famList) > 0) then 
          if (dualMesh) then 
             jBeg = BCData(mm)%jnBeg-1 ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg-1 ; iEnd = BCData(mm)%inEnd
          else
             jBeg = BCData(mm)%jnBeg; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd
          end if

          nNodes = nNodes + (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
          nCells = nCells + (iEnd - iBeg )*(jEnd - jBeg)
       end if famInclude
    end do
    
  end subroutine getWallSize

  subroutine getOSurfBufferSizes(famList, il, jl, kl, iSize, rSize, dualMesh)

    ! Subroutine to get the required buffer sizes. This one is pretty
    ! easy, but we use a routine to make it look the same as for hte
    ! oBlock. Note that these bufer sizes are over-estimates. They are
    ! the maximum possible amount of data to send. 

    use constants
    implicit none

    ! Input/OUtput
    integer(kind=intType), intent(in), dimension(:) :: famList
    integer(kind=intType), intent(in) :: il, jl ,kl
    logical, intent(in) :: dualMesh
    integer(kind=intType), intent(out) :: rSize, iSize

    ! Working
    integer(kind=intType) :: mm, nNodes, nCells, nBBox, nLeaves

    ! Initalization
    iSize = 3 ! For the block sizes
    iSize = iSize + 4 ! For the maxCells/nCells/nNodes variables
    rSize = 0

    call getWallSize(famList, nNodes, nCells, dualMesh)

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

       rSize = rSize + nNodes ! surface delta

       rSize = rSize + nBBox*6 ! Cell bounding boxes

       rSize = rSize + nLeaves*12 ! Bounding boxes for leaves
    end if
  end subroutine getOSurfBufferSizes

  subroutine packOSurf(famList, oSurf, dualMesh)

    use constants
    use overset, only : oversetWall

    implicit none

    ! Pack up the search coordines in this oSurf into its own buffer
    ! so we are ready to send it.

    ! Input/Output Parameters
    integer(kind=intType), intent(in), dimension(:) :: famList
    type(oversetWall), intent(inout) :: oSurf
    logical, intent(in) :: dualMesh
    ! Working paramters
    integer(kind=intType) :: rSize, iSize, mm, i, j, nNodes, nCells

    call getOSurfBufferSizes(famList, oSurf%il, oSurf%kl, oSurf%kl, isize, rSize, dualMesh)

    ! Allocate the buffers
    allocate(oSurf%rBuffer(rSize), oSurf%iBuffer(iSize))

    oSurf%iBuffer(1) = oSurf%il
    oSurf%iBuffer(2) = oSurf%jl
    oSurf%iBuffer(3) = oSurf%kl
    oSurf%iBuffer(4) = oSurf%nNodes
    oSurf%iBuffer(5) = oSurf%nCells
    oSurf%iBuffer(6) = oSurf%maxCells
    oSurf%iBuffer(7) = oSurf%cluster

    if (oSurf%nNodes > 0) then 
       iSize = 7
       do j=1, oSurf%nCells
          do i=1, 4
             iSize = iSize + 1
             oSurf%iBuffer(iSize) = oSurf%conn(i, j)
          end do
       end do

       do i=1, oSurf%maxCells
          iSize = iSize + 1
          oSurf%iBuffer(iSize) = oSurf%iBlank(i)
       end do

       do i=1, oSurf%nCells
          iSize = iSize + 1
          oSurf%iBuffer(iSize) = oSurf%cellPtr(i)
       end do

       do i=1, oSurf%ADT%nLeaves
          iSize = iSize + 1
          oSurf%iBuffer(iSize) = oSurf%ADT%ADTree(i)%children(1)
          iSize = iSize + 1
          oSurf%iBuffer(iSize) = oSurf%ADT%ADTree(i)%children(2)
       end do

       ! Done with the integer values, do the real ones
       rSize = 0

       do i=1, oSurf%nNodes
          do j=1, 3
             rSize = rSize + 1
             oSurf%rBuffer(rSize) = oSurf%x(j, i)
          end do
       end do

       do i=1, oSurf%nNodes
          rSize = rSize + 1
          oSurf%rBuffer(rSize) = oSurf%delta(i)
       end do

       do i=1, oSurf%ADT%nBboxes
          oSurf%rBuffer(rSize+1:rSize+6) = oSurf%ADT%xBBox(:, i)
          rSize = rSize + 6
       end do

       do i=1, oSurf%ADT%nLeaves
          oSurf%rBuffer(rSize+1:rSize+6) = oSurf%ADT%ADTree(i)%xMin(:)
          rSize = rSize + 6

          oSurf%rBuffer(rSize+1:rSize+6) = oSurf%ADT%ADTree(i)%xMax(:)
          rSize = rSize + 6
       end do
    end if
  end subroutine packOSurf

  subroutine unpackOSurf(oSurf)

    use constants
    use overset, only : oversetWall
    use kdtree2_module, only : kdtree2_create
    implicit none

    ! Input/Output Parameters
    type(oversetWall), intent(inout) :: oSurf

    ! Working paramters
    integer(kind=intType) :: rSize, iSize, idom, i,  j, k, n, iNode

    ! Set the sizes of this oSurf
    oSurf%il = oSurf%iBuffer(1)
    oSurf%jl = oSurf%iBuffer(2)
    oSurf%kl = oSurf%iBuffer(3)

    oSurf%nNodes = oSurf%iBuffer(4)
    oSurf%nCells = oSurf%iBuffer(5)
    oSurf%maxCells = oSurf%iBuffer(6)
    oSurf%cluster = oSurf%iBuffer(7)

    iSize = 7
    rSize = 0

    ! Allocate the arrays now that we know the sizes
    allocate(oSurf%x(3, oSurf%nNodes))
    allocate(oSurf%delta(oSurf%nNodes))
    allocate(oSurf%conn(4, oSurf%nCells))
    allocate(oSurf%iBlank(oSurf%maxCells))
    allocate(oSurf%cellPtr(oSurf%nCells))
    allocate(oSurf%nte(4, oSurf%nNodes))
    ! Once we know the sizes, allocate all the arrays in the
    ! ADTree. Since we are not going to call the *actual* build routine
    ! for the ADT, we need to set all the information ourselves. This
    ! essentially does the same thing as buildSerialHex.
    oSurf%ADT%adtType = adtSurfaceADT
    oSurf%ADT%nNodes = oSurf%nNodes
    oSurf%ADT%nTetra = 0
    oSurf%ADT%nPyra = 0
    oSurf%ADT%nPrisms = 0
    oSurf%ADT%nTria = 0
    oSurf%ADT%nQuads = oSurf%nCells
    oSurf%ADT%coor => oSurf%x
    oSurf%ADT%quadsConn => oSurf%conn
    nullify(oSurf%ADT%triaConn)
    oSurf%ADT%nBBoxes = oSurf%nCells
    allocate(oSurf%ADT%xBBOX(6, oSurf%nCells))
    allocate(oSurf%ADT%elementType(oSurf%nCells))
    allocate(oSurf%ADT%elementID(oSurf%nCells))
    oSurf%ADT%comm = MPI_COMM_SELF
    oSurf%ADT%nProcs = 1
    oSurf%ADT%myID = 0

    ! All hexas
    oSurf%ADT%elementType = adtQuadrilateral

    do i=1, oSurf%nCells
       oSurf%ADT%elementID(i) = i
    end do

    oSurf%ADT%nLeaves = oSurf%ADT%nBBoxes - 1
    if(oSurf%ADT%nBBoxes <= 1) oSurf%ADT%nLeaves = oSurf%ADT%nLeaves + 1

    allocate(oSurf%ADT%ADTree(oSurf%ADT%nLeaves))

    ! Now continue copying out the values if necessary:
    if (oSurf%nNodes > 0) then 
       do j=1, oSurf%nCells
          do i=1, 4
             iSize = iSize + 1
             oSurf%conn(i, j) = oSurf%iBuffer(iSize)
          end do
       end do

       do i=1, oSurf%maxCells
          iSize = iSize + 1
          oSurf%iBlank(i) = oSurf%iBuffer(iSize)
       end do

       do i=1, oSurf%nCells
          iSize = iSize + 1
          oSurf%cellPtr(i) = oSurf%iBuffer(iSize)
       end do

       do i=1, oSurf%ADT%nLeaves
          iSize = iSize + 1
          oSurf%ADT%ADTree(i)%children(1) = oSurf%iBuffer(iSize)
          iSize = iSize + 1
          oSurf%ADT%ADTree(i)%children(2) = oSurf%iBuffer(iSize)
       end do

       ! Done with the integer values, do the real ones
       rSize = 0

       do i=1, oSurf%nNodes
          do j=1, 3
             rSize = rSize + 1
             oSurf%x(j, i) = oSurf%rBuffer(rSize)
          end do
       end do

       do i=1, oSurf%nNodes
          rSize = rSize + 1
          oSurf%delta(i) = oSurf%rBuffer(rSize)
       end do

       do i=1, oSurf%ADT%nBboxes
          oSurf%ADT%xBBox(:, i) = oSurf%rBuffer(rSize+1:rSize+6)
          rSize = rSize + 6
       end do

       do i=1, oSurf%ADT%nLeaves
          oSurf%ADT%ADTree(i)%xMin(:) = oSurf%rBuffer(rSize+1:rSize+6)
          rSize = rSize + 6

          oSurf%ADT%ADTree(i)%xMax(:) = oSurf%rBuffer(rSize+1:rSize+6)
          rSize = rSize + 6
       end do
    end if

    ! Build the KDTree
    if (oSurf%nNodes > 0) then 
       oSurf%tree => kdtree2_create(oSurf%x)
    end if

    ! Build the inverse of the connectivity, the nodeToElem array. 
    oSurf%nte = 0
    do i=1, oSurf%nCells
       do j=1, 4
          n = oSurf%conn(j, i)
          inner:do k=1,4
             if (oSurf%nte(k, n) == 0) then 
                oSurf%nte(k, n) = i
                exit inner
             end if
          end do inner
       end do
    end do

    ! Flag this oSurf as being allocated:
    oSurf%allocated = .True.
    deallocate(oSurf%rBuffer, oSurf%iBuffer)
  end subroutine unpackOSurf
end module oversetPackingRoutines
