module oversetUtilities
contains
  subroutine emptyFringe(fringe)
    use constants
    use block, only : fringeType
    implicit none
    ! Input/Output
    type(fringeType), intent(inout) :: fringe

    ! Initialize data in empty fringe
    fringe%quality = large
    fringe%origQuality = large
    fringe%donorProc = -1
    fringe%donorBlock = -1
    fringe%dI = -1
    fringe%dJ = -1
    fringe%dK = -1
    fringe%myBlock = -1
    fringe%myI = -1
    fringe%myJ = -1
    fringe%myK = -1
    fringe%donorFrac = -one
    fringe%gInd = -1
    fringe%status = 0
    call setIsCompute(fringe%status, .True.)
  end subroutine emptyFringe

  subroutine printOverlapMatrix(overlap)

    ! This is a debugging routine to print out the overlap matrix.
    use constants
    use communication, only : myid
    use overset, only : CSRMatrix
    implicit none

    ! Input/output
    type(CSRMatrix), intent(in) :: overlap

    ! Working
    integer(kind=intType) :: i, jj

    if (myid == 0) then 
       ! Now dump out who owns what:
       do i=1, overlap%nrow
          write(*, "(a,I4, a)", advance='no'), 'Row:', i, "   "
          do jj=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
             write(*, "(a,I2, a, e10.5)", advance='no'), "(", overlap%colInd(jj), ")", overlap%data(jj)
          end do
          write(*, *) " "
       end do

       print *, '--------------------------------------'
       ! Now dump out who owns what:
       do i=1, overlap%nRow
          write(*, "(a,I4, a)", advance='no'), 'Row:', i, "   "
          do jj=overlap%rowPtr(i), overlap%rowPtr(i+1)-1
             write(*, "(a,I2, a, I8)", advance='no'), "(", overlap%colInd(jj), ")", int(overlap%assignedProc(jj))
          end do
          write(*, *) " "
       end do
    end if
  end subroutine printOverlapMatrix

  subroutine getCumulativeForm(sizeArray, n, cumArray)

    use constants
    implicit none

    ! Input/Output
    integer(kind=intType), dimension(n), intent(in) :: sizeArray
    integer(kind=intType), intent(in) :: n
    integer(kind=intType), dimension(0:n), intent(out) :: cumArray

    ! Working
    integer(kind=intType) :: i

    cumArray(0) = 0
    do i=1, n
       cumArray(i) = cumArray(i-1) + sizeArray(i)
    end do

  end subroutine getCumulativeForm

  subroutine transposeOverlap(A, B)

    ! Create the matrix Create the matrix transpose. 
    ! Inspired by: https://people.sc.fsu.edu/~jburkardt/f_src/sparsekit/sparsekit.f90
    use constants
    use overset, only : CSRMatrix
    implicit none

    ! Input/Output
    type(CSRMatrix), intent(in) :: A
    type(CSRMatrix), intent(inout) :: B

    ! Working
    integer(kind=intType) :: col, colp1, i, k, next

    ! A CSR matrix is the same as a CSC matrix of the transpose. So
    ! essentially the algorithm is convert A as a CSC matrix to B
    ! (a CSR matrix)

    ! Allocate space for everything in B
    B%nnz = A%nnz
    B%nRow = A%nCol
    B%nCol = A%nRow
    allocate(B%data(B%nnz), B%colInd(B%nnz), &
         B%assignedProc(B%nnz), B%rowPtr(B%nRow + 1))
    B%allocated = .True.
    !  Compute lengths of rows of B (ie the columns of A)

    B%rowPtr = 0

    do i = 1, A%nRow
       do k = A%rowPTr(i), A%rowPtr(i+1)-1
          colp1 = A%colInd(k) +1
          B%rowPtr(colp1) = B%rowPtr(colp1) + 1

       end do
    end do
    !
    !  Compute pointers from lengths.
    !
    B%rowPtr(1) = 1
    do i = 1, A%nRow
       B%rowPtr(i+1) = B%rowPtr(i) + B%rowPtr(i+1)
    end do
    !
    !  Do the actual copying.
    !
    do i = 1, A%nRow
       do k = A%rowPtr(i), A%rowPtr(i+1)-1
          col = A%colInd(k)
          next = B%rowPtr(col)

          B%data(next) = A%data(k)
          B%assignedProc(next) = A%assignedProc(k) 

          B%colInd(next) = i 
          B%rowPtr(col) = next + 1
       end do
    end do

    !  Reshift B%rowPtr
    do i = A%nRow, 1, -1
       B%rowPtr(i+1) = B%rowPtr(i)
    end do
    B%rowPtr(1) = 1 

  end subroutine transposeOverlap

  subroutine deallocateCSRMatrix(mat1)

    use constants
    use overset, only : CSRMatrix
    implicit none

    type(CSRMatrix), intent(inout) :: mat1

    if (mat1%allocated) then 
       deallocate(&
            mat1%data, &
            mat1%colInd, &
            mat1%rowPtr, &
            mat1%assignedProc)
    end if

  end subroutine deallocateCSRMatrix

  subroutine computeFringeProcArray(fringes, n, fringeProc, cumFringeProc, nFringeProc)
    ! Compute the breakpoints "cumFringeProc" for a list of sorted n
    ! fringes "fringes". nFringeProc is the total number of unique
    ! processors. fringeProc is the processor number for each section.
    use constants
    use block, only : fringeType
    use overset, only : CSRMatrix
    use communication, only : nProc

    implicit none

    ! Input/Output
    type(fringeType), intent(in), dimension(n) :: fringes
    integer(kind=intType), intent(in) :: n
    integer(kind=intType), intent(out) :: nFringeProc
    integer(kind=intType), intent(out) :: fringeProc(nProc), cumFringeProc(1:nProc+1)

    ! Working
    integer(kind=intType) :: i, currentProc

    fringeProc = -1
    nFringeProc = 0
    cumFringeProc(1) = 1
    currentProc = -1

    do i=1, n
       if (currentProc /= fringes(i)%donorProc) then 
          nFringeProc = nFringeProc + 1
          cumFringeProc(nFringeProc) = i
          fringeProc(nFringeProc) = fringes(i)%donorProc
          currentProc = fringes(i)%donorProc
       end if
    end do

    ! Finally, the nFringeProc+1 entry is always n+1
    cumFringeProc(nFringeProc+1) = n + 1

  end subroutine computeFringeProcArray

  subroutine fracToWeights(frac, weights)
    use constants
    implicit none
    real(kind=realType), intent(in), dimension(3) :: frac
    real(kind=realType), intent(out), dimension(8) :: weights

    weights(1) = (one-frac(1))*(one-frac(2))*(one-frac(3))
    weights(2) = (    frac(1))*(one-frac(2))*(one-frac(3))
    weights(3) = (one-frac(1))*(    frac(2))*(one-frac(3))
    weights(4) = (    frac(1))*(    frac(2))*(one-frac(3))
    weights(5) = (one-frac(1))*(one-frac(2))*(    frac(3))
    weights(6) = (    frac(1))*(one-frac(2))*(    frac(3))
    weights(7) = (one-frac(1))*(    frac(2))*(    frac(3))
    weights(8) = (    frac(1))*(    frac(2))*(    frac(3))
  end subroutine fracToWeights


  subroutine fracToWeights2(frac, weights)
    use constants
    implicit none
    real(kind=realType), intent(in), dimension(3) :: frac
    real(kind=realType), intent(out), dimension(8) :: weights

    weights(1) = (one-frac(1))*(one-frac(2))*(one-frac(3))
    weights(2) = (    frac(1))*(one-frac(2))*(one-frac(3))
    weights(3) = (    frac(1))*(    frac(2))*(one-frac(3))
    weights(4) = (one-frac(1))*(    frac(2))*(one-frac(3))

    weights(5) = (one-frac(1))*(one-frac(2))*(    frac(3))
    weights(6) = (    frac(1))*(one-frac(2))*(    frac(3))
    weights(7) = (    frac(1))*(    frac(2))*(    frac(3))
    weights(8) = (one-frac(1))*(    frac(2))*(    frac(3))

  end subroutine fracToWeights2


  subroutine deallocateOBlocks(oBlocks, n)

    ! This subroutine deallocates all data stores in a list of oBlocks
    use constants
    use adtBuild, only : destroySerialHex
    use overset, only : oversetBlock
    implicit none

    ! Input Params
    type(oversetBlock), dimension(n), intent(inout) :: oBLocks
    integer(kind=intType) :: n

    ! Working Parameters
    integer(kind=intType) :: i

    do i=1, n

       ! oBlock:
       if (oblocks(i)%allocated) then 
          deallocate(&
               oBlocks(i)%hexaConn, &
               oBlocks(i)%globalCell, &
               oBLocks(i)%nearWall, &
               oBLocks(i)%invalidDonor, &
               oBlocks(i)%qualDonor, &
               oBlocks(i)%xADT)
          if (allocated(oblocks(i)%rbuffer)) then 
             deallocate(oBlocks(i)%rBuffer, &
                  oBlocks(i)%iBuffer)
          end if
          call destroySerialHex(oBlocks(i)%ADT)
       end if
    end do
  end subroutine deallocateOBlocks

  subroutine deallocateOFringes(oFringes, n)

    ! This subroutine deallocates all data stores in a list of oFringes
    use constants
    use overset, only : oversetFringe
    implicit none

    ! Input Params
    type(oversetFringe), dimension(n), intent(inout) :: oFringes
    integer(kind=intType) :: n

    ! Working Parameters
    integer(kind=intType) :: i

    do i=1, n
       if (oFringes(i)%allocated) then 
          deallocate(&
               oFringes(i)%x, &
               oFringes(i)%quality, &
               oFringes(i)%origQuality, &
               oFringes(i)%myBlock, &
               oFringes(i)%myIndex, &
               oFringes(i)%donorProc, &
               oFringes(i)%donorBlock, &
               oFringes(i)%dI, &
               oFringes(i)%dJ, &
               oFringes(i)%dK, &
               oFringes(i)%donorFrac, &
               oFringes(i)%gInd, &
               oFringes(i)%isWall, &
               oFringes(i)%xSeed, &
               oFringes(i)%wallInd)
          if (allocated(oFringes(i)%rbuffer)) then 
             deallocate(oFringes(i)%rBuffer,&
                  oFringes(i)%iBuffer)
          end if
       end if
       oFringes(i)%allocated = .False. 
    end do
  end subroutine deallocateOFringes

  subroutine deallocateOSurfs(oSurfs, n)

    ! This subroutine deallocates all data stores in a list of oSurfs

    use constants
    use adtBuild, only : destroySerialQuad
    use overset, only : oversetWall
    use kdtree2_module, only : kdtree2destroy
    implicit none

    ! Input Params
    type(oversetWall), dimension(n), intent(inout) :: oSurfs
    integer(kind=intType) :: n

    ! Working Parameters
    integer(kind=intType) :: i

    do i=1, n
       if (oSurfs(i)%allocated) then 
          deallocate(&
               oSurfs(i)%x, &
               oSurfs(i)%conn, &
               oSurfs(i)%iblank, &
               oSurfs(i)%cellPtr)
          call destroySerialQuad(oSurfs(i)%ADT)
          if (oSurfs(i)%nNodes > 0) then 
             call kdtree2destroy(oSurfs(i)%tree)
          end if
       end if
       oSurfs(i)%allocated = .False.
    end do
  end subroutine deallocateOSurfs

  subroutine wallsOnBlock(wallsPresent) 

    use constants
    use blockPointers, only : nBkGlobal
    use cgnsGrid, only : cgnsDoms
    implicit none

    logical, intent(out) :: wallsPresent
    integer(kind=intType) :: mm
    wallsPresent = .False.
    ! Check THE ORIGINAL CGNS blocks for BCs, because the block may have
    ! been split. 
    do mm=1, cgnsDoms(nbkGlobal)%nBocos
       if (cgnsDoms(nbkGlobal)%bocoInfo(mm)%BCType == NSWallAdiabatic .or. &
            cgnsDoms(nbkGlobal)%bocoInfo(mm)%BCType == NSWallIsothermal .or. &
            cgnsDoms(nbkGlobal)%bocoInfo(mm)%BCType == EulerWall) then
          wallsPresent = .True.
       end if
    end do
  end subroutine wallsOnBlock

  subroutine flagForcedReceivers(tmp)

    use constants
    use blockPointers, only : nx, ny, nz, ie, je, ke, BCData, BCFaceID, nBocos, BCType
    implicit none

    ! This is generic routine for filling up a 3D array of 1st level halos
    ! cells (1:ie, 1:je, 1:ke) indicating cells that are forced
    ! receivers. BlockPointers must have already been set.

    integer(kind=intType), intent(out), dimension(1:ie, 1:je, 1:ke) :: tmp
    integer(kind=intType) :: i, j, k, mm, iStart, iEnd, jStart, jEnd, kStart, kEnd

    tmp = 0
    do mm=1,nBocos
       ! Just record the ranges necessary and we'll add in a generic
       ! loop. Why is it the first three? Well, the first level of halos
       ! off of an overset outer bound is completely
       ! meaningless. Essentially we ignore those. So the outer two
       ! layers of cells are indices 2 and 3. Therefore the first 3 on
       ! either side need to be flagged as invalid.

       select case (BCFaceID(mm))
       case (iMin)
          iStart=1; iEnd=3;
          jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (iMax)
          iStart=nx; iEnd=ie;
          jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (jMin)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=1; jEnd=3;
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (jMax)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=ny; jEnd=je;
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (kMin)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
          kStart=1; kEnd=3;
       case (kMax)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
          kStart=nz; kEnd=ke;
       end select

       if (BCType(mm) == OversetOuterBound) then
          do k=kStart, kEnd
             do j=jStart, jEnd
                do i=iStart, iEnd
                   tmp(i, j, k) = 1
                end do
             end do
          end do
       end if
    end do
  end subroutine flagForcedReceivers


  ! Utility function for unpacking/accessing the status variable

  function isDonor(i)
    use constants
    implicit none
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
    integer(kind=intType), intent(in) :: i
    call getStatus(i, isDonor, isDonor, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  end function isDonor

  function isHole(i)
    use constants
    implicit none
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
    integer(kind=intType), intent(in) :: i
    call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  end function isHole

  function isCompute(i)
    use constants
    implicit none
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
    integer(kind=intType), intent(in) :: i
    call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  end function isCompute

  function isFloodSeed(i)
    use constants
    implicit none
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
    integer(kind=intType), intent(in) :: i
    call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  end function isFloodSeed

  function isFlooded(i)
    use constants
    implicit none
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
    integer(kind=intType), intent(in) :: i
    call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  end function isFlooded

  function isWall(i)
    use constants
    implicit none
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
    integer(kind=intType), intent(in) :: i
    call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  end function isWall

  function isWallDonor(i)
    use constants
    implicit none
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
    integer(kind=intType), intent(in) :: i
    call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  end function isWallDonor

  subroutine setIsDonor(i, flag)
    use constants
    implicit none
    integer(kind=intType), intent(inout) :: i
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor, flag
    call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
    call setStatus(i, flag   , isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  end subroutine setIsDonor

  subroutine setIsHole(i, flag)
    use constants
    implicit none
    integer(kind=intType), intent(inout) :: i
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor, flag
    call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
    call setStatus(i, isDonor, flag  , isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
  end subroutine setIsHole

  subroutine setIsCompute(i, flag)
    use constants
    implicit none
    integer(kind=intType), intent(inout) :: i
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor, flag
    call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
    call setStatus(i, isDonor, isHole, flag     , isFloodSeed, isFlooded, isWall, isWallDonor)
  end subroutine setIsCompute

  subroutine setIsFloodSeed(i, flag)
    use constants
    implicit none
    integer(kind=intType), intent(inout) :: i
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor, flag
    call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
    call setStatus(i, isDonor, isHole, isCompute, flag       , isFlooded, isWall, isWallDonor)
  end subroutine setIsFloodSeed

  subroutine setIsFlooded(i, flag)
    use constants
    implicit none
    integer(kind=intType), intent(inout) :: i
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor, flag
    call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
    call setStatus(i, isDonor, isHole, isCompute, isFloodSeed, flag     , isWall, isWallDonor)
  end subroutine setIsFlooded

  subroutine setIsWall(i, flag)
    use constants
    implicit none
    integer(kind=intType), intent(inout) :: i
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor, flag
    call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
    call setStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, flag  , isWallDonor)
  end subroutine setIsWall

  subroutine setIsWallDonor(i, flag)
    use constants
    implicit none
    integer(kind=intType), intent(inout) :: i
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor, flag
    call getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)
    call setStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, flag)
  end subroutine setIsWallDonor

  subroutine setStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)

    use constants
    implicit none
    integer(kind=intType), intent(out) :: i
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
    i = 0

    if (isDonor  )   i = i + 1
    if (isHole   )   i = i + 2
    if (isCompute)   i = i + 4
    if (isFloodSeed) i = i + 8
    if (isFlooded  ) i = i + 16
    if (isWall     ) i = i + 32
    if (isWallDonor) i = i + 64

  end subroutine setStatus

  subroutine getStatus(i, isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor)

    use constants
    implicit none
    logical :: isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall, isWallDonor
    integer(kind=intType) :: i, j
    j = i

    isDonor = .False.
    isHole = .False.
    isCompute = .False.
    isFloodSeed = .False.
    isFlooded = .False.
    isWall = .False.
    isWallDonor = .False.

    if (j/64 > 0) then 
       isWallDonor = .True. 
       j = j - 64
    end if

    if (j/32 > 0) then 
       isWall = .True. 
       j = j - 32
    end if

    if (j/16 > 0) then 
       isFlooded = .True. 
       j = j - 16
    end if

    if (j/8 > 0) then 
       isFloodSeed = .True. 
       j = j - 8
    end if

    if (j/4 > 0) then 
       isCompute = .True. 
       j = j - 4
    end if

    if (j/2 > 0) then 
       isHole = .True. 
       j = j - 2
    end if

    if (j/1 > 0) then 
       isDonor = .True. 
       j = j - 1
    end if
  end subroutine getStatus

  !
  subroutine binSearchNodes(arr, searchNode, nn, searchInd)

    !  binSearchNodes does binary search for a node 'searchNode'  
    !  in arr(1:nn) and returns index 'searchInd' where           
    !  'searchNode' lies in arr. searchInd = -1 if not found.     

    use constants
    implicit none

    ! Input parameters
    integer(kind=intType), intent(in) :: nn, searchNode
    integer(kind=intType), intent(out) :: searchInd
    integer(kind=intType), intent(in) :: arr(nn)

    ! Local variables
    integer(kind=intType) :: first, last, middle

    first = 1
    last = nn

    middle = (first+last)/2

    do while (first <= last)
       if (arr(middle) < searchNode) then
          first = middle + 1
       else if (arr(middle) == searchNode) then
          searchInd = middle
          exit
       else
          last = middle -1
       end if

       middle = (first+last)/2
    end do !while

    if (first > last) then
       searchInd = -1
       print*, ' binSearchNode fails for searchNode ',searchNode
       STOP 
    end if
  end subroutine binSearchNodes

  subroutine binSearchPocketEdgeType(arr, search, nn, searchInd)

    !  binSearchPocketEdgeType does binary searche for            
    !  pocketEdgeType 'search' edge and returns index 'searchInd' 
    !  where 'search' lies in arr.                                

    use constants
    use overset ! cannot use only becuase of <= operator
    implicit none

    ! Input parameters
    integer(kind=intType), intent(in) :: nn
    integer(kind=intType), intent(out) :: searchInd

    type(pocketEdge), intent(in) :: search
    type(pocketEdge), dimension(*), intent(in) :: arr

    ! Local variables
    integer(kind=intType) :: first, last, middle

    first = 1
    last = nn

    middle = (first+last)/2

    do while (first <= last)
       if (arr(middle) < search) then
          first = middle + 1
       else if (arr(middle) == search) then
          searchInd = middle
          exit
       else
          last = middle -1
       end if

       middle = (first+last)/2
    end do !while

    if (first > last) then
       print*, ' binSearchPocketEdgeType fails for Edge with nodes ',&
            search%n1, search%n2
       STOP 
    end if
  end subroutine binSearchPocketEdgeType

  !
  subroutine qsortEdgeType(arr, nn)
    !
    !       qsortEdgeType sorts the given number of oversetString master   
    !       Edges in increasing order based on the <= operator for this    
    !       derived data type.                                             
    !       (Generously copied from qsortFringeType.F90)                   
    !
    use constants
    use overset ! cannot use only becuase of <= operator
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn

    type(oversetEdge), dimension(*), intent(inout) :: arr
    !
    !      Local variables.
    !
    integer(kind=intType), parameter :: m = 7

    integer(kind=intType) :: nStack
    integer(kind=intType) :: i, j, k, r, l, jStack, ii

    integer :: ierr

    type(oversetEdge) :: a, tmp

    integer(kind=intType), allocatable, dimension(:) :: stack
    integer(kind=intType), allocatable, dimension(:) :: tmpStack

    ! Allocate the memory for stack.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                         &
         call terminate("qsortEdgeType", &
         "Memory allocation failure for stack")

    ! Initialize the variables that control the sorting.

    jStack = 0
    l      = 1
    r      = nn

    ! Start of the algorithm

    do

       ! Check for the size of the subarray.

       if((r-l) < m) then

          ! Perform insertion sort

          do j=l+1,r
             a = arr(j)
             do i=(j-1),l,-1
                if(arr(i) <= a) exit
                arr(i+1) = arr(i)
             enddo
             arr(i+1) = a
          enddo

          ! In case there are no more elements on the stack, exit from
          ! the outermost do-loop. Algorithm has finished.

          if(jStack == 0) exit

          ! Pop stack and begin a new round of partitioning.

          r = stack(jStack)
          l = stack(jStack-1)
          jStack = jStack - 2

       else

          ! Subarray is larger than the threshold for a linear sort.
          ! Choose median of left, center and right elements as
          ! partitioning element a.
          ! Also rearrange so that (l) <= (l+1) <= (r).

          k = (l+r)/2
          tmp      = arr(k)             ! Swap the elements
          arr(k)   = arr(l+1)           ! k and l+1.
          arr(l+1) = tmp

          if(arr(r) < arr(l)) then
             tmp    = arr(l)             ! Swap the elements
             arr(l) = arr(r)             ! r and l.
             arr(r) = tmp
          endif

          if(arr(r) < arr(l+1)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(r)           ! r and l+1.
             arr(r)   = tmp
          endif

          if(arr(l+1) < arr(l)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(l)           ! l and l+1.
             arr(l)   = tmp
          endif

          ! Initialize the pointers for partitioning.

          i = l+1
          j = r
          a = arr(l+1)

          ! The innermost loop

          do

             ! Scan up to find element >= a.
             do
                i = i+1
                if(a <= arr(i)) exit
             enddo

             ! Scan down to find element <= a.
             do
                j = j-1
                if(arr(j) <= a) exit
             enddo

             ! Exit the loop in case the pointers i and j crossed.

             if(j < i) exit

             ! Swap the element i and j.

             tmp    = arr(i)
             arr(i) = arr(j)
             arr(j) = tmp
          enddo

          ! Swap the entries j and l+1. Remember that a equals
          ! arr(l+1).

          arr(l+1) = arr(j)
          arr(j)   = a

          ! Push pointers to larger subarray on stack,
          ! process smaller subarray immediately.

          jStack = jStack + 2
          if(jStack > nStack) then

             ! Storage of the stack is too small. Reallocate.

             allocate(tmpStack(nStack), stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortEdgeType", &
                  "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortEdgeType", &
                  "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortEdgeType", &
                  "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortEdgeType", &
                  "Deallocation error for tmpStack")
          endif

          if((r-i+1) >= (j-l)) then
             stack(jStack)   = r
             r               = j-1
             stack(jStack-1) = j
          else
             stack(jStack)   = j-1
             stack(jStack-1) = l
             l               = j
          endif

       endif
    enddo

    ! Release the memory of stack.

    deallocate(stack, stat=ierr)
    if(ierr /= 0)                         &
         call terminate("qsortEdgeType", &
         "Deallocation error for stack")

    ! Check in debug mode whether the array is really sorted.

    if( debug ) then
       do i=1,(nn-1)
          if(arr(i+1) < arr(i))                 &
               call terminate("qsortEdgeType", &
               "Array is not sorted correctly")
       enddo
    endif

  end subroutine qsortEdgeType

  subroutine qsortFringeType(arr, nn)
    !
    !       qsortFringeListTy sorts the given number of fringes            
    !       increasing order based on the <= operator for this derived     
    !       data type.                                                     
    !
    use constants
    use block ! Cannot use-only becuase of <= operator
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn

    type(fringeType), dimension(*), intent(inout) :: arr
    !
    !      Local variables.
    !
    integer(kind=intType), parameter :: m = 7

    integer(kind=intType) :: nStack
    integer(kind=intType) :: i, j, k, r, l, jStack, ii

    integer :: ierr

    type(fringeType) :: a, tmp

    integer(kind=intType), allocatable, dimension(:) :: stack
    integer(kind=intType), allocatable, dimension(:) :: tmpStack

    ! Allocate the memory for stack.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                         &
         call terminate("qsortfringeType", &
         "Memory allocation failure for stack")

    ! Initialize the variables that control the sorting.

    jStack = 0
    l      = 1
    r      = nn

    ! Start of the algorithm

    do

       ! Check for the size of the subarray.

       if((r-l) < m) then

          ! Perform insertion sort

          do j=l+1,r
             a = arr(j)
             do i=(j-1),l,-1
                if(arr(i) <= a) exit
                arr(i+1) = arr(i)
             enddo
             arr(i+1) = a
          enddo

          ! In case there are no more elements on the stack, exit from
          ! the outermost do-loop. Algorithm has finished.

          if(jStack == 0) exit

          ! Pop stack and begin a new round of partitioning.

          r = stack(jStack)
          l = stack(jStack-1)
          jStack = jStack - 2

       else

          ! Subarray is larger than the threshold for a linear sort.
          ! Choose median of left, center and right elements as
          ! partitioning element a.
          ! Also rearrange so that (l) <= (l+1) <= (r).

          k = (l+r)/2
          tmp      = arr(k)             ! Swap the elements
          arr(k)   = arr(l+1)           ! k and l+1.
          arr(l+1) = tmp

          if(arr(r) < arr(l)) then
             tmp    = arr(l)             ! Swap the elements
             arr(l) = arr(r)             ! r and l.
             arr(r) = tmp
          endif

          if(arr(r) < arr(l+1)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(r)           ! r and l+1.
             arr(r)   = tmp
          endif

          if(arr(l+1) < arr(l)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(l)           ! l and l+1.
             arr(l)   = tmp
          endif

          ! Initialize the pointers for partitioning.

          i = l+1
          j = r
          a = arr(l+1)

          ! The innermost loop

          do

             ! Scan up to find element >= a.
             do
                i = i+1
                if(a <= arr(i)) exit
             enddo

             ! Scan down to find element <= a.
             do
                j = j-1
                if(arr(j) <= a) exit
             enddo

             ! Exit the loop in case the pointers i and j crossed.

             if(j < i) exit

             ! Swap the element i and j.

             tmp    = arr(i)
             arr(i) = arr(j)
             arr(j) = tmp
          enddo

          ! Swap the entries j and l+1. Remember that a equals
          ! arr(l+1).

          arr(l+1) = arr(j)
          arr(j)   = a

          ! Push pointers to larger subarray on stack,
          ! process smaller subarray immediately.

          jStack = jStack + 2
          if(jStack > nStack) then

             ! Storage of the stack is too small. Reallocate.

             allocate(tmpStack(nStack), stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortfringeType", &
                  "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortfringeType", &
                  "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortfringeType", &
                  "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortfringeType", &
                  "Deallocation error for tmpStack")
          endif

          if((r-i+1) >= (j-l)) then
             stack(jStack)   = r
             r               = j-1
             stack(jStack-1) = j
          else
             stack(jStack)   = j-1
             stack(jStack-1) = l
             l               = j
          endif

       endif
    enddo

    ! Release the memory of stack.

    deallocate(stack, stat=ierr)
    if(ierr /= 0)                         &
         call terminate("qsortfringeType", &
         "Deallocation error for stack")

    ! Check in debug mode whether the array is really sorted.

    if( debug ) then
       do i=1,(nn-1)
          if(arr(i+1) < arr(i))                 &
               call terminate("qsortfringeType", &
               "Array is not sorted correctly")
       enddo
    endif

  end subroutine qsortFringeType

  !
  subroutine qsortPocketEdgeType(arr, nn)
    !
    !       qsortPocketEdgeType sorts the given number of oversetString    
    !       master Edges in increasing order based on the <= operator for  
    !       this derived data type.                                        
    !       (Generously copied from qsortFringeType.F90)                   
    !
    use constants
    use overset ! Cannot use-only becuase of <= operator
    use utils, onlY : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn

    type(pocketEdge), dimension(*), intent(inout) :: arr
    !
    !      Local variables.
    !
    integer(kind=intType), parameter :: m = 7

    integer(kind=intType) :: nStack
    integer(kind=intType) :: i, j, k, r, l, jStack, ii

    integer :: ierr

    type(pocketEdge) :: a, tmp

    integer(kind=intType), allocatable, dimension(:) :: stack
    integer(kind=intType), allocatable, dimension(:) :: tmpStack

    ! Allocate the memory for stack.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                         &
         call terminate("qsortEdgeType", &
         "Memory allocation failure for stack")

    ! Initialize the variables that control the sorting.

    jStack = 0
    l      = 1
    r      = nn

    ! Start of the algorithm

    do

       ! Check for the size of the subarray.

       if((r-l) < m) then

          ! Perform insertion sort

          do j=l+1,r
             a = arr(j)
             do i=(j-1),l,-1
                if(arr(i) <= a) exit
                arr(i+1) = arr(i)
             enddo
             arr(i+1) = a
          enddo

          ! In case there are no more elements on the stack, exit from
          ! the outermost do-loop. Algorithm has finished.

          if(jStack == 0) exit

          ! Pop stack and begin a new round of partitioning.

          r = stack(jStack)
          l = stack(jStack-1)
          jStack = jStack - 2

       else

          ! Subarray is larger than the threshold for a linear sort.
          ! Choose median of left, center and right elements as
          ! partitioning element a.
          ! Also rearrange so that (l) <= (l+1) <= (r).

          k = (l+r)/2
          tmp      = arr(k)             ! Swap the elements
          arr(k)   = arr(l+1)           ! k and l+1.
          arr(l+1) = tmp

          if(arr(r) < arr(l)) then
             tmp    = arr(l)             ! Swap the elements
             arr(l) = arr(r)             ! r and l.
             arr(r) = tmp
          endif

          if(arr(r) < arr(l+1)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(r)           ! r and l+1.
             arr(r)   = tmp
          endif

          if(arr(l+1) < arr(l)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(l)           ! l and l+1.
             arr(l)   = tmp
          endif

          ! Initialize the pointers for partitioning.

          i = l+1
          j = r
          a = arr(l+1)

          ! The innermost loop

          do

             ! Scan up to find element >= a.
             do
                i = i+1
                if(a <= arr(i)) exit
             enddo

             ! Scan down to find element <= a.
             do
                j = j-1
                if(arr(j) <= a) exit
             enddo

             ! Exit the loop in case the pointers i and j crossed.

             if(j < i) exit

             ! Swap the element i and j.

             tmp    = arr(i)
             arr(i) = arr(j)
             arr(j) = tmp
          enddo

          ! Swap the entries j and l+1. Remember that a equals
          ! arr(l+1).

          arr(l+1) = arr(j)
          arr(j)   = a

          ! Push pointers to larger subarray on stack,
          ! process smaller subarray immediately.

          jStack = jStack + 2
          if(jStack > nStack) then

             ! Storage of the stack is too small. Reallocate.

             allocate(tmpStack(nStack), stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortEdgeType", &
                  "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortEdgeType", &
                  "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortEdgeType", &
                  "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortEdgeType", &
                  "Deallocation error for tmpStack")
          endif

          if((r-i+1) >= (j-l)) then
             stack(jStack)   = r
             r               = j-1
             stack(jStack-1) = j
          else
             stack(jStack)   = j-1
             stack(jStack-1) = l
             l               = j
          endif

       endif
    enddo

    ! Release the memory of stack.

    deallocate(stack, stat=ierr)
    if(ierr /= 0)                         &
         call terminate("qsortEdgeType", &
         "Deallocation error for stack")

    ! Check in debug mode whether the array is really sorted.

    if( debug ) then
       do i=1,(nn-1)
          if(arr(i+1) < arr(i))                 &
               call terminate("qsortEdgeType", &
               "Array is not sorted correctly")
       enddo
    endif

  end subroutine qsortPocketEdgeType


  subroutine checkOverset (level, sps, totalOrphans)

    !
    !       CheckOverset checks the integrity of the overset connectivity  
    !       and holes. For every comptue cell (iblank = 1) it checks that  
    !       every cell in its stencil are not blanked. If even 1 cell is   
    !      * found with an incomplete stencil it is a fatal error. 
    use constants
    use blockPointers, only : il, jl, kl, iblank, flowDoms, nDom, orphans, &
         iBegOR, jBegOr, kBegOr, nbkGlobal, nOrphans
    use stencils, only : visc_drdw_stencil, N_visc_drdw
    use communication, only : myid, adflow_comm_world
    use utils, only : setPointers, EChk
    implicit none

    ! Input/Output
    integer(kind=intType), intent(in) :: level, sps
    integer(kind=intType), intent(out) :: totalOrphans

    ! Working
    integer(kind=intType) :: i, j, k, nn, ii, jj, kk, n, ierr
    integer(kind=intType) :: magic, localOrphans, i_stencil

    magic = 33
    localOrphans = 0
    do nn=1, nDom
       call setPointers(nn, level, sps)

       ! On the first pass count up the total number of orphans for this block
       n = 0
       do k=2, kl
          do j=2, jl
             do i=2, il
                if (iblank(i,j,k) == 0 .or. iblank(i,j,k)==-2 .or. iblank(i,j,k)==-3) then 

                   stencilLoop: do i_stencil=1, N_visc_drdw
                      ii = visc_drdw_stencil(i_stencil, 1) + i
                      jj = visc_drdw_stencil(i_stencil, 2) + j
                      kk = visc_drdw_stencil(i_stencil, 3) + k

                      if (ii >= 2 .and. jj >= 2 .and. kk>=2 .and. &
                           ii <= il .and. jj <= jl .and. kk <= kl) then 

                         if (iBlank(ii, jj, kk) == 1) then 
                            ! This cell is an orphan:
                            print *,'Error in connectivity: ', nbkglobal, i+iBegOR-1, j+jBegOr-1, k+kBegOr-1
                            n = n + 1
                         end if
                      end if
                   end do stencilLoop
                end if
             end do
          end do
       end do

       localOrphans = localOrphans + n

       ! Remove any existing orphans 
       if (associated(flowDoms(nn, level, sps)%orphans)) then 
          deallocate(flowDoms(nn, level, sps)%orphans)
       end if
       allocate(flowDoms(nn, level, sps)%orphans(3, n))

       ! Save the total number of orphans on this block
       flowDoms(nn, level, sps)%nOrphans = n

       ! Manual set information from blockPointers that would be set
       ! with setPointers()
       orphans => flowDoms(nn, level, sps)%orphans
       nOrphans = n

       ! On the first pass count up the total number of orphans for this block
       n = 0
       do k=2, kl
          do j=2, jl
             do i=2, il
                if (iblank(i,j,k) == 0 .or. iblank(i,j,k)==-2 .or. iblank(i,j,k)==-3) then 

                   stencilLoop2: do i_stencil=1, N_visc_drdw
                      ii = visc_drdw_stencil(i_stencil, 1) + i
                      jj = visc_drdw_stencil(i_stencil, 2) + j
                      kk = visc_drdw_stencil(i_stencil, 3) + k

                      if (ii >= 2 .and. jj >= 2 .and. kk>=2 .and. &
                           ii <= il .and. jj <= jl .and. kk <= kl) then 

                         if (iBlank(ii, jj, kk) == 1) then 
                            ! This cell is an orphan:
                            n = n + 1
                            orphans(:, n) = (/ii, jj, kk/)
                         end if
                      end if
                   end do stencilLoop2
                end if
             end do
          end do
       end do
    end do

    ! Determine the total number of overset orphans
    call mpi_allreduce(localOrphans, totalOrphans, 1, adflow_integer, MPI_SUM, &
         adflow_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    if (myid == 0) then 
       print *, 'Total number of orphans:', totalOrphans
    end if


  end subroutine checkOverset

  !
  !       Finds quadratic uvw weights [-1:1] for interpolant point       
  !       xSearch by solving for its co-ord in donor element defined by  
  !       xElem using newton-raphson iteration.                          
  !       Fractions uvwQuadratic come with initial guess from ADT search 
  !       used for linear interpolation.                                 

  subroutine computeQuadraticWeights(xSearch,xElem,uvwQuadratic)

    use constants

    implicit none

    ! Input variables
    real(kind=realType), intent(in)    :: xSearch(3), xElem(3, -1:1, -1:1, -1:1)
    real(kind=realType), intent(inout) :: uvwQuadratic(3)

    ! Working variables
    ! -----------------------------------------------------------------
    ! newton related
    integer(kind=intType) :: n, niter
    real(kind=realType)   :: resid,residtol
    real(kind=realType)   :: B(3,4)

    ! others
    integer(kind=intType) :: j, l, iii, jjj, kkk
    real(kind=realType)   :: ff(27), shp(3,3), psi(3)
    real(kind=realType)   :: dff(27,3),dshp(3,3,3) !-> differentials wrt 3 psi dirs
    logical               :: ok_flag
    ! -----------------------------------------------------------------

    ! Initialize newton parameters
    niter = 0 
    resid = 1.0e10
    residtol = 1.0e-15


    ! Begin newton iterations to solve for uvw weights wrt quadratic element


    newton_loop: do while (niter < 10 .and. resid > residtol)

       !step 1: find weights
       !--------------------

       ! Initialize weights 
       psi(1:3) = uvwQuadratic(1:3)

       ! Precopute the FE shape functions for each j-direction
       do j=1,3
          shp(1, j) = half*psi(j)*(psi(j) - one)
          shp(2, j) = -(psi(j)**2-1)
          shp(3, j) = half*psi(j)*(psi(j) + one)
       end do

       ! These are the 27 quadratic weights
       ff(1 )   = shp(1, 1)*shp(1, 2)*shp(1, 3)
       ff(2 )   = shp(2, 1)*shp(1, 2)*shp(1, 3)
       ff(3 )   = shp(3, 1)*shp(1, 2)*shp(1, 3)

       ff(4 )   = shp(1, 1)*shp(2, 2)*shp(1, 3)
       ff(5 )   = shp(2, 1)*shp(2, 2)*shp(1, 3)
       ff(6 )   = shp(3, 1)*shp(2, 2)*shp(1, 3)

       ff(7 )   = shp(1, 1)*shp(3, 2)*shp(1, 3)
       ff(8 )   = shp(2, 1)*shp(3, 2)*shp(1, 3)
       ff(9 )   = shp(3, 1)*shp(3, 2)*shp(1, 3)

       ff(10)   = shp(1, 1)*shp(1, 2)*shp(2, 3)
       ff(11)   = shp(2, 1)*shp(1, 2)*shp(2, 3)
       ff(12)   = shp(3, 1)*shp(1, 2)*shp(2, 3)

       ff(13)   = shp(1, 1)*shp(2, 2)*shp(2, 3)
       ff(14)   = shp(2, 1)*shp(2, 2)*shp(2, 3)
       ff(15)   = shp(3, 1)*shp(2, 2)*shp(2, 3)

       ff(16)   = shp(1, 1)*shp(3, 2)*shp(2, 3)
       ff(17)   = shp(2, 1)*shp(3, 2)*shp(2, 3)
       ff(18)   = shp(3, 1)*shp(3, 2)*shp(2, 3)

       ff(19)   = shp(1, 1)*shp(1, 2)*shp(3, 3)
       ff(20)   = shp(2, 1)*shp(1, 2)*shp(3, 3)
       ff(21)   = shp(3, 1)*shp(1, 2)*shp(3, 3)

       ff(22)   = shp(1, 1)*shp(2, 2)*shp(3, 3)
       ff(23)   = shp(2, 1)*shp(2, 2)*shp(3, 3)
       ff(24)   = shp(3, 1)*shp(2, 2)*shp(3, 3)

       ff(25)   = shp(1, 1)*shp(3, 2)*shp(3, 3)
       ff(26)   = shp(2, 1)*shp(3, 2)*shp(3, 3)
       ff(27)   = shp(3, 1)*shp(3, 2)*shp(3, 3)

       !step 2: find differentials of weights
       !-------------------------------------

       ! Linearize the FE shape functions wrt each psi(j) direction
       ! Note: only derivatives wrt psi(j) for any j are non-zero, rest are zero

       dshp(:, :, :) = 0.d0
       do j=1,3
          dshp(1, j, j) = half*(psi(j) - one) + half*psi(j)
          dshp(2, j, j) = -2.d0*psi(j)
          dshp(3, j, j) = half*(psi(j) + one) + half*psi(j)
       end do

       ! Linearize 27 quadratic weights wrt each psi(j) dir, build from dshp

       loop_psi_j: do j=1,3
          dff(1, j) =  dshp(1, 1, j)* shp(1, 2   )* shp(1, 3   ) &
               +  shp(1, 1   )*dshp(1, 2, j)* shp(1, 3   ) &
               +  shp(1, 1   )* shp(1, 2   )*dshp(1, 3, j)

          dff(2, j) =  dshp(2, 1, j)* shp(1, 2   )* shp(1, 3   ) &
               +  shp(2, 1   )*dshp(1, 2, j)* shp(1, 3   ) &
               +  shp(2, 1   )* shp(1, 2   )*dshp(1, 3, j)

          dff(3, j) =  dshp(3, 1, j)* shp(1, 2   )* shp(1, 3   ) &
               +  shp(3, 1   )*dshp(1, 2, j)* shp(1, 3   ) &
               +  shp(3, 1   )* shp(1, 2   )*dshp(1, 3, j)

          dff(4, j) =  dshp(1, 1, j)* shp(2, 2   )* shp(1, 3   ) &
               +  shp(1, 1   )*dshp(2, 2, j)* shp(1, 3   ) &
               +  shp(1, 1   )* shp(2, 2   )*dshp(1, 3, j)

          dff(5, j) =  dshp(2, 1, j)* shp(2, 2   )* shp(1, 3   ) &
               +  shp(2, 1   )*dshp(2, 2, j)* shp(1, 3   ) &
               +  shp(2, 1   )* shp(2, 2   )*dshp(1, 3, j)

          dff(6, j) =  dshp(3, 1, j)* shp(2, 2   )* shp(1, 3   ) &
               +  shp(3, 1   )*dshp(2, 2, j)* shp(1, 3   ) &
               +  shp(3, 1   )* shp(2, 2   )*dshp(1, 3, j)

          dff(7, j) =  dshp(1, 1, j)* shp(3, 2   )* shp(1, 3   ) &
               +  shp(1, 1   )*dshp(3, 2, j)* shp(1, 3   ) &
               +  shp(1, 1   )* shp(3, 2   )*dshp(1, 3, j)

          dff(8, j) =  dshp(2, 1, j)* shp(3, 2   )* shp(1, 3   ) &
               +  shp(2, 1   )*dshp(3, 2, j)* shp(1, 3   ) &
               +  shp(2, 1   )* shp(3, 2   )*dshp(1, 3, j)

          dff(9, j) =  dshp(3, 1, j)* shp(3, 2   )* shp(1, 3   ) &
               +  shp(3, 1   )*dshp(3, 2, j)* shp(1, 3   ) &
               +  shp(3, 1   )* shp(3, 2   )*dshp(1, 3, j)

          dff(10, j) =  dshp(1, 1, j)* shp(1, 2   )* shp(2, 3   ) &
               +  shp(1, 1   )*dshp(1, 2, j)* shp(2, 3   ) &
               +  shp(1, 1   )* shp(1, 2   )*dshp(2, 3, j)

          dff(11, j) =  dshp(2, 1, j)* shp(1, 2   )* shp(2, 3   ) &
               +  shp(2, 1   )*dshp(1, 2, j)* shp(2, 3   ) &
               +  shp(2, 1   )* shp(1, 2   )*dshp(2, 3, j)

          dff(12, j) =  dshp(3, 1, j)* shp(1, 2   )* shp(2, 3   ) &
               +  shp(3, 1   )*dshp(1, 2, j)* shp(2, 3   ) &
               +  shp(3, 1   )* shp(1, 2   )*dshp(2, 3, j)

          dff(13, j) =  dshp(1, 1, j)* shp(2, 2   )* shp(2, 3   ) &
               +  shp(1, 1   )*dshp(2, 2, j)* shp(2, 3   ) &
               +  shp(1, 1   )* shp(2, 2   )*dshp(2, 3, j)

          dff(14, j) =  dshp(2, 1, j)* shp(2, 2   )* shp(2, 3   ) &
               +  shp(2, 1   )*dshp(2, 2, j)* shp(2, 3   ) &
               +  shp(2, 1   )* shp(2, 2   )*dshp(2, 3, j)

          dff(15, j) =  dshp(3, 1, j)* shp(2, 2   )* shp(2, 3   ) &
               +  shp(3, 1   )*dshp(2, 2, j)* shp(2, 3   ) &
               +  shp(3, 1   )* shp(2, 2   )*dshp(2, 3, j)

          dff(16, j) =  dshp(1, 1, j)* shp(3, 2   )* shp(2, 3   ) &
               +  shp(1, 1   )*dshp(3, 2, j)* shp(2, 3   ) &
               +  shp(1, 1   )* shp(3, 2   )*dshp(2, 3, j)

          dff(17, j) =  dshp(2, 1, j)* shp(3, 2   )* shp(2, 3   ) &
               +  shp(2, 1   )*dshp(3, 2, j)* shp(2, 3   ) &
               +  shp(2, 1   )* shp(3, 2   )*dshp(2, 3, j)

          dff(18, j) =  dshp(3, 1, j)* shp(3, 2   )* shp(2, 3   ) &
               +  shp(3, 1   )*dshp(3, 2, j)* shp(2, 3   ) &
               +  shp(3, 1   )* shp(3, 2   )*dshp(2, 3, j)

          dff(19, j) =  dshp(1, 1, j)* shp(1, 2   )* shp(3, 3   ) &
               +  shp(1, 1   )*dshp(1, 2, j)* shp(3, 3   ) &
               +  shp(1, 1   )* shp(1, 2   )*dshp(3, 3, j)

          dff(20, j) =  dshp(2, 1, j)* shp(1, 2   )* shp(3, 3   ) &
               +  shp(2, 1   )*dshp(1, 2, j)* shp(3, 3   ) &
               +  shp(2, 1   )* shp(1, 2   )*dshp(3, 3, j)

          dff(21, j) =  dshp(3, 1, j)* shp(1, 2   )* shp(3, 3   ) &
               +  shp(3, 1   )*dshp(1, 2, j)* shp(3, 3   ) &
               +  shp(3, 1   )* shp(1, 2   )*dshp(3, 3, j)

          dff(22, j) =  dshp(1, 1, j)* shp(2, 2   )* shp(3, 3   ) &
               +  shp(1, 1   )*dshp(2, 2, j)* shp(3, 3   ) &
               +  shp(1, 1   )* shp(2, 2   )*dshp(3, 3, j)

          dff(23, j) =  dshp(2, 1, j)* shp(2, 2   )* shp(3, 3   ) &
               +  shp(2, 1   )*dshp(2, 2, j)* shp(3, 3   ) &
               +  shp(2, 1   )* shp(2, 2   )*dshp(3, 3, j)

          dff(24, j) =  dshp(3, 1, j)* shp(2, 2   )* shp(3, 3   ) &
               +  shp(3, 1   )*dshp(2, 2, j)* shp(3, 3   ) &
               +  shp(3, 1   )* shp(2, 2   )*dshp(3, 3, j)

          dff(25, j) =  dshp(1, 1, j)* shp(3, 2   )* shp(3, 3   ) &
               +  shp(1, 1   )*dshp(3, 2, j)* shp(3, 3   ) &
               +  shp(1, 1   )* shp(3, 2   )*dshp(3, 3, j)

          dff(26, j) =  dshp(2, 1, j)* shp(3, 2   )* shp(3, 3   ) &
               +  shp(2, 1   )*dshp(3, 2, j)* shp(3, 3   ) &
               +  shp(2, 1   )* shp(3, 2   )*dshp(3, 3, j)

          dff(27, j) =  dshp(3, 1, j)* shp(3, 2   )* shp(3, 3   ) &
               +  shp(3, 1   )*dshp(3, 2, j)* shp(3, 3   ) &
               +  shp(3, 1   )* shp(3, 2   )*dshp(3, 3, j)

       end do loop_psi_j


       ! Step 3: construct Jacobian d(R(psi(:))/d(psi(:)) and residue R(psi(:))
       ! ----------------------------------------------------------------------

       ! loop over x,y,z dirs, stored row-wise
       loop_xyz: do n=1,3 

          ! construct LHS -d(R(psi(:))/d(psi(:))

          ! loop over each psi dir j
          do j=1,3

             ! loop over nodes
             l = 0
             b(n, j) = 0.d0

             do kkk=-1,1
                do jjj=-1,1
                   do iii=-1,1
                      l = l + 1
                      b(n,j) = b(n,j) + dff(l,j) * xElem(n, iii, jjj, kkk)
                   end do
                end do
             end do
          end do !j

          ! construct RHS (Xp - sum_i(ff_i * X_i), i =1,27) for nth row (x,y or z)
          l = 0
          b(n, 4) = xSearch(n) 

          ! loop over nodes
          do kkk=-1,1
             do jjj=-1,1
                do iii=-1,1
                   l = l + 1
                   b(n,4) = b(n,4) - ff(l) * xElem(n, iii, jjj, kkk)
                end do
             end do
          end do

       end do loop_xyz

       ! Get d(uvwQuadratic) weights
       ! invert 3x3 matrix returns solution in b(:,4) 
       call matrixinv3by3(b, ok_flag)
       if (.not. ok_flag) stop 'Can not invert B in computeQuadraticWeights'

       ! update uvwQuadratic weights
       do n=1,3
          uvwQuadratic(n) = uvwQuadratic(n) + b(n,4)

          ! sanity check: if weights lie outside [-1:1] reset to 0.5
          if( (uvwQuadratic(n)-1)*(uvwQuadratic(n)+1) > 0.d0) uvwQuadratic(n) = half
       end do

       ! L2-norm of residue
       resid = 0.d0
       do n=1,3
          resid = resid + b(n,4)**2
       end do

       niter = niter +1

    end do newton_loop
  end subroutine computeQuadraticWeights
  !
  !      Plain invert of A to solve Ax=B                                
  !      a(3,3) --- matrix to be inverted                               
  !      a(:,4) --- contains the right hand side vector, also stores    
  !                 the final solution                                  

  subroutine matrixinv3by3(a, ok_flag)

    use constants
    implicit none

    ! Input variables
    real(kind=realType), intent(inout) :: a(3, 4)
    logical,             intent(out)   :: ok_flag   

    ! Working variables
    integer(kind=intType) :: n
    real(kind=realType) :: det, cofactor(3,3), ainv(3, 3), rin(3,1), rout(3,1), myeps

    myeps = 1.0e-10
    ainv = 0.d0

    det =   a(1,1)*a(2,2)*a(3,3)  &
         - a(1,1)*a(2,3)*a(3,2)  &
         - a(1,2)*a(2,1)*a(3,3)  &
         + a(1,2)*a(2,3)*a(3,1)  &
         + a(1,3)*a(2,1)*a(3,2)  &
         - a(1,3)*a(2,2)*a(3,1)


    if (abs(det) <= myeps) then
       ainv = 0.0d0
       ok_flag = .false.
       return
    end if

    cofactor(1,1) = +(a(2,2)*a(3,3)-a(2,3)*a(3,2))
    cofactor(1,2) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
    cofactor(1,3) = +(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    cofactor(2,1) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
    cofactor(2,2) = +(a(1,1)*a(3,3)-a(1,3)*a(3,1))
    cofactor(2,3) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
    cofactor(3,1) = +(a(1,2)*a(2,3)-a(1,3)*a(2,2))
    cofactor(3,2) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
    cofactor(3,3) = +(a(1,1)*a(2,2)-a(1,2)*a(2,1))

    ainv = transpose(cofactor) / det

    ok_flag = .true.

    do n=1, 3
       rin(n,1) = a(n, 4)
    end do

    !save solution on a(:,4)
    rout = matmul(ainv,rin)

    do n=1, 3
       a(n, 4) = rout(n,1)
    end do

  end subroutine matrixinv3by3

  subroutine fringeReduction(level, sps)

    use constants
    use blockPointers, only : nDom, il, jl, kl, fringes
    use stencils, only : visc_drdw_stencil, n_visc_drdw
    use utils, only : setPointers
    implicit none

    ! Input/Output
    integer(kind=intType), intent(in) :: level, sps

    ! Working
    integer(kind=intType) :: i, j, k, nn, ii, jj, kk, i_stencil
    logical :: computeCellFound

    do nn=1, nDom
       call setPointers(nn, level, sps)

       do k=2, kl
          do j=2, jl
             do i=2, il

                ! Check if this cell is a fringe:
                if (fringes(i, j, k)%donorProc /= -1) then 

                   computeCellFound = .False.

                   stencilLoop2: do i_stencil=1, N_visc_drdw
                      ii = visc_drdw_stencil(i_stencil, 1) + i
                      jj = visc_drdw_stencil(i_stencil, 2) + j
                      kk = visc_drdw_stencil(i_stencil, 3) + k

                      if (isCompute(fringes(ii, jj, kk)%status)) then 
                         ! This is a compute cell
                         computeCellFound = .True.
                      end if
                   end do stencilLoop2

                   if (.not. computeCellFound) then 
                      ! This cell is a hole no compute cell
                      ! surrounding a fringe, we can hard iblank it.
                      call emptyFringe(fringes(i, j, k))
                      call setIsHole(fringes(i, j, k)%status, .True.)
                      call setIsCompute(fringes(i, j, k)%status, .False.)
                   end if
                end if
             end do
          end do
       end do
    end do

  end subroutine fringeReduction

  subroutine irregularCellCorrection(level, sps)

    use constants
    use blockPointers, only : nDom, il, jl, kl, fringes
    use utils, only : setPointers
    implicit none

    ! Input/Output
    integer(kind=intType), intent(in) :: level, sps

    ! Working
    integer(kind=intType) :: i, j, k, nn

    do nn=1, nDom
       call setPointers(nn, level, sps)

       do k=2, kl
          do j=2, jl
             do i=2, il
                if (isDonor(fringes(i, j, k)%status) .and. &
                     fringes(i, j, k)%donorProc /= -1) then 
                   ! Clear this fringe
                   call emptyFringe(fringes(i, j, k))
                end if
             end do
          end do
       end do
    end do


  end subroutine irregularCellCorrection


  function checkOversetPresent()

    ! This routine determines if there are any oveset boundaries
    ! present in the mesh.

    use constants
    use blockPointers, only : nDom, nBocos, BCType
    use communication, only : adflow_comm_world
    use utils, only : setPointers, EChk
    implicit none

    ! Function
    logical :: checkOversetPresent, local

    ! Working
    integer(Kind=intType) :: nn, mm, ierr

    local = .False.
    do nn=1, nDom
       call setPointers(nn, 1_intType, 1_intType)

       do mm=1, nBocos
          if (BCType(mm) == OversetOuterBound) then
             local = .True.
          end if
       end do
    end do

    call mpi_allreduce(local, checkOversetPresent, 1, MPI_LOGICAL, MPI_LOR, ADflow_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

  end function checkOversetPresent

  subroutine setIblankArray(level, sps)

    use constants
    use blockPointers, only : nDom, il, jl, kl, fringes, iblank, flowDoms
    use communication, only : myid, commPatternCell_2nd, internalCell_2nd,&
         adflow_comm_world
    use utils, only : setPointers, EChk
    use haloExchange, only : whalo1to1IntGeneric
    implicit none

    ! Input/Output
    integer(kind=intType), intent(in) :: level, sps

    ! Working
    integer(kind=intType) :: i, j, k, nn
    integer(kind=intType) :: nCompute, nFringe, nBlank, nFloodSeed, nFlooded
    integer(kind=intType) :: counts(5), ierr
    nCompute = 0
    nFringe = 0
    nBlank = 0
    nFloodSeed = 0
    nFlooded = 0

    do nn=1, nDom
       call setPointers(nn, level, sps)
       iBlank = 1
       do k=2, kl
          do j=2, jl
             do i=2, il

                if (fringes(i, j, k)%donorProc /= -1) then
                   iblank(i, j, k) = -1
                   nFringe = nFringe + 1

                else if (isFloodSeed(fringes(i, j, k)%status)) then
                   iBlank(i, j, k) = -3
                   nFloodSeed = nFloodSeed + 1

                else if (isFlooded(fringes(i, j, k)%status)) then
                   iBlank(i, j, k) = -2
                   nFlooded = nFlooded + 1

                else if (isHole(fringes(i, j, k)%status)) then 
                   iBlank(i, j, k) = 0
                   nBlank = nBlank + 1

                else
                   ! Compute cell
                   nCompute = nCompute + 1
                end if

             end do
          end do
       end do
    end do

    ! Update the iblank info. 
    domainLoop:do nn=1, nDom
       flowDoms(nn, level, sps)%intCommVars(1)%var => &
            flowDoms(nn, level, sps)%iblank(:, :, :)
    end do domainLoop

    ! Run the generic integer exchange
    call wHalo1to1IntGeneric(1, level, sps, commPatternCell_2nd, internalCell_2nd)

    call mpi_reduce((/nCompute, nFringe, nBlank, nFlooded, nFloodSeed/), &
         counts, 5, adflow_integer, MPI_SUM, 0, adflow_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    if (myid == 0) then 
       print *, '+--------------------------------+'
       print *, '| Compute   Cells:', counts(1)
       print *, '| Fringe    Cells:', counts(2)
       print *, '| Blanked   Cells:', counts(3)
       print *, '| Flooded   Cells:', counts(4)
       print *, '| FloodSeed Cells:', counts(5)
       print *, '+--------------------------------+' 
    end if
  end subroutine setIblankArray

  subroutine dumpIblank(level, sps)

    use constants
    use blockPointers, only: il, jl, kl, x, nDom, iBlank
    use communication, only : myID
    use utils, only : setPointers
    implicit none

    ! Input/Output
    integer(kind=intType), intent(in) :: level, sps

    ! Working
    integer(kind=intType) :: i, j, k, nn
    real(kind=realType) :: xp(3)
    character(80) :: fileName

    write (fileName,"(a,I2.2,a)") "proc_", myid, ".dat"
    open(unit=19,file=trim(fileName),form='formatted')

    do nn=1, nDom
       call setPointers(nn, level, sps)
       do k=2, kl
          do j=2, jl
             do i=2, il

                ! Compute the cell center:
                xp = eighth*(&
                     x(i-1, j-1, k-1, :) + &
                     x(i  , j-1, k-1, :) + &
                     x(i-1, j  , k-1, :) + &
                     x(i  , j  , k-1, :) + &
                     x(i-1, j-1, k  , :) + &
                     x(i  , j-1, k  , :) + &
                     x(i-1, j  , k  , :) + &
                     x(i  , j  , k  , :))

                write(19, "(E18.10, E18.10, E18.10, I3)"), xp(1), xp(2), xp(3), iblank(i, j, k)
             end do
          end do
       end do
    end do

    close(19)

  end subroutine dumpIblank

  subroutine getWorkArray(overlap, work)

    use constants
    use communication, only : myid
    use overset, only : CSRMatrix, nDomTotal
    implicit none

    ! Input/Output
    type(CSRMatrix), intent(in) :: overlap
    integer(kind=intType), dimension(:,:), allocatable :: work

    ! Local variables
    integer(kind=intType) :: nWork, jj, iDom, jDom

    nWork = 0
    do jj=1,overlap%nnz
       if (overlap%assignedProc(jj) == myid) then 
          nWork = nWork + 1
       end if
    end do
    allocate(work(4, nWork))
    
    nWork = 0
    do iDom=1, nDomTotal
       do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
          jDom = overlap%colInd(jj)
          if (overlap%assignedProc(jj) == myID) then 
             nWork = nWork + 1
             work(1, nWork) = iDom
             work(2, nWork) = jDom
             work(3, nWork) = jj
             work(4, nWork) = 0
          end if
       end do
    end do
  end subroutine getWorkArray

end module oversetUtilities
