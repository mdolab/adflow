module oversetInitialization

contains

  subroutine initializeStatus(level, sps)

    ! This subroutine initializes the status variable for use in the
    ! overset process

    use constants
    use blockPointers, only : nDom, il, jl, kl, ib, jb, kb, status
    use oversetUtilities, only : setIsCompute
    use utils, only : setPointers
    implicit none

    ! Input params
    integer(kind=intType), intent(in) :: level, sps

    ! Working parameters
    integer(kind=intType) :: nn, i, j, k

    do nn=1, nDom
       call setPointers(nn, level, sps)

       ! Now loop over the owned cells and set the isCompute flag to
       ! true. 
       
       do k=2, kl
          do j=2, jl
             do i=2, il
                status(i,j,k) = 0
                call setIsCompute(status(i, j, k), .True. )
             end do
          end do
       end do
    end do

  end subroutine initializeStatus
  
  subroutine reInitializeStatus(level, sps)

    ! This subroutine reinitializes the status variable. However, if
    ! cell is a wallDonor, that information is kept. 

    use constants
    use blockPointers, only : nDom, il, jl, kl, ib, jb, kb, status
    use oversetUtilities, only : setIsCompute, isWallDonor, setIsWallDonor
    use utils, only : setPointers
    implicit none

    ! Input params
    integer(kind=intType), intent(in) :: level, sps

    ! Working parameters
    integer(kind=intType) :: nn, i, j, k
    logical :: wDonor
    do nn=1, nDom
       call setPointers(nn, level, sps)

       ! Now loop over the owned cells and set the isCompute flag to
       ! true. 
       
       do k=2, kl
          do j=2, jl
             do i=2, il
                wDonor = isWallDonor(status(i, j, k))
                status(i, j, k) = 0
                call setIsCompute(status(i, j, k), .True. )
                call setIsWallDonor(status(i, j, k), wDonor)
             end do
          end do
       end do
    end do
  end subroutine reInitializeStatus


  subroutine initializeOBlock(oBlock, nn, level, sps)

    ! This routine allocates the data for the supplied oBlock using the
    !  data currently in blockPointers
    use constants
    use overset, only : oversetBlock, clusters, cumDomProc
    use inputOverset, only : backgroundVolScale, nearWallDist
    use blockPointers, only : x, globalCell, il, jl, kl, ib, jb, kb, &
         ie, je, ke, vol, iBlank, xSeed, forcedRecv
    use adtBuild, only : buildSerialHex
    use communication, only : myID
    use stencils, only : visc_drdw_stencil, n_visc_drdw
    use utils, only : mynorm2
    use oversetUtilities, only :  wallsOnBlock
    implicit none 

    ! Input Params
    type(oversetBlock), intent(inout) :: oBlock
    integer(kind=intType) :: nn, level, sps, kk

    ! Working paramters
    integer(kind=intType) :: i, j, k, mm, nADT, nHexa, planeOffset
    integer(kind=intType) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
    real(kind=realType) :: factor, frac,  dist, xp(3)
    integer(kind=intType) :: i_stencil, ii, jj, iii
    logical :: wallsPresent
    logical, allocatable, dimension(:, :, :)  :: nearWallTmp

    ! Set all the sizes for this block.
    oBlock%il = il
    oBlock%jl = jl
    oBlock%kl = kl

    oBlock%proc = myID
    oBlock%block = nn
    oBlock%cluster = clusters(cumDomProc(myid) + nn)
    call wallsOnBlock(wallsPresent)

    ! Do the reset of the allocs
    allocate( &
         oBlock%qualDonor(1, ie*je*ke), &
         oBlock%globalCell(0:ib, 0:jb, 0:kb), &
         oBlock%invalidDonor(1:ie, 1:je, 1:ke))

    ! Invalid Donor array is simply if the cell is a forced receiver
    ! or not.
    do k=1,ke
       do j=1,je
          do i=1,ie
             oBlock%invalidDonor(i,j,k) = forcedRecv(i,j,k)
          end do
       end do
    end do

    ! Compute the qualDonor depending on if we have a wall block or not. 
    mm = 0
    do k=1,ke
       do j=1,je
          do i=1,ie
             mm = mm + 1
             if (wallsPresent) then 
                oBlock%qualDonor(1, mm) = vol(i, j, k)**third 
             else
                oBlock%qualDonor(1, mm) = (backGroundVolScale*vol(i, j, k))**third
             end if
          end do
       end do
    end do

    !Copy over global cell
    oBlock%globalCell = globalCell

    ! Now setup the data for the ADT
    nHexa = il * jl * kl
    nADT = ie * je * ke

    allocate(oBlock%xADT(3, nADT), oBlock%hexaConn(8, nHexa))
    ! Fill up the xADT using cell centers (dual mesh)
    mm = 0

    ! Allocate the nearWall
    allocate(oBlock%nearWall(1:il, 1:jl, 1:kl))
    oBlock%nearWall = 0

    allocate(nearWallTmp(1:ie, 1:je, 1:ke))
    nearWallTmp = .False.

    do k=1, ke
       do j=1, je
          do i=1, ie
             mm = mm + 1
             xp = eighth*(&
                  x(i-1, j-1, k-1, :) + &
                  x(i  , j-1, k-1, :) + &
                  x(i-1, j  , k-1, :) + &
                  x(i  , j  , k-1, :) + &
                  x(i-1, j-1, k  , :) + &
                  x(i  , j-1, k  , :) + &
                  x(i-1, j  , k  , :) + &
                  x(i  , j  , k  , :))
             oBlock%xADT(:, mm) = xp

             ! Determine if this point is near wall. Note that the
             ! boundary halos sill have xSeed as "large" so these won't
             ! be flagged as nearWall. We will account for this below. 
             dist = mynorm2(xp - xSeed(i, j, k, :))
             if (dist < nearWallDist) then 
                nearWallTmp(i, j, k) = .True. 
             end if
          end do
       end do
    end do

    ! Now finally set the nearwall for the dual mesh cells. It is
    ! considered a near wall if all "nodes" of the dual mesh cell are
    ! also near wall. Have to be carful not to count boundary halos
    ! since they do not have nearWallTmp Values.

    do k=1, kl
       do j=1, jl
          do i=1, il
             if (&
                  (nearWallTmp(i  , j  , k  ) .or. globalCell(i  , j,   k  ) < 0) .and. &
                  (nearWallTmp(i+1, j  , k  ) .or. globalCell(i+1, j,   k  ) < 0) .and. &
                  (nearWallTmp(i  , j+1, k  ) .or. globalCell(i  , j+1, k  ) < 0) .and. &
                  (nearWallTmp(i+1, j+1, k  ) .or. globalCell(i+1, j+1, k  ) < 0) .and. &
                  (nearWallTmp(i  , j  , k+1) .or. globalCell(i  , j,   k+1) < 0) .and. &
                  (nearWallTmp(i+1, j  , k+1) .or. globalCell(i+1, j,   k+1) < 0) .and. &
                  (nearWallTmp(i  , j+1, k+1) .or. globalCell(i  , j+1, k+1) < 0) .and. &
                  (nearWallTmp(i+1, j+1, k+1) .or. globalCell(i+1, j+1, k+1) < 0)) then 
                oBlock%nearWall(i, j, k) = 1
             end if
          end do
       end do
    end do

    deallocate(nearWallTmp)
    mm = 0
    ! These are the 'elements' of the dual mesh.
    planeOffset = ie * je
    do k=2, ke
       do j=2, je
          do i=2, ie
             mm = mm + 1
             oBlock%hexaConn(1, mm) = (k-2)*planeOffset + (j-2)*ie + (i-2) + 1
             oBlock%hexaConn(2, mm) = oBlock%hexaConn(1, mm) + 1 
             oBlock%hexaConn(3, mm) = oBlock%hexaConn(2, mm) + ie
             oBlock%hexaConn(4, mm) = oBlock%hexaConn(3, mm) - 1 

             oBlock%hexaConn(5, mm) = oBlock%hexaConn(1, mm) + planeOffset
             oBlock%hexaConn(6, mm) = oBlock%hexaConn(2, mm) + planeOffset
             oBlock%hexaConn(7, mm) = oBlock%hexaConn(3, mm) + planeOffset
             oBlock%hexaConn(8, mm) = oBlock%hexaConn(4, mm) + planeOffset
          end do
       end do
    end do

    ! Call the custom build routine -- Serial only, only Hexa volumes,
    ! we supply our own ADT Type

    call buildSerialHex(nHexa, nADT, oBlock%xADT, oBlock%hexaConn, oBlock%ADT)

    ! Flag this block as being allocated
    oBlock%allocated = .True.

  end subroutine initializeOBlock

  subroutine initializeOFringes(oFringe, nn)

    ! This subroutine initializes the fringe information for the given
    ! block, level and spectral instance. It is assumed that
    ! blockPointers are already set. 
    use constants
    use communication, only : myID
    use blockPointers
    use overset, only : oversetFringe, clusters, cumDomProc
    use stencils, only : visc_drdw_stencil, N_visc_drdw
    use inputOverset, only :  backgroundVolScale
    use utils, only : isWallType
    use oversetUtilities, only :  wallsOnBlock, windIndex, unwindindex
    implicit none

    ! Input Params
    type(oversetFringe), intent(inout) :: oFringe
    integer(kind=intType), intent(in) :: nn

    ! Working Params
    integer(kind=intTYpe) :: i, j, k, mm, iDim, ii, jj, kk, iii, jjj, myI, myJ, myK
    integer(kind=intTYpe) :: iStart, iEnd, jStart, jEnd, kStart, kEnd
    logical :: wallsPresent
    integer(kind=intType) :: i_stencil
    real(kind=realType) :: dist, frac, xp(3)
    ! Check if we have walls:
    call wallsOnBLock(wallsPresent)

    ! Set the sizes for the oFringe and allocate the required space. 
    oFringe%il = il
    oFringe%jl = jl
    oFringe%kl = kl
    oFringe%nx = nx
    oFringe%ny = ny
    oFringe%nz = nz
    oFringe%block = nn
    oFringe%cluster = clusters(cumDomProc(myid) + nn) 
    oFringe%proc = myid

    mm = nx*ny*nz
    allocate(oFringe%x(3, mm))
    allocate(oFringe%xSeed(3, mm))
    allocate(oFringe%wallInd(mm))
    allocate(oFringe%isWall(mm))
    oFringe%isWall = 0
    oFringe%xSeed = large
    oFringe%wallInd = 0

    ! Assume each cell will get just one donor. It's just a guess, it
    ! will be expanded if necessary so the exact value doesn't matter.
    allocate(oFringe%fringeIntBuffer(5, mm), ofringe%fringeRealBuffer(4, mm))
    
    oFringe%nDonor = 0
    ! Now loop over the actual compute cells, setting the cell center
    ! value 'x', the volume and flag these cells as compute
    ii = 0
    do k=2, kl
       do j=2, jl
          do i=2, il
             ii = ii + 1
             do iDim=1, 3
                oFringe%x(iDim, ii) = eighth*(&
                     x(i-1, j-1, k-1, iDim) + &
                     x(i  , j-1, k-1, iDim) + &
                     x(i-1, j  , k-1, iDim) + &
                     x(i  , j  , k-1, iDim) + &
                     x(i-1, j-1, k  , iDim) + &
                     x(i  , j-1, k  , iDim) + &
                     x(i-1, j  , k  , iDim) + &
                     x(i  , j  , k  , iDim))
             end do
             oFringe%xSeed(:, ii) = xSeed(i, j, k, :)
             oFringe%wallInd(ii) = wallInd(i, j, k)
             oFringe%fringeIntBuffer(4, ii) = nn 
             oFringe%fringeIntBuffer(5, ii) = windIndex(i, j, k, il, jl, kl)

          end do
       end do
    end do

    ! We also need to flag a single layer of cells next a wall
    ! boundary condition as being "isWall". This information is
    ! necessary to be able to determine the "wall donors" which are
    ! the flood seeds.

    do mm=1, nBocos
       select case(BCFaceID(mm))
       case (iMin)
          iStart=2; iEnd=2;
          jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (iMax)
          iStart=il; iEnd=il;
          jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (jMin)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=2; jEnd=2;
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (jMax)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=jl; jEnd=jl;
          kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
       case (kMin)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
          kStart=2; kEnd=2;
       case (kMax)
          iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
          jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
          kStart=kl; kEnd=kl;
       end select

       if (isWallType(BCType(mm)) .or. BCType(mm) == SubsonicOutflow & 
            .or. BCType(mm) == SubSonicInflow) then 
          do k=kStart, kEnd
             do j=jStart, jEnd
                do i=iStart, iEnd
                   ! Recompute the index
                   ii = (k-2)*nx*ny + (j-2)*nx + (i-2) + 1
                   oFringe%isWall(ii) = bcFaceID(mm)
                end do
             end do
          end do
       end if
    end do ! BocoLoop

    ! Flag this set of fringes as being allocated
    oFringe%allocated = .True.
  end subroutine initializeOFringes

  subroutine initializeOSurf(famList, oSurf, dualMesh, cluster)

    ! This routine builds the ADT tree for any wall surfaces for the
    ! block currently being pointed to by block Pointers.
    use constants
    use overset, only : oversetWall
    use blockPointers, only : nBocos, BCData, BCFaceID, il, jl, kl, &
         BCFaceID, x, BCType, rightHanded, nbklocal
    use adtBuild, only : buildSerialQuad
    use kdtree2_module, onlY : kdtree2_create
    use oversetPackingRoutines, only : getWallSize
    use sorting, only : bsearchIntegers
    implicit none 

    ! Input Params
    integer(kind=intType), intent(in), dimension(:) :: famList
    type(oversetWall), intent(inout) :: oSurf
    logical, intent(in) :: dualMesh 
    integer(kind=intType), intent(in) :: cluster

    ! Working paramters
    integer(kind=intType) :: i, j, k, n, ii, jj, jjj, mm, ni, nj, nodeCount
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, nNodes, maxCells, nCells, iNode
    logical :: regularOrdering

    ! Set all the sizes for this block.
    oSurf%il = il
    oSurf%jl = jl
    oSurf%kl = kl

    call getWallSize(famList, nNodes, maxCells, dualMesh)

    oSurf%nNodes = nNodes
    oSurf%maxCells = maxCells
    oSurf%cluster = cluster
    ! Allocate space for the x array and connectivity array. cellPtr is
    ! larger than necessary.
    allocate(oSurf%x(3, nNodes), oSurf%conn(4, maxCells), &
         oSurf%cellPtr(maxCells), oSurf%iBlank(maxCells), &
         oSurf%delta(nNodes), oSurf%nte(4, nNodes))

    ii = 0 ! Cumulative node counter
    jj = 0 ! Cumulative cell counter (with iblanks)
    jjj = 0 !Cumulative cell counter (without iblanks)
    nodeCount = 0

    do mm=1, nBocos
       famInclude: if (bsearchIntegers(BCData(mm)%famID, famList) > 0) then 

          select case(BCFaceID(mm))
          case(iMin, jMax, kMin)
             regularOrdering = .True.
          case default
             regularOrdering = .False.
          end select

          ! Now this can be reversed *again* if we have a block that
          ! is left handed. 
          if (.not. rightHanded) then 
             regularOrdering = .not. (regularOrdering)
          end if
          
          ! THIS IS SUPER IMPORTANT: It is absolutely critical that the
          ! wall be built *FROM THE DUAL MESH!!* IT WILL NOT WORK IF YOU
          ! USE THE PRIMAL MESH! The -1 for the node ranges below gives
          ! the extra '1' node for the mesh formed from the dual cells. 

          dualCheck: if (dualMesh) then 
             jBeg = BCData(mm)%jnBeg-1 ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg-1 ; iEnd = BCData(mm)%inEnd
             ! Now fill up the point array
             do j=jBeg, jEnd
                do i=iBeg, iEnd 
                   ii = ii +1
                   select case(BCFaceID(mm))
                   case(imin)
                      oSurf%x(:,ii) = fourth*(x(1, i, j, :) + x(1, i+1, j, :) + &
                           x(1, i, j+1, :) + x(1, i+1, j+1, :))
                   case(imax)
                      oSurf%x(:,ii) = fourth*(x(il, i, j, :) + x(il, i+1, j, :) + &
                           x(il, i, j+1, :) + x(il, i+1, j+1, :))
                   case(jmin) 
                      oSurf%x(:,ii) = fourth*(x(i, 1, j, :) + x(i+1, 1, j, :) + &
                           x(i, 1, j+1, :) + x(i+1, 1, j+1, :))
                   case(jmax) 
                      oSurf%x(:,ii) = fourth*(x(i, jl, j, :) + x(i+1, jl, j, :) + &
                           x(i, jl, j+1, :) + x(i+1, jl, j+1, :))
                   case(kmin) 
                      oSurf%x(:,ii) = fourth*(x(i, j, 1, :) + x(i+1, j, 1, :) + &
                           x(i, j+1, 1, :) + x(i+1, j+1, 1, :))
                   case(kmax) 
                      oSurf%x(:,ii) = fourth*(x(i, j, kl, :) + x(i+1, j, kl, :) + &
                           x(i, j+1, kl, :) + x(i+1, j+1, kl, :))
                   end select
                end do
             end do

             ! Fill up the conn array. Note that don't take the
             ! surface`normal direction (in or out) or the cell
             ! handed-ness into account...it is not necessary since we
             ! are just getting distance to the wall, which is
             ! independent of the orientation.

             ni = iEnd - iBeg + 1
             nj = jEnd - jBeg + 1
             do j=0, nj-2
                do i=0, ni-2
                   jj = jj + 1
                   if (regularOrdering) then 
                      oSurf%conn(1, jj) = nodeCount + (j  )*ni + i + 1 ! n1
                      oSurf%conn(2, jj) = nodeCount + (j  )*ni + i + 2 ! n2
                      oSurf%conn(3, jj) = nodeCount + (j+1)*ni + i + 2 ! n3
                      oSurf%conn(4, jj) = nodeCount + (j+1)*ni + i + 1 ! n4
                   else
                      oSurf%conn(1, jj) = nodeCount + (j  )*ni + i + 1 ! n1
                      oSurf%conn(2, jj) = nodeCount + (j+1)*ni + i + 1 ! n4
                      oSurf%conn(3, jj) = nodeCount + (j+1)*ni + i + 2 ! n3
                      oSurf%conn(4, jj) = nodeCount + (j  )*ni + i + 2 ! n2
                   end if
                end do
             end do
             nodeCount = nodeCount + ni*nj

             ! We don't care about iBlank, cellPtr or delta for the dual
             ! mesh
             oSurf%iBlank = 1
             oSurf%cellPtr = 0
             oSurf%delta = zero
          else ! Using the primal mesh
             jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd

             ! Now fill up the point array. Owned node loop.
             do j=jBeg, jEnd
                do i=iBeg, iEnd 
                   ii = ii +1
                   select case(BCFaceID(mm))
                   case(imin)
                      oSurf%x(:,ii) = x(1, i, j, :)
                   case(imax)
                      oSurf%x(:,ii) = x(il, i, j, :) 
                   case(jmin) 
                      oSurf%x(:,ii) = x(i, 1, j, :) 
                   case(jmax) 
                      oSurf%x(:,ii) = x(i, jl, j, :) 
                   case(kmin) 
                      oSurf%x(:,ii) = x(i, j, 1, :) 
                   case(kmax) 
                      oSurf%x(:,ii) = x(i, j, kl, :) 
                   end select
                   oSurf%delta(ii) = BCData(mm)%deltaNode(i, j)
                end do
             end do

             ! Fill up the conn array being careful to *only* adding
             ! cells that are not already blanked. 
             ni = iEnd - iBeg + 1
             nj = jEnd - jBeg + 1
             do j=0, nj-2
                do i=0, ni-2
                   jjj = jjj + 1
                   oSurf%iBlank(jjj) = BCData(mm)%iblank(iBeg+i+1,jBeg+j+1)
                   if (oSurf%iBlank(jjj) == 1) then 
                      jj = jj + 1
                      if (regularOrdering) then 
                         oSurf%conn(1, jj) = nodeCount + (j  )*ni + i + 1 ! n1
                         oSurf%conn(2, jj) = nodeCount + (j  )*ni + i + 2 ! n2
                         oSurf%conn(3, jj) = nodeCount + (j+1)*ni + i + 2 ! n3
                         oSurf%conn(4, jj) = nodeCount + (j+1)*ni + i + 1 ! n4
                      else
                         oSurf%conn(1, jj) = nodeCount + (j  )*ni + i + 1 ! n1
                         oSurf%conn(2, jj) = nodeCount + (j+1)*ni + i + 1 ! n4
                         oSurf%conn(3, jj) = nodeCount + (j+1)*ni + i + 2 ! n3
                         oSurf%conn(4, jj) = nodeCount + (j  )*ni + i + 2 ! n2
                      end if
                      oSurf%cellPtr(jj) = jjj
                   end if
                end do
             end do
             nodeCount = nodeCount + ni*nj
          end if dualCheck
       end if famInclude
    end do

    ! Set the actual number of cells
    oSurf%nCells = jj

    ! Build the tree itself.
    call buildSerialQuad(oSurf%nCells, nNodes, oSurf%x, oSurf%conn, oSurf%ADT)

    ! Build the KDTree
    if (oSurf%nNodes > 0) then 
       oSurf%tree => kdtree2_create(oSurf%x)
    end if

    ! Build the inverse of the connectivity, the nodeToElem array. 
    oSurf%nte = 0
    do i=1, oSurf%nCells
       do j=1, 4
          n = oSurf%conn(j, i)
          inner:do k=1, 4
             if (oSurf%nte(k, n) == 0) then 
                oSurf%nte(k, n) = i
                exit inner
             end if
          end do inner
       end do
    end do

    ! Flag this wall as being allocated
    oSurf%allocated = .True.

  end subroutine initializeOSurf

end module oversetInitialization
