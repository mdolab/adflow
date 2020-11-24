module wallSearches
contains
  subroutine wallSearch(aSurf, bSurf)

    use constants
    use oversetData, only : oversetWall, clusterAreas
    use inputOverset, only : nearWallDist
    use adtLocalSearch, only : minDistanceTreeSearchSinglePoint
    use adtData, only : adtBBoxTargetType, adtLeafType
    use adtUtils, only : stack
    use utils, only : mynorm2
    implicit none

    ! Input/Output
    type(oversetWall), intent(inout) :: aSurf, bSurf

    ! Working Varaibles
    integer(kind=intType) :: i, jj, k, iElem, maxLevels, nNeighbours, nOtherElem, iOther, otherElem
    integer(kind=intType) :: nInterpol, elemID, intInfo(3), factor, jelem, otherElems(4)
    real(kind=realType) :: uvw(5), xx(4), dist, q1(3, 4), q2(3, 4), delta, radius1, radius2

    ! Variables we have to pass the ADT search routine
    integer(kind=intType), dimension(:), pointer :: frontLeaves
    integer(kind=intType), dimension(:), pointer :: frontLeavesNew
    type(adtBBoxTargetType), dimension(:), pointer :: BB
    real(kind=realType),   dimension(3,2) :: dummy
    integer(kind=intType), dimension(:), allocatable :: tmpCellArr
    integer(kind=intType), dimension(:), allocatable :: tmpNodeElem
    integer(kind=intType), dimension(:), allocatable :: mask
    type(adtLeafType), dimension(:), pointer :: ADTree

    logical :: overlapped

    if (aSurf%nCells == 0 .or. bSurf%nCells == 0) then
       ! Either block doesn't have walls, so there is nothing do but just
       ! return.
       return
    end if

    if (clusterAreas(bSurf%cluster) <= clusterAreas(aSurf%cluster)) then
       ! B is smaller so we don't need to do anything
       return
    end if

    nInterpol = 0
    ! Allocate the (pointer) memory that may be resized as necessary for
    ! the singlePoint search routine.
    allocate(BB(10), frontLeaves(25), frontLeavesNew(25), stack(100))

    ! Basically what we are doing it looping all of our bSurf NODES. We
    ! use a special "surface containment search". Essentially all we are
    ! looking for is if a point it inside of of an actual element
    ! BBox. If it isn't inside any BBox then we know it it can't
    ! overlap. This is essentialy fast cull of the majority of panels so
    ! we can later just focus on the ones that may actually overlap.
    ! node as being blanked

    ! Start with a max 10 layers (each with an unreduced 8 cells)
    maxLevels = 1
    allocate(tmpCellArr(3*3), mask(aSurf%nCells), tmpNodeElem(bSurf%nNodes))
    tmpNodeElem(:) = 0
    mask = 0
    ADTree => aSurf%ADT%ADTree

    do i=1, bSurf%nNodes

       xx(1:3) = bSurf%x(:, i)
       xx(4) = large

       ! Just check if it is inside the root bounding box..ie the full
       ! bounding box of the surface. This is would appear conservative,
       ! but isn't good enough. We need to expand by nearWallDist since
       ! it is possible a overlap occurs right at the edge of the
       ! bounding box.
       if(xx(1) >= ADTree(1)%xMin(1) - nearWallDist .and. &
            xx(1) <= ADTree(1)%xMax(4)  + nearWallDist .and. &
            xx(2) >= ADTree(1)%xMin(2) - nearWallDist .and. &
            xx(2) <= ADTree(1)%xMax(5) + nearWallDist .and. &
            xx(3) >= ADTree(1)%xMin(3) - nearWallDist .and. &
            xx(3) <= ADTree(1)%xMax(6) + nearWallDist) then

          ! Now find the closest element on the other mesh for this
          ! node. This is the regular (expensive) closest point search

          call minDistanceTreeSearchSinglePoint(aSurf%ADT, xx, intInfo, uvw, &
               dummy, nInterpol, BB, frontLeaves, frontLeavesNew)
          elemID = intInfo(3)
          tmpNodeElem(i) = elemID
       end if
    end do

    ! Loop over the cells now since this is eventually want we need to blank out:
    cellLoop: do i=1,bSurf%nCells

       ! Extract out the elems found on the other mesh for the 4 nodes
       ! on my element. There could be none, 1 or up to 4 other elements.

       nOtherElem = 0
       do jj=1, 4
          otherElem = tmpNodeElem(bSurf%conn(jj, i))
          if (otherElem /= 0) then
             nOtherElem = nOtherElem + 1
             otherElems(nOtherElem) = otherElem
          end if
       end do

       ! Get the coordinates of my quad
       do jj=1,4
          q1(:, jj) = bSurf%x(:, bSurf%conn(jj, i))
       end do

       do iOther=1, nOtherElem

          elemID = otherElems(iOther)

          ! Get coordinates of the other (found) quad
          do jj=1,4
             q2(:, jj) = aSurf%x(:, aSurf%conn(jj, elemID))
          end do

          ! Do a quick check of the cell itself. If it overlaps,
          ! we're done and don't need to deal with neighbor cells at
          ! all.
          call quadOverlap(q1, q2, overlapped)
          if (overlapped) then
             bSurf%iBlank(bSurf%cellPtr(i)) = -2

             ! No need to do anything else
             cycle CellLoop
          end if

          ! Otherwise, we need to do more work.
          radius1 = getCellRadius(q1)
          radius2 = getCellRadius(q2)

          ! We technically only should only need to add 1 here, but
          ! to be safer, we'll have at least two layers to check.
          factor = int(radius1/radius2) + 2

          if (factor > maxLevels) then
             deallocate(tmpCellArr)
             maxLevels = factor
             allocate(tmpCellArr((1+2*maxLevels)**2))
          end if

          ! This is where it gets interesing: We can determine the
          ! number of recursive radiating layers we need to check
          ! based on the relative size of the the two quads.

          ! Now for the fun part: Recursion!
          nNeighbours = 0
          call getNeighbourCells(aSurf, mask, elemID, factor, tmpCellArr, nNeighbours)

          ! Now just blindly check them until we run out of find an overlapped one:

          do iElem = 1, nNeighbours
             elemID = tmpCellArr(iElem)

             ! Return the mask for this elem back to 0
             mask(elemID) = 0

             ! Get coordinates of the other quad
             do jj=1,4
                q2(:, jj) = aSurf%x(:, aSurf%conn(jj, elemID))
             end do

             ! Do the actual overlap calc for the found cell:
             call quadOverlap(q1, q2, overlapped)

             if (overlapped) then
                bSurf%iBlank(bSurf%cellPtr(i)) = -2

                ! No need to do anything else, but we we do need to
                ! flip all the mask elements back for the next
                ! iteration of cellLoop

                do jElem=iElem+1, nNeighbours
                   mask(tmpCellArr(jElem)) = 0
                end do

                cycle CellLoop

             end if
          end do
       end do
    end do cellLoop

    deallocate(BB, frontLeaves, frontLeavesNew, stack, tmpCellArr, mask)

  end subroutine wallSearch

  recursive subroutine getNeighbourCells(aSurf, mask, baseElemID, layers, elemList, nElemFound)

    ! This routine recursively assembles a list all neighbours within
    ! "layers" of the the baseELemID. The elemList is sorted such that
    ! there are no duplicates:

    use constants
    use oversetData, only : oversetWall
    implicit none

    ! Input/Output
    type(oversetWall), intent(inout) :: aSurf
    integer(kind=intType), intent(inout), dimension(:) :: mask, elemList
    integer(kind=intType), intent(in) :: baseElemID, layers
    integer(kind=intType), intent(inout) :: nElemFound

    ! Working
    integer(kind=intType) :: i, iNode, iCell, curElem

    ! The recusive chain ends when layers == 0
    if (layers == 0) then
       return
    end if

    ! Loop over the nodes of the given quad:
    do i=1, 4
       iNode = aSurf%conn(i, baseElemID)

       ! Loop over the (up to 4) cells surrounding this node use the
       ! node->elem (nte) array
       do iCell=1,4
          curElem = aSurf%nte(iCell, iNode)
          if (curElem /= 0) then
             ! This is a real cell:

             if (mask(curElem) /= baseElemID .and. mask(curElem) == 0) then
                ! we know we don't need to add the baseElemID and if its
                ! already in the mask we don't have to do anything either

                nElemFound = nElemFound + 1
                elemList(nElemFound) = curElem
                mask(curElem) = 1
                ! Now recursively call again, with the baseElement of curElem and 1 fewer levels
                call getNeighbourCells(aSurf, mask, curElem, layers-1, elemList, nElemFound)
             end if
          end if
       end do
    end do
  end subroutine getNeighbourCells

  subroutine quadOverlap(q1, q2, overlapped)
    ! Given two quad in *3D* determine if they overlap using the
    ! separation axis theorem after projecting onto the plane defined by
    ! the cell normal. Check both normals from each quad.

    use constants
    use utils, only : mynorm2, cross_prod

    implicit none

    ! input/output
    real(kind=realType), dimension(3, 4), intent(in) :: q1, q2
    logical , intent(out) :: overlapped

    ! Working
    integer(kind=intType) :: ii, jj
    real(kind=realType), dimension(2, 4) :: qq1, qq2
    real(kind=realType), dimension(3) :: axis1, axis2, n1, n2, normal, v1, v2, c1, c2
    real(kind=realType) :: e1, e2
    ! Check distance between cell centers
    c1 = zero
    c2 = zero
    do ii = 1,4
       c1 = c1 + fourth*q1(:,ii)
       c2 = c2 + fourth*q2(:,ii)
    end do

    ! Get get max distance between center and node:
    e1 = zero
    e2 = zero
    do ii=1,4
       e1 = max(e1, mynorm2(c1 - q1(:, ii)))
       e2 = max(e2, mynorm2(c2 - q2(:, ii)))
    end do

    ! Check if distance between cell center sid beyond the threshold
    if (mynorm2(c1-c2) .ge. (e1 + e2)) then
       overlapped = .False.
       return
    end if

    ! The two quads *may* be overlapped. We have to do it hard way.

    ! Normal of first quad
    v1 = q1(:, 3) - q1(:, 1)
    v2 = q1(:, 4) - q1(:, 2)
    call cross_prod(v1, v2, n1)
    n1 = n1 / mynorm2(n1)

    ! Normal of second quad
    v1 = q2(:, 3) - q2(:, 1)
    v2 = q2(:, 4) - q2(:, 2)
    call cross_prod(v1, v2, n2)
    n2 = n2/mynorm2(n2)

    ! f the normals are not in the same direction, must be a thin
    ! surface.
    if (dot_product(n1, n2) < zero) then
       overlapped = .False.
       return
    end if

    do ii=1, 2
       if (ii == 1) then
          normal = n1
          axis1 = q1(:, 2) - q1(:, 1)
       else
          normal = n2
          axis1 = q2(:, 2) - q2(:, 1)
       end if

       ! Project axis1 onto the plane and normalize
       axis1 = axis1 - dot_product(axis1, normal)*normal
       axis1 = axis1/mynorm2(axis1)

       ! Axis 2 is now the normal cross axis1
       call cross_prod(normal, axis1, axis2)
       axis2 = axis2/mynorm2(axis2)

       do jj=1, 4
          qq1(1, jj) = dot_product(axis1, q1(:, jj))
          qq1(2, jj) = dot_product(axis2, q1(:, jj))

          qq2(1, jj) = dot_product(axis1, q2(:, jj))
          qq2(2, jj) = dot_product(axis2, q2(:, jj))
       end do
       call quadOverlap2D(qq1, qq2, overlapped)

       if (overlapped) then
          return
       end if
    end do

  end subroutine quadOverlap

  subroutine quadOverlap2D(q1, q2, overlapped)
    ! Given two quad in *2D* determine if they overlap using the
    ! separation axis theorem

    use constants
    implicit none

    ! input/output
    real(kind=realType), dimension(2, 4), intent(in) :: q1, q2
    logical , intent(out) :: overlapped

    ! Working
    real(kind=realType), dimension(4) :: tmp1, tmp2
    integer(kind=intType) :: ii, jj, kk, jjp1
    real(kind=realType), dimension(2) :: axis, p0
    real(kind=realType) :: min1, max1, min2, max2
    overlapped = .True.
    tmp1 = zero
    tmp2 = zero
    quadLoop: do ii=1, 2 ! Loop over the two quads
       edgeLoop: do jj=1, 4 ! Loop over the edges of each quad
          jjp1 = mod(jj, 4)+1

          if (ii == 1) then
             axis = q1(:, jjp1) - q1(:, jj)
             p0 = q1(:, jj)
          else
             axis = q2(:, jjp1) - q2(:, jj)
             p0 = q2(:, jj)
          end if

          ! Take the axis normal
          axis = (/axis(2), -axis(1)/)

          ! Take the dot products
          do kk=1,4
             tmp1(kk) = dot_product(axis, q1(:, kk) - p0)
             tmp2(kk) = dot_product(axis, q2(:, kk) - p0)
          end do

          min1 = minval(tmp1)
          max1 = maxval(tmp1)

          min2 = minval(tmp2)
          max2 = maxval(tmp2)

          if (max1 < min2 .or. max2 < min1) then
             overlapped = .False.
             ! We can just jump right out since we know they cannot
             ! overlap.
             exit quadLoop
          end if
       end do edgeLoop
    end do quadLoop
  end subroutine quadOverlap2D

  function getCellRadius(q)

    use constants
    use utils, only : mynorm2
    implicit none
    ! Input
    real(kind=realType), dimension(3, 4) :: q
    real(kind=realType) :: getCellRadius

    ! Working
    real(kind=realType) :: c(3)
    integer(kind=intType) :: ii

    c = zero
    do ii=1, 4
       c = c + fourth*q(:, ii)
    end do


    getCellRadius = zero
    do ii=1, 4
       getCellRadius = max(getCellRadius, mynorm2(c - q(:, ii)))
    end do

  end function getCellRadius

end module wallSearches
