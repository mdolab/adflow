subroutine wallSearch(aWall, bWall)

  use constants
  use overset
  use inputOverset
  use adtLocalSearch
  use communication
  implicit none

  ! Input/Output
  type(oversetWall), intent(inout) :: aWall, bWall

  ! Working Varaibles
  integer(kind=intType) :: i, j, elem, jj, n, ii, k
  logical :: found
  integer(kind=intType) :: nInterpol, elemID, intInfo(3)
  real(kind=realType) :: uvw(5), xx(4), dist, q1(3, 4), q2(3, 4), delta

  ! Variables we have to pass the ADT search routine
  integer(kind=intType), dimension(:), pointer :: frontLeaves
  integer(kind=intType), dimension(:), pointer :: frontLeavesNew
  type(adtBBoxTargetType), dimension(:), pointer :: BB
  real(kind=realType),   dimension(3,2) :: dummy
  integer(kind=intType), dimension(:), allocatable :: tmpNodeElem, tmpCellElem
  type(adtLeafType), dimension(:), pointer :: ADTree

  logical :: overlapped, debugit
  if (aWall%nNodes == 0 .or. bWall%nNodes == 0) then 
     ! Either block doesn't have walls, so there is nothing do but just
     ! return. 
     return
  end if

  if (clusterAreas(bWall%cluster) <= clusterAreas(aWall%cluster)) then 
     ! B is smaller so we don't need to do anything
     return
  end if

  nInterpol = 0
  ! Allocate the (pointer) memory that may be resized as necessary for
  ! the singlePoint search routine. 
  allocate(BB(10), frontLeaves(25), frontLeavesNew(25), stack(100))

  ! Basically what we are doing it looping all of our bWall NODES. We
  ! use a special "surface containment search". Essentially all we are
  ! looking for is if a point it inside of of an actual element
  ! BBox. If it isn't inside any BBox then we know it it can't
  ! overlap. This is essentialy fast cull of the majority of panels so
  ! we can later just focus on the ones that may actually overlap.
  ! node as being blanked

  allocate(tmpNodeElem(bWall%nNodes), tmpCellElem(bWall%nCells))
  tmpNodeElem(:) = 0
  tmpCellElem(:) = 0

  ADTree => aWall%ADT%ADTree



  do i=1, bWall%nNodes

     xx(1:3) = bWall%x(:, i)
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

        call minDistanceTreeSearchSinglePoint(aWall%ADT, xx, intInfo, uvw, &
             dummy, nInterpol, BB, frontLeaves, frontLeavesNew)

        ! Don't accept the element just yet. Check that that it is
        ! within a factor of our node tolernace. We have to check both
        ! the node on bWall and the nodes on the cell we found becuase
        ! one could be bigger. 
        dist = sqrt(uvw(4))
        elemID = intInfo(3)

       delta = bWall%delta(i)

        do k=1,4
           delta = max(delta, aWall%delta(aWall%conn(k, elemID)))
        end do

        if (dist < max(nearWallDist, 10*delta)) then 
           ! Store the closest element for this node
           tmpNodeElem(i) = elemID
        end if
     end if
  end do

  ! Also check the cell centers
  do i=1,bWall%nCells

     ! Compute the cell center
     xx(1:3) = zero
     do jj=1,4
        xx(1:3) = xx(1:3) + fourth*bWall%x(:, bwall%conn(jj, i))
     end do

     xx(4) = large

     ! Just check if it is in the first bounding box:
     if(xx(1) >= ADTree(1)%xMin(1) .and. &
          xx(1) <= ADTree(1)%xMax(4) .and. &
          xx(2) >= ADTree(1)%xMin(2) .and. &
          xx(2) <= ADTree(1)%xMax(5) .and. &
          xx(3) >= ADTree(1)%xMin(3) .and. &
          xx(3) <= ADTree(1)%xMax(6)) then

        ! Now find the closest element on the other mesh for this
        ! node. This is the regular (expensive) closest point search

        call minDistanceTreeSearchSinglePoint(aWall%ADT, xx, intInfo, uvw, &
             dummy, nInterpol, BB, frontLeaves, frontLeavesNew)

        ! Don't accept the element just yet. Check that that it is
        ! within a factor of our node tolernace. 
        dist = sqrt(uvw(4))
        elemID = intInfo(3)

        delta = zero
        do k=1, 4
           delta = max(delta, bWall%delta(bWall%conn(k, i)))
           delta = max(delta, aWall%delta(aWall%conn(k, elemID)))
        end do

        if (dist < max(nearWallDist, 10*delta)) then 
        ! Store the closest element for this node
           tmpCellElem(i) = elemID
        end if
     end if
  end do

  ! On the next pass loop over the *cells*
  do i=1, bWall%nCells

     ! Get my coordiantes for my (3D) quad
     do jj=1,4
        q1(:, jj) = bWall%x(:, bWall%conn(jj, i))
     end do

     ! We first check overlap with the 4 elements found from
     ! projecting the nodes.
     do j=1, 4
        n = bWall%conn(j, i)
        elem = tmpNodeElem(n)

        if (elem > 0) then 

           ! Get coordinates of the other quad
           do jj=1,4
              q2(:, jj) = aWall%x(:, aWall%conn(jj, elem))
           end do

           call  quadOverlap(q1, q2, overlapped)
           
           if (overlapped) then 
              bWall%iBlank(bWall%cellPtr(i)) = -2
           end if
        end if
     end do

     ! And check the cell center
     elem = tmpCellElem(i)
     if (elem > 0) then 
        ! Get coordinates of the other quad

        do jj=1,4
           q2(:, jj) = aWall%x(:, aWall%conn(jj, elem))
        end do

        call  quadOverlap(q1, q2, overlapped)
        if (overlapped) then 
           bWall%iBlank(bWall%cellPtr(i)) = -2 ! -2 means it was overlapped and got blanked
        end if
     end if
  end do

  ! ============= This is quite inefficient. ==============
  !
  ! Now do the all the searches in reverse. Probably don't need the
  ! cell center check anymore since the other search will take care of
  ! it. 

  deallocate(tmpNodeElem, tmpCellElem)
  allocate(tmpNodeElem(aWall%nNodes), tmpCellElem(aWall%nCells))
  tmpNodeElem(:) = 0
  tmpCellElem(:) = 0

  ADTree => bWall%ADT%ADTree
  do i=1, aWall%nNodes

     xx(1:3) = aWall%x(:, i)
     xx(4) = large

     ! Just check if it is inside the root bounding box..ie the full
     ! bounding box of the surface. This is pretty conservative.
     if(xx(1) >= ADTree(1)%xMin(1) .and. &
          xx(1) <= ADTree(1)%xMax(4) .and. &
          xx(2) >= ADTree(1)%xMin(2) .and. &
          xx(2) <= ADTree(1)%xMax(5) .and. &
          xx(3) >= ADTree(1)%xMin(3) .and. &
          xx(3) <= ADTree(1)%xMax(6)) then

        ! Now find the closest element on the other mesh for this
        ! node. This is the regular (expensive) closest point search

        call minDistanceTreeSearchSinglePoint(bWall%ADT, xx, intInfo, uvw, &
             dummy, nInterpol, BB, frontLeaves, frontLeavesNew)

        ! Don't accept the element just yet. Check that that it is
        ! within a factor of our node tolernace. We have to check both
        ! the node on bWall and the nodes on the cell we found becuase
        ! one could be bigger. 
        dist = sqrt(uvw(4))
        elemID = intInfo(3)

        delta = aWall%delta(i)

        do k=1,4
           delta = max(delta, bWall%delta(bWall%conn(k, elemID)))
        end do

        if (dist < max(nearWallDist, 10*delta)) then 
           ! Store the closest element for this node
           tmpNodeElem(i) = elemID
        end if
     end if
  end do

  ! Also check the cell centers
  do i=1,aWall%nCells

     ! Compute the cell center
     xx(1:3) = zero
     do jj=1,4
        xx(1:3) = xx(1:3) + fourth*aWall%x(:, awall%conn(jj, i))
     end do

     xx(4) = large

     ! Just check if it is in the first bounding box:
     if(xx(1) >= ADTree(1)%xMin(1) .and. &
          xx(1) <= ADTree(1)%xMax(4) .and. &
          xx(2) >= ADTree(1)%xMin(2) .and. &
          xx(2) <= ADTree(1)%xMax(5) .and. &
          xx(3) >= ADTree(1)%xMin(3) .and. &
          xx(3) <= ADTree(1)%xMax(6)) then

        ! Now find the closest element on the other mesh for this
        ! node. This is the regular (expensive) closest point search

        call minDistanceTreeSearchSinglePoint(bWall%ADT, xx, intInfo, uvw, &
             dummy, nInterpol, BB, frontLeaves, frontLeavesNew)

        ! Don't accept the element just yet. Check that that it is
        ! within a factor of our node tolernace. 
        dist = sqrt(uvw(4))
        elemID = intInfo(3)

        delta = zero
        do k=1, 4
           delta = max(delta, aWall%delta(aWall%conn(k, i)))
           delta = max(delta, bWall%delta(bWall%conn(k, elemID)))
        end do

        if (dist < max(nearWallDist, 10*delta)) then 
        ! Store the closest element for this node
           tmpCellElem(i) = elemID
        end if
     end if
  end do

  ! On the next pass loop over the *cells*
  do i=1, aWall%nCells

     ! Get my coordiantes for my (3D) quad
     do jj=1,4
        q1(:, jj) = aWall%x(:, aWall%conn(jj, i))
     end do

     ! We first check overlap with the 4 elements found from
     ! projecting the nodes.
     do j=1, 4
        n = aWall%conn(j, i)
        elem = tmpNodeElem(n)

        if (elem > 0) then 

           ! Get coordinates of the other quad
           do jj=1,4
              q2(:, jj) = bWall%x(:, bWall%conn(jj, elem))
           end do

           call  quadOverlap(q1, q2, overlapped)
           
           if (overlapped) then 
              bWall%iBlank(bWall%cellPtr(elem)) = -2
           end if
        end if
     end do

     ! And check the cell center
     elem = tmpCellElem(i)
     if (elem > 0) then 
        ! Get coordinates of the other quad

        do jj=1,4
           q2(:, jj) = bWall%x(:, bWall%conn(jj, elem))
        end do

        call  quadOverlap(q1, q2, overlapped)
        if (overlapped) then 
           ! bWall is still the one blanked
           bWall%iBlank(bWall%cellPtr(elem)) = -2 ! -2 means it was overlapped and got blanked
        end if
     end if
  end do


  deallocate(BB, frontLeaves, frontLeavesNew, stack)

end subroutine wallSearch

subroutine quadOverlap(q1, q2, overlapped) 
  ! Given two quad in *3D* determine if they overlap using the
  ! separation axis theorem after projecting onto the plane defined by
  ! the cell normal. Check both normals from each quad.

  use constants
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
     e1 = max(e1, norm2(c1 - q1(:, ii)))
     e2 = max(e2, norm2(c2 - q2(:, ii)))
  end do

  ! Check if distance between cell center sid beyond the threshold
  if (norm2(c1-c2) .ge. (e1 + e2)) then
     overlapped = .False.
     return
  end if

  ! The two quads *may* be overlapped. We have to do it hard way. 

  ! Normal of first quad
  v1 = q1(:, 3) - q1(:, 1)
  v2 = q1(:, 4) - q1(:, 2)
  call cross_prod(v1, v2, n1)
  n1 = n1 / norm2(n1)

  ! Normal of second quad
  v1 = q2(:, 3) - q2(:, 1)
  v2 = q2(:, 4) - q2(:, 2)
  call cross_prod(v1, v2, n2)
  n2 = n2/norm2(n2)

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
     axis1 = axis1/norm2(axis1)

     ! Axis 2 is now the normal cross axis1
     call cross_prod(normal, axis1, axis2)
     axis2 = axis2/norm2(axis2)

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

subroutine cross_prod(a,b,c)

  use precision

  ! Inputs
  real(kind=realType), dimension(3), intent(in) :: a,b

  ! Outputs
  real(kind=realType), dimension(3), intent(out) :: c

  c(1) = a(2) * b(3) - a(3) * b(2)
  c(2) = a(3) * b(1) - a(1) * b(3)
  c(3) = a(1) * b(2) - a(2) * b(1)

end subroutine cross_prod

