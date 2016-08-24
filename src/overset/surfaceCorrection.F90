subroutine surfaceCorrection(oBlock, oFringe, offset, n)

  use overset
  use adtAPI
  use kdtree2_module
  use inputOverset
  use sorting, only : unique
  implicit none

  ! Input/Output
  type(oversetBlock), intent(inout) :: oBlock
  type(oversetFringe), intent(inout) :: oFringe
  integer(kind=intType), intent(in) :: n
  real(kind=realType), intent(out), dimension(3, n) :: offset

  ! Working 
  integer(kind=intType) :: i, j, k, ii, jj, nInterpol
  integer(kind=intType) :: cellID, idx, nUnique

  integer(kind=intType), dimension(3) :: intInfoF, intInfoB
  integer(kind=intType), dimension(4) :: nodesB, nodesF
  real(kind=realType), dimension(3, 2) :: dummy
  real(kind=realType), dimension(5) :: uvwF, uvwB
  real(kind=realType), dimension(3) :: ptB, ptF, yy
  real(kind=realType), dimension(4) :: weightsF, weightsB, xx
  real(kind=realType) :: ratio, fact, distY, q1(3, 4), q2(3, 4), dB
  type(kdtree2_result) :: results(1)
  logical :: overlapped1, overlapped2, overlapped

  integer(kind=intType), dimension(:), allocatable :: link, tmp
  real(kind=realType), dimension(:, :), allocatable :: uniqueWallPts, masterOffset
  ! Variables we have to pass the ADT search routine
  integer(kind=intType), dimension(:), pointer :: frontLeaves
  integer(kind=intType), dimension(:), pointer :: frontLeavesNew
  type(adtBBoxTargetType), dimension(:), pointer :: BB
  type(oversetWall), pointer :: bWall, fWall

  ! Set pointers to walls (even if they are empty) for the cluster
  ! containing the search block (bWall) and the cluster containing the
  ! fringe pts (fWall)
  bWall => clusterWalls(oBlock%cluster)
  fWall => clusterWalls(oFringe%cluster)

  ! Determine if we can make a quick exit:
  if (bWall%nNodes == 0 .or. fWall%nNodes == 0) then 
     ! oBlock cluster or fringeCluster do not have walls. Cannot have
     !  a surface-surface overlap!
     return
  end if

  
  ! Allocate the (pointer) memory that may be resized as necessary for
  ! the singlePoint search routine. 
  allocate(BB(10), frontLeaves(25), frontLeavesNew(25))

  nInterpol = 0

  ! Determine the surface points from our list of fringes. We have the
  ! the list global indices, which essentially just serves to
  ! determine a compact list of nodes to search for. We need to make a
  ! tmp array since unique overwrites the array.
  allocate(tmp(n), link(n))
  tmp = oFringe%wallInd
  call unique(tmp, n, nUnique, link)
  allocate(uniqueWallPts(3, nUnique), masterOffset(3, nUnique))
  masterOffset = zero

  do i=1, n
     uniqueWallPts(:, link(i)) = oFringe%xSeed(:, i)
  end do

  masterLoop: do ii=1, nUnique

     ! The search point we are dealing with:
     xx(1:3) = uniqueWallPts(:, ii)
     xx(4) = large

     ! Project the point onto the oBlock
     call minDistanceTreeSearchSinglePoint(bWall%ADT, &
          xx, intInfoB, uvwB, dummy, nInterpol, BB, frontLeaves, frontLeavesNew)
     dB = sqrt(uvwB(4))

     if (uvwB(1) > zero .and. uvwB(1) < one .and. &
          uvwB(2) > zero .and. uvwB(2) < one) then 

        ! Extract the 4 nodes for this quad element
        do k=1, 4
           q1(:, k) = bWall%x(:, bWall%conn(k, intInfoB(3)))
        end do

        ! This is a little inefficient...what we want to do is
        ! determine the 4 quads surrounding the point I'm looking
        ! for. Use the KDTree to determine the index of the node in
        ! question, then use the nToElem pointer to get the 4 quads
        ! surrounding my node.  
        call kdtree2_n_nearest(fWall%tree, xx(1:3), 1, results)

        idx = results(1)%idx ! Node index on fWall

        overlapped1 = .False.
        overlapped2 = .False.
        overlapped = .False.
        ! Now loop over (up to 4) of the quads surrounding this node:
        quadLoop: do j=1, 4

           cellID = fWall%nte(j, idx)
           if (cellID > 0) then 

              do k=1, 4
                 q2(:, k) = fWall%x(:, fWall%conn(k, cellID))
              end do

              ! Now see if the two quads overlap in the flat sense
              call quadOverlap(q1, q2, overlapped1)
              overlapped2 = .False. 
              if (dB < nearWallDist) then 
                 overlapped2 = .True. 
              end if

              if (overlapped1 .and. overlapped2) then 
                 overlapped = .True.
                 exit quadLoop
              end if
           end if
        end do quadLoop

        if (overlapped) then 

           nodesB = bWall%conn(:, intInfoB(3))
           call getWeights(uvwB(1:2), weightsB)

           ptB = zero
           do j=1,4
              ptB = ptB + weightsB(j)*bWall%x(:, nodesB(j))
           end do

           ! Now set the offset for the wall. 
           masterOffset(:, ii) =  ptB - xx(1:3)
        end if
     end if
  end do masterLoop

  ! Now that we've determined the number of surface offsets, we can
  ! loop back throught he actual nodes and set the wall offset. 

  do ii=1, n
     ! Attenuate the offset over nearWallDist
     distY = norm2(oFringe%x(:, ii) - uniqueWallPts(:, link(ii)))
     ratio = distY / nearWallDist
     fact = max(one - ratio**3, zero)
     offset(:, ii) = offset(:, ii) + fact*masterOffset(:, link(ii))
  end do

  ! Make sure to clean up the pointer allocations
  deallocate(BB, frontLeaves, frontLeavesNew)
  deallocate(uniqueWallPts, masterOffset, tmp, link)

contains
  subroutine getWeights(uv, weights)
    use constants
    implicit none

    real(kind=realType), intent(in) :: uv(2)
    real(kind=realType), intent(out) :: weights(4)
    weights(1) = (one - uv(1))*(one - uv(2))
    weights(2) = (      uv(1))*(one - uv(2))
    weights(3) = (      uv(1))*(      uv(2))
    weights(4) = (one - uv(1))*(      uv(2))
  end subroutine getWeights
end subroutine surfaceCorrection
