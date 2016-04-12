subroutine surfaceCorrection(oBlock, oFringe, bWall, fWall, offset, n)

  use overset
  use adtAPI
  use BCTypes
  use kdtree2_module
  use inputOverset
  implicit none

  ! Input/Output
  type(oversetBlock), intent(inout) :: oBlock
  type(oversetFringe), intent(inout) :: oFringe
  type(oversetWall), intent(inout) :: bWall, fWall
  integer(kind=intType), intent(in) :: n
  real(kind=realType), intent(out), dimension(3, n) :: offset

  ! Working 
  integer(kind=intType) :: i, j, k, ii, jj, nx, ny, nz, myIndex, nInterpol
  integer(kind=intType) :: cellID, idx
  integer(kind=intType) :: iStart, iEnd, iInc
  integer(kind=intType) :: jStart, jEnd, jInc
  integer(kind=intType) :: kStart, kEnd, kInc

  integer(kind=intType), dimension(3) :: intInfoF, intInfoB
  integer(kind=intType), dimension(4) :: nodesB, nodesF
  real(kind=realType), dimension(3, 2) :: dummy
  real(kind=realType), dimension(5) :: uvwF, uvwB
  real(kind=realType), dimension(3) :: ptB, ptF, yy, masterOffset
  real(kind=realType), dimension(4) :: weightsF, weightsB, xx
  real(kind=realType) :: ratio, dB, dF, fact, distY, q1(3, 4), q2(3, 4)
  type(kdtree2_result) :: results(1)
  logical :: overlapped1, overlapped2, overlapped

  ! Variables we have to pass the ADT search routine
  integer(kind=intType), dimension(:), pointer :: frontLeaves
  integer(kind=intType), dimension(:), pointer :: frontLeavesNew
  type(adtBBoxTargetType), dimension(:), pointer :: BB

  ! Allocate the (pointer) memory that may be resized as necessary for
  ! the singlePoint search routine. 
  allocate(BB(10), frontLeaves(25), frontLeavesNew(25))

  nInterpol = 0
  nx = oFringe%il -1 
  ny = oFringe%jl -1 
  nz = oFringe%kl -1 

  masterLoop: do ii=1, n

     if (oFringe%isWall(ii) > 0) then 

        ! The search point we are dealing with:
        xx(1:3) = oFringe%x(:, ii)
        xx(4) = large

      
        ! Project the point onto the oBlock
        call minDistanceTreeSearchSinglePoint(bWall%ADT, xx, intInfoB, uvwB, &
             dummy, nInterpol, BB, frontLeaves, frontLeavesNew)
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
               
                 if (overlapped1 .and.  overlapped2) then 
                    overlapped = .True.
                    exit quadLoop
                 end if
              end if
           end do quadLoop
            
           if (overlapped) then 
              ! Remember to re-initialize the distance. This determines
              ! the distance to my wall. Essentially half the off-wall distance
              xx(4) = large
              call minDistanceTreeSearchSinglePoint(fWall%ADT, xx, intInfoF, uvwF, &
                   dummy, nInterpol, BB, frontLeaves, frontLeavesNew)
              
              ! This is now close the point is to my own wall. Essentialy
              ! the is just the offwall spacing. 
              dF = sqrt(uvwF(4))
              
              ! Now compute the locations on the quad of each
              ! projection
              nodesB = bWall%conn(:, intInfoB(3))
              nodesF = fWall%conn(:, intInfoF(3))
              call getWeights(uvwB(1:2), weightsB)
              call getWeights(uvwF(1:2), weightsF)

              ptB = zero
              ptF = zero
              do j=1,4
                 ptB = ptB + weightsB(j)*bWall%x(:, nodesB(j))
                 ptF = ptF + weightsF(j)*fWall%x(:, nodesF(j))
              end do

              ! Now set the offset for the wall. 
              masterOffset =  ptB - ptF
          
              ! Last thing we need to do is add an attenuating
              ! offset for the nodes in the off-wall direction. 

              ! Back out the i,j,k index of this node from myIndex
              myIndex = oFringe%myIndex(ii)-1
              i = mod(myIndex, nx) + 2
              j = mod(myIndex/nx, ny) + 2
              k = myIndex/(nx*ny) + 2

              iStart = i; iEnd = i
              jStart = j; jEnd = j
              kStart = k; kEnd = k

              select case(oFringe%isWall(ii))
              case(iMin, iMax)
                 iStart=2; iEnd=nx+1
              case(jMin, jMax)
                 jStart=2; jEnd=ny+1
              case(kMin, kMax)
                 kStart=2; kEnd=nz+1
              end select
              do k=kStart, kEnd
                 do j=jStart, jEnd
                    do i=iStart, iEnd
                       ! Recompute the index
                       jj = (k-2)*nx*ny + (j-2)*nx + (i-2) + 1

                       ! Extract the curent point
                       yy = oFringe%x(:, jj)

                       ! Get this distance to the "wall". We just
                       ! use the parametric position we found on
                       ! our fringe.

                       distY = sqrt((yy(1)-ptF(1))**2 + (yy(2)-ptF(2))**2 + (yy(3)-ptF(3))**2)

                       ! Now we can finally compute the normalize ratio
                       ratio = (distY - dF) / dF / 500
                       fact = max(one-ratio**3, zero)
                       offset(:, jj) = offset(:, jj) + fact*masterOffset
                    end do
                 end do
              end do
           end if
        end if
     end if
  end do masterLoop

  ! Make sure to clean up the pointer allocations
  deallocate(BB, frontLeaves, frontLeavesNew)

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
