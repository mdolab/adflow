subroutine surfaceCorrection(oBlock, oFringe, bWall, fWall, offset, n)

  use overset
  use adtAPI
  use BCTypes
  implicit none

  ! Input/Output
  type(oversetBlock), intent(inout) :: oBlock
  type(oversetFringe), intent(inout) :: oFringe
  type(oversetWall), intent(inout) :: bWall, fWall
  integer(kind=intType), intent(in) :: n
  real(kind=realType), intent(out), dimension(3, n) :: offset

  ! Working 
  integer(kind=intType) :: i, j, k, ii, jj, nx, ny, nz, myIndex, nInterpol
  integer(kind=intType) :: iStart, iEnd, iInc
  integer(kind=intType) :: jStart, jEnd, jInc
  integer(kind=intType) :: kStart, kEnd, kInc

  integer(kind=intType), dimension(3) :: intInfoF, intInfoB, intInfoC
  integer(kind=intType), dimension(4) :: nodesB, nodesF, nodesC
  real(kind=realType), dimension(3, 2) :: dummy
  real(kind=realType), dimension(5) :: uvwF, uvwB, uvwC
  real(kind=realType), dimension(3) :: vecF, vecB, vecC, normF, normB, &
       ptB, ptF,ptC, yy, masterOffset,v1, v2, sss, normalB, normalF
  real(kind=realType), dimension(4) :: weightsF, weightsB, weightsC,xx, xxx
  real(kind=realType) :: ratio, dB, dF, dp, fact, distY

  ! Variables we have to pass the ADT search routine
  integer(kind=intType), dimension(:), pointer :: frontLeaves
  integer(kind=intType), dimension(:), pointer :: frontLeavesNew
  type(adtBBoxTargetType), dimension(:), pointer :: BB

  ! Allocate the (pointer) memory that may be resized as necessary for
  ! the singlePoint search routine. 
  allocate(BB(10), frontLeaves(25), frontLeavesNew(25), stack(100))

  ! Basic algorithm is:

  ! For each wall node (isWall > 0) in 'oFringe'
  !
  !    Project onto its own wall 'fWall'
  !
  !    Project onto its donor's wall 'bWall'
  !
  !    if (u,v) in [0,1] for bWall projection AND distances are comparable
  !
  !        Compute offset vector entry. 
  !
  !        Loop over all other nodes from oFringe along ray:
  !
  !            Set offset by attenuating the wall offset. 

  nInterpol = 0
  nx = oFringe%il -1 
  ny = oFringe%jl -1 
  nz = oFringe%kl -1 
  offset = zero
  masterLoop: do ii=1, n

     if (oFringe%isWall(ii) > 0) then 

        ! The search point we are dealing with:
        xx(1:3) = oFringe%x(:, ii)
        xx(4) = large

        ! Project onto the oBlock *first* since we may be able to
        ! short out out of this doesn't find a solution in the 0-1
        ! range.
        call minDistanceTreeSearchSinglePoint(bWall%ADT, xx, intInfoB, uvwB, &
             dummy, nInterpol, BB, frontLeaves, frontLeavesNew)

        if (uvwB(1) > zero .and. uvwB(1) < one .and. &
             uvwB(2) > zero .and. uvwB(2) < one) then 
           
           ! Ok, so we found a contained solution one the other
           ! wall. Now lets check our own. We almost don't need to do
           ! this, since it essentially has to be the fact below our
           ! point. 

           call minDistanceTreeSearchSinglePoint(fWall%ADT, xx, intInfoF, uvwF, &
                dummy, nInterpol, BB, frontLeaves, frontLeavesNew)
                   
           ! We only continue if the distances are *close*. It could
           ! happen that the first constrained solution is very far
           ! away. Essentialy what this second check does is determine
           ! what is "close to a wall"; it gives us a scaling for the problem

           dB = sqrt(uvwB(4))
           dF = sqrt(uvwF(4))
           ratio = dB/dF
           
           ! We area pretty generours with the distance check. This is
           ! essentially anything within about 500 y+. 

           if (ratio > 1/500_realType .and. ratio < 500_realType) then

              ! Now compute the actual surface points of each projection
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
     
              v1 = bWall%x(:, nodesB(3)) - bWall%x(:, nodesB(1))
              v2 = bWall%x(:, nodesB(4)) - bWall%x(:, nodesB(2))
              normalB(1) = (v1(2)*v2(3) - v1(3)*v2(2))
              normalB(2) = (v1(3)*v2(1) - v1(1)*v2(3))
              normalB(3) = (v1(1)*v2(2) - v1(2)*v2(1))
              
              v1 = fWall%x(:, nodesF(3)) - fWall%x(:, nodesF(1))
              v2 = fWall%x(:, nodesF(4)) - fWall%x(:, nodesF(2))
              normalF(1) = (v1(2)*v2(3) - v1(3)*v2(2))
              normalF(2) = (v1(3)*v2(1) - v1(1)*v2(3))
              normalF(3) = (v1(1)*v2(2) - v1(2)*v2(1))
              

              vecB = xx(1:3) - ptB
              vecF = xx(1:3) - ptF

              ! We do one last check before we commit to a correction:
              ! The dot product of the normalized vectors but put
              ! substantially in the same direction.

              normB = vecB / dB
              normF = vecF / dF

              ! Check our dot product
              normalB = normalB / sqrt(normalB(1)**2 + normalB(2)**2 + normalB(3)**2)
              normalF = normalF / sqrt(normalF(1)**2 + normalF(2)**2 + normalF(3)**2)

              dp = normB(1)*normF(1) + normB(2)*normF(2) + normB(3)*normF(3)

              if (abs(dp) > 0.98_realType) then 
                 
                 ! Now set the offset for the wall. 
                 masterOffset = vecF - vecB

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
                 
                 ! Note the that 'Max' faces have the start at the
                 ! wall so iEnd
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
                          ratio = (distY - dF) / dF 
                          fact = max(one-ratio**3, zero)
                          offset(:, jj) = offset(:, jj) + fact*masterOffset
                       end do
                    end do
                 end do
              end if
           end if
        end if
     end if
  end do masterLoop

  ! Make sure to clean up the pointer allocations
  deallocate(BB, frontLeaves, frontLeavesNew, stack)

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
