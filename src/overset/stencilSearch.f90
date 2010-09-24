!
!      ******************************************************************
!      *                                                                *
!      * File:          stencilSearch.f90                               *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 04-02-2005                                      *
!      * Last modified: 08-23-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine stencilSearch(level, sps, blk, ind, xp, weights, flag)
!
!      ******************************************************************
!      *                                                                *
!      * stencilSearch finds the interpolation stencil and weights for  *
!      * coordinates xp on the given block, level, and time spectral    *
!      * solution. The algorithm used is a stencil jumping (or gradient *
!      * search) based on the 3 parametric interpolants, which are      *
!      * determined via Newton iterations. The starting indices are     *
!      * given by the vector ind. Flag output possibilities are:        *
!      * - flag >= 0:  search terminated on a 1-to-1 block boundary.    *
!      *               New restart indices, local block, and processor  *
!      *               are given by ind, blk, and flag respectively.    *
!      * - StopAtBoco: search terminiated on block boundary with last   *
!      *               indices given by ind.                            *
!      * - HitMaxIter: maximum number of jumping iterations exceeded    *
!      *               with last indices given by ind.                  *
!      * - Success:   stencil found successfully with indices and 3     *
!      *              parametric interpolants given by ind and s.       *
!      *                                                                *
!      ******************************************************************
!
       use block
       use searchMod
       implicit none
!
!      Local parameters
!
       integer(kind=intType), parameter :: maxJump        = 4
       integer(kind=intType), parameter :: itermaxNewton  = 15
       integer(kind=intType), parameter :: itermaxJumping = 50

       real(kind=realType), parameter :: oscilTol = 0.001_realType
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in)    :: level, sps
       integer(kind=intType), intent(inout) :: blk, ind(3)
       integer(kind=intType), intent(out)   :: flag

       real(kind=realType), intent(in)  :: xp(3)
       real(kind=realType), intent(out) :: weights(nInterp)
!
!      Local variables
!
       integer(kind=intType) :: i, j, k, n, m, c, dind(3)
 
       integer(kind=intType), dimension(3,iterMaxJumping) :: indSave

       real(kind=realType) :: val, s12, s13, s23, s123, ds(3)
       real(kind=realType) :: s(3), a(3,3), f(3), xc(3,8), q(3,8), x1(3)

       logical :: converged
!
!      Function definitions.
!
       logical :: badDonorQuality
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Some initializations.

       s = half
       converged = .false.

       ! Start the search loop.
 
       jumping: do m = 1,itermaxJumping

         ! Store the current set of indices.

         indSave(:,m) = ind
 
         ! Compute the centroids of all 8 cells in the current stencil.
 
         do c = 1,3
           n = 0
           do k = ind(3),ind(3)+1
             do j = ind(2),ind(2)+1
               do i = ind(1),ind(1)+1
                 n = n + 1
                 xc(c,n) = eighth &
                         * sum(flowDoms(blk,level,sps)%x(i-1:i, &
                                                         j-1:j, &
                                                         k-1:k,c))
               end do
             end do
           end do
         end do
 
         ! Calculate the vector q of relative positions.

         q(:,1) =  xc(:,1)
         q(:,2) = -xc(:,1) + xc(:,2)
         q(:,3) = -xc(:,1) + xc(:,3)
         q(:,4) = -xc(:,1) + xc(:,5)
         q(:,5) = -q(:,2)  + xc(:,4) - xc(:,3)
         q(:,6) = -q(:,2)  - xc(:,5) + xc(:,6)
         q(:,7) = -q(:,3)  - xc(:,5) + xc(:,7)
         q(:,8) = -q(:,5)  + xc(:,5) - xc(:,6) + xc(:,8) - xc(:,7)

         ! Compute xp relative to node 1.
 
         x1 = xp - q(:,1)

         ! Start the newton iterations. This will determine the values
         ! of s that will generate the point xp via tri-linear
         ! interpolation from the current stencil.
 
         newton: do n = 1,itermaxNewton
 
           ! Compute the rhs, xp - x(s).

           s12  = s(1)*s(2)
           s13  = s(1)*s(3)
           s23  = s(2)*s(3)
           s123 = s12*s(3)

           f = x1 - (q(:,2)*s(1) + q(:,3)*s(2) + q(:,4)*s(3) &
                  +  q(:,5)*s12  + q(:,6)*s13  + q(:,7)*s23  &
                  +  q(:,8)*s123)

           ! Calculate the jacobian a. This is simply dx(s)/ds.
 
           a(:,1) = q(:,2) + q(:,5)*s(2) + q(:,6)*s(3) + q(:,8)*s23
           a(:,2) = q(:,3) + q(:,5)*s(1) + q(:,7)*s(3) + q(:,8)*s13
           a(:,3) = q(:,4) + q(:,6)*s(1) + q(:,7)*s(2) + q(:,8)*s12
 
           ! Solve the system a*ds = xp - x(s) = f using Cramer's rule.
           ! First compute the determinant. The cutoff is needed to
           ! safeguard against degenerate elements. Then compute ds.
 
           val = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) &
               + a(2,1)*(a(1,3)*a(3,2) - a(1,2)*a(3,3)) &
               + a(3,1)*(a(1,2)*a(2,3) - a(1,3)*a(2,2))

           val = sign(one,val)/max(abs(val),eps)

           ds(1) = val*((a(2,2)*a(3,3) - a(2,3)*a(3,2))*f(1) &
                 +      (a(1,3)*a(3,2) - a(1,2)*a(3,3))*f(2) &
                 +      (a(1,2)*a(2,3) - a(1,3)*a(2,2))*f(3))
           ds(2) = val*((a(2,3)*a(3,1) - a(2,1)*a(3,3))*f(1) &
                 +      (a(1,1)*a(3,3) - a(1,3)*a(3,1))*f(2) &
                 +      (a(1,3)*a(2,1) - a(1,1)*a(2,3))*f(3))
           ds(3) = val*((a(2,1)*a(3,2) - a(2,2)*a(3,1))*f(1) &
                 +      (a(1,2)*a(3,1) - a(1,1)*a(3,2))*f(2) &
                 +      (a(1,1)*a(2,2) - a(1,2)*a(2,1))*f(3))

           ! Update the values of s.

           s = s + ds
 
           ! Check for convergence.
 
           if(sqrt(ds(1)**2 + ds(2)**2 + ds(3)**2) < thresholdReal) then
             converged = .true.
             exit
           end if
 
         end do newton
 
         ! Compute the steps to take in computational space and apply
         ! the jump limit. The limit is applied to s as a real value as
         ! opposed to floor(s) to guard against an element of s being
         ! outside the bounds of default integer.
 
         s    = min(abs(s), real(maxJump, realType))*sign(one, s)
         dind = floor(s)
 
         ! If dind is 0 in all directions, then the 3 interpolants s
         ! are between 0 and 1 and we have found the correct stencil.
         ! If the Newton iterations didn't converge though, cycle to
         ! try again.
 
         if (sum(abs(dind)) == 0) then
           if (converged) then

             ! Well done - set the flag for success, compute the
             ! weights, and check the quality of the donor stencil.

             flag = Success
             call getWeights(s, weights)
             if (badDonorQuality(level, sps, blk, ind, weights)) &
               flag = BadDonor
             return

           else
             cycle
           end if
         end if

         ! Add the delta indices to the current ones.

         ind = ind + dind

         ! Restrict the indices such that they remain within the block
         ! boundaries. Note that the stencil is allowed to straddle a
         ! boundary but not cross it completely (i.e. origin cell index
         ! must satisfy 1 <= i <= il, same for j and k).

         ind(1) = max(1_intType, min(ind(1), flowDoms(blk,level,1)%il))
         ind(2) = max(1_intType, min(ind(2), flowDoms(blk,level,1)%jl))
         ind(3) = max(1_intType, min(ind(3), flowDoms(blk,level,1)%kl))

         ! Check if an oscillation between 2 cells is occurring.

         if (m >= 3) then
           if (all(ind          == indSave(:,m-1)) .and. &
               all(indSave(:,m) == indSave(:,m-2))) then

             ! Make sure the oscillation isn't occurring because the 
             ! Newton iterations did not converge.

             if (.not. converged) cycle

             ! Check if the current stencil is actually withing the 
             ! tolerance, in which case the oscillation is occurring due
             ! to precision issues, so accept the current stencil.

             if (minval(s) >= -oscilTol .and. &
                 maxval(s) <=  oscilTol+one) then

               ! Well done - set the flag for success, compute the
               ! weights, and check the quality of the donor stencil.

               flag = Success
               call getWeights(s, weights)
               if (badDonorQuality(level, sps, blk, ind, weights)) &
                 flag = BadDonor
               return

             end if
           end if
         end if

         ! Calculate the actual change in indices with the restriction.

         i = sum(abs(ind - indSave(:,m)))

         ! If the actual change totals 0 (and from above the jump needed
         ! is not 0), then we are stuck on a block boundary and it's time
         ! to set the output and stop.

         boundaryHalt: if (i == 0) then

           ! Restrict the indices to an owned cell because one or more
           ! indices may be equal to 1.

           indSave(:,m) = max(2_intType, ind)

           ! Set dind to the delta to get from indSave to a first level
           ! direct face halo. This is done only if the desired jump
           ! given by dind wants to move off that face.

           ! I-direction.

           if (ind(1) == flowDoms(blk,level,1)%il .and. dind(1) > 0) then
             dind(1) = 1
           else if (ind(1) == 1 .and. dind(1) < 0) then
             dind(1) = -1
           else
             dind(1) = 0
           end if

           ! J-direction.

           if (ind(2) == flowDoms(blk,level,1)%jl .and. dind(2) > 0) then
             dind(2) = 1
           else if (ind(2) == 1 .and. dind(2) < 0) then
             dind(2) = -1
           else
             dind(2) = 0
           end if

           ! K-direction.

           if (ind(3) == flowDoms(blk,level,1)%kl .and. dind(3) > 0) then
             dind(3) = 1
           else if (ind(3) == 1 .and. dind(3) < 0) then
             dind(3) = -1
           else
             dind(3) = 0
           end if

           ! Loop over the deltas just set. If they are non-zero then the
           ! jump desired was off that face. Multiple options can occur
           ! when we're stuck on a corner or edge of the block. In those
           ! cases, preference is given to jumping across a 1-to-1 face
           ! so that all possible jumping is exploited before declaring
           ! the point an orphan.

           do k = 1,3
             if (dind(k) /= 0) then

               ! Reset to the owned indices, then adjust the kth one so
               ! the indices of a direct halo are obtained.

               ind    = indSave(:,m)
               ind(k) = ind(k) + dind(k)

               ! Do a quick check of the iblank value. If it's 2 then
               ! we know its a physical boundary so skip the next step.

               if (flowDoms(blk,level,sps)%iblank(ind(1),ind(2),ind(3)) &
                   == 2) cycle

               ! Loop over the 1-to-1 subfaces for this block and check
               ! if the cell given by ind is part of it.

               do j=1,flowDoms(blk,level,sps)%n1to1
                 i = flowDoms(blk,level,sps)%nBocos + j

                 if (ind(1) <= flowDoms(blk,level,sps)%icend(i) .and. &
                     ind(1) >= flowDoms(blk,level,sps)%icbeg(i) .and. &
                     ind(2) <= flowDoms(blk,level,sps)%jcend(i) .and. &
                     ind(2) >= flowDoms(blk,level,sps)%jcbeg(i) .and. &
                     ind(3) <= flowDoms(blk,level,sps)%kcend(i) .and. &
                     ind(3) >= flowDoms(blk,level,sps)%kcbeg(i)) then

                   ! Set the block pointers and transform the indices to
                   ! the connecting block.

                   call setPointers(blk, level, sps)
                   call transformCellIndices(i, ind)

                   ! Set the new processor and block.

                   flag = flowDoms(blk,level,1)%neighProc(i)
                   blk  = flowDoms(blk,level,1)%neighBlock(i)
                   return

                 end if
               end do
             end if
           end do

           ! If the routine hasn't returned yet then it is stuck on a
           ! non-1-to-1 boundary, so set the flag and exit.

           flag = StopAtBoco
           return

         end if boundaryHalt

         ! Reset the values of s just to be safe. The Newton iterations
         ! should converge in 3-5 iterations regardless of this. Also
         ! reset the value of converged since a jump has occurred.

         s = half
         converged = .false.

       end do jumping
 
       ! Jumping iterations were exceeded. Set the flag accordingly.
 
       flag = HitMaxIter
 
       end subroutine stencilSearch
