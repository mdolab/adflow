!
!      ******************************************************************
!      *                                                                *
!      * File:          badDonorQuality.f90                             *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 08-11-2005                                      *
!      * Last modified: 08-19-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       logical function badDonorQuality(level, sps, blk, ind, weights)
!
!      ******************************************************************
!      *                                                                *
!      * badDonorQuality determines if a donor stencil on the given     *
!      * level, spectral mode, and local block is acceptable or not.    *
!      * The indices vector gives the i,j,k position of the base cell   *
!      * for the interpolation type, and the weights of the donor       *
!      * members must also be provided. A donor stencil is considered   *
!      * bad if it contains any holes or if one minus the sum of the    *
!      * weights of any fringe members is less that the input allowable *
!      * quality.                                                       *
!      *                                                                *
!      ******************************************************************
!
       use block
       use inputOverset
       use searchMod
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: level, sps, blk, ind(3)

       real(kind=realType), intent(in) :: weights(nInterp)
!
!      Local variables
!
       integer(kind=intType) :: i, j, k, nn

       real(kind=realType) :: quality
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the result and the quality.

       badDonorQuality = .false.
       quality         = one

       ! Check the stencil based on the interpolation type.

       select case (interpolationType)
         case (TriLinear)

           ! A trilinear stecil includes 8 cells with i,j,k being the
           ! minimum indices (i.e. the origin).

           nn = 0

           do k = ind(3),ind(3)+1
             do j = ind(2),ind(2)+1
               do i = ind(1),ind(1)+1

                 nn = nn + 1

                 if (flowDoms(blk,level,sps)%iblank(i,j,k) <=  0) then
                   badDonorQuality = .true.
                   return
                 else if (flowDoms(blk,level,sps)%iblank(i,j,k) >= 9) then
                   badDonorQuality = .true.
                   quality = quality - weights(nn)
                 end if

               end do
             end do
           end do

       end select

       ! Return if no fringe was detected in the stencil or if the 
       ! allowable input quality is being ignored.

       if (ignoreCutoffQuality .or. .not. badDonorQuality) return

       ! Compare the cutoff with the computed quality and flip the
       ! result if needed.

       if (quality >= allowableDonorQuality) badDonorQuality = .false.

       end function badDonorQuality
