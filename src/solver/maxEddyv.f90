!
!      ******************************************************************
!      *                                                                *
!      * File:          maxEddyv.f90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-25-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine maxEddyv(eddyvisMax)
!
!      ******************************************************************
!      *                                                                *
!      * maxEddyv determines the maximum value of the eddy viscosity    *
!      * ratio of the block given by the pointers in blockPointes.      *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), intent(out) :: eddyvisMax
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k

       real(kind=realType) :: eddyvis
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the maximum value to zero and return immediately if
       ! not an eddy viscosity model is used.

       eddyvisMax = zero
       if(.not. eddyModel) return

       ! Loop over the owned cells of this block.

       do k=2,kl
         do j=2,jl
           do i=2,il

             ! Compute the local viscosity ratio and take the maximum
             ! with the currently stored value.

             eddyvis = rev(i,j,k)/rlv(i,j,k)
             eddyvisMax = max(eddyvisMax, eddyvis)

           enddo
         enddo
       enddo

       end subroutine maxEddyv
