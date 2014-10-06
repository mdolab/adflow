!
!      ******************************************************************
!      *                                                                *
!      * File:          resetBCPointers.f90                             *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 10-02-2014                                      *
!      * Last modified: 10-02-2014                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine resetBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                                rev1, rev2, offset)
!
!      ******************************************************************
!      *                                                                *
!      * setBCPointers sets the pointers needed for the boundary        *
!      * condition treatment on a general face, such that the boundary  *
!      * routines are only implemented once instead of 6 times.         *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn, offset

       real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
       real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
       real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
       real(kind=realType), dimension(:,:),   pointer :: rev1, rev2
!
!      Local variables
!
       integer(kind=intType) :: id, ih, ierr, i, j, k
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! nullify all pointers
       nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)

       end subroutine resetBCPointers
