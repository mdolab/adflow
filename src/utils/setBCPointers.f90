!
!      ******************************************************************
!      *                                                                *
!      * File:          setBcPointers.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-17-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
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
       integer(kind=intType) :: id, ih
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the face id on which the subface is located and set
       ! the pointers accordinly.

       select case (BCFaceID(nn))

         case (iMin)
           id = 2 + offset;   ih = 1 - offset
           ww1 => w(ih,1:,1:,:); ww2 => w(id,1:,1:,:)
           pp1 => p(ih,1:,1:);   pp2 => p(id,1:,1:)

           if( viscous ) then
             rlv1 => rlv(ih,1:,1:); rlv2 => rlv(id,1:,1:)
           endif

           if( eddyModel ) then
             rev1 => rev(ih,1:,1:); rev2 => rev(id,1:,1:)
           endif

         !===============================================================

         case (iMax)
           id = il - offset;  ih = ie + offset
           ww1 => w(ih,1:,1:,:); ww2 => w(id,1:,1:,:)
           pp1 => p(ih,1:,1:);   pp2 => p(id,1:,1:)

           if( viscous ) then
             rlv1 => rlv(ih,1:,1:); rlv2 => rlv(id,1:,1:)
           endif

           if( eddyModel ) then
             rev1 => rev(ih,1:,1:); rev2 => rev(id,1:,1:)
           endif

         !===============================================================

         case (jMin)
           id = 2 + offset;   ih = 1 - offset
           ww1 => w(1:,ih,1:,:); ww2 => w(1:,id,1:,:)
           pp1 => p(1:,ih,1:);   pp2 => p(1:,id,1:)

           if( viscous ) then
             rlv1 => rlv(1:,ih,1:); rlv2 => rlv(1:,id,1:)
           endif

           if( eddyModel ) then
             rev1 => rev(1:,ih,1:); rev2 => rev(1:,id,1:)
           endif

         !===============================================================

         case (jMax)
           id = jl - offset;  ih = je + offset
           ww1 => w(1:,ih,1:,:); ww2 => w(1:,id,1:,:)
           pp1 => p(1:,ih,1:);   pp2 => p(1:,id,1:)

           if( viscous ) then
             rlv1 => rlv(1:,ih,1:); rlv2 => rlv(1:,id,1:)
           endif

           if( eddyModel ) then
             rev1 => rev(1:,ih,1:); rev2 => rev(1:,id,1:)
           endif

         !===============================================================

         case (kMin)
           id = 2 + offset;   ih = 1 - offset
           ww1 => w(1:,1:,ih,:); ww2 => w(1:,1:,id,:)
           pp1 => p(1:,1:,ih);   pp2 => p(1:,1:,id)

           if( viscous ) then
             rlv1 => rlv(1:,1:,ih); rlv2 => rlv(1:,1:,id)
           endif

           if( eddyModel ) then
             rev1 => rev(1:,1:,ih); rev2 => rev(1:,1:,id)
           endif

         !===============================================================

         case (kMax)
           id = kl - offset;  ih = ke + offset
           ww1 => w(1:,1:,ih,:); ww2 => w(1:,1:,id,:)
           pp1 => p(1:,1:,ih);   pp2 => p(1:,1:,id)

           if( viscous ) then
             rlv1 => rlv(1:,1:,ih); rlv2 => rlv(1:,1:,id)
           endif

           if( eddyModel ) then
             rev1 => rev(1:,1:,ih); rev2 => rev(1:,1:,id)
           endif

       end select

       end subroutine setBCPointers
