!
!      ******************************************************************
!      *                                                                *
!      * File:          initBleedsFromHalos.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-11-2007                                      *
!      * Last modified: 11-30-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initBleedsFromHalos
!
!      ******************************************************************
!      *                                                                *
!      * initBleedsFromHalos sets the prescribed states in the bleed    *
!      * regions to the current data in the halo cells. This is         *
!      * typically to accomplish a consistent restart from a previously *
!      * computed solution.                                             *
!      *                                                                *
!      ******************************************************************
!
       use bleedFlows
       use blockPointers
       use BCTypes
       use flowVarRefState
       use inputTimeSpectral
       implicit none
!
!      Local parameter.
!
       real(kind=realType), parameter :: twothird = two*third
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, sps, i, j

       real(kind=realType), dimension(:,:),   pointer :: pp1, ps
       real(kind=realType), dimension(:,:,:), pointer :: ww1
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if no outflow bleeds are present.

       if(nInflowBleeds == 0 .and. nOutflowBleeds == 0) return

       ! Loop over the number of time instances and domains.

       do sps=1,nTimeIntervalsSpectral
         do nn=1,nDom

           ! Set the pointers to this block on the finest grid level.
           ! This routine is only called when restarting, hence the
           ! finest grid must be used.

           call setPointers(nn, 1_intType, sps)

           ! Loop over the boundary conditions.

           do mm=1,nBocos

             ! Check for an inflow bleed.

             if(BCType(mm) == MassBleedInflow) then

               call terminate("initBleedsFromHalos", &
                              "No inflow bleeds yet")

             endif

             ! Check for an outflow bleed.

             if(BCType(mm) == MassBleedOutflow) then

               ! Determine the block face on which the subface is
               ! located and set the pointer for pp1 accordingly.

               select case (BCFaceID(mm))
                 case (iMin)
                   pp1 => p(1,1:,1:);  ww1 => w(1,1:,1:,1:)
                 case (iMax)
                   pp1 => p(ie,1:,1:); ww1 => w(ie,1:,1:,1:)
                 case (jMin)
                   pp1 => p(1:,1,1:);  ww1 => w(1:,1,1:,1:)
                 case (jMax)
                   pp1 => p(1:,je,1:); ww1 => w(1:,je,1:,1:)
                 case (kMin)
                   pp1 => p(1:,1:,1);  ww1 => w(1:,1:,1,1:)
                 case (kMax)
                   pp1 => p(1:,1:,ke); ww1 => w(1:,1:,ke,1:)
               end select

               ! Loop over the range of the subface and copy the data
               ! for the static pressure. If a k-equation is present
               ! a correction must be performed to obtain the true
               ! static pressure.

               ps => BCData(mm)%ps

               if( kPresent ) then
                 do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                     ps(i,j) = pp1(i,j) &
                             - twothird*ww1(i,j,irho)*ww1(i,j,itu1)
                   enddo
                 enddo
               else
                 do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                     ps(i,j) = pp1(i,j)
                   enddo
                 enddo
               endif

             endif

           enddo
         enddo
       enddo

       end subroutine initBleedsFromHalos
