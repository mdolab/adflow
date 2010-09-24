!
!      ******************************************************************
!      *                                                                *
!      * File:          mixingPlaneComm.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-27-2005                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mixingPlaneComm(level, color)
!
!      ******************************************************************
!      *                                                                *
!      * mixingPlaneComm determines the communication/interpolation     *
!      * pattern for the mixing plane interface corresponding to the    *
!      * given color on the given grid level.                           *
!      *                                                                *
!      ******************************************************************
!
       use commMixing
       use interfaceGroups
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, color
!
!      Local variables.
!
       integer(kind=intType) :: slideID
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check if this processor actually participates in a search.
       ! If not, return immediately.

       if(.not. myInterfaces(color)%procContributes) return

       ! Store the slide id a bit easier.

       slideID = myInterfaces(color)%globalSlideID

       ! Determine the interpolation pattern for both parts of the
       ! interface. The distinction between part 1 and part 2 is the
       ! sign of the group number. For part 1 the group number is
       ! minus the slide id and for part 2 the group number has a
       ! positive sign. This means that both parts can be built with
       ! the identical routine, except that the slide id's of the
       ! donor and halo must be reversed.

       ! Part 1.

       call mixingPlaneInterpol(level, color, -slideID, slideID, &
                                myInterfaces(color)%nSlices1,    &
                                commPatternMixing(level,color,1))

       ! Part 2.

       call mixingPlaneInterpol(level, color, slideID, -slideID, &
                                myInterfaces(color)%nSlices2,    &
                                commPatternMixing(level,color,2))

       end subroutine mixingPlaneComm
