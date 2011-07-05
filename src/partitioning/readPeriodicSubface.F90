!
!      ******************************************************************
!      *                                                                *
!      * File:          readPeriodicSubface.F90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-07-2005                                      *
!      * Last modified: 11-07-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readPeriodicSubface(cgnsInd, cgnsBase, zone, conn,  &
                                      connectName, periodic,          &
                                      rotationCenter, rotationAngles, &
                                      translation)
!
!      ******************************************************************
!      *                                                                *
!      * readPeriodicSubface reads the possible periodic info for the   *
!      * given general subface connectivity.                            *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use su_cgns
       implicit none
!
!      Subroutine arguments
!
       integer,          intent(in) :: cgnsInd, cgnsBase, zone, conn
       character(len=*), intent(in) :: connectName

       logical, intent(out) :: periodic

       real(kind=realType), dimension(3), intent(out) :: rotationCenter
       real(kind=realType), dimension(3), intent(out) :: rotationAngles
       real(kind=realType), dimension(3), intent(out) :: translation
#ifdef USE_NO_CGNS

       call terminate("readPeriodicSubface", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
!
!      Local variables.
!
       integer :: ierr
       integer :: jj
       integer :: mass, len, time, temp, angle

       real(kind=cgnsPerType), dimension(3) :: rotCenter, rotAngles
       real(kind=cgnsPerType), dimension(3) :: tlation

       real(kind=realType) :: mult, trans
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check if this is a periodic boundary.

       call cg_conn_periodic_read_f(cgnsInd, cgnsBase, zone, conn, &
                                    rotCenter, rotAngles, tlation, ierr)

       testPeriodic: if(ierr == all_ok) then

         ! Subface is a periodic boundary. Check if the unit for
         ! the rotation angles is specified.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", zone, &
                        "ZoneGridConnectivity_t", 1,             &
                        "GridConnectivity_t", conn,              &
                        "GridConnectivityProperty_t", 1,         &
                        "Periodic_t", 1, "DataArray_t", 2, "end")
         if(ierr /= all_ok)                      &
           call terminate("readPeriodicSubface", &
                          "Something wrong when calling cg_goto_f")

         call cg_units_read_f(mass, len, time, temp, angle, ierr)
         if(ierr == error)                   &
           call terminate("readPeriodicSubface", &
                          "Something wrong when calling cg_units_read_f")

         ! Check if the angle dimensions were specified.

         if(ierr == all_ok .and. angle /= Null) then

           ! Determine the conversion factor to radians.

           call siAngle(angle, mult, trans)

         else

           ! Angle units not specified. Assume radians.
           ! Processor 0 writes a warning to stdout.

           if(myID == 0) then

             print "(a)", "#"
             print "(a)", "#                      Warning"
             print 100, trim(cgnsDoms(zone)%zonename), trim(connectName)
             print "(a)", "#"

 100         format("# Zone",1X,A,", General connectivity",1X,A, &
                    ": No unit specified for periodic angles, &
                    &assuming radians.")

           endif

           ! Set mult to one.

           mult = one

         endif

         ! Store the info. Convert the rotation center and the
         ! translation vector to meters.

         periodic = .true.

         rotationCenter = rotCenter*cgnsDoms(zone)%LRef
         rotationAngles = rotAngles*mult
         translation    = tlation*cgnsDoms(zone)%LRef

         ! Make sure that the rotation angles are such that it
         ! corresponds to an integer value of the number of
         ! sections per wheel.

         mult = sqrt(rotationAngles(1)**2 &
              +      rotationAngles(2)**2 &
              +      rotationAngles(3)**2)

         if(mult > eps) then

           ! Nonzero angle specified. Determine the number of
           ! sections for the full wheel, which is an integer.

           jj = nint(two*pi/mult)

           ! Store the correction factor for the angles in
           ! mult and correct the periodic angles accordingly.

           mult = two*pi/(jj*mult)

           rotationAngles = mult*rotationAngles

         endif

       else testPeriodic

         ! Subface is a normal boundary. Set periodic to .false. and
         ! initialize the periodic data to zero to avoid possible
         ! problems due to uninitialized data.

         periodic = .false.

         rotationCenter = zero
         rotationAngles = zero
         translation    = zero

       endif testPeriodic

#endif

       end subroutine readPeriodicSubface
