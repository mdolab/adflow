!
!      ******************************************************************
!      *                                                                *
!      * File:          timePeriodSpectral.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-23-2004                                      *
!      * Last modified: 06-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine timePeriodSpectral
!
!      ******************************************************************
!      *                                                                *
!      * timePeriodSpectral determines the time of one period for the   *
!      * time spectral method. It is possible that sections have        *
!      * different periodic times.                                      *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       use inputMotion
       use inputPhysics
       use section
       implicit none
!
!      Local parameter.
!
       real(kind=realType), parameter :: tol  = 1.e-5_realType
!
!      Local variables.
!
       integer :: ierr
       integer(kind=intType) :: nn

       real(kind=realType) :: tt, omega
       real(kind=realType) :: timePeriod

       logical :: timeDetermined
!
!      Function definition.
!
       real(kind=realType) :: commonTimeSpectral
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! This routine is only used for the spectral solutions. Return
       ! immediately if a different mode is solved.

       if(equationMode /= timeSpectral) return

       ! First check if a rotational frequency has been specified.
       ! Only for external flows.

       timeDetermined = .false.

       externalTest: if(flowType == externalFlow) then

         ! X-rotation.

         if(degreeFourXRot > 0) then
          timePeriod     = two*pi/omegaFourXRot
          timeDetermined = .true.
         endif

         ! Y-rotation.

         if(degreeFourYRot > 0) then
           tt = two*pi/omegaFourYRot

           ! Check if a time period was already determined. If so, try
           ! to determine a common time. Otherwise just copy the data.

           if( timeDetermined ) then
             timePeriod = commonTimeSpectral(timePeriod, tt)
           else
             timePeriod     = tt
             timeDetermined = .true.
           endif
         endif

         ! Z-rotation.

         if(degreeFourZRot > 0) then
           tt = two*pi/omegaFourZRot

           ! Check if a time period was already determined. If so, try
           ! to determine a common time. Otherwise just copy the data.

           if( timeDetermined ) then
             timePeriod = commonTimeSpectral(timePeriod, tt)
           else
             timePeriod     = tt
             timeDetermined = .true.
           endif
         endif

       endif externalTest

       ! If it was possible to determine the time, copy it to the
       ! sections and return.

       if( timeDetermined ) then
         do nn=1,nSections
           sections(nn)%timePeriod = timePeriod/sections(nn)%nSlices
         enddo
         return
       endif

       ! Try to determine the periodic time via the rotation rate of the
       ! sections and its number of slices.

       sectionLoop: do nn=1,nSections

         ! Test if the section is rotating, because only for rotating
         ! sections the periodic time can be determined.

         testRotating: if( sections(nn)%rotating ) then

           ! Determine the magnitude of the rotation rate and the
           ! corresponding periodic time period.

           omega = sqrt(sections(nn)%rotRate(1)**2 &
                 +      sections(nn)%rotRate(2)**2 &
                 +      sections(nn)%rotRate(3)**2)

           tt = two*pi/omega

           ! If a time period was already determined, check if this is
           ! identical to tt. If not print an error message and exit.

           if( timeDetermined ) then

             tt = abs(tt-timePeriod)/timePeriod
             if(tt > tol) then
               if(myID == 0)                          &
                 call terminate("timePeriodSpectral", &
                                "Rotational frequencies of the rotating &
                                 &sections are not identical.")
               call mpi_barrier(SUmb_comm_world, ierr)
             endif

           else

             ! Just copy the data.

             timePeriod     = tt
             timeDetermined = .true.

           endif

         endif testRotating
       enddo sectionLoop

       ! Divide the periodic time by the number of slices to get the
       ! characteristic time for every section.

       do nn=1,nSections
         sections(nn)%timePeriod = timePeriod/sections(nn)%nSlices
       enddo

       ! Return if it was possible to determine the time.

       if( timeDetermined ) return

       ! Periodic time could not be determined. Print an error
       ! message and exit.

       if(myID == 0)                          &
         call terminate("timePeriodSpectral", &
                        "Not possible to determine the periodic time &
                        &for the time spectral method")
       call mpi_barrier(SUmb_comm_world, ierr)

       end subroutine timePeriodSpectral

!      ==================================================================

       function commonTimeSpectral(t1, t2)
!
!      ******************************************************************
!      *                                                                *
!      * The function commonTimeSpectral determines the smallest        *
!      * possible common time between t1 and t2, such that              *
!      * tcommon = n1*t1 = n2*t2 and n1, n2 integers.                   *
!      *                                                                *
!      ******************************************************************
!
       use communication
       implicit none
!
!      Function definition
!
       real(kind=realType) :: commonTimeSpectral
!
!      Function arguments.
!
       real(kind=realType), intent(in) :: t1, t2
!
!      Local parameters.
!
       integer(kind=intType), parameter :: nMax = 100
       real(kind=realType),   parameter :: tol  = 1.e-5_realType
!
!      Local variables.
!
       integer                :: ierr
       integer(kind=intType) :: n1, n2
       real(kind=realType)   :: tt1, tt2, tt, ratio
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the largest time in tt1 and the smallest in tt2 and
       ! compute the ratio tt1/tt2, which is >= 1

       tt1   = max(t1, t2)
       tt2   = min(t1, t2)
       ratio = tt1/tt2

       ! Loop to find the smallest integer values of n1 and n2, such
       ! that n1*tt1 = n2*tt2. Note that due to the previous definition
       ! n2 >= n1.

       do n1=1,nMax
         tt = n1*ratio
         n2 = nint(tt)
         tt = abs(tt-n2)
         if(tt <= tol) exit
       enddo

       ! Check if a common time was found

       if(n1 > nMax) then
         if(myID == 0)                          &
           call terminate("commonTimeSpectral", &
                          "No common time periodic time found. &
                          &Problem may not be periodic")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Set the common time.

       commonTimeSpectral = n1*tt1

       end function commonTimeSpectral
