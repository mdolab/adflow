!
!      ******************************************************************
!      *                                                                *
!      * File:          defaultMonitor.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-25-2003                                      *
!      * Last modified: 04-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine defaultMonitor
!
!      ******************************************************************
!      *                                                                *
!      * defaultMonitor sets the default set of variables to be         *
!      * monitored during the convergence. This set depends on the      *
!      * governing equations to be solved.                              *
!      *                                                                *
!      ******************************************************************
!
       use inputPhysics
       use monitor
       use cgnsNames
       implicit none
!
!      Local variables.
!
       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! CPU time is written to stdout.

       showCPU = .true.

       ! Determine the governing equations to be solved.

       select case (equations)
         case (EulerEquations)

           ! Set the number of summation and maximum monitor variables
           ! and allocate the memory for the monitoring names.
           ! A distinction is made between internal and external flows,
           ! because cl and cd do not make a lot of sense for the former.

           if(flowType == internalFlow) then

             ! Internal flow; only the density residual is monitored.

             nMonSum = 1; nMonMax = 0; nMon = 1
             allocate(monNames(nMon), stat=ierr)
             if(ierr /= 0)                       &
               call terminate("defaultMonitor", &
                              "Memory allocation failure for monNames")

             ! Set the names for the variables to be monitored.

             monNames(1) = cgnsL2resRho

           else

             ! External; also lift and drag is monitored.

             nMonSum = 3; nMonMax = 0; nMon = 3
             allocate(monNames(nMon), stat=ierr)
             if(ierr /= 0)                       &
               call terminate("defaultMonitor", &
                              "Memory allocation failure for monNames")

             ! Set the names for the variables to be monitored.

             monNames(1) = cgnsL2resRho
             monNames(2) = cgnsCl
             monNames(3) = cgnsCd

           endif

         case (NSEquations)

           ! Set the number of summation and maximum monitor variables
           ! and allocate the memory for the monitoring names.
           ! A distinction is made between internal and external flows,
           ! because cl and cd do not make a lot of sense for the former.

           if(flowType == internalFlow) then

             ! Internal flow; only the density residual is monitored.

             nMonSum = 1; nMonMax = 0; nMon = 1
             allocate(monNames(nMon), stat=ierr)
             if(ierr /= 0)                       &
               call terminate("defaultMonitor", &
                              "Memory allocation failure for monNames")

             ! Set the names for the variables to be monitored.

             monNames(1) = cgnsL2resRho

           else

             ! External; also lift and drag (total and viscous)
             ! is monitored.

             nMonSum = 4; nMonMax = 0; nMon = 4
             allocate(monNames(nMon), stat=ierr)
             if(ierr /= 0)                       &
               call terminate("defaultMonitor", &
                              "Memory allocation failure for monNames")

             ! Set the names for the variables to be monitored.

             monNames(1) = cgnsL2resRho
             monNames(2) = cgnsCl
             monNames(3) = cgnsCd
             monNames(4) = cgnsCdv

           endif

         case (RANSEquations)

           ! Set the number of summation and maximum monitor variables
           ! and allocate the memory for the monitoring names.
           ! A distinction is made between internal and external flows,
           ! because cl and cd do not make a lot of sense for the former.

           if(flowType == internalFlow) then

             ! Internal flow; the density residual as well as the
             ! maximum values of yplus and the eddy viscosity ration
             ! are monitored.

             nMonSum = 1; nMonMax = 2; nMon = 3
             allocate(monNames(nMon), stat=ierr)
             if(ierr /= 0)                       &
               call terminate("defaultMonitor", &
                              "Memory allocation failure for monNames")

             ! Set the names for the variables to be monitored.

             monNames(1) = cgnsL2resRho
             monNames(2) = cgnsYplusMax
             monNames(3) = cgnsEddyMax

           else

             ! External; also lift and drag (total and viscous)
             ! is monitored.


             nMonSum = 4; nMonMax = 2; nMon = 6
             allocate(monNames(nMon), stat=ierr)
             if(ierr /= 0)                       &
               call terminate("defaultMonitor", &
                              "Memory allocation failure for monNames")

             ! Set the names for the variables to be monitored.

             monNames(1) = cgnsL2resRho
             monNames(2) = cgnsCl
             monNames(3) = cgnsCd
             monNames(4) = cgnsCdv
             monNames(5) = cgnsYplusMax
             monNames(6) = cgnsEddyMax

           endif

       end select

       end subroutine defaultMonitor
