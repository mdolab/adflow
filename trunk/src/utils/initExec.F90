!
!      ******************************************************************
!      *                                                                *
!      * File:          initExec.F90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-10-2002                                      *
!      * Last modified: 10-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initExec(execName, lenExec, paramName, &
                           lenParam, nArgs, sizeCInt)
!
!      ******************************************************************
!      *                                                                *
!      * InitExec checks the number of command line arguments,          *
!      * set the name of the parameter file and determines the          *
!      * communicator for the given color.                              *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use inputIO
       use iteration
       implicit none
!
!      Subroutine arguments
!
       integer, intent(in) :: lenExec, lenParam, nArgs
       integer, intent(in) :: sizeCInt

       character(len=lenExec),  intent(in) :: execName
       character(len=lenParam), intent(in) :: paramName
!
!      Local variables
!
       integer :: ierr

       character(len=maxStringLen) :: errorMessage
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set standAloneMode to .true. to indicate that the stand alone
       ! mode of the executable is built and deformingGrid to .false.
       ! to indicate that a rigid grid problem is solved here.
       ! Furthermore, set the communicator to mpi_comm_world.
       ! These initializations can be hard coded like this, because this
       ! routine is only called for the stand alone mode. When the code
       ! is run in multi-disciplinary mode via python, the python
       ! equivalent of this routine is used where it is possible to
       ! specify a different communicator than mpi_comm_world.

       standAloneMode  = .true.
       deforming_Grid  = .false.
       SUmb_comm_world = mpi_comm_world

       ! Determine the rank and number of processors inside the group
       ! defined by SUmb_comm_world.

       call mpi_comm_rank(SUmb_comm_world, myID,  ierr)
       call mpi_comm_size(SUmb_comm_world, nProc, ierr)

       ! Check if the program is called correctly.

       if(nArgs == 1) then

         ! No parameter file specified. Processor 0 prints an error
         ! message while the others wait to get killed.

         if(myID == 0) then

           ! Write an error message. This depends whether a parallel
           ! or a sequential code is created.

           if( SU_MPI_isSequential ) then
             write(errorMessage,*) "Usage: ", execName(:lenExec), &
                                   " <parameter file>"
           else
             write(errorMessage,*) "Usage: mpirun -np #procs ", &
                                   execName(:lenExec),          &
                                   "<parameter file>"
           endif

           call terminate("initExec", errorMessage)

         endif

         call mpi_barrier(SUmb_comm_world, ierr)

       endif

       ! Check if a valid executable has been built, i.e. some form
       ! of IO must be possible in parallel mode. Either CGNS or
       ! parallel IO must have been enabled at compile time.

#ifdef USE_NO_CGNS
       if((.not. SU_MPI_isSequential) .and. SU_MPI_noMPIO) then
         if(myID == 0) &
           call terminate("initExec", &
                          "Both CGNS and parallel IO were disabled at &
                          &compile time. Not possible to perform IO")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif
#endif

       ! Set the name of the parameter file.

       paramFile = paramName(:lenParam)
       paramFile = adjustl(paramFile)
       paramFile = trim(paramFile)

       ! Write a message to stdout with information how the
       ! executable was built.

       call writeIntroMessage

       ! Check if the size of the C and fortran integers match.
       ! If not, processor 0 writes an error message and the
       ! program terminates.

       if(sizeCInt /= sizeOfInteger) then

         write(errorMessage,101) sizeCInt, sizeOfInteger
 101     format("Inconsistent size of C and Fortran integers: ", &
                i2," vs ", i2," bytes.")

         if(myID == 0) call terminate("initExec", errorMessage)
         call mpi_barrier(SUmb_comm_world, ierr)

       endif

       ! Allocate the memory for sendRequests and recvRequests.

       allocate(sendRequests(nProc), &
                recvRequests(nProc), stat=ierr)
       if(ierr /= 0)                &
         call terminate("initExec", &
                        "Memory allocation failure for sendRequests &
                        &and recvRequests")

       end subroutine initExec

!========================================================================

       subroutine writeIntroMessage
!
!      ******************************************************************
!      *                                                                *
!      * writeIntroMessage writes a message to stdout with              *
!      * information how the executable was built, e.g. whether single  *
!      * or double precision is used for the integers and reals, etc.   *
!      * To avoid a messy output only processor 0 prints this info.     *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       implicit none
!
!      Local variables
!
       character(len=7) :: integerString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return if this is not processor 0.

       if(myID > 0) return

       ! I'm processor 0. Write the info to stdout.

       print "(a)", "#"
       print "(a)", "# SUmb, multiblock structured flow solver"
       print "(a)", "#"
       print "(a)", "# This code solves the 3D RANS, laminar NS or &
                    &Euler equations"
       print "(a)", "# on multiblock structured hexahedral grids."

       if( SU_MPI_isSequential ) then
         print "(a)", "# This is a sequential executable compiled &
                      &with the following options:"
       else
         write(integerString,"(i7)") nProc
         integerString = adjustl(integerString)
         print "(3a)", "# This is a parallel executable running on ", &
                       trim(integerString), " processors."
         print "(a)", "# It has been compiled with the &
                      &following options:"
       endif

       if( debug ) then
         print "(a)", "# - Debug mode."
       else
         print "(a)", "# - Optimized mode."
       endif

#ifdef USE_LONG_INT
       print "(a)", "# - Size of standard integers: 8 bytes."
#else
       print "(a)", "# - Size of standard integers: 4 bytes."
#endif

#ifdef USE_SINGLE_PRECISION
       print "(a)", "# - Size of standard floating point types: &
                    &4 bytes."

#elif  USE_QUADRUPLE_PRECISION
       print "(a)", "# - Size of standard floating point types: &
                    &16 bytes."
#else
       print "(a)", "# - Size of standard floating point types: &
                    &8 bytes."
#endif

#ifdef USE_PV3
       print "(a)", "# - With pV3 support"
#else
       print "(a)", "# - Without pV3 support"
#endif

#ifdef USE_NO_CGNS
       print "(a)", "# - Without cgns support"
#else
       print "(a)", "# - With cgns support"
#endif

       if(.not. SU_MPI_isSequential) then
         if( SU_MPI_noMPIO ) then
           print "(a)", "# - Without parallel IO support"
         else
           print "(a)", "# - With parallel IO support"
         endif
       endif

#ifdef USE_NO_SIGNALS
       print "(a)", "# - Without support for signals."
#else
       print "(a)", "# - With support for signals."
#endif

       print "(a)", "#"

       end subroutine writeIntroMessage
