!
!      ******************************************************************
!      *                                                                *
!      * File:          killFunctions.F90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-12-2003                                      *
!      * Last modified: 10-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine set_signal_write
!
!      ******************************************************************
!      *                                                                *
!      * set_signal_write sets the localSignal to signalWrite. On the   *
!      * finest mesh this means that after the current iteration a      *
!      * solution is written. On the coarser grids this signal will be  *
!      * ignored.                                                       *
!      *                                                                *
!      * This routine is only compiled when signalling is supported.    *
!      *                                                                *
!      ******************************************************************
!
#ifndef USE_NO_SIGNALS

       use communication
       use constants
       use inputPhysics
       use killSignals
       use iteration
       implicit none
!
!      Local variables.
!
       character(len=7) :: integerString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution.                                               *
!      *                                                                *
!      ******************************************************************
!
       ! The user signals must be reconnected again, because the
       ! connection is lost after a signal has been given.

       call connect_signals

       ! Test for a sequential executable.

       testSequential: if( SU_MPI_isSequential ) then

         ! Sequential executable.  Write a general message that a write
         ! signal has been received.

         print "(a)", "#"
         print 100
 100     format("# Received write signal.")

       else testSequential

         ! Parallel executable.  Write a general message that this
         ! processor has received a write signal.

         write(integerString,'(i7)') myID
         integerString = adjustl(integerString)
         integerString = trim(integerString)
         print "(a)", "#"
         print 101, integerString(:len_trim(integerString))
 101     format("# Processor",1X,A,": Received write signal.")

       endif testSequential

       ! Check if a signal was set previously.

       testPrevious: if(localSignal /= noSignal) then

         ! A signal was set before. Set it back to noSignal
         ! and write a message that this has been done.

         localSignal = noSignal

         print 102
 102     format("# Signal was set previously and will now be &
                &overwritten to no signal")

       else testPrevious

         ! No signal was set yet.
         ! Determine the situation and act accordingly.

         if(groundLevel == 1) then

           ! Finest grid. A solution file will be written after either
           ! this multigrid cycle or this time step. This depends whether
           ! this is a steady or an unsteady computation.
           ! In both cases localSignal must be set to signalWrite.

           localSignal = signalWrite

           select case (equationMode)

             case (steady, timeSpectral)
               print 103
 103           format("# Solution will be written after this multigrid &
                      &cycle")

             case (unsteady)
               print 104
 104           format("# Solution will be written after this time step")
           end select

         else

           ! Coarser grid. Signal information will be ignored.

           print 105
           print 106
 105       format("# Solver is still on a coarse grid and therefore &
                  &the signal is ignored.")
 106       format("# Use kill -USR2 if you want to go to the next &
                  &finer grid level.")

         endif

       endif testPrevious

       ! Write a blank line, such that the message is pretty clear.

       print "(a)", "#"

#endif /* USE_NO_SIGNALS */

       end subroutine set_signal_write

!      ==================================================================

       subroutine set_signal_write_quit
!
!      ******************************************************************
!      *                                                                *
!      * set_signal_write_quit sets the localSignal to                  *
!      * signalWriteQuit. On the finest mesh this means that after      *
!      * the current iteration the solution is written and the          *
!      * computation is stopped. On the coarser grids the solution is   *
!      * transferred to the next finer level and the computation is     *
!      * continued there.                                               *
!      *                                                                *
!      * This routine is only compiled when signalling is supported.    *
!      *                                                                *
!      ******************************************************************
!
#ifndef USE_NO_SIGNALS

       use communication
       use constants
       use inputPhysics
       use killSignals
       use iteration
       implicit none
!
!      Local variables.
!
       character(len=7) :: integerString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! The user signals must be reconnected again, because the
       ! connection is lost after a signal has been given.
       ! Only when signalling is wanted.

       call connect_signals

       ! Test for a sequential executable.

       testSequential: if( SU_MPI_isSequential ) then

         ! Sequential executable. Write a general message that a write and
         ! quit signal has been received.

         print "(a)", "#"
         print 200
 200     format("# Received write and quit signal.")

       else testSequential

         ! Parallel executable. Write a general message that this
         ! processor has received a write and quit signal.

         write(integerString,'(i7)') myID
         integerString = adjustl(integerString)
         integerString = trim(integerString)
         print "(a)", "#"
         print 201, integerString(:len_trim(integerString))
 201     format("# Processor",1X,A,": Received write and quit signal.")

       endif testSequential

       ! Check if a signal was set previously.

       testPrevious: if(localSignal /= noSignal) then

         ! A signal was set before. Set it back to noSignal
         ! and write a message that this has been done.

         localSignal = noSignal

         print 202
 202     format("# Signal was set previously and will now be &
                &overwritten to no signal")

       else testPrevious

         ! Set localSignal to signalWriteQuit.

         localSignal = signalWriteQuit

         ! Determine the situation and act accordingly.

         if(groundLevel == 1) then

           ! Finest grid. A solution file will be written and the
           ! computation stopped after either this multigrid cycle or
           ! this time step. This depends whether this is a steady or an
           ! unsteady computation.

           select case (equationMode)
             case (steady, timeSpectral)
               print 203
 203           format("# Solution will be written and computation &
                      &stopped after this multigrid cycle")

             case (unsteady)
               print 204
 204           format("# Solution will be written and computation &
                      &stopped after this time step")
           end select

         else

           ! Coarser grid. The solution is transferred to the next finer
           ! grid level. This happens either after this multigrid cycle
           ! or after this time step depending on whether this is a
           ! steady or an unsteady computation.

           print 205
 205       format("# Solver is still on the coarse grid.")

           select case (equationMode)
             case (steady, timeSpectral)
               print 206
 206           format("# Solution will be transferred to the next &
                      &finer grid after this multigrid cycle.")

             case (unsteady)
               print 207
 207           format("# Solution will be transferred to the next &
                      &finer grid after this time step.")
           end select

         endif

       endif testPrevious

       ! Write a blank line, such that the message is pretty clear.

       print "(a)", "#"

#endif /* USE_NO_SIGNALS */

       end subroutine set_signal_write_quit
