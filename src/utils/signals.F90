  subroutine set_signal_write
      !
      !       set_signal_write sets the localSignal to signalWrite. On the
      !       finest mesh this means that after the current iteration a
      !       solution is written. On the coarser grids this signal will be
      !       ignored.
      !       This routine is only compiled when signalling is supported.
      !
#ifndef USE_NO_SIGNALS

      use constants
      use communication, only: myID
      use inputPhysics, only: equationMode
      use killSignals, only: localSignal, noSignal, signalWrite
      use iteration, only: groundLevel
      use commonFormats, only: strings
      implicit none
      !
      !      Local variables.
      !
      character(len=7) :: integerString

      ! The user signals must be reconnected again, because the
      ! connection is lost after a signal has been given.

      call connect_signals

      ! Parallel executable.  Write a general message that this
      ! processor has received a write signal.

      write (integerString, '(i7)') myID
      integerString = adjustl(integerString)
      integerString = trim(integerString)
      print "(a)", "#"
      print strings, "# Processor ", integerString(:len_trim(integerString)), ": Received write signal."

      ! Check if a signal was set previously.

      testPrevious: if (localSignal /= noSignal) then

          ! A signal was set before. Set it back to noSignal
          ! and write a message that this has been done.

          localSignal = noSignal

          print "(a)", "# Signal was set previously and will now be overwritten to no signal"

          else testPrevious

          ! No signal was set yet.
          ! Determine the situation and act accordingly.

          if (groundLevel == 1) then

              ! Finest grid. A solution file will be written after either
              ! this multigrid cycle or this time step. This depends whether
              ! this is a steady or an unsteady computation.
              ! In both cases localSignal must be set to signalWrite.

              localSignal = signalWrite

              select case (equationMode)

              case (steady, timeSpectral)
                  print "(a)", "# Solution will be written after this iteration"

              case (unsteady)
                  print "(a)", "# Solution will be written after this time step"
              end select

          else

              ! Coarser grid. Signal information will be ignored.

              print "(a)", "# Solver is still on a coarse grid and therefore the signal is ignored."
              print "(a)", "# Use kill -USR2 if you want to go to the next finer grid level."

          end if

      end if testPrevious

      ! Write a blank line, such that the message is pretty clear.

      print "(a)", "#"

#endif /* USE_NO_SIGNALS */

  end subroutine set_signal_write

  !      ==================================================================

  subroutine set_signal_write_quit
      !
      !       set_signal_write_quit sets the localSignal to
      !       signalWriteQuit. On the finest mesh this means that after
      !       the current iteration the solution is written and the
      !       computation is stopped. On the coarser grids the solution is
      !       transferred to the next finer level and the computation is
      !       continued there.
      !       This routine is only compiled when signalling is supported.
      !
#ifndef USE_NO_SIGNALS

      use constants
      use communication, only: myID
      use inputPhysics, only: equationMode
      use killSignals, only: localSignal, noSignal, signalWriteQuit
      use iteration, only: groundLevel
      use commonFormats, only: strings
      implicit none
      !
      !      Local variables.
      !
      character(len=7) :: integerString

      ! The user signals must be reconnected again, because the
      ! connection is lost after a signal has been given.
      ! Only when signalling is wanted.

      call connect_signals

      ! Parallel executable. Write a general message that this
      ! processor has received a write and quit signal.

      write (integerString, '(i7)') myID
      integerString = adjustl(integerString)
      integerString = trim(integerString)
      print "(a)", "#"
      print strings, "# Processor ", integerString(:len_trim(integerString)), ": Received write and quit signal."

      ! Check if a signal was set previously.

      testPrevious: if (localSignal /= noSignal) then

          ! A signal was set before. Set it back to noSignal
          ! and write a message that this has been done.

          localSignal = noSignal

          print "(a)", "# Signal was set previously and will now be overwritten to no signal"

          else testPrevious

          ! Set localSignal to signalWriteQuit.

          localSignal = signalWriteQuit

          ! Determine the situation and act accordingly.

          if (groundLevel == 1) then

              ! Finest grid. A solution file will be written and the
              ! computation stopped after either this multigrid cycle or
              ! this time step. This depends whether this is a steady or an
              ! unsteady computation.

              select case (equationMode)
              case (steady, timeSpectral)
                  print "(a)", "# Solution will be written and computation stopped after this multigrid cycle"

              case (unsteady)
                  print "(a)", "# Solution will be written and computation stopped after this time step"
              end select

          else

              ! Coarser grid. The solution is transferred to the next finer
              ! grid level. This happens either after this multigrid cycle
              ! or after this time step depending on whether this is a
              ! steady or an unsteady computation.

              print "(a)", "# Solver is still on the coarse grid."

              select case (equationMode)
              case (steady, timeSpectral)
                  print "(a)", "# Solution will be transferred to the next finer grid after this multigrid cycle."

              case (unsteady)
                  print "(a)", "# Solution will be transferred to the next finer grid after this time step."
              end select

          end if

      end if testPrevious

      ! Write a blank line, such that the message is pretty clear.

      print "(a)", "#"

#endif /* USE_NO_SIGNALS */

  end subroutine set_signal_write_quit
