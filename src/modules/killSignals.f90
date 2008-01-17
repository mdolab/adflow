!
!      ******************************************************************
!      *                                                                *
!      * File:          killSignals.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-12-2003                                      *
!      * Last modified: 08-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module killSignals
!
!      ******************************************************************
!      *                                                                *
!      * This module contains the variables used to handle the          *
!      * kill signals from the user. The user can send two signals that *
!      * can be processed by the code, namely kill -USR1 and kill -USR2.*
!      * The former signal will cause the program to dump a solution    *
!      * File after the current iteration, the latter will dump a       *
!      * solution file and the computation will be stopped. The         *
!      * definition of iteration is different for steady and unsteady.  *
!      * For steady it means after the current multigrid iteration, for *
!      * unsteady after the current time step.                          *
!      *                                                                *
!      * The handling can be switched off at compile time using the     *
!      * compiler flag -DUSE_NO_SIGNALS.                                *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

       ! Definition of some constants for the signalling. As in the
       ! reduce functions a maximum is determined, it is important that
       ! noSignal is less than signalWrite and signalWriteQuit.

       integer(kind=intType), parameter :: noSignal        = 0
       integer(kind=intType), parameter :: signalWrite     = 1
       integer(kind=intType), parameter :: signalWriteQuit = 2

       ! localSignal:  Signal stored on this processor.
       ! globalSignal: Maximum of the local signals.

       integer(kind=intType) :: localSignal, globalSignal

       end module killSignals
