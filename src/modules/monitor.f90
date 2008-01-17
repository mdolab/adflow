!
!      ******************************************************************
!      *                                                                *
!      * File:          monitor.f90                                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-25-2003                                      *
!      * Last modified: 10-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       module monitor
!
!      ******************************************************************
!      *                                                                *
!      * This module contains the variables to be monitored during the  *
!      * convergence as well as the arrays to store the convergence.    *
!      * The latter are only allocated by processor 0.                  *
!      * The default variables to be monitored depend on the governing  *
!      * equations to be solved.                                        *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
       save
!
!      Parameters for the format to write the convergence to stdout.
!
       integer, parameter :: fieldWidth   = 12
       integer, parameter :: decimalWidth =  5
!
!      ******************************************************************
!      *                                                                *
!      * Variables to compute the convergence info.                     *
!      *                                                                *
!      ******************************************************************
!
       ! nMonSum: Number of monitoring variables for which the sum over
       !          all processors must be taken. Note that this is an
       !          integer, because of MPI.
       ! nMonMax: Number of monitoring variables for which the maximum
       !          over all processors must be taken. Note that this is
       !          an integer, because of MPI.
       ! nMon:    The sum of nmonSum and nmonMax

       integer :: nMonSum, nMonMax, nMon

       ! monLoc(nMon):  Array for the local summation/maximum of the
       !                monitoring variables.
       ! monGlob(nMon): Idem, but for the global summation/maximum.
       ! monRef(nMon):  Idem, but for the reference values needed
       !                for an unsteady computation.

       real(kind=realType), dimension(:), allocatable :: monLoc
       real(kind=realType), dimension(:), allocatable :: monGlob
       real(kind=realType), dimension(:), allocatable :: monRef

       ! monNames(nMon): The names of the variables to be monitored.

       character(len=maxCGNSNameLen), dimension(:), allocatable :: &
                                                            monNames

       ! monMachOrHMax: Whether or not the maximum value of the Mach
       !                number and total enthalpy difference must be
       !                monitored and thus computed.
       ! showCPU:       Whether or not the CPU time must be shown
       !                in the output to screen.

       logical :: monMachOrHMax
       logical :: showCPU

       ! monMassSliding:  Whether or not to monitor the mass flow of
       !                  the sliding interfaces.
       ! monMassFamilies: Whether or not the mass flow of at least
       !                  one family must be monitored.

       logical :: monMassSliding
       logical :: monMassFamilies
!
!      ******************************************************************
!      *                                                                *
!      * Variables to store the convergence info.                       *
!      *                                                                *
!      ******************************************************************
!
       ! nIterOld: Number of iterations done by a previous computation,
       !           i.e. before the restart. This value is 0 when started
       !           from scratch. NiterOld is an integer, because of
       !           cgns.
       ! nIterCur: Current number of iterations. Also niterCur is an
       !           integer, because of cgns.

       integer :: nIterOld, nIterCur

       ! convArray(0:nIterMax,nsps,nmon): 3D array to store the
       !                                  convergence histories.

       real(kind=cgnsRealType), dimension(:,:,:), allocatable :: &
                                                               convArray
!
!      ******************************************************************
!      *                                                                *
!      * Variables to store the time accurate history.                  *
!      * Only allocated for a time accurate computation.                *
!      *                                                                *
!      ******************************************************************
!
       ! nTimeStepsRestart:   Number of time steps taken in an earlier
       !                      unsteady computation from which a restart
       !                      is performed.
       ! timeStepUnsteady:    The current unsteady time step number;
       !                      restart is not taken into account.
       ! timeUnsteady:        Amount of physical time of the current
       !                      simulation; only relevant in unsteady mode.
       ! timeUnsteadyRestart: Amount of physical time from a previous
       !                      simulation from which a restart is
       !                      performed.

       integer(kind=intType) :: nTimeStepsRestart, timeStepUnsteady
       real(kind=realType)   :: timeUnsteady, timeUnsteadyRestart

       ! timeArray(nTimeMax):          Array to store the values of the
       !                               time at every time step.
       ! timeDataArray(nTimeMax,nMon): 2D array to store the variables
       !                               to be monitored at every time
       !                               step. No need to store a spectral
       !                               index here.

       real(kind=cgnsRealType), dimension(:), allocatable :: timeArray

       real(kind=cgnsRealType), dimension(:,:), allocatable :: &
                                                          timeDataArray

       ! writeGrid:    Whether or not a grid file must be written.
       ! writeVolume:  Idem for a volume solution file.
       ! writeSurface: Idem for a surface solution file.

       logical :: writeGrid, writeVolume, writeSurface

       end module monitor
