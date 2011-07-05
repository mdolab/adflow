!
!      ******************************************************************
!      *                                                                *
!      * File:          localModules.F90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-11-2002                                      *
!      * Last modified: 03-23-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module localMG
!
!      ******************************************************************
!      *                                                                *
!      * Locally used module for the storage of the string describing   *
!      * the cycle strategy. Inside the code this info is converted to  *
!      * an equivalent integer array. The string is only used to be     *
!      * able to define predefined cycling strategies, like sg, 2v, 4w, *
!      * etc.                                                           *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the string, which stores the multigrid cycling   *
!      * strategy.                                                      *
!      *                                                                *
!      ******************************************************************
!
       character (len=maxStringLen) :: mgDescription

       end module localMG

!      ==================================================================

       module allInputParam
!
!      ******************************************************************
!      *                                                                *
!      * Locally used module which includes all input parameter modules *
!      * and some logicals to see whether or not monitoring variables,  *
!      * surface output variables and extra volume output variables     *
!      * were specified.                                                *
!      *                                                                *
!      ******************************************************************
!
       use inputDiscretization
       use inputIO
       use inputIteration
       use inputMotion
       use inputOverset
       use inputParallel
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use inputVisualization
       implicit none
       save

       ! Set the parameter none, which is used as a check to see
       ! whether or not some key parameters were specified.

       integer(kind=intType), parameter :: none = 0

       ! monDturb:             Whether or not the turbulent residuals
       !                        must be monitored. This must be done via
       !                        this construction, because during the
       !                        reading of the monitoring variables the
       !                        turbulence model might not be known.
       ! monitorSpecified:     Whether or not the monitoring variables
       !                        were specified.
       ! surfaceOutSpecified: Whether or not the surface output
       !                        variables were specified.
       ! volumeOutSpecified:  Whether or not the volume output
       !                        variables were specified.

       logical :: monDturb
       logical :: monitorSpecified
       logical :: surfaceOutSpecified
       logical :: volumeOutSpecified

       ! liftDirSpecified: Whether or not the lift direction was
       !                     specified.

       logical :: liftDirSpecified

       end module allInputParam
