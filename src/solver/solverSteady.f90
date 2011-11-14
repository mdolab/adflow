!
!      ******************************************************************
!      *                                                                *
!      * File:          solverSteady.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-12-2003                                      *
!      * Last modified: 03-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine solverSteady
!
!      ******************************************************************
!      *                                                                *
!      * solverSteady solves the steady equations for the multigrid     *
!      * level groundLevel, which can be found in the module            *
!      * iteration.                                                     *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       call solveState

       end subroutine solverSteady
