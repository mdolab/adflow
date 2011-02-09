!
!      ******************************************************************
!      *                                                                *
!      * File:          initCurveFitDataSae.f90                         *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 08-20-2004                                      *
!      * Last modified: 04-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initCurveFitDataSae
!
!      ******************************************************************
!      *                                                                *
!      * initCurveFitDataSae contains the curve fit constants for       *
!      * the wall function data for the Spalart-Allmaras (Edwards       *
!      * modification) turbulence model.                                *
!      *                                                                *
!      ******************************************************************
!
       use flowVarRefState
       use paramTurb
       implicit none
!
!      Local variables.
!
   !   integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       call terminate("initCurveFitDataSae", &
                      "Not implemented yet")

       end subroutine initCurveFitDataSae
