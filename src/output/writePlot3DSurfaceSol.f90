!
!      ******************************************************************
!      *                                                                *
!      * File:          writePlot3DSurfaceSol.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-30-2005                                      *
!      * Last modified: 11-01-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writePlot3DSurfaceSol
!
!      ******************************************************************
!      *                                                                *
!      * writePlot3DSurfaceSol writes a surface solution file when the  *
!      * plot3D format is used for the grids. The solution file written *
!      * in this routine is not a plot3D solution file. It is just an   *
!      * internally used format.                                        *
!      *                                                                *
!      ******************************************************************
!
       use communication
       implicit none
!
!      Local variables.
!

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number and names of the solution files.

       call surfSolFileNamesWrite



       if(myID == 0) then
         print "(a)", "#"
         print "(a)", "#                 Warning"
         print "(a)", "# writePlot3DSurfaceSol is not implemented yet"
         print "(a)", "#"
       endif

       end subroutine writePlot3DSurfaceSol
