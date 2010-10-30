!
!      ******************************************************************
!      *                                                                *
!      * File:          readDensityPlot3D.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-21-2005                                      *
!      * Last modified: 10-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readDensityPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readDensityPlot3D reads the density from the restart file when *
!      * Plot3D format is used. It checks if the density is actual      *
!      * present, determines the position and calls the general reading *
!      * routine for Plot3D solution files.                             *
!      *                                                                *
!      ******************************************************************
!
       use cgnsNames
       use communication
       use IOModule
       use restartMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, p
!
!      Function definitions.
!
       integer(kind=intType) :: bsearchStrings
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Find the index of the density in the solution file. If not
       ! found print an error message and exit.

       p  = nVar
       nn = bsearchStrings(cgnsDensity, varNames, p)
       if(nn == 0) then
         if(myID == 0)                         &
           call terminate("readDensityPlot3D", &
                          "Density not found in restart file")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       nn = sorted2Or(nn)

       ! Determine the offset from the beginning of the file where this
       ! solution record starts and call the general reading routine.

       P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol
       p = irho

       call readPlot3DVar(solID, rhoScale, p)

       end subroutine readDensityPlot3D
