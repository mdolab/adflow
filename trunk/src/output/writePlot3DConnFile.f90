!
!      ******************************************************************
!      *                                                                *
!      * File:          writePlot3DConnFile.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-29-2006                                      *
!      * Last modified: 03-29-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writePlot3DConnFile
!
!      ******************************************************************
!      *                                                                *
!      * writePlot3DConnFile writes the Plot3D connectivity file, if    *
!      * desired. Only processor 0 does the writing.                    *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use inputIO
       use IOModule
       use outputMod
       implicit none
!
!      Local parameter.
!
       integer, parameter :: writeunit = 40
!
!      Local variables.
!
       integer :: ierr

       character(len=512) :: string
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if the connectivity file should not be
       ! written or if I'm not processor 0.

       if(.not. writePlot3DConn) return
       if(myID /= 0)             return

       call terminate("writePlot3DConnFile", "Not implemented yet")

       ! Open the plot3D connectivity file for writing and check if it
       ! went okay.

       open(unit=writeunit, file=plot3DConnFile, status="unknown", &
            action="write", iostat=ierr)
       if(ierr /= 0) then
         write(string,101) trim(plot3DConnFile)
         call terminate("writePlot3DConnFile", string)
       endif




       ! Format statements.

 101   format("Plot3D connectivity file",1x,a,1x, &
              "could not be opened.")

       end subroutine writePlot3DConnFile
