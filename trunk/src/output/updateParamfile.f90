!
!      ******************************************************************
!      *                                                                *
!      * File:          updateParamfile.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-11-2003                                      *
!      * Last modified: 03-29-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine updateParamfile
!
!      ******************************************************************
!      *                                                                *
!      * updateParamfile updates the parameter file, such that a        *
!      * restart can be done automatically. This is only done if        *
!      * autoParameterUpdate is .true.. Only processor 0 performs       *
!      * this task.                                                     *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use inputIO
       use inputIteration
       use inputPhysics
       use iteration
       use monitor
       use outputMod
       implicit none
!
!      Local parameter
!
       integer, parameter :: writeUnit = 35
!
!      Local variables
!
       integer :: ios

       character(len=2*maxStringLen) :: errorMessage
       character(len=8)              :: integerString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if no grid and solution file has
       ! been written.

       if(.not. (writeGrid .or. writeVolume) ) return

       ! Check if the parameter file must be updated and if I am
       ! processor 0; only processor 0 performs the update.

       testUpdate: if(autoParameterUpdate .and. myID == 0) then

         ! Open the parameter file for appending.

         open(unit=writeUnit, file=paramFile, status="old", &
              action="write", position="append", iostat=ios)
         if(ios /= 0) then
           write(errorMessage,*) "Parameter file ", trim(paramFile), &
                                 " could not be updated"
           call terminate("updateParamfile", errorMessage)
         endif

         ! Write a new line, such that the output looks nicer.

         write(writeUnit,*)

         ! First write a comment line if a volume file is written.
         ! This depends on the mode of the equations being solved.

         if( writeVolume ) then

           select case (equationMode)
             case (steady, timeSpectral)
               write(integerString,"(i8)") iterTot + nIterOld
               integerString = adjustl(integerString)
               write(writeUnit,100) trim(integerString)
 100           format(" # Restart after",1X,A,1X,"iterations")

             case (unsteady)
               write(integerString,"(i8)") timeStepUnsteady + &
                                           nTimeStepsRestart
               integerString = adjustl(integerString)
               write(writeUnit,101) trim(integerString)
 101           format(" # Restart after",1X,A,1X,"time steps")
           end select

         endif

         ! Write the format, if needed.

         if( writeFormatInParam ) then
           select case(fileFormatWrite)
             case (cgnsFormat)
               write(writeUnit,*) "      File format read: CGNS"

             case (plot3DFormat)
               write(writeUnit,*) "      File format read: PLOT3D"
           end select
         endif

         ! In case a grid file has been written, write the name
         ! of the new grid file. If multiple files were written
         ! only the first name is relevant for the restart.

         if( writeGrid ) &
           write(writeUnit,*) "             Grid file: ", &
                              trim(gridFileNames(1))

         ! Write the name of the new restart file, which is the
         ! current solution file. Only if a solution was written.

         if( writeVolume ) then
           write(writeUnit,*) "          Restart file: ", &
                              trim(volSolFileNames(1))

           ! Write the MG start level and restart.

           write(writeUnit,*) " Multigrid start level: 1"
           write(writeUnit,*) "               Restart: yes"
         endif

         ! Close the file.

         close(unit=writeUnit)

       endif testUpdate

       end subroutine updateParamfile
