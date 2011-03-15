!
!      ******************************************************************
!      *                                                                *
!      * File:          readParamFile.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-11-2002                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readParamFile
!
!      ******************************************************************
!      *                                                                *
!      * readParamFile sets the values of the input parameters.         *
!      * First the default values are set whenever possible, then the   *
!      * input parameters are read from the given parameter file and    *
!      * finally a check is made to be sure that all desired parameters *
!      * have been specified.                                           *
!      *                                                                *
!      * All processors read the parameter file. Due to the small size  *
!      * this causes no overhead.                                       *
!      *                                                                *
!      * If the parameter file cannot be found, a template is created.  *
!      * This file must be adjusted by the user to create a valid       *
!      * parameter file.                                                *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       use allInputParam
       implicit none
!
!      Local variables
!
       integer, parameter :: readUnit = 32

       integer :: ios, ierr

       character (len=2*maxStringLen) :: errorMessage
       character (len=512)            :: string
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the input parameters to their default values

       call setDefaultValues

       ! Open the given file for reading

       open(unit=readUnit, file=paramFile, status="old", &
            action="read", iostat=ios)

       ! It is possible that in case the parameter file was not specified,
       ! processor 0 creates the default file and in the mean time other
       ! processors can open this default file. To avoid this a barrier
       ! is implemented here.

       call mpi_barrier(SUmb_comm_world, ierr)

       ! In case the file could not be opened, a template will be
       ! created. As only 1 file is needed, processor 0 performs this
       ! task.

       if(ios /= 0) then

         if(myID == 0) then

           ! Creat the template parameter file.

           call createTemplate

           ! Create the error message and terminate the program.

           write(errorMessage,*) "Parameter file ", trim(paramFile), &
                                 " not found. A template file has &
                                 &been created"
           call terminate("readParamFile", errorMessage)

         endif

         ! As only processor 0 writes the error message, the other
         ! processors should wait before they are killed.

         call mpi_barrier(SUmb_comm_world, ierr)

       endif

       ! Loop to read the data

       do

         ! Read a string from the file. In case the end of the file
         ! has been reached, exit the loop.

         read(unit=readUnit, fmt="(a512)", iostat=ios) string
         if(ios /= 0) exit

         ! Analyze the string and extract the possible info

         call analyzeString(string)

       enddo

       ! Close the file

       close(unit=readUnit)

       ! Check if all the desired input parameters were specified and
       ! print warnings if some irrelevant ones are specified.

       call checkInputParam

       ! Read the cp curve fits from file if a variable cp model
       ! must be used.

       if(cpModel == cpTempCurveFits) call readCpTempCurveFits

       ! Determine the number of governing equations and set the
       ! corresponding parameters accordingly.

       call setEquationParameters

       ! Extract the multigrid info from the string.

       call extractMgInfo

       ! Set some iteration parameters, depending on the governing
       ! equations to be solved and the iterative strategy to be used.

       call setIterationParam

       ! If no monitoring variables were specified, set the default set.
       ! Idem for the surface and volume output variables.

       if(.not. monitorSpecified)    call defaultMonitor
       if(.not. surfaceOutSpecified) call defaultSurfaceOut
       if(.not. volumeOutSpecified)  call defaultVolumeOut

       ! Check the monitoring and output variables.

       call checkMonitor
       call checkOutput

       end subroutine readParamFile
