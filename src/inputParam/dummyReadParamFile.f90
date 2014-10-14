
       subroutine dummyreadParamFile
!
!      ******************************************************************
!      *                                                                *
!      * This subroutine is the same as readParamFile EXCEPT it does not*
!      * read the actual file. Values are set diectly from python for   *
!      * all the options and then this file is run.                     *
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
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!

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
       if(.not. isoOutSpecified)     call defaultIsoOut

       ! Check the monitoring and output variables.

       call checkMonitor
       call checkOutput


     end subroutine dummyreadParamFile
