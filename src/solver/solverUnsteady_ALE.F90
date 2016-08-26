!
!       File:          solverUnsteadyBDF.F90                           
!       Author:        Edwin van der Weide                             
!       Starting date: 03-13-2003                                      
!       Last modified: 11-21-2007                                      
!
subroutine solverUnsteady_ALE(alecallback_python)
  !
  !       solverUnsteadyMD  solves the unsteady equations using the BDF  
  !       schemes for the multigrid level groundLevel. However, it       
  !       this version includes a call-back to python for aero-elastic   
  !       simulations                                                    
  !       Modified for ALE scheme, by HDN                                
  !
  use inputIteration
  use inputUnsteady
  use blockPointers
  use inputDiscretization
  use iteration
  use killSignals
  use monitor
  use communication 
  use section
  use inputMotion
  use inputTimeSpectral
  implicit none
  !
  !      Local variables.
  
  integer(kind=intType) :: iter, nTimeSteps
  external alecallback_python

  !       Begin execution                                                

  call solverUnsteadyWrapBegin

  ! Set time steps
  nTimeSteps = nTimeStepsCoarse
  if(groundLevel == 1) nTimeSteps = nTimeStepsFine

  timeLoop: do iter=1,nTimeSteps
     
     timeStepUnsteady = timeStepUnsteady + 1
     timeUnsteady     = timeUnsteady     + deltaT

     ! if(myID == 0) call unsteadyHeader

     call alecallback_python(timeUnsteady)

     call solverUnsteadyWrapInLoop

  enddo timeLoop

  call solverUnsteadyWrapEnd

end subroutine solverUnsteady_ALE
