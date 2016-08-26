subroutine initFlowRestart
  !
  !       initFlowRestart loads restart information from the restart     
  !       file into the state variables.                                 
  !
  use constants
  use IOModule, only : IOVar
  use restartMod, only : solFiles, copySpectral, halosRead, interpolSpectral
  use utils, only : terminate
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType) :: ierr

  ! Initialize halosRead to .false. This will be overwritten
  ! if halo values are read during the restart.

  halosRead = .false.


  ! Determine the number and names of the solution files and
  ! determine the contents of IOVar, the data structure to
  ! generalize the IO.
  
  call determineSolFileNames
  call setIOVar
  
  ! Determine the format of the files and read them.
  ! Note that halosRead is possibly overwritten in the
  ! folloing select case statement below
  
  call readRestartFile()
  
  ! Copy or interpolate the spectral solution, if needed.
  
  if( copySpectral )     call copySpectralSolution
  if( interpolSpectral ) call interpolateSpectralSolution
  
  ! Release the memory of the file names and IOVar.

  deallocate(solFiles, IOVar, stat=ierr)
  if(ierr /= 0)                &
       call terminate("initFlow", &
       "Deallocation failure for solFiles and IOVar")
  
  ! At the moment the pressure is stored at the location of the
  ! total energy. Copy the pressure to its own arrays and
  ! compute the total energy.
  
  call setPressureAndComputeEnergy(halosRead)

  ! Initialize the halo cells if a restart is performed and the
  ! entire flow field if this is not the case.

  call initializeHalos(halosRead)
  
  ! Initialize the dependent flow variables and the halo values.
  call initDepvarAndHalos(halosRead)

end subroutine initFlowRestart
