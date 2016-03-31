!
!      ******************************************************************
!      *                                                                *
!      * File:          initFlowRestart.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-06-2003                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine initFlowRestart
  !
  !      ******************************************************************
  !      *                                                                *
  !      * initFlowRestart loads restart information from the restart     *
  !      * file into the state variables.                                 *
  !      *                                                                *
  !      ******************************************************************
  !
  use bleedFlows
  use communication
  use inputIO
  use IOModule
  use restartMod
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
  
  select case (fileFormatRead)
  case (cgnsFormat)
     call readRestartFile()
     
  case (plot3DFormat)
     call readRestartFilePlot3D()
  end select

  ! Copy or interpolate the spectral solution, if needed.
  
  if( copySpectral )     call copySpectralSolution
  if( interpolSpectral ) call interpolateSpectralSolution
  
  ! Release the memory of the file names and IOVar.

  deallocate(solFiles, IOVar, stat=ierr)
  if(ierr /= 0)                &
       call returnFail("initFlow", &
       "Deallocation failure for solFiles and IOVar")
  
  ! At the moment the pressure is stored at the location of the
  ! total energy. Copy the pressure to its own arrays and
  ! compute the total energy.
  
  call setPressureAndComputeEnergy(halosRead)

  ! If no halo values were read and inflow or outflow bleed
  ! regions are present, print a warning that the restart will
  ! not be perfect.
  
  if((.not. halosRead) .and. myID == 0) then
     if(nInflowBleeds > 0 .or. nOutflowBleeds > 0) then
        print "(a)", "#"
        print "(a)", "#                 Warning"
        print "(a)", "# Halo values could not read during the &
             &restart."
        print "(a)", "# The restart is not consistent, because &
             &the halo info is needed for a"
        print "(a)", "# consistent initialization of the bleed &
             &regions."
        print "(a)", "# A constant extrapolation is used instead."
        print "(a)", "#"
     endif
  endif
  
  ! Initialize the bleed regions from the halos if the halos
  ! were read.

  if( halosRead ) call initBleedsFromHalos

  ! Initialize the halo cells if a restart is performed and the
  ! entire flow field if this is not the case.

  call initializeHalos(halosRead)
  
  ! Initialize the dependent flow variables and the halo values.
  call initDepvarAndHalos(halosRead)

end subroutine initFlowRestart
