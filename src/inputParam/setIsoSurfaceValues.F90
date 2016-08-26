subroutine initializeIsoSurfaceVariables(values, nValues)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * isoVariables extracts from the given string the extra          *
  !      * iso surface variables to be written to the solution file.      *
  !      *                                                                *
  !      ******************************************************************
  !
  use constants
  use extraOutput, only : isoValues, isoSurfaceNames, nIsoSurface
  implicit none
  !
  !      Subroutine arguments.
  !
  real(kind=realType), dimension(nValues), intent(in) :: values
  integer(kind=intType), intent(in) :: nValues

  ! Basically just copy into module
  if (allocated(isoValues)) then
     deallocate(isoValues)
  end if

  if (allocated(isoSurfaceNames)) then
     deallocate(isoSurfaceNames)
  end if

  nIsoSurface = nValues
  allocate(isoValues(nIsoSurface))
  allocate(isoSurfaceNames(nIsoSurface))

  isoValues = values

end subroutine initializeIsoSurfaceVariables

subroutine setIsoSurfaceVariable(variable, iVar)

  ! Set variable to iVar. initializeIsoSurfaceVariables MUST be called
  ! first with the desired number of values to set. 

  use constants
  use cgnsNames
  use extraOutput
  use communication, only : myID
  use utils, only : EChk
  implicit none
  !
  !      Subroutine arguments.
  !
  character(len=*), intent(in):: variable
  integer(kind=intType) :: iVar

  select case (variable)
  case("rho") 
     isoSurfaceNames(iVar) = cgnsDensity                 
  case("vx")  
     isoSurfaceNames(iVar) = cgnsVelX
  case("vy")
     isoSurfaceNames(iVar) = cgnsVelY
  case("vz")
     isoSurfaceNames(iVar) = cgnsVelZ
  case("P")
     isoSurfaceNames(iVar) = cgnsPressure
  case ("mx")
     isoSurfaceNames(iVar) = cgnsMomX
  case ("my")
     isoSurfaceNames(iVar) = cgnsMomY
  case ("mz")
     isoSurfaceNames(iVar) = cgnsMomZ
  case ("rvx")
     isoSurfaceNames(iVar) = cgnsRelVelX
  case ("rvy")
     isoSurfaceNames(iVar) = cgnsRelVelY
  case ("rvz")
     isoSurfaceNames(iVar) = cgnsRelVelZ
  case ("rhoe")
     isoSurfaceNames(iVar) = cgnsEnergy
  case ("temp")
     isoSurfaceNames(iVar) = cgnsTemp
  case ("vort")
     isoSurfaceNames(iVar) = cgnsVortMagn
  case ("vortx")
     isoSurfaceNames(iVar) = cgnsVortX
  case ("vorty")
     isoSurfaceNames(iVar) = cgnsVortY
  case ("vortz")
     isoSurfaceNames(iVar) = cgnsVortZ
  case ("cp")
     isoSurfaceNames(iVar) = cgnsCp
  case ("mach")
     isoSurfaceNames(iVar) = cgnsMach
  case ("rmach")
     isoSurfaceNames(iVar) = cgnsRelMach
  case ("macht")
     isoSurfaceNames(iVar) = cgnsMachTurb
  case ("ptloss")
     isoSurfaceNames(iVar) = cgnsPTotLoss
  case ("eddy")
     isoSurfaceNames(iVar) = cgnsEddy
  case ("eddyratio")
     isoSurfaceNames(iVar) = cgnsEddyRatio
  case ("dist")
     isoSurfaceNames(iVar) = cgnsWallDist
  case ("resrho")
     isoSurfaceNames(iVar) = cgnsResRho
  case("shock")
     isoSurfaceNames(iVar) = cgnsShock
  case("filteredShock")
     isoSurfaceNames(iVar) = cgnsFilteredShock


  case default

     if(myID == 0) Then
        print *,'Error: ', variable, 'cannot be used as an isoSurface'
     end if
     call EChk(-99, __FILE__, __LINE__)
  end select
end subroutine setIsoSurfaceVariable

