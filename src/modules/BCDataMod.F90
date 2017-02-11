module BCDataMod
  !
  !       This local module contains the variables and subroutine to     
  !       handle the prescribed boundary data.                           
  !
  use cgnsGrid
  implicit none
  save

  ! nbcVarMax: Parameter, which defines the maximum number of
  !              prescribed variables for a boundary.

  integer, parameter :: nbcVarMax = 21

  ! mass(nbcVarMax):   Unit of mass for the prescribed data.
  ! length(nbcVarMax): Unit of length for the prescribed data.
  ! time(nbcVarMax):   Unit of time for the prescribed data.
  ! temp(nbcVarMax):   Unit of temperature for the prescribed data.
  ! angle(nbcVarMax):  Unit of angle for the prescribed data.

  integer, dimension(nbcVarMax) :: mass, length, time, temp, angle

  ! nDataSet: Number of data sets present for the active face.
  ! cgnsBoco: The corresponding boundary face in the cgns block.
  ! nbcVar:   Theoretically possible number of variables
  !           prescribed for this face.

  integer(kind=intType) :: nDataSet, cgnsBoco, nbcVar

  ! xf(:,:,:): pointer to the coordinates of the block face to
  !            which the active subface belongs.

  real(kind=realType), dimension(:,:,:), pointer :: xf

  ! axis(3):    Axial unit vector in the local cylindrical
  !             coordinate system.
  ! radVec1(3): First radial unit vector in the local cylindrical
  !             coordinate system.
  ! radVec2(3): Second radial unit vector.

  real(kind=realType), dimension(3) :: axis, radVec1, radVec2

  ! bcVarNames(nbcVarMax): The cgns names of the possible
  !                        prescribed variables.

  character(len=maxCGNSNameLen), dimension(nbcVarMax) :: bcVarNames

  ! axAssumed:               Whether or not the x-axis is assumed
  !                          to be the axial direction.
  ! massflowPrescribed:      Whether or not subsonic inflow boundaries
  !                          are present with prescribed massflow.
  ! bcVarPresent(nbcVarMax): Whether or not the possible
  !                          variables are actually prescribed.

  logical :: axAssumed, massflowPrescribed
  logical, dimension(nbcVarMax) :: bcVarPresent

  ! dataSet: Pointer for the data sets of the corresponding cgns
  !          boundary face.

  type(cgnsBcDatasetType), pointer, dimension(:) :: dataSet
#ifndef USE_TAPENADE
  type(cgnsBcDatasetType), pointer, dimension(:) :: dataSetd
  real(kind=realType), dimension(:,:,:), allocatable :: bcVarArrayd
#endif 
end module BCDataMod
