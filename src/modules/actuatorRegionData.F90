module actuatorRegionData
  use constants

  type actuatorRegionType

     character(len=maxStringLen) :: famName
     integer(kind=intType) :: famID
     ! The block indexes of the cells included in this region
     integer(kind=intType), dimension(:, :), pointer :: cellIDs

     ! The total number of cells included this proc has
     integer(kind=intType) :: nCellIDs

     ! F is the total Force to be applied on this region
     real(kind=realType) :: F(3)

     ! T is the total torque to be applied on this regoin
     real(kind=realType) :: T
     real(kind=realType), dimension(3) :: axisVec

     ! Q is the total heat flux to be added on this regoin
     real(kind=realType) :: Q

     ! Volume is the total integrated volume of all cells (on all
     ! procs) included in this region
     real(kind=realType) :: volume

     integer(kind=intType), dimension(:), allocatable :: blkPtr

     ! Set the defaults for solution relaxation
     real(kind=realType) :: relaxStart = -one
     real(kind=realType) :: relaxEnd = -one
  end type actuatorRegionType

  integer(kind=intType), parameter :: nActuatorRegionsMax=10
  type(actuatorRegionType), dimension(nActuatorRegionsMax), target :: actuatorRegions
  integer(kind=intTYpe) :: nActuatorRegions=0

#ifndef USE_TAPENADE
  type(actuatorRegionType), dimension(nActuatorRegionsMax), target :: actuatorRegionsd
#endif
end module actuatorRegionData
