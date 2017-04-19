module actuatorRegionData
  use constants

  type actuatorRegionType

     character(len=maxStringLen) :: famName
     integer(kind=intType) :: famID
     ! The block indexes of the cells included in this region
     integer(kind=intType), dimension(:, :), pointer :: cellIDs
     ! The total number of cells included this proc has
     integer(kind=intType) :: nCellIDs

     ! The scaling factor to account for the edges of the region. 
     real(kind=realType), dimension(:), pointer :: Factor

     ! F is the total Force to be applied on this region
     real(kind=realType) :: F(3)

     ! T is the total torque to be applied on this regoin
     real(kind=realType) :: T

     ! Volume is the total integrated volume of all cells (on all
     ! procs) inluleded in this region
     real(kind=realType) :: volume

     ! The physical distance 
     real(kind=realType) :: smoothDistance
     integer(kind=intType), dimension(:), allocatable :: blkPtr
  end type actuatorRegionType

  integer(kind=intType), parameter :: nActuatorRegionsMax=10
  type(actuatorRegionType), dimension(nActuatorRegionsMax), target :: actuatorRegions
  integer(kind=intTYpe) :: nActuatorRegions=0
end module actuatorRegionData
