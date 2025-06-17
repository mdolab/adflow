module actuatorRegionData
    use constants

    type actuatorRegionType

        character(len=maxStringLen) :: famName
        integer(kind=intType) :: famID
        ! The block indexes of the cells included in this region
        integer(kind=intType), dimension(:, :), pointer :: cellIDs

        ! The total number of cells included this proc has
        integer(kind=intType) :: nCellIDs

        real(kind=realType), dimension(:, :), pointer :: force
        real(kind=realType), dimension(:), pointer :: heat

        integer(kind=intType), dimension(:), allocatable :: blkPtr

        ! Set the defaults for solution relaxation
        real(kind=realType) :: relaxStart = -one
        real(kind=realType) :: relaxEnd = -one
    end type actuatorRegionType

    integer(kind=intType), parameter :: nActuatorRegionsMax = 10
    type(actuatorRegionType), dimension(nActuatorRegionsMax), target :: actuatorRegions
    integer(kind=intTYpe) :: nActuatorRegions = 0

#ifndef USE_TAPENADE
    type(actuatorRegionType), dimension(nActuatorRegionsMax), target :: actuatorRegionsd
#endif
end module actuatorRegionData
