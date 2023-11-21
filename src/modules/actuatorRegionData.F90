module actuatorRegionData
    use constants

    type actuatorRegionType

        character(len=maxStringLen) :: famName
        integer(kind=intType) :: famID
        ! The block indexes of the cells included in this region
        integer(kind=intType), dimension(:, :), pointer :: cellIDs

        ! The total number of cells included this proc has
        integer(kind=intType) :: nCellIDs

        ! the force vector to be applied on this region
        ! this is equal to torque * axisVec
        real(kind=realType) :: force(3)

        ! magnitude of the total torque to be applied on this region
        real(kind=realType) :: torque
        ! vector that determines the direction of the applied torque
        real(kind=realType), dimension(3) :: axisVec

        ! total heat flux to be added on this regoin
        real(kind=realType) :: heat

        ! Volume is the total integrated volume of all cells (on all
        ! procs) included in this region
        real(kind=realType) :: volume
        real(kind=realType) :: volLocal

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
