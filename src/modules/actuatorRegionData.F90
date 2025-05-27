module actuatorRegionData
    use constants

    type actuatorRegionType

        character(len=maxStringLen) :: famName, actType
        integer(kind=intType) :: famID
        ! The block indexes of the cells included in this region
        integer(kind=intType), dimension(:, :), pointer :: cellIDs

        ! The tangent unit-vectors of the cells included in this region
        real(kind=realType), dimension(:, :), pointer :: cellTangentials
        ! The radii of the cells included in this region
        real(kind=realType), dimension(:), pointer :: cellRadii

        ! The total number of cells included this proc has
        integer(kind=intType) :: nCellIDs

        ! the force vector to be applied on this region
        ! this is equal to torque * axisVec
        real(kind=realType) :: force(3)

        ! thrust is the total thrust magnitude to be applied on this region
        real(kind=realType) :: thrust
        real(kind=realType) :: swirlFact
        real(kind=realType) :: mDistribParam
        real(kind=realType) :: nDistribParam
        real(kind=realType) :: distribPDfactor
        real(kind=realType) :: innerZeroThrustRadius
        real(kind=realType) :: spinnerRadius
        real(kind=realType) :: rootDragFactor

        ! magnitude of the total torque to be applied on this region
        real(kind=realType) :: torque
        ! vector that determines the direction of the applied torque and force
        real(kind=realType), dimension(3) :: axisVec
        real(kind=realType), dimension(:, :), pointer :: thrustVec
        real(kind=realType), dimension(:, :), pointer :: swirlVec

        ! total heat flux to be added on this regoin
        real(kind=realType) :: heat

        ! Volume is the total integrated volume of all cells (on all
        ! procs) included in this region
        real(kind=realType) :: volume
        real(kind=realType) :: volLocal

        real(kind=realType) :: totalThrustSum
        real(kind=realType) :: totalSwirlSum

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
