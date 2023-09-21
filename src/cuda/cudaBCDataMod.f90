! @File    :   cudaBCDataMod.f90
! @Desc    :   CUDA FORTRAN version of BCDataMod.F90. See that file for variable names and descriptions.

module cudaBCDataMod
    ! Local module with variables and subroutines to handle boundary data
    ! --- CUDA FORTRAN module---
    use cudafor    
    ! --- ADflow module ---
    use cgnsGrid
    use BCDataMod, only: nbcVarMax, nDataSet, cgnsBoco, nbcVar
    
    implicit none
    save
 
    integer(kind=intType), dimension(nbcVarMax), constant :: mass
    integer(kind=intType), dimension(nbcVarMax), constant :: length
    integer(kind=intType), dimension(nbcVarMax), constant :: time
    integer(kind=intType), dimension(nbcVarMax), constant :: temp
    integer(kind=intType), dimension(nbcVarMax), constant :: angle
    real(kind=realType), dimension(:, :, :), pointer, device :: xf
    real(kind=realType), dimension(3), constant :: axis
    real(kind=realType), dimension(3), constant :: radVec1
    real(kind=realType), dimension(3), constant :: radVec2
    logical, device :: axAssumed, massflowPrescribed
    logical, dimension(nbcVarMax), constant :: bcVarPresent
    type(cgnsBcDatasetType), pointer, dimension(:), device :: dataSet

! #ifndef USE_TAPENADE
!     type(cgnsBcDatasetType), pointer, dimension(:), device :: dataSetd
!     real(kind=realType), dimension(:, :, :), allocatable, device :: bcVarArrayd
! #endif

    contains

    subroutine copycudaBCDataMod
        ! Pointer assignment to have CPU vars known by the 'h_' variable name here
        use BCDataMod, only: h_mass => mass, h_length => length, h_time => time, h_temp => temp, h_angle => angle, &
        h_xf => xf, h_axis => axis, h_radVec1 => radVec1, h_radVec2 => radVec2, &
        h_axAssumed => axAssumed, h_massflowPrescribed => massflowPrescribed, &
        h_bcVarPresent => bcVarPresent, h_dataSet => dataSet
        ! h_dataSetd => dataSetd, &
        ! h_bcVarArrayd => bcVarArrayd

        implicit none

        ! Copy host (CPU) vars to GPU vars
        mass = h_mass;  length = h_length; time = h_time; temp = h_temp; angle = h_angle
        xf = h_xf
        axis = h_axis; radVec1 = h_radVec1; radVec2 = h_radVec2
        axAssumed = h_axAssumed; massflowPrescribed = h_massflowPrescribed
        bcVarPresent = h_bcVarPresent
        dataSet = h_dataSet
        ! dataSetd = h_dataSetd
        ! bcVarArrayd = h_bcVarArrayd

    end subroutine

end module cudaBCDataMod