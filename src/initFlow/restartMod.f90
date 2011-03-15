!
!      ******************************************************************
!      *                                                                *
!      * File:          restartMod.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-20-2003                                      *
!      * Last modified: 10-11-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module restartMod
!
!      ******************************************************************
!      *                                                                *
!      * This local module contains variables used when reading the     *
!      * restart file(s).                                               *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      *                 Variables in this module.                      *
!      *                                                                *
!      *           Variables used by both CGNS and Plot3D.              *
!      *                                                                *
!      ******************************************************************
!
       ! nVar:     Number of variables stored in the solution file.
       ! solID:    Loop variables for the number of solutions to be read.
       ! varNames: Variable names, sorted in increasing order,
       !           of the variables.

       integer :: nVar
       integer(kind=intType) :: solID

       character(len=maxCGNSNameLen), allocatable, dimension(:) :: &
                                                                 varNames

       ! rhoScale: Scale factor for the density.
       ! velScale: Scale factor for the velocity.
       ! pScale:   Scale factor for the pressure.
       ! muScale:  Scale factor for the molecular viscosity.

       real(kind=realType) :: rhoScale, velScale, pScale, muScale

       ! nSolsRead:            Number of solution files to read.
       ! solFiles(nSolsRead):  Names of the solution files to be read.
       ! interpolSpectral:     Whether or not to interpolate the 
       !                       coordinates/solutions for the time
       !                       spectral mode.
       ! copySpectral:         Whether or not to copy the solutions
       !                       for the time spectral mode.

       integer(kind=intType) :: nSolsRead
       logical ::               interpolSpectral, copySpectral

       character(len=maxStringLen), dimension(:), allocatable :: solFiles
!
!      ******************************************************************
!      *                                                                *
!      *                 Variables only used by CGNS.                   *
!      *                                                                *
!      ******************************************************************
!
       ! rangeMin(3):    Lower index in i, j and k direction of the
       !                 range to be read.
       ! rangeMax(3):    Upper index in i, j and k direction of the
       !                 range to be read.
       ! varTypes(nVar): Variable types of the variables stored.

       integer, dimension(3)              :: rangeMin, rangeMax
       integer, allocatable, dimension(:) :: varTypes

       ! cgnsInd:  File index of the CGNS file.
       ! cgnsBase: Base of the CGNS file, always set to 1.
       ! cgnsZone: Zone ID in the CGNS file.
       ! cgnsSol:  Solution ID in the zone ID of the CGNS file.
       ! location: Location where the variables are stored in CGNS.
       !           Supported possibilities are Vertex and CellCentered.

       integer :: cgnsInd, cgnsBase, cgnsZone, cgnsSol, location

       ! zoneNumbers: Corresponding zoneNumbers of the sorted
       !              zoneNames.
       ! zoneNames:   Zone names, sorted in increasing order, of the
       !              zones in the CGNS restart file.

       integer(kind=intType), allocatable, dimension(:) :: zoneNumbers

       character(len=maxCGNSNameLen), allocatable, dimension(:) :: &
                                                               zoneNames

       ! buffer(2:il,2:jl,2:kl): Buffer to read and store the cell
       !                         centered values.
       ! bufferVertex(:):        Additional buffer needed to read
       !                         vertex data and transform them into
       !                         cell centered values.

       real(kind=cgnsRealType), dimension(:,:,:), allocatable :: buffer
       real(kind=cgnsRealType), dimension(:,:,:), allocatable :: &
                                                            bufferVertex
!
!      ******************************************************************
!      *                                                                *
!      *               Variables only used by Plot3D.                   *
!      *                                                                *
!      ******************************************************************
!
       ! sorted2Or: Mapping of the sorted names to the original
       !            names.

       integer(kind=intType), dimension(:), allocatable :: sorted2Or

       ! sizeVolumeSol:   Size in bytes to store one cell centered
       !                  variable for the entire grid.
       ! sizeConvHistory: Size in bytes of the convergence history.
       ! sizeHeader:      Size in bytes of the file header, i.e. 
       !                  everything before the solution variables.

       integer(kind=mpi_offset_kind) :: sizeHeader
       integer(kind=mpi_offset_kind) :: sizeVolumeSol, sizeConvHistory

       end module restartMod
