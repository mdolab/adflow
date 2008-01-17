!
!      ******************************************************************
!      *                                                                *
!      * File:          BCDataMod.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-24-2004                                      *
!      * Last modified: 04-13-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module BCDataMod
!
!      ******************************************************************
!      *                                                                *
!      * This local module contains the variables used to handle the    *
!      * the prescribed boundary data.                                  *
!      *                                                                *
!      ******************************************************************
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
       ! iBeg:     Starting index i-direction.
       ! jBeg:     Starting index j-direction.
       ! iEnd:     Ending index i-direction.
       ! jEnd:     Ending index j-direction.

       integer(kind=intType) :: nDataSet, cgnsBoco, nbcVar
       integer(kind=intType) :: iBeg, jBeg, iEnd, jEnd

       ! nFreestreamSubfaces:     Number of supersonic inflow subfaces
       !                          for which the flow field variables must
       !                          be set to the free stream values.
       ! freestreamSubfaces(:,3): the corresponding block, subface and
       !                          spectral solution.

       integer(kind=intType) :: nFreestreamSubfaces
       integer(kind=intType), dimension(:,:), allocatable :: &
                                              freestreamSubfaces

       ! nTurbFreestreamSubfaces:     Number of subfaces for which the
       !                              turbulence variables must be set
       !                              to the free stream values.
       ! turbFreestreamSubfaces(:,3): the corresponding block, subface
       !                              and spectral solution.

       integer(kind=intType) :: nTurbFreestreamSubfaces
       integer(kind=intType), dimension(:,:), allocatable :: &
                                              turbFreestreamSubfaces

       ! bcVarArray(:,:,nbcVar): array to store the interpolated
       !                         values for the active face.

       real(kind=realType), dimension(:,:,:), allocatable :: bcVarArray

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

       end module BCDataMod
