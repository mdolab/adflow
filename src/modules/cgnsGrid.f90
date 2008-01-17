!
!      ******************************************************************
!      *                                                                *
!      * File:          cgnsGrid.f90                                    *
!      * Author:        Edwin van der Weide, Steve Repsher,             *
!      *                Seonghyeon Hahn                                 *
!      * Starting date: 12-17-2002                                      *
!      * Last modified: 10-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       module cgnsGrid
!
!      ******************************************************************
!      *                                                                *
!      * This module contains the derived data type for storing the     *
!      * information of the original cgns grid file. Information stored *
!      * is number of blocks, block sizes, zone names, etc. this info   *
!      * is needed again when the solution is written to file. Remember *
!      * that the original blocks may be split to obtain a better       *
!      * load balance. Note that this info is stored on all processors. *
!      *                                                                *
!      * Apart from the derived data type for the cgns blocks, this     *
!      * module also contains the name of the base and the physical     *
!      * dimensions of the problem.                                     *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived datatype to store the actual     *
!      * data of the boundary conditions.                               *
!      *                                                                *
!      ******************************************************************
!
       type cgnsBCDataArray

         ! The units in which the data is specified.

         integer :: mass, len, time, temp, angle

         ! The name of the array.

         character(len=maxCGNSNameLen) :: arrayName

         ! Number of dimensions for which the data is specified.

         integer :: nDimensions

         ! Number of data points of every dimensions. upper limit is
         ! three, although for BC data the maximum is usually 2.

         integer(kind=intType), dimension(3) :: dataDim

         ! The actual data. Assumed is that only floating point data
         ! is prescribed and not integer or character data. Note that
         ! dataArr is a 1D array even if the data is multi-dimensional.

         real(kind=realType), pointer, dimension(:) :: dataArr

       end type cgnsBCDataArray
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived datatype to store the prescribed *
!      * boundary data for a boundary subface.                          *
!      *                                                                *
!      ******************************************************************
!
       type cgnsBCDatasetType

         ! Name of the dataset.

         character(len=maxCGNSNameLen) :: datasetName

         ! Boundary condition type.

         integer :: BCType

         ! The number of Dirichlet arrays in the data set.

         integer(kind=intType) :: nDirichletArrays

         ! The number of Neumann arrays in the data set.

         integer(kind=intType) :: nNeumannArrays

         ! The Dirichlet arrays.

         type(cgnsBCDataArray), pointer, dimension(:) :: dirichletArrays

         ! The Neumann arrays.

         type(cgnsBCDataArray), pointer, dimension(:) :: neumannArrays

       end type cgnsBCDatasetType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type to store cgns 1 to 1   *
!      * block to block, i.e. continuous grid lines across block        *
!      * boundaries, connectivities.                                    *
!      *                                                                *
!      ******************************************************************
!
       type cgns1to1ConnType

         ! Name of the interface.

         character(len=maxCGNSNameLen) :: connectName

         ! Name of the zone/block interfacing with the current zone/block.

         character(len=maxCGNSNameLen) :: donorName

         ! Zone/block ID of the zone/block interfacing with the current
         ! zone/block.

         integer(kind=intType) :: donorBlock

         ! Range of points of this subface.

         integer(kind=intType) :: iBeg, jBeg, kBeg
         integer(kind=intType) :: iEnd, jEnd, kEnd

         ! Range of points for the donor block.

         integer(kind=intType) :: diBeg, djBeg, dkBeg
         integer(kind=intType) :: diEnd, djEnd, dkEnd

         ! Short hand notation defining the relative orientation of the
         ! two zones.

         integer(kind=intType) :: l1, l2, l3

         ! Whether or not the subface is a periodic boundary.

         logical :: periodic

         ! The center of rotation for a periodic boundary.

         real(kind=realType), dimension(3) :: rotationCenter

         ! The rotation angles for a periodic boundary.

         real(kind=realType), dimension(3) :: rotationAngles

         ! The translation vector for a periodic boundary.

         real(kind=realType), dimension(3) :: translation

       end type cgns1to1ConnType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type to store cgns          *
!      * non-matching abutting block to block connectivities.           *
!      *                                                                *
!      ******************************************************************
!
       type cgnsNonMatchAbuttingConnType

         ! Number of donor blocks. It is possible that the subface
         ! abuts multiple donor blocks.

         integer(kind=intType) :: nDonorBlocks

         ! Names of the interfaces. Dimension [nDonorBlocks].

         character(len=maxCGNSNameLen), pointer, dimension(:) :: &
                                                            connectNames

         ! Names of the zone/block interfacing with the current
         ! zone/block. Dimension [nDonorBlocks].

         character(len=maxCGNSNameLen), pointer, dimension(:) :: &
                                                           donorNames

         ! Zone/block IDs of the zones/blocks interfacing with the
         ! current zone/block. Dimension [nDonorBlocks].

         integer(kind=intType), pointer, dimension(:) :: donorBlocks

         ! Range of points of this subface.

         integer(kind=intType) :: iBeg, jBeg, kBeg
         integer(kind=intType) :: iEnd, jEnd, kEnd

         ! Block face IDs of the donor blocks, which abut this subface.
         ! Dimension [nDonorBlocks].

         integer(kind=intType), pointer, dimension(:) :: donorFaceIDs

         ! Whether or not the subface is a periodic boundary.

         logical :: periodic

         ! The center of rotation for a periodic boundary.

         real(kind=realType), dimension(3) :: rotationCenter

         ! The rotation angles for a periodic boundary.

         real(kind=realType), dimension(3) :: rotationAngles

         ! The translation vector for a periodic boundary.

         real(kind=realType), dimension(3) :: translation

       end type cgnsNonMatchAbuttingConnType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type to store cgns overset  *
!      * connectivity (i.e. overlapping grids to be handled via the     *
!      * chimera approach).                                             *
!      *                                                                *
!      ******************************************************************
!
       type cgnsOversetConnType

         ! Name of the interface.

         character(len=maxCGNSNameLen) :: connectName

         ! Name of the zone/block interfacing with the current zone/block.

         character(len=maxCGNSNameLen) :: donorName

         ! Zone/block id of the zone/block interfacing with the current
         ! zone/block.

         integer(kind=intType) :: donorBlock

         ! Number of points to be interpolated (should equal the
         ! number of points in the donor list.

         integer(kind=intType) :: npnts

         ! Indices for this block to be interpolated.
         ! [dimension(3,npnts)]

         integer(kind=intType), pointer, dimension(:,:) :: ibndry

         ! Indices for donor block that provide information.
         ! [dimension(3,npnts)]

         integer(kind=intType), pointer, dimension(:,:) :: idonor

         ! Interpolation weights for the donor stencil
         ! [dimension(3,npnts)]

         real(kind=realType), pointer, dimension(:,:) :: interp

       end type cgnsOversetConnType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type to store cgns overset  *
!      * holes. (these points are grouped together and ignored when     *
!      * calculating residuals).                                        *
!      *                                                                *
!      ******************************************************************
!
       type cgnsHolesType

         ! Name of the interface.

         character(len=maxCGNSNameLen) :: holeName

         ! Number of points in this hole set.

         integer(kind=intType) :: npnts

         ! Indices for the hole points.
         ! [dimension(3,npnts)]

         integer(kind=intType), pointer, dimension(:,:) :: indices

       end type cgnsHolesType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type to store cgns block    *
!      * boundary conditions.                                           *
!      *                                                                *
!      ******************************************************************
!
       type cgnsBocoType

         ! Name of the boundary condition.

         character(len=maxCGNSNameLen) :: bocoName

         ! CGNS and internal boundary condition type.

         integer :: BCTypeCGNS
         integer(kind=intType) :: BCType

         ! Name of the CGNS user defined data node if the CGNS
         ! boundary condition is UserDefined.

         character(len=maxCGNSNameLen) :: userDefinedName

         ! The way the boundary condition faces are specified; either
         ! a point range or an individual set of points.

         integer :: ptSetType

         ! Number of points in the boundary condition set defining this
         ! boundary region. For a point range this is 2.

         integer(kind=intType) :: nPnts

         ! Index vector indicating the computational coordinate
         ! direction of the boundary condition patch normal.

         integer :: normalIndex

         ! A flag indicating whether or not boundary normals are defined.
         ! normalListFlag == 0: normals are not defined.
         ! normalListFlag == 1: normals are defined.

         integer :: normalListFlag

         ! Data type used for the definition of the normals. Admissible
         ! types are realSingle and realDouble.

         integer :: normalDataType

         ! Corresponding family number. If the face does not belong to
         ! a family this value is 0.

         integer(kind=intType) :: familyID

         ! The number of the sliding mesh interface of which this
         ! boco is part. 0 means that this family is not part of a
         ! sliding mesh interface. This value can be positive and
         ! negative in order to distinguish between the two sides of the
         ! interface. The absolute value is the actual ID of the
         ! interface.

         integer(kind=intType) :: slidingID

         ! Number of boundary condition datasets for the current
         ! boundary condition.

         integer(kind=intType) :: nDataSet

         ! The actual boundary condition data sets.

         type(cgnsBCDatasetType), pointer, dimension(:) :: dataSet

         ! Whether or not I actually allocated the memory for data_set.
         ! It is possible that data_set points to corresponding entry
         ! of a family.

         logical :: dataSetAllocated

         ! The rotation center and rotation rate of the boundary face.
         ! It is possible that this differs from the rotation rate of
         ! the corresponding block, e.g. for a casing in a
         ! turbomachinery problem.

         real(kind=realType), dimension(3) :: rotCenter, rotRate

         ! Range of points of this subface.

         integer(kind=intType) :: iBeg, jBeg, kBeg
         integer(kind=intType) :: iEnd, jEnd, kEnd

         ! Whether or not this subface is an actual face. Some mesh
         ! generators (such as ICEM CFD hexa) include edges and points
         ! as boundary conditions. These should not be considered by
         ! the flow solver. in those cases, actual_face is .false.

         logical :: actualFace

       end type cgnsBocoType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type to store the data of a *
!      * cgns block.                                                    *
!      *                                                                *
!      ******************************************************************
!
       type cgnsBlockInfoType
!
!        ****************************************************************
!        *                                                              *
!        * Information read from the cgns file.                         *
!        *                                                              *
!        ****************************************************************
!
         ! The type of the zone. Should be structured. Note that this
         ! is an integer and not integer(kind=intType).

         integer :: zoneType

         ! Zone name for this block.

         character(len=maxCGNSNameLen) :: zoneName

         ! The number or subblocks and the processor ID's on which they
         ! are stored. Due to the possibility of splitting the block
         ! during runtime, multiple processors could store a part of
         ! the block.

         integer                        :: nSubBlocks
         integer, dimension(:), pointer :: procStored

         ! The local block ID's of the subblocks.

         integer, dimension(:), pointer :: localBlockID

         ! The corresponding nodal ranges of the subblocks.

         integer, dimension(:), pointer :: iBegOr, jBegOr, kBegOr
         integer, dimension(:), pointer :: iEndOr, jEndOr, kEndOr

         ! The units in which the grid is specified.

         integer :: mass, len, time, temp, angle

         ! Whether or not grid units are specified.

         logical :: gridUnitsSpecified

         ! The conversion factor to meters for this block.

         real(kind=realType) :: LRef

         ! Corresponding family number. If the block does not belong to
         ! a family this value is 0.

         integer(kind=intType) :: familyID

         ! Nodal block dimensions.

         integer(kind=intType) :: il, jl, kl

         ! Cell block dimensions.

         integer(kind=intType) :: nx, ny, nz

         ! Total number of 1 to 1 block to block connectivities, i.e.
         ! continous grid lines, for this block. Also the number of
         ! 1 to 1 connectivities stored in general connectivity nodes
         ! is incorporated in n1to1.

         integer(kind=intType) :: n1to1

         ! Number of 1 to 1 block to block connectivities stored in
         ! general connectivities.

         integer(kind=intType) :: n1to1General

         ! Array of 1 to 1 block to block connectivities.

         type(cgns1to1ConnType), pointer, dimension(:) :: conn1to1

         ! Number of non-matching abutting block to block
         ! connectivities.

         integer(kind=intType) :: nNonMatchAbutting

         ! Array of non-matching abutting block to block connectivities.

         type(cgnsNonMatchAbuttingConnType), pointer, dimension(:) :: &
                                                   connNonMatchAbutting

         ! Number of overset block to block connectivities, and the
         ! total number of overset boundary cells.

         integer(kind=intType) :: nOverset, nCellsOverset

         ! Array of overset block to block connectivities.

         type(cgnsOversetConnType), pointer, dimension(:) :: connOver

         ! Number of overset hole sets for this block.

         integer(kind=intType) :: nHoles

         ! Array of overset hole sets.

         type(cgnsHolesType), pointer, dimension(:) :: hole

         ! Number of boundary conditions for this block.

         integer(kind=intType) :: nBocos

         ! Array of boundary conditions.

         type(cgnsBocoType), pointer, dimension(:) :: bocoInfo

         ! Whether or not a rotating frame is specified.

         logical :: rotatingFrameSpecified

         ! The corresponding rotation center and rotation rate.

         real(kind=realType), dimension(3) :: rotCenter, rotRate

       end type cgnsBlockInfoType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type to store the data of a *
!      * cgns family.                                                   *
!      *                                                                *
!      ******************************************************************
!
       type cgnsFamilyType

         ! Name of the family.

         character(len=maxCGNSNameLen) :: familyName

         ! Type of the boundary condition and family BC name.

         integer :: BCTypeCGNS
         integer(kind=intType) :: BCType

         character(len=maxCGNSNameLen) :: bcName

         ! Name of the CGNS user defined data node if the CGNS
         ! boundary condition is UserDefined.

         character(len=maxCGNSNameLen) :: userDefinedName

         ! The number of the sliding mesh interface of which this
         ! family is part. 0 means that this family is not part of a
         ! sliding mesh interface. This value can be positive and
         ! negative in order to distinguish between the two sides of the
         ! interface. The absolute value is the actual ID of the
         ! interface.

         integer(kind=intType) :: slidingID

         ! The number of the bleed flow region of which this family is
         ! part. 0 means that this family does not belong to a bleed
         ! flow region. There is no need to distinguish between an
         ! inflow and an outflow bleed, because they have different
         ! boundary conditions.

         integer(kind=intType) :: bleedRegionID

         ! Whether or not the mass flow must be monitored for this family.

         logical :: monitorMassFlow

         ! Number of boundary condition datasets for this family.

         integer(kind=intType) :: nDataSet

         ! The actual boundary condition data sets.

         type(cgnsBCDatasetType), pointer, dimension(:) :: dataSet

         ! Whether or not a rotating frame is specified.

         logical :: rotatingFrameSpecified

         ! The corresponding rotation center and rotation rate.

         real(kind=realType), dimension(3) :: rotCenter, rotRate

       end type cgnsFamilyType
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the variables stored in this module.             *
!      *                                                                *
!      ******************************************************************
!
       ! Dimensions of the cell and of the physical dimensions.
       ! Both should be 3 for this code. Note that these are integers
       ! and not integers(kind=intType).

       integer :: cgnsCellDim, cgnsPhysDim

       ! Number of blocks (zones) in the cgns grid.

       integer(kind=intType) :: cgnsNDom

       ! Array of cgns blocks.

       type(cgnsBlockInfoType), allocatable, dimension(:) :: cgnsDoms

       ! Number of families in the cgns grid.

       integer(kind=intType) :: cgnsNFamilies

       ! Array of families.

       type(cgnsFamilyType), allocatable, dimension(:) :: cgnsFamilies

       ! Number of sliding mesh interfaces in the grid.

       integer(kind=intType) :: cgnsNSliding

       ! The corresponding family ID's of the sliding interfaces.

       integer(kind=intType), allocatable, dimension(:,:) :: famIDsSliding

       ! Number of domain interfaces, i.e. interfaces with other CFD
       ! codes, in the grid.

       integer(kind=intType) :: cgnsNDomainInterfaces

       ! The family and BC ID's of the domain interfaces.

       integer(kind=intType), allocatable, dimension(:) :: &
                                            famIDsDomainInterfaces
       integer(kind=intType), allocatable, dimension(:,:) :: &
                                            bcIDsDomainInterfaces

       ! Name of the cgns base.

       character(len=maxCGNSNameLen) :: cgnsBaseName

       ! Whether or not there are overset grids present.

       logical :: oversetPresent

       ! massFlowFamilyInv(:,:):  Array to store the local contributions
       !                          from the central part of the flux to
       !                          the mass flow of a family and the
       !                          sliding mesh interfaces. Dimension is
       !                          (0:nn,nTimeIntervalsSpectral, where
       !                          nn is the number of families for which
       !                          the mass flow must be monitored plus
       !                          2*cgnsNSliding (if the mass flow through
       !                          the sliding interfaces must be monitored).
       !                          The reason for 2*cgnsNSliding is that each
       !                          side of a sliding interface is monitored.
       !                          The first index starts at 0 to store
       !                          all the faces that are not on a
       !                          sliding interface.
       ! massFlowFamilyDiss(:,:): Idem for the dissipative part.

       real(kind=realType), allocatable, dimension(:,:) :: massFlowFamilyInv
       real(kind=realType), allocatable, dimension(:,:) :: massFlowFamilyDiss

       end module cgnsGrid
