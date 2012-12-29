!
!     ******************************************************************
!     *                                                                *
!     * File:          adtData.f90                                     *
!     * Author:        Edwin van der Weide                             *
!     * Starting date: 11-27-2004                                      *
!     * Last modified: 02-17-2005                                      *
!     *                                                                *
!     ******************************************************************
!
      module adtData
!
!     ******************************************************************
!     *                                                                *
!     * Module, which defines the derived data types and the arrays to *
!     * store multiple ADT's. An array is chosen to store multiple     *
!     * ATD's rather than a linked list, because this is more          *
!     * convenient during the search. When the ADT's are built there   *
!     * is some additional work due to reallocation. However this is   *
!     * negligible due to the usage of pointers.                       *
!     *                                                                *
!     ******************************************************************
!
      use precision
#ifdef USE_COMPLEX
      use complexify
#endif
      implicit none
      save
!
!     ******************************************************************
!     *                                                                *
!     * Define the functions needed for the sorting of the derived     *
!     * data types to be private, i.e. they can only be accessed       *
!     * within this module.                                            *
!     *                                                                *
!     ******************************************************************
!
      public
      private :: adtBBoxTargetTypeLessEqual
      private :: adtBBoxTargetTypeLess
      private :: adtTypeAssign
!
!     ******************************************************************
!     *                                                                *
!     * Definition of some constants.                                  *
!     *                                                                *
!     ******************************************************************
!
      real(kind=realType), parameter :: adtZero   = 0.0_realType
      real(kind=realType), parameter :: adtFourth = 0.25_realType
      real(kind=realType), parameter :: adtHalf   = 0.5_realType
      real(kind=realType), parameter :: adtOne    = 1.0_realType
!
!     ******************************************************************
!     *                                                                *
!     * Definition of the parameters, which define whether the ADT     *
!     * corresponds to surface or volume elements and the parameter,   *
!     * which defines the maximum number of coordinates an ADT can     *
!     * handle in one interpolation round.                             *
!     *                                                                *
!     ******************************************************************
!
      integer, parameter :: adtSurfaceADT = 1
      integer, parameter :: adtVolumeADT  = 2

      integer(kind=intType), parameter :: nCoorMaxLowerLimit = 100000
!
!     ******************************************************************
!     *                                                                *
!     * Definition of the parameters, which define the supported       *
!     * element types. To save memory these parameters and the arrays  *
!     * containing the data are of a different integer type.           *
!     *                                                                *
!     ******************************************************************
!
      integer(kind=adtElementType), parameter :: adtTriangle      = 1
      integer(kind=adtElementType), parameter :: adtQuadrilateral = 2
      integer(kind=adtElementType), parameter :: adtTetrahedron   = 3
      integer(kind=adtElementType), parameter :: adtPyramid       = 4
      integer(kind=adtElementType), parameter :: adtPrism         = 5
      integer(kind=adtElementType), parameter :: adtHexahedron    = 6
!
!     ******************************************************************
!     *                                                                *
!     * Definition of the derived data type store a leaf of an ADT.    *
!     * The ADT itself is an array of these leaves.                    *
!     *                                                                *
!     ******************************************************************
!
      type adtLeafType

        ! children(2): Children of the parent. If negative it means that
        !              it is a terminal leaf and the absolute values
        !              indicate the bounding box id's. Note that it is
        !              allowed that 1 child is negative and the other
        !              positive.
        ! xMin(6):     The minimum coordinates of the leaf.
        ! xMax(6):     The maximum coordinates of the leaf.

        integer(kind=intType), dimension(2) :: children
#ifdef USE_COMPLEX
        complex(kind=realType), dimension(6)  ::   xMin, xMax
#else
        real(kind=realType),   dimension(6) :: xMin, xMax
#endif
      end type adtLeafType
!
!     ******************************************************************
!     *                                                                *
!     * The definition of adtBBoxTargetType, which stores the data of  *
!     * a possible bounding box which minimizes the distances to the   *
!     * given coordinate.                                              *
!     *                                                                *
!     ******************************************************************
!
      type adtBBoxTargetType

        ! ID:       The id of the bounding box in the list.
        ! posDist2: the possible minimum distance squared to the active
        !           coordinate.

        integer(kind=intType) :: ID
#ifdef USE_COMPLEX
        complex(kind=realType)   :: posDist2
#else
        real(kind=realType)   :: posDist2
#endif
      end type adtBBoxTargetType

      ! Interfaces for the extension of the operators <= and <.
      ! These are needed for the sorting of BBoxTargetType. Note
      ! that the = operator does not need to be defined, because
      ! BBoxTargetType only contains primitive types.

      interface operator(<=)
        module procedure adtBBoxTargetTypeLessEqual
      end interface

      interface operator(<)
        module procedure adtBBoxTargetTypeLess
      end interface
!
!     ******************************************************************
!     *                                                                *
!     * Definition of the derived data type to store an ADT.           *
!     *                                                                *
!     ******************************************************************
!
      type adtType

        ! comm:   The communicator of this ADT.
        ! nProcs: The number of processors which participate in this
        !         ADT.
        ! myID:   My processor ID in the processor group of comm.

        integer :: comm, nProcs, myID

        ! adtType:  Type of ADT. Possibilities are adtSurfaceADT and
        !           adtVolumeADT.
        ! adtID:    The given ID of the ADT.
        ! isActive: Whether or not the ADT is active. If not, this
        !           entry could be used during a reallocation.

        integer           :: adtType
        character(len=64) :: adtID
        logical           :: isActive

        ! nNodes:  Number of local nodes in the given grid.
        ! nTria:   Number of local triangles in the given grid.
        ! nQuads:  Idem for the quadrilaterals.
        ! nTetra:  Idem for the tetrahedra.
        ! nPyra:   Idem for the pyramids.
        ! nPrisms: Idem for the prisms.
        ! nHexa:   Idem for the hexahedra.

        integer(kind=intType) :: nNodes, nTria, nQuads
        integer(kind=intType) :: nTetra, nPyra, nPrisms, nHexa

        ! coor(3,nNodes): Nodal coordinates of the local grid.
        !                 To save memory this pointer is not
        !                 allocated, but set to the data given.
#ifdef USE_COMPLEX
        complex(kind=realType), dimension(:,:), pointer :: coor
#else
        real(kind=realType), dimension(:,:), pointer :: coor
#endif
        ! triaConn(3,nTria):     Local connectivity of the triangles.
        !                        To save memory this pointer is not
        !                        allocated, but set to the data given.
        ! quadsConn(4,nQuads):   Idem for the quadrilaterals.
        ! tetraConn(4,nTetra):   Idem for the tetrahedra.
        ! pyraConn(5,nPyra):     Idem for the pyramids.
        ! prismsConn(6,nPrisms): Idem for the prisms.
        ! hexaConn(8,nHexa):     Idem for the hexahedra.

        integer(kind=intType), dimension(:,:), pointer :: triaConn
        integer(kind=intType), dimension(:,:), pointer :: quadsConn
        integer(kind=intType), dimension(:,:), pointer :: tetraConn
        integer(kind=intType), dimension(:,:), pointer :: pyraConn
        integer(kind=intType), dimension(:,:), pointer :: prismsConn
        integer(kind=intType), dimension(:,:), pointer :: hexaConn

        ! nRootLeaves:        Number of non-empty root leaves.
        !                     This number is of course less than or
        !                     equal to nProcs.
        ! myEntryInRootProcs: If the local tree is not empty, this
        !                     contains the entry in rootLeavesProcs.
        ! rootLeavesProcs(:): The corresponding processor ID's.
        ! rootBBoxes(3,2,:):  The 3D bounding boxes of the non-empty
        !                     root leaves.

        integer :: nRootLeaves, myEntryInRootProcs
        integer, dimension(:), pointer :: rootLeavesProcs
#ifdef USE_COMPLEX
        complex(kind=realType), dimension(:,:,:), pointer :: rootBBoxes
#else
        real(kind=realType), dimension(:,:,:), pointer :: rootBBoxes
#endif

        ! nBBoxes:              Number of bounding boxes stored in
        !                       the ADT.
        ! elementType(nBBoxes): The corresponding element type of the
        !                       bounding box.
        ! elementID(nBBoxes):   The corresponding entry in the element
        !                       connectivity of the bounding box.
        ! xBBox(6,nBBoxes):     The coordinates of the bounding boxes
        !                       of the elements stored in the ADT.

        integer(kind=intType) :: nBBoxes

        integer(kind=adtElementType), dimension(:), pointer :: elementType
        integer(kind=intType),     dimension(:), pointer :: elementID
#ifdef USE_COMPLEX
        complex(kind=realType),     dimension(:,:), pointer :: xBBox
#else
        real(kind=realType),     dimension(:,:), pointer :: xBBox
#endif

        ! nLeaves:         Number of present in the ADT. Due to the
        !                  variable splitting the tree is optimally
        !                  balanced and therefore nLeaves = nBBoxes -1.
        ! ADTree(nLeaves): The alternating digital tree.

        integer(kind=intType) :: nLeaves
        type(adtLeafType), dimension(:), pointer :: ADTree

      end type adtType

      ! Interface for the extension of the operator =.

      interface assignment(=)
        module procedure adtTypeAssign
      end interface
!
!     ******************************************************************
!     *                                                                *
!     *           Variables stored in this module.                     *
!     *                                                                *
!     ******************************************************************
!
      ! ADTs(:): The array to store the different ADT's.

      type(adtType), dimension(:), allocatable :: ADTs

      ! nProcRecv:      Number of processors from which I receive
      !                 coordinates that must be searched in my ADT.
      ! nCoorMax:       Maximum number of coordinates that can be
      !                 searched during an interpolation round.
      ! nRounds:        Number of rounds in the outer loop of the search
      !                 algorithm.
      ! nLocalInterpol: Number of local coordinates that must be
      !                 searched in the locally stored tree.

      integer :: nProcRecv

      integer(kind=intType) :: nCoorMax
      integer(kind=intType) :: nRounds
      integer(kind=intType) :: nLocalInterpol


      ! procRecv:         Processor ID's from which I will receive
      !                   coordinates.
      ! nCoorProcRecv:    Number of coordinates I must receive from the
      !                   processors which send coordinates to me.
      ! nCoorPerRootLeaf: Number of coordinates, which may be searched
      !                   in each of the local ADT's. The array is in
      !                   cumulative storage format.
      ! mCoorPerRootLeaf: Idem, but its contents changes during the
      !                   iterative algorithm.
      ! coorPerRootLeaf:  The ID's of the corresponding coordinates.

      integer, dimension(:), allocatable :: procRecv

      integer(kind=intType), dimension(:), allocatable :: nCoorProcRecv
      integer(kind=intType), dimension(:), allocatable :: nCoorPerRootLeaf
      integer(kind=intType), dimension(:), allocatable :: mCoorPerRootLeaf
      integer(kind=intType), dimension(:), allocatable :: coorPerRootLeaf

      !=================================================================

      contains

        !===============================================================

        logical function adtBBoxTargetTypeLessEqual(g1,g2)
!
!       ****************************************************************
!       *                                                              *
!       * This function returns .true. if g1 <= g2. The comparison is  *
!       * firstly based on the possible minimum distance such that the *
!       * most likely candidates are treated first. In case of ties    *
!       * the boundary box ID is considered.                           *
!       *                                                              *
!       * Function intent(in) arguments.                               *
!       * ------------------------------                               *
!       * g1, g2: The two instances of the derived datatype that most  *
!       *         be compared.                                         *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Function arguments.
!
        type(adtBBoxTargetType), intent(in) :: g1, g2
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Compare the possible minimum distances.

        if(g1%posDist2 < g2%posDist2) then
          adtBBoxTargetTypeLessEqual = .true.
          return
        else if(g1%posDist2 > g2%posDist2) then
          adtBBoxTargetTypeLessEqual = .false.
          return
        endif

        ! Compare the bounding box ID's.

        if(g1%ID < g2%ID) then
          adtBBoxTargetTypeLessEqual = .true.
          return
        else if(g1%ID > g2%ID) then
          adtBBoxTargetTypeLessEqual = .false.
          return
        endif

        ! g1 and g2 are identical. Return .true.

        adtBBoxTargetTypeLessEqual = .true.

        end function adtBBoxTargetTypeLessEqual

        !===============================================================

        logical function adtBBoxTargetTypeLess(g1,g2)
!
!       ****************************************************************
!       *                                                              *
!       * This function returns .true. if g1 < g2. The comparison is   *
!       * firstly based on the possible minimum distance such that the *
!       * most likely candidates are treated first. In case of ties    *
!       * the boundary box ID is considered.                           *
!       *                                                              *
!       * Function intent(in) arguments.                               *
!       * ------------------------------                               *
!       * g1, g2: The two instances of the derived datatype that most  *
!       *         be compared.                                         *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Function arguments.
!
        type(adtBBoxTargetType), intent(in) :: g1, g2
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Compare the possible minimum distances.

        if(g1%posDist2 < g2%posDist2) then
          adtBBoxTargetTypeLess = .true.
          return
        else if(g1%posDist2 > g2%posDist2) then
          adtBBoxTargetTypeLess = .false.
          return
        endif

        ! Compare the bounding box ID's.

        if(g1%ID < g2%ID) then
          adtBBoxTargetTypeLess = .true.
          return
        else if(g1%ID > g2%ID) then
          adtBBoxTargetTypeLess = .false.
          return
        endif

        ! g1 and g2 are identical. Return .false.

        adtBBoxTargetTypeLess = .false.

        end function adtBBoxTargetTypeLess

        !===============================================================

        subroutine adtTypeAssign(g1, g2)
!
!       ****************************************************************
!       *                                                              *
!       * This subroutine defines the generic assignment operator for  *
!       * the derived datatype adtType. The contents of g1 is copied   *
!       * into g2, where pointers just point to the other pointers,    *
!       * i.e. no additional allocation takes place.                   *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * g2: Entity whose data should be copied.                      *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * g1: Entity whose data should be assigned.                    *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        type(adtType), intent(in)  :: g2
        type(adtType), intent(out) :: g1
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        g1%comm     = g2%comm
        g1%nProcs   = g2%nProcs
        g1%myID     = g2%myID
        g1%adtType  = g2%adtType
        g1%adtID    = g2%adtID
        g1%isActive = g2%isActive

        g1%nNodes  = g2%nNodes
        g1%nTria   = g2%nTria
        g1%nQuads  = g2%nQuads
        g1%nTetra  = g2%nTetra
        g1%nPyra   = g2%nPyra
        g1%nPrisms = g2%nPrisms
        g1%nHexa   = g2%nHexa

        g1%coor       => g2%coor
        g1%triaConn   => g2%triaConn
        g1%quadsConn  => g2%quadsConn
        g1%tetraConn  => g2%tetraConn
        g1%pyraConn   => g2%pyraConn
        g1%prismsConn => g2%prismsConn
        g1%hexaConn   => g2%hexaConn

        g1%nRootLeaves        = g2%nRootLeaves
        g1%myEntryInRootProcs = g2%myEntryInRootProcs
        g1%rootLeavesProcs   => g2%rootLeavesProcs
        g1%rootBBoxes        => g2%rootBBoxes

        g1%nBBoxes     =  g2%nBBoxes
        g1%elementType => g2%elementType
        g1%elementID   => g2%elementID
        g1%xBBox       => g2%xBBox

        g1%nLeaves = g2%nLeaves
        g1%ADTree => g2%ADTree

        end subroutine adtTypeAssign

      end module adtData
