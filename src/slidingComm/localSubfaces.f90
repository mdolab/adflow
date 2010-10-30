!
!      ******************************************************************
!      *                                                                *
!      * File:          localSubfaces.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-17-2003                                      *
!      * Last modified: 04-07-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module localSubfacesMod
!
!      ******************************************************************
!      *                                                                *
!      * Local module to store the halo information for the cells as    *
!      * well as the subfaces and its quadrilateral faces whose values  *
!      * must be interpolated.                                          *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the derived data type to store the information   *
!      * of the quadrilateral faces for a subface on the interface.     *
!      *                                                                *
!      ******************************************************************
!
       integer(kind=porType), parameter :: normalQuad    = 0_porType
       integer(kind=porType), parameter :: piCrossed     = 1_porType
       integer(kind=porType), parameter :: polarSingular = 2_porType

       type localSubfaceType

         ! nQuad:    Number of quadrilateral faces on this subface.
         ! nNode:    Number of nodes on this subface.
         ! nDual:    Number of dual quadrilaterals.
         ! blockID:  Local block ID.
         ! faceID:   Block face on which the subface is located.

         integer(kind=intType) :: nQuad, nNode, nDual
         integer(kind=intType) :: blockID, faceID

         ! zeroCrossedQuad: Whether or not the subface crosses the line
         !                  theta == zero for the primary quadrilateral
         !                  surface grid.
         ! piCrossedQuad:   Idem for the line theta == pi.
         ! zeroCrossedDual: Idem for the dual quadrilateral surface
         !                  grid for theta == 0.
         ! piCrossedDual:   Idem for the line theta == pi.

         logical :: zeroCrossedQuad, piCrossedQuad
         logical :: zeroCrossedDual, piCrossedDual

         ! searchQuad(nQuad): Whether or not a donor must be searched
         !                    for the quadrilateral.
         ! searchNode(nNode): Idem for the nodes

         logical, dimension(:), pointer :: searchQuad, searchNode

         ! statusQuad(nQuad): What kind of quadrilateral face.
         !                    Possibilities are: - normalQuad
         !                                       - zeroCrossed
         !                                       - piCrossed
         !                                       - polarSingular
         !                    It is present to indicate whether or not
         !                    the face intersects the lines theta = 0
         !                    and theta = pi. Especially in the latter
         !                    case something special must be done.
         !                    polarSingular is a special case to
         !                    indicate that the cell center has a
         !                    smaller radius than all the nodes. Can
         !                    happen when the cell center coincides with
         !                    the polar singularity.
         ! statusDual(nQuad): Idem for the dual face.

         integer(kind=porType), dimension(:), pointer :: statusQuad
         integer(kind=porType), dimension(:), pointer :: statusDual

         ! indHalo1(nQuad,3): Indices of the corresponding 1st halo cell.
         ! indHalo2(nQuad,3): Indices of the corresponding 2nd halo cell.
         !                    A -1 for indHalo2(..,1) indicates that
         !                    this cell is not a direct halo.
         ! indHaloN(nNode,3): Indices of the corresponding halo node.

         integer(kind=intType), dimension(:,:), pointer :: indHalo1
         integer(kind=intType), dimension(:,:), pointer :: indHalo2
         integer(kind=intType), dimension(:,:), pointer :: indHaloN

         ! donorDualQ(nQuad): Index of the donor dual in the arrays to
         !                    store the surface grid if this quad is
         !                    to be searched.

         integer(kind=intType), dimension(:), pointer :: donorDualQ

         ! connQuad(nQuad,4): Local connectivity of the quadrilateral
         !                    grid faces.
         ! connDual(nDual,4): Idem for the dual faces.
         ! storeQuad(nQuad):  Whether or not the quadrilateral must
         !                    be stored in the entire surface mesh
         !                    of the sliding interface.
         ! storeDual(nDual):  Idem for the dual faces.

         integer(kind=intType), dimension(:,:), pointer :: connQuad
         integer(kind=intType), dimension(:,:), pointer :: connDual

         logical, dimension(:), pointer :: storeQuad
         logical, dimension(:), pointer :: storeDual

         ! nRotations(nQuad): Number of rotations needed needed to make
         !                    the quadrilateral face match a donor.
         ! donorProc(nQuad):  Processor ID in the communication group of
         !                    the sliding interface where the donor is
         !                    stored.

         integer(kind=intType), dimension(:), pointer :: nRotations
         integer(kind=intType), dimension(:), pointer :: donorProc

         ! coorQuad(nQuad,3): Cylindrical coordinates of the cell
         !                    centers.
         ! u(nQuad):          Interpolation weight in u-direction.
         ! v(nQuad):          Interpolation weight in v-direction.

         real(kind=realType), dimension(:,:), pointer :: coorQuad
         real(kind=realType), dimension(:),   pointer :: u, v

         ! coorN(nNode,3):    Cylindrical coordinates of the nodes.
         ! coorNInt(nNode,3): Idem, but now for the node one index
         !                    into the block.

         real(kind=realType), dimension(:,:), pointer :: coorN
         real(kind=realType), dimension(:,:), pointer :: coorNInt

       end type localSubfaceType
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the derived data type to store for the local     *
!      * blocks the type information of the cells, including the 1st    *
!      * halos. See the parameters below for the possibilities; another *
!      * possibility is a sliding mesh halo, but then the haloInfo gets *
!      * the number of the sliding interface.                           *
!      *                                                                *
!      ******************************************************************
!
       integer(kind=intType), parameter :: boundaryHalo = -3_intType
       integer(kind=intType), parameter :: unownedDonor = -2_intType
       integer(kind=intType), parameter :: ownedDonor   = -1_intType
       integer(kind=intType), parameter :: internalCell =  0_intType

       type blockHaloType

         ! haloInfo(ie,je,ke): The halo type for the cells.

         integer(kind=intType), dimension(:,:,:), pointer :: haloInfo

       end type blockHaloType
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the variables stored in this module.             *
!      *                                                                *
!      ******************************************************************
!
       ! nMySubfaces1:              Number of subfaces locally stored
       !                            for part 1 of the interface.
       ! mySubfaces1(nMySubfaces1): The corresponding interfaces.

       integer(kind=intType) :: nMySubfaces1
       type(localSubfaceType), dimension(:), allocatable :: mySubfaces1

       ! nMySubfaces2:              Number of subfaces locally stored
       !                            for part 2 of the interface.
       ! mySubfaces2(nMySubfaces2): The corresponding interfaces.

       integer(kind=intType) :: nMySubfaces2
       type(localSubfaceType), dimension(:), allocatable :: mySubfaces2

       ! donorDoms(nDom): Array of this type to store the true
       !                  donor info.

       type(blockHaloType), dimension(:), allocatable :: donorDoms

       end module localSubfacesMod
