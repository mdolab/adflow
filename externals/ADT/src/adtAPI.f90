!
!     ******************************************************************
!     *                                                                *
!     * File:          adtAPI.f90                                      *
!     * Author:        Edwin van der Weide                             *
!     * Starting date: 02-09-2006                                      *
!     * Last modified: 02-09-2006                                      *
!     *                                                                *
!     ******************************************************************
!
      module adtAPI
!
!     ******************************************************************
!     *                                                                *
!     * Module, which defines the API of the ADT routines. It is       *
!     * included in a module, such that an explicit interface is       *
!     * present.                                                       *
!     *                                                                *
!     ******************************************************************
!
      use adtBuild
      use adtSearch
      use adtUtils
      implicit none

      !=================================================================

      contains

        !===============================================================

        subroutine adtBuildSurfaceADT(nTria, nQuads,   nNodes,    &
                                      coor,  triaConn, quadsConn, &
                                      BBox,  useBBox,  comm,      &
                                      adtID)
!
!       ****************************************************************
!       *                                                              *
!       * This routine builds the 6 dimensional ADT, which stores the  *
!       * given surface grid. The memory intensive part of these       *
!       * arguments, the arrays with the coordinates and               *
!       * connectivities, are not copied. Instead pointers are set to  *
!       * these arrays. It is therefore the responsibility of the user *
!       * not to deallocate this memory before all the searches have   *
!       * been performed.                                              *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nNodes:    Number of local nodes in the given grid.          *
!       * nTria:     Idem for the triangles.                           *
!       * nQuads:    Idem for the quadrilaterals.                      *
!       * BBox(3,2): The possible bounding box. Only elements within   *
!       *            this box will be stored in the ADT.               *
!       * useBBox:   Whether or not to use the bounding box.           *
!       * comm:      MPI-communicator for the global ADT.              *
!       * adtID:     The ID of the ADT.                                *
!       *                                                              *
!       * Subroutine intent(in), target arguments.                     *
!       * ----------------------------------------                     *
!       * coor(3,nNodes):      Nodal coordinates of the local grid.    *
!       * triaConn(3,nTria):   Local connectivity of the triangles.    *
!       * quadsConn(4,nQuads): Idem for the quadrilaterals.            *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer, intent(in)          :: comm
        character(len=*), intent(in) :: adtID

        integer(kind=adtIntType), intent(in) :: nTria
        integer(kind=adtIntType), intent(in) :: nQuads
        integer(kind=adtIntType), intent(in) :: nNodes

        logical, intent(in) :: useBBox

        integer(kind=adtIntType), dimension(:,:), intent(in) :: triaConn
        integer(kind=adtIntType), dimension(:,:), intent(in) :: quadsConn

        real(kind=adtRealType), dimension(3,2), intent(in) :: BBox

        real(kind=adtRealType), dimension(:,:), intent(in) :: coor

        !===============================================================

        ! Call the subroutine buildSurfaceADT to do the actual work.

        call buildSurfaceADT(nTria,    nQuads,    nNodes, coor,    &
                             triaConn, quadsConn, BBox,   useBBox, &
                             comm,     adtID)

        end subroutine adtBuildSurfaceADT

        !***************************************************************
        !***************************************************************

        subroutine adtBuildVolumeADT(nTetra,    nPyra,    nPrisms,    &
                                     nHexa,     nNodes,   coor,       &
                                     tetraConn, pyraConn, prismsConn, &
                                     hexaConn,  BBox,     useBBox,    &
                                     comm,      adtID)
!
!       ****************************************************************
!       *                                                              *
!       * This routine builds the 6 dimensional ADT, which stores the  *
!       * given volume grid. The memory intensive part of these        *
!       * arguments, the arrays with the coordinates and               *
!       * connectivities, are not copied. Instead pointers are set to  *
!       * these arrays. It is therefore the responsibility of the user *
!       * not to deallocate this memory before all the searches have   *
!       * been performed.                                              *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nNodes:    Number of local nodes in the given grid.          *
!       * nTetra:    Idem for the tetrahedra.                          *
!       * nPyra:     Idem for the pyramids.                            *
!       * nPrisms:   Idem for the prisms.                              *
!       * nHexa:     Idem for the hexahedra.                           *
!       * BBox(3,2): The possible bounding box. Only elements within   *
!       *            this box will be stored in the ADT.               *
!       * useBBox:   Whether or not to use the bounding box.           *
!       * comm:      MPI-communicator for the global ADT.              *
!       * adtID:     The ID of the ADT.                                *
!       *                                                              *
!       * Subroutine intent(in), target arguments.                     *
!       * ----------------------------------------                     *
!       * coor(3,nNodes):        Nodal coordinates of the local grid.  *
!       * tetraConn(4,nTetra):   Local connectivity of the tetrahedra. *
!       * pyraConn(5,nPyra):     Idem for the pyramids.                *
!       * prismsConn(6,nPrisms): Idem for the prisms.                  *
!       * hexaConn(8,nHexa):     Idem for the hexahedra.               *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer, intent(in)          :: comm
        character(len=*), intent(in) :: adtID

        integer(kind=adtIntType), intent(in) :: nTetra
        integer(kind=adtIntType), intent(in) :: nPyra
        integer(kind=adtIntType), intent(in) :: nPrisms
        integer(kind=adtIntType), intent(in) :: nHexa
        integer(kind=adtIntType), intent(in) :: nNodes

        logical, intent(in) :: useBBox

        integer(kind=adtIntType), dimension(:,:), intent(in) :: tetraConn
        integer(kind=adtIntType), dimension(:,:), intent(in) :: pyraConn
        integer(kind=adtIntType), dimension(:,:), intent(in) :: prismsConn
        integer(kind=adtIntType), dimension(:,:), intent(in) :: hexaConn

        real(kind=adtRealType), dimension(3,2), intent(in) :: BBox

        real(kind=adtRealType), dimension(:,:), intent(in) :: coor

        !===============================================================

        ! Call the subroutine buildVolumeADT to do the actual work.

        call buildVolumeADT(nTetra,     nPyra,     nPrisms,   nHexa,    &
                            nNodes,     coor,      tetraConn, pyraConn, &
                            prismsConn, hexaConn,  BBox,      useBBox,  &
                            comm,       adtID)

        end subroutine adtBuildVolumeADT

        !***************************************************************
        !***************************************************************

        subroutine adtContainmentSearch(nCoor,  coor,        adtID,     &
                                        procID, elementType, elementID, &
                                        uvw,    nInterpol,   arrDonor,  &
                                        arrInterpol)
!
!       ****************************************************************
!       *                                                              *
!       * This routine attempts for every coordinate to find the       *
!       * element in the given ADT, which contains that coordinate.    *
!       * If no element is found the corresponding entry in procID is  *
!       * set to -1 to indicate failure.                               *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nCoor:     Number of coordinates for which the element must  *
!       *            be determined.                                    *
!       * coor:      The coordinates of these points.                  *
!       * adtID:     The ADT to be searched.                           *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * procID:      The ID of the processor in the group of the ADT *
!       *              where the element containing the point is       *
!       *              stored. If no element is found for a given      *
!       *              point the corresponding entry in procID is set  *
!       *              to -1 to indicate failure. Remember that the    *
!       *              processor ID's start at 0 and not at 1.         *
!       * elementType: The type of element which contains the point.   *
!       * elementID:   The entry in the connectivity of this element   *
!       *              which contains the point.                       *
!       * uvw:         The parametric coordinates of the point in the  *
!       *              transformed element; this transformation is     *
!       *              such that every element is transformed into a   *
!       *              standard element in parametric space. The u, v  *
!       *              and w coordinates can be used to determine the  *
!       *              actual interpolation weights.                   *
!       * arrInterpol: Array with the interpolated data.               *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=adtIntType), intent(in) :: nCoor, nInterpol
        character(len=*),         intent(in) :: adtID

        real(kind=adtRealType), dimension(:,:), intent(in) :: coor
        real(kind=adtRealType), dimension(:,:), intent(in) :: arrDonor

        integer,                  dimension(:), intent(out) :: procID
        integer(kind=adtIntType), dimension(:), intent(out) :: elementID

        integer(kind=adtElementType), dimension(:), intent(out) :: &
                                                            elementType
        real(kind=adtRealType), dimension(:,:), intent(out) :: uvw
        real(kind=adtRealType), dimension(:,:), intent(out) :: arrInterpol

        !===============================================================

        ! Call the subroutine containmentSearch to do the actual work.

        call containmentSearch(nCoor,       coor,      adtID, procID,    &
                               elementType, elementID, uvw,   nInterpol, &
                               arrDonor,    arrInterpol)

        end subroutine adtContainmentSearch

        !***************************************************************
        !***************************************************************

        subroutine adtDeallocateADTs(adtID)
!
!       ****************************************************************
!       *                                                              *
!       * This routine deallocates the memory for the given entry in   *
!       * the array ADTs and it tries to reallocate ADTs itself        *
!       * accordingly.                                                 *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * adtID: The entry in ADTs to be deallocated.                  *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        character(len=*), intent(in) :: adtID

        !===============================================================

        ! Call the subroutine deallocateADTs to do the actual work.

        call deallocateADTs(adtID)

        end subroutine adtDeallocateADTs

        !***************************************************************
        !***************************************************************

        subroutine adtFailSafeSearch(nCoor,    coor,        adtID,     &
                                     procID,   elementType, elementID, &
                                     uvw,      dist2,       nInterpol, &
                                     arrDonor, arrInterpol)
!
!       ****************************************************************
!       *                                                              *
!       * This routine attempts for every coordinate to find the       *
!       * element in the given ADT, which contains that coordinate.    *
!       * If no element is found a minimum distance search is          *
!       * performed, such that always an interpolation can be          *
!       * performed. To indicate that the element does not contain the *
!       * point the element ID is negated.                             *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nCoor: Number of coordinates for which the element must be   *
!       *        determined.                                           *
!       * coor:  The coordinates of these points.                      *
!       * adtID: The ADT to be searched.                               *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * ---------------------------------                            *
!       * procID:      The ID of the processor in the group of the ADT *
!       *              where the element containing the point is       *
!       *              stored. If no element is found for a given      *
!       *              point the corresponding entry in procID is set  *
!       *              to -1 to indicate failure. Remember that the    *
!       *              processor ID's start at 0 and not at 1.         *
!       * elementType: The type of element which contains the point.   *
!       * elementID:   The entry in the connectivity of this element   *
!       *              which contains the point. The ID is negative if *
!       *              the coordinate is outside the element, i.e. if  *
!       *              a minimum distance search had to be used.       *
!       * uvw:         The parametric coordinates of the point in the  *
!       *              transformed element; this transformation is     *
!       *              such that every element is transformed into a   *
!       *              standard element in parametric space. The u, v  *
!       *              and w coordinates can be used to determine the  *
!       *              actual interpolation weights.                   *
!       * arrInterpol: Array with the interpolated data.               *
!       *                                                              *
!       * Subroutine intent(inout) arguments.                          *
!       * -----------------------------------                          *
!       * dist2: Minimum distance squared of the coordinates to the    *
!       *        elements of the ADT. On input it should be            *
!       *        initialized by the calling program, possibly to a     *
!       *        large value. In this way it is possible to handle     *
!       *        periodic problems as efficiently as possible.         *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=adtIntType), intent(in) :: nCoor, nInterpol
        character(len=*),         intent(in) :: adtID

        real(kind=adtRealType), dimension(:,:), intent(in) :: coor
        real(kind=adtRealType), dimension(:,:), intent(in) :: arrDonor

        integer,                  dimension(:), intent(out) :: procID
        integer(kind=adtIntType), dimension(:), intent(out) :: elementID

        integer(kind=adtElementType), dimension(:), intent(out) :: &
                                                              elementType

        real(kind=adtRealType), dimension(:,:), intent(out) :: uvw
        real(kind=adtRealType), dimension(:,:), intent(out) :: arrInterpol

        real(kind=adtRealType), dimension(:), intent(inout) :: dist2

        !===============================================================

        ! Call the subroutine failSafeSearch to do the actual work.

        call failSafeSearch(nCoor,       coor,      adtID,    procID, &
                            elementType, elementID, uvw,      dist2,  &
                            nInterpol,   arrDonor, arrInterpol)

        end subroutine adtFailSafeSearch

        !***************************************************************
        !***************************************************************

        subroutine adtMinDistanceSearch(nCoor,    coor,        adtID,     &
                                        procID,   elementType, elementID, &
                                        uvw,      dist2,       nInterpol, &
                                        arrDonor, arrInterpol)
!
!       ****************************************************************
!       *                                                              *
!       * This routine attempts for every coordinate to find the       *
!       * element in the given ADT which minimizes the distance to     *
!       * this point.                                                  *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nCoor: Number of coordinates for which the element must be   *
!       *        determined.                                           *
!       * coor:  The coordinates of these points.                      *
!       * adtID: The ADT to be searched.                               *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * procID:      The ID of the processor in the group of the ADT *
!       *              where the element containing the point is       *
!       *              stored. If no element is found for a given      *
!       *              point the corresponding entry in procID is set  *
!       *              to -1 to indicate failure. Remember that the    *
!       *              processor ID's start at 0 and not at 1.         *
!       * elementType: The type of element which contains the point.   *
!       * elementID:   The entry in the connectivity of this element   *
!       *              which contains the point. The ID is negative if *
!       *              the coordinate is outside the element.          *
!       * uvw:         The parametric coordinates of the point in the  *
!       *              transformed element; this transformation is     *
!       *              such that every element is transformed into a   *
!       *              standard element in parametric space. The u, v  *
!       *              and w coordinates can be used to determine the  *
!       *              actual interpolation weights. If the tree       *
!       *              corresponds to a surface mesh the third entry   *
!       *              of this array will not be filled.               *
!       * arrInterpol: Array with the interpolated data.               *
!       *                                                              *
!       * Subroutine intent(inout) arguments.                          *
!       * -----------------------------------                          *
!       * dist2: Minimum distance squared of the coordinates to the    *
!       *        elements of the ADT. On input it should be            *
!       *        initialized by the calling program, possibly to a     *
!       *        large value. In this way it is possible to handle     *
!       *        periodic problems as efficiently as possible.         *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=adtIntType), intent(in) :: nCoor, nInterpol
        character(len=*),         intent(in) :: adtID

        real(kind=adtRealType), dimension(:,:), intent(in) :: coor
        real(kind=adtRealType), dimension(:,:), intent(in) :: arrDonor

        integer, dimension(:), intent(out) :: procID
        integer(kind=adtIntType), dimension(:), intent(out) :: elementID

        integer(kind=adtElementType), dimension(:), intent(out) :: &
                                                            elementType

        real(kind=adtRealType), dimension(:,:), intent(out) :: uvw
        real(kind=adtRealType), dimension(:,:), intent(out) :: arrInterpol

        real(kind=adtRealType), dimension(:), intent(inout) :: dist2

        !===============================================================

        ! Call the subroutine minDistanceSearch to do the actual work.

        call minDistanceSearch(nCoor,       coor,      adtID,    procID, &
                               elementType, elementID, uvw,      dist2,  &
                               nInterpol,   arrDonor,  arrInterpol)

        end subroutine adtMinDistanceSearch

      end module adtAPI
