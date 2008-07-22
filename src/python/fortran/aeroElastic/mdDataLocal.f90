!
!      ******************************************************************
!      *                                                                *
!      * File:          mdDataLocal.f90                                 *
!      * Author:        C.A.(Sandy) Mader, Edwin van der Weide          *
!      * Starting date: 10-24-2007                                      *
!      * Last modified: 10-24-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       module mdDataLocal
!
!      ******************************************************************
!      *                                                                *
!      * This module contains arrays to store data needed for           *
!      * multi-disciplinary problems. This data is meant to be accessed *
!      * by a controlling python script, which transfers it to a        *
!      * different discipline. In this module, tha data is stored only  *
!      * for the local processor                                        *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
       save

       ! mdNFamilies:   The number of families present in the grid.
       ! mdFamilyNames: The cgns names of these families.

       integer(kind=intType) :: mdNFamiliesLocal
       character(len=maxCGNSNameLen), dimension(:), allocatable :: &
                                                         mdFamilyNamesLocal

       ! mdNSurfNodes(0:nProc,mdNFamilies): Number of surface nodes
       !                                    per processor per family.
       !                                    Stored in cumulative
       !                                    storage format.

       integer(kind=intType), dimension(:), allocatable :: mdNSurfNodesLocal

       ! mdNSurfPatches(0:nProc,mdNFamilies): Number of surface patches
       !                                      per processor per family
       !                                      in cumulative storage
       !                                      format.
       ! mdPatchDimensions(3,nPatches):       The nodal dimensions of
       !                                      these patches.

       integer(kind=intType), dimension(:,:), allocatable :: &
                                                        mdNSurfPatchesLocal
       integer(kind=intType), dimension(:,:), allocatable :: &
                                                     mdPatchDimensionsLocal

       ! mdSurfInd(4,mdNSurfNodes):   Surface indices of all nodes on the
       !                              surface; index 4 is the block ID.
       ! mdSurfxx(3,mdNSurfNodes):    Surface coordinates of all
       !                              these nodes.
       ! mdSurfForce(3,mdNSurfNodes): Aerodynamic forces in these nodes.
       ! mdSurfVal(mdNSurfNodes):     Storage for the value of a
       !                              scalar variable at the surface.

       integer(kind=intType), dimension(:,:), allocatable :: mdSurfIndLocal

       real(kind=realType), dimension(:,:), allocatable :: mdSurfxxLocal
       real(kind=realType), dimension(:,:), allocatable :: mdSurfForceLocal
       real(kind=realType), dimension(:),   allocatable :: mdSurfValLocal

     end module mdDataLocal
