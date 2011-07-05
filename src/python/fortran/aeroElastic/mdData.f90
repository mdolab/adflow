!
!      ******************************************************************
!      *                                                                *
!      * File:          mdData.f90                                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-16-2004                                      *
!      * Last modified: 03-24-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module mdData
!
!      ******************************************************************
!      *                                                                *
!      * This module contains arrays to store data needed for           *
!      * multi-disciplinary problems. This data is meant to be accessed *
!      * by a controlling python script, which transfers it to a        *
!      * different discipline.                                          *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
       save

       ! mdNFamilies:   The number of families present in the grid.
       ! mdFamilyNames: The cgns names of these families.

       integer(kind=intType) :: mdNFamilies
       character(len=maxCGNSNameLen), dimension(:), allocatable :: &
                                                         mdFamilyNames

       ! mdNSurfNodes(0:nProc,mdNFamilies): Number of surface nodes
       !                                    per processor per family.
       !                                    Stored in cumulative
       !                                    storage format.

       integer(kind=intType), dimension(:,:), allocatable :: mdNSurfNodes

       ! mdNSurfPatches(0:nProc,mdNFamilies): Number of surface patches
       !                                      per processor per family
       !                                      in cumulative storage
       !                                      format.
       ! mdPatchDimensions(3,nPatches):       The nodal dimensions of
       !                                      these patches.

       integer(kind=intType), dimension(:,:), allocatable :: &
                                                        mdNSurfPatches
       integer(kind=intType), dimension(:,:), allocatable :: &
                                                     mdPatchDimensions

       ! mdSurfInd(4,mdNSurfNodes):   Surface indices of all nodes on the
       !                              surface; index 4 is the block ID.
       ! mdSurfxx(3,mdNSurfNodes):    Surface coordinates of all
       !                              these nodes.
       ! mdSurfForce(3,mdNSurfNodes): Aerodynamic forces in these nodes.
       ! mdSurfVal(mdNSurfNodes):     Storage for the value of a
       !                              scalar variable at the surface.

       integer(kind=intType), dimension(:,:), allocatable :: mdSurfInd

       real(kind=realType), dimension(:,:), allocatable :: mdSurfxx
       real(kind=realType), dimension(:,:), allocatable :: mdSurfForce
       real(kind=realType), dimension(:),   allocatable :: mdSurfVal

       end module mdData
