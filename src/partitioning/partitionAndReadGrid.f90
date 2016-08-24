!
!      ******************************************************************
!      *                                                                *
!      * File:          partitionAndReadGrid.f90                        *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 12-17-2002                                      *
!      * Last modified: 11-30-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine partitionAndReadGrid(partitionOnly)
!
!      ******************************************************************
!      *                                                                *
!      * partitionAndReadGrid determines the partitioning of the        *
!      * multiblock grid over the processors and reads the grid of the  *
!      * blocks (or block parts) assigned to this processor. Other      *
!      * preprocessing activities, such as the proper setup of the halo *
!      * communication structure, creation of coarse grids and wall     *
!      * distance computation, are performed in the preprocessing       *
!      * library.                                                       *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use blockPointers
       use cgnsGrid
       use inputIO
       use IOModule
       use partitionMod
       use utils, only : terminate
       implicit none

       logical, intent(in) :: partitionOnly
!
!      Local variables
!
       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of grid files that must be read,
       ! as well as the corresponding file names.

       call determineGridFileNames

       ! Read the number of blocks and the block sizes of the grid stored
       ! in the cgns grid file. This info is stored on all processors

       call readBlockSizes
       call determineNeighborIDs
       
       ! Some extra work is required to determine the IDs for the 
       ! sliding mesh and domain interfaces, which are specified via
       ! UserDefined boundary conditions in CGNS. The IDs for other
       ! code interfaces (domain interfaces) are stored in the module
       ! cgnsGrid. For sliding meshes, this routine pairs the two
       ! sides of the slide.

       call determineInterfaceIDs

       ! If we are just doing a partition test, return
       if (partitionOnly) then 
          return
       end if

       ! Determine the number of blocks to be stored on this processor
       ! and the relation to the original grid. Remember that blocks
       ! can be split for load balancing reasons.

       call loadBalance

      ! Initialize the iblank array for the fine grid domains

       call initFineGridIblank

       ! Allocate the coordinates of the fine grid level and the
       ! derived data type used to read them.

       call allocCoorFineGrid

       ! Read the grid of the blocks (block parts) to be stored
       ! on this processor. 

       call readGrid

       ! Determine the number of colors of the sliding mesh interfaces,
       ! such that the computation of the communication pattern of the
       ! sliding meshes is as load balanced as possible.

       call slidingMeshGroups

       ! Determine for the time spectral mode the time of one period,
       ! the rotation matrices for the velocity components and
       ! create the fine grid coordinates of all time spectral locations.

       call timePeriodSpectral
       call timeRotMatricesSpectral
       call fineGridSpectralCoor

       ! Release the memory of fileIDs, gridFiles and IOVar.
       ! They are not needed anymore.

       deallocate(fileIDs, gridFiles, IOVar, stat=ierr)
       if(ierr /= 0)                            &
            call terminate("partitionAndReadGrid", &
                        "Deallocation failure for fileIDs, gridFiles &
                        &and IOVar")

       ! Check if for all faces on the block boundaries either a
       ! physical boundary condition or a connectivity has been
       ! specified and check if the 1 to 1 subfaces match.

       call checkFaces
       call check1to1Subfaces

       end subroutine partitionAndReadGrid
