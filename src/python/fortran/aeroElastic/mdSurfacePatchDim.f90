!
!      ******************************************************************
!      *                                                                *
!      * File:          mdSurfacePatchDim.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-26-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mdSurfacePatchDim(famID,startInd,endInd)
!
!      ******************************************************************
!      *                                                                *
!      * mdSurfacePatchDim determines the dimensions of the surface     *
!      * patches for the given family ID.                               *
!      * On return the variables startInd and endInd are set to the     *
!      * range this family occupies in the array of all families.       *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use communication
       use mdData
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)  :: famID
       integer(kind=intType), intent(out) :: startInd, endInd
!
!      Local variables
!
       integer :: ierr, size
       integer, dimension(nProc) :: recvcounts, displs

       integer(kind=intType) :: nn, mm, ii, jj
       integer(kind=intType) :: modFamID, nSurfPatchesLoc

       integer(kind=intType), dimension(:,:), allocatable :: dimLoc

       logical :: storeSubface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Perform a check to see if this routine is called correctly.
       ! If not, terminate the program.

       if(famID == 0 .and. cgnsNfamilies > 0) then
         if(myID == 0) &
           call terminate("mdSurfacePatchDim", &
                          "Family ID 0 is only allowed when no family &
                          &info is present in the grid")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Check if the number of patches per family is already stored.
       ! If not, determine it.

       if(.not. allocated(mdNSurfPatches) ) call mdCreateNPatch

       ! Allocate the memory for the local surface coordinates.
       ! ModFamID is introduced to take famID == 0 into account.

       modFamID = max(famID, 1_intType)
       nSurfPatchesLoc = mdNSurfPatches(myID+1,modFamID) &
                       - mdNSurfPatches(myID,  modFamID)

       allocate(dimLoc(3,nSurfPatchesLoc), stat=ierr)
       if(ierr /= 0)                         &
         call terminate("mdSurfacePatchDim", &
                        "Memory allocation failure for dimLoc")

       ! Store the coordinates of the local patch dimensions.

       ii = 0
       domains: do nn=1,nDom

         ! Have the pointers in blockPointers point to this block
         ! on the finest mg level. There is no need to distinguish
         ! between the different spectral solutions, because this
         ! info is the same for all of them.

         call setPointers(nn,1_intType,1_intType)

         ! Loop over the number of boundary subfaces of this block.

         bocos: do mm=1,nBocos

           ! Check if the data of this subface must be stored.

           storeSubface = .false.

           if(famID == 0) then

             ! No family info present; all solid wall points are stored.

             if(BCType(mm) == EulerWall       .or. &
                BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) storeSubface = .true.

           else

             ! Family info is present. Check if this subface belongs
             ! to the given familyID.

             jj = cgnsSubface(mm)
             jj = cgnsDoms(nbkGlobal)%bocoInfo(jj)%familyID
             if(jj == famID) storeSubface = .true.

           endif

           ! Store the number of nodes in the 3 coordinate directions
           ! of the subface if needed.

           if( storeSubface ) then
             ii = ii + 1

             dimLoc(1,ii) = inEnd(mm) - inBeg(mm) + 1
             dimLoc(2,ii) = jnEnd(mm) - jnBeg(mm) + 1
             dimLoc(3,ii) = knEnd(mm) - knBeg(mm) + 1
           endif

         enddo bocos
       enddo domains

       ! Test if the memory of mdPatchDimensions has already been
       ! allocated. If not, allocate it.

       if(.not. allocated(mdPatchDimensions) ) then
         jj = mdNSurfPatches(nProc,max(cgnsNfamilies,1_intType))
         allocate(mdPatchDimensions(3,jj), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("mdSurfacePatchDim", &
                          "Memory allocation failure for &
                          &mdPatchDimensions")
       endif

       ! Construct the arrays recvcounts and displs needed for
       ! allgatherv.

       do nn=1,nProc
         recvcounts(nn) = 3*(mdNSurfPatches(nn,  modFamID) &
                        -    mdNSurfPatches(nn-1,modFamID))
         displs(nn)     = 3* mdNSurfPatches(nn-1,modFamID)
       enddo

       ! Call allgatherv to gather the data.

       size = 3*nSurfPatchesLoc
       call mpi_allgatherv(dimLoc, size, sumb_integer,            &
                           mdPatchDimensions, recvcounts, displs, &
                           sumb_integer, SUmb_comm_world, ierr)

       ! Release the memory of dimLoc.

       deallocate(dimLoc, stat=ierr)
       if(ierr /= 0)                         &
         call terminate("mdSurfacePatchDim", &
                        "Deallocation failure for dimLoc")

       ! Set the values of startInd and endInd to the range of the
       ! given family ID.

       startInd = mdNSurfPatches(0,    modFamID) + 1
       endInd   = mdNSurfPatches(nProc,modFamID)

       end subroutine mdSurfacePatchDim

!      ==================================================================

       subroutine mdDeleteSurfacePatchDim
!
!      ******************************************************************
!      *                                                                *
!      * mdDeleteSurfacePatchDim deallocates the memory of              *
!      * mdNSurfPatches and mdPatchDimensions, if allocated.            *
!      *                                                                *
!      ******************************************************************
!
       use mdData
       implicit none
!
!      Local variables.
!
       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       if( allocated(mdNSurfPatches) ) then

         deallocate(mdNSurfPatches, stat=ierr)
         if(ierr /= 0)                               &
           call terminate("mdDeleteSurfacePatchDim", &
                          "Deallocation error for mdNSurfPatches")
       endif

       if( allocated(mdPatchDimensions) ) then

         deallocate(mdPatchDimensions, stat=ierr)
         if(ierr /= 0)                               &
           call terminate("mdDeleteSurfacePatchDim", &
                          "Deallocation error for mdPatchDimensions")
       endif

       end subroutine mdDeleteSurfacePatchDim

!      ==================================================================

       subroutine mdCreateNPatch
!
!      ******************************************************************
!      *                                                                *
!      * mdCreateNPatch determines the total number of surface          *
!      * patches per processor per family.                              *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use communication
       use mdData
       implicit none
!
!      Local variables.
!
       integer :: ierr, size

       integer(kind=intType) :: jj, nn, mm

       integer(kind=intType), &
          dimension(max(cgnsNfamilies,1_intType)) :: nPatchesLoc
       integer(kind=intType), &
          dimension(max(cgnsNfamilies,1_intType),nProc) :: nPatches
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the local number of patches per family.

       nPatchesLoc = 0
       domains: do nn=1,nDom

         ! Have the pointers in blockPointers point to this block
         ! on the finest mg level. There is no need to distinguish
         ! between the different spectral solutions, because this
         ! info is the same for all of them.

         call setPointers(nn,1_intType,1_intType)

         ! Loop over the number of boundary subfaces of this block.

         bocos: do mm=1,nBocos

           ! Determine the family ID of the subface and check if it is
           ! a positive integer. In that case it belongs to a family.
           ! Otherwise it may still be stored if there are no families
           ! in the grid.

           jj = cgnsSubface(mm)
           jj = cgnsDoms(nbkGlobal)%bocoInfo(jj)%familyID

           if(jj > 0) then
             nPatchesLoc(jj) = nPatchesLoc(jj) + 1
           else
             if(cgnsNfamilies == 0 .and.             &
                (BCType(mm) == EulerWall        .or. &
                 BCType(mm) == NSWallAdiabatic  .or. &
                 BCType(mm) == NSWallIsothermal)) then
               nPatchesLoc(1) = nPatchesLoc(1) + 1
             endif
           endif

         enddo bocos
       enddo domains

       ! Gather the number of surface patches on all processors.

       size = max(cgnsNfamilies,1_intType)
       call mpi_allgather(nPatchesloc, size, sumb_integer, &
                          nPatches,    size, sumb_integer, &
                          SUmb_comm_world, ierr)

       ! Allocate the memory of mdNSurfPatches if it has not been
       ! allocated before.

       if(.not. allocated(mdNSurfPatches) ) then
         allocate(mdNSurfPatches(0:nProc,size), stat=ierr)
         if(ierr /= 0)                      &
           call terminate("mdCreateNPatch", &
                          "Memory allocation failure for &
                          &mdNSurfPatches")
       endif

       ! Create mdNSurfPatches, which is a cumulative storage version
       ! of nPatches. It is cumulative also for the families,
       ! such that one big array can be used to store the data of the
       ! different families.

       do nn=1,size
         if(nn == 1) then
           mdNSurfPatches(0,nn) = 0
         else
           mdNSurfPatches(0,nn) = mdNSurfPatches(nProc,nn-1)
         endif

         do mm=1,nProc
           mdNSurfPatches(mm,nn) = mdNSurfPatches(mm-1,nn) &
                                 + nPatches(nn,mm)
         enddo
       enddo

       end subroutine mdCreateNPatch
