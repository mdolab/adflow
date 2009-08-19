!
!      ******************************************************************
!      *                                                                *
!      * File:          mdNSurfNodesLocal.f90                           *
!      * Author:        C.A.(Sandy) Mader,Edwin van der Weide           *
!      * Starting date: 10-24-2007                                      *
!      * Last modified: 10-24-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mdCreateNsurfNodesLocal
!
!      ******************************************************************
!      *                                                                *
!      * mdCreateNsurfNodes determines the number of surface nodes      *
!      * nodes for all the families in the grid. If no family info is   *
!      * present, it will determine the total number of points on the   *
!      * surfaces.                                                      *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use communication
       use mdDataLocal
       implicit none
!
!      Local variables.
!
       integer :: ierr, size

       integer(kind=intType) :: jj, mm, nn, nni, nnj, nnk

       integer(kind=intType), &
         dimension(max(cgnsNfamilies,1_intType)) :: nNodesLoc
       integer(kind=intType), &
         dimension(max(cgnsNfamilies,1_intType),nProc) :: nNodes
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!

       ! Determine the local number of surface nodes per family.

       nNodesLoc = 0
       domains: do nn=1,nDom
          !print *,'ndoms',nn
         ! Have the pointers in blockPointers point to the 1st spectral
         ! solution of this block on the finest mg level. As the number
         ! of surface nodes is the same for all spectral modes, this
         ! is okay.

         call setPointers(nn,1_intType,1_intType)

         ! Loop over the number of boundary subfaces of this block.

         bocos: do mm=1,nBocos
            !print *,'bocos',mm
            !check to see whether family boundary conditions are present
            if (cgnsDoms(nbkGlobal)%BCFamilies==.true.)then
               !BC families are present
               
               ! Determine the family id of the subface and check if
               ! it is a positive integer. In that case it belongs to
               ! a family.
               
               jj = cgnsSubface(mm)
               jj = cgnsDoms(nbkGlobal)%bocoInfo(jj)%familyID
               
               familyTest: if(jj > 0) then
                  
                  ! Determine the number of nodes in every coordinate
                  ! direction and update mdNSurfNodesLoc.
                  
                  nni = inEnd(mm) - inBeg(mm) + 1
                  nnj = jnEnd(mm) - jnBeg(mm) + 1
                  nnk = knEnd(mm) - knBeg(mm) + 1
                  
                  nNodesLoc(jj) = nNodesLoc(jj) + nni*nnj*nnk

               else familyTest
                  if(myID == 0)                            &
                       call terminate("mdCreateSurfCoorList", &
                       "Family ID 0 is only allowed when no family &
                       &info is present in the grid")
               endif familyTest
            else
               ! Subface does not belong to a family. It is possible that
               ! no family info is present in the grid. In that case all
               ! solid wall boundary points are accumulated.

               if((BCType(mm) == EulerWall       .or. &
                 BCType(mm) == NSWallAdiabatic .or. &
                 BCType(mm) == NSWallIsothermal)) then          
!!$               if(cgnsNfamilies == 0 .and.            &
!!$                (BCType(mm) == EulerWall       .or. &
!!$                 BCType(mm) == NSWallAdiabatic .or. &
!!$                 BCType(mm) == NSWallIsothermal)) then

                  nni = inEnd(mm) - inBeg(mm) + 1
                  nnj = jnEnd(mm) - jnBeg(mm) + 1
                  nnk = knEnd(mm) - knBeg(mm) + 1
                  
                  nNodesLoc(1) = nNodesLoc(1) + nni*nnj*nnk
                  
               endif

            endif
         enddo bocos
      enddo domains

       ! Gather the number of surface nodes on all processors.

       size = max(cgnsNfamilies,1_intType)
      
!!$       call mpi_allgather(nNodesLoc, size, sumb_integer, &
!!$                          nNodes,    size, sumb_integer, &
!!$                          SUmb_comm_world, ierr)

       ! Test if the memory of mdNSurfNodes has already been allocated.
       ! If not, allocate it.

       if(.not. allocated(mdNSurfNodesLocal) ) then
         allocate(mdNSurfNodesLocal(size), stat=ierr)
         if(ierr /= 0)                          &
           call terminate("mdCreateNsurfNodesLocal", &
                          "Memory allocation failure for mdNSurfNodesLocal")
       endif

       ! Create mdNSurfNodesLocal, which is a cumulative storage format
       ! version of nNodes. It is cumulative also for the families,
       ! such that one big array can be used to store the data of the
       ! different families.

       do nn=1,size
          !print *,'size',nn,size,nNodesLoc(nn)
          mdNSurfNodesLocal(nn) =  nNodesLoc(nn)
!!$          if (nn==1) then
!!$             mdNSurfNodesLocal(2,nn) =  0
!!$          else
!!$             mdNSurfNodesLocal(2,nn) =  mdNSurfNodesLocal(2,nn-1)+&
!!$                  nNodesLoc(nn-1)
!!$          endif
       enddo

     end subroutine mdCreateNsurfNodesLocal

!      ==================================================================

       subroutine mdDeleteNsurfNodesLocal
!
!      ******************************************************************
!      *                                                                *
!      * mdDeleteNsurfNodes deallocates the memory of mdNSurfNodes.     *
!      *                                                                *
!      ******************************************************************
!
       use mdDataLocal
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
       ! Deallocate the memory of mdNSurfNodes
       ! if it has been allocated.

       if( allocated(mdNSurfNodesLocal) ) then

         deallocate(mdNSurfNodesLocal, stat=ierr)
         if(ierr /= 0)                          &
           call terminate("mdDeleteNsurfNodesLocal", &
                          "Deallocation error for mdNSurfNodesLocal")
       endif

     end subroutine mdDeleteNsurfNodesLocal
