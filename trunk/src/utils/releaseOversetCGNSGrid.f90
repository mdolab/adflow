!
!      ******************************************************************
!      *                                                                *
!      * File:          releaseOversetCGNSGrid.f90                      *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 08-02-2005                                      *
!      * Last modified: 08-30-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine releaseOversetCGNSGrid
!
!      ******************************************************************
!      *                                                                *
!      * releaseOversetCGNSGrid deallocates the memory for all of the   *
!      * overset connectivities and holes stored in the CGNS domains    *
!      * for the global blocks.                                         *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the cgns domains.

       do nn=1,cgnsNDom
         call releaseOversetCGNSZone(nn)
       end do

       end subroutine releaseOversetCGNSGrid
!
!      ==================================================================
!
       subroutine releaseOversetCGNSZone(nn)
!
!      ******************************************************************
!      *                                                                *
!      * releaseOversetCGNSZone deallocates the memory for all of the   *
!      * overset connectivities and holes for CGNS domain nn.           *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
!
!      Local variables.
!
       integer :: ierr
       integer(kind=intType) :: i
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Deallocate all of the overset connectivities.

       do i = 1,cgnsDoms(nn)%nOverset
         deallocate(cgnsDoms(nn)%connOver(i)%ibndry, &
                    cgnsDoms(nn)%connOver(i)%idonor, &
                    cgnsDoms(nn)%connOver(i)%interp, stat=ierr)
         if(ierr /= 0) &
           call terminate("releaseOversetCGNSZone", &
                          "Deallocation failure for connectivity data")
       end do

       deallocate(cgnsDoms(nn)%connOver, stat=ierr)
       if(ierr /= 0) &
         call terminate("releaseOversetCGNSZone", &
                        "Deallocation failure for connectivities")

       ! Deallocate all of the hole sets.

       do i = 1,cgnsDoms(nn)%nHoles
         deallocate(cgnsDoms(nn)%hole(i)%indices, stat=ierr)
         if(ierr /= 0) &
           call terminate("releaseOversetCGNSZone", &
                          "Deallocation failure for hole indices")
       end do

       deallocate(cgnsDoms(nn)%hole, stat=ierr)
       if(ierr /= 0) &
         call terminate("releaseOversetCGNSZone", &
                        "Deallocation failure for hole sets")

       ! Set some variables to 0 to be consistent.

       cgnsDoms(nn)%nOverset      = 0
       cgnsDoms(nn)%nCellsOverset = 0
       cgnsDoms(nn)%nHoles        = 0

       end subroutine releaseOversetCGNSZone
