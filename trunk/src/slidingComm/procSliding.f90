!
!      ******************************************************************
!      *                                                                *
!      * File:          procSliding.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-26-2003                                      *
!      * Last modified: 03-01-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine procSliding(distrSliding)
!
!      ******************************************************************
!      *                                                                *
!      * procSliding determines for each sliding mesh interface the     *
!      * processors which contribute.                                   *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use communication
       use interfaceGroups
       use tmpSliding
       implicit none
!
!      Subroutine arguments.
!
       type(tmpSlidingType), dimension(*), intent(out) :: distrSliding
!
!      Local variables.
!
       integer :: ierr, ns

       integer, dimension(nProc)   :: recvcounts
       integer, dimension(nProc+1) :: displs

       integer(kind=intType) :: i, j, k, nn, mm, ii

       integer(kind=intType), dimension(max(cgnsNSliding,1_intType),2) :: &
                                                             nFaces, tmp
       integer(kind=intType), dimension(max(cgnsNSliding,1_intType)) :: &
                                                                  buffer

       integer(kind=intType), dimension(:), allocatable :: bufferRecv
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of faces of each sliding mesh interface
       ! on this processor. Both parts of the interface are considered.

       do nn=1,cgnsNSliding
         nFaces(nn,1) = 0
         nFaces(nn,2) = 0
       enddo

       do nn=1,nDom

         ! Set the pointers for this block; note that only the 1st
         ! time spectral interval needs to be considered.

         call setPointers(nn,1_intType,1_intType)

         ! Loop over the boundary subfaces and find the subfaces
         ! on the sliding interfaces.

         do mm=1,nBocos
           if(BCType(mm) == slidingInterface) then

             i = max(inEnd(mm) - inBeg(mm), 1_intType)
             j = max(jnEnd(mm) - jnBeg(mm), 1_intType)
             k = max(knEnd(mm) - knBeg(mm), 1_intType)

             ii = groupNum(mm)
             if(ii < 0) then
               nFaces(-ii,1) = nFaces(-ii,1) + k*j*i
             else
               nFaces(ii,2)  = nFaces(ii,2)  + k*j*i
             endif
           endif
         enddo
       enddo

       ! First determine for every interface the maximum number of faces
       ! stored on one processor.

       do nn=1,cgnsNSliding
         tmp(nn,1) = nFaces(nn,1) + nFaces(nn,2)
       enddo

       ns = cgnsNSliding
       call mpi_allreduce(tmp, buffer, ns, sumb_integer, mpi_max, &
                          SUmb_comm_world, ierr)

       do nn=1,cgnsNSliding
         distrSliding(nn)%nFaceMax = buffer(nn)
       enddo

       ! Determine for every sliding mesh interface the number of
       ! processors for both sides of the interface. Also determine
       ! the sum, because this is needed in mpi_allgatherv.

       do nn=1,cgnsNSliding
         if(nFaces(nn,1) > 0) nFaces(nn,1) = 1
         if(nFaces(nn,2) > 0) nFaces(nn,2) = 1
       enddo

       call mpi_allreduce(nFaces, tmp, 2*ns, sumb_integer, mpi_sum, &
                          SUmb_comm_world, ierr)

       i = 0
       j = 0
       do nn=1,cgnsNSliding
         distrSliding(nn)%nProcs1 = tmp(nn,1)
         distrSliding(nn)%nProcs2 = tmp(nn,2)

         i = i + tmp(nn,1)
         j = j + tmp(nn,2)

         ! Allocate the memory for the processor id's.

         allocate(distrSliding(nn)%procs1(tmp(nn,1)), &
                  distrSliding(nn)%procs2(tmp(nn,2)), stat=ierr)

         if(ierr /= 0)                   &
           call terminate("procSliding", &
                          "Memory allocation failure for procs1 &
                          &and procs2")
       enddo

       ! Allocate the memory for the receive buffer in mpi_allgatherv.

       i = max(i,j)
       allocate(bufferRecv(i), stat=ierr)
       if(ierr /= 0)                   &
         call terminate("procSliding", &
                        "Memory allocation failure for bufferRecv")

       ! Loop over the two sides of a sliding mesh interface.

       sides: do mm=1,2

         ! Determine the number of sides mm of the sliding mesh
         ! interfaces this processor contributes to. Note that integers
         ! are used here, because of later usage of mpi_allgatherv.
         ! Furthermore nFaces(..,1) is used as a buffer to send the
         ! data. This is to avoid a crossing of array boundaries in
         ! the call to mpi_allgatherv when no sliding mesh interfaces
         ! are present.

         ns = 0
         do nn=1,cgnsNSliding
           if(nFaces(nn,mm) > 0) then
             ns = ns + 1
             nFaces(ns,1) = nn
           endif
         enddo

         ! Determine the number of integers which will be received from
         ! every processor in the mpi_allgatherv later on. Note the usage
         ! of MPI_integer and not sumb_integer.

         call mpi_allgather(ns, 1, mpi_integer, recvcounts, 1, &
                            mpi_integer, SUmb_comm_world, ierr)

         ! Create the array displs, needed in mpi_allgatherv.
         ! The last entry is not needed in mpi_allgatherv, but is
         ! stored for later use.

         displs(1) = 0
         do nn=1,nProc
           displs(nn+1) = displs(nn) + recvcounts(nn)
         enddo

         ! Use mpi_allgatherv to gather the data from all processors.

         call mpi_allgatherv(nFaces, ns, sumb_integer, bufferRecv, &
                             recvcounts, displs, sumb_integer,     &
                             SUmb_comm_world, ierr)

         ! Now the data in bufferRecv must be stored in the
         ! correct place in distrSliding. Use buffer as a counter,
         ! which must must initialized first.

         do nn=1,cgnsNSliding
           buffer(nn) = 0
         enddo

         do nn=1,nProc
           do ii=(displs(nn)+1),displs(nn+1)

             ! Store the sliding mesh interface in i and the new
             ! index in j.

             i         = bufferRecv(ii)
             buffer(i) = buffer(i) + 1
             j         = buffer(i)

             ! Store the appropriate processor in the correct
             ! place in distrSliding. As the processor numbering
             ! starts at 0, 1 must be substracted from nn.

             if(mm == 1) then
               distrSliding(i)%procs1(j) = nn - 1
             else
               distrSliding(i)%procs2(j) = nn - 1
             endif

           enddo
         enddo

       enddo sides

       ! Release the memory of bufferRecv.

       deallocate(bufferRecv, stat=ierr)
       if(ierr /= 0)                   &
         call terminate("procSliding", &
                        "Deallocation failure for bufferRecv")

       end subroutine procSliding
