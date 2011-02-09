!
!      ******************************************************************
!      *                                                                *
!      * File:          checkFaces.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-10-2004                                      *
!      * Last modified: 06-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine checkFaces
!
!      ******************************************************************
!      *                                                                *
!      * checkFaces determines whether or not a boundary condition or   *
!      * a connectivity has been specified for all block faces. If this *
!      * is not the case, the corresponding blocks are printed and the  *
!      * code terminates.                                               *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use communication
       use inputPhysics
       use inputTimeSpectral
       implicit none
!
!      Local variables.
!
       integer :: ierr, size

       integer, dimension(nProc) :: recvcounts, displs

       integer(kind=intType) :: mm, nn, sps, multiple
       integer(kind=intType) :: nBad, faceID, nBadGlobal

       integer(kind=intType), dimension(nProc) :: counts
       integer(kind=intType), &
               dimension(4,nDom*nTimeIntervalsSpectral) :: bad

       integer(kind=intType), dimension(:,:), allocatable :: badGlobal

       real(kind=realType) :: dummy

       logical :: blockIsBad

       character(len=7) :: intString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * Determine the local bad blocks.                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of spectral solutions to be checked and
       ! the number of local blocks.

       nBad = 0
       do sps=1,nTimeIntervalsSpectral
         do nn=1,nDom

           ! Check if the current block is okay.

           call setPointers(nn, 1_intType, sps)
           call checkFacesBlock(blockIsBad, faceID, multiple)

           ! If the block is bad, update nBad and store the info in bad.

           if( blockIsBad ) then
             nBad        = nBad + 1
             bad(1,nBad) = nbkGlobal
             bad(2,nBad) = faceID
             bad(3,nBad) = multiple
             bad(4,nBad) = sps
           endif

         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Determine the global number of bad blocks and gather this      *
!      * information.                                                   *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of bad blocks per processor.

       call mpi_allgather(nBad, 1, sumb_integer, counts, 1, &
                          sumb_integer, SUmb_comm_world, ierr)

       ! Determine the global number of bad blocks and the arrays
       ! recvcounts and displs needed for the call to allgatherv.

       nBadGlobal    = counts(1)
       recvcounts(1) = 4*counts(1)
       displs(1)     = 0

       do nn=2,nProc
         nBadGlobal     = nBadGlobal + counts(nn)
         recvcounts(nn) = 4*counts(nn)
         displs(nn)     = displs(nn-1) + recvcounts(nn-1)
       enddo

       ! Allocate the memory to store all the bad blocks.

       allocate(badGlobal(4,nBadGlobal), stat=ierr)
       if(ierr /= 0)                  &
         call terminate("checkFaces", &
                        "Memory allocation failure for badGlobal")

       ! Gather the data.

       size = 4*nBad
       call mpi_allgatherv(bad, size, sumb_integer, badGlobal, &
                           recvcounts, displs, sumb_integer,   &
                           SUmb_comm_world, ierr)

       ! Sort the bad blocks and get rid of the multiple entries.
       ! The last argument is .false. to indicate that only the
       ! integers must be sorted; dummy is passed for consistency.

       call sortBadEntities(nBadGlobal, badGlobal, dummy, .false.)

       ! If bad blocks are present, print them to stdout and terminate.

       testBadPresent: if(nBadGlobal > 0) then

         ! The data is only written by processor 0.

         testRootProc: if(myID == 0) then

           ! Write a header.

           write(intString,"(i6)") nBadGlobal
           intString = adjustl(intString)

           print "(a)", "#"
           print 101,  trim(intString)
           print "(a)", "# is not correct. &
                        &Here is the list of bad blocks"
           print "(a)", "#"

           ! Write the bad blocks.

           do mm=1,nBadGlobal

             ! Abbreviate the contents of badGlobal a bit easier.

             nn       = badGlobal(1,mm)
             faceID   = badGlobal(2,mm)
             multiple = badGlobal(3,mm)
             sps      = badGlobal(4,mm)

             ! If multiple spectral solutions are read, write the info
             ! about it. Otherwise just start the line with # sign.

             if(nTimeIntervalsSpectral > 1) then
               write(intString,"(i6)") nBadGlobal
               intString = adjustl(intString)
               write(*,102,advance="no") trim(intString)
             else
               write(*,"(a)",advance="no") "# "
             endif

             ! Write the error message, depending on the value of
             ! faceID and multiple.


             if( multiple == 0) then

               select case (faceID)
                 case (iMin)
                   print 103, trim(cgnsDoms(nn)%zoneName)
                 case (iMax)
                   print 104, trim(cgnsDoms(nn)%zoneName)
                 case (jMin)
                   print 105, trim(cgnsDoms(nn)%zoneName)
                 case (jMax)
                   print 106, trim(cgnsDoms(nn)%zoneName)
                 case (kMin)
                   print 107, trim(cgnsDoms(nn)%zoneName)
                 case (kMax)
                   print 108, trim(cgnsDoms(nn)%zoneName)
                 case default
                   print 109, trim(cgnsDoms(nn)%zoneName)
               end select

             else

               select case (faceID)
                 case (iMin)
                   print 110, trim(cgnsDoms(nn)%zoneName)
                 case (iMax)
                   print 111, trim(cgnsDoms(nn)%zoneName)
                 case (jMin)
                   print 112, trim(cgnsDoms(nn)%zoneName)
                 case (jMax)
                   print 113, trim(cgnsDoms(nn)%zoneName)
                 case (kMin)
                   print 114, trim(cgnsDoms(nn)%zoneName)
                 case (kMax)
                   print 115, trim(cgnsDoms(nn)%zoneName)
                 case default
                   print 116, trim(cgnsDoms(nn)%zoneName)
               end select

             endif

           enddo

 101       format("# Found ",a," blocks for which the boundary or &
                   &connectivity information")
 102       format("# Spectral grid ",a, ", ")
 103       format("Zone ",a, ": iMin block face not fully described")
 104       format("Zone ",a, ": iMax block face not fully described")
 105       format("Zone ",a, ": jMin block face not fully described")
 106       format("Zone ",a, ": jMax block face not fully described")
 107       format("Zone ",a, ": kMin block face not fully described")
 108       format("Zone ",a, ": kMax block face not fully described")
 109       format("Zone ",a, ": Multiple block faces not fully &
                  &described")
 110       format("Zone ",a, ": iMin block face multiple described")
 111       format("Zone ",a, ": iMax block face multiple described")
 112       format("Zone ",a, ": jMin block face multiple described")
 113       format("Zone ",a, ": jMax block face multiple described")
 114       format("Zone ",a, ": kMin block face multiple described")
 115       format("Zone ",a, ": kMax block face multiple described")
 116       format("Zone ",a, ": Multiple block faces multiple described")

           ! Terminate.

           call terminate("checkFaces", &
                          "Wrong block boundary info found")

         endif testRootProc

         ! The other processors wait to get killed.

         call mpi_barrier(SUmb_comm_world, ierr)

       endif testBadPresent

       ! Deallocate the memory of badGlobal.

       deallocate(badGlobal, stat=ierr)
       if(ierr /= 0)                   &
         call terminate("checkFaces", &
                        "Deallocation failure for badGlobal")

       end subroutine checkFaces

!      ==================================================================

       subroutine checkFacesBlock(blockIsBad, faceID, multiple)
!
!      ******************************************************************
!      *                                                                *
!      * checkFacesBlock checks if for the currently active block all   *
!      * the necessary boundary and connectivity info is specified.     *
!      * If not, blockIsBad is set to .true. and the block face ID      *
!      * which is bad, is returned as well.                             *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(out) :: faceID, multiple
       logical,               intent(out) :: blockIsBad
!
!      Local variables.
!
       integer, dimension(2:jl,2:kl), target :: iMinFace, iMaxFace
       integer, dimension(2:il,2:kl), target :: jMinFace, jMaxFace
       integer, dimension(2:il,2:jl), target :: kMinFace, kMaxFace

       integer, dimension(:,:), pointer :: face

       integer(kind=intType) :: nn, i, j, k, nBad
       integer(kind=intType) :: istart, iEnd, jStart, jEnd

       logical :: iMinBad, iMaxBad, jMinBad, jMaxBad
       logical :: kMinBad, kMaxBad
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize iMinFace, etc to 0. These arrays will contain the
       ! number of specified connectivities and BC's for each face on
       ! the block boundary.

       iMinFace = 0; iMaxFace = 0
       jMinFace = 0; jMaxFace = 0
       kMinFace = 0; kMaxFace = 0

       ! Loop over the number of subfaces of this block.

       do nn=1,nSubface

         ! Determine the block face on which the subface is located and
         ! set a couple of variables accordingly. Note that the nodal
         ! ranges are used to determine the correct range of faces,
         ! because the actual cell range can contain halo cells.

         select case(BCFaceID(nn))
           case (iMin)
             face => iMinFace; istart = min(jnBeg(nn),jnEnd(nn)) + 1
                               iEnd   = max(jnBeg(nn),jnEnd(nn))
                               jStart = min(knBeg(nn),knEnd(nn)) + 1
                               jEnd   = max(knBeg(nn),knEnd(nn))
           case (iMax)
             face => iMaxFace; istart = min(jnBeg(nn),jnEnd(nn)) + 1
                               iEnd   = max(jnBeg(nn),jnEnd(nn))
                               jStart = min(knBeg(nn),knEnd(nn)) + 1
                               jEnd   = max(knBeg(nn),knEnd(nn))
           case (jMin)
             face => jMinFace; istart = min(inBeg(nn),inEnd(nn)) + 1
                               iEnd   = max(inBeg(nn),inEnd(nn))
                               jStart = min(knBeg(nn),knEnd(nn)) + 1
                               jEnd   = max(knBeg(nn),knEnd(nn))
           case (jMax)
             face => jMaxFace; istart = min(inBeg(nn),inEnd(nn)) + 1
                               iEnd   = max(inBeg(nn),inEnd(nn))
                               jStart = min(knBeg(nn),knEnd(nn)) + 1
                               jEnd   = max(knBeg(nn),knEnd(nn))
           case (kMin)
             face => kMinFace; istart = min(inBeg(nn),inEnd(nn)) + 1
                               iEnd   = max(inBeg(nn),inEnd(nn))
                               jStart = min(jnBeg(nn),jnEnd(nn)) + 1
                               jEnd   = max(jnBeg(nn),jnEnd(nn))
           case (kMax)
             face => kMaxFace; istart = min(inBeg(nn),inEnd(nn)) + 1
                               iEnd   = max(inBeg(nn),inEnd(nn))
                               jStart = min(jnBeg(nn),jnEnd(nn)) + 1
                               jEnd   = max(jnBeg(nn),jnEnd(nn))
         end select

         ! Loop over the faces and update their counter.

         do j=jStart,jEnd
           do i=istart,iEnd
             face(i,j) = face(i,j) + 1
           enddo
         enddo

       enddo

       ! Determine the bad block faces.

       multiple = 0
       iMinBad  = .false.
       iMaxBad  = .false.

       ! iMin and iMax face.

       do k=2,kl
         do j=2,jl
           if(iMinFace(j,k) /= 1) then
             multiple = iMinFace(j,k)
             faceID   = iMin
             iMinBad  = .true.
           endif

           if(iMaxFace(j,k) /= 1) then
             multiple = iMaxFace(j,k)
             faceID   = iMax
             iMaxBad  = .true.
           endif
         enddo
       enddo

       ! jMin and jMax face.

       jMinBad = .false.
       jMaxBad = .false.

       do k=2,kl
         do i=2,il
           if(jMinFace(i,k) /= 1) then
             multiple = jMinFace(i,k)
             faceID   = jMin
             jMinBad  = .true.
           endif

           if(jMaxFace(i,k) /= 1) then
             multiple = jMaxFace(i,k)
             faceID   = jMax
             jMaxBad  = .true.
           endif
         enddo
       enddo

       ! kMin and kMax face.

       kMinBad = .false.
       kMaxBad = .false.

       do j=2,jl
         do i=2,il
           if(kMinFace(i,j) /= 1) then
             multiple = kMinFace(i,j)
             faceID   = kMin
             kMinBad  = .true.
           endif

           if(kMaxFace(i,j) /= 1) then
             multiple = kMaxFace(i,j)
             faceID   = kMax
             kMaxBad  = .true.
           endif
         enddo
       enddo

       ! Determine the number of bad block faces.

       nBad = 0
       if(iMinBad) nBad = nBad + 1
       if(iMaxBad) nBad = nBad + 1
       if(jMinBad) nBad = nBad + 1
       if(jMaxBad) nBad = nBad + 1
       if(kMinBad) nBad = nBad + 1
       if(kMaxBad) nBad = nBad + 1

       ! Set blockIsBad if bad faces are present and correct the face
       ! id to something weird if multiple bad faces are present.

       blockIsBad = .false.
       if(nBad > 0) blockIsBad = .true.
       if(nBad > 1) faceID = huge(faceID)

       end subroutine checkFacesBlock
