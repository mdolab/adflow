!
!      ******************************************************************
!      *                                                                *
!      * File:          setHelpVariablesWriting.f90                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-10-2005                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setHelpVariablesWriting
!
!      ******************************************************************
!      *                                                                *
!      * setHelpVariablesWriting determines the variables, which are    *
!      * needed to write the CGNS files.                                *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsGrid
       use communication
       use monitor
       use outputMod
       implicit none
!
!      Local variables.
!
       integer :: ierr, nSend
       integer, dimension(nProc) :: recvCounts, displs

       integer(kind=intType) :: i, nn

       integer(kind=intType), dimension(cgnsNDom) :: tmp
       integer(kind=intType), dimension(4,nDom)   :: buffer
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine for each CGNS block how many (sub) blocks are stored
       ! on this processor. Note that this info is the same for all
       ! spectral solutions, so the 1st is fine.

       allocate(nBlocksCGNSblock(0:cgnsNDom), blocksCGNSblock(nDom), &
                stat=ierr)
       if(ierr /= 0)                               &
         call terminate("setHelpVariablesWriting", &
                        "Memory allocation failure for &
                        &nBlocksCGNSblock and blocksCGNSblock.")

       nBlocksCGNSblock = 0
       do nn=1,nDom
         i = flowDoms(nn,1,1)%cgnsBlockID
         nBlocksCGNSblock(i) = nBlocksCGNSblock(i) + 1
       enddo

       ! Put nBlocksCGNSblock in cumulative storage format.
       ! Store this accumulated value in tmp, which serves as
       ! a counter later on.

       do i=1,cgnsNDom
         tmp(i)              = nBlocksCGNSblock(i-1)
         nBlocksCGNSblock(i) = nBlocksCGNSblock(i) + tmp(i)
       enddo

       ! Determine the values for blocksCGNSblock.

       do nn=1,nDom
         i = flowDoms(nn,1,1)%cgnsBlockID
         tmp(i) = tmp(i) + 1
         blocksCGNSblock(tmp(i)) = nn
       enddo

       ! If a grid file is being written and there are overset
       ! connectivities present, then build the mapping arrays needed to
       ! convert local indices and block numbers for donors.

       if(writeGrid .and. oversetPresent) then

         ! Gather the number of domains on each processor, put it in
         ! cumulative storage format, and build the arrays needed for
         ! the gathering.

         allocate(nDomPerProc(0:nProc), stat=ierr)
         if(ierr /= 0)                               &
           call terminate("setHelpVariablesWriting", &
                          "Memory allocation failure for nDomPerProc")

         nDomPerProc(0) = 0
         call mpi_allgather(nDom, 1, sumb_integer, nDomPerProc(1), 1, &
                            sumb_integer, SUmb_comm_world, ierr)

         recvCounts = 4*nDomPerProc(1:nProc)
         do nn=2,nProc
           nDomPerProc(nn) = nDomPerProc(nn) + nDomPerProc(nn-1)
         end do
         displs = 4*nDomPerProc(0:nProc-1)

         ! Allocate memory for the mapping array, fill the buffer, and
         ! gather the data from all processors.

         allocate(IDsBegOrAllDoms(4,nDomPerProc(nProc)), stat=ierr)
         if(ierr /= 0)                               &
           call terminate("setHelpVariablesWriting", &
                          "Memory allocation failure for &
                          &IDsBegOrAllDoms")

         do nn=1,nDom
           buffer(1,nn) = flowDoms(nn,1,1)%cgnsBlockID
           buffer(2,nn) = flowDoms(nn,1,1)%iBegOr
           buffer(3,nn) = flowDoms(nn,1,1)%jBegOr
           buffer(4,nn) = flowDoms(nn,1,1)%kBegOr
         enddo

         nSend = 4*nDom
         call mpi_allgatherv(buffer, nSend, sumb_integer,         &
                             IDsBegOrAllDoms, recvCounts, displs, &
                             sumb_integer, SUmb_comm_world, ierr)
       endif

       end subroutine setHelpVariablesWriting
