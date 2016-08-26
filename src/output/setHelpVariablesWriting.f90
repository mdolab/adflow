       subroutine setHelpVariablesWriting
!
!       setHelpVariablesWriting determines the variables, which are    
!       needed to write the CGNS files.                                
!
       use block
       use cgnsGrid
       use communication
       use monitor
       use outputMod
       use utils, only : terminate
       implicit none
!
!      Local variables.
!
       integer :: ierr, nSend
       integer, dimension(nProc) :: recvCounts, displs

       integer(kind=intType) :: i, nn

       integer(kind=intType), dimension(cgnsNDom) :: tmp
       integer(kind=intType), dimension(4,nDom)   :: buffer

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

       end subroutine setHelpVariablesWriting
