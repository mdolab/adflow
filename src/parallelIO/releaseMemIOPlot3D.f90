!
!      ******************************************************************
!      *                                                                *
!      * File:          releaseMemIOPlot3D.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-28-2005                                      *
!      * Last modified: 10-31-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine releaseMemIOPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * releaseMemIOPlot3D releases the memory of the variables needed *
!      * to perform the parallel IO. These variables are stored in the  *
!      * module IOModule.                                               *
!      *                                                                *
!      ******************************************************************
!
       use IOModule
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of IO parts and release the memory of
       ! P3D_IOParts. After that release P3D_IOParts itself.

       do nn=1,P3D_nIOParts
         deallocate(P3D_IOParts(nn)%blockID,  &
                    P3D_IOParts(nn)%indexW,   &
                    P3D_IOParts(nn)%posLocal, &
                    P3D_IOParts(nn)%indices,  &
                    P3D_IOParts(nn)%posComm,  &
                    P3D_IOParts(nn)%posNonLocal, stat=ierr)
         if(ierr /= 0)                          &
           call terminate("releaseMemIOPlot3D", &
                          "Deallocation failure for the member &
                          &variables of P3D_IOParts(nn)")
       enddo

       deallocate(P3D_IOParts, stat=ierr)
       if(ierr /= 0)                          &
         call terminate("releaseMemIOPlot3D", &
                        "Deallocation failure for P3D_IOParts")

       ! Release the memory of P3D_commPart. Note that posComm and
       ! posNonLocal is not allocated for P3D_commPart.

       deallocate(P3D_commPart%blockID,  &
                  P3D_commPart%indexW,   &
                  P3D_commPart%posLocal, &
                  P3D_commPart%indices, stat=ierr)
       if(ierr /= 0)                          &
         call terminate("releaseMemIOPlot3D", &
                        "Deallocation failure for the member &
                        &variables of P3D_commPart")

       ! Release the memory of the communication arrays.

       deallocate(P3D_procSend, P3D_procRecv, &
                  P3D_sendSize, P3D_recvSize, stat=ierr)
       if(ierr /= 0)                          &
         call terminate("releaseMemIOPlot3D", &
                        "Deallocation failure for the communication &
                        &arrays")

       ! Release the memory of the variables needed to write the
       ! record integers. Only if these have been allocated.

       if( allocated(P3D_recordIntegersWrite) ) then
         deallocate(P3D_recordIntegersWrite, &
                    P3D_recordPosition, stat=ierr)
         if(ierr /= 0)                          &
           call terminate("releaseMemIOPlot3D", &
                          "Deallocation failure for &
                          &P3D_recordIntegersWrite and &
                          &P3D_recordPosition")
       endif

       end subroutine releaseMemIOPlot3D
