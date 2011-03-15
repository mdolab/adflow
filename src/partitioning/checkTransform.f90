!
!      ******************************************************************
!      *                                                                *
!      * File:          checkTransform.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-04-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine checkTransform(transform, nZone, n1to1, printWarning)
!
!      ******************************************************************
!      *                                                                *
!      * checkTransform checks the transformation matrix between this   *
!      * zone and the donor for the given subrange. In case an error is *
!      * found it is tried to correct this.                             *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use cgnsGrid
       use communication
       implicit none
!
!      Subroutine arguments.
!
       integer, intent(in)                  :: nZone, n1to1
       integer, dimension(3), intent(inout) :: transform
       logical, intent(in)                  :: printWarning
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nDirFace, nDirDonor
       integer(kind=intType) :: sumTransform
       integer(kind=intType) :: l1, L2, l3

       integer(kind=intType), dimension(3) :: haloDir, donorDir
       integer(kind=intType), dimension(3,2) :: zoneRange, donorRange
       integer(kind=intType), dimension(3,3) :: tMat

       character(len=maxCGNSNameLen) :: zoneName, connectName
       character(len=2*maxStringLen)  :: errorMessage
!
!      Function definitions.
!
       integer(kind=intType) :: delta
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Copy the zoneName and connectName, just for readability later.

       zoneName    = cgnsDoms(nZone)%zoneName
       connectName = cgnsDoms(nZone)%conn1to1(n1to1)%connectName

       ! Copy the zone and donor range into zoneRange and donorRange.

       zoneRange(1,1) = cgnsDoms(nZone)%conn1to1(n1to1)%iBeg
       zoneRange(2,1) = cgnsDoms(nZone)%conn1to1(n1to1)%jBeg
       zoneRange(3,1) = cgnsDoms(nZone)%conn1to1(n1to1)%kBeg

       zoneRange(1,2) = cgnsDoms(nZone)%conn1to1(n1to1)%iEnd
       zoneRange(2,2) = cgnsDoms(nZone)%conn1to1(n1to1)%jEnd
       zoneRange(3,2) = cgnsDoms(nZone)%conn1to1(n1to1)%kEnd

       donorRange(1,1) = cgnsDoms(nZone)%conn1to1(n1to1)%diBeg
       donorRange(2,1) = cgnsDoms(nZone)%conn1to1(n1to1)%djBeg
       donorRange(3,1) = cgnsDoms(nZone)%conn1to1(n1to1)%dkEnd

       donorRange(1,2) = cgnsDoms(nZone)%conn1to1(n1to1)%diEnd
       donorRange(2,2) = cgnsDoms(nZone)%conn1to1(n1to1)%djEnd
       donorRange(3,2) = cgnsDoms(nZone)%conn1to1(n1to1)%dkEnd

       ! Determine the normal direction for the subface and do a trivial
       ! check to see if there is one.

       do nDirFace=1,3
         if(zoneRange(nDirFace,1) == zoneRange(nDirFace,2)) exit
       enddo

       if(nDirFace > 3) then
         if(myID == 0) then
           write(errorMessage,100) trim(connectName), trim(zoneName)
 100       format("1 to 1 subface",1X,A,1X,"of zone",1X,A, &
                  ": No constant index found")
           call terminate("checkTransform", errorMessage)
         endif

         ! Make sure that other processors wait until they are killed.

         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Idem for the donor subface.

       do nDirDonor=1,3
         if(donorRange(nDirDonor,1) == donorRange(nDirDonor,2)) exit
       enddo

       if(nDirDonor > 3) then
         if(myID == 0) then
           write(errorMessage,110) trim(connectName), trim(zoneName)
 110       format("1 to 1 subface",1X,A,1X,"of zone",1X,A, &
                  ": No constant index found for donor")
           call terminate("checkTransform", errorMessage)
         endif

         ! Make sure that other processors wait until they are killed.

         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Check if the sum of the absolute values of transform equals 6.
       ! If not, assume that the normal direction is not set correctly.

       sumTransform = abs(transform(1)) + abs(transform(2)) &
                    + abs(transform(3))
       if(sumTransform /= 6) then

         ! Change the normal direction of transform and check the
         ! sum again.

         transform(nDirFace) = nDirDonor
         sumTransform = abs(transform(1)) + abs(transform(2)) &
                      + abs(transform(3))

         if(sumTransform /= 6) then

           ! Something seriously wrong. I cannot repair this.

           if(myID == 0) then
             write(errorMessage,120) trim(connectName), trim(zoneName)
 120         format("1 to 1 subface",1X,A,1X,"of zone",1X,A,  &
                    ": Something seriously wrong with the &
                    &transformation matrix")
             call terminate("checkTransform", errorMessage)
           endif

           ! Make sure that other processors wait until they are killed.

           call mpi_barrier(SUmb_comm_world, ierr)
         else
           ! Repair successful, although the orientation might be wrong.
           ! This will be checked later. Anyway print a warning message
           ! if desired.

           if(myID == 0 .and. printWarning) then
             print "(a)", "#"
             print "(a)", "#                          Warning"
             print 130, trim(connectName), trim(zoneName)
 130         format("# 1 to 1 subface",1X,A,1X,"of zone",1X,A, &
                    ": Normal component of the transformation&
                    & matrix successfully corrected.")
             print "(a)", "#"
           endif
         endif
       endif

       ! Create the halo vector for the current zone. This vector is
       ! pointing outwards, i.e. in direction of the donor block.

       haloDir = 0
       haloDir(nDirFace) = 1
       if(zoneRange(nDirFace,1) == 1) haloDir(nDirFace) = -1

       ! Idem for the donor. Also this vector points from the current
       ! block to the donor block, albeit in donor block coordinates.

       donorDir = 0
       donorDir(nDirDonor) = -1
       if(donorRange(nDirDonor,1) == 1) donorDir(nDirDonor) = 1

       ! Determine the full transformation matrix.

       l1 = transform(1)
       L2 = transform(2)
       l3 = transform(3)

       tMat(1,1) = sign(1_intType,l1) * delta(l1,1_intType)
       tMat(2,1) = sign(1_intType,l1) * delta(l1,2_intType)
       tMat(3,1) = sign(1_intType,l1) * delta(l1,3_intType)

       tMat(1,2) = sign(1_intType,l2) * delta(l2,1_intType)
       tMat(2,2) = sign(1_intType,l2) * delta(l2,2_intType)
       tMat(3,2) = sign(1_intType,l2) * delta(l2,3_intType)

       tMat(1,3) = sign(1_intType,l3) * delta(l3,1_intType)
       tMat(2,3) = sign(1_intType,l3) * delta(l3,2_intType)
       tMat(3,3) = sign(1_intType,l3) * delta(l3,3_intType)

       ! Apply the transformation matrix to haloDir.

       l1 = haloDir(1)
       L2 = haloDir(2)
       l3 = haloDir(3)

       haloDir(1) = tMat(1,1)*l1 + tMat(1,2)*l2 + tMat(1,3)*l3
       haloDir(2) = tMat(2,1)*l1 + tMat(2,2)*l2 + tMat(2,3)*l3
       haloDir(3) = tMat(3,1)*l1 + tMat(3,2)*l2 + tMat(3,3)*l3

       ! If the transformation matrix is correct haloDir == donorDir.
       ! If this is not the case, there are two possibilities. Either
       ! the directions are just reversed, which means that the
       ! corresponding element of transform must be reversed, or they
       ! are really different. In the latter case it cannot be
       ! corrected here and the grid file must be adapted.

       if(haloDir(nDirDonor) == 0) then

         ! Something seriously wrong. Exit the program.

         if(myID == 0) then
           write(errorMessage,140) trim(connectName), trim(zoneName)
 140         format("1 to 1 subface",1X,A,1X,"of zone",1X,A,  &
                    ": Something seriously wrong with the &
                    &transformation matrix")
           call terminate("checkTransform", errorMessage)
         endif

         ! Make sure that other processors wait until they are killed.

         call mpi_barrier(SUmb_comm_world, ierr)

       else if(haloDir(nDirDonor)*donorDir(nDirDonor) < 0) then

         ! Simply reverse the sign of the corresponding entry in
         ! transform. Processor 0 prints a warning message if desired.

         transform(nDirFace) = -transform(nDirFace)

         if(myID == 0 .and. printWarning) then
           print "(a)", "#"
           print "(a)", "#                            Warning"
           print 150, trim(connectName), trim(zoneName)
 150       format("# 1 to 1 subface",1X,A,1X,"of zone",1X,A, &
                  ": Normal component of the transformation&
                  & matrix reversed")
           print "(a)", "#"
         endif

       endif

       end subroutine checkTransform
