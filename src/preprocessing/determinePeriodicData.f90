!
!      ******************************************************************
!      *                                                                *
!      * File:          determinePeriodicData.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-29-2007                                      *
!      * Last modified: 11-30-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determinePeriodicData(entityHalo,   nHalo, &
                                        externalComm, internalComm)
!
!      ******************************************************************
!      *                                                                *
!      * determinePeriodicData determines the periodic transformation   *
!      * for both the external and the internal communication patterns. *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use haloList
       use periodicInfo
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)            :: nHalo
       type(haloListType), dimension(:), intent(in) :: entityHalo

       type(commType),         intent(inout) :: externalComm
       type(internalCommType), intent(inout) :: internalComm
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, nn
       integer(kind=intType) :: nPerHalos, nInternal, nExternal

       integer(kind=intType), dimension(5) :: tmp

       type(periodicSubfacesHaloType), dimension(:), allocatable :: &
                                                                 periodic
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the halo's and determine the number of halo's
       ! for which periodic transformations are present.

       nn = 0
       do i=1,nHalo
         if(entityHalo(i)%nPeriodicSubfaces > 0) nn = nn + 1
       enddo

       ! Allocate the memory for periodic and copy the relevant data
       ! from entityHalo. Note that a shallow copy is made, i.e. the
       ! array for periodicSubfaces is not allocated. It just points
       ! to the entry in entityHalo.

       nPerHalos = nn
       allocate(periodic(nn), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("determinePeriodicData", &
                        "Memory allocation failure for periodic")

       nn = 0
       do i=1,nHalo
         if(entityHalo(i)%nPeriodicSubfaces > 0) then
           nn = nn + 1
           periodic(nn)%indexInHaloList   = i
           periodic(nn)%nPeriodicSubfaces = entityHalo(i)%nPeriodicSubfaces
           periodic(nn)%periodicSubfaces => entityHalo(i)%periodicSubfaces

           if(entityHalo(i)%donorProc == myID) then
             periodic(nn)%internalHalo = .true.
           else
             periodic(nn)%internalHalo = .false.
           endif
         endif
       enddo

       ! Make sure that the numbers of the periodic subfaces are
       ! sorted in increasing order. This is important for the sorting
       ! of the derived datatype periodicSubfacesHaloType.

       do i=1,nPerHalos
         do nn=1,periodic(i)%nPeriodicSubfaces
           tmp(nn) = periodic(i)%periodicSubfaces(nn)
         enddo

         call qsortIntegers(tmp, periodic(i)%nPeriodicSubfaces)

         do nn=1,periodic(i)%nPeriodicSubfaces
           periodic(i)%periodicSubfaces(nn) = tmp(nn)
         enddo
       enddo

       ! Sort periodic in increasing order, such that the number
       ! of different transformations can be determined.

       call qsortPeriodicSubfacesHaloType(periodic, nPerHalos)

       ! Determine the number of periodic halo's, which are also
       ! internal. These are now numbered first in periodic.

       nInternal = 0
       do i=1,nPerHalos
         if( periodic(i)%internalHalo ) nInternal = nInternal + 1
       enddo

       nExternal = nPerHalos - nInternal

       ! Call the internal subroutine to do the work for the
       ! internal and external communication pattern.

       nullify(internalComm%periodicData)
       if(nInternal == 0) then
         internalComm%nPeriodic = 0
       else
         call setPeriodicData(periodic, nInternal,    &
                              internalComm%nPeriodic, &
                              internalComm%periodicData)
       endif

       nullify(externalComm%periodicData)
       if(nExternal == 0) then
         externalComm%nPeriodic = 0
       else
         call setPeriodicData(periodic(nInternal+1:), nExternal, &
                              externalComm%nPeriodic,            &
                              externalComm%periodicData)
       endif

       ! Release the memory of period again. As the pointer was not
       ! allocated there is no need to release it either.

       deallocate(periodic, stat=ierr)
       if(ierr /= 0)                             &
         call terminate("determinePeriodicData", &
                        "Deallocation failure for periodic")

       !=================================================================

       contains

         !===============================================================

         subroutine setPeriodicData(perHalo, nPerHalo, nPeriodic, &
                                    periodicData)
!
!        ****************************************************************
!        *                                                              *
!        * setPeriodicData stores the periodic transformations and the  *
!        * corresponding halo's to which it must be applied in          *
!        * nPeriodic and periodicData. These variables are part of      *
!        * either the internal or external communication pattern.       *
!        *                                                              *
!        ****************************************************************
!
         use cgnsGrid
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=intType), intent(in) :: nPerHalo
         type(periodicSubfacesHaloType), dimension(:), intent(in) :: &
                                                                 perHalo

         integer(kind=intType), intent(out) :: nPeriodic
         type(periodicDataType), dimension(:), pointer :: periodicData
!
!        Local variables.
!
         integer(kind=intType) :: nn, mm, ll, kk
         integer(kind=intType) :: cgnsBlock, cgnsSub

         integer(kind=intType), dimension(0:nPerHalo) :: nHaloPerTransform

         integer(kind=intType), pointer, dimension(:)   :: blockID
         integer(kind=intType), pointer, dimension(:,:) :: indices

         real(kind=realType) :: theta, phi, psi
         real(kind=realType) :: cosTheta, cosPhi, cosPsi
         real(kind=realType) :: sinTheta, sinPhi, sinPsi

         real(kind=realType), dimension(3,3) :: rotMat, rotMatrix, tmpMat
         real(kind=realType), dimension(3)   :: trans, translation
         real(kind=realType), dimension(3)   :: rotCenter
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Determine the number of different periodic transformations
         ! as well as the number of halo's per transformation.
         ! Note that the operator == of periodicSubfacesHaloType only
         ! compares the periodic subface information of the halos.

         nHaloPerTransform(0) = 0

         nn = min(1_intType, nPerHalo)
         nHaloPerTransform(nn) = nn

         do i=2,nPerHalo
           if(perHalo(i) == perHalo(i-1)) then
             nHaloPerTransform(nn) = nHaloPerTransform(nn) + 1
           else
             nn = nn + 1
             nHaloPerTransform(nn) = nHaloPerTransform(nn-1) + 1
           endif
         enddo

         ! Set nPeriodic and allocate the memory for periodicData.
         ! If there are no periodic transformations return after
         ! nPeriodic is set.

         nPeriodic = nn
         if(nPeriodic == 0) return

         allocate(periodicData(nn), stat=ierr)
         if(ierr /= 0)                       &
           call terminate("setPeriodicData", &
                          "Memory allocation failure for periodicData")

         ! Loop over the number of periodic transformations.

         periodicLoop: do nn=1,nPeriodic

           ! Determine the number of halo's for this transformation
           ! and allocate the memory for the block ID and indices.
           ! Set pointers for readability.

           mm = nHaloPerTransform(nn) - nHaloPerTransform(nn-1)

           periodicData(nn)%nHalos = mm
           allocate(periodicData(nn)%block(mm), &
                    periodicData(nn)%indices(mm,3), stat=ierr)
           if(ierr /= 0)                       &
             call terminate("setPeriodicData", &
                            "Memory allocation failure for block &
                            &and indices.")

           blockID => periodicData(nn)%block
           indices => periodicData(nn)%indices

           ! Loop over the number of periodic halo's and set the
           ! blockID and indices.

           kk = 0
           do mm=(nHaloPerTransform(nn-1)+1),nHaloPerTransform(nn)
             kk = kk + 1
             ll = perHalo(mm)%indexInHaloList

             blockID(kk)   = entityHalo(ll)%myBlock
             indices(kk,1) = entityHalo(ll)%myI
             indices(kk,2) = entityHalo(ll)%myJ
             indices(kk,3) = entityHalo(ll)%myK
           enddo

           ! Determine the rotation matrix, rotation center and
           ! translation vector of this periodic transformation.
           ! The sorting has been such that all the halo's set above
           ! have the same transformation, so it is enough to consider
           ! the first element.

           kk = nHaloPerTransform(nn-1) + 1

           ! Set the rotation center to the rotation center of the
           ! first subface. This should be the same for all of them.
           ! Initialize the rotation matrix to the identity and the
           ! translation vector to zero.

           ll        = perHalo(kk)%periodicSubfaces(1)
           cgnsBlock = periodicGlobal(ll)%cgnsBlock
           cgnsSub   = periodicGlobal(ll)%cgnsSubface

           rotCenter = cgnsDoms(cgnsBlock)%conn1to1(cgnsSub)%rotationCenter

           rotMat(1,1) =  one; rotMat(1,2) = zero; rotMat(1,3) = zero
           rotMat(2,1) = zero; rotMat(2,2) =  one; rotMat(2,3) = zero
           rotMat(3,1) = zero; rotMat(3,2) = zero; rotMat(3,3) =  one

           trans = zero

           ! Loop over the number of periodic subface for the total
           ! periodic transformation.

           subfaceLoop: do mm=1,perHalo(kk)%nPeriodicSubfaces

             ll        = perHalo(kk)%periodicSubfaces(1)
             cgnsBlock = periodicGlobal(ll)%cgnsBlock
             cgnsSub   = periodicGlobal(ll)%cgnsSubface

             ! Store the data for the transformation of this subface
             ! a bit easier.

             translation = cgnsDoms(cgnsBlock)%conn1to1(cgnsSub)%translation

             theta = cgnsDoms(cgnsBlock)%conn1to1(cgnsSub)%rotationAngles(1)
             phi   = cgnsDoms(cgnsBlock)%conn1to1(cgnsSub)%rotationAngles(2)
             psi   = cgnsDoms(cgnsBlock)%conn1to1(cgnsSub)%rotationAngles(3)

             ! Construct the rotation matrix for this subface. Actually
             ! the inverse (== transpose) is constructed, because the
             ! cgns data is for the transformation from the current face
             ! to the donor and here the inverse is needed. Note
             ! furthermore that the sequence of rotation is first
             ! rotation around the x-axis, followed by rotation around
             ! the y-axis and finally rotation around the z-axis.

             cosTheta = cos(theta); cosPhi = cos(phi); cosPsi = cos(psi)
             sinTheta = sin(theta); sinPhi = sin(phi); sinPsi = sin(psi)

             rotMatrix(1,1) =  cosPhi*cosPsi
             rotMatrix(1,2) =  cosPhi*sinPsi
             rotMatrix(1,3) = -sinPhi

             rotMatrix(2,1) = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi
             rotMatrix(2,2) = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi
             rotMatrix(2,3) = sinTheta*cosPhi

             rotMatrix(3,1) = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi
             rotMatrix(3,2) = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi
             rotMatrix(3,3) = cosTheta*cosPhi

             ! Create the product of rotMat and rotMatrix. Copy the
             ! result back into rotMat.

             tmpMat(:,1) = rotMat(:,1)*rotMatrix(1,1) &
                         + rotMat(:,2)*rotMatrix(2,1) &
                         + rotMat(:,3)*rotMatrix(3,1)
             tmpMat(:,2) = rotMat(:,1)*rotMatrix(1,2) &
                         + rotMat(:,2)*rotMatrix(2,2) &
                         + rotMat(:,3)*rotMatrix(3,2)
             tmpMat(:,3) = rotMat(:,1)*rotMatrix(1,3) &
                         + rotMat(:,2)*rotMatrix(2,3) &
                         + rotMat(:,3)*rotMatrix(3,3)

             rotMat = tmpMat

             ! Update the translation vector. As we need the inverse
             ! of the transformation translation must be premultiplied
             ! by the rotation matrix and negated.

             trans = trans - rotMatrix(:,1)*translation(1) &
                           - rotMatrix(:,2)*translation(2) &
                           - rotMatrix(:,3)*translation(3)

           enddo subfaceLoop

           ! Store the results in periodicData(nn).

           periodicData(nn)%rotMatrix   = rotMat
           periodicData(nn)%rotCenter   = rotCenter
           periodicData(nn)%translation = trans

         enddo periodicLoop

         end subroutine setPeriodicData

       end subroutine determinePeriodicData
