module pointMatchedCommPattern


contains

  subroutine determineCommPattern(level)
    !
    !       determineCommPattern determines the communication pattern
    !       for the indicated grid level from the given block distribution
    !       and corresponding halo info. Both the first and second level
    !       cell halo communication pattern as well as the first level
    !       nodal halo communication pattern is determined.
    !       A recursive algorithm is used. First the face halo's are
    !       determined and from those the indirect halo's can be obtained
    !       by looping over the level of indirectness.
    !       This routine controls the creation of the communication
    !       pattern and basically contains the function calls to the
    !       subtasks.
    !
    use block
    use communication
    use haloList
    use periodicInfo
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer               :: ierr
    integer(kind=intType) :: i

    ! Determine the number of periodic faces.

    call determinePeriodicFaces

    ! Determine the amount of 1st and 2nd level cell halo's and 1st
    ! level node halo's.

    call determineNumberOfHalos(level)

    ! Allocate the memory for the help variables needed to
    ! determine the communication pattern.

    call allocMemHaloList(level)

    ! Determine the face halo's for both the 1st level cell and the
    ! 1st level node halo's.

    call determineFaceHalos(level)

    ! Determine the indirect 1st level node halo's.

    call determineIndirectHalos(nNodeHalo1st, iinode1st,    &
         nodeHalo1st, transformNode, &
         nodeIndex,   1_intType,     &
         1_intType,   0_intType, level)

    ! Release the memory of transformNode and nodeIndex.

    do i=1,nDom
       deallocate(nodeIndex(i)%entryList, stat=ierr)
       if(ierr /= 0)                           &
            call terminate("determineCommPattern", &
            "Deallocation error for &
            &nodeIndex(i)%entryList")
    enddo

    deallocate(nodeIndex, transformNode, stat=ierr)
    if(ierr /= 0)                           &
         call terminate("determineCommPattern", &
         "Deallocation error for nodeIndex &
         &and transformNode")

    ! Determine the indirect 1st level cell halo's.

    call determineIndirectHalos(nCellHalo1st, iicell1st,    &
         cellHalo1st, transformCell, &
         cellIndex,   2_intType,     &
         1_intType,    0_intType, level)

    ! Initialize the 2nd level cell halo list. Basically the 1st
    ! level cell halo list is copied.

    call init2ndLevelCellHalos

    ! Determine the indirect 2nd level cell halo's. As the indirect
    ! 1st level halo's are already treated an offset of 1 is passed.

    call determineIndirectHalos(nCellHalo2nd, iicell2nd,    &
         cellHalo2nd, transformCell, &
         cellIndex,   2_intType,     &
         2_intType,   1_intType, level)

    ! Release the memory of transformCell and cellIndex.

    do i=1,nDom
       deallocate(cellIndex(i)%entryList, stat=ierr)
       if(ierr /= 0)                             &
            call terminate("determineCommPattern", &
            "Deallocation error for &
            &cellIndex(i)%entryList")
    enddo

    deallocate(cellIndex, transformCell, stat=ierr)
    if(ierr /= 0)                           &
         call terminate("determineCommPattern", &
         "Deallocation error for cellIndex &
         &and transformCell")

    ! Sort the three lists in increasing order.

    call qsortHaloListType(nodeHalo1st, nNodeHalo1st)
    call qsortHaloListType(cellHalo1st, nCellHalo1st)
    call qsortHaloListType(cellHalo2nd, nCellHalo2nd)

    ! Determine the final communication data structures to store the
    ! halo info.

    call finalCommStructures(nodeHalo1st, nNodeHalo1st,  &
         commPatternNode_1st(level), &
         internalNode_1st(level), 0_intType)

    call finalCommStructures(cellHalo1st, nCellHalo1st,  &
         commPatternCell_1st(level), &
         internalCell_1st(level), 0_intType)

    call finalCommStructures(cellHalo2nd, nCellHalo2nd,  &
         commPatternCell_2nd(level), &
         internalCell_2nd(level), 0_intType)

    ! Determine the transformation for periodic halo's.

    call determinePeriodicData(nodeHalo1st, nNodeHalo1st,  &
         commPatternNode_1st(level), &
         internalNode_1st(level))

    call determinePeriodicData(cellHalo1st, nCellHalo1st,  &
         commPatternCell_1st(level), &
         internalCell_1st(level))

    call determinePeriodicData(cellHalo2nd, nCellHalo2nd,  &
         commPatternCell_2nd(level), &
         internalCell_2nd(level))

    ! Deallocate the memory for the 3 halo lists.

    call deallocatePointersHaloList(nodeHalo1st, nNodeHalo1st)
    call deallocatePointersHaloList(cellHalo1st, nCellHalo1st)
    call deallocatePointersHaloList(cellHalo2nd, nCellHalo2nd)

    deallocate(nodeHalo1st, cellHalo1st, cellHalo2nd, stat=ierr)
    if(ierr /= 0)                            &
         call terminate("determineCommPattern", &
         "Deallocation error for nodeHalo1st, &
         &cellHalo1st and cellHalo2nd")

    ! Deallocate the memory of periodicGlobal.

    deallocate(periodicGlobal, stat=ierr)
    if(ierr /= 0)                            &
         call terminate("determineCommPattern", &
         "Deallocation error for periodicGlobal")

  end subroutine determineCommPattern

  !-----------------------------------------------------------------

  subroutine deallocatePointersHaloList(entityHalo, nHalo)
    !
    !       deallocatePointersHaloList deallocates the memory of the
    !       pointer variables of entityHalo.
    !
    use haloList
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in)               :: nHalo
    type(haloListType), dimension(*), intent(inout) :: entityHalo
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i

    ! Loop over the number of halo's and deallocate the memory
    ! of the pointer variables, if allocated.

    do i=1,nHalo
       if( associated(entityHalo(i)%interp) ) then
          deallocate(entityHalo(i)%interp, stat=ierr)
          if(ierr /= 0)                                  &
               call terminate("deallocatePointersHaloList", &
               "Deallocation failure for interp")
       endif

       if( associated(entityHalo(i)%periodicSubfaces) ) then
          deallocate(entityHalo(i)%periodicSubfaces, stat=ierr)
          if(ierr /= 0)                                  &
               call terminate("deallocatePointersHaloList", &
               "Deallocation failure for periodicSubfaces")
       endif
    enddo

  end subroutine deallocatePointersHaloList

  subroutine determineIndirectHalos(nHalo, iihalo, entityHalo,  &
       transform, entityIndex,     &
       start, nLevel, offset,      &
       gridLevel)
    !
    !       determineIndirectHalos determines the indirect halo's via a
    !       recursive algorithm.
    !       Step 1.
    !       =======
    !       Determine for every indirect halo the closest face halo.
    !       If several options exist choose the one that does not
    !       correspond to a boundary halo, if possible. If several
    !       non-boundary halo's exist, just pick one. If all the closest
    !       face halo's are boundary halo's then there is no corresponding
    !       halo and the state is determined by the boundary conditions
    !       and/or extrapolation.
    !       Store the direction from the face halo to the indirect halo.
    !       Step 2.
    !       =======
    !       Determine the level of indirectness of every indirect halo.
    !       This is the sum of the absolute values of the elements of the
    !       direction vector. For 1st level halo's the maximum level of
    !       of indirectness is 2; for 2nd level halo's it is 5. These
    !       numbers are for 3 space dimensions.
    !       Step 3.
    !       =======
    !       Loop over the number of indirect levels.
    !         For every halo of the current level of indirectness do:
    !           - apply the transformation matrix of its corresponding
    !             face halo to the direction vector.
    !           - start in the donor cell of the face halo and travel in
    !             the direction of the transformed direction vector.
    !           - you either end up in an internal cell/node or in a halo.
    !             Case internal: you're done. Internal cell/node is the
    !                            donor.
    !             Case halo: this is guarenteed to be a halo of at least
    !                        one level of indirectness less than the
    !                        current level. Thus the donor is known and
    !                        you're done too. It is possible that it is a
    !                        boundary halo, but this is allowed.
    !       End loop over the number of indirect levels.
    !
    use haloList
    use indirectHalo
    use communication
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nHalo, start, nLevel, offset
    integer(kind=intType), intent(in) :: gridLevel
    integer(kind=intType), intent(inout) :: iihalo

    integer(kind=intType), dimension(:,:), intent(in) :: transform

    type(haloListType), dimension(:), intent(inout) :: entityHalo
    type(indexListType), dimension(:), intent(inout) :: entityIndex
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i

    ! Determine for every indirect halo the closest direct halo.

    call closestDirectHalos(entityHalo, entityIndex, start, &
         nLevel, offset, gridLevel)

    ! Sort the indirect halo's.

    call qsortIndHaloType(indHalo, nIndHalo)

    ! Abbreviate the number of indirect levels a bit easier, allocate
    ! the memory for nHaloPerLev and nHaloPerProc, and determine
    ! the values of nHaloPerLev.

    nLevOfInd = indHalo(nIndHalo)%levOfInd
    allocate(nHaloPerLev(0:nlevOfInd), nHaloPerProc(0:nProc), &
         stat=ierr)
    if(ierr /= 0)                                          &
         call terminate("determineIndirectHalos",             &
         "Allocation error for nHaloPerLev and &
         &nHaloPerProc")

    nHaloPerLev = 0
    do i=1,nIndHalo
       nHaloPerLev(indHalo(i)%levOfInd) = &
            nHaloPerLev(indHalo(i)%levOfInd) + 1
    enddo

    ! Put nHaloPerLev in cumulative storage format.

    do i=1,nlevOfInd
       nHaloPerLev(i) = nHaloPerLev(i) + nHaloPerLev(i-1)
    enddo

    ! Loop over the number of levels of indirectness.

    do i=1,nLevOfInd

       ! Determine the halo's for this level of indirectness.

       call indirectHalosPerLevel(i, iihalo, entityHalo, transform, &
            entityIndex)

       ! Synchronize the processors to avoid possible problems.

       call mpi_barrier(ADflow_comm_world, ierr)
    enddo

    ! Release the memory of the module indirectHalo.

    deallocate(indHalo, nHaloPerLev, nHaloPerProc, stat=ierr)
    if(ierr /= 0) &
         call terminate("determineIndirectHalos",        &
         "Deallocation error for indHalo, &
         &nHaloPerLev and nHaloPerProc")

    ! Check in debug mode if iihalo equals nHalo, as it should be.

    if( debug ) then
       if(iihalo /= nHalo)                        &
            call terminate("determineIndirectHalos", &
            "iihalo differs from nHalo")
    endif

  end subroutine determineIndirectHalos


  subroutine determinePeriodicData(entityHalo,   nHalo, &
       externalComm, internalComm)
    !
    !       determinePeriodicData determines the periodic transformation
    !       for both the external and the internal communication patterns.
    !
    use constants
    use communication
    use haloList
    use periodicInfo
    use utils, only : terminate
    use sorting, only :qsortIntegers

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
      !         setPeriodicData stores the periodic transformations and the
      !         corresponding halo's to which it must be applied in
      !         nPeriodic and periodicData. These variables are part of
      !         either the internal or external communication pattern.
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


  subroutine allocMemHaloList(level)
    !
    !       allocMemHaloList allocates the memory for the variables
    !       needed to construct the communication lists and the periodic
    !       information. These variables are located in the module
    !       haloList. Only the 1st level halo variables are allocated here
    !       to avoid unnecessary memory usage. The 2nd level cell halo
    !       are allocated later on in init2ndLevelCellHalos.
    !
    use block
    use haloList
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i
    integer(kind=intType) :: ie, je, ke, ib, jb, kb

    ! Allocate the memory for the 1st level cell and node halo lists.

    allocate(cellHalo1st(nCellHalo1st), &
         nodeHalo1st(nNodeHalo1st), stat=ierr)
    if(ierr /= 0)                                               &
         call terminate("allocMemHaloList",                        &
         "Memory allocation failure for cellHalo1st &
         &and nodeHalo1st")

    ! Initialize the level of indirectness to 0. This will only be
    ! overwritten for indirect boundary halo's. Also initialize the
    ! number of periodid subfaces to 0.

    do i=1,nCellHalo1st
       cellHalo1st(i)%levOfInd          = 0
       cellHalo1st(i)%nPeriodicSubfaces = 0
       nullify(cellHalo1st(i)%periodicSubfaces)
       nullify(cellHalo1st(i)%interp)
    enddo

    do i=1,nNodeHalo1st
       nodeHalo1st(i)%levOfInd          = 0
       nodeHalo1st(i)%nPeriodicSubfaces = 0
       nullify(nodeHalo1st(i)%periodicSubfaces)
       nullify(nodeHalo1st(i)%interp)
    enddo

    ! Allocate the memory to store the short hand of the transformation
    ! matrix for both cell and nodal halo's. In principle this matrix
    ! is only needed for the face (i.e. direct) halo's. However the
    ! difference between the number of 1st level halo's and the cell
    ! halo's is not so large and at this time of the program you should
    ! not worry too much about memory, because the metrics as well as
    ! the solution variables have not been allocated yet.

    allocate(transformCell(nCellHalo1st,3), &
         transformNode(nNodeHalo1st,3), stat=ierr)
    if(ierr /= 0)                                                 &
         call terminate("allocMemHaloList"  ,                        &
         "Memory allocation failure for transformCell &
         &and transformNode")

    ! Allocate the memory for nodeIndex and cellIndex, which will
    ! store the indices per block in the lists above.

    allocate(nodeIndex(nDom), cellIndex(nDom), stat=ierr)
    if(ierr /= 0)                                             &
         call terminate("allocMemHaloList",                      &
         "Memory allocation failure for nodeIndex &
         &and cellIndex")

    ! Loop over the number of blocks to allocate and initialize
    ! the elements of nodeIndex and cellIndex.

    do i=1,nDom

       ! Store the upper indices in the allocation a bit easier.

       ie = flowDoms(i,level,1)%ie
       je = flowDoms(i,level,1)%je
       ke = flowDoms(i,level,1)%ke

       ib = flowDoms(i,level,1)%ib
       jb = flowDoms(i,level,1)%jb
       kb = flowDoms(i,level,1)%kb

       ! Allocate the memory for entryList.

       allocate(nodeIndex(i)%entryList(0:ie,0:je,0:ke), &
            cellIndex(i)%entryList(0:ib,0:jb,0:kb), stat=ierr)
       if(ierr /= 0)                        &
            call terminate("allocMemHaloList", &
            "Memory allocation failure for entryList")

       ! Initialize entryList to zero. This serves as a check later
       ! on. Cell halo's are uniquely defined via the 1 to 1 block
       ! connectivity, but for node halo's (on the boundary of a
       ! subface) several possibilities exist.

       nodeIndex(i)%entryList = 0
       cellIndex(i)%entryList = 0

    enddo

  end subroutine allocMemHaloList

  subroutine determineFaceHalos(level)
    !
    !       determineFaceHalos determines the 1st level direct cell and
    !       node halo's. Direct halo means that at least one of the
    !       neighboring cell/nodes belongs is owned by the block.
    !       Consequently the halo can be found using the 1 to 1 block
    !       connectivity.
    !
    use constants
    use blockPointers
    use cgnsGrid
    use haloList
    use periodicInfo
    use utils, only : delta, terminate, setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j, k
    integer(kind=intType) :: ii, jj, kk, ll, mm, nn
    integer(kind=intType) :: iH, jH, kH, iD, jD, kD
    integer(kind=intType) :: indexPeriodic

    integer(kind=intType), dimension(3) :: myOffset, donorOffset
    integer(kind=intType), dimension(3) :: step

    integer(kind=intType), dimension(3,3) :: trMat
    integer(kind=intType), dimension(3,2) :: myCellRange
    integer(kind=intType), dimension(3,2) :: myNodeRange
    integer(kind=intType), dimension(3,2) :: donorCellRange

    type(cgnsPeriodicType) :: key

    ! Initialize the counter variables for the 1st level halo's to 0.

    iicell1st = 0
    iinode1st = 0

    ! Determine the 1st level cell and node halo lists by looping over
    ! the boundary faces of the blocks stored on this processor.

    domains: do nn=1,nDom

       ! Set the pointers for this block on the given level.

       call setPointers(nn,level,1_intType)
       !
       !         Loop over the boundary halo's first. The reason is that a
       !         node could belong to both a boundary and an internal subface.
       !         By looping first over the boundaries, the internal subfaces
       !         overwrite earlier set values by the boundary subface, which
       !         is desirable.
       !
       bocos: do mm=1,nBocos

          ! Determine the cell and nodal range for the halo's of this
          ! subface as well as the direction normal to the subface.

          call haloRanges(mm, myOffset, myCellRange, myNodeRange, step)
          !
          !           First treat the nodes on the subface.
          !
          ! Loop over the nodes of the boundary subface and store
          ! the halo info. For the edges of the subface it is possible
          ! that the node is already stored. This must be checked;
          ! otherwise this node will occur more than once in the list,
          ! which is not correct.

          do k=myNodeRange(3,1),myNodeRange(3,2),step(3)
             do j=myNodeRange(2,1),myNodeRange(2,2),step(2)
                do i=myNodeRange(1,1),myNodeRange(1,2),step(1)

                   ! Determine the indices of my nodal halo node.

                   iH = i + myOffset(1)
                   jH = j + myOffset(2)
                   kH = k + myOffset(3)

                   ! Store the halo, if it has not been stored yet.

                   if(nodeIndex(nn)%entryList(iH,jH,kH) == 0) then

                      ! Update the counter iinode1st. Store it in ii for
                      ! convenience. Store the index in %entryList.

                      iinode1st = iinode1st +1
                      ii        = iinode1st

                      nodeIndex(nn)%entryList(iH,jH,kH) = ii

                      ! Store the info in nodeHalo1st. Note that
                      ! donorBlock contains the boundary condition, although
                      ! this is not really important for node halo's.
                      ! Furthermore donorProc is set to -1 to indicate a
                      ! boundary halo and the donor indices are set to i,j,k.
                      ! In this way the extrapolated coordinates can be
                      ! obtained easily later on.

                      nodeHalo1st(ii)%myBlock = nn
                      nodeHalo1st(ii)%myI     = iH
                      nodeHalo1st(ii)%myJ     = jH
                      nodeHalo1st(ii)%myK     = kH

                      nodeHalo1st(ii)%donorProc  = -1
                      nodeHalo1st(ii)%donorBlock = BCType(mm)

                      nodeHalo1st(ii)%dI = i
                      nodeHalo1st(ii)%dJ = j
                      nodeHalo1st(ii)%dK = k

                   endif

                enddo
             enddo
          enddo
          !
          !           The cell halo's belonging to this subface. Direct cell
          !           halo's are unique and therefore info cannot already be
          !           written earlier.
          !
          ! Loop over the halo cells located adjacent to the subface.

          do k=myCellRange(3,1),myCellRange(3,2),step(3)
             do j=myCellRange(2,1),myCellRange(2,2),step(2)
                do i=myCellRange(1,1),myCellRange(1,2),step(1)

                   ! Check in debug mode whether this halo is already
                   ! stored. This should not be the case.

                   if( debug ) then
                      if(cellIndex(nn)%entryList(i,j,k) /= 0) &
                           call terminate("determineFaceHalos",  &
                           "boundary cell halo already stored")
                   endif

                   ! Update the counter iicell1st and store its value a
                   ! bit easier in ii and set entryList accordingly.

                   iicell1st = iicell1st +1
                   ii        = iicell1st

                   cellIndex(nn)%entryList(i,j,k) = ii

                   ! Store the info in cellHalo1st. Note that donorBlock
                   ! contains the boundary condition and donorProc is set
                   ! to -1 to indicate a boundary halo. The donor indices are
                   ! set to the owned cell on the other side of the subface.

                   cellHalo1st(ii)%myBlock = nn
                   cellHalo1st(ii)%myI     = i
                   cellHalo1st(ii)%myJ     = j
                   cellHalo1st(ii)%myK     = k

                   cellHalo1st(ii)%donorProc  = -1
                   cellHalo1st(ii)%donorBlock = BCType(mm)

                   cellHalo1st(ii)%dI = i - myOffset(1)
                   cellHalo1st(ii)%dJ = j - myOffset(2)
                   cellHalo1st(ii)%dK = k - myOffset(3)

                enddo
             enddo
          enddo

       enddo bocos
       !
       !         Loop over the 1 to 1 block to block boundaries.
       !
       n1to1Loop: do ll=1,n1to1

          ! Store the correct index for this subface, i.e. add the
          ! offset from the boundary subfaces.

          mm = nBocos + ll

          ! Check if the original subface is a periodic subface.
          ! Subfaces created by internal block splitting are indicated
          ! by 0 and are certainly not periodic. This must be tested
          ! first to avoid array overflow.

          indexPeriodic = 0

          kk = cgnsSubface(mm)
          if(kk > 0) then
             if(cgnsDoms(nbkGlobal)%conn1to1(kk)%periodic) then

                ! Determine the corresponding index in periodicGlobal.

                key%cgnsBlock   = nbkGlobal
                key%cgnsSubface = kk

                indexPeriodic = bsearchCGNSPeriodicType(key, periodicGlobal)

                if( debug ) then
                   if(indexPeriodic == 0)                 &
                        call terminate("determineFaceHalos", &
                        "Entry not found in periodicGlobal")
                endif
             endif
          endif

          ! Determine the cell and nodal range for the halo's of this
          ! subface as well as the direction normal to the subface.

          call haloRanges(mm, myOffset, myCellRange, &
               myNodeRange, step)

          ! Determine the complete transformation matrix from the
          ! given shorthand.

          trMat(1,1) = sign(1_intType,l1(mm)) * delta(l1(mm),1_intType)
          trMat(2,1) = sign(1_intType,l1(mm)) * delta(l1(mm),2_intType)
          trMat(3,1) = sign(1_intType,l1(mm)) * delta(l1(mm),3_intType)

          trMat(1,2) = sign(1_intType,l2(mm)) * delta(l2(mm),1_intType)
          trMat(2,2) = sign(1_intType,l2(mm)) * delta(l2(mm),2_intType)
          trMat(3,2) = sign(1_intType,l2(mm)) * delta(l2(mm),3_intType)

          trMat(1,3) = sign(1_intType,l3(mm)) * delta(l3(mm),1_intType)
          trMat(2,3) = sign(1_intType,l3(mm)) * delta(l3(mm),2_intType)
          trMat(3,3) = sign(1_intType,l3(mm)) * delta(l3(mm),3_intType)

          ! Determine the offset of the donor block.

          donorOffset(1) = trMat(1,1)*myOffset(1) &
               + trMat(1,2)*myOffset(2) &
               + trMat(1,3)*myOffset(3)
          donorOffset(2) = trMat(2,1)*myOffset(1) &
               + trMat(2,2)*myOffset(2) &
               + trMat(2,3)*myOffset(3)
          donorOffset(3) = trMat(3,1)*myOffset(1) &
               + trMat(3,2)*myOffset(2) &
               + trMat(3,3)*myOffset(3)
          !
          !           First treat the nodes on the subface.
          !
          ! Loop over the nodal range for this subface.

          do k=myNodeRange(3,1),myNodeRange(3,2),step(3)
             do j=myNodeRange(2,1),myNodeRange(2,2),step(2)
                do i=myNodeRange(1,1),myNodeRange(1,2),step(1)

                   ! Determine the donor indices by applying the
                   ! transformation matrix to i,j,k and adding the
                   ! offset to obtain the halo.

                   ii = i - myNodeRange(1,1)
                   jj = j - myNodeRange(2,1)
                   kk = k - myNodeRange(3,1)

                   iD = donorOffset(1) + dinBeg(mm) &
                        + trMat(1,1)*ii + trMat(1,2)*jj + trMat(1,3)*kk
                   jD = donorOffset(2) + djnBeg(mm) &
                        + trMat(2,1)*ii + trMat(2,2)*jj + trMat(2,3)*kk
                   kD = donorOffset(3) + dknBeg(mm) &
                        + trMat(3,1)*ii + trMat(3,2)*jj + trMat(3,3)*kk

                   ! Determine the indices of my nodal halo node.

                   iH = i + myOffset(1)
                   jH = j + myOffset(2)
                   kH = k + myOffset(3)

                   ! It is possible that this halo is already stored,
                   ! either as a boundary or as an internal halo. In the
                   ! former case it should be overwritten; in the latter
                   ! this is not strictly necessary, but it does not hurt.
                   ! Therefore simply overwrite the old index. If the
                   ! halo has not been stored yet, update iinode1st.
                   ! The index to store the info will be ii.

                   if(nodeIndex(nn)%entryList(iH,jH,kH) == 0) then
                      iinode1st = iinode1st +1
                      ii        = iinode1st

                      nodeIndex(nn)%entryList(iH,jH,kH) = ii
                   else
                      ii = nodeIndex(nn)%entryList(iH,jH,kH)
                   endif

                   ! Store the info in the correct place in nodeHalo1st.

                   nodeHalo1st(ii)%myBlock = nn
                   nodeHalo1st(ii)%myI     = iH
                   nodeHalo1st(ii)%myJ     = jH
                   nodeHalo1st(ii)%myK     = kH

                   nodeHalo1st(ii)%donorProc  = neighProc(mm)
                   nodeHalo1st(ii)%donorBlock = neighBlock(mm)

                   nodeHalo1st(ii)%dI = iD
                   nodeHalo1st(ii)%dJ = jD
                   nodeHalo1st(ii)%dK = kD

                   ! Store the short hand of the transformation matrix
                   ! for this halo.

                   transformNode(ii,1) = l1(mm)
                   transformNode(ii,2) = l2(mm)
                   transformNode(ii,3) = l3(mm)

                   ! It is possible that ii is treated earlier and hence
                   ! periodic info may have been stored. Remove this.

                   if(nodeHalo1st(ii)%nPeriodicSubfaces > 0) then
                      deallocate(nodeHalo1st(ii)%periodicSubfaces, stat=ierr)
                      if(ierr /= 0)                              &
                           call terminate("determineFaceHalos",     &
                           "Deallocation failure for &
                           &periodicSubfaces")
                      nullify(nodeHalo1st(ii)%periodicSubfaces)
                      nodeHalo1st(ii)%nPeriodicSubfaces = 0
                   endif

                   ! If the subface is periodic store the periodic info.

                   if(indexPeriodic > 0) then
                      nodeHalo1st(ii)%nPeriodicSubfaces = 1
                      allocate(nodeHalo1st(ii)%periodicSubfaces(1), &
                           stat=ierr)
                      if(ierr /= 0)                                   &
                           call terminate("determineFaceHalos",          &
                           "Memory allocation failure for &
                           &periodicSubfaces")
                      nodeHalo1st(ii)%periodicSubfaces(1) = indexPeriodic
                   endif

                enddo
             enddo
          enddo
          !
          !           The cell halo's belonging to this subface. Direct cell
          !           halo's are unique and therefore info cannot already be
          !           written earlier.
          !
          ! First determine the cell range of the donor block on
          ! the subface. This equals the nodal range, except that 1 is
          ! added to the smallest index. As it is possible that the
          ! index is running negatively, this should be taken into account.

          donorCellRange(1,1) = dinBeg(mm)
          donorCellRange(2,1) = djnBeg(mm)
          donorCellRange(3,1) = dknBeg(mm)

          donorCellRange(1,2) = dinEnd(mm)
          donorCellRange(2,2) = djnEnd(mm)
          donorCellRange(3,2) = dknEnd(mm)

          ! The loop to add 1 to the lowest index and to correct the
          ! index corresponding to the face we are on.

          do i=1,3
             if(donorCellRange(i,1) == donorCellRange(i,2)) then

                ! If the face corresponds to a min face, indicated by
                ! donorCellRange(i,1) == 1 then 1 must be added;
                ! otherwise nothing needs to be done.

                if(donorCellRange(i,1) == 1) then
                   donorCellRange(i,1) = 2
                   donorCellRange(i,2) = 2
                endif

             else if(donorCellRange(i,1) > donorCellRange(i,2)) then
                donorCellRange(i,2) = donorCellRange(i,2) + 1
             else
                donorCellRange(i,1) = donorCellRange(i,1) + 1
             endif
          enddo

          ! Loop over the halo cells located adjacent to the subface.

          do k=myCellRange(3,1),myCellRange(3,2),step(3)
             do j=myCellRange(2,1),myCellRange(2,2),step(2)
                do i=myCellRange(1,1),myCellRange(1,2),step(1)

                   ! Check in debug mode whether this halo is already
                   ! stored. This should not be the case.

                   if( debug ) then
                      if(cellIndex(nn)%entryList(i,j,k) /= 0) &
                           call terminate("determineFaceHalos",  &
                           "internal cell halo already stored")
                   endif

                   ! Determine the indices of the donor point by applying
                   ! the transformation matrix to i,j,k.

                   ii = i - myCellRange(1,1)
                   jj = j - myCellRange(2,1)
                   kk = k - myCellRange(3,1)

                   iD = donorCellRange(1,1) &
                        + trMat(1,1)*ii + trMat(1,2)*jj + trMat(1,3)*kk
                   jD = donorCellRange(2,1) &
                        + trMat(2,1)*ii + trMat(2,2)*jj + trMat(2,3)*kk
                   kD = donorCellRange(3,1) &
                        + trMat(3,1)*ii + trMat(3,2)*jj + trMat(3,3)*kk

                   ! Update the counter iicell1st and store its value a
                   ! bit easier in ii and set entryList accordingly.

                   iicell1st = iicell1st +1
                   ii        = iicell1st

                   cellIndex(nn)%entryList(i,j,k) = ii

                   ! Store the info in the correct place in cellHalo1st.

                   cellHalo1st(ii)%myBlock = nn
                   cellHalo1st(ii)%myI     = i
                   cellHalo1st(ii)%myJ     = j
                   cellHalo1st(ii)%myK     = k

                   cellHalo1st(ii)%donorProc  = neighProc(mm)
                   cellHalo1st(ii)%donorBlock = neighBlock(mm)

                   cellHalo1st(ii)%dI = iD
                   cellHalo1st(ii)%dJ = jD
                   cellHalo1st(ii)%dK = kD

                   ! Store the short hand of the transformation matrix
                   ! for this halo.

                   transformCell(ii,1) = l1(mm)
                   transformCell(ii,2) = l2(mm)
                   transformCell(ii,3) = l3(mm)

                   ! If the subface is periodic store the periodic info.
                   ! Note that for the cells it is not needed to check
                   ! if a previous transformation was already stored,
                   ! because cell halo's are unique.

                   if(indexPeriodic > 0) then
                      cellHalo1st(ii)%nPeriodicSubfaces = 1
                      allocate(cellHalo1st(ii)%periodicSubfaces(1), &
                           stat=ierr)
                      if(ierr /= 0)                                   &
                           call terminate("determineFaceHalos",          &
                           "Memory allocation failure for &
                           &periodicSubfaces")
                      cellHalo1st(ii)%periodicSubfaces(1) = indexPeriodic
                   endif

                enddo
             enddo
          enddo

       enddo n1to1Loop

    enddo domains

  end subroutine determineFaceHalos

  !      ==================================================================

  subroutine haloRanges(mm, offset, cellRange, nodeRange, step)
    !
    !       haloRanges determines the cell and nodal ranges for the given
    !       subface as well as the direction normal to the subface,
    !       pointing outwards. In case of negative running indices of the
    !       subface, step is set to -1; otherwise it is 1.

    use constants
    use blockPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: mm
    integer(kind=intType), dimension(3), intent(out) :: offset, step
    integer(kind=intType), dimension(3,2), intent(out) :: cellRange
    integer(kind=intType), dimension(3,2), intent(out) :: nodeRange
    !
    !      Local variables.
    !
    integer(kind=intType) :: i
    integer(kind=intType) :: cellHaloInd, cellHaloID

    ! Determine the offset in i, j and k direction, depending on the
    ! faceID. This can be interpreted as the outward pointing normal
    ! in the index domain. Furthermore store cellHaloInd and
    ! cellHaloID. This info is needed to construct the cell halo
    ! range correctly.

    offset = 0

    select case (BCFaceID(mm))
    case (iMin)
       offset(1)   = -1
       cellHaloInd =  1
       cellHaloID  =  1
    case (iMax)
       offset(1)   =  1
       cellHaloInd =  1
       cellHaloID  =  ie
    case (jMin)
       offset(2)   = -1
       cellHaloInd =  2
       cellHaloID  =  1
    case (jMax)
       offset(2)   =  1
       cellHaloInd =  2
       cellHaloID  =  je
    case (kMin)
       offset(3)   = -1
       cellHaloInd =  3
       cellHaloID  =  1
    case (kMax)
       offset(3)   =  1
       cellHaloInd =  3
       cellHaloID  =  ke
    end select

    ! Copy the nodal range.

    nodeRange(1,1) = inBeg(mm)
    nodeRange(2,1) = jnBeg(mm)
    nodeRange(3,1) = knBeg(mm)

    nodeRange(1,2) = inEnd(mm)
    nodeRange(2,2) = jnEnd(mm)
    nodeRange(3,2) = knEnd(mm)

    ! Determine the cell range. The cell numbering of a block starts
    ! at index 2, i.e. 1 higher than the node numbering. Consequently
    ! 1 must be added to the smallest indices of the nodal range.
    ! Take negative running indices into account and set step
    ! accordingly.

    do i=1,3
       if(nodeRange(i,1) > nodeRange(i,2)) then

          ! Negative running index.

          step(i)        = -1
          cellRange(i,1) = nodeRange(i,1)
          cellRange(i,2) = nodeRange(i,2) + 1
       else

          ! Positive running index.

          step(i)        = 1
          cellRange(i,1) = nodeRange(i,1) + 1
          cellRange(i,2) = nodeRange(i,2)
       endif
    enddo

    ! Correct the cell range for the index corresponding to the face
    ! we are on.

    cellRange(cellHaloInd,1) = cellHaloID
    cellRange(cellHaloInd,2) = cellHaloID

  end subroutine haloRanges


  subroutine init2ndLevelCellHalos
    !
    !       init2ndLevelCellHalos initializes the 2nd level cell halo
    !       list. Basically the 1st level cell halo list is copied and the
    !       counter iicell2nd is set to nCellHalo1st. This means that
    !       the 2nd level cell halo's are appended to the first level
    !       halo's. They are stored in a separate list, because the
    !       communication pattern of the 2nd level halo's is separate from
    !       the 1st level halo's. Efficiency is the reason to do this; it
    !       is more efficient to send one big message than two smaller
    !       ones.
    !
    use haloList
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j, jj

    ! Allocate the memory for the 2nd level halo list.

    allocate(cellHalo2nd(nCellHalo2nd), stat=ierr)
    if(ierr /= 0)                             &
         call terminate("init2ndLevelCellHalos", &
         "Memory allocation failure for cellHalo2nd")

    ! Initialize iicell2nd to nCellHalo1st.

    iiCell2nd = nCellHalo1st

    ! Copy the information from the 1st level cell halo's into the
    ! second level cell halo list. Make sure to make a deep copy.

    do i=1,nCellHalo1st
       cellHalo2nd(i)%myBlock    = cellHalo1st(i)%myBlock
       cellHalo2nd(i)%myI        = cellHalo1st(i)%myI
       cellHalo2nd(i)%myJ        = cellHalo1st(i)%myJ
       cellHalo2nd(i)%myK        = cellHalo1st(i)%myK
       cellHalo2nd(i)%donorProc  = cellHalo1st(i)%donorProc
       cellHalo2nd(i)%donorBlock = cellHalo1st(i)%donorBlock
       cellHalo2nd(i)%dI         = cellHalo1st(i)%dI
       cellHalo2nd(i)%dJ         = cellHalo1st(i)%dJ
       cellHalo2nd(i)%dK         = cellHalo1st(i)%dK
       cellHalo2nd(i)%levOfInd   = cellHalo1st(i)%levOfInd

       nullify(cellHalo2nd(i)%interp)

       cellHalo2nd(i)%nPeriodicSubfaces = &
            cellHalo1st(i)%nPeriodicSubfaces

       if(cellHalo2nd(i)%nPeriodicSubfaces > 0) then
          jj = cellHalo2nd(i)%nPeriodicSubfaces
          allocate(cellHalo2nd(i)%periodicSubfaces(jj), stat=ierr)
          if(ierr /= 0)                             &
               call terminate("init2ndLevelCellHalos", &
               "Memory allocation failure for &
               &periodicSubfaces")
          do j=1,jj
             cellHalo2nd(i)%periodicSubfaces(j) = &
                  cellHalo1st(i)%periodicSubfaces(j)
          enddo
       else
          nullify(cellHalo2nd(i)%periodicSubfaces)
       endif
    enddo

    ! Initialize the level of indirectness for the rest of the list
    ! to 0 and initialize the periodic data to 0 as well.

    do i=(nCellHalo1st+1), nCellHalo2nd
       cellHalo2nd(i)%levOfInd          = 0
       cellHalo2nd(i)%nPeriodicSubfaces = 0
       nullify(cellHalo2nd(i)%periodicSubfaces)
       nullify(cellHalo2nd(i)%interp)
    enddo

  end subroutine init2ndLevelCellHalos

  subroutine qsortHaloListType(arr, nn)
    !
    !       qsortHaloListType sorts the given number of halo's in
    !       increasing order based on the <= operator for this derived
    !       data type.
    !
    use haloList
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn

    type(haloListType), dimension(*), intent(inout) :: arr
    !
    !      Local variables.
    !
    integer(kind=intType), parameter :: m = 7

    integer(kind=intType) :: nStack
    integer(kind=intType) :: i, j, k, r, l, jStack, ii

    integer :: ierr

    type(haloListType) :: a, tmp

    integer(kind=intType), allocatable, dimension(:) :: stack
    integer(kind=intType), allocatable, dimension(:) :: tmpStack

    ! Allocate the memory for stack.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                         &
         call terminate("qsortHaloListType", &
         "Memory allocation failure for stack")

    ! Initialize the variables that control the sorting.

    jStack = 0
    l      = 1
    r      = nn

    ! Start of the algorithm

    do

       ! Check for the size of the subarray.

       if((r-l) < m) then

          ! Perform insertion sort

          do j=l+1,r
             a = arr(j)
             do i=(j-1),l,-1
                if(arr(i) <= a) exit
                arr(i+1) = arr(i)
             enddo
             arr(i+1) = a
          enddo

          ! In case there are no more elements on the stack, exit from
          ! the outermost do-loop. Algorithm has finished.

          if(jStack == 0) exit

          ! Pop stack and begin a new round of partitioning.

          r = stack(jStack)
          l = stack(jStack-1)
          jStack = jStack - 2

       else

          ! Subarray is larger than the threshold for a linear sort.
          ! Choose median of left, center and right elements as
          ! partitioning element a.
          ! Also rearrange so that (l) <= (l+1) <= (r).

          k = (l+r)/2
          tmp      = arr(k)             ! Swap the elements
          arr(k)   = arr(l+1)           ! k and l+1.
          arr(l+1) = tmp

          if(arr(r) < arr(l)) then
             tmp    = arr(l)             ! Swap the elements
             arr(l) = arr(r)             ! r and l.
             arr(r) = tmp
          endif

          if(arr(r) < arr(l+1)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(r)           ! r and l+1.
             arr(r)   = tmp
          endif

          if(arr(l+1) < arr(l)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(l)           ! l and l+1.
             arr(l)   = tmp
          endif

          ! Initialize the pointers for partitioning.

          i = l+1
          j = r
          a = arr(l+1)

          ! The innermost loop

          do

             ! Scan up to find element >= a.
             do
                i = i+1
                if(a <= arr(i)) exit
             enddo

             ! Scan down to find element <= a.
             do
                j = j-1
                if(arr(j) <= a) exit
             enddo

             ! Exit the loop in case the pointers i and j crossed.

             if(j < i) exit

             ! Swap the element i and j.

             tmp    = arr(i)
             arr(i) = arr(j)
             arr(j) = tmp
          enddo

          ! Swap the entries j and l+1. Remember that a equals
          ! arr(l+1).

          arr(l+1) = arr(j)
          arr(j)   = a

          ! Push pointers to larger subarray on stack,
          ! process smaller subarray immediately.

          jStack = jStack + 2
          if(jStack > nStack) then

             ! Storage of the stack is too small. Reallocate.

             allocate(tmpStack(nStack), stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortHaloListType", &
                  "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortHaloListType", &
                  "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortHaloListType", &
                  "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortHaloListType", &
                  "Deallocation error for tmpStack")
          endif

          if((r-i+1) >= (j-l)) then
             stack(jStack)   = r
             r               = j-1
             stack(jStack-1) = j
          else
             stack(jStack)   = j-1
             stack(jStack-1) = l
             l               = j
          endif

       endif
    enddo

    ! Release the memory of stack.

    deallocate(stack, stat=ierr)
    if(ierr /= 0)                         &
         call terminate("qsortHaloListType", &
         "Deallocation error for stack")

    ! Check in debug mode whether the array is really sorted.

    if( debug ) then
       do i=1,(nn-1)
          if(arr(i+1) < arr(i))                 &
               call terminate("qsortHaloListType", &
               "Array is not sorted correctly")
       enddo
    endif

  end subroutine qsortHaloListType

  subroutine finalCommStructures(entityHalo, nHalo, commPattern, &
       internalComm, nInterp)
    !
    !       FinalCommStructures determines the communication data
    !       structures used in the flow solver, commPattern and
    !       internalComm, from the given haloList, entityHalo.
    !
    use communication
    use haloList
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in)            :: nHalo
    type(haloListType), dimension(*), intent(in) :: entityHalo

    type(commType),         intent(out) :: commPattern
    type(internalCommType), intent(out) :: internalComm

    integer(kind=intType), intent(in) :: nInterp
    !
    !      Local variables.
    !
    integer :: ierr, count, dest, source, sizeRecv

    integer, dimension(mpi_status_size)  :: mpiStatus

    integer(kind=intType) :: i, j
    integer(kind=intType) :: ii, jj, kk, ll, mm, nn, pp
    integer(kind=intType) :: sizeBuffer

    integer(kind=intType), dimension(:), allocatable :: nHaloPerProc
    integer(kind=intType), dimension(:), allocatable :: sendInfo

    integer(kind=intType), dimension(:,:), allocatable :: buffer
    integer(kind=intType), dimension(:,:), allocatable :: recvBuf

    real(kind=realType), dimension(:,:), allocatable :: bufInt
    real(kind=realType), dimension(:,:), allocatable :: recvBufInt

    ! Allocate the memory for nHaloPerProc and initialize its
    ! values to 0. Boundary halo's are stored with a processor id -1
    ! and the array will be put in cumulative storage format later,
    ! which explains the boundaries in the allocation.
    ! Also allocate the memory for sendInfo, which is needed in the
    ! all to all communication call.

    allocate(nHaloPerProc(0:nProc), sendInfo(0:nProc-1), stat=ierr)
    if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
         "Memory allocation error for nHaloPerProc &
         &and sendInfo")

    nHaloPerProc = 0

    ! Loop over the number of halo's and determine the nHaloPerProc.
    ! Note that the boundary conditions are stored in entityHalo with
    ! a processor id -1.

    do i=1,nHalo
       ii = entityHalo(i)%donorProc + 1
       nHaloPerProc(ii) = nHaloPerProc(ii) + 1
    enddo

    ! Perform an all to all communication, such that the processors
    ! know how many messages will be received as well as their size.

    call mpi_alltoall(nHaloPerProc(1), 1,  adflow_integer, &
         sendInfo(0),      1, adflow_integer, &
         ADflow_comm_world, ierr)

    ! Allocate the memory for indexRecvProc, the index in the
    ! receive info for a certain processor. Initialize these value
    ! to 0, indicating that nothing is received from a processor.

    allocate(commPattern%indexRecvProc(0:nProc-1), stat=ierr)
    if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
         "Memory allocation error for indexRecvProc")

    do i=0,(nProc-1)
       commPattern%indexRecvProc(i) = 0
    enddo

    ! Determine the number of processors from which i receive data.
    ! Receive data here means when an actual exchange takes place.
    ! Furthermore determine the size of the buffer for the
    ! nonblocking sends.

    commPattern%nProcRecv = 0
    sizeBuffer             = 0
    do i=0,(nProc-1)
       ii = i+1
       if(nHaloPerProc(ii) > 0 .and. i /= myId) then
          commPattern%nProcRecv         = commPattern%nProcRecv + 1
          commPattern%indexRecvProc(i) = commPattern%nProcRecv

          sizeBuffer = sizeBuffer + nHaloPerProc(ii)
       endif
    enddo

    ! Allocate the memory for recvProc, nrecv, nrecvCum
    ! and recvList.

    ii = commPattern%nProcRecv
    allocate(commPattern%recvProc(ii),   &
         commPattern%nrecv(ii),       &
         commPattern%nrecvCum(0:ii), &
         commPattern%recvList(ii), stat=ierr)
    if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
         "Memory allocation error for recvProc, etc")

    ! Allocate memory for buffers, needed for the nonblocking sends.

    allocate(buffer(4,sizeBuffer), bufInt(ninterp,sizeBuffer), &
         stat=ierr)
    if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
         "Memory allocation error for buffer, bufInt.")

    ! Repeat the loop over the processors, but now store receive info.
    ! Jj is the counter for the current index in the receive processors
    ! and mm the counter in the buffer used to send the halo info.

    jj = 0
    mm = 1
    commPattern%nrecvCum(0) = 0
    do i=0,(nProc-1)
       ii = i+1
       if(nHaloPerProc(ii) > 0 .and. i /= myId) then

          ! Update the counter jj, store the number of halo's a bit
          ! easier in kk and set recvProc and nrecv.

          jj = jj + 1
          kk = nHaloPerProc(ii)

          commPattern%recvProc(jj) = i
          commPattern%nrecv(jj)     = kk
          commPattern%nrecvCum(jj) = commPattern%nrecvCum(jj-1) + kk

          ! Allocate the memory for the receive list.

          allocate(commPattern%recvList(jj)%block(kk),     &
               commPattern%recvList(jj)%indices(kk,3), &
               stat=ierr)
          if(ierr /= 0)                             &
               call terminate("finalCommStructures", &
               "Memory allocation failure for block and &
               &indices of recvlist")

          ! Copy the halo info in the receive list of commPattern.

          do j=1,kk
             ll = j + nHaloPerProc(i)

             commPattern%recvList(jj)%block(j) = entityHalo(ll)%myBlock

             commPattern%recvList(jj)%indices(j,1) = entityHalo(ll)%myI
             commPattern%recvList(jj)%indices(j,2) = entityHalo(ll)%myJ
             commPattern%recvList(jj)%indices(j,3) = entityHalo(ll)%myK
          enddo

          ! Copy the donor info in the buffer and send it to the
          ! appropriate processor.

          nn = mm - 1
          do j=1,kk
             ll = j + nHaloPerProc(i)
             nn = nn +1

             buffer(1,nn) = entityHalo(ll)%donorBlock
             buffer(2,nn) = entityHalo(ll)%dI
             buffer(3,nn) = entityHalo(ll)%dJ
             buffer(4,nn) = entityHalo(ll)%dK

             do pp = 1,ninterp
                bufInt(pp,nn) = entityHalo(ll)%interp(pp)
             end do
          enddo

          ! Copy some values to be sure integer data is used in the
          ! MPI call.

          count = 4*kk
          dest  = i

          ! And send the stuff.

          call mpi_isend(buffer(1,mm), count, adflow_integer, &
               dest, dest+2, ADflow_comm_world,     &
               sendRequests(jj), ierr)

          ! Now send the interpolants, if any are being exchanged.

          if(ninterp > 0) then
             count = nInterp*kk
             call mpi_isend(bufInt(1,mm), count, adflow_real, &
                  dest, dest+3, ADflow_comm_world,  &
                  recvRequests(jj), ierr)
          end if

          ! Update mm to the index in buffer for the next processor.

          mm = mm + kk
       endif

       ! Put nHaloPerProc in cumulative storage format.

       nHaloPerProc(ii) = nHaloPerProc(ii) + nHaloPerProc(i)
    enddo

    ! Determine the number of internal memory to memory copies and
    ! allocate the memory for it.

    ii = nHaloPerProc(myId+1) - nHaloPerProc(myId)
    internalComm%ncopy = ii

    allocate(internalComm%donorBlock(ii),          &
         internalComm%donorIndices(ii,3),      &
         internalComm%donorInterp(ii,ninterp), &
         internalComm%haloBlock(ii),           &
         internalComm%haloIndices(ii,3), stat=ierr)
    if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
         "Memory allocation failure for internalComm")

    ! Copy the info from the halo list.

    do i=1,internalComm%ncopy
       ii = i + nHaloPerProc(myId)

       internalComm%donorBlock(i) = entityHalo(ii)%donorBlock
       internalComm%haloBlock(i)  = entityHalo(ii)%myBlock

       internalComm%donorIndices(i,1) = entityHalo(ii)%dI
       internalComm%donorIndices(i,2) = entityHalo(ii)%dJ
       internalComm%donorIndices(i,3) = entityHalo(ii)%dK

       do pp = 1,ninterp
          internalComm%donorInterp(i,pp) = entityHalo(ii)%interp(pp)
       end do

       internalComm%haloIndices(i,1) = entityHalo(ii)%myI
       internalComm%haloIndices(i,2) = entityHalo(ii)%myJ
       internalComm%haloIndices(i,3) = entityHalo(ii)%myK
    enddo

    ! Determine from the earlier MPI_alltoall call the processors
    ! to whom i must send data in a normal data exchange. Also
    ! determine the amount i must send.

    ! First determine the number of processors to which i must send
    ! data and the size of the buffer for receiving data.
    ! Allocate and determine indexSendProc as well.

    allocate(commPattern%indexSendProc(0:nProc-1), stat=ierr)
    if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
         "Memory allocation error for indexSendProc")

    do i=0,(nProc-1)
       commPattern%indexSendProc(i) = 0
    enddo

    commPattern%nProcSend = 0
    sizeBuffer = 0

    do i=0,(nProc-1)
       if(sendInfo(i) > 0 .and. i /= myId) then
          commPattern%nProcSend         = commPattern%nProcSend + 1
          commPattern%indexSendProc(i) = commPattern%nProcSend

          sizeBuffer = max(sizeBuffer,sendInfo(i))
       endif
    enddo

    ! Allocate the memory for sendProc, nsend, nsendCum
    ! and sendList

    ii = commPattern%nProcSend
    allocate(commPattern%sendProc(ii),   &
         commPattern%nsend(ii),       &
         commPattern%nsendCum(0:ii), &
         commPattern%sendList(ii),   stat=ierr)
    if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
         "Memory allocation error for sendProc, etc")

    ! Repeat the loop over the number of processors, but now store
    ! the sending info. Use ii as a counter.

    ii = 0
    commPattern%nsendCum(0) = 0
    do i=0,(nProc-1)
       if(sendInfo(i) > 0 .and. i /= myId) then
          ii = ii + 1

          commPattern%sendProc(ii) = i
          commPattern%nsend(ii)     = sendInfo(i)
          commPattern%nsendCum(ii) = commPattern%nsendCum(ii-1) &
               + sendInfo(i)
       endif
    enddo

    ! Allocate the memory for the send lists.

    do i=1,commPattern%nProcSend
       ii = commPattern%nsend(i)
       allocate(commPattern%sendList(i)%block(ii),          &
            commPattern%sendList(i)%indices(ii,3),      &
            commPattern%sendList(i)%interp(ii,ninterp), &
            stat=ierr)
       if(ierr /= 0)                             &
            call terminate("finalCommStructures", &
            "Memory allocation failure for block, interp, &
            &and indices of sendlist")
    enddo

    ! Allocate the memory for the receiving buffers.

    allocate(recvBuf(4,sizeBuffer), &
         recvBufInt(ninterp,sizeBuffer), stat=ierr)
    if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
         "Memory allocation failure for recvBuffers")

    ! Loop over the number of processors to receive my send info.

    recvSendInfo: do i=1,commPattern%nProcSend

       ! Block until a message arrives.

       call mpi_probe(mpi_any_source, myId+2, ADflow_comm_world, &
            mpiStatus, ierr)

       ! Store the source processor a bit easier and determine the
       ! index in the send data structure.

       source = mpiStatus(mpi_source)
       ii     = commPattern%indexSendProc(source)

       ! Perform some checks in debug mode.

       if( debug ) then
          if(ii <= 0 .or. ii > commPattern%nProcSend) &
               call terminate("finalCommStructures",     &
               "Send processor not in the list")

          call mpi_get_count(mpiStatus, adflow_integer, sizeRecv, ierr)
          if(sizeRecv /= 4*commPattern%nsend(ii)) &
               call terminate("finalCommStructures",     &
               "Unexpected size of message")
       endif

       ! Receive the message. As it has already arrived a blocking
       ! receive can be used.

       sizeRecv = 4*commPattern%nsend(ii)

       call mpi_recv(recvBuf, sizeRecv, adflow_integer, source, &
            myId+2, ADflow_comm_world, mpiStatus, ierr)

       ! Now receive the interpolants, if any.

       if(ninterp > 0) then
          sizeRecv = nInterp*commPattern%nsend(ii)
          call mpi_recv(recvBufInt, sizeRecv, adflow_real, source, &
               myId+3, ADflow_comm_world, mpiStatus, ierr)
       end if

       ! Store the info I must send to this processor in a normal
       ! data exchange.

       do j=1,commPattern%nsend(ii)
          commPattern%sendList(ii)%block(j) = recvBuf(1,j)

          commPattern%sendList(ii)%indices(j,1) = recvBuf(2,j)
          commPattern%sendList(ii)%indices(j,2) = recvBuf(3,j)
          commPattern%sendList(ii)%indices(j,3) = recvBuf(4,j)

          do pp = 1,ninterp
             commPattern%sendList(ii)%interp(j,pp) = recvBufInt(pp,j)
          end do
       enddo

    enddo recvSendInfo

    ! Complete the nonblocking sends. Dest is used as a temporary
    ! storage such that an integer is passed to the MPI call

    dest = commPattern%nProcRecv
    do i=1,commPattern%nProcRecv
       call mpi_waitany(dest, sendRequests, count, mpiStatus, ierr)
    enddo

    ! Repeat the call if any interpolants were exchanged.

    if(ninterp > 0) then
       do i=1,commPattern%nProcRecv
          call mpi_waitany(dest, recvRequests, count, mpiStatus, ierr)
       enddo
    end if

    ! Release the memory of the help variables allocated in this
    ! subroutine.

    deallocate(nHaloPerProc, sendInfo, buffer, recvBuf, &
         bufInt, recvBufInt, stat=ierr)
    if(ierr /= 0) &
         call terminate("finalCommStructures", &
         "Deallocation error for nHaloPerProc, &
         &sendInfo, and temporary buffers")

  end subroutine finalCommStructures

  subroutine closestDirectHalos(entityHalo, entityIndex, &
       start, nLevel, offset, gridLevel)
    !
    !       closestDirectHalos determines the number of indirect halo's
    !       to be treated and its corresponding direct halo.
    !
    use block
    use bcHalo
    use haloList
    use indirectHalo
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: start, nLevel, offset
    integer(kind=intType), intent(in) :: gridLevel

    type(haloListType), dimension(:), intent(in) :: entityHalo
    type(indexListType), dimension(:), intent(in) :: entityIndex
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j, k, jj, nn, mm
    integer(kind=intType) :: il, jl, kl
    integer(kind=intType) :: iStart, jStart, kStart, iEnd, jEnd, kEnd
    integer(kind=intType) :: nHaloDirI, nHaloDirJ, nHaloDirK
    integer(kind=intType) :: iindHalo, levOfInd

    integer(kind=intType), dimension(3) :: dir, haloDir, ii

    type(bcHaloType), dimension(3) :: bcHalos
    !

    ! Determine the number of indirect halo's for which the donor must
    ! be determined. Initialize iindHalo, the actual counter, to 0.

    nIndHalo = getNumberIndirectHalos(start, nLevel, offset, &
         gridLevel)
    iindHalo = 0

    ! Allocate the memory for indHalo

    allocate(indHalo(nIndHalo), stat=ierr)
    if(ierr /= 0)                         &
         call terminate("closestDirectHalos", &
         "Memory allocation failure for indHalo")

    ! Determine the lower bound for the block with halo's.
    ! This is identical for all the blocks.

    iStart = start - nLevel
    jStart = iStart
    kStart = iStart

    ! Loop over the blocks on this processor.

    domains: do nn=1,nDom

       ! Store the number of nodes in the three directions of the
       ! current block a bit easier.

       il = flowDoms(nn,gridLevel,1)%il
       jl = flowDoms(nn,gridLevel,1)%jl
       kl = flowDoms(nn,gridLevel,1)%kl

       ! Determine the upper boundaries of this block with halo's.

       iEnd = il + nLevel
       jEnd = jl + nLevel
       kEnd = kl + nLevel

       ! Loop over all the entities of the full halo block.

       iLoop: do i=iStart,iEnd

          ! Determine the halo in i-direction.

          dir(1) = min(i-start,max(i-il,0_intType))

          if(dir(1) == 0) then
             nHaloDirI = 0
          else
             nHaloDirI = 1
             haloDir(1) = 1
          endif

          jLoop: do j=jStart,jEnd

             ! Determine the halo in j-direction.

             dir(2) = min(j-start,max(j-jl,0_intType))

             if(dir(2) == 0) then
                nHaloDirJ = nHaloDirI
             else
                nHaloDirJ          = nHaloDirI + 1
                haloDir(nHaloDirJ) = 2
             endif

             kLoop: do k=kStart,kEnd

                ! Determine the halo in k-direction.

                dir(3) = min(k-start,max(k-kl,0_intType))

                if(dir(3) == 0) then
                   nHaloDirK = nHaloDirJ
                else
                   nHaloDirK          = nHaloDirJ + 1
                   haloDir(nHaloDirK) = 3
                endif

                ! Determine the level of indirectness relative to the
                ! block boundary. This is the sum of the absolute values
                ! of dir.

                levOfInd = abs(dir(1)) + abs(dir(2)) + abs(dir(3))

                ! If levOfInd <= 1, this means that we are either dealing
                ! with an owned entity of the block or a direct halo.
                ! As this is not an indirect halo continue with the next
                ! entity.

                if(levOfInd <= 1) cycle

                ! A distinction must now be made between 1st level halo's
                ! and higher level halo's. In this code it is done such
                ! that all first level halo's are determined first,
                ! followed by the higher level halo's. The reasons is that
                ! the communication pattern of the first level halo's is
                ! stored separately. The distinction between both cases is
                ! the value of offset. Offset == 0 means 1st level halo's;
                ! offset > 0 higher level halo's. Anyway, the level of
                ! indirect halo's for 1st level halo's is either 2 or 3.
                ! They can be distinguished by higher level halo's by the
                ! fact that their number of halo directions equal their
                ! level of indirectness. So continue with the next halo
                ! if we are dealing with such a case here.

                if(offset > 0 .and. nHaloDirK == levOfInd) cycle

                ! This indirect halo should be stored in the indHalo.
                ! Independent of the choice of the corresponding direct
                ! halo quite a bit of info can already be set. Update the
                ! counter iindHalo and do this.
                ! Note that 1 is substracted from levOfInd, because
                ! in the other routines levOfInd is the level of
                ! indirectness of the halo's and not the distance to
                ! the nearest owned entity. The difference is 1.

                iindHalo = iindHalo +1

                indHalo(iindHalo)%myBlock = nn

                indHalo(iindHalo)%myI = i
                indHalo(iindHalo)%myJ = j
                indHalo(iindHalo)%myK = k

                indHalo(iindHalo)%levOfInd = levOfInd -1

                ! The number of possibilities for the corresponding direct
                ! halo are nHaloDirK. Therefore loop over the number of
                ! possibilities and take an internal halo if possible. If
                ! several options exist, just take one. If all nearest
                ! direct halo's are boundary halo's then this indirect halo
                ! will also be a boundary halo.

                do mm=1,nHaloDirK

                   ! Determine the corresponding direct halo. Note that in
                   ! dir the direction from the nearest owned entity is stored.
                   ! Consequently 1 must be added/substracted from this
                   ! direction in the coordinate direction specified by
                   ! haloDir(mm) to obtain the direct halo.

                   ii(1) = i - dir(1)
                   ii(2) = j - dir(2)
                   ii(3) = k - dir(3)

                   ii(haloDir(mm)) = ii(haloDir(mm)) &
                        + sign(1_intType, dir(haloDir(mm)))

                   ! Store the index of the direct halo in entityHalo
                   ! a bit easier.

                   jj = entityIndex(nn)%entryList(ii(1),ii(2),ii(3))

                   ! Copy the data in bcHalos for possible later use.
                   ! Remember that donorBlock contains the boundary
                   ! condition in case jj is a boundary halo.

                   bcHalos(mm)%directHalo = jj
                   bcHalos(mm)%bc         = entityHalo(jj)%donorBlock

                   ! Check if jj > 0. If not this means that the halo information
                   ! in the cgns file is not correct.

                   if(jj == 0)                            &
                        call terminate("closestDirectHalos", &
                        "Closest direct halo not in halo &
                        &list. Something wrong with BC info?")

                   if(entityHalo(jj)%donorProc >= 0) then

                      ! Face halo is an internal block boundary halo. Store
                      ! the rest of the info in indHalo.

                      indHalo(iindHalo)%myDirectHalo = jj

                      indHalo(iindHalo)%donorProc = entityHalo(jj)%donorProc

                      ! Direct halo has been set. Exit the loop.

                      exit

                   endif
                enddo

                ! Check if a halo has been set. If not this is a boundary
                ! halo. Set the direct halo with the most important
                ! boundary condition.

                if(mm > nHaloDirK) then

                   ! Sort the boundary halo's in increasing order. The
                   ! most important will be in position nHaloDirK after
                   ! the sorting.

                   call sortBCHaloType(bcHalos, nHaloDirK)

                   ! Set the value of myDirectHalo and set donorProc
                   ! to -1 to indicate a boundary halo.

                   indHalo(iindHalo)%myDirectHalo = &
                        bcHalos(nHaloDirK)%directHalo
                   indHalo(iindHalo)%donorProc    = -1
                endif

             enddo kLoop
          enddo jLoop
       enddo iLoop

    enddo domains

  end subroutine closestDirectHalos


  function getNumberIndirectHalos(start, nLevel, offset, gridLevel)
    !
    !       getNumberIndirectHalos determines the number of indirect
    !       halo's for which the donor must be determined.
    !
    use block
    implicit none
    !
    !      Function type
    !
    integer(kind=intType) :: getNumberIndirectHalos
    !
    !      Function arguments
    !
    integer(kind=intType), intent(in) :: start, nLevel, offset
    integer(kind=intType), intent(in) :: gridLevel
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, il, jl, kl
    integer(kind=intType) :: i, j, k, ii, jj, kk

    ! Initialize getNumberIndirectHalos and loop over the blocks.

    getNumberIndirectHalos = 0

    domains: do nn=1,nDom

       ! Store the number of nodes in the three directions of the
       ! current block a bit easier.

       il = flowDoms(nn,gridLevel,1)%il
       jl = flowDoms(nn,gridLevel,1)%jl
       kl = flowDoms(nn,gridLevel,1)%kl

       ! Determine the size in every coordinate direction with and
       ! without halo's

       i = il - start + 2*offset + 1
       j = jl - start + 2*offset + 1
       k = kl - start + 2*offset + 1

       ii = il - start + 2*nLevel + 1
       jj = jl - start + 2*nLevel + 1
       kk = kl - start + 2*nLevel + 1

       ! Add the total number of halo's to getNumberIndirectHalos.

       getNumberIndirectHalos = getNumberIndirectHalos &
            + ii*jj*kk - i*j*k

       ! Substract the direct halo's. The convention is such that in
       ! a first loop the first level halo's are determined,
       ! offset == 0, and the higher level halo's in a next loop,
       ! offset > 0. Only in the former case the direct halo's must
       ! be sustracted.

       if(offset == 0) then

          ! Store the number of owned number of entities in i, j, k

          i = il - start + 1
          j = jl - start + 1
          k = kl - start + 1

          ! Substract the direct halo's

          getNumberIndirectHalos = getNumberIndirectHalos &
               - 2*(i*j + i*k + j*k)
       endif

    enddo domains

  end function getNumberIndirectHalos

  subroutine determineNumberOfHalos(level)
    !
    !       determineNumberOfHalos determines the amount of 1st and 2nd
    !       level cell halo's as well as the number of 1st level node
    !       halo's stored on this processor.
    !
    use blockPointers
    use haloList
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, nn

    ! Initialize the amount of halo cells and nodes to 0.

    nCellHalo1st = 0
    nCellHalo2nd = 0
    nNodeHalo1st = 0

    ! Loop over the number of blocks stored on this processor.

    do i=1,nDom

       ! Set the pointers for this block. The halo construction is
       ! the same for all time spectral solutions, so only the 1st
       ! needs to be considered.

       call setPointers(i,level,1_intType)

       ! Determine the number of 1st level halo cells for this block
       ! and add it to nCellHalo1st. Note the ie == nx + 2, etc.

       nn = nx*ny*nz

       nCellHalo1st = nCellHalo1st - nn + ie*je*ke

       ! Idem for the second level halo's. However there is no variable
       ! which stores nx + 4, so it is computed.

       nCellHalo2nd = nCellHalo2nd - nn + (nx+4)*(ny+4)*(nz+4)

       ! Idem for the 1st level node halo's. Use is made of the fact
       ! that ib == il + 2, etc.

       nNodeHalo1st = nNodeHalo1st + ib*jb*kb - il*jl*kl
    enddo

  end subroutine determineNumberOfHalos

  subroutine indirectHalosPerLevel(level, iihalo, entityHalo, &
       transform, entityIndex)
    !
    !       indirectHalosPerLevel determines the donor cells for the
    !       halo's of the given level of indirectness. From the known
    !       appropriate direct halo and its donor, the corresponding cell
    !       in the donor block is determined.
    !
    use haloList
    use indirectHalo
    use communication
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in)    :: level
    integer(kind=intType), intent(inout) :: iihalo

    integer(kind=intType), dimension(:,:), intent(in) :: transform

    type(haloListType),  dimension(:), intent(inout) :: entityHalo
    type(indexListType), dimension(:), intent(inout) :: entityIndex
    !
    !      Local variables.
    !
    integer :: ierr
    integer :: count, dest, source, sizeRecv

    integer, dimension(mpi_status_size)  :: mpiStatus

    integer, allocatable, dimension(:) :: sizeMessage

    integer(kind=intType) :: i, ii, nn, mm, ms
    integer(kind=intType) :: l1, l2, l3, proc
    integer(kind=intType) :: start, eend, nProcsSend, nProcsRecv
    integer(kind=intType) :: sizeSendBuf, sizeRecvBuf
    integer(kind=intType) :: nLocalHalos
    integer(kind=intType) :: nItemReturn, nItemSend, nItemAlloc

    integer(kind=intType), dimension(2) :: tmpBuf

    integer(kind=intType), allocatable, dimension(:) :: counter
    integer(kind=intType), allocatable, dimension(:) :: sendBuf
    integer(kind=intType), allocatable, dimension(:) :: recvBuf
    integer(kind=intType), allocatable, dimension(:) :: localHalos

    ! Determine the number of items per entity that is returned from
    ! the requested processor. This is 5 + level, because of the
    ! storage of possible periodic transformations.
    ! Also set nItemSend to 7 and determine the number of items that
    ! should be allocated in the communication buffers.

    nItemReturn = 5 + level
    nItemSend   = 7
    nItemAlloc  = max(nItemReturn, nItemSend)

    ! Abbreviate the start and ending index for this level
    ! a bit easier.

    start = nHaloPerLev(level-1) +1
    eend = nHaloPerLev(level)

    ! Allocate the memory for sizeMessage and counter, which will
    ! be used in reduceScatter to determine the amount of messages
    ! and the total size.

    allocate(sizeMessage(nProc), counter(2*nProc), stat=ierr)
    if(ierr /= 0)                             &
         call terminate("indirectHalosPerLevel", &
         "Memory allocation failure for sizeMessage &
         &and counter")

    ! Determine the amount of halo's per processors. Here use is made
    ! of the fact that boundary halo's are flagged with a processor
    ! ID of -1.

    nHaloPerProc = 0
    do i=start,eend
       ii               = indHalo(i)%donorProc+1
       nHaloPerProc(ii) = nHaloPerProc(ii) +1
    enddo

    ! Determine the amount of messages as well as the total number of
    ! halos I have to send and put nHaloPerProc in cumulative
    ! storage format.

    nProcsSend  = 0
    sizeSendBuf = 0

    do i=1,nProc
       if(nHaloPerProc(i) > 0 .and. i /= (myID+1)) then

          ! Something must be send to this proc. Update nProcsSend
          ! and sizeSendBuf and set the two elements of counter
          ! appropriately.

          nProcsSend     = nProcsSend +1
          sizeSendBuf    = sizeSendBuf + nHaloPerProc(i)
          counter(2*i-1) = 1
          counter(2*i)   = nHaloPerProc(i)

       else

          ! Nothing needs to be sent. Set the two elements of
          ! counter to 0.

          counter(2*i-1) = 0
          counter(2*i)   = 0

       endif

       ! Set the size of the message to 2 and put nHaloPerProc in
       ! cumulative storage format.

       sizeMessage(i)  = 2
       nHaloPerProc(i) = nHaloPerProc(i) + nHaloPerProc(i-1)

    enddo

    ! Call reduceScatter to determine the number of processors
    ! from which I receive data and the total amount of data.
    ! tmpBuf is used as a temporary buffer to receive the data.

    call mpi_reduce_scatter(counter, tmpBuf, sizeMessage,           &
         adflow_integer, mpi_sum, ADflow_comm_world, &
         ierr)

    nProcsRecv  = tmpBuf(1)
    sizeRecvBuf = tmpBuf(2)

    ! Allocate the memory for the send buffer.

    allocate(sendBuf(nItemAlloc*sizeSendBuf), stat=ierr)
    if(ierr /= 0)                                           &
         call terminate("indirectHalosPerLevel",               &
         "Memory allocation error for sendBuf.")

    ! Send the data I must send. Use nonblocking sends to
    ! avoid deadlock. nn is used as a counter for the current
    ! message to be sent and ms as starting index for the
    ! active message in sendBuf.

    nn = 0
    ms = 1
    sendproc: do i=1,nProc
       if(counter(2*i) > 0) then

          ! Something must be send to this processor. Update nn.

          nn = nn +1

          ! Fill this part of the send buffer

          proc = i-1
          call fillSendBuf(sendBuf(ms:), proc, entityHalo, &
               transform, level, mm)

          ! Send this buffer. Make sure that integers are used for
          ! count, destination and tag.

          count = nItemSend*mm
          dest  = i-1
          call mpi_isend(sendBuf(ms), count, adflow_integer, dest, dest, &
               ADflow_comm_world, sendRequests(nn), ierr)

          ! Update ms to the starting index in buffer for the next
          ! message to be sent.

          ms = ms + nItemAlloc*mm

       endif
    enddo sendproc

    ! Loop over the boundary halos for this level. The value of start
    ! defined earlier is still okay; only end needs to be redefined.

    eend = nHaloPerLev(level-1) + nHaloPerProc(0)

    bocos: do i=start,eend

       ! Determine the corresponding direct boundary halo.

       ii = indHalo(i)%myDirectHalo

       ! Determine the vector from the direct halo to its donor.
       ! Store this in l1, l2 and l3.

       l1 = entityHalo(ii)%dI - entityHalo(ii)%myI
       l2 = entityHalo(ii)%dJ - entityHalo(ii)%myJ
       l3 = entityHalo(ii)%dK - entityHalo(ii)%myK

       ! Store the info for this boundary halo. Set the boundary
       ! condition, i.e. donorBlock, to the boundary condition of the
       ! closest direct halo. This has been chosen to be the most
       ! important boundary condition.
       ! DonorProc is set to -1 to indicate a boundary halo.

       iihalo = iihalo +1
       entityHalo(iihalo)%myBlock = indHalo(i)%myBlock
       entityHalo(iihalo)%myI     = indHalo(i)%myI
       entityHalo(iihalo)%myJ     = indHalo(i)%myJ
       entityHalo(iihalo)%myK     = indHalo(i)%myK

       entityHalo(iihalo)%donorProc  = -1
       entityHalo(iihalo)%donorBlock = entityHalo(ii)%donorBlock

       ! Determine the vector from my indices to the indices of my
       ! donor. This is the vector from me to the donor indices of
       ! the closest halo, projected on the vector (l1,l2,l3).
       ! Although in general this will not give to indices of my true
       ! donor, this info is stored, because the halo list could be
       ! both a node or a cell halo list. The true donor will be
       ! extracted later.
       ! This implies that my donor could also be a halo, but its
       ! level of indirectness is guaranteed one less than mine.

       ii = l1*(entityHalo(ii)%dI - entityHalo(iihalo)%myI) &
            + l2*(entityHalo(ii)%dJ - entityHalo(iihalo)%myJ) &
            + l3*(entityHalo(ii)%dK - entityHalo(iihalo)%myK)

       entityHalo(iihalo)%dI = ii*l1 + indHalo(i)%myI
       entityHalo(iihalo)%dJ = ii*l2 + indHalo(i)%myJ
       entityHalo(iihalo)%dK = ii*l3 + indHalo(i)%myK

       ! Copy the level of indirectness.

       entityHalo(iihalo)%levOfInd = indHalo(i)%levOfInd

       ! Store the entry of entityHalo in the i,j,k indices
       ! of in entityIndex.

       ii = indHalo(i)%myBlock
       l1 = indHalo(i)%myI
       l2 = indHalo(i)%myJ
       l3 = indHalo(i)%myK

       entityIndex(ii)%entryList(l1,l2,l3) = iihalo

    enddo bocos

    ! Treat the internal block boundary halo's, whose donor is stored
    ! on the same processor. Store this data in localHalos, for
    ! which the memory must be allocated. First determine the number
    ! of local halo's.

    nLocalHalos = nHaloPerProc(myID+1) - nHaloPerProc(myID)
    allocate(localHalos(nItemAlloc*nLocalHalos), stat=ierr)
    if(ierr /= 0)                             &
         call terminate("indirectHalosPerLevel", &
         "Memory allocation error for localHalos")

    proc = myID
    call fillSendBuf(localHalos, proc, entityHalo, transform, &
         level, mm)

    ! Receive the messages in arbitrary sequence, determine the
    ! corresponding halo info and send the message back.

    allocate(recvBuf(nItemAlloc*sizeRecvBuf),  stat=ierr)
    if(ierr /= 0)                             &
         call terminate("indirectHalosPerLevel", &
         "Memory allocation error for recvBuf.")

    ! Initialize ms to 1 and start the loop over the number of
    ! messages to be received.

    ms = 1
    recvproc: do i=1,nProcsRecv

       ! Block until a message arrives.

       call mpi_probe(mpi_any_source, myID, ADflow_comm_world, &
            mpiStatus, ierr)

       ! Find the source and size of the message.

       source = mpiStatus(mpi_source)
       call mpi_get_count(mpiStatus, adflow_integer, sizeRecv, ierr)

       ! Check in debug mode that the incoming message is of
       ! correct size.

       if( debug ) then
          if(sizeRecv == mpi_undefined .or. &
               mod(sizeRecv,nItemSend) /= 0)  &
               call terminate("indirectHalosPerLevel", &
               "Unexpected size of message")
       endif

       ! Receive the message. As it has already arrived a blocking
       ! receive can be used.

       call mpi_recv(recvBuf(ms), sizeRecv, adflow_integer, &
            source, myID, ADflow_comm_world, mpiStatus, ierr)

       ! Determine the number of halo's in the receive buffer.

       mm = sizeRecv/nItemSend

       ! Find the donors for the halo's in the receive buffer as
       ! well as the periodic info.

       call findDonorsRecvBuffer(recvBuf(ms:), mm, entityHalo, &
            entityIndex, level, nItemReturn)

       ! Send the modified receive buffer back to the source processor.

       count = nItemReturn*mm
       call mpi_isend(recvBuf(ms), count, adflow_integer, source, &
            source+1, ADflow_comm_world, recvRequests(i), ierr)

       ! Update the starting index ms for the next message.

       ms = ms + nItemAlloc*mm

    enddo recvproc

    ! Find the donors for the locally stored halo's and store them
    ! in the list.

    call findDonorsRecvBuffer(localHalos, nLocalHalos, entityHalo, &
         entityIndex, level, nItemReturn)

    proc = myID
    call storeHalosInList(localHalos, nLocalHalos, proc, level, &
         nItemReturn, entityHalo, entityIndex, &
         iihalo)

    ! Complete the 1st series of nonblocking sends.

    do i=1,nProcsSend
       call mpi_waitany(nProcsSend, sendRequests, count, mpiStatus, ierr)
    enddo

    ! Loop over the processors to which I sent data to find out the
    ! halo information. Now these messages must be received.

    secondRecv: do i=1,nProcsSend

       ! Block until a message arrives.

       call mpi_probe(mpi_any_source, myID+1, ADflow_comm_world, &
            mpiStatus, ierr)

       ! Find the source and size of the message.

       source = mpiStatus(mpi_source)
       call mpi_get_count(mpiStatus, adflow_integer, sizeRecv, ierr)

       ! Check in debug mode that the incoming message is of
       ! correct size.

       if( debug ) then
          if(sizeRecv == mpi_undefined .or.         &
               mod(sizeRecv,nItemReturn) /= 0)        &
               call terminate("indirectHalosPerLevel", &
               "Unexpected size of message")
       endif

       ! Store the number of halo's in mm.

       mm = sizeRecv/nItemReturn

       ! Receive the message. Use a blocking receive, as the message
       ! has already arrived.

       call mpi_recv(sendBuf, sizeRecv, adflow_integer, source, &
            myID+1, ADflow_comm_world, mpiStatus, ierr)

       ! Store the donors in the list.

       proc = source
       call storeHalosInList(sendBuf, mm, proc, level, nItemReturn, &
            entityHalo, entityIndex, iihalo)
    enddo secondRecv

    ! Complete the second series of nonblocking sends.

    do i=1,nProcsRecv
       call mpi_waitany(nProcsRecv, recvRequests, count, mpiStatus, ierr)
    enddo

    ! Release the memory allocated in this subroutine.

    deallocate(sizeMessage, counter, sendBuf, localHalos, &
         recvBuf, stat=ierr)
    if(ierr /= 0)                             &
         call terminate("indirectHalosPerLevel", &
         "Deallocation error for sizeMessage, etc.")


  end subroutine indirectHalosPerLevel

  !=================================================================

  subroutine fillSendBuf(sendBuf, proc, entityHalo, transform, &
       level, mm)
    !
    !       fillSendBuf fills the buffer, which must be sent to the
    !       given processor.
    !
    use haloList
    use indirectHalo
    use utils, only : delta
    implicit none
    !
    !      Subroutine arguments
    !
    integer(kind=intType), intent(in)  :: proc, level
    integer(kind=intType), intent(out) :: mm

    integer(kind=intType), dimension(:), intent(out) :: sendBuf
    integer(kind=intType), dimension(:,:), intent(in) :: transform

    type(haloListType), dimension(:), intent(in) :: entityHalo
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, ii, jj, m
    integer(kind=intType) :: l1, l2, l3

    integer(kind=intType), dimension(3,3) :: trMat


    ! Initialize m to 0.

    m = 0

    ! Loop over the number of halo's for this processor.

    do i=nHaloPerProc(proc)+1,nHaloPerProc(proc+1)

       ! Abbreviate the current entry in indHalo in ii and the
       ! entry of its corresponding direct halo in jj.

       ii = i + nHaloPerLev(level-1)
       jj = indHalo(ii)%myDirectHalo

       ! Determine the transformation matrix between the block of
       ! the direct halo and its donor block.

       l1 = transform(jj,1)
       l2 = transform(jj,2)
       l3 = transform(jj,3)

       trMat(1,1) = sign(1_intType,l1) * delta(l1,1_intType)
       trMat(2,1) = sign(1_intType,l1) * delta(l1,2_intType)
       trMat(3,1) = sign(1_intType,l1) * delta(l1,3_intType)

       trMat(1,2) = sign(1_intType,l2) * delta(l2,1_intType)
       trMat(2,2) = sign(1_intType,l2) * delta(l2,2_intType)
       trMat(3,2) = sign(1_intType,l2) * delta(l2,3_intType)

       trMat(1,3) = sign(1_intType,l3) * delta(l3,1_intType)
       trMat(2,3) = sign(1_intType,l3) * delta(l3,2_intType)
       trMat(3,3) = sign(1_intType,l3) * delta(l3,3_intType)

       ! Store the direction from the direct to the indirect
       ! halo in l1, l2 and l3

       l1 = indHalo(ii)%myI - entityHalo(jj)%myI
       l2 = indHalo(ii)%myJ - entityHalo(jj)%myJ
       l3 = indHalo(ii)%myK - entityHalo(jj)%myK

       ! Fill the send buffer with block ID and i,j and k
       ! indices of the donor of the direct halo and the
       ! transformed direction to reach the donor of the
       ! indirect halo. For these last three values the
       ! transformation matrix must be applied to l1, l2, l3.

       m = m+1; sendBuf(m) = entityHalo(jj)%donorBlock
       m = m+1; sendBuf(m) = entityHalo(jj)%dI
       m = m+1; sendBuf(m) = entityHalo(jj)%dJ
       m = m+1; sendBuf(m) = entityHalo(jj)%dK

       m = m+1; sendBuf(m) = trMat(1,1)*l1 + trMat(1,2)*l2 + trMat(1,3)*l3
       m = m+1; sendBuf(m) = trMat(2,1)*l1 + trMat(2,2)*l2 + trMat(2,3)*l3
       m = m+1; sendBuf(m) = trMat(3,1)*l1 + trMat(3,2)*l2 + trMat(3,3)*l3

    enddo

    ! Set the return variable mm to the number of halos stored in
    ! the send buffer.

    mm = m/7

  end subroutine fillSendBuf

  !=================================================================

  subroutine findDonorsRecvBuffer(recvBuf, nHalos, entityHalo, &
       entityIndex, level, nItemReturn)
    !
    !       findDonorsRecvBuffer finds the donor cells for the halo
    !       information stored in recvBuf. On return recvBuf contains
    !       for every halo the following information: processor ID,
    !       block ID, the i,j,k indices of the donor cell and periodic
    !       information. The number of periodic subfaces stored is level,
    !       where a 0 indicates that the subface is not periodic.
    !
    use haloList
    use communication
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nHalos, level, nItemReturn

    integer(kind=intType), dimension(:), intent(inout) :: recvBuf

    type(haloListType),  dimension(:), intent(in) :: entityHalo
    type(indexListType), dimension(:), intent(in) :: entityIndex
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j, k, ii, jj, kk, mm, nn
    integer(kind=intType) :: db, l1, L2, l3
    integer(kind=intType) :: nPeriodic

    integer(kind=intType), dimension(:), allocatable :: tmpBuf

    ! Allocate the memory for tmpBuf to the size needed to store the
    ! return information.

    allocate(tmpBuf(nItemReturn*nHalos), stat=ierr)
    if(ierr /= 0)                            &
         call terminate("findDonorsRecvBuffer", &
         "Memory allocation failure for tmpBuf");

    ! Initialize nn and mm to 0. nn is the counter for the incoming
    ! halo information (7 per halo) and mm for the outgoing info
    ! (nItemReturn per halo).

    nn = 0
    mm = 0

    ! Loop over the number of halos in the receive buffer.

    do ii=1,nHalos

       ! Store the incoming information a bit easier.

       db = recvBuf(nn+1)
       i  = recvBuf(nn+2)
       j  = recvBuf(nn+3)
       k  = recvBuf(nn+4)
       l1 = recvBuf(nn+5)
       l2 = recvBuf(nn+6)
       l3 = recvBuf(nn+7)

       ! At the moment i,j,k are the indices of the donor of the direct
       ! halo and l1,l2,l3 the path from the direct halo to the
       ! indirect halo. Add l1, L2, l3 to i,j,k such that they store
       ! the indices of the donor of the indirect halo.

       i = i + l1
       j = j + l2
       k = k + l3

       ! Store the entry (if there is one) in entityHalo in jj.

       jj = entityIndex(db)%entryList(i,j,k)

       ! Now determine the situation we are dealing with.

       if(jj == 0) then

          ! Donor is an owned entity of the block. Store the
          ! appropriate info in tmpBuf. There are no periodic
          ! subfaces for this donor.

          tmpBuf(mm+1) = myID
          tmpBuf(mm+2) = db
          tmpBuf(mm+3) = i
          tmpBuf(mm+4) = j
          tmpBuf(mm+5) = k
          nPeriodic    = 0

       else

          ! Donor is also a halo, but its level of indirectness is at
          ! least one less than the one for which information is to
          ! be found. Still there are two possibilities. Either this
          ! is an internal block boundary halo or a physical boundary
          ! halo. In the former case the corresponding halo is stored,
          ! in the latter case the indices of the boundary halo are
          ! returned.

          if(entityHalo(jj)%donorProc == -1) then

             ! Physical boundary halo. Store the appropriate
             ! info in tmpBuf. No periodic subfaces for this donor.

             tmpBuf(mm+1) = myID
             tmpBuf(mm+2) = db
             tmpBuf(mm+3) = i
             tmpBuf(mm+4) = j
             tmpBuf(mm+5) = k
             nPeriodic    = 0

          else

             ! Internal block boundary halo. Store the appropriate
             ! info in tmpBuf, including the possible periodic
             ! subfaces.

             tmpBuf(mm+1) = entityHalo(jj)%donorProc
             tmpBuf(mm+2) = entityHalo(jj)%donorBlock
             tmpBuf(mm+3) = entityHalo(jj)%dI
             tmpBuf(mm+4) = entityHalo(jj)%dJ
             tmpBuf(mm+5) = entityHalo(jj)%dK

             nPeriodic = entityHalo(jj)%nPeriodicSubfaces
             do kk=1,nPeriodic
                tmpBuf(mm+5+kk) = entityHalo(jj)%periodicSubfaces(kk)
             enddo

          endif

       endif

       ! Fill the remaining part reserved for the periodic subfaces
       ! with 0's. A 0 indicates that no periodic subface is crossed.

       do kk=(nPeriodic+1),level
          tmpBuf(mm+5+kk) = 0
       enddo

       ! Update nn and mm for the next halo.

       nn = nn + 7
       mm = mm + nItemReturn

    enddo

    ! Copy the data from tmpBuf into recvBuf and delete tmpBuf
    ! afterwards.

    nn = nItemReturn*nHalos
    do i=1,nn
       recvBuf(i) = tmpBuf(i)
    enddo

    deallocate(tmpBuf, stat=ierr)
    if(ierr /= 0)                            &
         call terminate("findDonorsRecvBuffer", &
         "Deallocation failure for tmpBuf")

  end subroutine findDonorsRecvBuffer

  !=================================================================

  subroutine storeHalosInList(buffer, bufSize, proc, level, &
       nItemReturn, entityHalo,      &
       entityIndex, iihalo)
    !
    !       storeHalosInList stores the halo info present in buf, which
    !       has been retreived from the given processor, in the correct
    !       place in entityHalo and entityIndex.
    !
    use haloList
    use indirectHalo
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in)    :: bufSize, proc
    integer(kind=intType), intent(in)    :: level, nItemReturn
    integer(kind=intType), intent(inout) :: iihalo

    integer(kind=intType), dimension(:), intent(in) :: buffer

    type(haloListType),  dimension(:), intent(inout) :: entityHalo
    type(indexListType), dimension(:), intent(inout) :: entityIndex
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j, k
    integer(kind=intType) :: ii, nn, blockID, iii
    integer(kind=intType) :: nPeriodic

    ! Store the start index (-1) for this processor in indHalo in nn.

    nn = nHaloPerLev(level-1) + nHaloPerProc(proc)

    ! Loop over the number of halo's stored in the buffer.

    do ii=1,bufSize

       iii = (ii-1)*nItemReturn

       ! Update the counters iihalo and nn.

       iihalo = iihalo +1
       nn     = nn +1

       ! Store the i,j,k indices and the block ID of the current
       ! halo a bit easier.

       blockID = indHalo(nn)%myBlock
       i       = indHalo(nn)%myI
       j       = indHalo(nn)%myJ
       k       = indHalo(nn)%myK

       ! Store the entry of entityHalo in the i,j,k indices
       ! of in entityIndex.

       entityIndex(blockID)%entryList(i,j,k) = iihalo

       ! Store the info of the current halo in entityHalo.

       entityHalo(iihalo)%myBlock = blockID

       entityHalo(iihalo)%myI = i
       entityHalo(iihalo)%myJ = j
       entityHalo(iihalo)%myK = k

       entityHalo(iihalo)%donorProc  = buffer(iii+1)
       entityHalo(iihalo)%donorBlock = buffer(iii+2)

       entityHalo(iihalo)%dI = buffer(iii+3)
       entityHalo(iihalo)%dJ = buffer(iii+4)
       entityHalo(iihalo)%dK = buffer(iii+5)

       ! Determine the number of periodic subfaces in the buffer.

       nPeriodic = 0
       do i=6,nItemReturn
          if(buffer(iii+i) > 0) nPeriodic = nPeriodic + 1
       enddo

       ! Check if the corresponding direct halo borders a periodic
       ! subface and update nPeriodic accordingly.

       j = indHalo(nn)%myDirectHalo
       nPeriodic = nPeriodic + entityHalo(j)%nPeriodicSubfaces

       ! If periodic subfaces are present for this halo, allocate
       ! the memory for periodicSubfaces and copy the data from
       ! both the buffer and the direct halo.

       if(nPeriodic > 0) then
          entityHalo(iihalo)%nPeriodicSubfaces = nPeriodic
          allocate(entityHalo(iihalo)%periodicSubfaces(nPeriodic), &
               stat=ierr)
          if(ierr /= 0)                        &
               call terminate("storeHalosInList", &
               "Memory allocation failure for &
               &periodicSubfaces")
          nPeriodic = 0
          do i=6,nItemReturn
             if(buffer(iii+i) > 0) then
                nPeriodic = nPeriodic + 1
                entityHalo(iihalo)%periodicSubfaces(nPeriodic) = &
                     buffer(iii+i)
             endif
          enddo

          do i=1,entityHalo(j)%nPeriodicSubfaces
             nPeriodic = nPeriodic + 1
             entityHalo(iihalo)%periodicSubfaces(nPeriodic) = &
                  entityHalo(j)%periodicSubfaces(i)
          enddo
       endif

    enddo

    ! Check in debug mode if the buffer size was correct.

    if( debug ) then
       if(nn /= nHaloPerLev(level-1) + nHaloPerProc(proc+1)) &
            call terminate("storeHalosInList",                  &
            "Something wrong with buffer size")
    endif

  end subroutine storeHalosInList

  subroutine qsortPeriodicSubfacesHaloType(arr, nn)
    !
    !       qsortPeriodicSubfacesHaloType sorts the given number of halo's
    !       with periodic subfaces in increasing order based on the
    !       <= operator for this derived data type.
    !
    use periodicInfo
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn

    type(periodicSubfacesHaloType), dimension(*), intent(inout) :: arr
    !
    !      Local variables.
    !
    integer(kind=intType), parameter :: m = 7

    integer(kind=intType) :: nStack
    integer(kind=intType) :: i, j, k, r, l, jStack, ii

    integer :: ierr

    type(periodicSubfacesHaloType) :: a, tmp

    integer(kind=intType), allocatable, dimension(:) :: stack
    integer(kind=intType), allocatable, dimension(:) :: tmpStack

    ! Allocate the memory for stack.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                                     &
         call terminate("qsortPeriodicSubfacesHaloType", &
         "Memory allocation failure for stack")

    ! Initialize the variables that control the sorting.

    jStack = 0
    l      = 1
    r      = nn

    ! Start of the algorithm

    do

       ! Check for the size of the subarray.

       if((r-l) < m) then

          ! Perform insertion sort

          do j=l+1,r
             a = arr(j)
             do i=(j-1),l,-1
                if(arr(i) <= a) exit
                arr(i+1) = arr(i)
             enddo
             arr(i+1) = a
          enddo

          ! In case there are no more elements on the stack, exit from
          ! the outermost do-loop. Algorithm has finished.

          if(jStack == 0) exit

          ! Pop stack and begin a new round of partitioning.

          r = stack(jStack)
          l = stack(jStack-1)
          jStack = jStack - 2

       else

          ! Subarray is larger than the threshold for a linear sort.
          ! Choose median of left, center and right elements as
          ! partitioning element a.
          ! Also rearrange so that (l) <= (l+1) <= (r).

          k = (l+r)/2
          tmp      = arr(k)             ! Swap the elements
          arr(k)   = arr(l+1)           ! k and l+1.
          arr(l+1) = tmp

          if(arr(r) < arr(l)) then
             tmp    = arr(l)             ! Swap the elements
             arr(l) = arr(r)             ! r and l.
             arr(r) = tmp
          endif

          if(arr(r) < arr(l+1)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(r)           ! r and l+1.
             arr(r)   = tmp
          endif

          if(arr(l+1) < arr(l)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(l)           ! l and l+1.
             arr(l)   = tmp
          endif

          ! Initialize the pointers for partitioning.

          i = l+1
          j = r
          a = arr(l+1)

          ! The innermost loop

          do

             ! Scan up to find element >= a.
             do
                i = i+1
                if(a <= arr(i)) exit
             enddo

             ! Scan down to find element <= a.
             do
                j = j-1
                if(arr(j) <= a) exit
             enddo

             ! Exit the loop in case the pointers i and j crossed.

             if(j < i) exit

             ! Swap the element i and j.

             tmp    = arr(i)
             arr(i) = arr(j)
             arr(j) = tmp
          enddo

          ! Swap the entries j and l+1. Remember that a equals
          ! arr(l+1).

          arr(l+1) = arr(j)
          arr(j)   = a

          ! Push pointers to larger subarray on stack,
          ! process smaller subarray immediately.

          jStack = jStack + 2
          if(jStack > nStack) then

             ! Storage of the stack is too small. Reallocate.

             allocate(tmpStack(nStack), stat=ierr)
             if(ierr /= 0)                                     &
                  call terminate("qsortPeriodicSubfacesHaloType", &
                  "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                                     &
                  call terminate("qsortPeriodicSubfacesHaloType", &
                  "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                                     &
                  call terminate("qsortPeriodicSubfacesHaloType", &
                  "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                                     &
                  call terminate("qsortPeriodicSubfacesHaloType", &
                  "Deallocation error for tmpstack")
          endif

          if((r-i+1) >= (j-l)) then
             stack(jStack)   = r
             r               = j-1
             stack(jStack-1) = j
          else
             stack(jStack)   = j-1
             stack(jStack-1) = l
             l               = j
          endif

       endif
    enddo

    ! Release the memory of stack.

    deallocate(stack, stat=ierr)
    if(ierr /= 0)                                     &
         call terminate("qsortPeriodicSubfacesHaloType", &
         "Deallocation error for stack")

    ! Check in debug mode whether the array is really sorted.

    if( debug ) then
       do i=1,(nn-1)
          if(arr(i+1) < arr(i))                             &
               call terminate("qsortPeriodicSubfacesHaloType", &
               "Array is not sorted correctly")
       enddo
    endif

  end subroutine qsortPeriodicSubfacesHaloType

  subroutine determinePeriodicFaces
    !
    !       determinePeriodicFaces determines and stores the number of
    !       periodic faces present in the complete mesh. The sequence of
    !       storing the data is such that the array periodicGlobal is
    !       sorted with the definition of the < operator for this datatype.
    !
    use cgnsGrid
    use periodicInfo
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr
    integer(kind=intType) :: nn, ii, i

    ! Determine the number of periodic faces present in the cgns grid.

    nPeriodicGlobal = 0
    do nn=1,cgnsNDom
       do i=1,cgnsDoms(nn)%n1to1
          if( cgnsDoms(nn)%conn1to1(i)%periodic ) &
               nPeriodicGlobal = nPeriodicGlobal + 1
       enddo
    enddo

    ! Allocate the memory for periodicGlobal.

    allocate(periodicGlobal(nPeriodicGlobal), stat=ierr)
    if(ierr /= 0)                              &
         call terminate("determinePeriodicFaces", &
         "Memory allocation failure for periodicGlobal")

    ! Repeat the loop over the faces of the cgns grid and store the
    ! periodic faces.

    ii = 0
    do nn=1,cgnsNDom
       do i=1,cgnsDoms(nn)%n1to1
          if( cgnsDoms(nn)%conn1to1(i)%periodic ) then

             ii = ii + 1
             periodicGlobal(ii)%cgnsBlock   = nn
             periodicGlobal(ii)%cgnsSubface = i

          endif
       enddo
    enddo

  end subroutine determinePeriodicFaces

  function bsearchCGNSPeriodicType(key, base)
    !
    !       bsearchCGNSPeriodicType returns the index in base where key
    !       is stored. A binary search algorithm is used here, so it is
    !       assumed that base is sorted in increasing order. In case key
    !       appears more than once in base, the result is arbitrary.
    !       If key is not found, a zero is returned.
    !
    use periodicInfo
    implicit none
    !
    !      Function type
    !
    integer(kind=intType) :: bsearchCGNSPeriodicType
    !
    !      Function arguments.
    !
    type(cgnsPeriodicType), intent(in)               :: key
    type(cgnsPeriodicType), dimension(:), intent(in) :: base
    integer(kind=intType)                            :: nn
    !
    !      Local variables.
    !
    integer(kind=intType) :: ii, pos, start
    logical               :: entryFound

    ! Initialize some values.

    start       = 1
    ii          = size(base)
    entryFound = .false.

    ! Binary search to find key.

    do
       ! Condition for breaking the loop

       if(ii == 0) exit

       ! Determine the position in the array to compare.

       pos = start + ii/2

       ! In case this is the entry, break the search loop.

       if(base(pos) == key) then
          entryFound = .true.
          exit
       endif

       ! In case the search key is larger than the current position,
       ! only parts to the right must be searched. Remember that base
       ! is sorted in increasing order. Nothing needs to be done if the
       ! key is smaller than the current element.

       if(base(pos) < key) then
          start = pos +1
          ii    = ii -1
       endif

       ! Modify ii for the next branch to search.

       ii = ii/2
    enddo

    ! Set bsearchCGNSPeriodicType.
    ! This depends whether the key was found.

    if( entryFound ) then
       bsearchCGNSPeriodicType = pos
    else
       bsearchCGNSPeriodicType = 0
    endif

  end function bsearchCGNSPeriodicType
  subroutine qsortIndHaloType(arr, nn)
    !
    !       qsortIndHaloType sorts the given number of indirect halo's
    !       in increasing order based on the <= operator for this derived
    !       data type.
    !
    use indirectHalo
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn

    type(indirectHaloType), dimension(*), intent(inout) :: arr
    !
    !      Local variables.
    !
    integer(kind=intType), parameter :: m = 7

    integer(kind=intType) :: nStack
    integer(kind=intType) :: i, j, k, r, l, jStack, ii

    integer :: ierr

    type(indirectHaloType) :: a, tmp

    integer(kind=intType), allocatable, dimension(:) :: stack
    integer(kind=intType), allocatable, dimension(:) :: tmpStack

    ! Allocate the memory for stack.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                        &
         call terminate("qsortIndHaloType", &
         "Memory allocation failure for stack")

    ! Initialize the variables that control the sorting.

    jStack = 0
    l      = 1
    r      = nn

    ! Start of the algorithm

    do

       ! Check for the size of the subarray.

       if((r-l) < m) then

          ! Perform insertion sort

          do j=l+1,r
             a = arr(j)
             do i=(j-1),l,-1
                if(arr(i) <= a) exit
                arr(i+1) = arr(i)
             enddo
             arr(i+1) = a
          enddo

          ! In case there are no more elements on the stack, exit from
          ! the outermost do-loop. Algorithm has finished.

          if(jStack == 0) exit

          ! Pop stack and begin a new round of partitioning.

          r = stack(jStack)
          l = stack(jStack-1)
          jStack = jStack - 2

       else

          ! Subarray is larger than the threshold for a linear sort.
          ! Choose median of left, center and right elements as
          ! partitioning element a.
          ! Also rearrange so that (l) <= (l+1) <= (r).

          k = (l+r)/2
          tmp      = arr(k)             ! Swap the elements
          arr(k)   = arr(l+1)           ! k and l+1.
          arr(l+1) = tmp

          if(arr(r) < arr(l)) then
             tmp    = arr(l)             ! Swap the elements
             arr(l) = arr(r)             ! r and l.
             arr(r) = tmp
          endif

          if(arr(r) < arr(l+1)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(r)           ! r and l+1.
             arr(r)   = tmp
          endif

          if(arr(l+1) < arr(l)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(l)           ! l and l+1.
             arr(l)   = tmp
          endif

          ! Initialize the pointers for partitioning.

          i = l+1
          j = r
          a = arr(l+1)

          ! The innermost loop

          do

             ! Scan up to find element >= a.
             do
                i = i+1
                if(a <= arr(i)) exit
             enddo

             ! Scan down to find element <= a.
             do
                j = j-1
                if(arr(j) <= a) exit
             enddo

             ! Exit the loop in case the pointers i and j crossed.

             if(j < i) exit

             ! Swap the element i and j.

             tmp    = arr(i)
             arr(i) = arr(j)
             arr(j) = tmp
          enddo

          ! Swap the entries j and l+1. Remember that a equals
          ! arr(l+1).

          arr(l+1) = arr(j)
          arr(j)   = a

          ! Push pointers to larger subarray on stack,
          ! process smaller subarray immediately.

          jStack = jStack + 2
          if(jStack > nStack) then

             ! Storage of the stack is too small. Reallocate.

             allocate(tmpStack(nStack), stat=ierr)
             if(ierr /= 0)                        &
                  call terminate("qsortIndHaloType", &
                  "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                        &
                  call terminate("qsortIndHaloType", &
                  "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                        &
                  call terminate("qsortIndHaloType", &
                  "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                        &
                  call terminate("qsortIndHaloType", &
                  "Deallocation error for tmpstack")
          endif

          if((r-i+1) >= (j-l)) then
             stack(jStack)   = r
             r               = j-1
             stack(jStack-1) = j
          else
             stack(jStack)   = j-1
             stack(jStack-1) = l
             l               = j
          endif

       endif
    enddo

    ! Release the memory of stack.

    deallocate(stack, stat=ierr)
    if(ierr /= 0)                        &
         call terminate("qsortIndHaloType", &
         "Deallocation error for stack")

    ! Check in debug mode whether the array is really sorted.

    if( debug ) then
       do i=1,(nn-1)
          if(arr(i+1) < arr(i))                &
               call terminate("qsortIndHaloType", &
               "Array is not sorted correctly")
       enddo
    endif

  end subroutine qsortIndHaloType

end module pointMatchedCommPattern
