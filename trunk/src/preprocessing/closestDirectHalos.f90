!
!      ******************************************************************
!      *                                                                *
!      * File:          closestDirectHalos.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-31-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine closestDirectHalos(entityHalo, entityIndex, &
                                     start, nLevel, offset, gridLevel)
!
!      ******************************************************************
!      *                                                                *
!      * closestDirectHalos determines the number of indirect halo's    *
!      * to be treated and its corresponding direct halo.               *
!      *                                                                *
!      ******************************************************************
!
       use block
       use bcHalo
       use haloList
       use indirectHalo
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
!      Function definitions.
!
       integer(kind=intType) :: getNumberIndirectHalos
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
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

       !=================================================================

       function getNumberIndirectHalos(start, nLevel, offset, gridLevel)
!
!      ******************************************************************
!      *                                                                *
!      * getNumberIndirectHalos determines the number of indirect       *
!      * halo's for which the donor must be determined.                 *
!      *                                                                *
!      ******************************************************************
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
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
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
