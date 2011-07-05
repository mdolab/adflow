!
!      ******************************************************************
!      *                                                                *
!      * File:          createCoarseBoundary.f90                        *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 02-06-2005                                      *
!      * Last modified: 10-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine createCoarseBoundary(level, sps, nHaloBndryAdd)
!
!      ******************************************************************
!      *                                                                *
!      * CreateCoarseBoundary creates the ibndry overset list for the   *
!      * given level and spectral solution from the next finer level.   *
!      * At the end of this routine the blocks will have a preliMinary  *
!      * boundary list and iblank array. It is preliMinary because      *
!      * adjacent blocks could have decided that some of their halos    *
!      * need to be part of the boundary, in which case they need to    *
!      * tell the connecting block to add them. This is not dealt with  *
!      * in this routine. Rather, the iblank on these halos are set to  *
!      * a particular value as a signal.                                *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use boundaryList
       use communication
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in)  :: level, sps
       integer(kind=intType), intent(out) :: nHaloBndryAdd
!
!      Local variables.
!
       integer               :: ierr
       integer(kind=intType) :: i, j, k, l, m, n, fineLevel
       integer(kind=intType) :: i1, i2, j1, j2, k1, k2, fact
       integer(kind=intType) :: del(3), cp(3), cm(3)
 
       integer(kind=intType), pointer :: iblankFine(:,:,:)

       logical :: doChecks
!
!      Function definitions.
!
       logical :: is1to1Halo
!
!      Interfaces.
!
       interface
         subroutine reallocateInteger(intArray, newSize, oldSize, &
                                       alwaysFreeMem)
           use precision
           implicit none

           integer(kind=intType), dimension(:), pointer :: intArray
           integer(kind=intType), intent(in) :: newSize, oldSize
           logical, intent(in) :: alwaysFreeMem
         end subroutine reallocateInteger

         !===============================================================

         subroutine reallocateInteger2(intArray,            &
                                        newSize1, newSize2, &
                                        oldSize1, oldSize2, &
                                        alwaysFreeMem)
           use precision
           implicit none

           integer(kind=intType), dimension(:,:), pointer :: intArray
           integer(kind=intType), intent(in) :: newSize1, newSize2, &
                                                   oldSize1, oldSize2
           logical, intent(in) :: alwaysFreeMem
         end subroutine reallocateInteger2
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Next finer level.
 
       fineLevel = level - 1
 
       ! Loop over the domains and use the finer level to initialize.

       domains1: do n = 1,nDom

         ! Allocate memory for the iblank array of this coarse block
         ! only if needed since this is only done once in preprocessing.

         if (.not. associated(flowDoms(n,level,sps)%iblank)) then
           i = flowDoms(n,level,sps)%ib
           j = flowDoms(n,level,sps)%jb
           k = flowDoms(n,level,sps)%kb

           allocate(flowDoms(n,level,sps)%iblank(0:i,0:j,0:k), &
                    stat=ierr)
           if (ierr /= 0)                           &
             call terminate("createCoarseBoundary", &
                            "Memory allocation failure for iblank")
         end if

         ! Allocate some memory for the coarse block's boundary indices.
         ! We assume the actual # of cells will not be greater than the
         ! fine #. Extra memory is released later in this routine.

         i = flowDoms(n,fineLevel,sps)%nCellsOverset &
           + flowDoms(n,fineLevel,sps)%nOrphans

         allocate(flowDoms(n,level,sps)%ibndry(3,i), stat=ierr)
         if (ierr /= 0)                           &
           call terminate("createCoarseBoundary", &
                          "Memory allocation failure for ibndry")

         ! Set pointers to make things readable.
 
         call setPointers(n, level, sps)

         ! Initialize iblank on the cell halos. See the description of
         ! the called routine. OversetOuterBoundaries are set to -1
         ! which serve as a check in the next procedure.

         call changeIblanks(.true., -1_intType)

         ! Set the pointer to the iblanks of the finer level.
 
         iblankFine => flowDoms(n,fineLevel,sps)%iblank
 
         ! Initialize the counter for overset cells and then loop over
         ! the interior cells and set the iblank.
 
         m = 0
         do k = 2,kl
           k1 = mgKFine(k,1)
           k2 = mgKFine(k,2)
           do j = 2,jl
             j1 = mgJFine(j,1)
             j2 = mgJFine(j,2)
             do i = 2,il
               i1 = mgIFine(i,1)
               i2 = mgIFine(i,2)
 
               ! Initialize the iblank based on a few simple rules
               ! designed to ensure that the holes contract rather than
               ! expand into the originial flow field of the fine level.
               ! (1) if the coarse cell contains any of the finer field
               !     then consider it also as field.
               ! (2) if the coarse cell contains any of the fine boundary
               !     then add it to the boundary only if it contains no
               !     no field cells (contracting hole).
               ! (3) otherwise it is a hole.
 
               if (any(iblankFine(i1:i2,j1:j2,k1:k2) == 1)) then
                 iblank(i,j,k) = 1
               else if (any(iblankFine(i1:i2,j1:j2,k1:k2) >= 9)) then
                 m = m + 1
                 iblank(i,j,k) = 9
                 ibndry(:,m) = (/ i, j, k /)
               else
                 iblank(i,j,k) = 0
               end if
 
             end do
           end do
         end do

       ! Loop over the boundary conditions which are specified as
       ! overset outer boundaries and perform the same type process as
       ! above. The reason is to allow the boundary to extend into all
       ! the halos so as to prevent an overlap problems on coarse grids.

         bocos: do l = 1,nBocos
           if (BCType(l) == OversetOuterBound) then

             do k = kcBeg(l),kcEnd(l)
               k1 = mgKFine(k,1)
               k2 = mgKFine(k,2)
               do j = jcBeg(l),jcEnd(l)
                 j1 = mgJFine(j,1)
                 j2 = mgJFine(j,2)
                 do i = icBeg(l),icEnd(l)
                   i1 = mgIFine(i,1)
                   i2 = mgIFine(i,2)

                   ! Skip this cell unless it is actually part of the
                   ! outer boundary. This prevents 1-to-1 halos from
                   ! being added to the fringe here.

                   if (iblank(i,j,k) /= -1) cycle

                   ! Check if any of the fine cells are part of the
                   ! boundary and if so, make the coarse cell part as
                   ! well. Note there is not 'else' part because the
                   ! iblanks here have already been initialized.

                   if (any(iblankFine(i1:i2,j1:j2,k1:k2) >= 9)) then
                     m = m + 1
                     iblank(i,j,k) = 9
                     ibndry(:,m) = (/ i, j, k /)
                   end if

                 end do
               end do
             end do

           end if
         end do bocos
 
         ! Copy the current # of boundary cells to the flowDoms.
 
         flowDoms(n,level,sps)%nCellsOverset = m
 
       end do domains1
 
       ! At this point we need to make an iblank exchange among the
       ! blocks. The reason is because the next few steps will add cells
       ! to the boundary where needed, and all additions are based on a
       ! current boundary being adjacent to a field cell. Since no
       ! field cells will change at this point, we make the exchange so
       ! the rest of the process goes smoothly.
 
       call exchangeIblanks(level, sps, commPatternCell_1st, &
                                        internalCell_1st)
 
       ! We need to check the fringe and make sure that a double fringe
       ! layer exists to maintain accuracy when a higher-order scheme
       ! is used. Note that any cells added here are entirely inside
       ! the hole of the finer block, thus we will also store the
       ! nearest neighbor from the current boundary to cells we add.
       ! This will be needed to find their estimated donor locations.
 
       ! Initialize the number of cells that this processor needs to
       ! tell other blocks to add to 0.
 
       nHaloBndryAdd = 0
 
       ! Loop over the local domains again.
 
       domains2: do n = 1,nDom

         ! Set pointers to make things readable.

         call setPointers(n, level, sps)

         ! Reset the counter for overset cells, and set the beginning
         ! index for boundary cells entirely over holes.
 
         m = nCellsOverset
         blockBndry(n)%istartOverHole = m + 1
 
         ! Allocate memory to store the nearest neighbors.
 
         i = blockBndry(n)%istartOverHole
         j = ubound(ibndry,2)
 
         allocate(blockBndry(n)%nearestBndry(i:j), stat=ierr)
         if (ierr /= 0)                           &
           call terminate("createCoarseBoundary", &
                          "Memory allocation failure for nearest")
 
         ! Loop over the boundary already created and do the checks for
         ! each coordinate direction k.
 
         do i = 1,nCellsOverset
           do k = 1,3
 
             ! Find the indices for the adjacent cells in this direction.
 
             del    = 0
             del(k) = 1
 
             cp = ibndry(:,i) + del
             cm = ibndry(:,i) - del
 
             ! If either of the adjacent cells is in the field then we
             ! will perform 2 checks. In both checks. We only add owned
             ! cells or direct halo cells that have their boco type
             ! specified as "oversetOuterBound" (i.e. only cells with
             ! iblanks currently < 0 will be added.
             !
             ! If a 1-to-1 halo cell is deemed to be added, we just set
             ! the iblank to 10 + the nearest neighbor index. This will
             ! allow us to later have the neighboring block make sure it
             ! is added to its boundary and also send it the initial
             ! donor estimate.
             !
             ! Note that cells with actual bcs are never added because
             ! their iblanks are > 0 here.

             if (iblank(cp(1),cp(2),cp(3)) == 1) then
               doChecks = .true.
               fact = -1
             else if (iblank(cm(1),cm(2),cm(3)) == 1) then
               doChecks = .true.
               fact = 1 
             else
               doChecks = .false.
             end if

             if (doChecks) then
 
               ! Search the 3 x 3 face normal to the current direction
               ! for holes. These cells are added to the boundary since
               ! they will be needed for computing the viscous residual
               ! of the adjacent field cell.
 
               do k1 = ibndry(3,i)-1+del(3), ibndry(3,i)+1-del(3)
                 do j1 = ibndry(2,i)-1+del(2), ibndry(2,i)+1-del(2)
                   do i1 = ibndry(1,i)-1+del(1), ibndry(1,i)+1-del(1)
 
                     if (iblank(i1,j1,k1) <= 0) then
 
                       if (is1to1Halo(i1,j1,k1)) then
                         iblank(i1,j1,k1) = 10 + i
                         nHaloBndryAdd  = nHaloBndryAdd + 1
                       else
                         m = m + 1
                         iblank(i1,j1,k1) = 9
                         ibndry(:,m) = (/ i1, j1, k1 /)
                         blockBndry(n)%nearestBndry(m) = i
                       end if
 
                     end if
 
                   end do
                 end do
               end do
 
               ! Loop over the size of the fringe in the current
               ! direction, away from the adjacent field cell. If one
               ! of the those cells is a hole then add it to the
               ! boundary since it will be needed for computing the
               ! inviscid residual for a higher-order scheme.

               do j = 1,fringeSize-1

                 del(k) = j*fact
                 cp = ibndry(:,i) + del

                 if (iblank(cp(1),cp(2),cp(3)) <= 0) then
 
                   if (is1to1Halo(cp(1),cp(2),cp(3))) then
                     iblank(cp(1),cp(2),cp(3)) = 10 + i
                     nHaloBndryAdd = nHaloBndryAdd + 1
                   else
                     m = m + 1
                     iblank(cp(1),cp(2),cp(3)) = 9
                     ibndry(:,m) = cp
                     blockBndry(n)%nearestBndry(m) = i
                   end if
 
                 end if
               end do
             end if
           end do
         end do
 
         ! Copy the current # of boundary cells to the flowDoms.
 
         flowDoms(n,level,sps)%nCellsOverset = m
 
       end do domains2
 
       ! Release the extra memory allocated for the boundary indices
       ! and the nearest neighbors.
 
       domains3: do n = 1,nDom
         i = flowDoms(n,     level,sps)%nCellsOverset
         j = ubound(flowDoms(n,level,sps)%ibndry,2)
         if (j == 0) cycle
         call reallocateInteger2(flowDoms(n,level,sps)%ibndry, &
                                 3_intType, i, 3_intType, j, .true.)
         i = i - blockBndry(n)%istartOverHole + 1
         j = size(blockBndry(n)%nearestBndry)
         call reallocateInteger(blockBndry(n)%nearestBndry, &
                                i, j, .true.)
       end do domains3
 
       end subroutine createCoarseBoundary
 
       !=================================================================

       logical function is1to1Halo(i,j,k)
!
!      ******************************************************************
!      *                                                                *
!      * Is1to1Halo determines whether a cell is a 1to1 connectivity    *
!      * halo, which is distinguished from a bc halo by the cell's      *
!      * iblank value.                                                  *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       implicit none
!
!      Function arguments
!
       integer(kind=intType), intent(in) :: i, j, k
!
!      Local variables.
!
       integer(kind=intType), dimension(3) :: ii
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize the result and if the iblank is < 0 immediately
       ! return - this is a bc halo specified as oversetOuterBound.
 
       is1to1Halo = .false.
       if (iblank(i,j,k) < 0) return
 
       ii(1) = max(2-i, i-il, 0_intType)
       ii(2) = max(2-j, j-jl, 0_intType)
       ii(3) = max(2-k, k-kl, 0_intType)
 
       if (sum(ii) > 0) is1to1Halo = .true.
 
       end function is1to1Halo
