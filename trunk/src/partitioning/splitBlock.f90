!
!      ******************************************************************
!      *                                                                *
!      * File:          splitBlock.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-23-2004                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine splitBlock(compBlock, nSub, nCells, ranges)
!
!      ******************************************************************
!      *                                                                *
!      * splitBlock tries to split the given computational block into   *
!      * the desired number of subblocks nSub. However it can happen    *
!      * that nSub is a strange number and a different splitting is     *
!      * performed. On return, nSub contains the actual number into the *
!      * block is split. This number is smaller or equal to nSub on     *
!      * entry. As it is possible that the computational block itself   *
!      * is a subblock of an original cgns block, on return ranges will *
!      * contain the nodal ranges of the subblocks in the original cgns *
!      * block. The splitting attempts to keep the needed multigrid     *
!      * capabilities as much as possible.                              *
!      * A recursive bisection algorithm is used.                       *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use constants
       use inputIteration
       use partitionMod
       implicit none
!
!      Subroutine arguments.
!
       type(distributionBlockType), intent(in)    :: compBlock
       integer(kind=intType),       intent(in)    :: nCells
       integer(kind=intType),       intent(inout) :: nSub

       integer(kind=intType), dimension(nSub,3,2), intent(out) :: ranges
!
!      Local variables.
!
       integer(kind=intType) :: nLevels, level, nn, mm, nTarget
       integer(kind=intType) :: nSplit, nSplitNew

       integer(kind=intType), dimension(nSub) :: nSubblocks

       logical, dimension(3) :: viscousDir
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the viscous directions of the block.

       viscousDir = .false.
       do nn=1,compBlock%nBocos

         ! Check for viscous boundary conditions and if found set the
         ! corresponding direction to viscous.

         if(compBlock%BCType(nn) == NSWallAdiabatic .or. &
            compBlock%BCType(nn) == NSWallIsothermal) then

           select case (compBlock%BCFaceID(nn))
             case (iMin,iMax)
               viscousDir(1) = .true.

             case (jMin,jMax)
               viscousDir(2) = .true.

             case (kMin,kMax)
               viscousDir(3) = .true.
           end select
         endif
       enddo

       ! Set viscousDir to .false. if an explicit smoother is used.

       if(smoother == RungeKutta) viscousDir = .false.

       ! Determine the number of levels in the bisection.

       nLevels = log(real(nSub,realType))/log(two)
       if(2**nLevels < nSub) nLevels = nLevels + 1

       ! Initialize the range of the first splittable block to the
       ! entire block.

       ranges(1,1,1) = 1; ranges(1,1,2) = compBlock%il
       ranges(1,2,1) = 1; ranges(1,2,2) = compBlock%jl
       ranges(1,3,1) = 1; ranges(1,3,2) = compBlock%kl

       ! Initialize the number of blocks to be split to 1 and the
       ! number of blocks it should be split into to nSub.

       nSplit = 1
       nSubblocks(1) = nSub

       ! Loop over the number of levels

       levelLoop: do level=1,nLevels

         ! Initialize nSplitNew to nSplit; it will be used to store the
         ! position of the new subblocks. At the end of this loop this
         ! value will be set to nSplit for the new level.

         nSplitNew = nSplit

         ! Loop over the number of blocks to be split.

         blockSplitLoop: do nn=1,nSplit

           ! Determine the situation we are having it.

           select case (nSubblocks(nn))

             case (2_intType,3_intType)

               ! Subblock must be split in either two or three. Store
               ! the number into the new subblocks must be in nSubblocks.
               ! As the routine split2block puts the subblock with the
               ! target number of cells in position nSplitNew, this
               ! block should not be split any further and thus gets a
               ! value of 1. Position nn gets the old value -1, which
               ! will be either 1 or 2

               nSplitNew             = nSplitNew + 1
               nSubblocks(nSplitNew) = 1
               nSubblocks(nn)        = nSubblocks(nn) - 1

               ! Split the block into two, where a block with nCells
               ! should be split off.

               call split2block(nSub, nn, nSplitNew, nCells, ranges, &
                                viscousDir)

             !===========================================================

             case (4_intType:)

               ! Subblock must be split in 4 or more. First determine the
               ! number of subblocks into the two new subblocks must be
               ! split in the next round. Save the old value of
               ! nSubblocks(nn) in mm.

               mm = nSubblocks(nn)

               nSplitNew             = nSplitNew + 1
               nSubblocks(nSplitNew) = nSubblocks(nn)/2
               nSubblocks(nn)        = nSubblocks(nn) &
                                     - nSubblocks(nSplitNew)

               ! Determine the number of cells for the new subblock
               ! nSplitNew.

               nTarget = (ranges(nn,1,2) - ranges(nn,1,1)) &
                       * (ranges(nn,2,2) - ranges(nn,2,1)) &
                       * (ranges(nn,3,2) - ranges(nn,3,1))
               nTarget = nTarget*(real(nSubblocks(nSplitNew),realType) &
                       /          real(mm,realType))

               ! Split the block into two, where one block should get
               ! approximately nTarget cells.

               call split2block(nSub, nn, nSplitNew, nTarget, ranges, &
                                viscousDir)

           end select

         enddo blockSplitLoop

         ! Set nSplit to nSplitNew for the next level of bisection.

         nSplit = nSplitNew

       enddo levelLoop

       ! Get rid of the possible zero cell ranges. In that case
       ! nSub will change as well. Add the offset of the computational
       ! block, as this might be a subblock itself.

       mm   = nSub
       nSub = 0
       do nn=1,mm

         ! Test whether this is a non-zero cell subblock.

         if(ranges(nn,1,2) > ranges(nn,1,1) .and. &
            ranges(nn,2,2) > ranges(nn,2,1) .and. &
            ranges(nn,3,2) > ranges(nn,3,1)) then

           ! This is a valid subblock. Update nSub and copy the range
           ! including the offset of the computational block.

           nSub = nSub + 1

           ranges(nSub,1,1) = ranges(nn,1,1) + compBlock%iBegor - 1
           ranges(nSub,2,1) = ranges(nn,2,1) + compBlock%jBegor - 1
           ranges(nSub,3,1) = ranges(nn,3,1) + compBlock%kBegor - 1

           ranges(nSub,1,2) = ranges(nn,1,2) + compBlock%iBegor - 1
           ranges(nSub,2,2) = ranges(nn,2,2) + compBlock%jBegor - 1
           ranges(nSub,3,2) = ranges(nn,3,2) + compBlock%kBegor - 1

         endif
       enddo

       end subroutine splitBlock

!========================================================================

       subroutine split2block(nSub, n1, n2, nTarget, ranges, &
                              viscousDir)
!
!      ******************************************************************
!      *                                                                *
!      * split2block splits the block stored in ranges(n1,:,:) into     *
!      * two. The new blocks are stored in ranges(n1,:,:) and           *
!      * ranges(n2,:,:), where the n2 will store the block with the     *
!      * number of cells closest to nTarget.
!      *                                                                *
!      ******************************************************************
!
       use inputIteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nSub, n1, n2, nTarget

       integer(kind=intType), dimension(nSub,3,2), intent(inout) :: ranges

       logical, dimension(3), intent(in) :: viscousDir
!
!      Local variables.
!
       integer(kind=intType) :: nMG, mm, ii, jj, i, iiOpt, iOpt
       integer(kind=intType) :: nCellOpt

       integer(kind=intType), dimension(3) :: nc, mc, c, nf, prefDir
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the maximum of the number of mg levels needed in the
       ! cycle and the mg start level. Store this value in nMG.

       nMG = max(nMGLevels, mgStartlevel)

       ! Store the number of cells of the original subblock in nc.

       nc(1) = ranges(n1,1,2) - ranges(n1,1,1)
       nc(2) = ranges(n1,2,2) - ranges(n1,2,1)
       nc(3) = ranges(n1,3,2) - ranges(n1,3,1)

       ! Determine its multigrid capabilities in the three directions.

       loopDir: do mm=1,3

         ! Initialize nf(mm), number of mg levels, to 1 and ii to 2.
         ! nf is used as a temporary buffer.

         nf(mm) = 1
         ii     = 2

         ! Loop to determine the number of mg levels.

         do
           ! Conditions to exit the loop. The first condition is the
           ! truncation to the maximum number of mg levels needed by the
           ! iterative algorithm of the solver. The second condition is
           ! simply that the number of cells do not allow more standard
           ! mg levels.

           if(nf(mm) == nMG) exit
           if(mod(nc(mm),ii) /= 0) exit

           ! Update nf(mm) and multiply ii by 2 to test for the
           ! next level.

           nf(mm) = nf(mm) +1
           ii     = 2*ii
         enddo

         ! Determine the constants c(mm) and mc(mm), such that
         ! nc(mm) = c(mm)*mc(mm). Mc(mm) is the number of cells in every
         ! direction per supercell. A supercell is the smallest possible
         ! block able to do multigridding, i.e. mc(mm) = 2**(nf(mm)-1)
         ! if the block allows for nMG multigrid levels.

         c(mm)  = 2*nc(mm)/ii
         mc(mm) = nc(mm)/c(mm)

       enddo loopDir

       ! Test whether the multigrid requirements are not too restrictive
       ! for this subblock. If they are, relax them, if possible.

       if(nMG > 1 .and. max(c(1),c(2),c(3)) == 1) then
          c(1) =  c(1)*2;  c(2) =  c(2)*2;  c(3) =  c(3)*2
         mc(1) = mc(1)/2; mc(2) = mc(2)/2; mc(3) = mc(3)/2
       endif

       ! Determine the number of faces on a slice in the
       ! three directions.

       nf(1) = nc(2)*nc(3)
       nf(2) = nc(1)*nc(3)
       nf(3) = nc(1)*nc(2)

       ! Determine the preferred split direction; this the direction
       ! which results in the least amount of additional faces. This is
       ! important, because the number of flux evaluations is
       ! proportional to the number of faces.

       if(nf(1) < nf(2)) then
         if(nf(1) < nf(3)) then
           prefDir(1) = 1
           if(nf(2) < nf(3)) then
             prefDir(2) = 2
             prefDir(3) = 3
           else
             prefDir(2) = 3
             prefDir(3) = 2
           endif
         else
           prefDir(1) = 3
           prefDir(2) = 1
           prefDir(3) = 2
         endif
       else if(nf(2) < nf(3)) then
         prefDir(1) = 2
         if(nf(1) < nf(3)) then
           prefDir(2) = 1
           prefDir(3) = 3
         else
           prefDir(2) = 3
           prefDir(3) = 1
         endif
       else
         prefDir(1) = 3
         prefDir(2) = 2
         prefDir(3) = 1
       endif

       ! Take viscousDir into account to determine the preference
       ! direction. For implicit methods it may not be a good idea to
       ! split the block parallel to a viscous wall.

       ! Check the 1st preference direction. If it is viscous swap
       ! it with the best inviscid direction, if available.

       if( viscousDir(prefDir(1)) ) then
         if(.not. viscousDir(prefDir(2)) ) then
           mm         = prefDir(1)
           prefDir(1) = prefDir(2)
           prefDir(2) = mm
         else if(.not. viscousDir(prefDir(3)) ) then
           mm         = prefDir(1)
           prefDir(1) = prefDir(3)
           prefDir(3) = mm
         endif
       endif

       ! Check the second preference direction. If it is viscous
       ! try to swap it whith the third direction.

       if(viscousDir(prefDir(2)) .and. &
          .not. viscousDir(prefDir(3)) ) then
         mm         = prefDir(2)
         prefDir(2) = prefDir(3)
         prefDir(3) = mm
       endif

       ! Determine the splitting which is best. Initialize nCellOpt
       ! to a ridiculously high number.

       nCellOpt = 10*nc(1)*nc(2)*nc(3)

       ! Loop over the three directions.

       do mm=1,3

         ! Store the current direction considered in ii and determine
         ! the index i, which gives a number of cells as close as
         ! possible to the desired number.

         ii = prefDir(mm)
         i  = nint(real(nTarget,realType) &
            /      real(mc(ii)*nf(ii),realType),intType)
         i  = max(i,1_intType)

         ! Determine whether the corresponding number of cells is
         ! closer to the target number than the currently stored value.
         ! If so, store the settings for this splitting.

         jj = i*mc(ii)*nf(ii)
         if(abs(jj - nTarget) < abs(nCellOpt - nTarget)) then
           nCellOpt = jj
           iiOpt    = ii
           iOpt     = i
         endif

       enddo

       ! Copy the range of subblock n1 into n2 as initialization.

       ranges(n2,1,1) = ranges(n1,1,1); ranges(n2,1,2) = ranges(n1,1,2)
       ranges(n2,2,1) = ranges(n1,2,1); ranges(n2,2,2) = ranges(n1,2,2)
       ranges(n2,3,1) = ranges(n1,3,1); ranges(n2,3,2) = ranges(n1,3,2)

       ! Determine the nodal index where the splitting takes place.

       jj = iOpt*mc(iiOpt) + ranges(n1,iiOpt,1)

       ! Adapt the corresponding indices of the new subblocks n1 and
       ! n2, such that it corresponds to the new situation; n2 should
       ! contain the subblock, which contains a number of cells as
       ! close as possible to nTarget.

       ranges(n2,iiOpt,2) = jj
       ranges(n1,iiOpt,1) = jj

       end subroutine split2block
