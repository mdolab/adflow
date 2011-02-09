!
!      ******************************************************************
!      *                                                                *
!      * File:          getIndices.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-26-2005                                      *
!      * Last modified: 10-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine getBegIndices(offset, sizeEntity, blockDim, &
                                iBeg, jBeg, kBeg, lBeg)
!
!      ******************************************************************
!      *                                                                *
!      * getBegIndices determines for the given offset in bytes from    *
!      * the beginning of the block, the size of an entity and the      *
!      * block dimensions the starting indices in i, j and k direction  *
!      * as well as the l-direction. This last dimension is not a real  *
!      * block dimension, but the number of physical variables          *
!      * read/written. Therefore it is treated slightly differently.    *
!      *                                                                *
!      ******************************************************************
!
       use IOModule
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(inout) :: offset
       integer(kind=intType), intent(in)    :: sizeEntity
       integer(kind=intType), dimension(3), intent(in) :: blockDim

       integer(kind=intType), intent(out) :: iBeg, jBeg, kBeg, lBeg
!
!      Local variables.
!
       integer(kind=intType) :: ii, nEntities
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the total number of entities in the block.

       nEntities = blockDim(1) * blockDim(2) * blockDim(3)

       ! Determine the situation we are having here.

       if(offset <= 0) then

         ! Lower bound corresponds to the very beginning of the block.
         ! Set the begin indices to 1.

         iBeg = 1; jBeg = 1; kBeg = 1; lBeg = 1

       else if(offset >= P3D_nVar*sizeEntity*nEntities) then

         ! Lower bound is at a position after the maximum indices of this
         ! block. This usually only happens when i-blanking is present.

         iBeg = 1; jBeg = 1; kBeg = 1; lBeg = P3D_nVar + 1

       else

         ! Lower bound is somewhere in the middle of the record.
         ! Determine the l-index.

         ii     = sizeEntity*nEntities
         lBeg   = offset/ii + 1
         offset = offset - (lBeg-1)*ii

         ! Determine the k-index.

         ii     = sizeEntity*blockDim(1)*blockDim(2)
         kBeg   = offset/ii + 1
         offset = offset - (kBeg-1)*ii

         ! Determine the j-index and the i-index.

         ii     = sizeEntity*blockDim(1)
         jBeg   = offset/ii + 1
         offset = offset - (jBeg-1)*ii
         iBeg   = offset/sizeEntity + 1

       endif

       ! Correct the indices for a possible index shift, because owned
       ! cell data starts at index 2. Only needed for cell centered
       ! data without halo's.

       if(P3D_DataStorage == cellDataNoHalo) then
         iBeg = iBeg + 1; jBeg = jBeg + 1; kBeg = kBeg + 1
       endif

       end subroutine getBegIndices

!      ==================================================================

       subroutine getEndIndices(offset, sizeEntity, blockDim, &
                                iEnd, jEnd, kEnd, lEnd)
!
!      ******************************************************************
!      *                                                                *
!      * getEndIndices determines for the given offset in bytes from    *
!      * the end of the block, the size of an entity and the block      *
!      * dimensions the ending indices in i, j and k direction as well  *
!      * as the l-direction. This last dimension is not a real block    *
!      * dimension, but the number of physical variables read/written.  *
!      * Therefore it is treated slightly differently.                  *
!      *                                                                *
!      ******************************************************************
!
       use IOModule
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(inout) :: offset
       integer(kind=intType), intent(in)    :: sizeEntity
       integer(kind=intType), dimension(3), intent(in) :: blockDim

       integer(kind=intType), intent(out) :: iEnd, jEnd, kEnd, lEnd
!
!      Local variables.
!
       integer(kind=intType) :: ii
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the situation we are having here.

       if(offset <= 0) then

         ! Upper bound corresponds to the end of the global block.
         ! Set the indices accordingly.

         iEnd = blockDim(1)
         jEnd = blockDim(2)
         kEnd = blockDim(3)
         lEnd = P3D_nVar

       else

         ! The upper bound corresponds to a position somewhere in
         ! the middle of the array. Find out the indices by
         ! subtracting then from the back. First the l-index.

         ii     = sizeEntity*blockDim(1)*blockDim(2)*blockDim(3)
         lEnd   = P3D_nVar - offset/ii
         offset = offset - (P3D_nVar - lEnd)*ii

         ! The k-index.

         ii     = sizeEntity*blockDim(1)*blockDim(2)
         kEnd   = blockDim(3) - offset/ii
         offset = offset - (blockDim(3) - kEnd)*ii

         ! The j and i indices.

         ii     = sizeEntity*blockDim(1)
         jEnd   = blockDim(2) - offset/ii
         offset = offset - (blockDim(2) - jEnd)*ii
         iEnd   = blockDim(1) - offset/sizeEntity

       endif

       ! Correct the indices for a possible index shift, because
       ! owned cell data starts at index 2. Only needed for
       ! cell centered data without halo's.

       if(P3D_DataStorage == cellDataNoHalo) then
         iEnd = iEnd + 1; jEnd = jEnd + 1; kEnd = kEnd + 1
       endif

       end subroutine getEndIndices

!      ==================================================================

       subroutine getMinMaxIndices(iBegOr, jBegOr, kBegOr,             &
                                   iEndOr, jEndOr, kEndOr, il, jl, kl, &
                                   dataIsRead, includeHalos,           &
                                   iMin, jMin, kMin, iMax, jMax, kMax)
!
!      ******************************************************************
!      *                                                                *
!      * getMinMaxIndices determines the range iMin-iMax, jMin-jMax,    *
!      * kMin-kMax of the subblock that is active in the IO.            *
!      *                                                                *
!      ******************************************************************
!
       use IOModule
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: iBegOr, jBegOr, kBegOr
       integer(kind=intType), intent(in) :: iEndOr, jEndOr, kEndOr
       integer(kind=intType), intent(in) :: il, jl, kl

       logical, intent(in) :: dataIsRead, includeHalos

       integer(kind=intType), intent(out) :: iMin, jMin, kMin
       integer(kind=intType), intent(out) :: iMax, jMax, kMax
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the situation we are having here.

       select case (P3D_DataStorage)
         case (cellDataNoHalo)
           iMin = iBegOr + 1
           iMax = iEndOr
           jMin = jBegOr + 1
           jMax = jEndOr
           kMin = kBegOr + 1
           kMax = kEndOr

         !===============================================================

         case (nodeData)
           iMin = iBegOr
           iMax = iEndOr
           jMin = jBegOr
           jMax = jEndOr
           kMin = kBegOr
           kMax = kEndOr

           if(.not. dataIsRead) then
             if(iMin > 1) iMin = iMin + 1
             if(jMin > 1) jMin = jMin + 1
             if(kMin > 1) kMin = kMin + 1
           endif

         !===============================================================

         case (cellDataPlusHalo)

           ! Make a distinction between reading and writing.

           if( dataIsRead ) then

             ! Make the distinction whether or not to read the halo data.

             if( includeHalos ) then

               ! Halo data should be read.

               iMin = iBegOr
               iMax = iEndOr + 1
               jMin = jBegOr
               jMax = jEndOr + 1
               kMin = kBegOr
               kMax = kEndOr + 1

             else

               ! Halo data should not be read.

               iMin = iBegOr + 1
               iMax = iEndOr
               jMin = jBegOr + 1
               jMax = jEndOr
               kMin = kBegOr + 1
               kMax = kEndOr

             endif

           else

             ! Data is written. Only write the halo data if my
             ! boundaries coincide with the boundaries of the
             ! original block.

             iMin = iBegOr + 1; if(iMin == 2)  iMin = 1
             iMax = iEndOr;     if(iMax == il) iMax = iMax + 1
             jMin = jBegOr + 1; if(jMin == 2)  jMin = 1
             jMax = jEndOr;     if(jMax == jl) jMax = jMax + 1
             kMin = kBegOr + 1; if(kMin == 2)  kMin = 1
             kMax = kEndOr;     if(kMax == kl) kMax = kMax + 1

           endif

       end select

       end subroutine getMinMaxIndices
