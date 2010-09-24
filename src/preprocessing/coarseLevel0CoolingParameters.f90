!
!      ******************************************************************
!      *                                                                *
!      * File:          coarseLevel0CoolingParameters.f90               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-30-2005                                      *
!      * Last modified: 08-03-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine coarseLevel0CoolingParameters(level)
!
!      ******************************************************************
!      *                                                                *
!      * coarseLevel0CoolingParameters creates parameters for the       *
!      * level 0 cooling model on the coarse grid from the known values *
!      * of the fine grid.                                              *
!      *                                                                *
!      * This model has been developed by Pratt and Whitney and should  *
!      * not be given to third parties.                                 *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use block
       use coolingModelLevel0
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: levm1, nn, mm
       integer(kind=intType) :: blockID, indexDir
       integer(kind=intType) :: indSol, indNorm, indX1, indX2
       integer(kind=intType) :: jcBeg, jcEnd, icBeg, icEnd
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the finer grid level in levm1

       levm1 = level -1

       ! Loop over the number of cooling planes specified.

       coolingLoop: do nn=1,nPlanesLevel0CoolingModel

         ! Copy the number of subfaces and the parameters of the model
         ! from the finer grid level.

         level0Cooling(nn,level)%nSubfaces = level0Cooling(nn,levm1)%nSubfaces
         level0Cooling(nn,level)%mDotRatio = level0Cooling(nn,levm1)%mDotRatio
         level0Cooling(nn,level)%dpLog     = level0Cooling(nn,levm1)%dpLog
         level0Cooling(nn,level)%dTLog     = level0Cooling(nn,levm1)%dTLog

         ! Allocate the memory for the subface info.

         mm = level0Cooling(nn,level)%nSubfaces

         allocate(level0Cooling(nn,level)%blockID(mm),  &
                  level0Cooling(nn,level)%indexDir(mm), &
                  level0Cooling(nn,level)%indSol(mm),   &
                  level0Cooling(nn,level)%indNorm(mm),  &
                  level0Cooling(nn,level)%indX1(mm),    &
                  level0Cooling(nn,level)%indX2(mm),    &
                  level0Cooling(nn,level)%jcBeg(mm),    &
                  level0Cooling(nn,level)%jcEnd(mm),    &
                  level0Cooling(nn,level)%icBeg(mm),    &
                  level0Cooling(nn,level)%icEnd(mm),    &
                  stat=ierr)
         if(ierr /= 0)                                     &
           call terminate("coarseLevel0CoolingParameters", &
                          "Memory allocation failure for subface data")

         ! Loop over the number of subfaces for this cooling plane.

         subfaceLoop: do mm=1,level0Cooling(nn,level)%nSubfaces

           ! Easier storage of some of the variables.

           blockID  = level0Cooling(nn,levm1)%blockID(mm)
           indexDir = level0Cooling(nn,levm1)%indexDir(mm)
           indSol   = level0Cooling(nn,levm1)%indSol(mm)
           indX1    = level0Cooling(nn,levm1)%indX1(mm)
           indX2    = level0Cooling(nn,levm1)%indX2(mm)

           icBeg = level0Cooling(nn,levm1)%icBeg(mm)
           icEnd = level0Cooling(nn,levm1)%icEnd(mm)
           jcBeg = level0Cooling(nn,levm1)%jcBeg(mm)
           jcEnd = level0Cooling(nn,levm1)%jcEnd(mm)

           ! Determine the index direction and determine the subface
           ! parameters accordingly.

           select case (indexDir)
             case (iMin, iMax)

               ! Index direction is the I-direction. Use the multigrid
               ! info to determine the coarse grid cell index downstream
               ! of the injection plane. From that info indNorm, indX1
               ! and indX2 can be extracted.

               indSol = flowDoms(blockID,levm1,1)%mgICoarse(indSol,1)

               if(indX2 > indX1) then   ! Positive index is downstream
                 indNorm = indSol - 1
                 indX1   = indNorm
                 indX2   = indX1 + 1
               else                     ! Negative index is downstream
                 indNorm = indSol
                 indX1   = indNorm
                 indX2   = indX1 - 1
               endif

               ! Determine the coarse grid cell range on the plane.

               icBeg = flowDoms(blockID,levm1,1)%mgJCoarse(icBeg,1)
               icEnd = flowDoms(blockID,levm1,1)%mgJCoarse(icEnd,1)
               jcBeg = flowDoms(blockID,levm1,1)%mgKCoarse(jcBeg,1)
               jcEnd = flowDoms(blockID,levm1,1)%mgKCoarse(jcEnd,1)

             !===========================================================

             case (jMin, jMax)

               ! Index direction is the J-direction. Use the multigrid
               ! info to determine the coarse grid cell index downstream
               ! of the injection plane. From that info indNorm, indX1
               ! and indX2 can be extracted.

               indSol = flowDoms(blockID,levm1,1)%mgJCoarse(indSol,1)

               if(indX2 > indX1) then   ! Positive index is downstream
                 indNorm = indSol - 1
                 indX1   = indNorm
                 indX2   = indX1 + 1
               else                     ! Negative index is downstream
                 indNorm = indSol
                 indX1   = indNorm
                 indX2   = indX1 - 1
               endif

               ! Determine the coarse grid cell range on the plane.

               icBeg = flowDoms(blockID,levm1,1)%mgICoarse(icBeg,1)
               icEnd = flowDoms(blockID,levm1,1)%mgICoarse(icEnd,1)
               jcBeg = flowDoms(blockID,levm1,1)%mgKCoarse(jcBeg,1)
               jcEnd = flowDoms(blockID,levm1,1)%mgKCoarse(jcEnd,1)

             !===========================================================

             case (kMin, kMax)

               ! Index direction is the K-direction. Use the multigrid
               ! info to determine the coarse grid cell index downstream
               ! of the injection plane. From that info indNorm, indX1
               ! and indX2 can be extracted.

               indSol = flowDoms(blockID,levm1,1)%mgKCoarse(indSol,1)

               if(indX2 > indX1) then   ! Positive index is downstream
                 indNorm = indSol - 1
                 indX1   = indNorm
                 indX2   = indX1 + 1
               else                     ! Negative index is downstream
                 indNorm = indSol
                 indX1   = indNorm
                 indX2   = indX1 - 1
               endif

               ! Determine the coarse grid cell range on the plane.

               icBeg = flowDoms(blockID,levm1,1)%mgICoarse(icBeg,1)
               icEnd = flowDoms(blockID,levm1,1)%mgICoarse(icEnd,1)
               jcBeg = flowDoms(blockID,levm1,1)%mgJCoarse(jcBeg,1)
               jcEnd = flowDoms(blockID,levm1,1)%mgJCoarse(jcEnd,1)

           end select

           ! Store the info in the coarse grid derived data type.

           level0Cooling(nn,level)%blockID(mm)  = blockID
           level0Cooling(nn,level)%indexDir(mm) = indexDir

           level0Cooling(nn,level)%indSol(mm)  = indSol
           level0Cooling(nn,level)%indNorm(mm) = indNorm
           level0Cooling(nn,level)%indX1(mm)   = indX1
           level0Cooling(nn,level)%indX2(mm)   = indX2

           level0Cooling(nn,level)%jcBeg(mm) = jcBeg
           level0Cooling(nn,level)%jcEnd(mm) = jcEnd
           level0Cooling(nn,level)%icBeg(mm) = icBeg
           level0Cooling(nn,level)%icEnd(mm) = icEnd

         enddo subfaceLoop

       enddo coolingLoop

       end subroutine coarseLevel0CoolingParameters
