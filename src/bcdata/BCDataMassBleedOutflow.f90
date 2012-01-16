!
!      ******************************************************************
!      *                                                                *
!      * File:          BCDataMassBleedOutflow.f90                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-06-2007                                      *
!      * Last modified: 11-30-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine BCDataMassBleedOutflow(initPressure, useInternalOnly)
!
!      ******************************************************************
!      *                                                                *
!      * BCDataMassBleedOutflow adapts the prescribed static pressure   *
!      * for bleed outflow boundaries, such that the mass flow leaving  *
!      * the domain is closer to the prescribed value. If needed, the   *
!      * prescribed pressure is initialized to the internal pressure.   *
!      *                                                                *
!      ******************************************************************
!
       use bleedFlows
       use blockPointers
       use BCTypes
       use inputIteration
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: initPressure, useInternalOnly
!
!      Local variables.
!
       integer(kind=intType) :: level, nn, mm, sps, i, j
       integer(kind=intType) :: nLevels, ii, if1, if2, jf1, jf2
       integer(kind=intType) :: iiMax, jjMax

       integer(kind=intType), dimension(:,:), pointer :: iFine, jFine

       real(kind=realType) :: fact

       real(kind=realType), dimension(:,:),  pointer :: pp2, ps, psFine
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if no outflow bleeds are present.

       if(nOutflowBleeds == 0) return
!
!      ******************************************************************
!      *                                                                *
!      * Initialize the prescribed pressure, if desired.                *
!      * This must be done on the multigrid level currentLevel.         *
!      *                                                                *
!      ******************************************************************
!
       if( initPressure ) then

         ! Loop over the time instances and domains.

         do sps=1,nTimeIntervalsSpectral
           do nn=1,nDom

             ! Set the pointers to this block on currentLevel.

             call setPointers(nn,currentLevel,sps)

             ! Loop over the boundary subfaces and check for outflow
             ! mass bleeds.

             do mm=1,nBocos
               if(BCType(mm) == MassBleedOutflow) then

                 ! Determine the block face on which the subface is
                 ! located and set the pointer for pp2 accordingly.

                 select case (BCFaceID(mm))
                   case (iMin)
                     pp2 => p(2,1:,1:)
                   case (iMax)
                     pp2 => p(il,1:,1:)
                   case (jMin)
                     pp2 => p(1:,2,1:)
                   case (jMax)
                     pp2 => p(1:,jl,1:)
                   case (kMin)
                     pp2 => p(1:,1:,2)
                   case (kMax)
                     pp2 => p(1:,1:,kl)
                 end select

                 ! Loop over the range of the subface and copy the data
                 ! for the static pressure.

                 ps => BCData(mm)%ps

                 do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                     ps(i,j) = pp2(i,j)
                   enddo
                 enddo

               endif
             enddo
           enddo
         enddo
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Computation of the prescribed pressures.                       *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of time instances.

       do sps=1,nTimeIntervalsSpectral

         ! Compute the global parameters for the bleeds.

         call bleedFlowParameters(sps, useInternalOnly)

         ! Loop over the blocks.

         do nn=1,nDom

           ! Set the pointers to this block for currentLevel.

           call setPointers(nn,currentLevel,sps)

           ! Loop over the number of boundary subfaces and check for an
           ! outflow bleed region.

           do mm=1,nBocos
             if(BCType(mm) == MassBleedOutflow) then

               ! Determine the relaxation factor for the pressure.

               ii   = groupNum(mm)
               fact = 1.0_realType                               &
                    + relaxBleeds*(outflowBleeds(ii)%curMassFlux &
                    -              outflowBleeds(ii)%massFlux)   &
                    /              outflowBleeds(ii)%massFlux

               ! Make sure that the factor is within reasonable limits.

               fact = max(0.9_realType,min(fact,1.1_realType))
               print 101, fact, outflowBleeds(ii)%massFlux, &
                          outflowBleeds(ii)%curMassFlux
 101           format("# Outflowbleeds: ", f9.6,2(1x,e12.5))

               ! Multiply the prescribed pressure by fact.

               do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                 do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                   BCData(mm)%ps(i,j) = fact*BCData(mm)%ps(i,j)
                 enddo
               enddo

             endif
           enddo

         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Modify the prescribed pressure on the coarser grid levels.     *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of grid levels.

       nLevels = ubound(flowDoms,2)

       do level=(currentLevel+1),nLevels

         ! Loop over the number of time instances and local blocks.

         do sps=1,nTimeIntervalsSpectral
           do nn=1,nDom

             ! Set the pointers to this block.

             call setPointers(nn,level,sps)

             ! Loop over the number of boundary subfaces and check for
             ! an outflow bleed region.

             do mm=1,nBocos
               if(BCType(mm) == MassBleedOutflow) then

                 ! Set the pointer for the prescribed pressure for this
                 ! level as well as for the finer grid level.

                 ps => BCData(mm)%ps
                 psFine => flowDoms(nn,level-1,sps)%BCData(mm)%ps

                 ! Determine the block face on which the subface is
                 ! located and set some multigrid variables accordingly.

                 select case (BCFaceID(mm))

                   case (iMin,iMax)
                     iiMax = jl; jjMax = kl
                     iFine => mgJFine; jFine => mgKFine

                   case (jMin,jMax)
                     iiMax = il; jjMax = kl
                     iFine => mgIFine; jFine => mgKFine

                   case (kMin,kMax)
                     iiMax = il; jjMax = jl
                     iFine => mgIFine; jFine => mgJFine

                 end select

                 ! Loop over the j-direction of this subface.

                 do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd

                   ! Determine the two children in this direction.
                   ! Take care of the halo's, as this info is only
                   ! available for owned cells.

                   if(j < 2) then
                     jf1 = 1; jf2 = 1
                   else if(j > jjMax) then
                     jf1 = jFine(jjMax,2) +1; jf2 = jf1
                   else
                     jf1 = jFine(j,1); jf2 = jFine(j,2)
                   endif

                   ! Loop in the i-direction.

                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                     ! Determine the two children in this direction.
                     ! Same story as in j-direction.

                     if(i < 2) then
                       if1 = 1; if2 = 1
                     else if(i > iiMax) then
                       if1 = iFine(iiMax,2) +1; if2 = if1
                     else
                       if1 = iFine(i,1); if2 = iFine(i,2)
                     endif

                     ! Compute the pressure.

                     ps(i,j) = fourth*(psFine(if1,jf1) + psFine(if2,jf1) &
                             +         psFine(if1,jf2) + psFine(if2,jf2))
                   enddo
                 enddo

               endif
             enddo
           enddo
         enddo
       enddo

       end subroutine BCDataMassBleedOutflow
