!
!      ******************************************************************
!      *                                                                *
!      * File:          orphanAverage.f90                               *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 09-14-2005                                      *
!      * Last modified: 10-16-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine orphanAverage(wstart, wend, calcPressure, calcGamma, &
                                calcLamVis, calcEddyVis)
!
!      ******************************************************************
!      *                                                                *
!      * orphanAverage uses the neighboring cells of an overset orphan  *
!      * to set the flow state for the orphan cell by a simple average. *
!      * This routine operates on the block given by the block pointers *
!      * so it is assumed they are set.                                 *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use flowVarRefState
       use inputPhysics
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: wstart, wend

       logical, intent(in) :: calcPressure, calcGamma, calcLamVis
       logical, intent(in) :: calcEddyVis
!
!      Local variables.
!
       integer(kind=intType) :: oi, oj, ok, ni, nj, nk, i, l, m, n, nAvg

       integer(kind=intType), dimension(3) :: del

       real(kind=realType) :: nAvgReal
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if there are no orphans for this block.

       if (nOrphans == 0) return

       ! Change the iblanks on the fringe so that the average can
       ! distinguish between hole nad fringe neighbors.

       call changeIblanks(.false., 1_intType)

       ! Loop over the number of orphans.

       orphanLoop: do n = 1,nOrphans

         ! Store the orphan indices easier.

         m = nCellsOverset + n

         oi = ibndry(1,m)
         oj = ibndry(2,m)
         ok = ibndry(3,m)

         ! Initialize the number of neighbors used to 0 and also set
         ! the flow variables to zero such that an average can be
         ! accumulated below.

         nAvg = 0

         do l=wstart,wend
           w(oi,oj,ok,l) = zero
         end do

         if (calcPressure) p(oi,oj,ok)     = zero
         if (calcGamma)    gamma(oi,oj,ok) = zero
         if (calcLamVis)   rlv(oi,oj,ok)   = zero
         if (calcEddyVis)  rev(oi,oj,ok)   = zero

         ! Loop over the 3 coordinate directions, and for both the
         ! positive and negative direction set the delta vector to be a
         ! unit vector in that direction.

         directionLoop: do m = 1,3
           plusMinusLoop: do i = -1,1,2

             del    = 0
             del(m) = i

             ! Compute the neighbor indices and skip if it is outside the
             ! boundaries of the block.

             ni = oi + del(1)
             nj = oj + del(2)
             nk = ok + del(3)

             if (ni < 0 .or. ni > ib .or. &
                 nj < 0 .or. nj > jb .or. &
                 nk < 0 .or. nk > kb) cycle

             ! If the neighboring iblank value indicates a cell that is
             ! part of either the field, fringe, or boundary condition,
             ! then use its flow state in the average.

             if (iblank(ni,nj,nk) >= 1) then

               ! Update the number of neighbors used in the average and 
               ! compute the flow variables for the given range.

               nAvg = nAvg + 1

               do l=wstart,wend
                 w(oi,oj,ok,l) = w(oi,oj,ok,l) + w(ni,nj,nk,l)
               end do

               ! Check if the pressure, specific heat ratio, laminar
               ! viscosity, and/or eddy viscosity needs to be computed.

               if (calcPressure) &
                     p(oi,oj,ok) = p(oi,oj,ok) + p(ni,nj,nk)
               if (calcGamma) &
                     gamma(oi,oj,ok) = gamma(oi,oj,ok) + gamma(ni,nj,nk)
               if (calcLamVis) &
                     rlv(oi,oj,ok) = rlv(oi,oj,ok) + rlv(ni,nj,nk)
               if (calcEddyVis) &
                     rev(oi,oj,ok) = rev(oi,oj,ok) + rev(ni,nj,nk)

             end if

           end do plusMinusLoop
         end do directionLoop

         ! Check to make sure that at least 1 suitable neighbeor was
         ! found to use in the average. 

         checkNoNeighbors: if (nAvg > 0) then

           ! Divide each of the variables being computed by the number
           ! of neighbors used in the average.

           nAvgReal = real(nAvg, realType)

           ! Average the flow variables for the given range.

           do l=wstart,wend
             w(oi,oj,ok,l) = w(oi,oj,ok,l)/nAvgReal
           end do

           ! Check if the pressure, specific heat ratio, laminar
           ! viscosity, and/or eddy viscosity needs to be averaged.

           if (calcPressure) p(oi,oj,ok)     = p(oi,oj,ok)/nAvgReal
           if (calcGamma)    gamma(oi,oj,ok) = gamma(oi,oj,ok)/nAvgReal
           if (calcLamVis)   rlv(oi,oj,ok)   = rlv(oi,oj,ok)/nAvgReal
           if (calcEddyVis)  rev(oi,oj,ok)   = rev(oi,oj,ok)/nAvgReal

         else checkNoNeighbors

           ! No suitable neighbors were found in order to compute an
           ! average. Set the variables back to the the freestream.

           do l=wstart,wend
             w(oi,oj,ok,l) = wInf(l)
           end do

           if (calcPressure) p(oi,oj,ok)     = pInfCorr
           if (calcGamma)    gamma(oi,oj,ok) = gammaInf
           if (calcLamVis)   rlv(oi,oj,ok)   = muInf
           if (calcEddyVis)  rev(oi,oj,ok)   = eddyVisInfRatio*muInf

         end if checkNoNeighbors

       end do orphanLoop

       ! Change the iblanks on the fringe back to 0.

       call changeIblanks(.false., 0_intType)

       end subroutine orphanAverage
