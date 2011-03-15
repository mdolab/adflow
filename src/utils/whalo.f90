!
!      ******************************************************************
!      *                                                                *
!      * File:          whalo.f90                                       *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 03-07-2003                                      *
!      * Last modified: 11-09-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine whalo1(level, start, end, commPressure, commGamma, &
                         commViscous)
!
!      ******************************************************************
!      *                                                                *
!      * whalo1 exchanges all the 1st level internal halo's for the     *
!      * cell centered variables.                                       *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use commSliding
       use communication
       use flowVarRefState
       use inputPhysics
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, start, end
       logical, intent(in) :: commPressure, commGamma, commViscous
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, ll

       logical :: correctForK, commLamVis, commEddyVis, commVarGamma
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the logicals whether or not to communicate the viscosities.

       commLamVis = .false.
       if(viscous .and. commViscous) commLamVis = .true.

       commEddyVis = .false.
       if(eddyModel .and. commViscous) commEddyVis = .true.

       ! Set the logical whether or not to communicate gamma.

       commVarGamma = .false.
       if(commGamma .and. (cpModel == cpTempCurveFits)) &
         commVarGamma = .true.

       ! Exchange the 1 to 1 matching 1st level cell halo's.

       call whalo1to1(level, start, end, commPressure, commVarGamma, &
                      commLamVis, commEddyVis, commPatternCell_1st,  &
                      internalCell_1st)

       ! Exchange the sliding mesh 1st level cell halo's.

       mm = ubound(commSlidingCell_1st,1)
       call whaloSliding(level, start, end, commPressure,       &
                         commVarGamma, commLamVis, commEddyVis, &
                         commSlidingCell_1st, intSlidingCell_1st, mm)

       ! Exchange the mixing plane 1st level cell halo's.

       call whaloMixing(level, start, end, commPressure, commVarGamma, &
                        commLamVis, commEddyVis, 1_intType)

       ! Exchange the overset boundary cell data.

       call wOverset(level, start, end, commPressure,       &
                     commVarGamma, commLamVis, commEddyVis, &
                     commPatternOverset, internalOverset, mm)

       ! Average any overset orphans.

       do ll=1,nTimeIntervalsSpectral
         do nn=1,nDom
           call setPointers(nn,level,ll)
           call orphanAverage(start, end, commPressure, commGamma, &
                              commLamVis, commEddyVis)
         end do
       end do

       ! If both the pressure and the total energy has been communicated
       ! compute the energy again. The reason is that both values are
       ! interpolated and consequently the values are not consistent.
       ! The energy depends quadratically on the velocity.

       bothPAndE: if(commPressure .and. start <= irhoE .and. &
                     end >= irhoE) then

         ! First determine whether or not the total energy must be
         ! corrected for the presence of the turbulent kinetic energy.

         if( kPresent ) then
           if((level <= groundLevel) .or. turbCoupled) then
             correctForK = .true.
           else
             correctForK = .false.
           endif
         else
           correctForK = .false.
         endif

         ! Loop over the blocks to find the sliding mesh subfaces.
         ! Use is made of the fact the boundary conditions are identical
         ! for all spectral solutions. So that loop can be inside the
         ! test for the sliding mesh subface.

         domains: do nn=1,nDom
           do mm=1,flowDoms(nn,level,1)%nBocos
             if(flowDoms(nn,level,1)%BCType(mm) == slidingInterface) then

               ! Loop over the number of spectral solutions.

               do ll=1,nTimeIntervalsSpectral

                 ! Set the pointers for this block and compute the energy
                 ! for the halo cells of this sliding interface subface.

                 call setPointers(nn,level,ll)
                 call computeEtot(icBeg(mm), icEnd(mm), &
                                  jcBeg(mm), jcEnd(mm), &
                                  kcBeg(mm), kcEnd(mm), correctForK)
               enddo
             endif
           enddo

           ! Now treat the overset boundary.

           do ll=1,nTimeIntervalsSpectral
             call setPointers(nn,level,ll)
             call computeEtotBndryList(correctForK)
           enddo

         enddo domains

       endif bothPAndE

       end subroutine whalo1

!      ==================================================================

       subroutine whalo2(level, start, end, commPressure, commGamma, &
                         commViscous)
!
!      ******************************************************************
!      *                                                                *
!      * whalo2 exchanges all the 2nd level internal halo's for the     *
!      * cell centered variables.                                       *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use commSliding
       use communication
       use flowVarRefState
       use inputPhysics
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, start, end
       logical, intent(in) :: commPressure, commGamma, commViscous
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, ll
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

       logical :: correctForK, commLamVis, commEddyVis, commVarGamma
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the logicals whether or not to communicate the viscosities.

       commLamVis = .false.
       if(viscous .and. commViscous) commLamVis = .true.

       commEddyVis = .false.
       if(eddyModel .and. commViscous) commEddyVis = .true.

       ! Set the logical whether or not to communicate gamma.

       commVarGamma = .false.
       if(commGamma .and. (cpModel == cpTempCurveFits)) &
         commVarGamma = .true.

       ! Exchange the 1 to 1 matching 2nd level cell halo's.

       call whalo1to1(level, start, end, commPressure, commVarGamma, &
                      commLamVis, commEddyVis, commPatternCell_2nd,  &
                      internalCell_2nd)

       ! Exchange the sliding mesh 2nd level cell halo's.

       mm = ubound(commSlidingCell_2nd,1)
       call whaloSliding(level, start, end, commPressure,       &
                         commVarGamma, commLamVis, commEddyVis, &
                         commSlidingCell_2nd, intSlidingCell_2nd, mm)

       ! Exchange the mixing plane 2nd level cell halo's.

       call whaloMixing(level, start, end, commPressure, commVarGamma, &
                        commLamVis, commEddyVis, 2_intType)

       ! Exchange the overset boundary cell data.

       call wOverset(level, start, end, commPressure,       &
                     commVarGamma, commLamVis, commEddyVis, &
                     commPatternOverset, internalOverset, mm)

       ! Average any overset orphans.

       do ll=1,nTimeIntervalsSpectral
         do nn=1,nDom
           call setPointers(nn,level,ll)
           call orphanAverage(start, end, commPressure, commGamma, &
                              commLamVis, commEddyVis)
         end do
       end do

       ! If both the pressure and the total energy has been communicated
       ! compute the energy again. The reason is that both values are
       ! interpolated and consequently the values are not consistent.
       ! The energy depends quadratically on the velocity.

       bothPAndE: if(commPressure .and. start <= irhoE .and. &
                     end >= irhoE) then

         ! First determine whether or not the total energy must be
         ! corrected for the presence of the turbulent kinetic energy.

         if( kPresent ) then
           if((level <= groundLevel) .or. turbCoupled) then
             correctForK = .true.
           else
             correctForK = .false.
           endif
         else
           correctForK = .false.
         endif

         ! Loop over the blocks to find the sliding mesh subfaces.
         ! Use is made of the fact the boundary conditions are identical
         ! for all spectral solutions. So that loop can be inside the
         ! test for the sliding mesh subface.

         domains: do nn=1,nDom
           do mm=1,flowDoms(nn,level,1)%nBocos
             if(flowDoms(nn,level,1)%BCType(mm) == slidingInterface) then

               ! Loop over the number of spectral solutions.

               do ll=1,nTimeIntervalsSpectral

                 ! Set the pointers for this block.

                 call setPointers(nn,level,ll)

                 ! Set the range, depending on the block face on which
                 ! the subface is located.

                 select case (bcFaceID(mm))
                   case (iMin)
                     iBeg = 0;         jBeg = jcBeg(mm); kBeg = kcBeg(mm)
                     iEnd = 1;         jEnd = jcEnd(mm); kEnd = kcEnd(mm)
                   case (iMax)
                     iBeg = ie;        jBeg = jcBeg(mm); kBeg = kcBeg(mm)
                     iEnd = ib;        jEnd = jcEnd(mm); kEnd = kcEnd(mm)
                   case (jMin)
                     iBeg = icBeg(mm); jBeg = 0;         kBeg = kcBeg(mm)
                     iEnd = icEnd(mm); jEnd = 1;         kEnd = kcEnd(mm)
                   case (jMax)
                     iBeg = icBeg(mm); jBeg = je;        kBeg = kcBeg(mm)
                     iEnd = icEnd(mm); jEnd = jb;        kEnd = kcEnd(mm)
                   case (kMin)
                     iBeg = icBeg(mm); jBeg = jcBeg(mm); kBeg = 0
                     iEnd = icEnd(mm); jEnd = jcEnd(mm); kEnd = 1
                   case (kMax)
                     iBeg = icBeg(mm); jBeg = jcBeg(mm); kBeg = ke
                     iEnd = icEnd(mm); jEnd = jcEnd(mm); kEnd = kb
                 end select

                 ! Compute the total energy for the sliding mesh halo's

                 call computeEtot(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd, &
                                  correctForK)
               enddo
             endif
           enddo

           ! Now treat the overset boundary.

           do ll=1,nTimeIntervalsSpectral
             call setPointers(nn,level,ll)
             call computeEtotBndryList(correctForK)
           enddo

         enddo domains

       endif bothPAndE

       end subroutine whalo2
