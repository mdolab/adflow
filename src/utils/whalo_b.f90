!
!      ******************************************************************
!      *                                                                *
!      * File:          whalo_b.f90                                     *
!      * Author:        Gaetan K. W. Kenway                             *
!      * Starting date: 01-22-2015                                      *
!      * Last modified: 01-22-2015                                      *
!      *                                                                *
!      ******************************************************************

       subroutine whalo2_b(level, start, end, commPressure, commGamma, &
                         commViscous)
!
!      ******************************************************************
!      *                                                                *
!      * whalo2_b exchanges all the 2nd level internal halo's for the   *
!      * cell centered variables IN REVERSE MODE                        *
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

       
       mm = ubound(commPatternOverset, 1)
       call wOverset_b(level, start, end, commPressure, commVarGamma, &
            commLamVis, commEddyVis, commPatternOverset, internalOverset, mm)

       ! Exchange the 1 to 1 matching 2nd level cell halo's.

       call whalo1to1_b(level, start, end, commPressure, commVarGamma, &
            commLamVis, commEddyVis, commPatternCell_2nd,  &
            internalCell_2nd)

       ! NOTE: Only the 1to1 halo exchange and overset is done. whalosliding,
       ! whalomixing, orphanAverage and PandE corrections
       ! calculation are NOT implementent. 

     end subroutine whalo2_b
