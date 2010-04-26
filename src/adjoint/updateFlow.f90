!
!     ******************************************************************
!     *                                                                *
!     * File:          updateFlow.f90                                  *
!     * Authors:       C.A(Sandy) Mader                                *
!     * Starting date: 24-02-2009                                      *
!     * Last modified: 24-02-2009                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine updateFlow
!
!     ******************************************************************
!     *                                                                *
!     *  reruns the initialization routines to update AOA and other    *
!     *  flow variables after a design change                          *
!     *                                                                *
!     ******************************************************************
!
       
        use inputTimeSpectral ! spaceDiscr,nTimeIntervalsSpectral
        use section
        use inputPhysics  !
        use monitor
        use iteration
        implicit none


        integer(kind=intType) ::mm,nnn
        real(kind=realType), dimension(nSections) :: t


        !begin execution
        !print *,'in update flow'
        call referenceState
        !print *,'reference state'
        call setFlowInfinityState
        !print *,'flow infinity state'
        groundlevel = 1
        do mm=1,nTimeIntervalsSpectral
           
           ! Compute the time, which corresponds to this spectral solution.
           ! For steady and unsteady mode this is simply the restart time;
           ! for the spectral mode the periodic time must be taken into
           ! account, which can be different for every section.
           !print *,'sections'
           t = timeUnsteadyRestart
           
           if(equationMode == timeSpectral) then
              do nnn=1,nSections
                 t(nnn) = t(nnn) + (mm-1)*sections(nnn)%timePeriod &
                      /         real(nTimeIntervalsSpectral,realType)
              enddo
           endif
           !print *,'grid velocities',t,mm
           call gridVelocitiesFineLevel(.false., t, mm)
           !print *,'grid velocities coarse'
           call gridVelocitiesCoarseLevels(mm)
           !print *,'normal velocities'
           call normalVelocitiesAllLevels(mm)
           !print *,'slip velocities'
           call slipVelocitiesFineLevel(.false., t, mm)
           call slipVelocitiesCoarseLevels(mm)
           
        enddo

      end subroutine updateFlow
