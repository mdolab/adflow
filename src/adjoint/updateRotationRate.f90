!
!     ******************************************************************
!     *                                                                *
!     * File:          updateRotationRate.f90                          *
!     * Authors:       C.A(Sandy) Mader                                *
!     * Starting date: 02-09-2010                                      *
!     * Last modified: 02-09-2010                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine updateRotationRate(xrot,yrot,zrot)
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
        use inputPhysics  !pointRef
        use inputMotion   !rotPoint
        use cgnsGrid      !cgnsdoms
        use monitor
        use iteration
        use blockpointers
        implicit none

        real(kind=realType),intent(in)::xrot,yrot,zrot

        integer(kind=intType) ::mm,nnn,nn
        real(kind=realType), dimension(nSections) :: t

        groundlevel = 1
        do mm=1,nTimeIntervalsSpectral
           do nn=1,nDom
              
              ! Set the pointers for this block.
              
              call setPointers(nn, groundLevel, mm)
              !lref is outside
              cgnsDoms(nbkglobal)%rotRate(1) = xrot
              cgnsDoms(nbkglobal)%rotRate(2) = yrot
              cgnsDoms(nbkglobal)%rotRate(3) = zrot
              
           enddo
        enddo
       
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

      end subroutine updateRotationRate
