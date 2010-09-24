!
!      ******************************************************************
!      *                                                                *
!      * File:          initDepvarAndHalos.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-03-2004                                      *
!      * Last modified: 10-31-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initDepvarAndHalos(halosRead)
!
!      ******************************************************************
!      *                                                                *
!      * InitDepvarAndHalos computes the dependent flow variables,      *
!      * like viscosities, and initializes the halo cells by applying   *
!      * the boundary conditions and exchanging the internal halo's.    *
!      * This is all done on the start level grid.                      *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use flowVarRefState
       use inputIO
       use inputIteration
       use inputPhysics
       use inputTimeSpectral
       use iteration
       use monitor
       use section
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: halosRead
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm
       real(kind=realType)   :: relaxBleedsOr

       real(kind=realType), dimension(nSections) :: t

       logical :: initBleeds
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the logical whether or not to initialize the prescribed
       ! data for the bleed regions. If the halos were read the bleeds
       ! have been initialized already and nothing needs to be done.

       if( halosRead ) then
         initBleeds = .false.
       else
         initBleeds = .true.
       endif

       ! Compute the face velocities and for viscous walls the slip
       ! velocities. This is done for all the mesh levels.

       currentLevel = 1
       groundLevel  = 1

       do mm=1,nTimeIntervalsSpectral

         ! Compute the time, which corresponds to this spectral solution.
         ! For steady and unsteady mode this is simply the restart time;
         ! for the spectral mode the periodic time must be taken into
         ! account, which can be different for every section.

         t = timeUnsteadyRestart

         if(equationMode == timeSpectral) then
           do nn=1,nSections
             t(nn) = t(nn) + (mm-1)*sections(nn)%timePeriod &
                   /         real(nTimeIntervalsSpectral,realType)
           enddo
         endif

         call gridVelocitiesFineLevel(.false., t, mm)
         call gridVelocitiesCoarseLevels(mm)
         call normalVelocitiesAllLevels(mm)

         call slipVelocitiesFineLevel(.false., t, mm)
         call slipVelocitiesCoarseLevels(mm)

       enddo

       ! Loop over the number of spectral solutions and blocks
       ! to compute the laminar viscosity.

       do mm=1,nTimeIntervalsSpectral
         do nn=1,nDom
           call setPointers(nn,mgStartlevel,mm)
           call computeLamViscosity
         enddo
       enddo

       ! Exchange the solution on the multigrid start level.
       ! It is possible that the halo values are needed for the boundary
       ! conditions. Viscosities are not exchanged.

       call whalo2(mgStartlevel, 1_intType, nw, .true., .true., &
                   .false.)

       ! Apply all flow boundary conditions to be sure that the halo's
       ! contain the correct values. These might be needed to compute
       ! the eddy-viscosity. Also the data for the outflow bleeds
       ! is determined. As BCDataMassBleedOutflow is called one more
       ! time in this routine, relaxBleeds is set to zero when this
       ! routine is called here, such that the pressure is not adapted
       ! twice.

       currentLevel = mgStartlevel
       groundLevel  = mgStartlevel

       relaxBleedsOr = relaxBleeds
       relaxBleeds   = zero
       call BCDataMassBleedOutflow(initBleeds, .false.)
       relaxBleeds = relaxBleedsOr

       call applyAllBC(.true.)

       ! Loop over the number of spectral solutions and blocks
       ! to compute the eddy viscosities.

       do mm=1,nTimeIntervalsSpectral
         do nn=1,nDom

           ! Set the pointers for this block.

           call setPointers(nn,mgStartlevel,mm)

           ! Compute the eddy viscosity for rans computations using
           ! an eddy viscosity model.

           call computeEddyViscosity

           ! In case of a rans computation and no restart, initialize
           ! the turbulent variables a bit better for some turbulence
           ! models.

           ! if(equations == RANSEquations .and. (.not. restart)) then
           !
           !   select case (turbModel)
           !
           !     case (komegaWilcox, komegaModified, menterSST)
           !       call initKOmega(0_intType)
           !       call computeEddyViscosity
           !
           !   end select
           !
           ! endif

         enddo
       enddo

       ! Exchange the laminar and eddy viscosities.

       call whalo2(mgStartlevel, 1_intType, 0_intType, .false., &
                   .false., .true.)

       ! Apply the turbulent boundary conditions for a segregated
       ! solver twice and redo the mean flow boundary conditions.
       ! Just to be sure that everything is initialized properly for
       ! all situations possible.

       if( turbSegregated ) call applyAllTurbBC(.true.)
       call applyAllBC(.true.)
       call BCDataMassBleedOutflow(initBleeds, .false.)
       call applyAllBC(.true.)
       if( turbSegregated ) call applyAllTurbBC(.true.)

       ! Exchange the solution for the second time to be sure that all
       ! halo's are initialized correctly. As this is the initialization
       ! phase, this is not critical.

       call whalo2(mgStartlevel, 1_intType, nw, .true., .true., &
                   .true.)

       end subroutine initDepvarAndHalos
