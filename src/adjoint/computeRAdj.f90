!
!      ******************************************************************
!      *                                                                *
!      * File:          computeRAdj.f90                                 *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 02-01-2008                                      *
!      * Last modified: 02-01-2001                                      *
!      *                                                                *
!      ******************************************************************
!

subroutine computeRAdjoint()

!      Set Use Modules




!      Set Passed in Variables




!      Set Local Variables




! *************************************************************************
!      Begin Execution
! *************************************************************************


!      Mimic the Residual calculation in the main code

       !Compute the Pressure in the stencil based on the current 
       !States

       ! replace with Compute Pressure Adjoint!
       call computePressureAdj(wAdj, pAdj)
       

       ! Apply all boundary conditions to stencil.
       ! In case of a full mg mode, and a segegated turbulent solver,
       ! first call the turbulent boundary conditions, such that the
       ! turbulent kinetic energy is properly initialized in the halo's.

!###! Ignore Viscous for now
!###!       if(turbSegregated .and. (.not. corrections)) &
!###!         call applyAllTurbBCAdj(secondHalo)

       ! Apply all boundary conditions of the mean flow.

       call applyAllBCAdj(secondHalo)

       ! If case this routine is called in full mg mode call the mean
       ! flow boundary conditions again such that the normal momentum
       ! boundary condition is treated correctly.

       if(.not. corrections) call applyAllBCAdj(secondHalo)

       !Leave out State exchanges for now. If there are discrepancies 
       !Later, this may be a source...
!!$       ! Exchange the solution. Either whalo1 or whalo2
!!$       ! must be called.
!!$
!!$       if( secondHalo ) then
!!$         call whalo2(currentLevel, 1_intType, nVarInt, .true., &
!!$                     .true., .true.)
!!$       else
!!$         call whalo1(currentLevel, 1_intType, nVarInt, .true., &
!!$                     .true., .true.)
!!$       endif

       ! For full multigrid mode the bleeds must be determined, the
       ! boundary conditions must be applied one more time and the
       ! solution must be exchanged again.

       if(.not. corrections) then
         call BCDataMassBleedOutflowAdj(.true., .true.)
         call applyAllBCAdj(secondHalo)

       !Leave out State exchanges for now. If there are discrepancies 
       !Later, this may be a source...
!!$
!!$         if( secondHalo ) then
!!$           call whalo2(currentLevel, 1_intType, nVarInt, .true., &
!!$                       .true., .true.)
!!$         else
!!$           call whalo1(currentLevel, 1_intType, nVarInt, .true., &
!!$                       .true., .true.)
!!$         endif
       endif



       ! Reset the values of rkStage and currentLevel, such that
       ! they correspond to a new iteration.

       rkStage = 0
       currentLevel = groundLevel

       ! Compute the latest values of the skin friction velocity.
       ! The currently stored values are of the previous iteration.

       call computeUtauAdj

       ! Apply an iteration to the turbulent transport equations in
       ! case these must be solved segregatedly.

       if( turbSegregated ) call turbSolveSegregatedAdj

       ! Compute the time step.

       call timeStepAdj(.false.)

       ! Compute the residual of the new solution on the ground level.

       if( turbCoupled ) then
         call initresAdj(nt1MG, nMGVar)
         call turbResidualAdj
       endif

       call initresAdj(1_intType, nwf)
       call residualAdj



     end subroutine computeRAdjoint
