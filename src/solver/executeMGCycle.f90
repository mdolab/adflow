!
!      ******************************************************************
!      *                                                                *
!      * File:          executeMGCycle.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-15-2003                                      *
!      * Last modified: 08-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine executeMGCycle
!
!      ******************************************************************
!      *                                                                *
!      * executeMGCycle performs a multigrid cycle defined by           *
!      * cycling, see the module iteration.                             *
!      *                                                                *
!      ******************************************************************
!
       use flowVarRefState
       use iteration
       use inputIteration
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize currentLevel to groundLevel, the ground level
       ! of this multigrid cycle.

       currentLevel = groundLevel

       ! Loop over the number of steps in cycling.

       mgLoop: do nn=1,nstepsCycling

         ! Determine what must be done.

         select case (cycling(nn))

           case (-1_intType)

             ! Set the new currentLevel and and prolongate the
             ! the corrections to this grid level.

             currentLevel = currentLevel -1
             call transferToFineGrid(.true.)

           case ( 0_intType)

             ! Perform a smoothing iteration on this grid level.
             ! First determine the current situation. If this is the
             ! first entry in cycling the residual already contains the
             ! correct values and therefore it needs not to be computed.

             if(nn > 1) then

               ! Compute the residual if the previous action was not a
               ! restriction. In that case the residual already contains
               ! the correct values.

               if(cycling(nn-1) /= 1_intType) then

                 ! Initialize and compute the residual.

                 rkStage = 0
                 call timeStep(.false.)

                 if( turbCoupled ) then
                   call initres(nt1MG, nMGVar)
                   call turbResidual
                 endif

                 call initres(1_intType, nwf)
                 call residual

               endif
             endif

             ! Perform a smoothing step. Determine the smoother to
             ! be used and call the appropriate routine.

             select case (smoother)
               case (RungeKutta)
                 call RungeKuttaSmoother
               case (nlLusgs)
                 call terminate("executeMGCycle", &
                                "nlLusgs smoother not implemented yet")
               case (nlLusgsLine)
                 call terminate("executeMGCycle", &
                                "nlLusgsLine smoother not implemented &
                                &yet")
             end select

           case ( 1_intType)

             ! Restrict the solution and residual to the next coarser
             ! grid level. Inside transferToCoarseGrid currentLevel
             ! is updated.

             call transferToCoarseGrid

         end select

       enddo mgLoop

       ! Reset the values of rkStage and currentLevel, such that
       ! they correspond to a new iteration.

       rkStage = 0
       currentLevel = groundLevel

       ! Compute the latest values of the skin friction velocity.
       ! The currently stored values are of the previous iteration.

       call computeUtau

       ! Apply an iteration to the turbulent transport equations in
       ! case these must be solved segregatedly.

       if( turbSegregated ) call turbSolveSegregated

       ! Compute the time step.

       call timeStep(.false.)

       ! Compute the residual of the new solution on the ground level.

       if( turbCoupled ) then
         call initres(nt1MG, nMGVar)
         call turbResidual
       endif

       call initres(1_intType, nwf)
       call residual

       end subroutine executeMGCycle
