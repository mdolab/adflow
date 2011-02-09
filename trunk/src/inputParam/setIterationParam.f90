!
!      ******************************************************************
!      *                                                                *
!      * File:          setIterationParam.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-30-2003                                      *
!      * Last modified: 07-11-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setIterationParam
!
!      ******************************************************************
!      *                                                                *
!      * setIterationParam sets the parameters in the module            *
!      * iteration, which control the iterative strategy to be used.    *
!      *                                                                *
!      ******************************************************************
!
       use inputPhysics
       use inputIteration
       use iteration
       use flowVarRefState
       implicit none
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the equations to be solved.

       select case (equations)
         case (EulerEquations, NSEquations)

           ! Euler or laminar navier-stokes equations. Set the parameters
           ! such that the turbulent stuff is not called.

           nMGVar = nw
           nt1MG  = nt1
           nt2MG  = nt2

           restrictEddyVis = .false.
           turbSegregated  = .false.
           turbCoupled     = .false.

         case (RANSEquations)

           ! Rans computation. Set the parameters, depending on the
           ! user input.

           select case (turbTreatment)

             case (segregated)

               ! Segregated treatment of the turbulent transport
               ! equations. Multigrid is only applied to the mean
               ! flow equations.

               nMGVar = nwf
               nt1MG  = nwf + 1
               nt2MG  = nwf

               restrictEddyVis = .false.
               if( eddyModel ) restrictEddyVis = .true.

               turbSegregated = .true.
               turbCoupled    = .false.

               ! Increment nSubIterTurb, because at the moment it
               ! contains the number of additional subiterations.
               ! Also protect against a bad input value of nSubIterTurb.

               nSubIterTurb = nSubIterTurb + 1
               nSubIterTurb = max(nSubIterTurb, 1_intType)

             case (coupled)

               ! Mean flow and turbulent transport equations are solved
               ! in a coupled manner. Multigrid must be applied to all
               ! variables.

               nMGVar = nw
               nt1MG  = nt1
               nt2MG  = nt2

               restrictEddyVis = .false.
               turbSegregated  = .false.
               turbCoupled     = .true.

               ! If additional subiterations must be used for the
               ! turbulence this happens with a frozen flow field.
               ! Set turbSegregated to .true.

               if(nSubIterTurb > 0) turbSegregated = .true.

           end select

           ! For the algebraic models, baldwin lomax, the parameters
           ! specified above must be altered such that the code works
           ! correctly.

           select case (turbModel)
             case (baldwinLomax)

               ! Set the number of multigrid variables to the number of
               ! flow equations (== total number of equations) and set
               ! the number of turbulent equations accordingly.

               nMGVar = nw
               nt1MG  = nt1
               nt2MG  = nt2

               ! Make sure that the eddy-viscosity is restricted to the
               ! coarser levels and that the eddy viscosity is not
               ! recomputed on these levels.

               restrictEddyVis = .true.
               turbSegregated   = .false.
               turbCoupled      = .false.

           end select

       end select

       end subroutine setIterationParam
