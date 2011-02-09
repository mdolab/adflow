!
!      ******************************************************************
!      *                                                                *
!      * File:          defaultSurfaceOut.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-26-2003                                      *
!      * Last modified: 03-23-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine defaultSurfaceOut
!
!      ******************************************************************
!      *                                                                *
!      * defaultSurfaceOut sets the default set of surface variables    *
!      * to be written to the solution file. This set depends on the    *
!      * governing equations to be solved.                              *
!      *                                                                *
!      ******************************************************************
!
       use inputPhysics
       use extraOutput
       implicit none
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! First set the variables, which are independent from the
       ! governing equations to be solved.

       surfWriteRho  = .true.
       surfWriteP    = .false.
       surfWriteTemp = .false.
       surfWriteVx   = .true.
       surfWriteVy   = .true.
       surfWriteVz   = .true.
       surfWriteCp   = .true.
       surfWriteMach = .true.

       ! Set the values which depend on the equations to be solved.

       select case (equations)
         case (EulerEquations)
           surfWritePtotloss = .true.
           surfWriteCf       = .false.
           surfWriteCh       = .false.
           surfWriteYplus    = .false.
           surfWriteCfx      = .false.
           surfWriteCfy      = .false.
           surfWriteCfz      = .false.

         case (NSEquations)
           surfWritePtotloss = .false.
           surfWriteCf       = .true.
           surfWriteCh       = .false.
           surfWriteYplus    = .false.
           surfWriteCfx      = .true.
           surfWriteCfy      = .true.
           surfWriteCfz      = .true.

         case (RANSEquations)
           surfWritePtotloss = .false.
           surfWriteCf       = .true.
           surfWriteCh       = .false.
           surfWriteYplus    = .true.
           surfWriteCfx      = .true.
           surfWriteCfy      = .true.
           surfWriteCfz      = .true.
       end select

       end subroutine defaultSurfaceOut
