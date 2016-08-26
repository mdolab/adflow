!
!      ******************************************************************
!      *                                                                *
!      * File:          defaultIsoOut.f90                               *
!      * Author:        Gaetan Kenway                                   *
!      * Starting date: 07-21-2013                                      *
!      * Last modified: 07-21-2013                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine defaultIsoOut
!
!      ******************************************************************
!      *                                                                *
!      * defaultIsoOut sets the default set of additional            *
!      * variables to be written to the solution file; the primitive    *
!      * variables are always written. This additional set depends on   *
!      * the governing equations to be solved.                          *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use extraOutput
       use inputPhysics, only : equations
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
       isoWriteRho      = .false.
       isoWriteVx       = .false.
       isoWriteVy       = .false.
       isoWriteVz       = .false.
       isoWriteP        = .false.
       isoWriteMx       = .false.
       isoWriteMy       = .false.
       isoWriteMz       = .false.
       isoWriteRhoe     = .false.
       isoWriteTemp     = .false.
       isoWriteCp       = .false.
       isoWriteMach     = .false.
       isoWriteMachTurb = .false.
       isoWriteDist     = .false.
       isoWriteVort     = .false.
       isoWriteVortx    = .false.
       isoWriteVorty    = .false.
       isoWriteVortz    = .false.
       isoWritePtotloss = .false.
       isoWriteResRho   = .false.
       isoWriteResMom   = .false.
       isoWriteResRhoe  = .false.
       isoWriteShock    = .false.
       isoWriteFilteredShock = .false.
       ! Set the values which depend on the equations to be solved.

       select case (equations)
         case (EulerEquations)
           isoWriteEddyVis      = .false.
           isoWriteRatioEddyVis = .false.
           isoWriteResTurb      = .false.

         case (NSEquations)
           isoWriteEddyVis      = .false.
           isoWriteRatioEddyVis = .false.
           isoWriteResTurb      = .false.

         case (RANSEquations)
           isoWriteEddyVis      = .false.
           isoWriteRatioEddyVis = .false.
           isoWriteResTurb      = .false.
       end select

       end subroutine defaultIsoOut
