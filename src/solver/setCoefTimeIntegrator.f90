!
!      ******************************************************************
!      *                                                                *
!      * File:          setCoefTimeIntegrator.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-05-2004                                      *
!      * Last modified: 03-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setCoefTimeIntegrator
!
!      ******************************************************************
!      *                                                                *
!      * setCoefTimeIntegrator determines the coefficients of the       *
!      * time integration scheme in unsteady mode. Normally these are   *
!      * equal to the coefficients corresponding to the specified       *
!      * accuracy. However during the initial phase there are not       *
!      * enough states in the past and the accuracy is reduced.         *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use inputUnsteady
       use iteration
       use monitor
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: nn, nLevelsSet
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine which time integrator must be used.

       select case (timeAccuracy)
         case (firstOrder)

           ! 1st order. No need to check the number of available
           ! states in the past. Set the two coefficients and
           ! nLevelsSet to 2.

           coefTime(0) =  1.0_realType
           coefTime(1) = -1.0_realType

           nLevelsSet = 2

           !=============================================================

         case (secondOrder)

           ! Second order time integrator. Determine the amount of
           ! available states and set the coefficients accordingly.

           select case (nOldSolAvail)

             case (1_intType)
               coefTime(0) =  1.0_realType
               coefTime(1) = -1.0_realType
               nLevelsSet  = 2

             case default   ! 2 or bigger.
               coefTime(0) =  1.5_realType
               coefTime(1) = -2.0_realType
               coefTime(2) =  0.5_realType
               nLevelsSet  = 3

           end select

           !=============================================================

         case (thirdOrder)

           ! Third order time integrator.  Determine the amount of
           ! available states and set the coefficients accordingly.

           select case (nOldSolAvail)

             case (1_intType)
               coefTime(0) =  1.0_realType
               coefTime(1) = -1.0_realType
               nLevelsSet  = 2

             case (2_intType)
               coefTime(0) =  1.5_realType
               coefTime(1) = -2.0_realType
               coefTime(2) =  0.5_realType
               nLevelsSet  = 3

             case default   ! 3 or bigger.
               coefTime(0) = 11.0_realType/6.0_realType
               coefTime(1) = -3.0_realType
               coefTime(2) =  1.5_realType
               coefTime(3) = -1.0_realType/3.0_realType
               nLevelsSet  = 4

           end select

       end select

       ! Set the rest of the coefficients to 0 if not enough states
       ! in the past are available.

       do nn=nLevelsSet,nOldLevels
         coefTime(nn) = zero
       enddo

       end subroutine setCoefTimeIntegrator
