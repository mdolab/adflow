!
!      ******************************************************************
!      *                                                                *
!      * File:          setStageCoeffExplicitRK.f90                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-12-2006                                      *
!      * Last modified: 08-13-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setStageCoeffExplicitRK
!
!      ******************************************************************
!      *                                                                *
!      * setStageCoeffExplicitRK determines the coefficients of the     *
!      * stages for the explicit Runge Kutta time integration schemes   *
!      * for unsteady problems.                                         *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use inputUnsteady
       implicit none
!
!      Local variables.
!
       integer :: ierr
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of Runge Kutta stages as a function of
       ! the accuracy.

       select case (timeAccuracy)
         case (firstOrder)
           nRKStagesUnsteady = 1

         case (secondOrder)
           nRKStagesUnsteady = 2

         case (thirdOrder)
           nRKStagesUnsteady = 3

         case default
           call terminate("setStageCoeffExplicitRK", &
                          "No higher order stuff yet")
       end select

       ! Allocate and determine betaRKUnsteady and gammaRKUnsteady.

       allocate(betaRKUnsteady(nRKStagesUnsteady,nRKStagesUnsteady), &
                gammaRKUnsteady(nRKStagesUnsteady), stat=ierr)
       if(ierr /= 0)                               &
         call terminate("setStageCoeffExplicitRK", &
                        "Memory allocation failure for betaRKUnsteady &
                        &and gammaRKUnsteady.")

       betaRKUnsteady = zero

       select case (timeAccuracy)
         case (firstOrder)

           ! Just the forward Euler time integration scheme.

           betaRKUnsteady(1,1) = 1.0_realType
           gammaRKUnsteady(1)  = 0.0_realType

         !==============================================================

         case (secondOrder)

           ! The TVD Runge Kutta scheme which allows for the maximum
           ! CFL number (1.0).

           betaRKUnsteady(1,1) =  1.0_realType
           betaRKUnsteady(2,1) = -0.5_realType
           betaRKUnsteady(2,2) =  0.5_realType

           gammaRKUnsteady(1)  = 0.0_realType
           gammaRKUnsteady(2)  = 1.0_realType

         !==============================================================

         case (thirdOrder)

           ! Low storage (although not exploited in this implemetation)
           ! 3 stage scheme of Le and Moin.

           betaRKUnsteady(1,1) =   8.0_realType/15.0_realType
           betaRKUnsteady(2,1) = -17.0_realType/60.0_realType
           betaRKUnsteady(2,2) =   5.0_realType/12.0_realType
           betaRKUnsteady(3,2) =  -5.0_realType/12.0_realType
           betaRKUnsteady(3,3) =   3.0_realType/ 4.0_realType

           gammaRKUnsteady(1)  = 0.0_realType
           gammaRKUnsteady(2)  = 8.0_realType/15.0_realType
           gammaRKUnsteady(3)  = 2.0_realType/ 3.0_realType

           ! The TVD Runge Kutta scheme which allows for the maximum
           ! CFL number (1.0).

         ! betaRKUnsteady(1,1) =  1.0_realType
         ! betaRKUnsteady(2,1) = -3.0_realType/ 4.0_realType
         ! betaRKUnsteady(2,2) =  1.0_realType/ 4.0_realType
         ! betaRKUnsteady(3,1) = -1.0_realType/12.0_realType
         ! betaRKUnsteady(3,2) = -1.0_realType/12.0_realType
         ! betaRKUnsteady(3,3) =  2.0_realType/ 3.0_realType

         ! gammaRKUnsteady(1)  = 0.0_realType
         ! gammaRKUnsteady(2)  = 1.0_realType
         ! gammaRKUnsteady(3)  = 0.5_realType

         !==============================================================

         case default
           call terminate("setStageCoeffExplicitRK", &
                          "No higher order stuff yet")
       end select

       end subroutine setStageCoeffExplicitRK
