!
!      ******************************************************************
!      *                                                                *
!      * File:          setEquationParameters.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-21-2002                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setEquationParameters
!
!      ******************************************************************
!      *                                                                *
!      * setEquationParameters sets the number of variables in the      *
!      * governing equations, the number of turbulent variables, etc.   *
!      *                                                                *
!      ******************************************************************
!
       use flowVarRefState
       use inputPhysics
       use paramTurb
       implicit none
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the number of flow variables to 5, nt1 to 6. This is valid
       ! for all governing equations. Furthermore initialize viscous,
       ! kPresent and eddyModel to .False., which indicates an inviscid
       ! computation. For ns and rans this will be corrected.

       nwf = 5
       nt1 = 6

       viscous    = .false.
       kPresent  = .false.
       eddyModel = .false.

       ! Determine the set of governing equations to solve for and set
       ! the parameters accordingly.

       select case (equations)
         case (EulerEquations)
           nw  = 5
           nt2 = 5

         !===============================================================

         case (NSEquations)
           nw  = 5
           nt2 = 5

           viscous = .true.

         !===============================================================

         case (RANSEquations)

           viscous = .true.

           select case(turbModel)
             case (baldwinLomax)
               nw  = 5
               nt2 = 5

               eddyModel = .true.

             !===========================================================

             case (spalartAllmaras)
               nw  = 6
               nt2 = 6

               eddyModel = .true.
               if( wallFunctions ) call initCurveFitDataSa

             !===========================================================

             case (spalartAllmarasEdwards)
               nw  = 6
               nt2 = 6

               eddyModel = .true.
               if( wallFunctions ) call initCurveFitDataSae

             !===========================================================

             case (komegaWilcox)
               nw  = 7
               nt2 = 7

               kPresent  = .true.
               eddyModel = .true.
               if( wallFunctions ) call initCurveFitDataKw

             !===========================================================

             case (komegaModified)
               nw  = 7
               nt2 = 7

               kPresent  = .true.
               eddyModel = .true.
               if( wallFunctions ) call initCurveFitDataKwMod

             !===========================================================

             case (menterSST)
               nw  = 7
               nt2 = 7

               kPresent  = .true.
               eddyModel = .true.
               if( wallFunctions ) call initCurveFitDataSST

             !===========================================================

             case (ktau)
               nw  = 7
               nt2 = 7

               kPresent  = .true.
               eddyModel = .true.
               if( wallFunctions ) call initCurveFitDataKtau

             !===========================================================

             case (v2f)
               nw  = 9
               nt2 = 9

               rvfLimitK = 1.e-25_realType
               rvfLimitE = 1.e-25_realType

               if(rvfN == 6) then
                 rvfCmu   = rvfN6Cmu
                 rvfCl    = rvfN6Cl
               else
                 rvfCmu   = rvfN1Cmu
                 rvfCl    = rvfN1Cl
               endif

               kPresent  = .true.
               eddyModel = .true.
               if( wallFunctions ) call initCurveFitDataVf

           end select

       end select

       ! Determine the number of turbulent variables.

       nwt = nw - nwf

       end subroutine setEquationParameters
