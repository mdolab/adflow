!
!      ******************************************************************
!      *                                                                *
!      * File:          describeScheme.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-23-2003                                      *
!      * Last modified: 08-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine describeScheme(string)
!
!      ******************************************************************
!      *                                                                *
!      * describeScheme gives a short description about the scheme      *
!      * used to obtain the solution. The description is stored in the  *
!      * character array string.                                        *
!      *                                                                *
!      ******************************************************************
!
       use inputDiscretization
       use inputPhysics
       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       character(len=*), intent(out) :: string
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Write the basic scheme info.

       select case(spaceDiscr)
         case (dissScalar)
           write(string,100) "Scalar dissipation scheme", vis2, vis4
         case (dissMatrix)
           write(string,100) "Matrix dissipation scheme", vis2, vis4
         case (dissCusp)
           write(string,100) "CUSP dissipation scheme", vis2, vis4
         case (upwind)
           select case (limiter)
             case (firstOrder)
               write(string,110) "First order upwind scheme."
             case (noLimiter)
               write(string,111) kappaCoef
             case (vanAlbeda)
               write(string,112) kappaCoef
             case (minmod)
               write(string,113) kappaCoef
           end select

           select case (riemann)
             case (Roe)
               write(string,130) trim(string), &
                                 "Roe's approximate Riemann Solver."
             case (vanLeer)
               write(string,130) trim(string), &
                                 "Van Leer flux vector splitting."
             case (ausmdv)
               write(string,130) trim(string), &
                                 "ausmdv flux vector splitting."
           end select
 
       end select

 100   format(a,1x,", k2 = ", e12.5, ", k4 = ", e12.5,".")
 110   format(a)
 111   format("Second order upwind scheme using linear reconstruction, &
              &i.e. no limiter, kappa =", 1x,f7.3,".")
 112   format("Second order upwind scheme with Van Albeda limiter, &
              &kappa =", 1x,f7.3,".")
 113   format("Second order upwind scheme with Minmod limiter, &
              &kappa =", 1x,f7.3,".")
 130   format(a,1x,a)

       ! In case of the scalar dissipation scheme, write whether or not
       ! directional scaling has been applied.

       if(spaceDiscr == dissScalar) then
         if( dirScaling ) then
           write(string,200) trim(string), adis
         else
           write(string,210) trim(string)
         endif
       endif

 200   format(A,1X,"Directional scaling of dissipation with exponent", &
              1x,e12.5, ".")
 210   format(A,1X,"No directional scaling of dissipation.")

       ! For the Euler equations, write the inviscid wall boundary
       ! condition treatment.

       if(equations == EulerEquations) then
         select case (wallBcTreatment)
           case (constantPressure)
             write(string,300) trim(string), &
                               "Zero normal pressure gradIent"
           case (linExtrapolPressure)
             write(string,300) trim(string), &
                               "Linear extrapolation of normal &
                               &pressure gradIent"
           case (quadExtrapolPressure)
             write(string,300) trim(string), &
                               "Quadratic extrapolation of normal &
                               &pressure gradIent"
           case (normalMomentum)
             write(string,300) trim(string), &
                               "Normal momentum equation used to &
                               &determine pressure gradIent"
         end select
       endif
 300   format(A,1X,A,1X,"for inviscid wall boundary conditions.")

       ! If preconditioning is used, write the preconditioner.

       select case(precond)
         case (Turkel)
           write(string,400) trim(string), &
                             "Turkel preconditioner for inviscid fluxes."
         case (ChoiMerkle)
           write(string,400) trim(string), &
                        "Choi Merkle preconditioner for inviscid fluxes."
       end select
 400   format(a,1x,a)

       ! For a viscous computation write that a central discretization
       ! is used for the viscous fluxes.

       if( viscous ) then
         write(string,400) trim(string), &
                           "Central discretization for viscous fluxes."
       endif

       end subroutine describeScheme
