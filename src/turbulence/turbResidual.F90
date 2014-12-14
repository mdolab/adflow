!
!      ******************************************************************
!      *                                                                *
!      * File:          turbResidual.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-01-2003                                      *
!      * Last modified: 04-11-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine turbResidual
!
!      ******************************************************************
!      *                                                                *
!      * turbResidual computes the residual of the residual of the      *
!      * turbulent transport equations on the current multigrid level.  *
!      *                                                                *
!      ******************************************************************
!
       use inputDiscretization
       use inputPhysics
       use iteration
       use turbMod
       implicit none
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine whether or not a second order discretization for the
       ! advective terms must be used.

       secondOrd = .false.
       if(currentLevel == 1_intType .and. &
          orderTurb    == secondOrder) secondOrd = .true.

       ! Compute the quantities for certain turbulence models that
       ! need to be communicated between blocks.

       select case (turbModel)

         case (menterSST)
           call f1SST

       end select

       ! Determine the turbulence model we have to solve for.

       select case (turbModel)

         case (baldwinLomax)
           call bl(.true.)

         case (spalartAllmaras)
           call sa(.true.)

         case (komegaWilcox, komegaModified)
           call kw(.true.)

         case (menterSST)
           call SST(.true.)

         case (ktau)
           call kt(.true.)

         case (v2f)
           call vf(.true.)

         case default
           call terminate("turbResidual", &
                          "Turbulence model not implemented yet")

       end select

       end subroutine turbResidual
