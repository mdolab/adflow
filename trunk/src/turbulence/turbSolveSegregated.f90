!
!      ******************************************************************
!      *                                                                *
!      * File:          turbSolveSegregated.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-01-2003                                      *
!      * Last modified: 07-21-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine turbSolveSegregated
!
!      ******************************************************************
!      *                                                                *
!      * turbSolveSegregated solves the turbulent transport equations   *
!      * segregatedly, i.e. the mean flow variables are kept constant   *
!      * and the turbulent variables are updated.                       *
!      *                                                                *
!      ******************************************************************
!
       use flowVarRefState
       use inputDiscretization
       use inputIteration
       use inputPhysics
       use iteration
       use turbMod
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: iter
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
       if(groundLevel == 1_intType .and. &
          orderTurb   == secondOrder) secondOrd = .true.

       ! Loop over the number of iterations for the turbulence.

       do iter=1,nSubIterTurb

         ! Compute the quantities for certain turbulence models that
         ! need to be communicated between blocks.

         select case (turbModel)

           case (menterSST)
             call f1SST

         end select

         ! Determine the turbulence model we have to solve for.
         ! Note that the eddy viscosity is updated inside these routines.

         select case (turbModel)

           case (baldwinLomax)
             call bl(.false.)

           case (spalartAllmaras)
             call sa(.false.)

           case (komegaWilcox, komegaModified)
             call kw(.false.)

           case (menterSST)
             call SST(.false.)

           case (ktau)
             call kt(.false.)

           case (v2f)
             call vf(.false.)

           case default
             call terminate("turbSolveSegregated", &
                            "Turbulence model not implemented yet")

         end select

         ! Exchange the halo data. As it is guaranteed that we are on the
         ! finest mesh, exchange both layers of halo's.

         call whalo2(groundLevel, nt1, nt2, .false., .false., .true.)

       enddo

       end subroutine turbSolveSegregated
