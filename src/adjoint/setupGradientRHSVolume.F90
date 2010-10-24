
!
!     ******************************************************************
!     *                                                                *
!     * File:          setupGradientRHSVolume.F90                      *
!     * Author:         C.A.(Sandy) Mader                              *
!     * Starting date: 11-27-2009                                      *
!     * Last modified: 02-02-2010                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine setupGradientRHSVolume(costFunction)
        !subroutine setupGradientRHS(level,sps,costFunction)
!
!     ******************************************************************
!     *                                                                *
!     * Compute the partial derivative for the discrete ADjoint problem*
!     * in question. Notice that this right hand side is problem /     *
!     * cost function J dependent.                                     *
!     *                                                                *
!     * The ordering of the unknowns in the ADjoint vector used here   *
!     * is based on the global node numbering.                         *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use ADjointVars
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) ::  costFunction
!
!     Local variables.
!
      integer(kind=intType)::sps=1,level=1
      real(kind=realType), dimension(2) :: time
      real(kind=realType)               :: timeAdjLocal, timeAdj
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
#ifndef USE_NO_PETSC

      ! Send some feedback to screen.

      if( PETScRank==0 ) &
        write(*,10) "Assembling Gradient RHS vector..."

      ! Get the initial time.
      level = 1
      call cpu_time(time(1))
!
!     ******************************************************************
!     *                                                                *
!     * Select case over different choices of cost functions.          *
!     * The following routine calls fill in the PETSc vector dJdW by   *
!     * making calls to the PETSc routine VecSetValuesBlocked.         *
!     *                                                                *
!     ******************************************************************
!
      select case (costFunction)

        case (costFuncLiftCoef, &
              costFuncDragCoef, &
              costFuncForceXCoef,&
              costFuncForceYCoef,&
              costFuncForceZCoef,&
              costFuncMomXCoef, &
              costFuncMomYCoef, &
              costFuncMomZCoef)

           !For now, hard code this to the first time instance for a time
           !spectral case. In the long run maybe we specify this from the
           !input file...
           sps = 1
           
           
           call setupGradientRHSSpatial(level,costFunction,sps) 
   

        case(costFuncCmzAlpha, &
             costFuncCm0,&
             costFuncClAlpha,&
             costFuncCl0,&
             costFuncCdAlpha,&
             costFuncCd0,&
             costFuncCmzAlphaDot,&
             costFuncCmzq)

           call setupGradientRHSStability(level,costFunction)
        case default
          write(*,*) "Invalid cost function ", costFunction
          stop

      end select
!

      call f77flush()
      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

      ! Output formats.

   10 format(a)
   20 format(a,1x,f8.2)

#endif

    end subroutine setupGradientRHSVolume
