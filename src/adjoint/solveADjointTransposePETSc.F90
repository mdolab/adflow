!
!     ******************************************************************
!     *                                                                *
!     * File:          solveADjointTransposePETSc.F90                  *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 05-07-2010                                      *
!     * Last modified: 05-14-2010                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine solveADjointTransposePETSc
!
!     ******************************************************************
!     *                                                                *
!     * Solve the linear discrete ADjoint system of equations          *
!     *                                                                *
!     *     [dR/dW]T . psi = {dJdW}                                    *
!     *                                                                *
!     * using preconditioned GMRES provided by PETSc. It also stores   *
!     * the residual convergence history and the number of iterations  *
!     * to convergence. The convergence is verified by computing the   *
!     * norm of the residual vector r =  {dJdW} - [dR/dW]T . psi       *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use inputADjoint
      implicit none
!
!     Local variables.
!
      ! norm - Norm of error

      real(kind=realType)   :: norm

      real(kind=realType), dimension(2) :: time
      real(kind=realType)               :: timeAdjLocal, timeAdj,l2abs,curRes

!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
#ifndef USE_NO_PETSC

      ! Send some feedback to screen.

      if( PETScRank==0 .and. printTiming)  &
           write(*,10) "Solving ADjoint Transpose with PETSc..."

      ! Get the initial time.

      call cpu_time(time(1))

      if (restartADjoint) then
         !The user wants to restart the adjoint from the last point. Set
         !initial guess non-zero to true instead of zeroing the vector
         call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,PETScIerr)

      else
         call VecSet(psi,PETScZero,PETScIerr)
         call EChk(PETScIerr,__file__,__line__)
      end if
!
!     ******************************************************************
!     *                                                                *
!     * Solve the linear system of equations dRdWT . psi = dIdW using  *
!     * preconditioned GMRES.                                          *
!     *                                                                *
!     ******************************************************************

      adjResHist = 0.0_realType

      call KSPSetResidualHistory(ksp, adjResHist, adjMaxIter, &
                                 PETSC_FALSE, PETScIerr)
      call EChk(PETScIerr,__file__,__line__)

      ! If the user is doing a MDO problem there may be an
      ! agumentation to the RHS. This is set in agumentRHS.F90 and
      ! results in a non-zero vector in adjointRHS. For an aero-only
      ! problem this vector should be zero at this point. We compute:
      ! adjointRHS = -adjointRHS + dJdw
      
      call VecAYPX(adjointRHS,-1.0,dJdw,PETScIerr)
      call EChk(PETScIerr,__file__,__line__)

      ! Get Current Residual
      call MatMult(dRdWT,psi,adjointRes,PETScIerr)
      call EChk(PETScIerr,__file__,__line__)
      call VecAXPY(adjointRes,PETScNegOne,adjointRHS,PETScIerr)
      call EChk(PETScIerr,__file__,__line__)
      call VecNorm(adjointRHS,NORM_2,curRes,PETScIerr)  
      call EChk(PETScIerr,__file__,__line__)

      ! We are only going to overwrite adjRelTol and adjAbsTol

      L2abs = curRes * adjreltolrel
      
      if (L2Abs < adjAbsTol) then
         L2abs = adjabstol
      end if

      ! Solve the adjoint system of equations [dR/dW]T psi = dJ/dW. 

      call KSPSetTolerances(ksp,adjRelTol,L2Abs,adjDivTol,adjMaxIter,PETScIerr)
      call EChk(PETScIerr,__file__,__line__)
      call KSPSolve(ksp,dJdW,psi,PETScIerr)
      call EChk(PETScIerr,__file__,__line__)

      ! Get new time and compute the elapsed time.
      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, SUMB_PETSC_COMM_WORLD, PETScIerr)
!
!     ******************************************************************
!     *                                                                *
!     * Check the solution.                                            *
!     *                                                                *
!     ******************************************************************
!
      call MatMult(dRdWT,psi,adjointRes,PETScIerr)
      call EChk(PETScIerr,__file__,__line__)
      
      call VecAXPY(adjointRes,PETScNegOne,adjointRHS,PETScIerr)
      call EChk(PETScIerr,__file__,__line__)
      
      call VecNorm(adjointRes,NORM_2,norm,PETScIerr)
      call EChk(PETScIerr,__file__,__line__)
      
      ! Finally we MUST zero adjointRHS
      call VecZeroEntries(adjointRHS,PETscIerr)
      call EChk(PETScIerr,__file__,__line__)

      call KSPGetIterationNumber(ksp,adjConvIts,PETScIerr)
      call EChk(PETScIerr,__file__,__line__)
      
      ! Use the root processor to display the output summary, such as
      ! the norm of error and the number of iterations

      if( PETScRank==0 .and. printTiming) then
        write(*,20) "Solving ADjoint Transpose with PETSc time (s) =", timeAdj
        write(*,30) "Norm of error =",norm,"Iterations =",adjConvIts
        write(*,*) "------------------------------------------------"
        if( adjConvIts.lt.0 ) then
          write(*,40) "PETSc solver diverged after", -adjConvIts, &
                     "iterations..."
        else
          write(*,40) "PETSc solver converged after", adjConvIts, &
                       "iterations."
        endif
        write(*,*) "------------------------------------------------"
      endif

      ! Flush the output buffer and synchronize the processors.

      call f77flush()

      ! Output formats.

   10 format(a)
   20 format(a,1x,f8.2)
   30 format(1x,a,1x,e10.4,4x,a,1x,i4)
   40 format(1x,a,1x,i5,1x,a)

#endif
      
    end subroutine solveADjointTransposePETSc
      
