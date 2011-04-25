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
      use killsignals
      use inputADjoint
      use communication
      implicit none
!
!     Local variables.
!
      ! norm - Norm of error

      real(kind=realType)   :: norm,temp

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

      if( myid ==0 .and. printTiming)  &
           write(*,10) "Solving ADjoint Transpose with PETSc..."

      ! Allocate resHist of not already done so
      if (.not. allocated(adjResHist))then
         allocate(adjResHist(adjMaxIter))
      endif

      call cpu_time(time(1))

      !Get the initial time.
      if (restartADjoint) then
         !The user wants to restart the adjoint from the last point. Set
         !initial guess non-zero to true instead of zeroing the vector
         call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,PETScIerr)
      else
         call VecSet(psi,PETScZero,PETScIerr)
         call EChk(PETScIerr,__FILE__,__LINE__)
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
      call EChk(PETScIerr,__FILE__,__LINE__)

      ! If the user is doing a MDO problem there may be an
      ! agumentation to the RHS. This is set in agumentRHS.F90 and
      ! results in a non-zero vector in adjointRHS. For an aero-only
      ! problem this vector should be zero at this point. We compute:
      ! adjointRHS = -adjointRHS + dJdw
      
      call VecAYPX(adjointRHS,-1.0,dJdw,PETScIerr)
      call EChk(PETScIerr,__FILE__,__LINE__)


      ! Get Current Residual
      call MatMult(dRdWT,psi,adjointRes,PETScIerr)
      call EChk(PETScIerr,__FILE__,__LINE__)
      
      ! AdjointRes = AdjointRes - adjointRHS
      call VecAXPY(adjointRes,-1.0,adjointRHS,PETScIerr)
      call EChk(PETScIerr,__FILE__,__LINE__)
      
      ! Norm of adjoint Residual
      call VecNorm(adjointRes,NORM_2,curRes,PETScIerr)
      call EChk(PETScIerr,__FILE__,__LINE__)

      ! L2Abs is used to stipulate an exit criteria for adjreltolrel
      L2abs = curRes * adjreltolrel
      
      ! If L2Abs is less that what we actually want as the absolute
      ! tolerance, clip it
      if (L2Abs < adjAbsTol) then
         L2abs = adjabstol
      end if

      ! Set the tolerances
      call KSPSetTolerances(ksp,adjRelTol,L2Abs,adjDivTol,adjMaxIter,PETScIerr)
      call EChk(PETScIerr,__FILE__,__LINE__)

    
      ! Solve the adjoint system of equations [dR/dW]T psi = adjointRHS
      call KSPSolve(ksp,adjointRHS,psi,PETScIerr)
      call EChk(PETScIerr,__FILE__,__LINE__)

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
      call EChk(PETScIerr,__FILE__,__LINE__)
      
      call VecAXPY(adjointRes,-1.0,adjointRHS,PETScIerr)
      call EChk(PETScIerr,__FILE__,__LINE__)
      
      call VecNorm(adjointRes,NORM_2,norm,PETScIerr)
      call EChk(PETScIerr,__FILE__,__LINE__)
      
      ! Finally we MUST zero adjointRHS
      call VecZeroEntries(adjointRHS,PETscIerr)
      call EChk(PETScIerr,__FILE__,__LINE__)

      call KSPGetIterationNumber(ksp,adjConvIts,PETScIerr)
      call EChk(PETScIerr,__FILE__,__LINE__)
      
      ! Use the root processor to display the output summary, such as
      ! the norm of error and the number of iterations

      if( myid ==0 .and. printTiming) then
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

      ! Get the petsc converged reason and set the fail flag

      call KSPGetConvergedReason(ksp, adjointConvergedReason,PETScIerr)
      call EChk(PETScIerr,__FILE__,__LINE__)

      if (adjointConvergedReason ==  KSP_CONVERGED_RTOL .or. &
          adjointConvergedReason ==  KSP_CONVERGED_ATOL .or. &
          adjointConvergedReason ==  KSP_CONVERGED_HAPPY_BREAKDOWN) then
         adjointFailed = .False.
      else
         adjointFailed = .True.
      end if


   
      ! Output formats.

   10 format(a)
   20 format(a,1x,f8.2)
   30 format(1x,a,1x,e10.4,4x,a,1x,i4)
   40 format(1x,a,1x,i5,1x,a)

#endif
    end subroutine solveADjointTransposePETSc
      
