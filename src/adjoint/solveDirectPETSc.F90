!
!     ******************************************************************
!     *                                                                *
!     * File:          solveDirectPETSc.F90                            *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 06-04-2008                                      *
!     * Last modified: 06-04-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine solveDirectPETSc
!
!     ******************************************************************
!     *                                                                *
!     * Solve the linear discrete ADjoint system of equations          *
!     *                                                                *
!     *     [dR/dW] . psi = {dRdx}                                     *
!     *                                                                *
!     * using preconditioned GMRES provided by PETSc. It also stores   *
!     * the residual convergence history and the number of iterations  *
!     * to convergence. The convergence is verified by computing the   *
!     * norm of the residual vector r =  {dRdx} - [dR/dW] . psi        *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc 
      use ADjointVars
      use flowVarRefState !nw
      use blockPointers   !il,jl,kl,globalcell
      use communication
      implicit none
!
!     Local variables.
!
      ! norm - Norm of error

      real(kind=realType)   :: norm,vecval

      real(kind=realType), dimension(2) :: time
      real(kind=realType)               :: timeAdjLocal, timeAdj

      character(len=2*maxStringLen) :: errorMessage

      integer       :: nDimW
      integer(kind=intType) :: i, j, k, iw,idxres,nn,ierr,localoffset
      integer(kind=intType) :: iFile = 8
      character*32:: testfile,file1
      write(*,*)'In solveDirect Petsc.....'
!________________________________________________________________________
!
!     Begin execution
!________________________________________________________________________
           
      open(unit=iFile,file='ADw.txt',status='replace',action='write')
      !open(unit=iFile+1,file='ADw2.txt',status='replace',action='write')

      call mpi_barrier(SUmb_comm_world, ierr)

      write(testfile,100) myid!12
!100   format ('testfile',i5)  
100   format (i5)  
      !testfile = 'testfile'
      print *,'file: ',testfile
      testfile=adjustl(testfile)
      print *,'file2: ',testfile
      write(file1,101) trim(testfile)!testfile
101   format("testfile",a,".out")
      print *,'file3: ',file1!trim(testfile)!testfile
      open(unit=iFile+1+myID,file=file1,status='replace',action='write')
      call mpi_barrier(SUmb_comm_world, ierr)


      nDimW = nw * nCellsLocal
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
        write(*,10) "Solving Direct with PETSc..."

      ! Get the initial time.

      call cpu_time(time(1))

      ! Reset the solution vector psi by assigning the value zero to all
      ! its components.

      ! VecSet - Sets all components of vector to a single scalar value.
      !
      ! Synopsis
      !
      ! #include "petscvec.h" 
      ! call VecSet(Vec x,PetscScalar alpha, PetscErrorCode ierr)
      !
      ! Collective on Vec
      !
      ! Input Parameters
      !   x     - the vector
      !   alpha	- the scalar
      !
      ! Output Parameter
      !   x     - the vector
      !
      ! Note
      ! For a vector of dimension n, VecSet() computes
      ! x[i] = alpha, for i=1,...,n,
      ! so that all vector entries then equal the identical scalar
      ! value, alpha. Use the more general routine VecSetValues() to
      ! set different vector entries.
      !
      ! You CANNOT call this after you have called VecSetValues() but
      ! before you call VecAssemblyBegin/End(). 
      !
      ! see .../petsc/docs/manualpages/Vec/VecSet.html
      ! or PETSc users manual, pp.36

      call VecSet(psi,PETScZero,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("solveADjointPETSc", "Error in VecSet")
!
!     ******************************************************************
!     *                                                                *
!     * Solve the linear system of equations dRdWT . psi = dIdW using  *
!     * preconditioned GMRES.                                          *
!     *                                                                *
!     ******************************************************************
!
      ! KSPSetOperators - Sets the matrix associated with the linear
      !                   system and a (possibly) different one
      !                   associated with the preconditioner.
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPSetOperators(KSP ksp,Mat Amat,Mat Pmat, &
      !                      MatStructure flag,PetscErrorCode ierr)
      !
      ! Collective on KSP and Mat
      !
      ! Input Parameters
      !   ksp  - the KSP context
      !   Amat - the matrix associated with the linear system
      !   Pmat - the matrix to be used in constructing the 
      !          preconditioner, usually the same as Amat.
      !   flag - flag indicating information about the preconditioner
      !          matrix structure during successive linear solves.
      !          This flag is ignored the first time a linear system is
      !          solved, and thus is irrelevant when solving just one
      !          linear system.
      !
      ! see .../petsc/docs/manualpages/KSP/KSPSetOperators.html
      ! or PETSc users manual, pp.63

      ! Here the matrix that defines the linear system
      ! also serves as the preconditioning matrix.

      !call KSPSetOperators(ksp,dRdW,dRdW, &
       !                    DIFFERENT_NONZERO_PATTERN,PETScIerr)

      !if( PETScIerr/=0 ) &
      !  call terminate("solveADjointPETSc", "Error in KSPSetOperators.")

      ! Save the residual history.

      ! KSPSetResidualHistory - Sets the array used to hold the residual
      !                         history. If set, this array will contain
      !                         the residual norms computed at each
      !                         iteration of the solver.
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPSetResidualHistory(KSP ksp,PetscReal a[],PetscInt na, &
      !                           PetscTruth reset, PetscErrorCode ierr)
      !
      ! Not Collective
      !
      ! Input Parameters
      !   ksp   - iterative context obtained from KSPCreate()
      !   a     - array to hold history
      !   na    - size of a
      !   reset - PETSC_TRUE indicates the history counter is reset to
      !           zero for each new linear solve
      !
      ! Notes: The array is NOT freed by PETSc so the user needs to keep
      !        track of it and destroy once the KSP object is destroyed.
      !
      ! If 'na' is PETSC_DECIDE or 'a' is PETSC_NULL, then a default
      !  array of length 1000 is allocated. 
      !
      ! see .../petsc/docs/manualpages/KSP/KSPSetResidualHistory.html

      adjResHist = 0.0_realType
      !print *,'asjhis',shape(adjResHist),adjmaxiter
      call KSPSetResidualHistory(ksp, adjResHist, adjMaxIter, &
                                 PETSC_FALSE, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("solveADjointPETSc", & 
                       "Error in KSPSetResidualHistory.")

      ! Solve the adjoint system of equations [dR/dW]T psi = dJ/dW.

      ! KSPSolveTranspose - Solves the transpose of a linear system.
      !
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPSolveTranspose(KSP ksp,Vec b,Vec x, PetscErrorCode ierr)
      !
      ! Collective on KSP
      !
      ! Input Parameter
      !   ksp - iterative context obtained from KSPCreate()
      !   b   - right hand side vector
      !   x   - solution vector
      !
      ! see .../petsc/docs/manualpages/KSP/KSPSolveTranspose.html
      ! or PETSc users manual, pp.64

      call KSPSolve(ksp,pvr,psi,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("solveADjointPETSc", &
                       "Error in KSPSolveTranspose.")

      !call VecView(psi,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
      !pause
      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, PETSC_COMM_WORLD, PETScIerr)

!     *****************************************************************
!     *                                                               *
!     * Write the Solution to a File                                  *
!     *                                                               *
!     *****************************************************************

      !Access solution and write to a file for comparison

      !Get local processor offset for global index
      call setPointersAdj(1,1,1)
      localoffset = globalcell(2,2,2)

      call VecGetArray(psi,x_array,i_x, PETScIerr)

      print *,'ix',myID,shape(x_array),i_x,localoffset
      !  Access first local entry in vector with

      do i = 1,nDimW
          vecval = x_array(i_x + i)

          write(iFile,11) -vecval
11        format(1x,'',f20.7)
      enddo

      do nn=1,nDom
         call setPointersAdj(nn,1,1)
         do k= 2, kl
            do j= 2, jl
               do i= 2, il
                  do iw = 1, nw
		    
                     idxres   = (globalCell(i,j,k)-localoffset)*nw+iw
!                     print *,'index',i,j,k,iw,nn,i_x + idxres,i_x,idxres
	             vecval = x_array(i_x + idxres)
                    ! write(iFile+1, 12) -vecval ,nn,i,j,k,iw,idxres
                     write(iFile+1+myID, 12) -vecval ,nn,i,j,k,iw,idxres+localoffset*nw
                  end do
               end do
            end do
         end do
      enddo
12        format(1x,'drdx',f20.9,6I8)



      ! ...... other code

      call VecRestoreArray(psi,x_array,i_x, PETScIerr)
      close(iFile)
      close(iFile+1+myID)

      call VecDot(psi,dJdw,PETScScalarV,PETScIerr)

      print *,'PETScScalarV',PETScScalarV
      stop
!
!     ******************************************************************
!     *                                                                *
!     * Check the solution.                                            *
!     *                                                                *
!     ******************************************************************
!
      ! Compute the residual vector r = Ax - b in two steps:
      ! step #1) r = Ax

      ! MatMultTranspose - Computes matrix transpose times a vector.
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! call MatMultTranspose(Mat mat, Vec x,Vec y, PetscErrorCode ierr)
      !
      ! Collective on Mat and Vec
      ! 
      ! Input Parameters
      !   mat - the matrix
      !   x   - the vector to be multilplied
      !
      ! Output Parameters
      !   y -the result
      !
      ! Notes
      ! The vectors x and y cannot be the same. I.e.,
      ! one cannot call MatMultTranspose(A,y,y).
      !
      ! see .../petsc/docs/manualpages/Mat/MatMultTranspose.html
      ! or PETSc users manual, pp.57

!      call MatMultTranspose(dRdW,psi,pvr,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("solveADjointPETSc","Error in MatMultTranspose.")

      ! step #2) r = r - b

      ! VecAXPY - Computes y = alpha x + y.
      !
      ! Synopsis
      !
      ! #include "petscvec.h" 
      ! call VecAXPY(Vec y,PetscScalar alpha,Vec x,PetscErrorCode ierr)
      !
      ! Collective on Vec
      !
      ! Input Parameters
      !   alpha - the scalar
      !   x, y  - the vectors
      !
      ! Output Parameter
      !   y -output vector 
      !
      ! see .../petsc/docs/manualpages/Vec/VecAXPY.html
      ! or PETSc users manual, pp.38

!      call VecAXPY(pvr,PETScNegOne,dJdW,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("solveADjointPETSc", "Error in VecAXPY.")

      ! Check the error (compute the norm of the error vector pvr)

      ! VecNorm - Computes the vector norm.
      !
      ! Synopsis
      !
      ! #include "petscvec.h" 
      ! call VecNorm(Vec x,NormType type,PetscReal *val, &
      !              PetscErrorCode ierr)
      !
      ! Collective on Vec
      !
      ! Input Parameters
      !   x    - the vector
      !   type - one of NORM_1, NORM_2, NORM_INFINITY. Also available
      !          NORM_1_AND_2, which computes both norms and stores
      !          them in a two element array.
      !
      ! Output Parameter
      !   val  - the norm 
      !
      ! see .../petsc/docs/manualpages/Vec/VecNorm.html
      ! or PETSc users manual, pp.38

!      call VecNorm(pvr,NORM_2,norm,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("solveADjointPETSc", "Error in VecNorm.")

      ! Get the iteration number when convergence or divergence
      ! was detected.

      ! KSPGetIterationNumber - Gets the current iteration number;
      !                         if the KSPSolve() is complete, returns
      !                         the number of iterations used.
      ! Synopsis
      !
      ! #include "petscksp.h" 
      ! call KSPGetIterationNumber(KSP ksp,PetscInt *its, &
      !                            PetscErrorCode ierr)
      !
      ! Not Collective
      !
      ! Input Parameters
      !   ksp - the iterative context
      !
      ! Output Parameters
      !   its - number of iterations
      !
      ! Notes: During the ith iteration this returns i-1
      !
      ! see .../petsc/docs/manualpages/KSP/KSPGetIterationNumber.html
      ! or PETSc users manual, pp.64

      call KSPGetIterationNumber(ksp,adjConvIts,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("solveADjointPETSc", & 
                       "Error in KSPGetIterationNumber.")

      ! Use the root processor to display the output summary, such as
      ! the norm of error and the number of iterations
      norm = 3
      if( PETScRank==0 ) then
        write(*,20) "Solving ADjoint with PETSc time (s) =", timeAdj
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
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)

      ! Output formats.

   10 format(a)
   20 format(a,1x,f8.2)
   30 format(1x,a,1x,e10.4,4x,a,1x,i4)
   40 format(1x,a,1x,i5,1x,a)

#endif

      end subroutine solveDirectPETSc
