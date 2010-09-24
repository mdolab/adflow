!
!     ******************************************************************
!     *                                                                *
!     * File:          solverADjoint.f90                               *
!     * Authors:       Juan J. Alonso, Joaquim R. R. A. Martins,       *
!     *                C.A.(Sandy) Mader                               *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 01-17-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine solverADjoint
!
!     ******************************************************************
!     *                                                                *
!     * This subroutine is the main driver for the discrete adjoint    *
!     *  solver for SUmb.  The approach to the setup and               *
!     *  solution of the adjoint problem is based on the use of        *
!     *  automatic differentiation techniques (Tapenade from INRIA)    *
!     *  applied to a version of the residual computation for a single *
!     *  cell in the SUmb mesh.  This residual computation is meant to *
!     *  be in all ways identical to that performed by the actual flow *
!     *  solver (no approximations), including all flux routines,      *
!     *  limiters (if used) and all boundary conditions.               *
!     *                                                                *
!     * This version of the residual computation is differentiated (in *
!     *  reverse mode) using Tapenade and is used to compute the exact *
!     *  matrix for the discrete adjoint problem.  No approximations   *
!     *  are made.  In principle, this approach is able to handle      *
!     *  flow solver discretizations of arbitrary complexity in an     *
!     *  automatic fashion and can therefore dramatically reduce the   *
!     *  time of development of adjoint solvers while ensuring that    *
!     *  the sensitivity information obtained is exactly consistent    *
!     *  with (infinitely accurate) finite differencing of the solver. *
!     *                                                                *
!     * For efficiency reasons, this version of the adjoint solver     *
!     *  stores the adjoint matrix.  This is a significant overhead    *
!     *  (approximately an order of magnitude) in memory with respect  *
!     *  to the original solver.  We have decided to put up with this  *
!     *  for the time being as the approach has the potential to be    *
!     *  used, exactly, for governing equations as complex as entire   *
!     *  multiple equation turbulence models and even complex MHD      *
!     *  formulations.                                                 *
!     *                                                                *
!     * The resulting linear system is solved using Krylov subspace    *
!     *  methods embedded in the PETSc toolkit from Argonne National   *
!     *  laboratories.                                                 *
!     *                                                                *
!     * This is an implementation on the cell-centered version of      *
!     *  the SUmb solver.  It is not intended to be comprehensive.     *
!     *  The intent is to support the entire SUmb vertex solver once   *
!     *  it is completed and operational.                              *
!     *                                                                *
!     * This driver assumes that the flow solution has already been    *
!     *  computed and is stored in the block module. It also assumes   *
!     *  that PETSc has already been initialized.                      *
!     *                                                                *
!     * The adjoint solver is applied only to the finest grid level    *
!     *  and first time instance.                                      *
!     *                                                                *
!     ******************************************************************
!
      use ADjointVars
      use ADjointPETSc    !petsc_comm_world,petsc_rank
      use communication   ! myID, nProc
      use cgnsGrid        ! cgnsNFamilies
      use iteration       ! groundLevel
      use inputTimeSpectral !nIntervalsTimeSpectral
      use inputTSStabDeriv !TSstability
      implicit none
!
!     Local variables.
!
      integer(kind=intType) :: level, sps ,ierr,famID
      integer(kind=intType) :: costFunction,cfstart,cfend

      real(kind=realType)   :: CL, CD, Cfx,Cfy,Cfz, CMx, CMy, CMz
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
     if(myID == 0) then
	print *,'In SolverADjoint...'
     endif
#ifdef USE_NO_PETSC

      ! Sources compiled without PETSc support. Since PETSc is required
      ! to solve the ADjoint system of equations, the adjoint solver
      ! cannot be run.

      call terminate("solverADjoint", &
                     "Routine should not be called if no PETSc support &
                     &is selected.")

#else
      ! Set the relevant grid level and time instance.

      level = 1 ! finest grid level
      sps   = 1 ! first time instance

      groundLevel = level
!
!     ******************************************************************
!     *                                                                *
!     * Cost function values: aerodynamic coefficients.                *
!     *                                                                *
!     ******************************************************************
!
      do sps  = 1,nTimeIntervalsSpectral
         call computeAeroCoef(CL,CD,Cfx,Cfy,Cfz,CMx,CMy,CMz,level,sps)

         ! Write the cost function values; only processor 0 does this.
         
         if(myID == 0) then
            write(*,*) "Cost function values: sps:",sps
            write(*,*) " CL  =", CL
            write(*,*) " CD  =", CD
            write(*,*) " CFx =", CFx   
            write(*,*) " CFy =", CFy
            write(*,*) " CFz =", CFz
            write(*,*) " CMx =", CMx
            write(*,*) " CMy =", CMy
            write(*,*) " CMz =", CMz
            write(*,*)
         endif
      enddo
!
!     ******************************************************************
!     *                                                                *
!     * Preprocessing tasks.                                           *
!     *                                                                *
!     ******************************************************************
!
      ! Assertion testing, memory allocation and global number indexing.

      call preprocessingADjoint(level)


      ! Determine the number of design variables, initialize the arrays
      ! that store the cost functions (values,names,gradients) and
      ! design variables (values,names,lower and upper bounds).

      call designInit


      ! Initialize PETSc.

      call initializePETSc

      !call PetscOptionsSetValue('-malloc_debug',PETSC_NULL_CHARACTER,PETScIerr)
      !initializeWarping
      !begin execution
      if(cgnsNfamilies > 0) then
         famID = 1
      else
         famID = 0
      endif
      call initializeWarping(famID)
    
      ! Create all the necessary PETSc objects.

      call createPETScVars


      ! Perform some verifications if in DEBUG mode.
      !moved after PETSc initialization because PETsc now included in debugging...

      !if( debug ) then

        ! Verify the node-based residual routine.

!      call verifyRAdj(level)
!stop
!      call verifyResiduals(level)
!return
!stop
!!$        ! Verify the node-based ADjoint residual routine.
!!$
!!$!	call verifydRdW(level,sps)
!      call verifydRdWFile(level,sps)
!!$!        call verifydRdwFileFD(level)
!     call verifydRdxFile(level)
!return
      call verifydRdxsFile
!!$!      call verifydRdxFileFD(level)
!!$!stop
return
!!$        ! Verify the dRdx routine
!!$
!!$! 	call verifydRdx(level,sps)
!!$!stop	
!      call verifydRdExtra(level)
!!$	call verifydRdExtraFDAD(level)	
!!$
!!$        ! Verify the ADjoint routine for the forces
!!$
!!$
!         call verifyForcesAdj(level) 
!         call verifyTSStabilityDerivAdj(level)
!return
!!$!stop	
!!$        
!!$	! Verify the force derivatives
!!$
!!$	
!call verifydCfdx(level)
!!$
!!$!stop
!!$	! Verify the force derivatives
!	call verifydCfdwfile(level)	
!!$!        call verifydCfdw(level)
!!$!stop
!      call verifydIdwfile(level)
!return
!!$	!verify the coupling derivatives
!!$	!call verifyForceCouplingAdj(level)
!!$	!call verifydSdw(level)
!!$	!call verifydSdx(level)
!!$
!!$!stop
!!$	!print *,'Going to call verifydSdx'
!!$	!call verifydSdx(level)
!!$	!print *,'Done VerifydSdx'
!!$	!stop
!!$	!print *,'Going to call verifydWdx'
!!$	!call verifydSdw(level)
!!$	!print *,'Done VerifydWdx'
!!$	!stop
!!$
!!$!stop
!!$!return


  !    endif


!
!     ******************************************************************
!     *                                                                *
!     * Set up the ADjoint matrix dR/dW using Tapenade AD routines.    *
!     * => dRdW(il,jl,kl,nw,nw)                                        *
!     *                                                                *
!     ******************************************************************
!
      !print *,'calling setupADjointMatrix'
      !call setupADjointMatrix(level)
      !call setupADjointMatrixTranspose(level)
      call setupAllResidualMatrices(level)
!return
!stop

      ! Reordered for ASM preconditioner
      ! Create the Krylov subspace linear solver context,
      ! the preconditioner context, and set their various options.

      !call createPETScKsp(level)
      call setupPETScKsp(level)
      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        print "(a)", "# ... Krylov subspace created;"



!
!     ******************************************************************
!     *                                                                *
!     * Set up the residual sensitivity w.r.t. design variables dR/da. *
!     * => dRda(il,jl,kl,nw,nDesign)                                   *
!     *                                                                *
!     ******************************************************************
!

     ! call setupGradientMatrixExtra(level)

     ! call setupGradientMatrixSpatial(level)

      call setupVolumeSurfaceDerivativesDV

!      call setupSurfaceDirect
!
!      call solveDirectPETSc
 !  return
!!$!stop
!
!     ******************************************************************
!     *                                                                *
!     * Compute the sensitivity of each cost function J.               *
!     *                                                                *
!     ******************************************************************
!
      if(TSStability)then
         cfstart = 9
         cfend = 16
      else
         cfstart =8
         cfend = 8!8
      end if
!      functionLoop: do costFunction = 1,1!10,10!1, 1!nCostFunction
      functionLoop: do costFunction = cfstart,cfend

        !***************************************************************
        !                                                              *
        ! Set up the RHS adjoint vector dJ/dW, ie, the cost function J *
        ! sensitivity w.r.t. W.                                        *
        ! => dJ/dW(il,jl,kl,nw)                                        *
        !                                                              *
        !***************************************************************
         !print *,'costfunction',costfunction 
        !call setupADjointRHS(level,sps,costFunction)
        call setupADjointRHS(level,costFunction)

        ! Solve the discrete ADjoint problem using PETSc's Krylov
        ! solver and preconditioner.
        ! => psi

        !call solveADjointPETSc
        call solveADjointTransposePETSc

        !***************************************************************
        !                                                              *
        ! Set up the cost function sensitivity w.r.t. design variables.*
        ! => dJ/da(nDesign)                                            *
        !                                                              *
        !***************************************************************

        call setupGradientRHSFlow(level,costFunction)
        call setupGradientRHSVolume(level,costFunction)

        ! Compute the total sensitivity dIda(nDesignExtra) and
        ! store it in the array functionGrad(nCostFunc,nDesignExtra)
        ! for the extra design variables in the root processor.

        call computeADjointGradientExtra(costFunction)

	! Compute the total sensitivity dIdx(ndesignSpatial) store in
	! it in the array spatialGrad(nCostFunc,nDesignSpatial) for the 
	! spatial design variables

        call computeADjointGradientSpatial(costFunction)
        !print *,'computing surface'
        call computeADjointGradientSurfaceDV(costFunction)
!!$
!!$        ! Write the adjoint field solution, the convergence history and
!!$        ! the cost function total sensitivity to file/screen.
!!$
!!$        call writeADjoint(level,sps,costFunction)
!!$
      enddo functionLoop

!
!     ******************************************************************
!
	
      ! Destroy all the PETSc objects.

      call destroyPETScVars
      
      ! Finalize PETSc.
      !print *,'finalizing petsc'

      call finalizePETSc
      
      !print *,'petsc Finalized'
      ! Release memory allocated in initDesign.
      
      if(allocated(functionName )) deallocate(functionName )
      if(allocated(functionValue)) deallocate(functionValue)
      if(allocated(functionGrad )) deallocate(functionGrad )
      
      !print *,' new additions'
      !new additions
      if(allocated(functionGradSurfaceDV)) deallocate(functionGradSurfaceDV)
      if(allocated(functionGradSurfaceDisp)) deallocate(functionGradSurfaceDisp)
      if(allocated(adjoint)) deallocate(adjoint)
      if(allocated(functionGradSpatial)) deallocate(functionGradSpatial)
      if(allocated(functionGradStruct)) deallocate(functionGradStruct)
      if(allocated(functionGradCoupling)) deallocate(functionGradCoupling)
      if(allocated(functionGradCouplingExp)) deallocate(functionGradCouplingExp)

      !print *,'xdesign'
      if(allocated(xDesignVarName )) deallocate(xDesignVarName )
      if(allocated(xDesignVar     )) deallocate(xDesignVar     )
      if(allocated(xDesignVarLower)) deallocate(xDesignVarLower)
      if(allocated(xDesignVarUpper)) deallocate(xDesignVarUpper)
   
      ! Release the memory of the adjoint variables.
      !print *,'releaseMemADjoint'
      !this function causes segmentation fault in block splitting mode
      do sps  = 1,nTimeIntervalsSpectral
         call releaseMemADjoint(level,sps)
      end do
      ! Output formats.

   10 format(/,3x,a)
   20 format(1x,i3,3x,e10.3)

#endif

      end subroutine solverADjoint
