
!
!     ******************************************************************
!     *                                                                *
!     * File:          setupADjointRHSStruct.f90                       *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 08-18-2008                                      *
!     * Last modified: 08-18-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine setupADjointRHSStruct(lenadjoint,structAdjoint)
!
!     ******************************************************************
!     *                                                                *
!     * Compute the right hand side Augmentation of the adjoint problem*
!     * contributed by the structures.                                 *
!     *                                                                *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use ADjointVars
      use mdData              !mdNSurfNodes
      use communication       ! procHalo(currentLevel)%nProcSend, myID

      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType),intent(in)::lenadjoint
      integer(kind=intType) :: modFamID
      real(kind=realType), intent(in),dimension(3,lenadjoint) :: structAdjoint
!      real(kind=realType), intent(in),dimension(3,nSurfNodesLocal) :: structAdjoint
!
!     Local variables.
!
      real(kind=realType), dimension(2) :: time
      real(kind=realType)               :: timeAdjLocal, timeAdj

      integer(kind=intType)::n,m,idxsurf
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
!#ifndef USE_NO_PETSC
!      !determine the number of surface nodes for coupling matrix
!      call mdCreateNSurfNodesLocal
      modFamID = max(0, 1_intType)


      ! Send some feedback to screen.

      if( PETScRank==0 ) &
        write(*,10) "Assembling ADjoint RHS Augmentation vector..."

      ! Get the initial time.

      call cpu_time(time(1))


      do m=1,nSurfNodesLocal
         do n=1,3
                 
            idxSurf = (m-1)*3+n + (mdNsurfNodes(myID,modFamID)*3)

            if (structAdjoint(n,m).ne.0.0)then
               call VecSetValues(phic, 1, idxSurf-1,  &
                    structAdjoint(n,m), INSERT_VALUES, PETScIerr)
               if( PETScIerr/=0 ) &
                    print *,'matrix setting error'
               
            endif
         enddo
      enddo
!
!     ******************************************************************
!     *                                                                *
!     * Complete the PETSc vector assembly process.                    *
!     *                                                                *
!     ******************************************************************
!
      ! VecAssemblyBegin - Begins assembling the vector. This routine
      ! should be called after completing all calls to VecSetValues().
      !
      ! Synopsis
      !
      ! #include "petscvec.h" 
      ! call VecAssemblyBegin(Vec vec, PetscErrorCode ierr)
      !
      ! Collective on Vec
      !
      ! Input Parameter
      !   vec -the vector 
      !
      ! see .../petsc/docs/manualpages/Vec/VecAssemblyBegin.html

      call VecAssemblyBegin(phic,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupADjointRHSStruct", "Error in VecAssemblyBegin")

      ! VecAssemblyEnd - Completes assembling the vector. This routine
      ! should be called after VecAssemblyBegin().
      !
      ! Synopsis
      !
      ! #include "petscvec.h" 
      ! call VecAssemblyEnd(Vec vec, PetscErrorCode ierr)
      !
      ! Collective on Vec
      !
      ! Input Parameter
      !   vec -the vector 
      !
      ! see .../petsc/docs/manualpages/Vec/VecAssemblyEnd.html

      call VecAssemblyEnd  (phic,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupADjointRHSStruct", "Error in VecAssemblyEnd")

      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Assembling ADjoint RHS Struct vector time (s) = ", timeAdj
!
!     ******************************************************************
!     *                                                                *
!     * Visualize the assembled vector.                                *
!     *                                                                *
!     ******************************************************************
!
      ! VecView - Views a vector object.
      !
      ! Synopsis
      !
      ! #include "petscvec.h" 
      ! PetscErrorCode PETSCVEC_DLLEXPORT VecView(Vec vec, &
      !                                              PetscViewer viewer)
      !
      ! Collective on Vec
      !
      ! Input Parameters
      !   v      - the vector
      !   viewer - an optional visualization context
      !
      ! Notes
      ! The available visualization contexts include
      !   PETSC_VIEWER_STDOUT_SELF  - standard output (default)
      !   PETSC_VIEWER_STDOUT_WORLD - synchronized standard output where
      !    only the first processor opens the file. All other processors
      !    send their data to the first processor to print.
      !
      ! see .../petsc/docs/manualpages/Vec/VecView.html
      ! or PETSc users manual, pp.36,148

      if( debug ) then
	call VecView(phic,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
	!call VecView(phic,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
        if( PETScIerr/=0 ) &
          call terminate("setupADjointRHSStruct", "Error in VecView")
        pause
      endif

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)

      !now multiply with the dSdw derivatives
      
      ! Compute the second contribution term phic^T dSdw, which requires
      ! a matrix-vector multiplication.

      call MatMultTranspose(dSdw,phic,dJcdw,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("SetupADjointRHSStruct", &
                       "Error in MatMultTranspose X")

      ! Add the first structural augmentation on to the base dJ/dw term.
      ! the negative sign will be taken into account in the total derivatives.

      call VecAXPY(dJdW,PETScOne,dJcdW,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupADjointRHSStruct", &
                       "Error in VecAXPY X")


      ! Output formats.

   10 format(a)
   20 format(a,1x,f8.2)

!#endif

    end subroutine setupADjointRHSStruct
