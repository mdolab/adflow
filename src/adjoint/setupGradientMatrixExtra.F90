!
!     ******************************************************************
!     *                                                                *
!     * File:          setupGradientMatrixExtra.F90                    *
!     * Author:        Andre C. Marta,C.A.(Sandy) Mader                *
!     * Starting date: 08-23-2005                                      *
!     * Last modified: 03-01-2007                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine setupGradientMatrixExtra(level)
!
!     ******************************************************************
!     *                                                                *
!     * Compute the residual partial sensitivity with respect to the   *
!     * extra design variables:                                        *
!     *                                                                *
!     *                                       dpartial R(i,j,k,nw)     *
!     *   dRda(il,jl,kl,nw,nDesignExtra) = --------------------------  *
!     *                                   dpartial alpha(nDesignExtra) *
!     *                                                                *
!     *   where nDesignExtra includes:                                 *
!     *   > angle of attack                                            *
!     *   > side slip angle                                            *
!     *   > Mach Number                                 .              *
!     *                                                                *
!     * This routine has to be consistent with the design variable     *
!     * ordering done in "setDesignVaribles".                          *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use ADjointVars
      use blockPointers   ! il,jl,kl, detJinvFlow, w, b0, dw, globalNode
      use communication   ! myID, nProc
      use flowVarRefState ! vInf
      use inputTimeSpectral ! nTimeIntervalsSpectral
      use inputDiscretization ! spaceDiscr
      use inputPhysics    ! velDirFreestream
      use iteration       ! magCoupled
     
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: level
!
!     Local variables.
!
      integer(kind=intType) :: iCell, jCell, kCell,sps
      integer(kind=intType) :: nn, m, n,liftIndex

      real(kind=realType), dimension(2) :: time
      real(kind=realType)               :: timeAdjLocal, timeAdj

      ! dR/da local block matrix at node (iNode,jNode,kNode)

      real(kind=realType), dimension(nw,ndesignextra) :: dRdaLocal

      ! idxmg - global row index
      ! idxng - global column index

      integer(kind=intType) :: idxmg(nw), idxng

      ! auxiliar variables to compute dR/dalpha and dR/dbeta

      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral) :: wAdj, wAdjB
      real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral)  :: xAdj, xAdjB

      real(kind=realType), dimension(nw,nTimeIntervalsSpectral) :: dwAdj, dwAdjB

      REAL(KIND=REALTYPE) :: machadj, machcoefadj, uinfadj, pinfcorradj
      REAL(KIND=REALTYPE) :: machadjb, machcoefadjb,machgridadj, machgridadjb
      REAL(KIND=REALTYPE) :: prefadj, rhorefadj
      REAL(KIND=REALTYPE) :: pinfdimadj, rhoinfdimadj
      REAL(KIND=REALTYPE) :: rhoinfadj, pinfadj
      REAL(KIND=REALTYPE) :: murefadj, timerefadj
      REAL(KIND=REALTYPE) :: alphaadj, betaadj
      REAL(KIND=REALTYPE) :: alphaadjb, betaadjb
      real(kind=realType), dimension(3) ::rotRateAdj,rotCenterAdj,rotrateadjb
      REAL(KIND=REALTYPE) :: xblockcorneradj(2, 2, 2, 3,nTimeIntervalsSpectral), xblockcorneradjb(2&
           &  , 2, 2, 3,nTimeIntervalsSpectral)

      logical :: secondHalo,exchangeTurb,correctfork,finegrid
      integer(kind=intType):: discr

      character(len=2*maxStringLen) :: errorMessage

      ! pvr row block
      real(kind=realType), dimension(nw) :: pvrlocal
      real(kind=realType) ::test
      integer(kind=intType) ::idxmgb
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
        write(*,10) "Assembling dR/da matrix..."

      ! Get the initial time.

      call cpu_time(time(1))

      ! Set the grid level of the current MG cycle, the value of the
      ! discretization and the logical correctForK.

      currentLevel = level
      discr        = spaceDiscr
      fineGrid     = .true.

      ! Determine whether or not the total energy must be corrected
      ! for the presence of the turbulent kinetic energy and whether
      ! or not turbulence variables should be exchanged.

      correctForK  = .false.
      exchangeTurb = .false.

     
      ! Set the value of secondHalo, depending on the situation.
      ! In the full MG (currentLevel < groundLevel) the second halo is
      ! always set; otherwise only on the finest mesh in the current mg
      ! cycle.

      if(currentLevel <= groundLevel) then
         secondHalo = .true.
      else
         secondHalo = .false.
      endif

!
!     ******************************************************************
!     *                                                                *
!     * Exchange halo data to make sure it is up-to-date.              *
!     * (originally called inside "rungeKuttaSmoother" subroutine).    *
!     *                                                                *
!     ******************************************************************
!
      ! Exchange the pressure if the pressure must be exchanged early.
      ! Only the first halo's are needed, thus whalo1 is called.
      ! Only on the fine grid.
      
      if(exchangePressureEarly .and. currentLevel <= groundLevel) &
           call whalo1(currentLevel, 1_intType, 0_intType, .true.,&
           .false., .false.)
      
      ! Apply all boundary conditions to all blocks on this level.
      
      call applyAllBC(secondHalo)
      
      ! Exchange the solution. Either whalo1 or whalo2
      ! must be called.
      
      if( secondHalo ) then
         call whalo2(currentLevel, 1_intType, nMGVar, .true., &
              .true., .true.)
      else
         call whalo1(currentLevel, 1_intType, nMGVar, .true., &
              .true., .true.)
      endif

!
!     ******************************************************************
!     *                                                                *
!     * Compute the  matrix dR/da using Tapenade's reverse mode        *
!     * of Automatic Differentiation.  NOTE: This is the reason I have *
!     * been writing the word "ADjoint" with A and D capitalized. A    *
!     * simple play with letter so that:                               *
!     *                                                                *
!     * ADjoint = Automatically Differentiated adjoint                 *
!     *                                                                *
!     ******************************************************************
!

      ! Get the initial time.

      call cpu_time(time(1))

      ! Loop over the number of local blocks.

      domainLoop: do nn=1,nDom

         ! Loop over the number of time instances for this block.

         spectralLoop: do sps=1,nTimeIntervalsSpectral

            call setPointersAdj(nn,level,sps)
           ! Loop over location of output (R) cell of residual
            do kCell = 2, kl
               do jCell = 2, jl
                  do iCell = 2, il
                     ! Copy the state w to the wAdj array in the stencil
                     call copyADjointStencil(wAdj, xAdj,xBlockCornerAdj,alphaAdj,&
           betaAdj,MachAdj,machCoefAdj,machGridAdj,iCell, jCell, kCell,&
           nn,level,sps,&
           prefAdj,rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
           rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
           murefAdj, timerefAdj,pInfCorrAdj,liftIndex)

                     mLoop: do m = 1, nw       
                        ! Loop over output cell residuals (R)

                        ! Initialize the seed for the reverse mode
                        dwAdjb(:,:) = 0.; dwAdjb(m,sps) = 1.
                        dwAdj(:,:)  = 0.
                        wAdjb(:,:,:,:,:)  = 0.  !dR(m)/dw
                        xAdjb(:,:,:,:,:)  = 0.  !dR(m)/dx
                        alphaAdjb = 0.
                        betaAdjb = 0.
                        MachAdjb = 0.
                        MachgridAdjb = 0.
		        rotrateadjb(:)=0.

                        ! Call reverse mode of residual computation
                        call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
&  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
&  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
&  icell, jcell, kcell, nn, level, sps, correctfork, secondhalo, prefadj&
&  , rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj&
&  , rotrateadjb, rotcenteradj, murefadj, timerefadj, pinfcorradj, &
&  liftindex)


                        dRdaLocal(m,nDesignAOA) =alphaAdjb
                        dRdaLocal(m,nDesignSSA) =betaAdjb
                        dRdaLocal(m,nDesignMach) =machAdjb
                        dRdaLocal(m,nDesignMachGrid) =machgridAdjb
			dRdaLocal(m,nDesignRotX) =rotrateadjb(1)
			dRdaLocal(m,nDesignRotY) =rotrateadjb(2)
			dRdaLocal(m,nDesignRotZ) =rotrateadjb(3)
                        pvrlocal(m) = machgridAdjb!machAdjb!rotrateadjb(3)*timeref
                     enddo mLoop

!!$!temporary code to setup Direct solve
!!$                     idxmgb = globalCell(icell,jcell,kcell)
!!$                     
!!$                     test = sum(pvrlocal(:))
!!$                     !print *,'test',test
!!$                     if ( test.ne.0 .and. idxmgb.ne.-5 .and. idxmgb>=0 .and. idxmgb<nCellsGlobal*nTimeIntervalsSpectral) then
!!$                        !print *,'setting PETSc Vector',sum(wAdjB(icell,jcell,kcell,:))
!!$                        !pvrlocal(:) = wFD2(iCell-1, jCell-1, kCell-1,:)
!!$                        
!!$                        !                call VecSetValuesBlocked(dJdW, 1, idxmgb, dJdWlocal, &
!!$                        !                                         INSERT_VALUES, PETScIerr)
!!$                        call VecSetValuesBlocked(pvr, 1, idxmgb, pvrlocal, &
!!$                             INSERT_VALUES, PETScIerr)
!!$                        
!!$                        if( PETScIerr/=0 ) then
!!$                           write(errorMessage,99) &
!!$                                "Error in VecSetValuesBlocked for global node", &
!!$                           idxmgb
!!$                           call terminate("setupADjointRHSAeroCoeff", &
!!$                                errorMessage)
!!$                        endif
!!$                     endif

              ! Transfer the block Jacobians to the global [dR/da]
              ! matrix by setting the corresponding block entries of
              ! the PETSc matrix dRda.
              !
              ! Global matrix column idxmg function of node indices.
              ! (note: index displaced by previous design variables)

   	      !Angle of Attack
              do m=1,nw
                idxmg(m) = globalCell(iCell,jCell,kCell) * nw + m - 1
              enddo
              idxng = nDesignAOA - 1
	      !print *,'index',idxmg,'n',idxng
              call MatSetValues(dRda, nw, idxmg, 1, idxng, &
                                dRdaLocal(:,nDesignAOA), INSERT_VALUES, PETScIerr)

              if( PETScIerr/=0 ) then
                write(errorMessage,99) &
                      "Error in MatSetValues for global column", idxng
                call terminate("setupGradientMatrixExtra", errorMessage)
              endif

              ! Side slip angle
              do m=1,nw
                idxmg(m) = globalCell(iCell,jCell,kCell) * nw + m - 1
              enddo
              idxng = nDesignSSA - 1

              call MatSetValues(dRda, nw, idxmg, 1, idxng, &
                                dRdaLocal(:,nDesignSSA), INSERT_VALUES, PETScIerr)

              if( PETScIerr/=0 ) then
                write(errorMessage,99) &
                      "Error in MatSetValues for global column", idxng
                call terminate("setupGradientMatrixExtra", errorMessage)
              endif

              !Mach Number
              do m=1,nw
                idxmg(m) = globalCell(iCell,jCell,kCell) * nw + m - 1
              enddo
              idxng = nDesignMach - 1

              call MatSetValues(dRda, nw, idxmg, 1, idxng, &
                                dRdaLocal(:,nDesignMach), INSERT_VALUES, PETScIerr)

              if( PETScIerr/=0 ) then
                write(errorMessage,99) &
                      "Error in MatSetValues for global column", idxng
                call terminate("setupGradientMatrixExtra", errorMessage)
              endif

	      !Mach Number Grid
              do m=1,nw
                idxmg(m) = globalCell(iCell,jCell,kCell) * nw + m - 1
              enddo
              idxng = nDesignMachGrid - 1

              call MatSetValues(dRda, nw, idxmg, 1, idxng, &
                                dRdaLocal(:,nDesignMachGrid), INSERT_VALUES, PETScIerr)

              if( PETScIerr/=0 ) then
                write(errorMessage,99) &
                      "Error in MatSetValues for global column", idxng
                call terminate("setupGradientMatrixExtra", errorMessage)
              endif

	      !X Rotation
              do m=1,nw
                idxmg(m) = globalCell(iCell,jCell,kCell) * nw + m - 1
              enddo
              idxng = nDesignRotX - 1

              call MatSetValues(dRda, nw, idxmg, 1, idxng, &
                                dRdaLocal(:,nDesignRotX), INSERT_VALUES, PETScIerr)

              if( PETScIerr/=0 ) then
                write(errorMessage,99) &
                      "Error in MatSetValues for global column", idxng
                call terminate("setupGradientMatrixExtra", errorMessage)
              endif

              !Y Rotation
              do m=1,nw
                idxmg(m) = globalCell(iCell,jCell,kCell) * nw + m - 1
              enddo
              idxng = nDesignRotY - 1

              call MatSetValues(dRda, nw, idxmg, 1, idxng, &
                                dRdaLocal(:,nDesignRotY), INSERT_VALUES, PETScIerr)

              if( PETScIerr/=0 ) then
                write(errorMessage,99) &
                      "Error in MatSetValues for global column", idxng
                call terminate("setupGradientMatrixExtra", errorMessage)
              endif

	      !Z Rotation
              do m=1,nw
                idxmg(m) = globalCell(iCell,jCell,kCell) * nw + m - 1
              enddo
              idxng = nDesignRotZ - 1

              call MatSetValues(dRda, nw, idxmg, 1, idxng, &
                                dRdaLocal(:,nDesignRotZ), INSERT_VALUES, PETScIerr)

              if( PETScIerr/=0 ) then
                write(errorMessage,99) &
                      "Error in MatSetValues for global column", idxng
                call terminate("setupGradientMatrixExtra", errorMessage)
              endif

           enddo
        enddo
     enddo
     
  enddo spectralLoop
enddo domainLoop


!!$call VecAssemblyBegin(pvr,PETScIerr)
!!$
!!$if( PETScIerr/=0 ) &
!!$     call terminate("setupASjointRHS", "Error in VecAssemblyBegin")  
!!$
!!$call VecAssemblyEnd  (pvr,PETScIerr)
!!$
!!$if( PETScIerr/=0 ) &
!!$     call terminate("setupADjointRHS", "Error in VecAssemblyEnd")
!!$
!!$if( debug ) then
!!$   call VecView(pvr,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
!!$   !call VecView(pvr,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
!!$   if( PETScIerr/=0 ) &
!!$        call terminate("setupADjointRHS", "Error in VecView")
!!$   pause
!!$endif

!
!     ******************************************************************
!     *                                                                *
!     * Complete the PETSc matrix assembly process.                    *
!     *                                                                *
!     ******************************************************************
!
      ! MatAssemblyBegin - Begins assembling the matrix. This routine
      !  should be called after completing all calls to MatSetValues().

      call MatAssemblyBegin(dRda,MAT_FINAL_ASSEMBLY,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupGradientMatrixExtra", &
                       "Error in MatAssemblyBegin")

      ! MatAssemblyEnd - Completes assembling the matrix. This routine
      !                  should be called after MatAssemblyBegin().

      call MatAssemblyEnd(dRda,MAT_FINAL_ASSEMBLY,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupGradientMatrixExtra", &
                       "Error in MatAssemblyEnd")

      ! Let PETSc know that the dRda matrix retains the same nonzero 
      ! pattern, in case the matrix is assembled again, as for a new
      ! point in the design space.

      ! MatSetOption - Sets a parameter option for a matrix.
      !   Some options may be specific to certain storage formats.
      !   Some options determine how values will be inserted (or added).
      !   Sorted,row-oriented input will generally assemble the fastest.
      !   The default is row-oriented, nonsorted input.

      call MatSetOption(dRda,MAT_NO_NEW_NONZERO_LOCATIONS,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupGradientMatrixExtra", &
                       "Error in MatSetOption")

      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Assembling dR/da matrix time (s) =", timeAdj
!
!     ******************************************************************
!     *                                                                *
!     * Visualize the assembled matrix.                                *
!     *                                                                *
!     ******************************************************************
!
      ! MatView - Visualizes a matrix object.

      if( debug ) then
        call MatView(dRda,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
        if( PETScIerr/=0 ) &
          call terminate("setupGradientMatrixExtra", "Error in MatView")
        pause
      endif

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)

      ! Output formats.

   10 format(a)
   20 format(a,1x,f8.2)
   99 format(a,1x,i6)

#endif

      end subroutine setupGradientMatrixExtra
