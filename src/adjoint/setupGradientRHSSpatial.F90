!
!     ******************************************************************
!     *                                                                *
!     * File:          setupGradientRHSSpatial.F90                     *
!     * Author:        C.A.(Sandy)Mader, Andre C. Marta                *
!     * Starting date: 02-06-2007                                      *
!     * Last modified: 06-15-2007                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine setupGradientRHSSpatial(level,costFunction,sps)
!
!     ******************************************************************
!     *                                                                *
!     * Compute the cost function partial sensitivity with respect to  *
!     * the grid coordinates:                                          *
!     *                                                                *
!     *                                   dpartial J                   *
!     *   dJdx(3*nNodesGlobal)  = --------------------------           *
!     *                           dpartial x(3*nNodesGlobal)           *
!     *                                                                *
!     * Supported cost functions J:                                    *
!     *   - lift coefficient Cl if "costFunction" = costFuncLiftCoef   *
!     *   - drag coefficient Cd if "costFunction" = costFuncDragCoef   *
!     *   - x-force coef.   Cfx if "costFunction" = costFuncForceXCoef *
!     *   - y-force coef.   Cfy if "costFunction" = costFuncForceYCoef *
!     *   - z-force coef.   Cfz if "costFunction" = costFuncForceZCoef *
!     *   - x-moment coef.  Cmx if "costFunction" = costFuncMomXCoef   *
!     *   - y-moment coef.  Cmy if "costFunction" = costFuncMomYCoef   *
!     *   - z-moment coef.  Cmz if "costFunction" = costFuncMomZCoef   *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      use ADjointVars
      use blockPointers       ! il,jl,kl, globalNode
      use cgnsGrid            ! cgnsDoms
      use communication       ! myID, nProc
      use inputPhysics        ! liftDirection, dragDirection
      use flowVarRefState     ! nw
      use inputDiscretization ! spaceDiscr, useCompactDiss
      use inputTimeSpectral   ! nTimeIntervalsSpectral
      use iteration           ! overset, currentLevel
      use monitor          ! timeUnsteadyRestart
      use section          ! nSections,section%
      use bcTypes             !imin,imax,jmin,jmax,kmin,kmax
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: level, costFunction
!
!     Local variables.
!
      integer(kind=intType) :: discr, nHalo,sps,nnn
      integer(kind=intType) :: iNode, jNode, kNode, mm, nn, m, n
      integer(kind=intType) :: ii, jj, kk, i1, j1, k1, i2, j2, k2

      integer(kind=intType) ::  i2Beg,  i2End,  j2Beg,  j2End
      integer(kind=intType) :: iiBeg, iiEnd, jjBeg, jjEnd
      integer(kind=intType) :: i,j,k,l,liftIndex

      logical :: correctForK

      real(kind=realType) :: Cl, Cd, Cfx, Cfy, Cfz, Cmx, Cmy, Cmz
      real(kind=realType),dimension(nTimeIntervalsSpectral) :: ClAdj,   CdAdj,            &
                             CfxAdj,  CfyAdj,  CfzAdj,  &
                             CmxAdj,  CmyAdj,  CmzAdj
      real(kind=realType),dimension(nTimeIntervalsSpectral) :: ClAdjB,  CdAdjB,           &
                             CfxAdjB, CfyAdjB, CfzAdjB, &
                             CmxAdjB, CmyAdjB, CmzAdjB 
  
      real(kind=realType), dimension(:,:,:,:), allocatable :: xAdj,xAdjB
      real(kind=realType), dimension(:,:,:,:), allocatable :: wAdj,wAdjB
      real(kind=realType), dimension(:,:,:), allocatable :: pAdj

      REAL(KIND=REALTYPE) :: machadj, machcoefadj, uinfadj, pinfcorradj
      REAL(KIND=REALTYPE) :: machadjb, machcoefadjb,machgridadj,machgridadjb
      REAL(KIND=REALTYPE) :: prefadj, rhorefadj
      REAL(KIND=REALTYPE) :: pinfdimadj, rhoinfdimadj
      REAL(KIND=REALTYPE) :: rhoinfadj, pinfadj
      REAL(KIND=REALTYPE) :: murefadj, timerefadj
      REAL(KIND=REALTYPE) :: alphaadj, betaadj
      REAL(KIND=REALTYPE) :: alphaadjb, betaadjb
      REAL(KIND=REALTYPE) :: rotcenteradj(3), rotrateadj(3), rotrateadjb(3)

      real(kind=realType), dimension(3) :: cFpAdj, cFvAdj
      real(kind=realType), dimension(3) :: cMpAdj, cMvAdj

      real(kind=realType), dimension(3) :: cFp, cFv
      real(kind=realType), dimension(3) :: cMp, cMv

      real(kind=realType) :: factI, factJ, factK, tmp,test

      integer(kind=intType), dimension(0:nProc-1) :: offsetRecv

      real(kind=realType), dimension(2) :: time
      real(kind=realType)               :: timeAdjLocal, timeAdj

      real(kind=realType), dimension(:,:,:), pointer :: norm
      real(kind=realType), dimension(:,:,:),allocatable:: normAdj
      real(kind=realType), dimension(3) :: refPoint
      real(kind=realType) :: yplusMax

      logical :: contributeToForce, viscousSubface,righthanded,secondHalo
      logical :: finegrid,exchangeturb

      real(kind=realType), dimension(nSections) :: t

      ! dJ/dx local vector at node (iNode,jNode,kNode)

      real(kind=realType),dimension(3) :: dJdxLocal

      ! idxmgb - global block row index

      integer(kind=intType) :: idxmgb

      character(len=2*maxStringLen) :: errorMessage

      integer :: ierr

      ! auxiliar variables to compute dJ/dx using FD

      real(kind=realType) :: xBase, xPert
      real(kind=realType) :: deltaFunc
      real(kind=realType) :: clBase, cdBase, cmxBase, cmyBase, cmzBase
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
        write(*,10) "Assembling dJ/dx vector..."

      ! Get the initial time.

      call cpu_time(time(1))

      ! Reset the RHS vector dJ/dx by assigning the value zero to all
      ! its components.

      ! VecSet - Sets all components of vector to a single scalar value.

      call VecSet(dJdx,PETScZero,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupGradientRHSSpatial", "Error in VecSet X")

!################## FD CODE - START ################################
      go to 999
!
!     ******************************************************************
!     *                                                                *
!     * Store the baseline values.                                     *
!     *                                                                *
!     ******************************************************************
!
      ! Compute the baseline cost function values.

      call computeAeroCoef(clBase,cdBase,cmxBase,cmyBase,cmzBase,&
                           level,sps)
!
!     ******************************************************************
!     *                                                                *
!     * Perturb X at node.                                             *
!     *                                                                *
!     ******************************************************************
!
      ! Loop over the number of local blocks.

      domainLoopFD: do nn=1,nDom

        ! Set some pointers to make the code more readable.

        call setPointers(nn,level,sps)

        ! Loop over the owned nodes of the block.

        kNodeLoop: do kNode=0,ke
          jNodeLoop: do jNode=0,je
            iNodeLoop: do iNode=0,ie

              do n = 1,3

                ! Store baseline value

                xBase = x(iNode,jNode,kNode,n)

                ! Perturb x (design variable).

                if( abs(xBase).gt.adjEpsFd ) then
                  xPert = ( one + adjRelFd ) * xBase
                else
                  xPert = adjAbsFd * sign(one, xBase)
                endif

                x(iNode,jNode,kNode,n) = xPert

                ! Recompute the metric.

                call metric(level, .true.)

                ! Compute the perturbed cost function values.

                call computeAeroCoef(cl,cd,cmx,cmy,cmz,level,sps)

                ! Compute and store the sensitivity value.
                ! Only root processor does this.

                if( myID==0 ) then

                  ! Compute the sensitivity by finite-differences.

                  select case (costFunction)
                    case(costFuncLiftCoef)
                      deltaFunc = cl - clBase
                    case(costFuncDragCoef)
                      deltaFunc = cd - cdBase
                    case(costFuncMomXCoef)
                      deltaFunc = cmx - cmxBase
                    case(costFuncMomYCoef)
                      deltaFunc = cmy - cmyBase
                    case(costFuncMomZCoef)
                      deltaFunc = cmz - cmzBase
                    case default
                      write(errorMessage,99) &
                         "Invalid cost function", costFunction
                      call terminate("setupGradientRHSExtra", &
                         errorMessage)
                  end select

                  if( abs(xBase).gt.adjEpsFd ) then
                    dJdxLocal = deltaFunc / ( adjRelFd * xBase )
                  else
                    dJdxLocal = deltaFunc / ( xPert - xBase )
                  endif

                  ! Transfer data to PETSc vector

                  idxmgb = globalNode(iNode,jNode,kNode)

!###                  select case (n)
!###                    case(1)
!###                      call VecSetValue(dJdx, idxmgb, dJdxLocal, &
!###                                       INSERT_VALUES, PETScIerr)
!###                    case(2)
!###                      call VecSetValue(dJdy, idxmgb, dJdxLocal, &
!###                                       INSERT_VALUES, PETScIerr)
!###                    case(3)
!###                      call VecSetValue(dJdz, idxmgb, dJdxLocal, &
!###                                       INSERT_VALUES, PETScIerr)
!###                    case default
!###                      stop
!###                  end select

                  if( PETScIerr/=0 ) then
                    write(errorMessage,99) &
                      "Error in VecSetValue for global node", idxmgb
                    call terminate("setupGradientRHSSpatial", &
                                   errorMessage)
                  endif

                endif

                ! Restore baseline x.

                x(iNode,jNode,kNode,n) = xBase

              enddo ! n

            enddo iNodeLoop
          enddo jNodeLoop
        enddo kNodeLoop

      enddo domainLoopFD

      ! Adjust the metric.

      call metric(level, .true.)

  999 continue
!################## FD CODE - END ################################


      ! Set the grid level of the current MG cycle, the value of the
      ! discretization and the logical correctForK.

      currentLevel = level
      discr        = spaceDiscr
      correctForK  = .false.
      fineGrid     = .true.
      exchangeTurb = .false.
      secondhalo = .true.


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

      call mpi_barrier(SUmb_comm_world, ierr)
      if( myID==0 ) call cpu_time(time(1))
!
!     ******************************************************************
!     *                                                                *
!     * Compute the d(forces)/dx (partial) using the tapenade routines.*
!     *                                                                *
!     ******************************************************************
         
      !*********************
      !from ForcesAndMoments.f90
      ! Determine the reference point for the moment computation in
      ! meters.

      refPoint(1) = LRef*pointRef(1)
      refPoint(2) = LRef*pointRef(2)
      refPoint(3) = LRef*pointRef(3)

      ! Initialize the force and moment coefficients to 0 as well as
      ! yplusMax.
       
      dJdxLocal(:) = zero
      ClAdj = 0
      CDAdj = 0
      CmxAdj = 0
      CmyAdj = 0
      CmzAdj = 0
      CfxAdj = 0
      CfyAdj = 0
      CfzAdj = 0

      yplusMax = zero

!for the spectral case can only use one at a time anyway..
!      spectralLoopAdj: do sps=1,nTimeIntervalsSpectral
         
         !print *,'zeroing force output'
         !zero the force components
         cFpAdj(1) = zero; cFpAdj(2) = zero; cFpAdj(3) = zero
         cFvAdj(1) = zero; cFvAdj(2) = zero; cFvAdj(3) = zero
         cMpAdj(1) = zero; cMpAdj(2) = zero; cMpAdj(3) = zero
         cMvAdj(1) = zero; cMvAdj(2) = zero; cMvAdj(3) = zero
         !***********************************
         
         ! Loop over the number of local blocks.
         
         domainLoop: do nn=1,nDom
      
	! Set some pointers to make the code more readable.

        call setPointersAdj(nn,level,sps)

        allocate(xAdj(0:ie,0:je,0:ke,3), stat=ierr)
        if(ierr /= 0)                              &
             call terminate("Memory allocation failure for xAdj.")

        allocate(xAdjB(0:ie,0:je,0:ke,3), stat=ierr)
        if(ierr /= 0)                              &
             call terminate("Memory allocation failure for xAdjB.")
        
        allocate(wAdj(0:ib,0:jb,0:kb,nw), stat=ierr)
        if(ierr /= 0)                              &
             call terminate("Memory allocation failure for wAdj.")

        allocate(wAdjB(0:ib,0:jb,0:kb,nw), stat=ierr)
        if(ierr /= 0)                              &
             call terminate("Memory allocation failure for wAdjB.")
        
        allocate(pAdj(0:ib,0:jb,0:kb), stat=ierr)
        if(ierr /= 0)                              &
             call terminate("Memory allocation failure for pAdj.")
        
        righthanded = flowDoms(nn,level,sps)%righthanded


        ! Copy the coordinates into xAdj and
        ! Compute the face normals on the subfaces
        call  copyADjointForcesStencil(wAdj,xAdj,alphaAdj,betaAdj,&
           MachAdj,machCoefAdj,machGridAdj,prefAdj,rhorefAdj, pinfdimAdj,&
           rhoinfdimAdj,rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,murefAdj,&
           timerefAdj,pInfCorrAdj,nn,level,sps,liftIndex)
                
        bocoLoop: do mm=1,nBocos

           ! Determine the range of cell indices of the owned cells
           ! Notice these are not the node indices
           iiBeg = BCData(mm)%icBeg
           iiEnd = BCData(mm)%icEnd
           jjBeg = BCData(mm)%jcBeg
           jjEnd = BCData(mm)%jcEnd
           
           i2Beg= BCData(mm)%inBeg+1; i2End = BCData(mm)%inEnd
           j2Beg= BCData(mm)%jnBeg+1; j2End = BCData(mm)%jnEnd
        
           ! Zero the seeds and then select one based on the cost function

           ClAdjB  = zero
           CDAdjB  = zero
           CfxAdjB = zero
           CfyAdjB = zero
           CfzAdjB = zero
           CmxAdjB = zero
           CmyAdjB = zero
           CmzAdjB = zero
           
           select case (costFunction)
           case (costFuncLiftCoef)
              ClAdjB(sps) = 1
           case (costFuncDragCoef)
              CdAdjB(sps) = 1
           case (costFuncForceXCoef)
              CfxAdjB(sps) = 1
           case (costFuncForceYCoef)
              CfyAdjB(sps) = 1
           case (costFuncForceZCoef)
              CfzAdjb(sps) = 1
           case (costFuncMomXCoef)
              CmxAdjB(sps) = 1
           case (costFuncMomYCoef)
              CmyAdjB(sps) = 1
           case (costFuncMomZCoef)
              CmzAdjB(sps) = 1 
           end select

           ! Initialize output.
           
           xAdjB(:,:,:,:) = zero ! > return dCf/dx
           !===============================================================
           ! Compute the force derivatives
           if(equationMode == timeSpectral) then
              do nnn=1,nSections
                 !t(nnn) = t(nnn) + (sps2-1)*sections(nnn)%timePeriod &
                 !     /         real(nTimeIntervalsSpectral,realType)
                 t(nnn) = t(nnn) + (sps-1)*sections(nnn)%timePeriod &
                      /         (nTimeIntervalsSpectral*1.0)!to make denomenator a real number...
              enddo
           endif

           call COMPUTEFORCESADJ_B(xadj, xadjb, wadj, wadjb, padj, iibeg, &
&  iiend, jjbeg, jjend, i2beg, i2end, j2beg, j2end, mm, cfxadj, cfxadjb&
&  , cfyadj, cfyadjb, cfzadj, cfzadjb, cmxadj, cmxadjb, cmyadj, cmyadjb&
&  , cmzadj, cmzadjb, yplusmax, refpoint, cladj, cladjb, cdadj, cdadjb, &
&  nn, level, sps, cfpadj, cmpadj, righthanded, secondhalo, alphaadj, &
&  alphaadjb, betaadj, betaadjb, machadj, machadjb, machcoefadj, &
&  machcoefadjb, machgridadj, machgridadjb, prefadj, rhorefadj, &
&  pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, murefadj, timerefadj, &
&  pinfcorradj, rotcenteradj, rotrateadj, rotrateadjb, liftindex, t)
 

          ! Loop over cells to store the jacobian

          do k = 0,ke
            do j = 0,je
              do i = 0,ie
                
                !*******************************************************    
                !                                                      *
                ! Store the Values in PETSc.                           *
                !                                                      *
                !*******************************************************

                dJdxLocal(:) = xAdjB(i,j,k,:)

                !*******************************************************
                !                                                      *
                ! Transfer the block Jacobians to the PETSc vector.    *
                !                                                      *
                !*******************************************************

                ! When using the block vector format, one can insert
                ! elements more efficiently using the block variant.

                ! VecSetValuesBlocked - Inserts or adds blocks of values
                !                    into certain locations of a vector.
                ! Synopsis
                ! 
                ! #include "petscvec.h" 
                ! call VecSetValuesBlocked(Vec x,PetscInt ni,          &
                !            const PetscInt ix[],const PetscScalar y[],&
                !            InsertMode iora,PetscErrorCode ierr) 
                !
                ! Not Collective
                !
                ! Input Parameters
                !   x    - vector to insert in
                !   ni   - number of blocks to add
                !   ix   - indices where to add in block count,
                !          rather than element count
                !   y    - array of values
                !   iora - either INSERT_VALUES or ADD_VALUES,
                !          where ADD_VALUES adds values to any existing
                !          entries, and INSERT_VALUES replaces existing
                !          entries with new values
                ! Notes
                ! VecSetValuesBlocked() sets x[bs*ix[i]+j] = y[bs*i+j],
                !   for j=0,...,bs, for i=0,...,ni-1. where bs was set
                !   with VecSetBlockSize().
                !
                ! Calls to VecSetValuesBlocked() with the INSERT_VALUES
                !   and ADD_VALUES options cannot be mixed without
                !   intervening calls to the assembly routines.
                !
                ! These values may be cached, so VecAssemblyBegin() and
                !   VecAssemblyEnd() MUST be called after all calls to
                !   VecSetValuesBlocked() have been completed.
                !
                ! VecSetValuesBlocked() uses 0-based indices in Fortran
                !   as well as in C. 
                !
                ! see .../petsc/docs/manualpages/Vec/VecSetValuesBlocked.html

                ! Global vector block row mgb function of node indices.
                !
                ! VecSetValuesBlocked() uses 0-based row numbers but the
                ! global node numbering already accounts for that since
                ! it starts at node 0.

                idxmgb = globalNode(i,j,k)
                !print *,'indices',idxmgb,nNodesGlobal
                if (idxmgb<nNodesGlobal.and.idxmgb>=0)then
                   test = sum(dJdxLocal(:))
                   if (test/=0)then
                      ! note: Add_values is used here to accumulate over the subfaces
                      ! possible to add an if != 0 clause to speed up????
                      
                      call VecSetValuesBlocked(dJdx, 1, idxmgb, dJdxLocal, &
                           ADD_VALUES, PETScIerr)
                      
                      if( PETScIerr/=0 ) then
                         write(errorMessage,99) &
                              "Error in VecSetValuesBlocked for global node", &
                              idxmgb
                         call terminate("setupGradientRHSSpatial", &
                              errorMessage)
                      endif
                   endif
                endif
              enddo 
            enddo 
          enddo 

	enddo bocoLoop

        !===============================================================
        
        ! Deallocate the xAdj.

                ! Deallocate the xAdj.
        deallocate(pAdj, stat=ierr)
        if(ierr /= 0)                              &
             call terminate("verifydCfdx", &
             "Deallocation failure for xAdj.")
             
        ! Deallocate the xAdj.
        deallocate(wAdj, stat=ierr)
        if(ierr /= 0)                              &
             call terminate("verifydCfdx", &
             "Deallocation failure for xAdj.") 

         deallocate(wAdjB, stat=ierr)
        if(ierr /= 0)                              &
             call terminate("verifydCfdx", &
             "Deallocation failure for xAdj.") 
        ! Deallocate the xAdj.
        deallocate(xAdj, stat=ierr)
        if(ierr /= 0)                              &
             call terminate("verifydCfdx", &
             "Deallocation failure for xAdj.") 

         deallocate(xAdjB, stat=ierr)
        if(ierr /= 0)                              &
             call terminate("verifydCfdx", &
             "Deallocation failure for xAdj.") 
 
     enddo domainLoop
!see above
!  enddo spectralLoopAdj
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

      call VecAssemblyBegin(dJdx,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupGradientRHSSpatial", &
                       "Error in VecAssemblyBegin")

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

      call VecAssemblyEnd  (dJdx,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupGradientRHSSpatial", &
                       "Error in VecAssemblyEnd X")

      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, SUMB_PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Assembling dJ/dspatial vector time (s) =", timeAdj
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
        call VecView(dJdx,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
        if( PETScIerr/=0 ) &
          call terminate("setupGradientRHSSpatial", "Error in VecView")
        pause
      endif

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

      ! Output formats.

   10 format(a,1x,i3,1x,a,1x,i3)
   20 format(a,1x,f8.2)
   99 format(a,1x,i6)

#endif

    end subroutine setupGradientRHSSpatial
