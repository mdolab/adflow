!
!     ******************************************************************
!     *                                                                *
!     * File:          verifyRAdj.f90                                  *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 04-17-2008                                      *
!     * Last modified: 04-17-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine verifyRAdj(level)
!
!     ******************************************************************
!     *                                                                *
!     * Computes the residual of the mean flow equations on the        *
!     * current level and specified time instance using the auxiliary  *
!     * routines based on the cell indices and compares it to the      *
!     * residual obtained with the flow solver routine.                *
!     *                                                                *
!     * This is only executed in debug mode.                           *
!     *                                                                *
!     ******************************************************************
!
!     Modules
      use blockPointers
      use communication !myID
      use iteration     !Currentlevel, GroundLevel
      use flowVarRefState !kPresent
      use inputTimeSpectral     !nTimeIntervalsSpectral
      use inputIO          !solfile,surfaceSolFile
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: level!, sps
!
!     Local variables.
!
      integer(kind=intType) :: iCell, jCell, kCell
      integer(kind=intType) :: mm, nn, ii,sps,liftIndex

      logical :: fineGrid, correctForK, exchangeTurb
      logical :: secondHalo!, correctForK
      

      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdj
      !real(kind=realType), dimension(-2:3,-2:3,-2:3,3)  :: xAdj
      real(kind=realType), dimension(-3:2,-3:2,-3:2,3)  :: xAdj


      real(kind=realType), dimension(nw)                :: dwAdj

      real(kind=realType) :: alphaAdj, betaAdj,MachAdj,machCoefAdj,machGridAdj
      REAL(KIND=REALTYPE) :: prefAdj, rhorefAdj,pInfCorrAdj
      REAL(KIND=REALTYPE) :: pinfdimAdj, rhoinfdimAdj
      REAL(KIND=REALTYPE) :: rhoinfAdj, pinfAdj
      REAL(KIND=REALTYPE) :: murefAdj, timerefAdj
      real(kind=realType), dimension(3) ::rotRateAdj,rotCenterAdj

      real(kind=realType), dimension(4) :: time
      real(kind=realType)               :: timeAdj, timeOri

      integer :: ierr
      real(kind=realType) :: differ

      character(len=maxStringLen) :: solFileBak, surfaceSolFileBak
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
      ! Get the initial time.

      call mpi_barrier(SUmb_comm_world, ierr)
      if(myID == 0) call cpu_time(time(1))


      if(myID == 0) then
        write(*,*) "Running verifyRAdj..."
        write(*,*) "proc block  i   j   k  sum(dwAdj)", &
                   "    sum(dw)  sum(diff)"
      endif

      ! Set the grid level of the current MG cycle, the value of the
      ! discretization and the logical fineGrid.

      currentLevel = level
!      discr        = spaceDiscr
      fineGrid     = .true.

      ! Determine whether or not the total energy must be corrected
      ! for the presence of the turbulent kinetic energy.

      if( kPresent ) then
         if((currentLevel <= groundLevel) .or. turbCoupled) then
            correctForK = .true.
         else
            correctForK = .false.
         endif
      else
         correctForK = .false.
      endif

      ! and whether or not turbulence variables should be exchanged
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
      
      ! Reset the values of rkStage and currentLevel, such that
      ! they correspond to a new iteration.

      rkStage = 0
      currentLevel = groundLevel

      ! Compute the latest values of the skin friction velocity.
      ! The currently stored values are of the previous iteration.
      
      call computeUtau
      
      ! Apply an iteration to the turbulent transport equations in
      ! case these must be solved segregatedly.
      
      if( turbSegregated ) call turbSolveSegregated

      ! Compute the time step.
      
      call timeStep(.false.)
      
      ! Compute the residual of the new solution on the ground level.
      
      if( turbCoupled ) then
         call initres(nt1MG, nMGVar)
         call turbResidual
      endif

      call initres(1_intType, nwf)
      call residual

      ! Get the final time for original routines.
      call mpi_barrier(SUmb_comm_world, ierr)
      if(myID == 0) then
        call cpu_time(time(2))
        timeOri = time(2)-time(1)
      endif

!
!     ******************************************************************
!     *                                                                *
!     * Compute the residual using the node-by-node based routine.     *
!     *                                                                *
!     ******************************************************************
!
      ! Get the initial ResidualAdj time.

      call mpi_barrier(SUmb_comm_world, ierr)
      if(myID == 0) call cpu_time(time(3))

      ! Loop over the number of local blocks.
      spectralLoop: do sps=1,nTimeIntervalsSpectral
         domainResidualLoop: do nn=1,nDom

            ! Set some pointers to make the code more readable.

            call setPointersAdj(nn,level,sps)

 
            !===============================================================
            ! Compute the residual for each cell.
            !print *,'aftersetpointeres',nn,il,jl,kl!,shape(w)
            do kCell=2,kl!2,kl
               !print *,'init il,jl,kl',il,jl,kl!,shape(w)
               do jCell=2,jl!2,jl
                  !print *,'interior il,jl,kl',il,jl,kl,nn!,shape(w)
                  do iCell=2,il
                     !print *,'interior2 il,jl,kl',il,jl,kl,nn!,shape(w)
                     ! Transfer the state w to the auxiliar array wAdj
                     ! and the coordinates x to the auxiliar array xAdj.
                     !print *,'copying Adjoint',nn,icell,il,jcell,jl,kcell,kl,shape(w)
                     call copyADjointStencil(wAdj, xAdj,alphaAdj,betaAdj,MachAdj,&
                          machCoefAdj,machGridAdj,iCell, jCell, kCell,prefAdj,&
                          rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                          rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
                          murefAdj, timerefAdj,pInfCorrAdj,liftIndex)
                     !print *,'rotrateadj,rotcenteradj',rotRateAdj,rotCenterAdj
                     !(wAdj, xAdj, iCell, jCell, kCell)
                     !print *,'interior3 il,jl,kl',il,jl,kl,nn!,shape(w)
                     ! Compute the total residual.
                     ! This includes inviscid and viscous fluxes, artificial
                     ! dissipation, and boundary conditions.                   
                     !print *,'Calling compute ADjoint'!,wadj(:,:,:,irho)!,wAdj,xAdj,dwAdj,alphaAdj,&

!                          betaAdj,MachAdj, MachCoefAdj,&
!                          iCell, jCell,  kCell, &
!                          nn,sps, correctForK,secondHalo,prefAdj,&
!                          rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
!                          rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
!                          murefAdj, timerefAdj,pInfCorrAdj,liftIndex
                     call computeRAdjoint(wAdj,xAdj,dwAdj,alphaAdj,betaAdj,MachAdj, &
                          MachCoefAdj,machGridAdj,iCell, jCell,  kCell, &
                          nn,sps, correctForK,secondHalo,prefAdj,&
                          rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                          rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
                          murefAdj, timerefAdj,pInfCorrAdj,liftIndex)

                     !print *,'interior4 il,jl,kl',il,jl,kl,nn!,shape(w)

!!$                     call computeRAdjoint(wAdj,        &
!!$                          xAdj,                       &
!!$                          dwAdj,                      &
!!$                          iCell, jCell,  kCell,       &
!!$                          nn,sps, correctForK,secondHalo)


                     ! Output the result sum(dw(:)) to debug.

                     differ = sum(dwAdj(:))-sum(dw(iCell,jCell,kCell,:))

!!$              if( differ > L2Conv ) &
!!$                write(*,10) myID, nn, iCell, jCell, kCell,               &
!!$                            sum(dwAdj(:)), sum(dw(iCell,jCell,kCell,:)), &
!!$                            differ
                     
!                     if( differ > 1e-14 ) &
                     if( abs(differ) > 1e-14 ) &
!                     if( abs(1) > 1e-14 ) &
                          write(*,10) myID, nn, iCell, jCell, kCell,               &
                          sum(dwAdj(:)), sum(dw(iCell,jCell,kCell,:)), &
                          differ,(dwAdj(2)), (dw(iCell,jCell,kCell,2)),(dwAdj(2))-(dw(iCell,jCell,kCell,2))

                     ! Store difference to output to volume solution file.
                     ! (resrho, resmom, resrhoe) have to be added to the volume
                     ! output variables in the parameter file.
                     
                     dw(iCell,jCell,kCell,:) = dw(iCell,jCell,kCell,:) - dwAdj(:)

                  enddo
               enddo
            enddo

         enddo domainResidualLoop
      enddo spectralLoop
      ! Get new time and compute the elapsed ResidualAdj time.

      call mpi_barrier(SUmb_comm_world, ierr)
      if(myID == 0) then
        call cpu_time(time(4))
        timeAdj = time(4)-time(3)
      endif

      ! Output elapsed time for the adjoint and FD computations.

      if( myID==0 ) then
        ! print *,'times',time
        print *, "====================================================="
        print *, " Time for original   residual =", timeOri
        print *, " Time for node-based residual =", timeAdj
        print *, " Factor                       =", timeAdj/timeOri
        print *, "====================================================="
      endif

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(SUmb_comm_world, ierr)
!
!     ******************************************************************
!
      ! Write the volume and surface solution of the residual.

      ! > store original filenames.

      solFileBak        = solFile
      surfaceSolFileBak = surfaceSolFile

      ! > determine the volume and surface file names based on
      !   the corresponding flow file names.

      write(solFile,99) "verifyRAdj", trim(solFile)
      write(surfaceSolFile,99) "verifyRAdj", trim(surfaceSolFile)

      ! > write residual to file.

      !call writeSol

      ! > restore original filenames.

      solFile        = solFileBak
      surfaceSolFile = surfaceSolFileBak

      ! Output formats.

  10  format(1x,i3,2x,i3,2x,i3,1x,i3,1x,i3,2x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3)
  99  format(2a)

      end subroutine verifyRAdj
