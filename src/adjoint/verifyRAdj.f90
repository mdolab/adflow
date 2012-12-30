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
  use adjointvars
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

  real(kind=realType), dimension(20) :: timings
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral) :: wAdj
  !real(kind=realType), dimension(-2:3,-2:3,-2:3,3)  :: xAdj
  real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral)  :: xAdj


  real(kind=realType), dimension(nw,nTimeIntervalsSpectral)                :: dwAdj

  real(kind=realType), dimension(2,2,2,3,nTimeIntervalsSpectral) ::xBlockCornerAdj
  real(kind=realType) :: alphaAdj, betaAdj,MachAdj,machCoefAdj,machGridAdj
  REAL(KIND=REALTYPE) :: prefAdj, rhorefAdj,pInfCorrAdj
  REAL(KIND=REALTYPE) :: pinfdimAdj, rhoinfdimAdj
  REAL(KIND=REALTYPE) :: rhoinfAdj, pinfAdj
  REAL(KIND=REALTYPE) :: murefAdj, timerefAdj
  real(kind=realType), dimension(3) ::rotRateAdj,rotCenterAdj
  real(kind=realType), dimension(3) ::rotPointAdj,pointRefAdj
  real(kind=realType), dimension(4) :: time
  real(kind=realType)               :: timeAdj, timeOri

  integer :: ierr
  real(kind=realType) :: differ,relerr

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
  if(myID == 0) time(1) = mpi_wtime()


  if(myID == 0) then
     write(*,*) "Running verifyRAdj..."
     write(*,*) "proc sps block  i   j   k  sum(dwAdj)", &
          "    sum(dw)  rel Err differ"
  endif

  groundLevel = 1
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

  ! Get the final time for original routines


  call mpi_barrier(SUmb_comm_world, ierr)
  time(2) = mpi_wtime()
  if(myID == 0) then
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
  timings(:) = 0.0
  time(3) = mpi_wtime()

  ! Loop over the number of local blocks.
  spectralLoop: do sps=1,nTimeIntervalsSpectral
     domainResidualLoop: do nn=1,nDom

        ! Set some pointers to make the code more readable.
        call setPointers(nn,level,sps)

        !===============================================================
        ! Compute the residual for each cell.
        do kCell=2,kl
           do jCell=2,jl
              do iCell=2,il

                 ! Transfer the state w to the auxiliar array wAdj
                 ! and the coordinates x to the auxiliar array xAdj.

                 call copyADjointStencil(wAdj, xAdj,xBlockCornerAdj,alphaAdj,&
                      betaAdj,MachAdj,machCoefAdj,machGridAdj,iCell, jCell, kCell,&
                      nn,level,sps,pointRefAdj,rotPointAdj,&
                      prefAdj,rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                      rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
                      murefAdj, timerefAdj,pInfCorrAdj,liftIndex)

                 ! Compute the total residual.
                 ! This includes inviscid and viscous fluxes, artificial
                 ! dissipation, and boundary conditions.                   

                 call computeRAdjoint(wAdj,xAdj,xBlockCornerAdj,dwAdj,alphaAdj,&
                          betaAdj,MachAdj, &
                          MachCoefAdj,machGridAdj,iCell, jCell,  kCell, &
                          nn,level,sps, correctForK,secondHalo,prefAdj,&
                          rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                          rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
                          pointRefAdj,rotPointAdj,&
                          murefAdj, timerefAdj,pInfCorrAdj,liftIndex)


                  differ = (sum(dwAdj(:,sps))-sum(dw(iCell,jCell,kCell,:)))
                  relerr = differ/(.5*sum(dw(iCell,jCell,kCell,:)) + .5*sum(dwadj(:,sps)))
                  if( abs(differ) > 1e-10 .and. abs(relerr) > 1e-8) &
                       write(*,10) myID,sps, nn, iCell, jCell, kCell,               &
                       sum(dwAdj(:,sps)), sum(dw(iCell,jCell,kCell,:)), &
                       relerr,differ
                 !,(dwAdj(2,sps)), (dw(iCell,jCell,kCell,2)),(dwAdj(1,sps))-(dw(iCell,jCell,kCell,1)),(dwAdj(2,sps))-(dw(iCell,jCell,kCell,2)),(dwAdj(3,sps))-(dw(iCell,jCell,kCell,3)),(dwAdj(4,sps))-(dw(iCell,jCell,kCell,4)),(dwAdj(5,sps))-(dw(iCell,jCell,kCell,5))
                 
!                  ! Store difference to output to volume solution file.
!                  ! (resrho, resmom, resrhoe) have to be added to the volume
!                  ! output variables in the parameter file.

                 !dw(iCell,jCell,kCell,:) = dw(iCell,jCell,kCell,:) - dwAdj(:,sps)

              enddo
           enddo
        enddo

     enddo domainResidualLoop
  enddo spectralLoop
  ! Get new time and compute the elapsed ResidualAdj time.
    

  call mpi_barrier(SUmb_comm_world, ierr)
  time(4) = mpi_wtime()
  if(myID == 0) then
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
     do icell=1,10
        print *,'Time:',icell,timings(icell)
     end do
  endif

  ! Output formats.

10 format(1x,i3,2x,i3,2x,i3,2x,i3,1x,i3,1x,i3,2x,e10.3,1x,e10.3,1x,&
        e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3)

end subroutine verifyRAdj
