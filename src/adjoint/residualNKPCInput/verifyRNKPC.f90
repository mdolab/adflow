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
subroutine verifyRNKPC
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
  !     Local variables.

  integer(kind=intType) :: iCell, jCell, kCell, nn, level, sps
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw, nTimeIntervalsSpectral) :: wAdj
  real(kind=realType), dimension(-3:2,-3:2,-3:2,3,&
       nTimeIntervalsSpectral) :: siAdj, sjAdj, skAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,&
       nTimeIntervalsSpectral) ::sFaceIAdj,sFaceJAdj,sFaceKAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,3,&
       nTimeIntervalsSpectral) :: sAdj
  real(kind=realType),dimension(nTimeIntervalsSpectral) :: volAdj  
  real(kind=realType),dimension(3) :: rotRateAdj
  real(kind=realType), dimension(nw,nTimeIntervalsSpectral)  :: dwAdj

  logical :: secondHalo,exchangeTurb,correctforK,ierr
  real(kind=realType) :: time(3),timeori,timeadj,differ,relerr
  
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !
  ! Get the initial time.
  time(1) = mpi_wtime()

  if(myID == 0) then
     write(*,*) "Running verifyRNKPC..."
     write(*,*) "proc sps block  i   j   k  sum(dwAdj)", &
          "    sum(dw)  rel Err differ"
  endif
  level= 1
  groundLevel = 1
  ! Set the grid level of the current MG cycle, the value of the
  ! discretization and the logical fineGrid.
 
  currentLevel = level

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

  time(2) = mpi_wtime()

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Compute the residual using the node-by-node based routine.     *
  !     *                                                                *
  !     ******************************************************************
  !
  ! Get the initial ResidualAdj time.
  
  timings(:) = 0.0
  level = 1
  ! Loop over the number of local blocks.
  spectralLoop: do sps=1,nTimeIntervalsSpectral
     domainResidualLoop: do nn=1,nDom

        call setPointersAdj(nn,level,sps)

        !===============================================================
        ! Compute the residual for each cell.
        do kCell=2,kl
           do jCell=2,jl
              do iCell=2,il
                 ! Copy the required stencil
                
                 call copyNKPCStencil(iCell, jCell, kCell, nn, level, sps, wAdj, &
                      siAdj, sjAdj, skAdj, sAdj, sfaceIAdj, sfaceJAdj, sfaceKAdj, rotRateAdj,&
                      voladj)
              
                 ! Compute the total residual.
                 ! This includes inviscid and viscous fluxes, artificial
                 ! dissipation, and boundary conditions.                   
                 call computeRNKPC(wAdj, dwAdj, siAdj, sjAdj, skAdj, sAdj, volAdj, &
                      sfaceIAdj, sfaceJAdj, sfaceKadj, rotRateAdj, iCell, jCell, kCell, &
                      nn, level, sps)
                 
                 differ = (sum(dwAdj(:,sps))-sum(dw(iCell,jCell,kCell,:)/vol(icell,jcell,kcell)))
                 relerr = differ/(.5*sum(dw(iCell,jCell,kCell,:)) + .5*sum(dwadj(:,sps)))

                 if( abs(differ) > 1e-12 .and. abs(relerr) > 1e-10) &
                      write(*,10) myID,sps, nn, iCell, jCell, kCell,               &
                      sum(dwAdj(:,sps)), sum(dw(iCell,jCell,kCell,:)), &
                      relerr,differ
              enddo
           enddo
        enddo

     enddo domainResidualLoop
  enddo spectralLoop
  ! Get new time and compute the elapsed ResidualAdj time.
    
  time(3) = mpi_wtime()

  call mpi_reduce(time(2)-time(1),timeOri,1,sumb_real,mpi_max,0,&
       SUmb_comm_world, ierr)

  call mpi_reduce(time(3)-time(2),timeAdj,1,sumb_real,mpi_max,0,&
       SUmb_comm_world, ierr)

  ! Output elapsed time for the adjoint and FD computations.

  if( myID==0 ) then
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

10 format(1x,i3,2x,i3,2x,i3,2x,i3,1x,i3,1x,i3,2x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3)

end subroutine verifyRNKPC
