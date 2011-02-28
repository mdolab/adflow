!
!     ******************************************************************
!     *                                                                *
!     * File:          setupADjointPCMatrixTranspose.F90               *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 05-07-2010                                      *
!     * Last modified: 05-14-2010                                      *
!     *                                                                *
!     ******************************************************************
subroutine setupADjointPCMatrixTranspose

  use ADjointVars
  use ADjointPetsc
  use blockPointers       ! i/j/kl/b/e, i/j/k/Min/MaxBoundaryStencil
  use communication       ! procHalo(currentLevel)%nProcSend
  use inputDiscretization ! spaceDiscr
  USE inputTimeSpectral   ! nTimeIntervalsSpectral
  use iteration           ! overset, currentLevel
  use flowVarRefState     ! nw
  use inputTimeSpectral   ! spaceDiscr
  use inputADjoint        !sigma
  implicit none

  !     Local variables.

  integer(kind=intType) :: iCell, jCell, kCell, nn, level, m
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw, nTimeIntervalsSpectral) :: wAdj,wadjb
  real(kind=realType), dimension(-3:2,-3:2,-3:2,3,&
       nTimeIntervalsSpectral) :: siAdj, sjAdj, skAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,&
       nTimeIntervalsSpectral) ::sFaceIAdj,sFaceJAdj,sFaceKAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,3,&
       nTimeIntervalsSpectral) :: sAdj
  real(kind=realType),dimension(nTimeIntervalsSpectral) :: volAdj  
  real(kind=realType),dimension(3) :: rotRateAdj
  real(kind=realType), dimension(nw,nTimeIntervalsSpectral)  :: dwAdj,dwadjb
  real(kind=realType), dimension(2) :: time
  real(kind=realType)               ::setupTime

  logical :: correctForK, secondHalo, exchangeTurb

  ! dR/dw stencil

  real(kind=realType), dimension(nw,nw,nTimeIntervalsSpectral) :: Aad, Bad, Cad, Dad, Ead, Fad, Gad
  real(kind=realType) ::eye(nw,nw)
  ! idxmgb - global block row index
  ! idxngb - global block column index

  integer(kind=intType) :: idxmgb, idxngb,ierr, sps, sps2,ilow,ihigh

  !Reference values of the dissipation coeff for the preconditioner
  real(kind=realType) :: vis2_ref, vis4_ref

  ! Set the grid level of the current MG cycle, the value of the
  ! discretization and the logical correctForK.
  level = 1_intType
  currentLevel = level
  time(1) = mpi_wtime()
  rkStage = 0
  currentLevel = groundLevel

  !     ******************************************************************
  !     *                                                                *
  !     * Compute the ADjoint matrix dR/dW using Tapenade's reverse mode *
  !     * of Automatic Differentiation.  NOTE: This is the reason I have *
  !     * been writing the word "ADjoint" with A and D capitalized. A    *
  !     * simple play with letter so that:                               *
  !     *                                                                *
  !     * ADjoint = Automatically Differentiated adjoint                 *
  !     *                                                                *
  !     ******************************************************************
  !
  ! Send some feedback to screen.

  if (myid == 0) then
     print * ,"Assembling NK KSP PC matrix..."
  end if
 
  !store the current values of vis2,vis4 and reset vis2 for preconditioner
  !method based on (Hicken and Zingg,2008) AIAA journal,vol46,no.11
  vis2_ref = vis2
  vis4_ref = vis4
  lumpedDiss=.True.
  
  call MatZeroEntries(dRdwPreT,ierr)
  call EChk(ierr,__file__,__line__)

  domainLoopAD: do nn=1,nDom

     ! Loop over the number of time instances for this block.
     spectralLoop: do sps=1,nTimeIntervalsSpectral
        call setPointersAdj(nn,level,sps)
        ! Loop over location of output (R) cell of residual
        do kCell = 2, kl
           do jCell = 2, jl
              do iCell = 2, il
                 ! Copy the state w to the wAdj array in the stencil

                 call copyNKPCStencil(iCell, jCell, kCell, nn, level, sps, wAdj, &
                      siAdj, sjAdj, skAdj, sAdj, sfaceIAdj, sfaceJAdj, sfaceKAdj, rotRateAdj,&
                      voladj)
                 Aad(:,:,:)  = zero
                 Bad(:,:,:)  = zero
                 Cad(:,:,:)  = zero
                 Dad(:,:,:)  = zero
                 Ead(:,:,:)  = zero
                 Fad(:,:,:)  = zero
                 Gad(:,:,:)  = zero

                 mLoop: do m = 1, nw      ! Loop over output cell residuals (R)
                    ! Initialize the seed for the reverse mode
                    dwAdjb(:,:) = 0.
                    dwAdjb(m,sps) = 1.
                    dwAdj(:,:)  = 0.
                    wAdjb(:,:,:,:,:)  = 0.  !dR(m)/dw

                    ! Call the reverse mode of residual computation.
                    !
                    !                          dR(iCell,jCell,kCell,l)
                    ! wAdjb(ii,jj,kk,n) = --------------------------------
                    !                     dW(iCell+ii,jCell+jj,kCell+kk,n)

                    ! Call reverse mode of residual computation

                    call COMPUTERNKPC_B(wadj, wadjb, dwadj, dwadjb, siadj, sjadj, &
                         &  skadj, sadj, voladj, sfaceiadj, sfacejadj, sfacekadj, rotrateadj, &
                         &  icell, jcell, kcell, nn, level, sps)

                    ! Store the block Jacobians (by rows).

                    Aad(m,:,:)  = wAdjB( 0, 0, 0,:,:)
                    Bad(m,:,:)  = wAdjB(-1, 0, 0,:,:)
                    Cad(m,:,:)  = wAdjB( 1, 0, 0,:,:)
                    Dad(m,:,:)  = wAdjB( 0,-1, 0,:,:)
                    Ead(m,:,:)  = wAdjB( 0, 1, 0,:,:)
                    Fad(m,:,:)  = wAdjB( 0, 0,-1,:,:)
                    Gad(m,:,:)  = wAdjB( 0, 0, 1,:,:)
                 enddo mLoop
                 ! Global matrix block row mgb function of node indices.

                 idxmgb = globalCell(iCell,jCell,kCell)

                 ! >>> center block A < W(i,j,k)
                 do sps2 = 1,nTimeIntervalsSpectral
                    idxngb = flowDoms(nn,level,sps2)%globalCell(iCell,jCell,kCell)
                    call MatSetValuesBlocked(dRdwPreT, 1, idxngb, 1, idxmgb,&
                         transpose(Aad(:,:,sps2)),ADD_VALUES,ierr)
                    call EChk(ierr,__file__,__line__)
                 enddo

                 ! >>> west block B < W(i-1,j,k)
                 if( (iCell-1) >= 0 ) then
                    idxngb = globalCell(iCell-1,jCell,kCell)
                    if (idxngb >=0 .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdwPreT,1,idxngb,1,idxmgb,&
                            transpose( Bad(:,:,sps)),ADD_VALUES,ierr)
                       call EChk(ierr,__file__,__line__)
                    endif
                 endif

                 ! >>> east block C < W(i+1,j,k)
                 if( (iCell+1) <= ib ) then
                    idxngb = globalCell(iCell+1,jCell,kCell)
                    if (idxngb<nCellsGlobal*nTimeIntervalsSpectral .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdwPreT, 1, idxngb, 1, idxmgb,&
                            transpose(Cad(:,:,sps)),ADD_VALUES,ierr)
                       call EChk(ierr,__file__,__line__)
                    endif
                 end if

                 ! >>> south block D < W(i,j-1,k)
                 if( (jCell-1) >= 0 ) then
                    idxngb = globalCell(iCell,jCell-1,kCell)
                    if (idxngb>=0 .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdwPreT, 1, idxngb, 1, idxmgb,&
                            transpose( Dad(:,:,sps)),ADD_VALUES,ierr)
                       call EChk(ierr,__file__,__line__)
                    endif
                 endif

                 ! >>> north block E < W(i,j+1,k)
                 if( (jCell+1) <= jb ) then
                    idxngb = globalCell(iCell,jCell+1,kCell)
                    if (idxngb<nCellsGlobal*nTimeIntervalsSpectral .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdwPreT, 1, idxngb, 1, idxmgb,&
                            transpose( Ead(:,:,sps)),ADD_VALUES,ierr)
                       call EChk(ierr,__file__,__line__)
                    endif
                 end if

                 ! >>> back block F < W(i,j,k-1)
                 if( (kCell-1) >= 0 ) then
                    idxngb = globalCell(iCell,jCell,kCell-1)
                    if (idxngb>=0 .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdwPreT, 1, idxngb, 1, idxmgb,&
                            transpose( Fad(:,:,sps)),ADD_VALUES,ierr)
                       call EChk(ierr,__file__,__line__)
                    endif
                 endif

                 ! >>> front block G < W(i,j,k+1)
                 if( (kCell+1) <= kb ) then
                    idxngb = globalCell(iCell,jCell,kCell+1)
                    if (idxngb<nCellsGlobal*nTimeIntervalsSpectral .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdwPreT, 1, idxngb, 1, idxmgb,&
                            transpose( Gad(:,:,sps)),ADD_VALUES,ierr)
                       call EChk(ierr,__file__,__line__)
                    endif
                 end if
              enddo
           enddo
        enddo
     enddo spectralLoop
  enddo domainLoopad

  !Return dissipation Parameters to normal
  vis2 = vis2_ref
  vis4 = vis4_ref

  call MatAssemblyBegin(dRdwPreT,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__file__,__line__)
  call MatAssemblyEnd  (dRdwPreT,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__file__,__line__)
#ifdef USE_PETSC_3
  call MatSetOption(dRdwPreT,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE,ierr)
  call EChk(ierr,__file__,__line__)
#else
  call MatSetOption(dRdwPreT,MAT_NO_NEW_NONZERO_LOCATIONS,ierr)
  call EChk(ierr,__file__,__line__)
#endif

  time(2) = mpi_wtime()
  call mpi_reduce(time(2)-time(1),setupTime,1,sumb_real,mpi_max,0,&
       SUmb_comm_world, ierr)

  if (myid == 0) then
     print *,'Done PC Assembly'
     print *,'Time:',setupTime
  end if
end subroutine setupADjointPCMatrixTranspose

















! subroutine setupADjointPCMatrixTranspose
!   !
!   !     ******************************************************************
!   !     *                                                                *
!   !     * Compute the transpose of the matrix dRdWPre of the discrete    *
!   !     * ADjoint problem for                                            *
!   !     * subsequent use by the solveADjointPETSc subroutine. The entries*
!   !     * in dRdwPre are obtained using the automatically differentiated *
!   !     * routines generated by Tapenade, with a modification to the     *
!   !     * dissipation parameters to improve the conditioning of the      *
!   !     * matrix.                                                        *
!   !     *                                                                *
!   !     * The ordering of the unknowns in the ADjoint matrix used here   *
!   !     * is based on the global node numbering and is consistent with   *
!   !     * the ordering used in the vector for the ADjoint problem        *
!   !     * assembled in setupADjointRHS.                                  *
!   !     *                                                                *
!   !     ******************************************************************
!   !
!   use ADjointPETSc
!   use ADjointVars ! nCellsGlobal, nCellsLocal, nOffsetLocal
!   use blockPointers       ! i/j/kl/b/e, i/j/k/Min/MaxBoundaryStencil
!   use cgnsGrid            ! cgnsDoms
!   use communication       ! procHalo(currentLevel)%nProcSend
!   use inputDiscretization ! spaceDiscr
!   USE inputTimeSpectral   ! nTimeIntervalsSpectral
!   use iteration           ! overset, currentLevel
!   use flowVarRefState     ! nw
!   use inputTimeSpectral   ! spaceDiscr
!   use inputADjoint        !sigma
!   implicit none
!   !
!   !     Local variables.
!   !
!   integer(kind=intType) :: discr, nHalo,level
!   integer(kind=intType) :: iCell, jCell, kCell
!   integer(kind=intType) :: mm, nn, m, n,idxstate,idxres
!   integer(kind=intType) :: ii, jj, kk, i, j, k,liftIndex,l

!   logical :: fineGrid, correctForK, exchangeTurb,secondhalo

!   real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral) :: wAdj, wAdjB
!   real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral)  :: xAdj, xAdjB

!   real(kind=realType), dimension(nw,nTimeIntervalsSpectral) :: dwAdj, dwAdjB

!   REAL(KIND=REALTYPE) :: machadj, machcoefadj, uinfadj, pinfcorradj
!   REAL(KIND=REALTYPE) :: machadjb, machcoefadjb,machgridadj, machgridadjb
!   REAL(KIND=REALTYPE) :: prefadj, rhorefadj
!   REAL(KIND=REALTYPE) :: pinfdimadj, rhoinfdimadj
!   REAL(KIND=REALTYPE) :: rhoinfadj, pinfadj
!   REAL(KIND=REALTYPE) :: murefadj, timerefadj
!   REAL(KIND=REALTYPE) :: alphaadj, betaadj
!   REAL(KIND=REALTYPE) :: alphaadjb, betaadjb
!   REAL(KIND=REALTYPE), DIMENSION(3) :: rotcenteradj
!   REAL(KIND=REALTYPE), DIMENSION(3) :: rotrateadj 
!   REAL(KIND=REALTYPE) :: rotcenteradjb(3), rotrateadjb(3)
!   REAL(KIND=REALTYPE) :: pointrefadj(3), pointrefadjb(3), rotpointadj(3)&
!        &  , rotpointadjb(3)
!   REAL(KIND=REALTYPE) :: xblockcorneradj(2, 2, 2, 3,nTimeIntervalsSpectral), xblockcorneradjb(2&
!        &  , 2, 2, 3,nTimeIntervalsSpectral)

!   integer(kind=intType), dimension(0:nProc-1) :: offsetRecv

!   real(kind=realType), dimension(2) :: time
!   real(kind=realType)               :: timeAdjLocal, timeAdj
!   real(kind=realType) :: prefadjb,rhorefadjb, pinfdimadjb,rhoinfdimadjb,&
!        rhoinfadjb,PINFADJB,MUREFADJB,PINFCORRADJB
!   ! dR/dw stencil

!   real(kind=realType), dimension(nw,nw,nTimeIntervalsSpectral) :: Aad, Bad, BBad, &
!        Cad, CCad, Dad, DDad, Ead, EEad, Fad, FFad, Gad, GGad

!   ! idxmgb - global block row index
!   ! idxngb - global block column index

!   integer(kind=intType) :: idxmgb, idxngb,ierr, sps,sps2

!   ! idxmgb - array of global row indices
!   ! idxngb - array of global column indices

!   integer(kind=intType), dimension(nw) :: idxmg, idxng

!   !Reference values of the dissipation coeff for the preconditioner
!   real(kind=realType) :: vis2_ref, vis4_ref

!   !
!   !     ******************************************************************
!   !     *                                                                *
!   !     * Begin execution.                                               *
!   !     *                                                                *
!   !     ******************************************************************
!   !
! #ifndef USE_NO_PETSC

!   ! Set the grid level of the current MG cycle, the value of the
!   ! discretization and the logical correctForK.
!   level = 1_intType
!   currentLevel = level
!   !discr        = spaceDiscr
!   fineGrid     = .true.

!   ! Determine whether or not the total energy must be corrected
!   ! for the presence of the turbulent kinetic energy and whether
!   ! or not turbulence variables should be exchanged.

!   if( kPresent ) then
!      if((currentLevel <= groundLevel) .or. turbCoupled) then
!         correctForK = .true.
!      else
!         correctForK = .false.
!      endif
!   else
!      correctForK = .false.
!   endif

!   exchangeTurb = .false.

!   ! Set the value of secondHalo, depending on the situation.
!   ! In the full MG (currentLevel < groundLevel) the second halo is
!   ! always set; otherwise only on the finest mesh in the current mg
!   ! cycle.

!   if(currentLevel <= groundLevel) then
!      secondHalo = .true.
!   else
!      secondHalo = .false.
!   endif

!   !
!   !     ******************************************************************
!   !     *                                                                *
!   !     * Exchange halo data to make sure it is up-to-date.              *
!   !     * (originally called inside "rungeKuttaSmoother" subroutine).    *
!   !     *                                                                *
!   !     ******************************************************************
!   !
!   ! Exchange the pressure if the pressure must be exchanged early.
!   ! Only the first halo's are needed, thus whalo1 is called.
!   ! Only on the fine grid.

!   if(exchangePressureEarly .and. currentLevel <= groundLevel) &
!        call whalo1(currentLevel, 1_intType, 0_intType, .true.,&
!        .false., .false.)

!   ! Apply all boundary conditions to all blocks on this level.

!   call applyAllBC(secondHalo)

!   ! Exchange the solution. Either whalo1 or whalo2
!   ! must be called.

!   if( secondHalo ) then
!      call whalo2(currentLevel, 1_intType, nMGVar, .true., &
!           .true., .true.)
!   else
!      call whalo1(currentLevel, 1_intType, nMGVar, .true., &
!           .true., .true.)
!   endif

!   ! Reset the values of rkStage and currentLevel, such that
!   ! they correspond to a new iteration.

!   rkStage = 0
!   currentLevel = groundLevel

!   !
!   !     ******************************************************************
!   !     *                                                                *
!   !     * Compute the ADjoint matrix dR/dW using Tapenade's reverse mode *
!   !     * of Automatic Differentiation.  NOTE: This is the reason I have *
!   !     * been writing the word "ADjoint" with A and D capitalized. A    *
!   !     * simple play with letter so that:                               *
!   !     *                                                                *
!   !     * ADjoint = Automatically Differentiated adjoint                 *
!   !     *                                                                *
!   !     ******************************************************************
!   !
!   ! Send some feedback to screen.

!   if( myid ==0 ) &
!        write(*,10) "Assembling ADjoint PC Transpose matrix..."

!   call mpi_barrier(SUmb_comm_world, ierr)
!   call EChk(PETScIerr,__file__,__line__)
!   call cpu_time(time(1))

!   !zero the matrix for dRdWPre Insert call
!   call MatZeroEntries(dRdwPreT,PETScIerr)
!   call EChk(PETScIerr,__file__,__line__)

!   !store the current values of vis2,vis4 and reset vis2 for preconditioner
!   !method based on (Hicken and Zingg,2008) AIAA journal,vol46,no.11
!   vis2_ref = vis2
!   vis4_ref = vis4
!   !sigma = 4 !setup as an input parameter....
!   !vis2 = vis2_ref+sigma*vis4_ref
!   !vis4 = 0.0
!   lumpedDiss=.True.


!   domainLoopAD: do nn=1,nDom

!      ! Loop over the number of time instances for this block.
!      spectralLoop: do sps=1,nTimeIntervalsSpectral
!         call setPointersAdj(nn,level,sps)

!         ! Loop over location of output (R) cell of residual
!         do kCell = 2, kl
!            do jCell = 2, jl
!               do iCell = 2, il
!                  ! Copy the state w to the wAdj array in the stencil
               
!                  call copyADjointStencil(wAdj, xAdj,xBlockCornerAdj,alphaAdj,&
!                       betaAdj,MachAdj,machCoefAdj,machGridAdj,iCell, jCell, kCell,&
!                       nn,level,sps,pointRefAdj,rotPointAdj,&
!                       prefAdj,rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
!                       rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
!                       murefAdj, timerefAdj,pInfCorrAdj,liftIndex)
             
!                  Aad(:,:,:)  = zero
!                  Bad(:,:,:)  = zero
!                  !BBad(:,:,:) = zero
!                  Cad(:,:,:)  = zero
!                  !CCad(:,:,:) = zero
!                  Dad(:,:,:)  = zero
!                  !DDad(:,:,:) = zero
!                  Ead(:,:,:)  = zero
!                  !EEad(:,:,:) = zero
!                  Fad(:,:,:)  = zero
!                  !FFad(:,:,:) = zero
!                  Gad(:,:,:)  = zero
!                  !GGad(:,:,:) = zero

!                  mLoop: do m = 1, nw      ! Loop over output cell residuals (R)
!                     ! Initialize the seed for the reverse mode
!                     dwAdjb(:,:) = 0.
!                     dwAdjb(m,sps) = 1.
!                     dwAdj(:,:)  = 0.
!                     wAdjb(:,:,:,:,:)  = 0.  !dR(m)/dw
!                     xadjb = 0.
!                     alphaadjb = 0.
!                     betaadjb = 0.
!                     machadjb = 0.
!                     rotrateadjb(:)=0.

!                     ! Call the reverse mode of residual computation.
!                     !
!                     !                          dR(iCell,jCell,kCell,l)
!                     ! wAdjb(ii,jj,kk,n) = --------------------------------
!                     !                     dW(iCell+ii,jCell+jj,kCell+kk,n)

!                     ! Call reverse mode of residual computation

!                     call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
!                          &  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
!                          &  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
!                          &  icell, jcell, kcell, nn, level, sps, correctfork, secondhalo, prefadj&
!                          &  , rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj&
!                          &  , rotrateadjb, rotcenteradj, rotcenteradjb, pointrefadj, pointrefadjb&
!                          &  , rotpointadj, rotpointadjb, murefadj, timerefadj, pinfcorradj, &
!                          &  liftindex)


!                     ! Store the block Jacobians (by rows).

!                     Aad(m,:,:)  = wAdjB( 0, 0, 0,:,:)
!                     Bad(m,:,:)  = wAdjB(-1, 0, 0,:,:)
!                     !BBad(m,:,:) = wAdjB(-2, 0, 0,:,:)
!                     Cad(m,:,:)  = wAdjB( 1, 0, 0,:,:)
!                     !CCad(m,:,:) = wAdjB( 2, 0, 0,:,:)
!                     Dad(m,:,:)  = wAdjB( 0,-1, 0,:,:)
!                     !DDad(m,:,:) = wAdjB( 0,-2, 0,:,:)
!                     Ead(m,:,:)  = wAdjB( 0, 1, 0,:,:)
!                     !EEad(m,:,:) = wAdjB( 0, 2, 0,:,:)
!                     Fad(m,:,:)  = wAdjB( 0, 0,-1,:,:)
!                     !FFad(m,:,:) = wAdjB( 0, 0,-2,:,:)
!                     Gad(m,:,:)  = wAdjB( 0, 0, 1,:,:)
!                     !GGad(m,:,:) = wAdjB( 0, 0, 2,:,:)



!                  enddo mLoop


!                  if(PETScBlockMatrix) then

!                     ! Global matrix block row mgb function of node indices.
!                     !
!                     ! MatSetValuesBlocked() uses 0-based row and column 
!                     ! numbers but the global node numbering already accounts
!                     ! for that since it starts at node 0.

                    
!                     idxmgb = globalCell(iCell,jCell,kCell)
!                     ! >>> center block A < W(i,j,k)
!                     do sps2 = 1,nTimeIntervalsSpectral
!                        idxngb = flowDoms(nn,level,sps2)%globalCell(iCell,jCell,kCell)!idxmgb
!                        call MatSetValuesBlocked(dRdWPret, 1, idxngb, 1, idxmgb, &
!                             transpose(Aad(:,:,sps2)), ADD_VALUES,PETScIerr)
!                        call EChk(PETScIerr,__file__,__line__)
!                     enddo

!                     ! >>> west block B < W(i-1,j,k)
!                     if( (iCell-1) >= 0 ) then
!                        idxngb = globalCell(iCell-1,jCell,kCell)!idxmgb
!                        if (idxngb >=0 .and. idxngb.ne.-5) then
!                           call MatSetValuesBlocked(dRdWPret, 1, idxngb, 1, idxmgb, &
!                                transpose(Bad(:,:,sps)), ADD_VALUES,PETScIerr)
!                           call EChk(PETScIerr,__file__,__line__)
!                        endif
!                     endif

!                     ! >>> east block C < W(i+1,j,k)
!                     if( (iCell+1) <= ib ) then
!                        idxngb = globalCell(iCell+1,jCell,kCell)!idxmgb
!                        if (idxngb<nCellsGlobal*nTimeIntervalsSpectral .and. idxngb.ne.-5) then
!                           call MatSetValuesBlocked(dRdWPret, 1, idxngb, 1, idxmgb, &
!                                transpose(Cad(:,:,sps)), ADD_VALUES,PETScIerr)
!                           call EChk(PETScIerr,__file__,__line__)
!                        endif
!                     end if

!                     ! >>> south block D < W(i,j-1,k)
!                     if( (jCell-1) >= 0 ) then
!                        idxngb = globalCell(iCell,jCell-1,kCell)!idxmgb
!                        if (idxngb>=0 .and. idxngb.ne.-5) then
!                           call MatSetValuesBlocked(dRdWPret, 1, idxngb, 1, idxmgb, &
!                                transpose(Dad(:,:,sps)), ADD_VALUES,PETScIerr)
!                           call EChk(PETScIerr,__file__,__line__)

!                        endif
!                     endif

!                     ! >>> north block E < W(i,j+1,k)
!                     if( (jCell+1) <= jb ) then
!                        idxngb = globalCell(iCell,jCell+1,kCell)!idxmgb
!                        if (idxngb<nCellsGlobal*nTimeIntervalsSpectral .and. idxngb.ne.-5) then
!                           call MatSetValuesBlocked(dRdWPret, 1, idxngb, 1, idxmgb, &
!                                transpose(Ead(:,:,sps)), ADD_VALUES,PETScIerr)
!                           call EChk(PETScIerr,__file__,__line__)
!                        endif
!                     end if

!                     ! >>> back block F < W(i,j,k-1)

!                     if( (kCell-1) >= 0 ) then
!                        idxngb = globalCell(iCell,jCell,kCell-1)!idxmgb
!                        if (idxngb>=0 .and. idxngb.ne.-5) then
!                           call MatSetValuesBlocked(dRdWPret, 1, idxngb, 1, idxmgb, &
!                                transpose(Fad(:,:,sps)), ADD_VALUES,PETScIerr)
!                           call EChk(PETScIerr,__file__,__line__)
!                        endif
!                     endif

!                     ! >>> front block G < W(i,j,k+1)
!                     if( (kCell+1) <= kb ) then
!                        idxngb = globalCell(iCell,jCell,kCell+1)!idxmgb
!                        if (idxngb<nCellsGlobal*nTimeIntervalsSpectral .and. idxngb.ne.-5) then
!                           call MatSetValuesBlocked(dRdWPret, 1, idxngb, 1, idxmgb, &
!                                transpose(Gad(:,:,sps)), ADD_VALUES,PETScIerr)
!                           call EChk(PETScIerr,__file__,__line__)
!                        endif
!                     end if
!                  else ! PETScBlockMatrix
!                     print *,'Non-Block Matrix not setup currently'
!                     stop
!                  endif ! PETScBlockMatrix
!                  !end do

!               enddo
!            enddo
!         enddo
!      enddo spectralLoop
!      !===============================================================

!   enddo domainLoopad

!   !Return dissipation Parameters to normal
!   vis2 = vis2_ref
!   vis4 = vis4_ref

!   call MatAssemblyBegin(dRdWPreT,MAT_FINAL_ASSEMBLY,PETScIerr)
!   call EChk(PETScIerr,__file__,__line__)
!   call MatAssemblyEnd  (dRdWPreT,MAT_FINAL_ASSEMBLY,PETScIerr)
!   call EChk(PETScIerr,__file__,__line__)
! #ifdef USE_PETSC_3
!   call MatSetOption(dRdWPreT,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE,PETScIerr)
!   call EChk(PETScIerr,__file__,__line__)
! #else
!   call MatSetOption(dRdWPreT,MAT_NO_NEW_NONZERO_LOCATIONS,PETScIerr)
!   call EChk(PETScIerr,__file__,__line__)
! #endif
!   ! Get new time and compute the elapsed time.

!   call cpu_time(time(2))
!   timeAdjLocal = time(2)-time(1)

!   ! Determine maximum time using MPI reduce
!   ! with operation mpi_max.

!   call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
!        mpi_max, 0, SUMB_PETSC_COMM_WORLD, PETScIerr)
!   call EChk(PETScIerr,__file__,__line__)
!   if( myid==0 ) &
!        write(*,20) "Assembling ADjoint PCT matrix time (s) = ", timeAdj

!   ! Output formats.

! 10 format(a)
! 20 format(a,1x,f8.2)

!   !=================================================================

! contains

!   !===============================================================

!   subroutine blockIndices(indexBlock, indexArray)
!     !
!     !       ****************************************************************
!     !       *                                                              *
!     !       * blockIndices fill the array containing the block matrix      *
!     !       * indices used when assembling the Jacobian matrix with the    *
!     !       * PETSc function MatSetValues().                               *
!     !       *                                                              *
!     !       ****************************************************************
!     !
!     implicit none
!     !
!     !       Subroutine arguments.
!     !
!     integer(kind=intType), intent(in) :: indexBlock
!     integer(kind=intType), dimension(nw), intent(out) :: indexArray
!     !
!     !       Local variables.
!     !
!     integer(kind=intType) :: idx
!     !
!     !       ****************************************************************
!     !       *                                                              *
!     !       * Begin execution.                                             *
!     !       *                                                              *
!     !       ****************************************************************
!     !
!     indexArray(1) = indexBlock * nw
!     do idx = 2, nw
!        indexArray(idx) = indexArray(idx-1) + 1
!     enddo

!   end subroutine blockIndices

!   !===============================================================
  
! #endif

! end subroutine setupADjointPCMatrixTranspose
