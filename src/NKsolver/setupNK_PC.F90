subroutine setupNK_PC(dRdwPre)

  !     ******************************************************************
  !     *                                                                *
  !     * Compute the dRdWPre matrix for the NK Solver                   *
  !     * This uses the reverse AD routines. The FD routines are MUCH    *
  !     * faster and SHOULD be used unless the problems warrants using   *
  !     * AD for the preconditioner.                                     *
  !     ******************************************************************
  !
  use ADjointVars
  use blockPointers       ! i/j/kl/b/e, i/j/k/Min/MaxBoundaryStencil
  use communication       ! procHalo(currentLevel)%nProcSend
  use inputDiscretization ! spaceDiscr
  USE inputTimeSpectral   ! nTimeIntervalsSpectral
  use iteration           ! overset, currentLevel
  use flowVarRefState     ! nw
  use inputTimeSpectral   ! spaceDiscr
  use inputDiscretization ! sigma,lumpedDiss
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  Mat dRdwPre
  Vec aVec,bVec
  !
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
  real(kind=realType)               ::setupTime,trace

  logical :: correctForK, secondHalo, exchangeTurb

  ! dR/dw stencil

  real(kind=realType), dimension(nw,nw,nTimeIntervalsSpectral) :: Aad, Bad, Cad, Dad, Ead, Fad, Gad
  real(kind=realType) ::eye(nw,nw)
  ! idxmgb - global block row index
  ! idxngb - global block column index

  integer(kind=intType) :: idxmgb, idxngb,ierr, sps, sps2,ilow,ihigh,i

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

  !store the current values of vis2,vis4 and reset vis2 for preconditioner
  !method based on (Hicken and Zingg,2008) AIAA journal,vol46,no.11
  vis2_ref = vis2
  vis4_ref = vis4
  lumpedDiss=.True.
  
  call MatZeroEntries(dRdwPre,ierr)
  call EChk(ierr,__FILE__,__LINE__)

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
                    call MatSetValuesBlocked(dRdWPre, 1, idxmgb, 1, idxngb,Aad(:,:,sps2), &
                         ADD_VALUES,ierr)
                    call EChk(ierr,__FILE__,__LINE__)
                 enddo

                 ! >>> west block B < W(i-1,j,k)
                 if( (iCell-1) >= 0 ) then
                    idxngb = globalCell(iCell-1,jCell,kCell)
                    if (idxngb >=0 .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdWPre, 1, idxmgb, 1, idxngb, Bad(:,:,sps), &
                            ADD_VALUES,ierr)
                       call EChk(ierr,__FILE__,__LINE__)
                    endif
                 endif

                 ! >>> east block C < W(i+1,j,k)
                 if( (iCell+1) <= ib ) then
                    idxngb = globalCell(iCell+1,jCell,kCell)
                    if (idxngb<nCellsGlobal*nTimeIntervalsSpectral .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdWPre, 1, idxmgb, 1, idxngb, Cad(:,:,sps), &
                            ADD_VALUES,ierr)
                       call EChk(ierr,__FILE__,__LINE__)
                    endif
                 end if

                 ! >>> south block D < W(i,j-1,k)
                 if( (jCell-1) >= 0 ) then
                    idxngb = globalCell(iCell,jCell-1,kCell)
                    if (idxngb>=0 .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdWPre, 1, idxmgb, 1, idxngb, Dad(:,:,sps), &
                            ADD_VALUES,ierr)
                       call EChk(ierr,__FILE__,__LINE__)
                    endif
                 endif

                 ! >>> north block E < W(i,j+1,k)
                 if( (jCell+1) <= jb ) then
                    idxngb = globalCell(iCell,jCell+1,kCell)
                    if (idxngb<nCellsGlobal*nTimeIntervalsSpectral .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdWPre, 1, idxmgb, 1, idxngb, Ead(:,:,sps), &
                            ADD_VALUES,ierr)
                       call EChk(ierr,__FILE__,__LINE__)
                    endif
                 end if

                 ! >>> back block F < W(i,j,k-1)
                 if( (kCell-1) >= 0 ) then
                    idxngb = globalCell(iCell,jCell,kCell-1)
                    if (idxngb>=0 .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdWPre, 1, idxmgb, 1, idxngb, Fad(:,:,sps), &
                            ADD_VALUES,ierr)
                       call EChk(ierr,__FILE__,__LINE__)
                    endif
                 endif

                 ! >>> front block G < W(i,j,k+1)
                 if( (kCell+1) <= kb ) then
                    idxngb = globalCell(iCell,jCell,kCell+1)
                    if (idxngb<nCellsGlobal*nTimeIntervalsSpectral .and. idxngb.ne.-5) then
                       call MatSetValuesBlocked(dRdWPre, 1, idxmgb, 1, idxngb, Gad(:,:,sps), &
                            ADD_VALUES,ierr)
                       call EChk(ierr,__FILE__,__LINE__)
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
  lumpedDiss=.False.
  call MatAssemblyBegin(dRdWPre,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call MatAssemblyEnd  (dRdWPre,MAT_FINAL_ASSEMBLY,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#ifdef USE_PETSC_3
  call MatSetOption(dRdWPre,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#else
  call MatSetOption(dRdWPre,MAT_NO_NEW_NONZERO_LOCATIONS,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif

  time(2) = mpi_wtime()
  call mpi_reduce(time(2)-time(1),setupTime,1,sumb_real,mpi_max,0,&
       SUmb_comm_world, ierr)
  if (myid == 0) then
     print *,'Assembly time:',setupTime
  end if

end subroutine setupNK_PC





