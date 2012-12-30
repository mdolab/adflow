!
!     ******************************************************************
!     *                                                                *
!     * File:          verifydRdxfile.f90                              *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 06-12-2009                                      *
!     * Last modified: 06-12-2009                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine verifydRdxfile(level)
!
!     ******************************************************************
!     *                                                                *
!     *  This subroutine computes the values for dRdw and compares     *
!     *  the Tapenade result to the finite-difference results.         *
!     *                                                                *
!     ******************************************************************

#ifndef USE_NO_PETSC
      use blockPointers ! block (nDoms,flowDoms), globalCell
      use flowvarrefstate
      use communication
      use iteration     ! groundLevel
      use inputTimeSpectral ! spaceDiscr,nTimeIntervalsSpectral
      use inputIO

      !from old verify routine
      use ADjointPETSc, only: drdx,insert_values,petscierr,&
           mat_final_assembly,petsc_viewer_draw_world,&
           petsc_viewer_stdout_world,add_values
      !use FDPETSc, only: DRDWFD
      use precision
      !use blockPointers
      !use flowvarrefstate
      !use iteration
      use inputIteration
      !implicit none
     

      implicit none

!
!     Subroutine arguments
      integer(kind=intType), intent(in) :: level


!
!     Local variables 
!
      integer(kind=intType) :: i, j, k, n, assembled,mmmm,nnnn
      integer(kind=intType) :: iCell,jCell,kCell,liftIndex
      !real(kind=realType), dimension(nw) :: dwL2
      !real(kind=realType), dimension(nx, ny, nz, nw) :: dwerr
      real(kind=realType), dimension(10) :: time
      real(kind=realType) :: timeRes
!      real(kind=realType), dimension(4) :: time
      real(kind=realType) ::  timeOri
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral) :: wAdj
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral) :: wAdjb
      !real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wFD
!!$      real(kind=realType), dimension(-2:3,-2:3,-2:3,3)  :: xAdj
!!$      real(kind=realType), dimension(-2:3,-2:3,-2:3,3)  :: xAdjb
!!$      real(kind=realType), dimension(-2:3,-2:3,-2:3,3)  :: xFD
      real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral)  :: xAdj
      real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral)  :: xAdjb
      !real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral)  :: xFD


      REAL(KIND=REALTYPE) :: machadj, machcoefadj, pinfcorradj,machgridadj
      REAL(KIND=REALTYPE) :: machadjb, machcoefadjb, pinfcorradjb,machgridadjb
      REAL(KIND=REALTYPE) :: prefadj, rhorefadj
      REAL(KIND=REALTYPE) :: pinfdimadj, rhoinfdimadj
      REAL(KIND=REALTYPE) :: rhoinfadj, pinfadj
      REAL(KIND=REALTYPE) :: rhoinfadjb
      REAL(KIND=REALTYPE) :: murefadj, timerefadj
      REAL(KIND=REALTYPE) :: alphaadj, betaadj
      REAL(KIND=REALTYPE) :: alphaadjb, betaadjb

      character fileName*32, dataName*32
      real(kind=realType), dimension(nw,nTimeIntervalsSpectral) :: dwAdj,dwAdjb,dwAdjRef
      !real(kind=realType), dimension(nw,nTimeIntervalsSpectral) :: dwAdjP, dwAdjM
      real(kind=realType) :: deltaw, wAdjRef
      real(kind=realType) :: timeAdj, timeFD, timeResAdj

      !integer(kind=intType), dimension(nDom) :: maxglobalcell
      integer(kind=intType) :: idx, ii, jj, kk, idxnode, idxres, m, l,nnn
      integer(kind=intType) :: sps, nTime, max_nTime, nHalo, nn, discr,sps2
      !real(kind=realType), allocatable, dimension(:,:,:,:) :: dRdwAdj,dRdwFD1,dRdwFD2,dRdwErr

      REAL(KIND=REALTYPE), DIMENSION(3) :: rotcenteradj
      REAL(KIND=REALTYPE), DIMENSION(3):: rotrateadj
      REAL(KIND=REALTYPE) :: rotcenteradjb(3), rotrateadjb(3)
      REAL(KIND=REALTYPE) :: pointrefadj(3), pointrefadjb(3), rotpointadj(3)&
           &  , rotpointadjb(3)
      REAL(KIND=REALTYPE) :: xblockcorneradj(2, 2, 2, 3,nTimeIntervalsSpectral), xblockcorneradjb(2&
           &  , 2, 2, 3,nTimeIntervalsSpectral)

      integer :: ierr, testnode
      logical :: fineGrid, correctForK, exchangeTurb,secondHalo
      REAL(KIND=REALTYPE) ::value1
 

      integer :: unitdRdx = 8,ierror
      character(len = 20)::outfile,testfile
      write(testfile,100) myid!12
100   format (i5)  
      testfile=adjustl(testfile)
      write(outfile,101) trim(testfile)!testfile
101   format("FDdRdxfile",a,".out")
      unitdrdx = 8+myID

      
      open (UNIT=unitdRdx,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifydRdxFile", &
           "Something wrong when &
           &calling open")

      !FULL FD Variables
      
!!$      !real(kind=realType), dimension(nx, ny, nz, nw,nx) :: dRdwErr, dRdwErrRel
!!$      real(kind=realType), dimension(nw) :: dRdwL2, dRdwL2Rel
!!$      real(kind=realtype), dimension(0:ib,0:jb,0:kb,1:nw)::dwp,dwm,dwtemp
!!$      real(kind=realType), dimension(0:ib,0:jb,0:kb,1:nw) :: wtemp
!!$      real(kind=realType), dimension(0:ib,0:jb,0:kb) :: ptemp
!!$      integer(kind=intType) :: istate, jstate, kstate,ires
!!$      real(kind=realType), dimension(-2:2,-2:2,-2:2) :: pAdjtemp
!!$      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdjtemp
!!$      real(kind=realType), dimension(nx,ny,nz,nw) :: wFD2
!!$      !real(kind=realType), dimension(ib*jb*kb*nw,ib*jb*kb*nw) :: dRdw
!!$      !real(kind=realType), dimension(ib*jb*kb*nw,ib*jb*kb*nw) :: dRdwFD

!     ******************************************************************
!     *                                                                *
!     *  Begin execution                                               *
!     *                                                                *
!     ******************************************************************

!!$      if( myID==0 ) write(*,*) "Running verifydRdxFile..."
!!$
!!$      currentLevel = level
!!$      !discr        = spaceDiscr
!!$      fineGrid     = .true.
!!$
!!$      ! Determine whether or not the total energy must be corrected
!!$      ! for the presence of the turbulent kinetic energy.
!!$
!!$      if( kPresent ) then
!!$         if((currentLevel <= groundLevel) .or. turbCoupled) then
!!$            correctForK = .true.
!!$         else
!!$            correctForK = .false.
!!$         endif
!!$      else
!!$         correctForK = .false.
!!$      endif
!!$
!!$      ! and whether or not turbulence variables should be exchanged
!!$      exchangeTurb = .false.
!!$
!!$      
!!$      ! Set the value of secondHalo, depending on the situation.
!!$      ! In the full MG (currentLevel < groundLevel) the second halo is
!!$      ! always set; otherwise only on the finest mesh in the current mg
!!$      ! cycle.
!!$
!!$      if(currentLevel <= groundLevel) then
!!$         secondHalo = .true.
!!$      else
!!$         secondHalo = .false.
!!$      endif
!!$
!!$
!!$!
!!$!     ******************************************************************
!!$!     *                                                                *
!!$!     * Exchange halo data to make sure it is up-to-date.              *
!!$!     * (originally called inside "rungeKuttaSmoother" subroutine).    *
!!$!     *                                                                *
!!$!     ******************************************************************
!!$!
!!$      ! Exchange the pressure if the pressure must be exchanged early.
!!$      ! Only the first halo's are needed, thus whalo1 is called.
!!$      ! Only on the fine grid.
!!$      
!!$      if(exchangePressureEarly .and. currentLevel <= groundLevel) &
!!$           call whalo1(currentLevel, 1_intType, 0_intType, .true.,&
!!$           .false., .false.)
!!$      
!!$      ! Apply all boundary conditions to all blocks on this level.
!!$      !print *,'second halo',secondhalo,exchangePressureEarly
!!$      call applyAllBC(secondHalo)
!!$      
!!$      ! Exchange the solution. Either whalo1 or whalo2
!!$      ! must be called.
!!$      
!!$      if( secondHalo ) then
!!$         call whalo2(currentLevel, 1_intType, nMGVar, .true., &
!!$              .true., .true.)
!!$      else
!!$         call whalo1(currentLevel, 1_intType, nMGVar, .true., &
!!$              .true., .true.)
!!$      endif
!!$
!!$     ! Reset the values of rkStage and currentLevel, such that
!!$      ! they correspond to a new iteration.
!!$
!!$      rkStage = 0
!!$      currentLevel = groundLevel
!!$
!      ! Compute the latest values of the skin friction velocity.
!      ! The currently stored values are of the previous iteration.
!      
!      call computeUtau
!      
!      ! Apply an iteration to the turbulent transport equations in
!      ! case these must be solved segregatedly.
!      
!      if( turbSegregated ) call turbSolveSegregated
!
!      ! Compute the time step.
!      
!      call timeStep(.false.)
!      
!      ! Compute the residual of the new solution on the ground level.
!      
!      if( turbCoupled ) then
!         call initres(nt1MG, nMGVar)
!         call turbResidual
!      endif!
!
!      call initres(1_intType, nwf)
!      call residual
!
!      ! Get the final time for original routines.
!      call mpi_barrier(SUmb_comm_world, ierr)
!      if(myID == 0) then
!        call cpu_time(time(2))
!        timeOri = time(2)-time(1)
!      endif
!!$
!!$
!!$!     ******************************************************************
!!$!     *                                                                *
!!$!     *  Compute d(dw)/d(w) using Tapenade                             *
!!$!     *                                                                *
!!$!     ******************************************************************
!!$
!!$      ! Get the initial AD time.

!!$      call mpi_barrier(SUmb_comm_world, ierr)
!!$      if( myID==0 ) call cpu_time(time(1))
!!$
!!$      !zero the matrix for dRdW Insert call
!!$      call MatZeroEntries(dRdx,PETScIerr)
!!$
!!$      if( PETScIerr/=0 ) &
!!$        call terminate("verifydRdxFile", "Error in MatZeroEntries drdx")
!!$      print *,'Entering Domain loop'
!!$      domainLoopAD: do nn=1,nDom
!!$         
!!$         ! Loop over the number of time instances for this block.
!!$
!!$         spectralLoop: do sps=1,nTimeIntervalsSpectral
!!$            !print *,'Setting Pointers',nn,level,sps
!!$            call setPointers(nn,level,sps)
!!$
!!$            ! Loop over location of output (R) cell of residual
!!$            do kCell = 2, kl
!!$               do jCell = 2, jl
!!$                  do iCell = 2, il
!!$                     !print *,'indices',icell,jcell,kcell
!!$                     ! Copy the state w to the wAdj array in the stencil
!!$!                     call copyADjointStencil(wAdj, xAdj, iCell, jCell, kCell)                  
!!$                     call copyADjointStencil(wAdj, xAdj,xBlockCornerAdj,alphaAdj,&
!!$           betaAdj,MachAdj,machCoefAdj,machGridAdj,iCell, jCell, kCell,&
!!$           nn,level,sps,pointRefAdj,rotPointAdj,&
!!$           prefAdj,rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
!!$           rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
!!$           murefAdj, timerefAdj,pInfCorrAdj,liftIndex)
!!$                     
!!$
!!$                     !print *,'Stencil Copied'
!!$!                     print *,'wadj',wadj
!!$                     mLoop: do m = 1, nw           ! Loop over output cell residuals (R)
!!$!                        print *,'initializing variables'
!!$                        ! Initialize the seed for the reverse mode
!!$                        dwAdjb(:,:) = 0.; dwAdjb(m,sps) = 1.
!!$                        dwAdj(:,:)  = 0.
!!$                        wAdjb(:,:,:,:,:)  = 0.  !dR(m)/dw
!!$                        xAdjb(:,:,:,:,:)  = 0.  !dR(m)/dw
!!$                        xblockcorneradjb  = 0.
!!$                        alphaadjb = 0.
!!$                        betaadjb = 0.
!!$                        machadjb = 0.
!!$!                        print *,'dwadjb',dwadjb,'wadjb',wadjb(0,0,0,:)
!!$
!!$  !                      print *,'calling reverse mode'
!!$!                        print *,'secondhalo',secondhalo
!!$                        	        ! Call the reverse mode of residual computation.
!!$                !
!!$                !                          dR(iCell,jCell,kCell,l)
!!$                ! wAdjb(ii,jj,kk,n) = --------------------------------
!!$                !                     dW(iCell+ii,jCell+jj,kCell+kk,n)
!!$ !                       print *,'input',wadj, wadjb, xadj, xadjb, dwadj, dwadjb, &
!!$ !                            &  alphaadj, alphaadjb, betaadj, betaadjb, machadj, machadjb, &
!!$ !                            &  machcoefadj, icell, jcell, kcell, nn, sps, correctfork, secondhalo, &
!!$ !                            &  prefadj, rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, &
!!$ !                            &  murefadj, timerefadj, pinfcorradj   
!!$                        !print *,'calling reverse mode'
!!$                        ! Call reverse mode of residual computation
!!$                        call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
!!$   &  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
!!$   &  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
!!$   &  icell, jcell, kcell, nn, level, sps, correctfork, secondhalo, prefadj&
!!$   &  , rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj&
!!$   &  , rotrateadjb, rotcenteradj, rotcenteradjb, pointrefadj, rotpointadj&
!!$   &  , rotpointadjb, murefadj, timerefadj, pinfcorradj, liftindex)
!!$
!!$                        ! Store the block Jacobians (by rows).
!!$                        !print *,'entering storage loop'
!!$                        do sps2 = 1,nTimeIntervalsSpectral
!!$                        do ii=-3,2!1,il-1
!!$                           do jj = -3,2!1,jl-1
!!$                              do kk = -3,2!1,kl-1
!!$                                 do l = 1,3
!!$                                    i = iCell + ii
!!$                                    j = jCell + jj
!!$                                    k = kCell + kk
!                                    !print *,'secondaryindicies',i,j,k,ii,jj,kk
!                                    if(i>zero .and. j>zero .and. k>zero .and. i<=il .and. j<=jl .and. k<=kl)then
!                                       idxstate = globalCell(i,j,k)*nw+l
!                                       idxres   = globalCell(iCell,jCell,kCell)*nw+m
!                                       !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
!                                       if( idxres>=0 .and. idxstate>=0) then
!                                          dRdwAdj(idxres,idxstate,nn,sps) = wAdjb(ii,jj,kk,l)!
!
!                                       endif
!                                    end if
!!$
!!$                                    if(i>=zero .and. j>=zero .and. k>=zero .and. i<=ie .and. j<=je .and. k<=ke)then
!!$                                       call setPointers(nn,level,sps2)
!!$                                       idxnode = globalnode(i,j,k)*3+l
!!$                                       call setPointers(nn,level,sps)
!!$                                       idxres   = globalCell(iCell,jCell,kCell)*nw+m
!!$                                       if( (idxres-1)>=0 .and. (idxnode-1)>=0) then
!!$                                          if (xAdjb(ii,jj,kk,l,sps2).ne.0.0)then
!                                          write(unitdrdx,13) idxnode,idxres,l,i,j,k,nn,m,kcell,jcell,icell,xAdjb(ii,jj,kk,l),xAdj(ii,jj,kk,l)!idxnode,idxres,m,icell,jcell,kcell,nn,n,k,j,i,nnn,value
!                                          !write(unitWarp,13) xderiv,i,j,k,n,nnn,nn,mm,ll
!13                                        format(1x,'drdx',11I8,2f18.10)
!                                             print 14,idxnode,idxres,l,i,j,k,nn,m,kcell,jcell,icell,xAdjb(ii,jj,kk,l)
!14                                           format(1x,'drdx',11I8,f18.10)
!!$                                             call MatSetValues(drdx, 1, idxres-1, 1, idxnode-1,   &
!!$                                                  xAdjb(ii,jj,kk,l,sps2), ADD_VALUES, PETScIerr)
!!$                                             if( PETScIerr/=0 ) &
!!$                                                  print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
!!$                                          endif
!!$                                       endif
!!$                                    endif                                   
!!$                                    
!!$                                 enddo !l
!!$                              enddo !kk
!!$                           enddo !jj
!!$                        enddo !ii
!!$                     
!!$                        !set values for symmtery plane normal derivatives
!!$                        do l = 1,3
!!$                           if (xblockcorneradjb(1,1,1,l,sps2).ne.0.0)then
!!$                              call setPointers(nn,level,sps2)
!!$                              idxnode = globalnode(1,1,1)*3+l
!!$                              call setPointers(nn,level,sps)
!!$                              call MatSetValues(drdx, 1, idxres-1, 1, idxnode-1,   &
!!$                                   xblockcorneradjb(1,1,1,l,sps2), ADD_VALUES, PETScIerr)
!!$                              if( PETScIerr/=0 ) &
!!$                                   print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
!!$                           endif
!!$                           if (xblockcorneradjb(2,1,1,l,sps2).ne.0.0)then
!!$                              call setPointers(nn,level,sps2)
!!$                              idxnode = globalnode(il,1,1)*3+l
!!$                              call setPointers(nn,level,sps)
!!$                              call MatSetValues(drdx, 1, idxres-1, 1, idxnode-1,   &
!!$                                   xblockcorneradjb(2,1,1,l,sps2), ADD_VALUES, PETScIerr)
!!$                              if( PETScIerr/=0 ) &
!!$                                   print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
!!$                           endif
!!$                           if (xblockcorneradjb(1,2,1,l,sps2).ne.0.0)then
!!$                              call setPointers(nn,level,sps2)
!!$                              idxnode = globalnode(1,jl,1)*3+l
!!$                              call setPointers(nn,level,sps)
!!$                              call MatSetValues(drdx, 1, idxres-1, 1, idxnode-1,   &
!!$                                   xblockcorneradjb(1,2,1,l,sps2), ADD_VALUES, PETScIerr)
!!$                              if( PETScIerr/=0 ) &
!!$                                   print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
!!$                           endif
!!$                           if (xblockcorneradjb(2,2,1,l,sps2).ne.0.0)then
!!$                              call setPointers(nn,level,sps2)
!!$                              idxnode = globalnode(il,jl,1)*3+l
!!$                              call setPointers(nn,level,sps)
!!$                              call MatSetValues(drdx, 1, idxres-1, 1, idxnode-1,   &
!!$                                   xblockcorneradjb(2,2,1,l,sps2), ADD_VALUES, PETScIerr)
!!$                              if( PETScIerr/=0 ) &
!!$                                   print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
!!$                           endif
!!$                           if (xblockcorneradjb(1,1,2,l,sps2).ne.0.0)then
!!$                              call setPointers(nn,level,sps2)
!!$                              idxnode = globalnode(1,1,kl)*3+l
!!$                              call setPointers(nn,level,sps)
!!$                              call MatSetValues(drdx, 1, idxres-1, 1, idxnode-1,   &
!!$                                   xblockcorneradjb(1,1,2,l,sps2), ADD_VALUES, PETScIerr)
!!$                              if( PETScIerr/=0 ) &
!!$                                   print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
!!$                           endif
!!$                           if (xblockcorneradjb(1,2,2,l,sps2).ne.0.0)then
!!$                              call setPointers(nn,level,sps2)
!!$                              idxnode = globalnode(1,jl,kl)*3+l
!!$                              call setPointers(nn,level,sps)
!!$                              call MatSetValues(drdx, 1, idxres-1, 1, idxnode-1,   &
!!$                                   xblockcorneradjb(1,2,2,l,sps2), ADD_VALUES, PETScIerr)
!!$                              if( PETScIerr/=0 ) &
!!$                                   print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
!!$                           endif
!!$                           if (xblockcorneradjb(2,1,2,l,sps2).ne.0.0)then
!!$                              call setPointers(nn,level,sps2)
!!$                              idxnode = globalnode(il,1,kl)*3+l
!!$                              call setPointers(nn,level,sps)
!!$                              call MatSetValues(drdx, 1, idxres-1, 1, idxnode-1,   &
!!$                                   xblockcorneradjb(2,1,2,l,sps2), ADD_VALUES, PETScIerr)
!!$                              if( PETScIerr/=0 ) &
!!$                                   print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
!!$                           endif
!!$                           if (xblockcorneradjb(2,2,2,l,sps2).ne.0.0)then
!!$                              call setPointers(nn,level,sps2)
!!$                              idxnode = globalnode(il,jl,kl)*3+l
!!$                              call setPointers(nn,level,sps)
!!$                              call MatSetValues(drdx, 1, idxres-1, 1, idxnode-1,   &
!!$                                   xblockcorneradjb(2,2,2,l,sps2), ADD_VALUES, PETScIerr)
!!$                              if( PETScIerr/=0 ) &
!!$                                   print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
!!$                           endif
!!$                        enddo
!!$
!!$                     end do
!!$
!!$                     end do mLoop
!!$                     !stop
!!$                  end do !iCell
!!$               end do !jCell
!!$            end do! kCell
!!$         end do spectralLoop
!!$      end do domainLoopAD
!!$      print *,'AD Completed'!,'indices',il,jl,kl
!!$!      stop
!!$      ! Get new time and compute the elapsed AD time.
!!$!      stop
!!$      call mpi_barrier(SUmb_comm_world, ierr)
!!$      if(myID == 0) then
!!$         call cpu_time(time(2))
!!$         timeAdj = time(2)-time(1)
!!$      endif
!!$      
!!$      call MatAssemblyBegin(drdx,MAT_FINAL_ASSEMBLY,PETScIerr)
!!$      
!!$      if( PETScIerr/=0 ) &
!!$           call terminate("verifydrdxFile","Error in MatAssemblyBegin")
!!$      
!!$      
!!$      call MatAssemblyEnd(drdx,MAT_FINAL_ASSEMBLY,PETScIerr)
!!$       
!!$       if( PETScIerr/=0 ) &
!!$            call terminate("verifydrdxFile","Error in MatAssemblyEnd")
!!$
!!$       !if( debug ) then
!!$          !call MatView(drdx,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
!!$          !call MatView(drdx,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
!!$          !if( PETScIerr/=0 ) &
!!$          !     call terminate("setupADjointMatrix", "Error in MatView")
!!$          !pause
!!$       !endif
!!$       
!       call cpu_time(time(4))
!       timeFD = time(4)-time(3)

      print *, 'in verifydRdxFile'
      call MatAssembled(drdx,assembled,ierr)
      call EChk(ierr,__FILE__,__LINE__)
      print *,'assembled spatial =',assembled
!      !now extract and write to a file
       do sps = 1,nTimeIntervalsSpectral
       do nn = 1,nDom
          call setPointers(nn,1,sps)
          do kCell = 1,kl!0, ke
             do jCell = 1,jl!0, je
                do iCell = 1,il!0, ie
                   !print *,'ie',icell,ie,jcell,je,kcell,ke
                   do m = 1, 3
                      do sps2 = 1,nTimeIntervalsSpectral
                      call setPointers(nn,1,sps)
                     idxnode   = globalnode(iCell,jCell,kCell)*3+m 
                     testnode = globalnode(iCell,jCell,kCell)
                     call setPointers(nn,1,sps2)
                     do nnn = 1,ndom
                        call setPointers(nnn,1,sps2)
                        DO I=2,Il
                           DO J=2,Jl
                              DO K=2,Kl
                                 do n = 1,nw
                                    idxres = globalCell(i,j,k)*nw+n
                                    if ((idxres-1)>=0 .and. (idxnode-1)>=0)then
                                       call MatGetValues(drdx,1,idxnode-1,1,idxres-1,value1,PETScIerr)
                                       !call MatGetOwnershipRange(drdx,mmmm,nnnn)
                                       !print *, mmmm,nnnn
                                       !print *,'idxres = ', idxres, 'idxnode =',idxnode
                                       !call EChk(PETScIerr,__FILE__,__LINE__)
                                       !print *,'value',value,i,j,k
                                       !if(value.ne.0)then
                                       if(abs(value1)>1e-10)then
                                          !write(unitWarp,12)ifaceptb,iedgeptb !'face',ifaceptb,'edge',iedgeptb
                                          !12                                     format(1x,'Face',6I2,'edge',12I2)
                                          write(unitdrdx,13) idxnode,idxres,m,icell,jcell,kcell,nn,sps,n,k,j,i,nnn,sps2,value1
                                          !print*,idxnode,idxres,m,icell,jcell,kcell,nn,sps,n,k,j,i,nnn,sps2,value1
                                          !write(unitWarp,13) xderiv,i,j,k,n,nnn,nn,mm,ll
13                                        format(1x,'drdx',14I8,f18.10)
                                       endif
                                    end if
                                 enddo
                              END DO
                           END DO
                        END DO
                        call setPointers(nn,1,sps2)
                     end do
                  end do
               enddo
            end do
         end do
      enddo
   end do
end do
!Print *,'barriercall',myID
call mpi_barrier(SUmb_comm_world, ierr)

 print *, 'closing unitdrdx'
 close(unitdrdx)

 call mpi_barrier(SUmb_comm_world, ierr)


 111  format(4I4, 3ES22.14)

#endif

end subroutine verifydRdxfile
