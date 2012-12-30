!
!     ******************************************************************
!     *                                                                *
!     * File:          verifydRdx.f90                                  *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 05-05-2008                                      *
!     * Last modified: 06-17-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine verifydRdx(level)
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
      use inputTimeSpectral ! spaceDiscr
      use inputIO

      !from old verify routine
      use ADjointPETSc, only: pvr,drdx,drdxfd,petscone,insert_values,petscierr,mat_final_assembly,petsc_viewer_draw_world,petsc_viewer_stdout_world
      use ADjointVars
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
      integer(kind=intType) :: i, j, k, n
      integer(kind=intType) :: iCell,jCell,kCell
      real(kind=realType), dimension(nw) :: dwL2
      real(kind=realType), dimension(nx, ny, nz, nw) :: dwerr
      real(kind=realType), dimension(10) :: time
      real(kind=realType) :: timeRes
!      real(kind=realType), dimension(4) :: time
      real(kind=realType) ::  timeOri
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdj
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdjb
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wFD
      real(kind=realType), dimension(-3:2,-3:2,-3:2,3)  :: xAdj
      real(kind=realType), dimension(-3:2,-3:2,-3:2,3)  :: xAdjb
      real(kind=realType), dimension(-3:2,-3:2,-3:2,3)  :: xFD

      REAL(KIND=REALTYPE) :: machadj, machcoefadj, pinfcorradj
      REAL(KIND=REALTYPE) :: machadjb, machcoefadjb, pinfcorradjb
      REAL(KIND=REALTYPE) :: prefadj, rhorefadj
      REAL(KIND=REALTYPE) :: pinfdimadj, rhoinfdimadj
      REAL(KIND=REALTYPE) :: rhoinfadj, pinfadj
      REAL(KIND=REALTYPE) :: rhoinfadjb
      REAL(KIND=REALTYPE) :: murefadj, timerefadj
      REAL(KIND=REALTYPE) :: alphaadj, betaadj
      REAL(KIND=REALTYPE) :: alphaadjb, betaadjb

      REAL(KIND=REALTYPE), DIMENSION(3) :: rotcenteradj
      REAL(KIND=REALTYPE), DIMENSION(3) :: rotrateadj
      REAL(KIND=REALTYPE) :: rotrateadjb(3)

      character fileName*32, dataName*32
      real(kind=realType), dimension(nw) :: dwAdj,dwAdjb,dwAdjRef
      real(kind=realType), dimension(nw) :: dwAdjP, dwAdjM
      real(kind=realType) :: deltax, xAdjRef
      real(kind=realType) :: timeAdj, timeFD, timeResAdj,test

      integer(kind=intType), dimension(nDom) :: maxglobalcell,maxglobalnode
      integer(kind=intType) :: idx, ii, jj, kk, idxnode, idxres, m, l,idxmgb
      integer(kind=intType) :: sps, nTime, max_nTime, nHalo, nn, discr
      real(kind=realType), allocatable, dimension(:,:,:,:) :: dRdxAdj,dRdxFD1,dRdxFD2,dRdxErr

      integer :: ierr,liftIndex
      logical :: fineGrid, correctForK, exchangeTurb,secondHalo

      !FULL FD Variables
      
      !real(kind=realType), dimension(nx, ny, nz, nw,nx) :: dRdwErr, dRdwErrRel
      real(kind=realType), dimension(3) :: dRdxL2, dRdxL2Rel
!!$      real(kind=realtype), dimension(0:ib,0:jb,0:kb,1:nw)::dwp,dwm,dwtemp
!!$!      real(kind=realType), dimension(0:ib,0:jb,0:kb,1:nw) :: wtemp
!!$      real(kind=realType), dimension(0:ie,0:je,0:ke,1:3) :: xtemp
!!$      real(kind=realType), dimension(0:ib,0:jb,0:kb) :: ptemp
      !real(kind=realtype),allocatable, dimension(:,:,:,:)::dwp,dwm,dwtemp
!      real(kind=realType), dimension(0:ib,0:jb,0:kb,1:nw) :: wtemp
      real(kind=realType),allocatable, dimension(:,:,:,:) :: xtemp
      real(kind=realType),allocatable, dimension(:,:,:) :: ptemp
      integer(kind=intType) :: istate, jstate, kstate,ires
      real(kind=realType), dimension(-2:2,-2:2,-2:2) :: pAdjtemp
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdjtemp
      real(kind=realType),allocatable, dimension(:,:,:,:) :: xFD2
      !real(kind=realType), dimension(nx,ny,nz,nw) :: xFD2
      !real(kind=realType), dimension(ib*jb*kb*nw,ib*jb*kb*nw) :: dRdw
      !real(kind=realType), dimension(ib*jb*kb*nw,ib*jb*kb*nw) :: dRdwFD

     ! pvr row block
      
      real(kind=realType), dimension(nw) :: pvrlocal
      character(len=2*maxStringLen) :: errorMessage


!     ******************************************************************
!     *                                                                *
!     *  Begin execution                                               *
!     *                                                                *
!     ******************************************************************

      if( myID==0 ) write(*,*) "Running verifydRdx..."

      currentLevel = level
      !discr        = spaceDiscr
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

      !allocate memory for error arrays
      max_nTime=-100
      do i = 1,nDom
         print *,'allocating i'
         maxglobalcell(i) = maxval(flowDoms(i,currentLevel,1)%globalCell(:,:,:))
         print *,'maxglobalcell',maxglobalcell(i) 
         maxglobalnode(i) = maxval(flowDoms(i,currentLevel,1)%globalnode(:,:,:))
         print *,'maxglobalnode',maxglobalnode(i)
         nTime     = nTimeIntervalsSpectral!sections(sectionID)%nTimeInstances
         if(nTime>=max_nTime) max_nTime = nTime
      enddo
      idx = maxval(maxglobalcell(:))
      idxnode = maxval(maxglobalnode(:))
      print *,'allocating',idx,nw*(idx+1),3*(idxnode+1),ndom,max_nTime
       allocate(dRdxErr(nw*(idx+1),3*(idxnode+1),1,max_nTime), &
               dRdxAdj(nw*(idx+1),3*(idxnode+1),1,max_nTime), &
               dRdxFD1(nw*(idx+1),3*(idxnode+1),1,max_nTime), &
               dRdxFD2(nw*(idx+1),3*(idxnode+1),1,max_nTime))
      print *,' allocated'

!           if(ierr /= 0)                       &
!                call terminate("memory?") 
      
      ! Initialize the temporary arrays.
      dRdxErr = 0
      dRdxAdj = 0
      dRdxFD1  = 0
      dRdxFD2  = 0

    !allocate memory for FD
      allocatedomains: do nn = 1,ndom
         print *,'domain',nn
         groundLevel = 1
         sps = 1
         call setPointers(nn,1,sps)
         allocate(flowDoms(nn,level,sps)%dwp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
         allocate(flowDoms(nn,level,sps)%dwm(0:ib,0:jb,0:kb,1:nw),stat=ierr)
         allocate(flowDoms(nn,level,sps)%dwtemp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
      end do allocatedomains


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
      print *,'second halo',secondhalo,exchangePressureEarly
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


!     ******************************************************************
!     *                                                                *
!     *  Compute d(dw)/d(w) using Tapenade                             *
!     *                                                                *
!     ******************************************************************

      ! Get the initial AD time.

      call mpi_barrier(SUmb_comm_world, ierr)
      if( myID==0 ) call cpu_time(time(1))
      
      print *,'Entering Domain loop'
      domainLoopAD: do nn=1,nDom
         
         ! Loop over the number of time instances for this block.

         spectralLoop: do sps=1,nTimeIntervalsSpectral
            print *,'Setting Pointers',nn,level,sps
            call setPointers(nn,level,sps)

            ! Loop over location of output (R) cell of residual
            do kCell = 2, kl
               do jCell = 2, jl
                  do iCell = 2, il
 !                    print *,'indices',icell,jcell,kcell
                     ! Copy the state w to the wAdj array in the stencil
!                     call copyADjointStencil(wAdj, xAdj, iCell, jCell, kCell)                  
                     call copyADjointStencil(wAdj, xAdj,alphaAdj,betaAdj,MachAdj,&
                          machCoefAdj,iCell, jCell, kCell,prefAdj,&
                          rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                          rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
                          murefAdj, timerefAdj,pInfCorrAdj,liftIndex)

!                     print *,'Stencil Copied'
!                     print *,'wadj',wadj
                     mLoop: do m = 1, nw           ! Loop over output cell residuals (R)
!                        print *,'initializing variables'
                        ! Initialize the seed for the reverse mode
                        dwAdjb(:) = 0.; dwAdjb(m) = 1.
                        dwAdj(:)  = 0.
                        wAdjb(:,:,:,:)  = 0.  !dR(m)/dw
                        xAdjb(:,:,:,:)  = 0.  !dR(m)/dx
                        alphaadjb = 0.
                        betaadjb = 0.
                        machadjb = 0.
!                        print *,'dwadjb',dwadjb,'wadjb',wadjb(0,0,0,:)

  !                      print *,'calling reverse mode'
!                        print *,'secondhalo',secondhalo
                        	        ! Call the reverse mode of residual computation.
                !
                !                          dR(iCell,jCell,kCell,l)
                ! wAdjb(ii,jj,kk,n) = --------------------------------
                !                     dW(iCell+ii,jCell+jj,kCell+kk,n)
 !                       print *,'input',wadj, wadjb, xadj, xadjb, dwadj, dwadjb, &
 !                            &  alphaadj, alphaadjb, betaadj, betaadjb, machadj, machadjb, &
 !                            &  machcoefadj, icell, jcell, kcell, nn, sps, correctfork, secondhalo, &
 !                            &  prefadj, rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, &
 !                            &  murefadj, timerefadj, pinfcorradj   
                        
                        ! Call reverse mode of residual computation
                        call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, dwadj, dwadjb, &
                             &  alphaadj, alphaadjb, betaadj, betaadjb, machadj, machadjb, &
                             &  machcoefadj, icell, jcell, kcell, nn, sps, correctfork, secondhalo, &
                             &  prefadj, rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, &
                             &  rotrateadj, rotrateadjb, rotcenteradj, murefadj, timerefadj, &
                             &  pinfcorradj, liftindex)
                       ! call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb,&
                       !      dwadj, dwadjb, icell, jcell, kcell, nn, sps,&
                       !      correctfork, secondhalo)
                        !print *,'wadjb',wAdjB( 0, 0, 0,:)
                        ! Store the block Jacobians (by rows).
 !                       print *,'entering storage loop'
                        do ii=-3,2!1,il-1
                           do jj = -3,2!1,jl-1
                              do kk = -3,2!1,kl-1
                                 do l = 1,3
                                    i = iCell + ii
                                    j = jCell + jj
                                    k = kCell + kk
                                    !print *,'secondaryindicies',i,j,k,ii,jj,kk
                                    if(i>zero .and. j>zero .and. k>zero .and. i<=ie .and. j<=je .and. k<=ke)then
                                       idxnode = globalNode(i,j,k)*3+l
                                       idxres   = globalCell(iCell,jCell,kCell)*nw+m
                                       !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
                                       if( idxres>=0 .and. idxnode>=0) then
                                          !dRdxAdj(idxres,idxnode,nn,sps) = xAdjb(ii,jj,kk,l)
                                          dRdxAdj(idxres,idxnode,1,1) = xAdjb(ii,jj,kk,l)
!!$                                          if (icell==13.and.jcell==5.and.kcell==5.and.ii==0.and.jj==-1.and.kk==0)then
!!$                                             print *,'xAdjb', xAdjb(ii,jj,kk,l),icell,jcell,kcell,m,ii,jj,kk,l,i,j,k
!!$                                             stop
!!$                                          endif
                                          
                                           !if ( xAdjb(ii,jj,kk,l).ne. 0.0) then
                                           !   print *,'xAdjb', xAdjb(ii,jj,kk,l),ii,jj,kk,l
!!$                                              if (i==12) then
!!$                                                 print *,'xAdjb', xAdjb(ii,jj,kk,l),icell,jcell,kcell,m,ii,jj,kk,l,i,j,k
!!$                                              endif
                                           !endif
                                       endif
                                    end if
                                 enddo !l
                              enddo !kk
                           enddo !jj
                        enddo !ii
                        
                     end do mLoop
                     !stop
                  end do !iCell
               end do !jCell
            end do! kCell
         end do spectralLoop
      end do domainLoopAD
      print *,'AD Completed'!,'indices',il,jl,kl
!      stop
      ! Get new time and compute the elapsed AD time.
!      stop
      call mpi_barrier(SUmb_comm_world, ierr)
      if(myID == 0) then
         call cpu_time(time(2))
         timeAdj = time(2)-time(1)
      endif

!
!     Compute d(dw)/d(x) using central finite differences
!_______________________________________________________
       
       deltax = 1.d-7
       print *, "deltax=", deltax
       xFD2(:,:,:,:) = 0. 
       call cpu_time(time(3))
       
       groundLevel = 1
       sps = 1
       nn=1

       print *,'Entering Domain loop'
       domainLoopFD1: do nn=1,1!nDom
          
          ! Loop over the number of time instances for this block.
          
          spectralLoop1: do sps=1,nTimeIntervalsSpectral
             print *,'Setting Pointers',nn,level,sps
             call setPointers(nn,level,sps)
             
             allocate(dwp(0:ib,0:jb,0:kb,1:nw),dwm(0:ib,0:jb,0:kb,1:nw),dwtemp(0:ib,0:jb,0:kb,1:nw))
             allocate(xtemp(0:ie,0:je,0:ke,1:3) ,ptemp(0:ib,0:jb,0:kb),xFD2(nx,ny,nz,nw))

             ! Loop over the inputs (x) 
             do n = 1, 3
                do kstate = 0, ke
                   print *, "iw, kstate =", n, kstate
                   do jstate = 0, je
                      do istate = 0, ie
                         !print *,'Setting Pointers loop',nn,level,sps
                         call setPointers(nn,level,sps)

                         !print *, "iw, istate =", n, istate
                         ! if (istate==icell .or. jstate==jcell .or. kstate==kcell) then
                         !Remember current values
                         dwtemp(:,:,:,:) = dw(:,:,:,:)
                         xtemp(:,:,:,:) = x(:,:,:,:)
                         ptemp(:,:,:) = p(:,:,:)
                         xAdjRef = x(istate,jstate,kstate,n)
                         
                         !print *,'x', x(istate,jstate,kstate,n), xAdjRef, deltax
                         !print *,'x12,2,2', x(12,2,2,n),x(11,2,2,n)
                         x(istate,jstate,kstate,n) = xAdjRef + deltax
!!$                   wAdjtemp = w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:)
!!$!                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$                   call computePressureAdj(wAdjtemp, pAdjtemp)
!!$                   p(istate, jstate, kstate) = pAdjtemp(0,0,0)
                         call xhalo(groundlevel)
                         call metric(groundlevel)
                         !call setPointers(nn,level,sps)
                         call checkSymmetry(groundlevel)
                         !call setPointers(nn,level,sps)
                         !print *, "Calling Residual ="
                         call initres(1_intType, nwf)
                         !call setPointers(nn,level,sps)
                         call applyAllBC(secondHalo)
                         !!call setPointers(nn,level,sps)
                         ! Exchange the solution. Either whalo1 or whalo2
                         ! must be called.
                         if( secondHalo ) then
                            call whalo2(currentLevel, 1_intType, nMGVar, .true., &
                                 .true., .true.)
                         else
                            call whalo1(currentLevel, 1_intType, nMGVar, .true., &
                                 .true., .true.)
                         endif
                         !call setPointers(nn,level,sps)
                         !if( secondHalo ) then
                         !   !  write(*,*)'2ndHalo..........'
                         !   call whalo2(currentLevel, 1_intType, nVarInt, .true., &
                         !        .true., .true.)
                         !else
                         !   ! write(*,*)'1stHalo..........'
                         !   call whalo1(currentLevel, 1_intType, nVarInt, .true., &
                         !        .true., .true.)
                         !endif
                         
                         call residual
                         call setPointers(nn,level,sps)
                         !print *, "Called Residual =", dw(istate,jstate,kstate,n)
                         
                         dwp(:,:,:,:) = dw(:,:,:,:)
                         dw(:,:,:,:) = dwtemp(:,:,:,:)
                         
                         x(istate,jstate,kstate,n) = xAdjRef - deltax
!!$                   wAdjtemp = w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:)
!!$!                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$                   call computePressureAdj(wAdjtemp, pAdjtemp)
!!$!                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$                   p(istate, jstate, kstate) = pAdjtemp(0,0,0)
                         call xhalo(groundlevel)
                         call metric(groundlevel)    
                         !call setPointers(nn,level,sps)
                         call checkSymmetry(groundlevel)
                         !call setPointers(nn,level,sps)
                         !print *, "Calling Residual ="
                         call initres(1_intType, nwf)
                        ! call setPointers(nn,level,sps)
                         call applyAllBC(secondHalo)
                         !call setPointers(nn,level,sps)
                         ! Exchange the solution. Either whalo1 or whalo2
                         ! must be called.
                         if( secondHalo ) then
                            call whalo2(currentLevel, 1_intType, nMGVar, .true., &
                                 .true., .true.)
                         else
                            call whalo1(currentLevel, 1_intType, nMGVar, .true., &
                                 .true., .true.)
                         endif
                         !call setPointers(nn,level,sps)
                         !if( secondHalo ) then
                         !   !write(*,*)'2ndHalo..........'
                         !   call whalo2(currentLevel, 1_intType, nVarInt, .true., &
                         !        .true., .true.)
                         !else
                         !   !write(*,*)'1stHalo..........'
                         !   call whalo1(currentLevel, 1_intType, nVarInt, .true., &
                         !        .true., .true.)
                         !endif
                         
                         call residual
                         call setPointers(nn,level,sps)
                         !print *, "Called Residual 2 =", dw(istate,jstate,kstate,n)
                         dwm(:,:,:,:) = dw(:,:,:,:)
                         dw(:,:,:,:) = dwtemp(:,:,:,:)
                         
                         !                             end if
                         
                         
                         !                   ival = 0
                         call setPointers(nn,level,sps)
                         
                         ! Loop over location of output (R) cell of residual
                         do kCell = 2, kl
                            do jCell = 2, jl
                               do iCell = 2, il
                                  ! Loop over output cell residuals (R)
                                  !print *,'ires',ires
                                  do iRes = 1, nw  
                                     !print *,'xfd2',shape(xfd2),icell,jcell,kcell,ires
                                     !column = (n-1) +(istate-2)*nw +(jstate-2)*nw*nx +(kstate-2)*nw*nx*ny
                                     xFD2(iCell-1, jCell-1, kCell-1, iRes) = (dwp(iCell,jCell,kCell,iRes)-dwm(iCell,jCell,kCell,iRes))/(2.0*deltax)
                                    
                                     !if (abs(wFD(iCell-1, jCell-1, kCell-1, iRes)) >1e-12) then
                                     !   ival = ival + 1
                                     !   row2(ival) = (ires-1) + (icell-2)*nw + (jcell-2)*nw*nx + (kcell-2)*nw*nx*ny
                                     !   col_values(ival)  = wFD(iCell-1, jCell-1, kCell-1, iRes)
                                     !   
                                     !  
                                     !  
                                     !end if
                                     !if(i>zero .and. j>zero .and. k>zero .and. i<il .and. j<jl .and. k<kl)then
                                     idxnode = globalNode(istate,jstate,kstate)*3+n
                                     idxres   = globalCell(iCell,jCell,kCell)*nw+ires
                                     !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
                                     if( idxres>=0 .and. idxnode>=0) then
                                        if (xFD2(iCell-1, jCell-1, kCell-1, iRes).ne. zero) then
                                           !print *,'xFD2',xFD2(iCell-1, jCell-1, kCell-1, iRes),icell,jcell,kcell
                                           !dRdxFD1(idxres,idxnode,nn,sps) = xFD2(iCell-1, jCell-1, kCell-1, iRes)
                                           dRdxFD1(idxres,idxnode,1,1) = xFD2(iCell-1, jCell-1, kCell-1, iRes)
                                           !if( istate ==1 .and. jstate==4 .and. kstate ==4 .and. n==3 .and. nn==1) then
                                           if(idxnode==flowdoms(1,1,1)%globalNode(1,4,4)*3+3) then
!                                           if( istate ==1 .and. jstate==1 .and. kstate ==1 .and. n==3 .and. nn==1) then
                                              if (  idxres<(nCellsGlobal(1_intType)*nw+nw)) then
                                                 print *,'ncells global',nCellsGlobal,nCellsGlobal(1_intType)*nw+nw,idxres,globalCell(iCell,jCell,kCell)
                                                 call VecSetValue(pvr, idxres-1, xFD2(iCell-1, jCell-1, kCell-1, iRes) ,INSERT_VALUES, PETScIerr)
                                                 if( PETScIerr/=0 ) then
                                                    write(errorMessage,99) &
                                                         "Error in VecSetValuesBlocked for global node", &
                                                         idxmgb
                                                    call terminate("setupADjointRHSAeroCoeff", &
                                                         errorMessage)
                                                 endif
                                              endif
                                           endif
                                           !Aad(m,:)  = wFD2(iCell-1, jCell-1, kCell-1, iRes)wAdjB( 0, 0, 0,:)
                                           !print *,'setting',dRdW, PETScOne, idxres, PETScOne, idxstate,   &
                                           !     dRdwFD(idxres,idxstate,nn,sps)
                                        endif
                                     endif
                                  enddo
!!$                                        call MatSetValues(dRdxFD, 1, idxres-1, 1, idxnode-1,   &
!!$                                             dRdxFD1(idxres,idxnode,nn,sps), INSERT_VALUES, PETScIerr)
!!$                                        if( PETScIerr/=0 ) &
!!$                                             print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
                                           
                                           
!!$                                        call MatSetValues(dRdW, 1, idxres-1, 1, idxstate-1,   &
!!$                                             dRdwFD(idxres,idxstate,nn,sps), INSERT_VALUES, PETScIerr)
!!$                                        if( PETScIerr/=0 ) &
!!$                                             print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")

!*************                                           
!!$                                  idxmgb = globalCell(icell,jcell,kcell)
!!$                                  
!!$                                  test = sum(xFD2(iCell-1, jCell-1, kCell-1,:))
!!$                                  if( istate ==1 .and. jstate==1 .and. kstate ==1 .and. n==3 .and. nn==1) then
!!$                                     !if( istate ==1 .and. jstate==1 .and. kstate ==1 .and. n==3) then
!!$                                     !print *,'test',test
!!$                                     if ( test.ne.0 .and. idxmgb.ne.-5 .and. idxmgb>=0 .and. idxmgb<nCellsGlobal) then
!!$                                        !print *,'setting PETSc Vector',sum(wAdjB(icell,jcell,kcell,:))
!!$                                        print *,'ncells global',nCellsGlobal,nCellsGlobal*nw+nw,idxres,idxmgb
!!$                                        pvrlocal(:) = xFD2(iCell-1, jCell-1, kCell-1,:)
!!$                                        
!!$                                        !                call VecSetValuesBlocked(dJdW, 1, idxmgb, dJdWlocal, &
!!$                                        !                                         INSERT_VALUES, PETScIerr)
!!$                                        call VecSetValuesBlocked(pvr, 1, idxmgb, pvrlocal, &
!!$                                             INSERT_VALUES, PETScIerr)
!!$                                        
!!$                                        if( PETScIerr/=0 ) then
!!$                                           write(errorMessage,99) &
!!$                                                "Error in VecSetValuesBlocked for global node", &
!!$                                                idxmgb
!!$                                           call terminate("setupADjointRHSAeroCoeff", &
!!$                                                errorMessage)
!!$                                        endif
!!$                                     endif
!!$                                  endif
!*****************                                  
                                  !idxmgb = globalCell(iCell,jCell,kCell)
                                  !print *,'globalcell',idxmgb,globalCell(iCell,jCell,kCell)
                                  ! >>> center block A < W(i,j,k)
                                           
                                  !idxngb = idxmgb
                                  !print *,'indicies0',idxmgb,idxngb
                                  !call MatSetValuesBlocked(dRdW, 1, idxmgb, 1, idxngb, &
                                  !     Aad, INSERT_VALUES,PETScIerr)
                                  !if( PETScIerr/=0 ) &
                                  !      print *,'matrix setting error'!call errAssemb("MatSetValuesBlocked", "Aad")
                                  
                                  ! if ( wAdjb(ii,jj,kk,l).ne. 0.0) then
                                  !    print *,'wAdjb', wAdjb(ii,jj,kk,l),i,j,k,l
                                  ! endif
                                  !endif
                                  !endif
                                  !end if
                                  
                               end do
                            end do
                         end do
                         
                         call setPointers(nn,level,sps)
                         !print *,'ressetting reference'
                         x(istate,jstate,kstate,n) = xAdjRef
                         !print *,'xreset', x(istate,jstate,kstate,n), xAdjRef, deltax
                         !print *,'x12,2,2 reset', x(12,2,2,n),x(11,2,2,n)
                         dw(:,:,:,:) = dwtemp(:,:,:,:)
                         x(:,:,:,:) = xtemp(:,:,:,:)
                         p(:,:,:) = ptemp(:,:,:)
!!$                   
!!$                   col_size = ival
!!$                   
!!$                   call MatSetValues(dRdWFD, col_size, row2, PetscOne, column, col_values, INSERT_VALUES, PETScIerr)
                      end do
                   end do
                end do
             end do
             deallocate(dwp,dwm,dwtemp)
             deallocate(xtemp ,ptemp,xFD2)
          enddo spectralLoop1
       enddo domainLoopFD1
!!$
!!$       call MatAssemblyBegin(dRdxFD,MAT_FINAL_ASSEMBLY,PETScIerr)
!!$
!!$       if( PETScIerr/=0 ) &
!!$            call terminate("setupADjointMatrix","Error in MatAssemblyBegin")
!!$       
!!$
!!$       call MatAssemblyEnd(dRdxFD,MAT_FINAL_ASSEMBLY,PETScIerr)
!!$       
!!$       if( PETScIerr/=0 ) &
!!$            call terminate("setupADjointMatrix","Error in MatAssemblyEnd")
!!$
!!$       !if( debug ) then
!!$          !call MatView(dRdWFD,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
!!$          !call MatView(dRdWFD,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
!!$          !if( PETScIerr/=0 ) &
!!$          !     call terminate("setupADjointMatrix", "Error in MatView")
!!$          !pause
!!$       !endif
       call VecAssemblyBegin(pvr,PETScIerr)
       
       if( PETScIerr/=0 ) &
            call terminate("setupASjointRHS", "Error in VecAssemblyBegin")  
       
       call VecAssemblyEnd  (pvr,PETScIerr)
       
       if( PETScIerr/=0 ) &
            call terminate("setupADjointRHS", "Error in VecAssemblyEnd")
       
       if( debug ) then
          call VecView(pvr,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
          !call VecView(pvr,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
          if( PETScIerr/=0 ) &
               call terminate("setupADjointRHS", "Error in VecView")
          pause
       endif
       
       call cpu_time(time(4))
       timeFD = time(4)-time(3)
       
       print *,'finished first FD calculation.....'

!     ******************************************************************
!     *                                                                *
!     *  Compute d(dw)/d(w) using central finite differences           *
!     *                                                                *
!     ******************************************************************

      ! Get the initial FD time.
      print *,'starting 2nd FD calculation'
      call mpi_barrier(SUmb_comm_world, ierr)
      if(myID == 0) call cpu_time(time(3))


      deltax = 1.d-7
      xFD(:,:,:,:) = 0. 
      call cpu_time(time(3))

!!$      !print *,'entering FD loop'
!!$      domainResidualLoopFDorig: do nn=1,nDom         
!!$
!!$         ! Determine the number of time instances for this block
!!$        
!!$         nTime     = nTimeIntervalsSpectral
!!$
!!$         ! Loop over the number of time instances for this block.
!!$
!!$         spectralLoop2: do sps=1,nTime
!!$!            print *,'setting Pointers'
!!$            call setPointers(nn,level,sps)
!!$
!!$
!!$               ! Loop over location of output (R) cell of residual
!!$               do kCell = 2, kl
!!$                  print *, "kCell =", kCell
!!$                  do jCell = 2, jl
!!$                     do iCell = 2, il
!!$                        
!!$                        do m = 1, nw           ! Loop over output cell residuals (R)
!!$!                           print *,'indices',icell,jcell,kcell,m
!!$                        ! Copy the state in the stencil
!!$                        ! actually no need to call again, but ...
!!$!                           call copyADjointStencil(wAdj,xAdj, iCell, jCell, kCell)
!!$                           call copyADjointStencil(wAdj, xAdj,alphaAdj,betaAdj,&
!!$                          MachAdj,MachCoefAdj,iCell, jCell, kCell,prefAdj,&
!!$                          rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
!!$                          rhoinfAdj, pinfAdj,&
!!$                          murefAdj, timerefAdj,pInfCorrAdj,liftIndex)
!!$!                           print *,'stencil copied'
!!$                           ! Loop over the inputs (w) 
!!$                           !print *,'secondhalo',secondhalo
!!$                        
!!$                           do ii = -3, 2
!!$                              do jj = -3, 2
!!$                                 do kk = -3, 2
!!$                                    do n = 1, 3
!!$                                       i = iCell+ii
!!$                                       j = jCell+jj
!!$                                       k = kCell+kk
!!$ !                                      print *,'secondary FD indices',i,j,k,n,ii,jj,kk
!!$!                                       if (ii==0 .or. jj==0 .or. kk==0) then
!!$                                       xAdjRef = xAdj(ii,jj,kk,n)
!!$                                       xAdj(ii,jj,kk,n) = xAdjRef + deltax
!!$                                       call computeRAdjoint(wAdj,xAdj,dwAdjP,alphaAdj,&
!!$                                            betaAdj,MachAdj, MachCoefAdj,&
!!$                                            iCell, jCell,  kCell, &
!!$                                            nn,sps, correctForK,secondHalo,prefAdj,&
!!$                                            rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
!!$                                            rhoinfAdj, pinfAdj,&
!!$                                            murefAdj, timerefAdj,pInfCorrAdj,liftIndex)
!!$                                          !call computeRAdjoint(wAdj,xAdj,dwAdjP,   &
!!$                                          !     iCell, jCell,  kCell, &
!!$                                          !     nn,sps, correctForK,secondHalo)
!!$                                          
!!$                                       xAdj(ii,jj,kk,n) = xAdjRef - deltax
!!$                                       call computeRAdjoint(wAdj,xAdj,dwAdjM,alphaAdj,&
!!$                                            betaAdj,MachAdj, MachCoefAdj,&
!!$                                            iCell, jCell,  kCell, &
!!$                                            nn,sps, correctForK,secondHalo,prefAdj,&
!!$                                            rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
!!$                                            rhoinfAdj, pinfAdj,&
!!$                                            murefAdj, timerefAdj,pInfCorrAdj,liftIndex)
!!$                                          !call computeRAdjoint(wAdj,xAdj,dwAdjM,   &
!!$                                          !     iCell, jCell,  kCell, &
!!$                                          !     nn,sps, correctForK,secondHalo)
!!$
!!$                                          
!!$                                       xFD(ii,jj,kk,n) = (dwAdjP(m)-dwAdjM(m))/(2.0*deltax)
!!$                                       xAdj(ii,jj,kk,n) = xAdjRef
!!$
!!$                                       if(i>zero .and. j>zero .and. k>zero .and. i<=ie .and. j<=je .and. k<=ke)then
!!$                                          !print *,'global FD indicies',idxstate,idxres
!!$                                          !print *,'indicies',i,j,k,n,shape(globalNode)
!!$                                          idxnode = globalNode(i,j,k)*3+n
!!$                                          !print *,'globalNode',globalNode(i,j,k),shape( dRdxFD2),shape(xFD)
!!$                                          idxres   = globalCell(iCell,jCell,kCell)*nw+m
!!$                                          if( idxres>=0 .and. idxnode>=0) then
!!$                                             !print *,'passed indicies',idxres,idxnode
!!$                                             !dRdxFD2(idxres,idxnode,nn,sps) = xFD(ii,jj,kk,n)
!!$                                             dRdxFD2(idxres,idxnode,1,1) = xFD(ii,jj,kk,n)
!!$                                          endif
!!$                                       end if
!!$!                                    end if
!!$                                    end do ! n
!!$                                 end do ! ii
!!$                              end do !jj
!!$                           end do !kk
!!$                        
!!$                        !print 111, iRes, iCell-1, jCell-1, kCell-1, dRdw(iCell-1, jCell-1, kCell-1, iRes), dRdwFD(iCell-1, jCell-1, kCell-1, iRes)
!!$                     end do !iRes
!!$                  end do !iCell
!!$               end do !jCell
!!$            end do !kCell
!!$         end do spectralLoop2
!!$      end do domainResidualLoopFDorig
!!$      
      call cpu_time(time(4))
      timeFD = time(4)-time(3)

      print *,'finished 2nd FD calculation'
!_______________________________________________________
!
!     Compute the errors in dR/dw
!_______________________________________________________

      print *, "Computing the the error..."
 outputDomainLoop: do nn = 1,1! nDom    
    call setpointers(nn,1,1)
    do m=1,nw
      do iCell=2,il
        do jCell=2,jl
          do kCell=2,kl
            do ii=-3,2
              do jj=-3,2
                do kk=-3,2
                 do n=1,3-1
                  i=iCell+ii; j=jCell+jj; k=kCell+kk
                     if(i>zero .and. j>zero .and. k>zero .and. i<=ie .and. j<=je .and. k<=ke)then
                        idxnode = globalNode(i,j,k)*3+n
                        idxres   = globalCell(iCell,jCell,kCell)*nw+m
                       
                        if( idxres>=0 .and. idxnode>=0) then
                           !write(*,*) iCell,jCell,kCell,m,ii,jj,kk,n,idxnode,idxres
                           dRdxErr(idxres, idxnode, 1, 1) = dRdxAdj(idxres, idxnode, 1, 1) - dRdxFD1(idxres, idxnode, 1, 1)
                           !dRdxErr(idxres, idxnode, 1, 1) = dRdxFD2(idxres, idxnode, 1, 1) - dRdxFD1(idxres, idxnode, 1, 1)
!                           if(dRdxFD(idxres,idxnode,1,1).ne.0 .and. dRdxErr(idxres, idxnode, 1, 1)>1e-10.and.dRdxAdj(idxres, idxnode, 1, 1)==0) then
!                           if(dRdxFD1(idxres,idxnode,1,1).ne.0 .and. dRdxErr(idxres, idxnode, 1, 1)>1e-10) then
                           if(dRdxErr(idxres, idxnode, 1, 1)>1e-9) then
                           !if(dRdxFD1(idxres,idxnode,1,1).ne.0 .and. dRdxErr(idxres, idxnode, 1, 1).ne.0) then
                           !if(dRdxFD1(idxres,idxnode,1,1).ne.0 .and. dRdxAdj(idxres, idxnode, 1, 1)==0.0) then
                           !if(dRdxAdj(idxres,idxnode,1,1).ne.0) then
                           !if(dRdxFD1(idxres,idxnode,1,1).ne.0) then
                           !if(dRdxFD(idxres,idxnode,1,1)>1e-12) then
                              !if((ii.ne.zero .and.jj.ne.zero).and.(ii.ne.zero .and.kk.ne.zero).and.(kk.ne.zero .and.jj.ne.zero))then
                              !print *,'test',(ii.ne.zero .and.jj.ne.zero).and.(ii.ne.zero .and.kk.ne.zero).and.(kk.ne.zero .and.jj.ne.zero),(ii.ne.zero .and.jj.ne.zero),(ii.ne.zero .and.kk.ne.zero),(kk.ne.zero .and.jj.ne.zero)
                                 write(*,*) nn,iCell,jCell,kCell,m,ii,jj,kk,n, idxres,idxnode,&
                                      dRdxAdj(idxres,idxnode,1,1), dRdxFD1(idxres,idxnode,1,1), dRdxErr(idxres,idxnode,1,1),dRdxFD2(idxres,idxnode,1,1)
                              !endif
                              !endif
                           end if
                        end if
                     endif
                end do
              end do
            end do
          end do
        end do
      end do
   end do
end do
enddo outputDomainLoop

   deallocatedomains: do nn = 1,ndom
         print *,'domain',nn
         groundLevel = 1
         sps = 1
         call setPointers(nn,1,sps)
         deallocate(flowDoms(nn,level,sps)%dwp)
         deallocate(flowDoms(nn,level,sps)%dwm)
         deallocate(flowDoms(nn,level,sps)%dwtemp)
      end do deallocatedomains

!      dRdwErrRel_1(:, :, :, :) = abs(dRdwErr_q(:, :, :, :)/dRdw_q(:,:,:,:))
!      print *
      print *, "max, min of dRdxErr =", maxval(dRdxErr(:,:,1,1)), minval(dRdxErr(:,:,1,1))
!      print *, "max, min of dRdxErrRel =", maxval(dRdxErrRel_q(:,:,:,1)), minval(dRdxErrRel_q(:,:,:,1))
!      print *

99    format(a,1x,i6)
 111  format(4I4, 3ES22.14)


#endif
  end subroutine verifydRdx
