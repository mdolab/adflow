!
!     ******************************************************************
!     *                                                                *
!     * File:          verifydRdw.f90                                  *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 05-05-2008                                      *
!     * Last modified: 05-06-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine verifydRdw(level)
!
!     ******************************************************************
!     *                                                                *
!     *  This subroutine computes the values for dRdw and compares     *
!     *  the Tapenade result to the finite-difference results.         *
!     *                                                                *
!     ******************************************************************


      use blockPointers ! block (nDoms,flowDoms), globalCell
      use flowvarrefstate
      use communication
      use iteration     ! groundLevel
      use inputTimeSpectral ! spaceDiscr
      use inputIO

      !from old verify routine
      use ADjointPETSc, only: drdw,drdwfd,petscone,insert_values,petscierr,mat_final_assembly,petsc_viewer_draw_world,petsc_viewer_stdout_world
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
      integer(kind=intType) :: iCell,jCell,kCell,liftIndex
      real(kind=realType), dimension(nw) :: dwL2
      real(kind=realType), dimension(nx, ny, nz, nw) :: dwerr
      real(kind=realType), dimension(10) :: time
      real(kind=realType) :: timeRes
!      real(kind=realType), dimension(4) :: time
      real(kind=realType) ::  timeOri
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdj
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdjb
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wFD
!!$      real(kind=realType), dimension(-2:3,-2:3,-2:3,3)  :: xAdj
!!$      real(kind=realType), dimension(-2:3,-2:3,-2:3,3)  :: xAdjb
!!$      real(kind=realType), dimension(-2:3,-2:3,-2:3,3)  :: xFD
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

      character fileName*32, dataName*32
      real(kind=realType), dimension(nw) :: dwAdj,dwAdjb,dwAdjRef
      real(kind=realType), dimension(nw) :: dwAdjP, dwAdjM
      real(kind=realType) :: deltaw, wAdjRef
      real(kind=realType) :: timeAdj, timeFD, timeResAdj

      integer(kind=intType), dimension(nDom) :: maxglobalcell
      integer(kind=intType) :: idx, ii, jj, kk, idxstate, idxres, m, l
      integer(kind=intType) :: sps, nTime, max_nTime, nHalo, nn, discr
      real(kind=realType), allocatable, dimension(:,:,:,:) :: dRdwAdj,dRdwFD1,dRdwFD2,dRdwErr

      REAL(KIND=REALTYPE), DIMENSION(3) :: rotcenteradj
      REAL(KIND=REALTYPE), DIMENSION(3):: rotrateadj
      REAL(KIND=REALTYPE) :: rotrateadjb(3)

      integer :: ierr
      logical :: fineGrid, correctForK, exchangeTurb,secondHalo

      !FULL FD Variables
      
      !real(kind=realType), dimension(nx, ny, nz, nw,nx) :: dRdwErr, dRdwErrRel
      real(kind=realType), dimension(nw) :: dRdwL2, dRdwL2Rel
      real(kind=realtype), dimension(0:ib,0:jb,0:kb,1:nw)::dwp,dwm,dwtemp
      real(kind=realType), dimension(0:ib,0:jb,0:kb,1:nw) :: wtemp
      real(kind=realType), dimension(0:ib,0:jb,0:kb) :: ptemp
      integer(kind=intType) :: istate, jstate, kstate,ires
      real(kind=realType), dimension(-2:2,-2:2,-2:2) :: pAdjtemp
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdjtemp
      real(kind=realType), dimension(nx,ny,nz,nw) :: wFD2
      !real(kind=realType), dimension(ib*jb*kb*nw,ib*jb*kb*nw) :: dRdw
      !real(kind=realType), dimension(ib*jb*kb*nw,ib*jb*kb*nw) :: dRdwFD

!     ******************************************************************
!     *                                                                *
!     *  Begin execution                                               *
!     *                                                                *
!     ******************************************************************

      if( myID==0 ) write(*,*) "Running verifydRdW..."

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
         nTime     = nTimeIntervalsSpectral!sections(sectionID)%nTimeInstances
         if(nTime>=max_nTime) max_nTime = nTime
      enddo
      idx = maxval(maxglobalcell(:))
      print *,'allocating',idx,nw*(idx+1),nw*(idx+1),ndom,max_nTime
      
      allocate(dRdwErr(nw*(idx+1),nw*(idx+1),ndom,max_nTime), &
               dRdwAdj(nw*(idx+1),nw*(idx+1),ndom,max_nTime), &
               dRdwFD1(nw*(idx+1),nw*(idx+1),ndom,max_nTime), &
               dRdwFD2(nw*(idx+1),nw*(idx+1),ndom,max_nTime))
      print *,' allocated'

!           if(ierr /= 0)                       &
!                call terminate("memory?") 
      
      ! Initialize the temporary arrays.
      dRdwErr = 0
      dRdwAdj = 0
      dRdwFD1  = 0


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
            call setPointersAdj(nn,level,sps)

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
                        do ii=-2,2!1,il-1
                           do jj = -2,2!1,jl-1
                              do kk = -2,2!1,kl-1
                                 do l = 1,nw
                                    i = iCell + ii
                                    j = jCell + jj
                                    k = kCell + kk
                                    !print *,'secondaryindicies',i,j,k,ii,jj,kk
                                    if(i>zero .and. j>zero .and. k>zero .and. i<=il .and. j<=jl .and. k<=kl)then
                                       idxstate = globalCell(i,j,k)*nw+l
                                       idxres   = globalCell(iCell,jCell,kCell)*nw+m
                                       !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
                                       if( idxres>=0 .and. idxstate>=0) then
                                          dRdwAdj(idxres,idxstate,nn,sps) = wAdjb(ii,jj,kk,l)
!!$                                          if (icell==13.and.jcell==5.and.kcell==5.and.ii==0.and.jj==-1.and.kk==0)then
!!$                                             print *,'wAdjb', wAdjb(ii,jj,kk,l),icell,jcell,kcell,m,ii,jj,kk,l,i,j,k
!!$                                             stop
!!$                                          endif
                                          
                                           !if ( wAdjb(ii,jj,kk,l).ne. 0.0) then
                                           !   print *,'wAdjb', wAdjb(ii,jj,kk,l),ii,jj,kk,l
!!$                                              if (i==12) then
!!$                                                 print *,'wAdjb', wAdjb(ii,jj,kk,l),icell,jcell,kcell,m,ii,jj,kk,l,i,j,k
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
!     Compute d(dw)/d(w) using central finite differences
!_______________________________________________________
       
       deltaw = 1.d-5
       print *, "deltaw=", deltaw
       wFD2(:,:,:,:) = 0. 
       call cpu_time(time(3))
       
       groundLevel = 1
       sps = 1
       nn=1
       
       ! Loop over the inputs (w) 
       do n = 1, nw
          do kstate = 2, kl
             print *, "iw, kstate =", n, kstate
             do jstate = 2, jl
                do istate = 2, il
                   
                   ! if (istate==icell .or. jstate==jcell .or. kstate==kcell) then
                   !Remember current values
                   dwtemp(:,:,:,:) = dw(:,:,:,:)
                   wtemp(:,:,:,:) = w(:,:,:,:)
                   ptemp(:,:,:) = p(:,:,:)
                   wAdjRef = w(istate,jstate,kstate,n)
                   
                   w(istate,jstate,kstate,n) = wAdjRef + deltaw
                   wAdjtemp = w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:)
!                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
                   call computePressureAdj(wAdjtemp, pAdjtemp)
                   p(istate, jstate, kstate) = pAdjtemp(0,0,0)
                   
                   !print *, "Calling Residual ="
                   call initres(1_intType, nwf)
                   
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
                   !print *, "Called Residual =", dw(istate,jstate,kstate,n)
                   
                   dwp(:,:,:,:) = dw(:,:,:,:)
                   dw(:,:,:,:) = dwtemp(:,:,:,:)
                   
                   w(istate,jstate,kstate,n) = wAdjRef - deltaw
                   wAdjtemp = w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:)
!                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
                   call computePressureAdj(wAdjtemp, pAdjtemp)
!                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
                   p(istate, jstate, kstate) = pAdjtemp(0,0,0)
                   
                   !print *, "Calling Residual ="
                   call initres(1_intType, nwf)
                   
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
                   !print *, "Called Residual 2 =", dw(istate,jstate,kstate,n)
                   dwm(:,:,:,:) = dw(:,:,:,:)
                   dw(:,:,:,:) = dwtemp(:,:,:,:)
                   
                   !                             end if

                   
!                   ival = 0
                   do iRes = 1, nw           ! Loop over output cell residuals (R)
                      ! Loop over location of output (R) cell of residual
                      do kCell = 2, kl
                         
                         do jCell = 2, jl
                            do iCell = 2, il
                               
                               !column = (n-1) +(istate-2)*nw +(jstate-2)*nw*nx +(kstate-2)*nw*nx*ny
                               wFD2(iCell-1, jCell-1, kCell-1, iRes) = (dwp(iCell,jCell,kCell,iRes)-dwm(iCell,jCell,kCell,iRes))/(2.0*deltaw)
                               
                               !if (abs(wFD(iCell-1, jCell-1, kCell-1, iRes)) >1e-12) then
                               !   ival = ival + 1
                               !   row2(ival) = (ires-1) + (icell-2)*nw + (jcell-2)*nw*nx + (kcell-2)*nw*nx*ny
                               !   col_values(ival)  = wFD(iCell-1, jCell-1, kCell-1, iRes)
                               !   
                                !  
                                !  
                               !end if
                               !if(i>zero .and. j>zero .and. k>zero .and. i<il .and. j<jl .and. k<kl)then
                                  idxstate = globalCell(istate,jstate,kstate)*nw+n
                                  idxres   = globalCell(iCell,jCell,kCell)*nw+ires
                                  !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
                                  if( idxres>=0 .and. idxstate>=0) then
                                     if (wFD2(iCell-1, jCell-1, kCell-1, iRes).ne. zero) then
                                        dRdwFD1(idxres,idxstate,nn,sps) = wFD2(iCell-1, jCell-1, kCell-1, iRes)
                                        !Aad(m,:)  = wFD2(iCell-1, jCell-1, kCell-1, iRes)wAdjB( 0, 0, 0,:)
                                        !print *,'setting',dRdW, PETScOne, idxres, PETScOne, idxstate,   &
                                        !     dRdwFD(idxres,idxstate,nn,sps)
                                        call MatSetValues(dRdWFD, 1, idxres-1, 1, idxstate-1,   &
                                             dRdwFD1(idxres,idxstate,nn,sps), INSERT_VALUES, PETScIerr)
                                        if( PETScIerr/=0 ) &
                                             print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")


!!$                                        call MatSetValues(dRdW, 1, idxres-1, 1, idxstate-1,   &
!!$                                             dRdwFD(idxres,idxstate,nn,sps), INSERT_VALUES, PETScIerr)
!!$                                        if( PETScIerr/=0 ) &
!!$                                             print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
                                        
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
                                     endif
                                  endif
                               !end if
                               
                            end do
                         end do
                      end do
                   end do
                   w(istate,jstate,kstate,n) = wAdjRef
                   dw(:,:,:,:) = dwtemp(:,:,:,:)
                   w(:,:,:,:) = wtemp(:,:,:,:)
                   p(:,:,:) = ptemp(:,:,:)
!!$                   
!!$                   col_size = ival
!!$                   
!!$                   call MatSetValues(dRdWFD, col_size, row2, PetscOne, column, col_values, INSERT_VALUES, PETScIerr)
                end do
             end do
          end do
       end do

       call MatAssemblyBegin(dRdWFD,MAT_FINAL_ASSEMBLY,PETScIerr)

       if( PETScIerr/=0 ) &
            call terminate("setupADjointMatrix","Error in MatAssemblyBegin")
       

       call MatAssemblyEnd(dRdWFD,MAT_FINAL_ASSEMBLY,PETScIerr)
       
       if( PETScIerr/=0 ) &
            call terminate("setupADjointMatrix","Error in MatAssemblyEnd")

       !if( debug ) then
          !call MatView(dRdWFD,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
          !call MatView(dRdWFD,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
          !if( PETScIerr/=0 ) &
          !     call terminate("setupADjointMatrix", "Error in MatView")
          !pause
       !endif
       
       call cpu_time(time(4))
       timeFD = time(4)-time(3)

!     ******************************************************************
!     *                                                                *
!     *  Compute d(dw)/d(w) using central finite differences           *
!     *                                                                *
!     ******************************************************************

      ! Get the initial FD time.
      !print *,'starting FD calculation'
      call mpi_barrier(SUmb_comm_world, ierr)
      if(myID == 0) call cpu_time(time(3))


      deltaw = 1.d-6
      wFD(:,:,:,:) = 0. 
      call cpu_time(time(3))

      !print *,'entering FD loop'
      domainResidualLoopFDorig: do nn=1,nDom         

         ! Determine the number of time instances for this block
        
         nTime     = nTimeIntervalsSpectral

         ! Loop over the number of time instances for this block.

         spectralLoop2: do sps=1,nTime
!            print *,'setting Pointers'
            call setPointersAdj(nn,level,sps)


               ! Loop over location of output (R) cell of residual
               do kCell = 2, kl
                  !print *, "iRes, kCell =", iRes, kCell
                  do jCell = 2, jl
                     do iCell = 2, il
                        
                        do m = 1, nw           ! Loop over output cell residuals (R)
!                           print *,'indices',icell,jcell,kcell,m
                        ! Copy the state in the stencil
                        ! actually no need to call again, but ...
!                           call copyADjointStencil(wAdj,xAdj, iCell, jCell, kCell)
                           call copyADjointStencil(wAdj, xAdj,alphaAdj,betaAdj,MachAdj,&
                                machCoefAdj,iCell, jCell, kCell,prefAdj,&
                                rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                                rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
                                murefAdj, timerefAdj,pInfCorrAdj,liftIndex)
!                           copyADjointStencil(wAdj, xAdj,alphaAdj,betaAdj,&
!                                MachAdj,MachCoefAdj,iCell, jCell, kCell,prefAdj,&
!                                rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
!                                rhoinfAdj, pinfAdj,&
!                                murefAdj, timerefAdj,pInfCorrAdj,liftindex)
!                           print *,'stencil copied'
                           ! Loop over the inputs (w) 
                           !print *,'secondhalo',secondhalo
                        
                           do ii = -2, 2
                              do jj = -2, 2
                                 do kk = -2, 2
                                    do n = 1, nw
                                       i = iCell+ii
                                       j = jCell+jj
                                       k = kCell+kk
 !                                      print *,'secondary FD indices',i,j,k,n,ii,jj,kk
                                       if (ii==0 .or. jj==0 .or. kk==0) then
                                          wAdjRef = wAdj(ii,jj,kk,n)
                                          wAdj(ii,jj,kk,n) = wAdjRef + deltaw
                                          call computeRAdjoint(wAdj,xAdj,dwAdj,alphaAdj,betaAdj,MachAdj, &
                                               MachCoefAdj,iCell, jCell,  kCell, &
                                               nn,sps, correctForK,secondHalo,prefAdj,&
                                               rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                                               rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
                                               murefAdj, timerefAdj,pInfCorrAdj,liftIndex)

                                          
                                          wAdj(ii,jj,kk,n) = wAdjRef - deltaw
                                          call computeRAdjoint(wAdj,xAdj,dwAdj,alphaAdj,betaAdj,MachAdj, &
                                               MachCoefAdj,iCell, jCell,  kCell, &
                                               nn,sps, correctForK,secondHalo,prefAdj,&
                                               rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                                               rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
                                               murefAdj, timerefAdj,pInfCorrAdj,liftIndex)

                                          !call computeRAdjoint(wAdj,xAdj,dwAdjM,   &
                                          !     iCell, jCell,  kCell, &
                                          !     nn,sps, correctForK,secondHalo)

                                          
                                          wFD(ii,jj,kk,n) = (dwAdjP(m)-dwAdjM(m))/(2.0*deltaw)
                                          wAdj(ii,jj,kk,n) = wAdjRef

                                          if(i>zero .and. j>zero .and. k>zero .and. i<=il .and. j<=jl .and. k<=kl)then
!                                             print *,'global FD indicies',idxstate,idxres
                                             idxstate = globalCell(i,j,k)*nw+n
                                             idxres   = globalCell(iCell,jCell,kCell)*nw+m
                                             if( idxres>=0 .and. idxstate>=0) then
                                                dRdwFD2(idxres,idxstate,nn,sps) = wFD(ii,jj,kk,n)
                                             endif
                                          end if
                                       end if
                                    end do ! n
                                 end do ! ii
                              end do !jj
                           end do !kk
                        
                        !print 111, iRes, iCell-1, jCell-1, kCell-1, dRdw(iCell-1, jCell-1, kCell-1, iRes), dRdwFD(iCell-1, jCell-1, kCell-1, iRes)
                     end do !iRes
                  end do !iCell
               end do !jCell
            end do !kCell
         end do spectralLoop2
      end do domainResidualLoopFDorig
      
      call cpu_time(time(4))
      timeFD = time(4)-time(3)

!_______________________________________________________
!
!     Compute the errors in dR/dw
!_______________________________________________________

      !print *, "Computing the the error..."
    
    do m=1,nw
      do iCell=2,il
        do jCell=2,jl
          do kCell=2,kl
            do ii=-2,2
              do jj=-2,2
                do kk=-2,2
                 do n=1,nw
                  i=iCell+ii; j=jCell+jj; k=kCell+kk
                     if(i>zero .and. j>zero .and. k>zero .and. i<=il .and. j<=jl .and. k<=kl)then
                        idxstate = globalCell(i,j,k)*nw+n
                        idxres   = globalCell(iCell,jCell,kCell)*nw+m
                        !write(*,*) iCell,jCell,kCell,m,ii,jj,kk,n,idxstate,idxres
                        if( idxres>=0 .and. idxstate>=0) then
                           dRdwErr(idxres, idxstate, 1, 1) = dRdwAdj(idxres, idxstate, 1, 1) - dRdwFD1(idxres, idxstate, 1, 1)
!                           if(dRdwFD(idxres,idxstate,1,1).ne.0 .and. dRdwErr(idxres, idxstate, 1, 1)>1e-10.and.dRdwAdj(idxres, idxstate, 1, 1)==0) then
                           if(dRdwFD1(idxres,idxstate,1,1).ne.0 .and. dRdwErr(idxres, idxstate, 1, 1)>1e-8) then
!                           if(dRdwFD1(idxres,idxstate,1,1).ne.0 .and. dRdwAdj(idxres, idxstate, 1, 1)==0.0) then
                           !if(dRdwFD(idxres,idxstate,1,1).ne.0) then
                           !if(dRdwFD(idxres,idxstate,1,1)>1e-12) then
                              !if((ii.ne.zero .and.jj.ne.zero).and.(ii.ne.zero .and.kk.ne.zero).and.(kk.ne.zero .and.jj.ne.zero))then
                              !print *,'test',(ii.ne.zero .and.jj.ne.zero).and.(ii.ne.zero .and.kk.ne.zero).and.(kk.ne.zero .and.jj.ne.zero),(ii.ne.zero .and.jj.ne.zero),(ii.ne.zero .and.kk.ne.zero),(kk.ne.zero .and.jj.ne.zero)
                                 write(*,*) iCell,jCell,kCell,ii,jj,kk, &
                                      dRdWAdj(idxres,idxstate,1,1), dRdwFD1(idxres,idxstate,1,1), dRdwErr(idxres,idxstate,1,1),dRdwFD2(idxres,idxstate,1,1)
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


!      dRdwErrRel_1(:, :, :, :) = abs(dRdwErr_q(:, :, :, :)/dRdw_q(:,:,:,:))
!      print *
      !print *, "max, min of dRdwErr =", maxval(dRdwErr(:,:,1,1)), minval(dRdwErr(:,:,1,1))
!      print *, "max, min of dRdwErrRel =", maxval(dRdwErrRel_q(:,:,:,1)), minval(dRdwErrRel_q(:,:,:,1))
!      print *


 111  format(4I4, 3ES22.14)



    end subroutine verifyDRdw
