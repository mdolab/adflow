!
!     ******************************************************************
!     *                                                                *
!     * File:          verifydRdM.f90                                  *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 06-03-2008                                      *
!     * Last modified: 06-03-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine verifydRdExtra(level)
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
      use ADjointPETSc, only: pvr,drda,petscone,insert_values,petscierr,mat_final_assembly,petsc_viewer_draw_world
      use ADjointVars
      !use FDPETSc, only: DRDWFD
      use precision
      !use blockPointers
      !use flowvarrefstate
      !use iteration
      use inputIteration
      use inputPhysics  !mach,machcoef
      use section
      use cgnsgrid
      use monitor
      !implicit none
     

      implicit none

!
!     Subroutine arguments
      integer(kind=intType), intent(in) :: level


!
!     Local variables 
!
      integer(kind=intType) :: i, j, k, n,mm
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
      REAL(KIND=REALTYPE) :: xblockcorneradj(2, 2, 2, 3), xblockcorneradjb(2&
&  , 2, 2, 3)
      REAL(KIND=REALTYPE) :: machadj, machcoefadj,machGridAdj, pinfcorradj
      REAL(KIND=REALTYPE) :: machadjb, machcoefadjb,machgridadjb, pinfcorradjb
      REAL(KIND=REALTYPE) :: prefadj, rhorefadj
      REAL(KIND=REALTYPE) :: pinfdimadj, rhoinfdimadj
      REAL(KIND=REALTYPE) :: rhoinfadj, pinfadj
      REAL(KIND=REALTYPE) :: rhoinfadjb
      REAL(KIND=REALTYPE) :: murefadj, timerefadj
      REAL(KIND=REALTYPE) :: alphaadj, betaadj
      REAL(KIND=REALTYPE) :: alphaadjb, betaadjb
      real(kind=realType), dimension(3) ::rotRateAdj,rotCenterAdj,rotrateadjb
      real(kind=realType), dimension(3) :: velDirFreestreamAdj
      real(kind=realType), dimension(3) :: liftDirectionAdj
      real(kind=realType), dimension(3) :: dragDirectionAdj

      character fileName*32, dataName*32
      real(kind=realType), dimension(nw) :: dwAdj,dwAdjb,dwAdjRef
      real(kind=realType), dimension(nw) :: dwAdjP, dwAdjM
      real(kind=realType) :: deltaw, ExtraAdjRef
      real(kind=realType) :: timeAdj, timeFD, timeResAdj,test

      integer(kind=intType), dimension(nDom) :: maxglobalcell
      integer(kind=intType) :: idx, ii, jj, kk, idxres, m, l,idxmgb,liftindex
      integer(kind=intType) :: sps, nTime, max_nTime, nHalo, nn, discr
      real(kind=realType), allocatable, dimension(:,:,:,:) :: dRdExtraAdj,dRdExtraFD,dRdExtraErr

      integer :: ierr
      logical :: fineGrid, correctForK, exchangeTurb,secondHalo

      !FULL FD Variables
      
      !real(kind=realType), dimension(nx, ny, nz, nw,nx) :: dRdwErr, dRdwErrRel
      real(kind=realType), dimension(nw) :: dRdwL2, dRdwL2Rel
      !real(kind=realtype), dimension(0:ib,0:jb,0:kb,1:nw)::dwp,dwm,dwtemp
      real(kind=realType), dimension(0:ib,0:jb,0:kb,1:nw) :: wtemp
      real(kind=realType), dimension(0:ib,0:jb,0:kb) :: ptemp
!      integer(kind=intType) :: istate, jstate, kstate,ires
      integer(kind=intType) :: ires
      real(kind=realType), dimension(-2:2,-2:2,-2:2) :: pAdjtemp
      real(kind=realType), dimension(-2:2,-2:2,-2:2,nw) :: wAdjtemp
      real(kind=realType), dimension(nx,ny,nz,nw) :: wFD2
      !real(kind=realType), dimension(ib*jb*kb*nw,ib*jb*kb*nw) :: dRdw
      !real(kind=realType), dimension(ib*jb*kb*nw,ib*jb*kb*nw) :: dRdwFD

      real(kind=realType), dimension(nSections) :: t
      
      ! pvr row block
      
      real(kind=realType), dimension(nw) :: pvrlocal
      character(len=2*maxStringLen) :: errorMessage
     !File Parameters
      integer :: unitM = 8,unitAoA = 30,unitSSA = 12,unitRotx = 13,ierror,unitRoty,unitRotz
      character(len = 16)::outfile,testfile

      write(testfile,100) myid!12
100   format (i5)  
      testfile=adjustl(testfile)
      write(outfile,101) trim(testfile)!testfile
101   format("CSMachfile",a,".out")
      !outfile = "CSMachfile.txt"
      unitM = 8+myID
      !outfile = "CSMachfile.txt"
      
      open (UNIT=unitM,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifydRdExtraCS", &
           "Something wrong when &
           &calling open")
      write(testfile,102) myid!12
102   format (i5)  
      testfile=adjustl(testfile)
      write(outfile,103) trim(testfile)!testfile
103   format("CSAoAfile",a,".out")
      unitAoA = 8+myID+nproc*1
      !outfile = "CSAOAfile.txt"
      
      open (UNIT=unitAoA,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifydRdExtraCS", &
           "Something wrong when &
           &calling open")
      write(testfile,104) myid!12
104   format (i5)  
      testfile=adjustl(testfile)
      write(outfile,105) trim(testfile)!testfile
105   format("CSSSAfile",a,".out")
      unitSSA = 8+myID+nproc*2
      !outfile = "CSSSAfile.txt"
      
      open (UNIT=unitSSA,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifydRdExtraCS", &
           "Something wrong when &
           &calling open")

      write(testfile,106) myid!12
106   format (i5)  
      testfile=adjustl(testfile)
      write(outfile,107) trim(testfile)!testfile
107   format("CSRotxfile",a,".out")
      unitrotx = 8+myID+nproc*3
      !outfile = "CSRotxfile.txt"
      
      open (UNIT=unitrotx,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifydRdExtraCS", &
           "Something wrong when &
           &calling open")

      write(testfile,108) myid!12
108   format (i5)  
      testfile=adjustl(testfile)
      write(outfile,109) trim(testfile)!testfile
109   format("CSRotyfile",a,".out")
      unitroty = 8+myID+nproc*4
      !outfile = "CSRotxfile.txt"
      
      open (UNIT=unitroty,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifydRdExtraCS", &
           "Something wrong when &
           &calling open")

      write(testfile,110) myid!12
110   format (i5)  
      testfile=adjustl(testfile)
      write(outfile,111) trim(testfile)!testfile
111   format("CSRotzfile",a,".out")
      unitrotz = 8+myID+nproc*5
      !outfile = "CSRotxfile.txt"
      
      open (UNIT=unitrotz,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifydRdExtraCS", &
           "Something wrong when &
           &calling open")


!     ******************************************************************
!     *                                                                *
!     *  Begin execution                                               *
!     *                                                                *
!     ******************************************************************

      if( myID==0 ) write(*,*) "Running verifydRdExtra..."

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
         !print *,'allocating i'
         maxglobalcell(i) = maxval(flowDoms(i,currentLevel,1)%globalCell(:,:,:))
         !print *,'maxglobalcell',maxglobalcell(i) 
         nTime     = nTimeIntervalsSpectral!sections(sectionID)%nTimeInstances
         if(nTime>=max_nTime) max_nTime = nTime
      enddo
      idx = maxval(maxglobalcell(:))
      !print *,'allocating',idx,nw*(idx+1),nw*(idx+1),ndom,max_nTime
      !allocate(dRdExtraErr(nw*(idx+1),nDesignExtra,ndom,max_nTime), &
      !         dRdExtraAdj(nw*(idx+1),nDesignExtra,ndom,max_nTime), &
      !         dRdExtraFD(nw*(idx+1),nDesignExtra,ndom,max_nTime))
      allocate(dRdExtraErr(nw*(idx+1),nDesignExtra,1,max_nTime), &
           dRdExtraAdj(nw*(idx+1),nDesignExtra,1,max_nTime), &
           dRdExtraFD(nw*(idx+1),nDesignExtra,1,max_nTime))
      !print *,' allocated'

    !allocate memory for FD
      allocatedomains: do nn = 1,ndom
         print *,'domain',nn
         groundLevel = 1
         sps = 1
         call setPointersAdj(nn,1,sps)
         allocate(flowDoms(nn,level,sps)%dwp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
         allocate(flowDoms(nn,level,sps)%dwm(0:ib,0:jb,0:kb,1:nw),stat=ierr)
         allocate(flowDoms(nn,level,sps)%dwtemp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
      end do allocatedomains

      print *,'domains allocated'
!           if(ierr /= 0)                       &
!                call terminate("memory?") 
      
      ! Initialize the temporary arrays.
      dRdExtraErr = 0
      dRdExtraAdj = 0
      dRdExtraFD  = 0


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


!     ******************************************************************
!     *                                                                *
!     *  Compute d(dw)/d(w) using Tapenade                             *
!     *                                                                *
!     ******************************************************************

      ! Get the initial AD time.

      call mpi_barrier(SUmb_comm_world, ierr)
      if( myID==0 ) call cpu_time(time(1))
      
      !print *,'Entering Domain loop'
      domainLoopAD: do nn=1,nDom!97,97!1,1!nDom
         
         ! Loop over the number of time instances for this block.

         spectralLoop: do sps=1,nTimeIntervalsSpectral
            !print *,'Setting Pointers'
            call setPointersAdj(nn,level,sps)

            ! Loop over location of output (R) cell of residual
            do kCell = 2, kl
               do jCell = 2, jl
                  do iCell = 2, il
 !                    print *,'indices',icell,jcell,kcell
                     ! Copy the state w to the wAdj array in the stencil
!                     call copyADjointStencil(wAdj, xAdj, iCell, jCell, kCell)                  
                     call copyADjointStencil(wAdj, xAdj,xBlockCornerAdj,alphaAdj,&
           betaAdj,MachAdj,&
           machCoefAdj,machGridAdj,iCell, jCell, kCell,prefAdj,&
           rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
           rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
           murefAdj, timerefAdj,pInfCorrAdj,liftIndex)
                     !print *,'liftiindex',liftindex
                     !print *,'rotcenteradj',nn,icell,jcell,kcell,rotcenteradj
!copyADjointStencil(wAdj, xAdj,alphaAdj,betaAdj,MachAdj,&
!                          machCoefAdj,iCell, jCell, kCell,prefAdj,&
!                          rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
!                          rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
!                          murefAdj, timerefAdj,pInfCorrAdj,liftIndex)
!                     print *,'Stencil Copied'

                     mLoop: do m = 1, nw           ! Loop over output cell residuals (R)
!                        print *,'initializing variables'
                        ! Initialize the seed for the reverse mode
                        dwAdjb(:) = 0.; dwAdjb(m) = 1.
                        dwAdj(:)  = 0.
                        wAdjb(:,:,:,:)  = 0.  !dR(m)/dw
                        machadjb = 0
                        

  !                      print *,'calling reverse mode'
!                        print *,'secondhalo',secondhalo
                        
                        ! Call reverse mode of residual computation
                        call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
&  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
&  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
&  icell, jcell, kcell, nn, sps, correctfork, secondhalo, prefadj, &
&  rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj, &
&  rotrateadjb, rotcenteradj, murefadj, timerefadj, pinfcorradj, &
&  liftindex)


                        ! Store the block Jacobians (by rows).
                        idxres   = globalCell(iCell,jCell,kCell)*nw+m
                        !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
                        if( idxres>=0) then
!!$                           dRdExtraAdj(idxres,nDesignMach,nn,sps) = machadjb
!!$                           pvrlocal(m) = machadjb
!!$                           dRdExtraAdj(idxres,nDesignAOA,nn,sps) = alphaadjb
!!$                           dRdExtraAdj(idxres,nDesignSSA,nn,sps) = betaadjb
!!$                           !print *,'results',machadjb,rotrateadjb,idxres,nDesignRotX,nn,sps
!!$                           dRdExtraAdj(idxres,nDesignRotX,nn,sps) = rotrateadjb(1)!*timeref
!!$                           dRdExtraAdj(idxres,nDesignRotY,nn,sps) = rotrateadjb(2)
!!$                           dRdExtraAdj(idxres,nDesignRotZ,nn,sps) = rotrateadjb(3)
                           dRdExtraAdj(idxres,nDesignMach,1,sps) = machadjb
                           pvrlocal(m) = machadjb
                           write(unitM,10) machadjb,nn,icell,jcell,kcell,m,idxres
10                         format(1x,'Mach ',f18.10,6I8)
                           dRdExtraAdj(idxres,nDesignMachGrid,1,sps) = machgridadjb
                           !pvrlocal(m) = machGridadjb
                           dRdExtraAdj(idxres,nDesignAOA,1,sps) = alphaadjb
                           write(unitAoA,11) alphaadjb,nn,icell,jcell,kcell,m,idxres
11                         format(1x,'AoA ',f18.10,6I8)
                           dRdExtraAdj(idxres,nDesignSSA,1,sps) = betaadjb
                           write(unitSSA,12) betaadjb,nn,icell,jcell,kcell,m,idxres
12                         format(1x,'SSA ',f18.10,6I8)
                           !print *,'results',machadjb,rotrateadjb,idxres,nDesignRotX,nn,sps
                           dRdExtraAdj(idxres,nDesignRotX,1,sps) = rotrateadjb(1)!*timeref
                           !pvrlocal(m) = rotrateadjb(1)*timeref
                           write(unitRotx,13) rotrateadjb(1),nn,icell,jcell,kcell,m,idxres
13                         format(1x,'Rotx ',f18.10,6I8)
                           write(unitRoty,14) rotrateadjb(2),nn,icell,jcell,kcell,m,idxres
14                         format(1x,'Roty ',f18.10,6I8)
                           write(unitRotz,15) rotrateadjb(3),nn,icell,jcell,kcell,m,idxres
15                         format(1x,'Rotz ',f18.10,6I8)
                           dRdExtraAdj(idxres,nDesignRotY,1,sps) = rotrateadjb(2)
                           dRdExtraAdj(idxres,nDesignRotZ,1,sps) = rotrateadjb(3)
                           !pvrlocal(m) = rotrateadjb(3)*timeref
                        endif
                                       
                     end do mLoop
                     idxmgb = globalCell(icell,jcell,kcell)
                
                     test = sum(pvrlocal(:))
                     !print *,'test',test
                     if ( test.ne.0 .and. idxmgb.ne.-5 .and. idxmgb>=0 .and. idxmgb<nCellsGlobal) then
                        !print *,'setting PETSc Vector',sum(wAdjB(icell,jcell,kcell,:))
                        !pvrlocal(:) = wFD2(iCell-1, jCell-1, kCell-1,:)
                        
                        !                call VecSetValuesBlocked(dJdW, 1, idxmgb, dJdWlocal, &
                        !                                         INSERT_VALUES, PETScIerr)
                        call VecSetValuesBlocked(pvr, 1, idxmgb, pvrlocal, &
                             INSERT_VALUES, PETScIerr)
                        
                        if( PETScIerr/=0 ) then
                           write(errorMessage,99) &
                                "Error in VecSetValuesBlocked for global node", &
                           idxmgb
                           call terminate("setupADjointRHSAeroCoeff", &
                                errorMessage)
                        endif
                     endif
                  end do !iCell
               end do !jCell
            end do! kCell
         end do spectralLoop
      end do domainLoopAD
      !print *,'AD Completed'
      ! Get new time and compute the elapsed AD time.
!      stop
      call mpi_barrier(SUmb_comm_world, ierr)
      if(myID == 0) then
         call cpu_time(time(2))
         timeAdj = time(2)-time(1)
      endif

!!$!
!!$!     Compute d(dw)/d(w) using central finite differences
!!$!_______________________________________________________
!!$       
!!$       deltaw = 1.d-5
!!$       !print *, "deltaw=", deltaw
!!$       wFD2(:,:,:,:) = 0. 
!!$       call cpu_time(time(3))
!!$       
!!$       groundLevel = 1
!!$       sps = 1
!!$       nn=97
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       ! if (istate==icell .or. jstate==jcell .or. kstate==kcell) then
!!$       !Remember current values
!!$       dwtemp(:,:,:,:) = dw(:,:,:,:)
!!$       wtemp(:,:,:,:) = w(:,:,:,:)
!!$       ptemp(:,:,:) = p(:,:,:)
!!$       ExtraAdjRef = mach
!!$       
!!$       Mach = ExtraAdjRef + deltaw
!!$       
!!$       !call checkInputParam
!!$       
!!$       call referenceState
!!$       
!!$       call setFlowInfinityState
!!$       
!!$!       wAdjtemp = w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:)
!!$!       !                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$!       call computePressureAdj(wAdjtemp, pAdjtemp)
!!$!       p(istate, jstate, kstate) = pAdjtemp(0,0,0)
!!$       
!!$       !print *, "Calling Residual ="
!!$       call initres(1_intType, nwf)
!!$       
!!$       call applyAllBC(secondHalo)
!!$       
!!$       ! Exchange the solution. Either whalo1 or whalo2
!!$       ! must be called.
!!$       if( secondHalo ) then
!!$          call whalo2(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       else
!!$          call whalo1(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       endif
!!$       !if( secondHalo ) then
!!$       !   !  write(*,*)'2ndHalo..........'
!!$       !   call whalo2(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !else
!!$       !   ! write(*,*)'1stHalo..........'
!!$       !   call whalo1(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !endif
!!$       
!!$       call residual
!!$       !print *, "Called Residual =", dw(istate,jstate,kstate,n)
!!$       nn=97
!!$       groundlevel = 1
!!$       sps = 1
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       dwp(:,:,:,:) = dw(:,:,:,:)
!!$       dw(:,:,:,:) = dwtemp(:,:,:,:)
!!$       
!!$       Mach = ExtraAdjRef - deltaw
!!$       
!!$       !call checkInputParam
!!$       
!!$       call referenceState
!!$       
!!$       call setFlowInfinityState
!!$       
!!$!       wAdjtemp = w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:)
!!$!       !                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$!       call computePressureAdj(wAdjtemp, pAdjtemp)
!!$       !                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$!       p(istate, jstate, kstate) = pAdjtemp(0,0,0)
!!$       
!!$       !print *, "Calling Residual ="
!!$       call initres(1_intType, nwf)
!!$       
!!$       call applyAllBC(secondHalo)
!!$       
!!$       ! Exchange the solution. Either whalo1 or whalo2
!!$       ! must be called.
!!$       if( secondHalo ) then
!!$          call whalo2(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       else
!!$          call whalo1(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       endif
!!$       !if( secondHalo ) then
!!$       !   !write(*,*)'2ndHalo..........'
!!$       !   call whalo2(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !else
!!$       !   !write(*,*)'1stHalo..........'
!!$       !   call whalo1(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !endif
!!$       
!!$       call residual
!!$       !print *, "Called Residual 2 =", dw(istate,jstate,kstate,n)
!!$       nn=97
!!$       groundlevel = 1
!!$       sps = 1
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       dwm(:,:,:,:) = dw(:,:,:,:)
!!$       dw(:,:,:,:) = dwtemp(:,:,:,:)
!!$       
!!$       !                             end if
!!$       
!!$       
!!$       !                   ival = 0
!!$       do kCell = 2, kl
!!$          do jCell = 2, jl
!!$             do iCell = 2, il
!!$                do iRes = 1, nw           ! Loop over output cell residuals (R)
!!$                   ! Loop over location of output (R) cell of residual
!!$                   
!!$                   !column = (n-1) +(istate-2)*nw +(jstate-2)*nw*nx +(kstate-2)*nw*nx*ny
!!$                   wFD2(iCell-1, jCell-1, kCell-1, iRes) = (dwp(iCell,jCell,kCell,iRes)-dwm(iCell,jCell,kCell,iRes))/(2.0*deltaw)
!!$                   
!!$                   idxres   = globalCell(iCell,jCell,kCell)*nw+ires
!!$                   !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
!!$                   if( idxres>=0) then
!!$                      if (wFD2(iCell-1, jCell-1, kCell-1, iRes).ne. zero) then
!!$                         !dRdExtraFD(idxres,nDesignMach,nn,sps) = wFD2(iCell-1, jCell-1, kCell-1, iRes)
!!$                         dRdExtraFD(idxres,nDesignMach,1,sps) = wFD2(iCell-1, jCell-1, kCell-1, iRes)
!!$                         !Aad(m,:)  = wFD2(iCell-1, jCell-1, kCell-1, iRes)wAdjB( 0, 0, 0,:)
!!$                         !print *,'setting',dRdW, PETScOne, idxres, PETScOne, idxstate,   &
!!$                         !     dRdwFD(idxres,idxstate,nn,sps)
!                         call MatSetValues(dRda, 1, idxres-1, 1, nDesignMach-1,   &
!                              dRdExtraFD(idxres,nDesignMach,nn,sps), INSERT_VALUES, PETScIerr)
!                         if( PETScIerr/=0 ) &
!                              print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
!!$                      end if
!!$                   endif
!!$                end do
!!$!!                idxmgb = globalCell(icell,jcell,kcell)
!!$                
!                test = sum(wFD2(iCell-1, jCell-1, kCell-1,:))
!                !print *,'test',test
!                if ( test.ne.0 .and. idxmgb.ne.-5 .and. idxmgb>=0 .and. idxmgb<nCellsGlobal) then
!                   !print *,'setting PETSc Vector',sum(wAdjB(icell,jcell,kcell,:))
!                   pvrlocal(:) = wFD2(iCell-1, jCell-1, kCell-1,:)
!                   
!                   !                call VecSetValuesBlocked(dJdW, 1, idxmgb, dJdWlocal, &
!                   !                                         INSERT_VALUES, PETScIerr)
!                   call VecSetValuesBlocked(pvr, 1, idxmgb, pvrlocal, &
!                        INSERT_VALUES, PETScIerr)
!                   
!                   if( PETScIerr/=0 ) then
!                      write(errorMessage,99) &
!                           "Error in VecSetValuesBlocked for global node", &
!                           idxmgb
!                      call terminate("setupADjointRHSAeroCoeff", &
!                           errorMessage)
!                   endif
!                endif
!!$                !idxmgb = globalCell(iCell,jCell,kCell)
!!$                !print *,'globalcell',idxmgb,globalCell(iCell,jCell,kCell)
!!$                ! >>> center block A < W(i,j,k)
!!$                
!!$                !idxngb = idxmgb
!!$                !print *,'indicies0',idxmgb,idxngb
!!$                !call MatSetValuesBlocked(dRdW, 1, idxmgb, 1, idxngb, &
!!$                !     Aad, INSERT_VALUES,PETScIerr)
!!$                !if( PETScIerr/=0 ) &
!!$                !      print *,'matrix setting error'!call errAssemb("MatSetValuesBlocked", "Aad")
!!$                
!!$                ! if ( wAdjb(ii,jj,kk,l).ne. 0.0) then
!!$                !    print *,'wAdjb', wAdjb(ii,jj,kk,l),i,j,k,l
!!$                ! endif
!!$                !end if
!!$                
!!$             end do
!!$          end do
!!$       end do
!!$!       w(istate,jstate,kstate,n) = wAdjRef
!!$       nn=97
!!$       groundlevel = 1
!!$       sps = 1
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       dw(:,:,:,:) = dwtemp(:,:,:,:)
!!$       w(:,:,:,:) = wtemp(:,:,:,:)
!!$       p(:,:,:) = ptemp(:,:,:)
!!$       mach = ExtraAdjRef
!!$
!       call MatAssemblyBegin(dRda,MAT_FINAL_ASSEMBLY,PETScIerr)!!
!
!       if( PETScIerr/=0 ) &
!            call terminate("setupADjointMatrix","Error in MatAssemblyBegin")
!       !!
!
!       call MatAssemblyEnd  (dRda,MAT_FINAL_ASSEMBLY,PETScIerr)
!       
!       if( PETScIerr/=0 ) &
!            call terminate("setupADjointMatrix","Error in MatAssemblyEnd")!
!
!       !if( debug ) then
!          call MatView(dRda,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
!          !call MatView(dRdW,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
!          if( PETScIerr/=0 ) &
!               call terminate("setupADjointMatrix", "Error in MatView")
!          !pause
!       !endif
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

!!$!     ******************************************************************
!!$!     *                                                                *
!!$!     *  Compute d(dw)/d(w) using central finite differences           *
!!$!     *                                                                *
!!$!     ******************************************************************
!!$
!      ! Get the initial FD time.
!      print *,'starting FD calculation'
!      call mpi_barrier(SUmb_comm_world, ierr)
!      if(myID == 0) call cpu_time(time(3))
!!
!
!      deltaw = 1.d-6
!      wFD(:,:,:,:) = 0. 
!      call cpu_time(time(3))!
!
!      print *,'entering FD loop'
!      domainResidualLoopFDorig: do nn=1,nDom         
!
!         ! Determine the number of time instances for this block
!        
!         nTime     = nTimeIntervalsSpectral!
!
!         ! Loop over the number of time instances for this block.!!
!
!         spectralLoop2: do sps=1,nTime
!!            print *,'setting Pointers'
!            call setPointersAdj(nn,level,sps)!
!
!
!               ! Loop over location of output (R) cell of residual
!               do kCell = 2, kl
!                  !print *, "iRes, kCell =", iRes, kCell
!                  do jCell = 2, jl
!                     do iCell = 2, il
!                        
!                        do m = 1, nw           ! Loop over output cell residuals (R)
!!                           print *,'indices',icell,jcell,kcell,m
!                        ! Copy the state in the stencil
!                        ! actually no need to call again, but ...
!!                           call copyADjointStencil(wAdj,xAdj, iCell, jCell, kCell)
!                           call copyADjointStencil(wAdj, xAdj,alphaAdj,betaAdj,&
!                          MachAdj,MachCoefAdj,iCell, jCell, kCell,prefAdj,&
!                          rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
!                          rhoinfAdj, pinfAdj,&
!                          murefAdj, timerefAdj,pInfCorrAdj)
!!                           print *,'stencil copied'
!                           ! Loop over the inputs (w) 
!                           !print *,'secondhalo',secondhalo
!                        
!                           do ii = -2, 2
!                              do jj = -2, 2
!                                 do kk = -2, 2
!                                    do n = 1, nw
!                                       i = iCell+ii
!                                       j = jCell+jj
!                                       k = kCell+kk
! !                                      print *,'secondary FD indices',i,j,k,n,ii,jj,kk
!                                       if (ii==0 .or. jj==0 .or. kk==0) then
!                                          wAdjRef = wAdj(ii,jj,kk,n)
!                                          wAdj(ii,jj,kk,n) = wAdjRef + deltaw
!                                          call computeRAdjoint(wAdj,xAdj,dwAdjP,alphaAdj,&
!                                               betaAdj,MachAdj, MachCoefAdj,&
!                                               iCell, jCell,  kCell, &
!                                               nn,sps, correctForK,secondHalo,prefAdj,&
 !                                              rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
 !                                              rhoinfAdj, pinfAdj,&
 !                                              murefAdj, timerefAdj,pInfCorrAdj)
 !                                         !call computeRAdjoint(wAdj,xAdj,dwAdjP,   &
 !                                         !     iCell, jCell,  kCell, &
!                                          !     nn,sps, correctForK,secondHalo)
!                                          
!                                          wAdj(ii,jj,kk,n) = wAdjRef - deltaw
!                                          call computeRAdjoint(wAdj,xAdj,dwAdjM,alphaAdj,&
!                                               betaAdj,MachAdj, MachCoefAdj,&
!                                               iCell, jCell,  kCell, &
!                                               nn,sps, correctForK,secondHalo,prefAdj,&
!                                               rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
!                                               rhoinfAdj, pinfAdj,&
!                                               murefAdj, timerefAdj,pInfCorrAdj)
!                                          !call computeRAdjoint(wAdj,xAdj,dwAdjM,   &
!                                          !     iCell, jCell,  kCell, &
!                                          !     nn,sps, correctForK,secondHalo)!!!
!
!                                          
!                                          wFD(ii,jj,kk,n) = (dwAdjP(m)-dwAdjM(m))/(2.0*deltaw)
!                                          wAdj(ii,jj,kk,n) = wAdjRef!
!
!                                          if(i>zero .and. j>zero .and. k>zero .and. i<il .and. j<jl .and. k<kl)then
!!                                             print *,'global FD indicies',idxstate,idxres
!                                             idxstate = globalCell(i,j,k)*nw+n
!                                             idxres   = globalCell(iCell,jCell,kCell)*nw+m
!                                             if( idxres>=0 .and. idxstate>=0) then
!                                                dRdwFD(idxres,idxstate,nn,sps) = wFD(ii,jj,kk,n)!
!                                             endif
!                                          end if
!                                       end if
!                                    end do ! n
!                                 end do ! ii
!                              end do !jj
!                           end do !kk
!                        
!                        !print 111, iRes, iCell-1, jCell-1, kCell-1, dRdw(iCell-1, jCell-1, kCell-1, iRes), dRdwFD(iCell-1, jCell-1, kCell-1, iRes)
!                     end do !iRes
!                  end do !iCell
!               end do !jCell
!            end do !kCell
!         end do spectralLoop2
!      end do domainResidualLoopFDorig
!      
!      call cpu_time(time(4))
!      timeFD = time(4)-time(3)
!!$
!!$!*******************************************
!!$!Now repeat for Alpha
!!$!*******************************************
!!$
!!$!
!!$!     Compute d(dw)/d(w) using central finite differences
!!$!_______________________________________________________
!!$       
!!$       deltaw = 1.d-5
!!$       !print *, "deltaw=", deltaw
!!$       wFD2(:,:,:,:) = 0. 
!!$       call cpu_time(time(3))
!!$       
!!$       groundLevel = 1
!!$       sps = 1
!!$       nn=97
!!$
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       ! if (istate==icell .or. jstate==jcell .or. kstate==kcell) then
!!$       !Remember current values
!!$       dwtemp(:,:,:,:) = dw(:,:,:,:)
!!$       wtemp(:,:,:,:) = w(:,:,:,:)
!!$       ptemp(:,:,:) = p(:,:,:)
!!$       ExtraAdjRef = alphaadj
!!$       
!!$       alphaadj = ExtraAdjRef + deltaw
!!$       
!!$       call adjustInflowAngleAdj(alphaAdj,betaAdj,velDirFreestreamAdj,&
!!$     liftDirectionAdj,dragDirectionAdj,liftIndex)
!!$       
!!$       liftDirection = liftDirectionAdj
!!$       !print *,'liftdirection 1',liftindex,liftDirection
!!$       veldirfreestream = veldirfreestreamAdj
!!$       dragdirection = dragdirectionAdj
!!$       !call checkInputParam
!!$       
!!$       call referenceState
!!$       
!!$       call setFlowInfinityState
!!$
!!$       do mm=1,nTimeIntervalsSpectral
!!$          
!!$          ! Compute the time, which corresponds to this spectral solution.
!!$          ! For steady and unsteady mode this is simply the restart time;
!!$          ! for the spectral mode the periodic time must be taken into
!!$          ! account, which can be different for every section.
!!$          
!!$          t = timeUnsteadyRestart
!!$          
!!$          if(equationMode == timeSpectral) then
!!$             do nn=1,nSections
!!$                t(nn) = t(nn) + (mm-1)*sections(nn)%timePeriod &
!!$                     /         real(nTimeIntervalsSpectral,realType)
!!$             enddo
!!$          endif
!!$          
!!$          call gridVelocitiesFineLevel(.false., t, mm)
!!$          call gridVelocitiesCoarseLevels(mm)
!!$          call normalVelocitiesAllLevels(mm)
!!$          
!!$          call slipVelocitiesFineLevel(.false., t, mm)
!!$          call slipVelocitiesCoarseLevels(mm)
!!$          
!!$       enddo
!!$       
!!$!       wAdjtemp = w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:)
!!$!       !                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$!       call computePressureAdj(wAdjtemp, pAdjtemp)
!!$!       p(istate, jstate, kstate) = pAdjtemp(0,0,0)
!!$       
!!$
!!$       call applyAllBC(secondHalo)
!!$       
!!$       ! Exchange the solution. Either whalo1 or whalo2
!!$       ! must be called.
!!$       if( secondHalo ) then
!!$          call whalo2(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       else
!!$          call whalo1(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       endif
!!$       !if( secondHalo ) then
!!$       !   !  write(*,*)'2ndHalo..........'
!!$       !   call whalo2(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !else
!!$       !   ! write(*,*)'1stHalo..........'
!!$       !   call whalo1(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !endif
!!$
!!$       !print *, "Calling Residual ="
!!$       call initres(1_intType, nwf)
!!$       
!!$       call residual
!!$       !print *, "Called Residual =", dw(istate,jstate,kstate,n)
!!$       groundLevel = 1
!!$       sps = 1
!!$       nn=97
!!$
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       dwp(:,:,:,:) = dw(:,:,:,:)
!!$       dw(:,:,:,:) = dwtemp(:,:,:,:)
!!$       
!!$       alphaadj = ExtraAdjRef - deltaw
!!$       
!!$       call adjustInflowAngleAdj(alphaAdj,betaAdj,velDirFreestreamAdj,&
!!$            liftDirectionAdj,dragDirectionAdj,liftIndex)
!!$       
!!$       liftDirection = liftDirectionAdj
!!$       !print *,'liftdirection 2',liftindex,liftDirection
!!$       veldirfreestream = veldirfreestreamAdj
!!$       dragdirection = dragdirectionAdj
!!$       !call checkInputParam
!!$       
!!$       call referenceState
!!$       
!!$       call setFlowInfinityState
!!$
!!$       do mm=1,nTimeIntervalsSpectral
!!$          
!!$          ! Compute the time, which corresponds to this spectral solution.
!!$          ! For steady and unsteady mode this is simply the restart time;
!!$          ! for the spectral mode the periodic time must be taken into
!!$          ! account, which can be different for every section.
!!$          
!!$          t = timeUnsteadyRestart
!!$          
!!$          if(equationMode == timeSpectral) then
!!$             do nn=1,nSections
!!$                t(nn) = t(nn) + (mm-1)*sections(nn)%timePeriod &
!!$                     /         real(nTimeIntervalsSpectral,realType)
!!$             enddo
!!$          endif
!!$          
!!$          call gridVelocitiesFineLevel(.false., t, mm)
!!$          call gridVelocitiesCoarseLevels(mm)
!!$          call normalVelocitiesAllLevels(mm)
!!$          
!!$          call slipVelocitiesFineLevel(.false., t, mm)
!!$          call slipVelocitiesCoarseLevels(mm)
!!$          
!!$       enddo
!!$       
!!$!       wAdjtemp = w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:)
!!$!       !                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$!       call computePressureAdj(wAdjtemp, pAdjtemp)
!!$       !                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$!       p(istate, jstate, kstate) = pAdjtemp(0,0,0)
!!$       
!!$      
!!$       call applyAllBC(secondHalo)
!!$       
!!$       ! Exchange the solution. Either whalo1 or whalo2
!!$       ! must be called.
!!$       if( secondHalo ) then
!!$          call whalo2(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       else
!!$          call whalo1(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       endif
!!$       !if( secondHalo ) then
!!$       !   !write(*,*)'2ndHalo..........'
!!$       !   call whalo2(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !else
!!$       !   !write(*,*)'1stHalo..........'
!!$       !   call whalo1(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !endif
!!$       
!!$       !print *, "Calling Residual ="
!!$       call initres(1_intType, nwf)
!!$       call residual
!!$       !print *, "Called Residual 2 =", dw(istate,jstate,kstate,n)
!!$       groundLevel = 1
!!$       sps = 1
!!$       nn=97
!!$
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       dwm(:,:,:,:) = dw(:,:,:,:)
!!$       dw(:,:,:,:) = dwtemp(:,:,:,:)
!!$       
!!$       !                             end if
!!$       
!!$       
!!$       !                   ival = 0
!!$       do kCell = 2, kl
!!$          do jCell = 2, jl
!!$             do iCell = 2, il
!!$                do iRes = 1, nw           
!!$                   ! Loop over location of output (R) cell of residual
!!$                                 
!!$                   wFD2(iCell-1, jCell-1, kCell-1, iRes) = (dwp(iCell,jCell,kCell,iRes)-dwm(iCell,jCell,kCell,iRes))/(2.0*deltaw)
!!$                   
!!$                   idxres   = globalCell(iCell,jCell,kCell)*nw+ires
!!$                   !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
!!$                   if( idxres>=0) then
!!$                      !print *,'fdresultaoa',wfd2(iCell-1, jCell-1, kCell-1, iRes),iCell-1, jCell-1, kCell-1, iRes
!!$                      if (wFD2(iCell-1, jCell-1, kCell-1, iRes).ne. zero) then
!!$                         !dRdExtraFD(idxres,nDesignAOA,nn,sps) = wFD2(iCell-1, jCell-1, kCell-1, iRes)
!!$                         dRdExtraFD(idxres,nDesignAOA,1,sps) = wFD2(iCell-1, jCell-1, kCell-1, iRes)
!                         call MatSetValues(dRda, 1, idxres-1, 1, nDesignMach-1,   &
!                              dRdExtraFD(idxres,nDesignMach,nn,sps), INSERT_VALUES, PETScIerr)
!                         if( PETScIerr/=0 ) &
!                              print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
!!$                      end if
!!$                   endif
!!$                end do
!!$
!!$                
!!$             end do
!!$          end do
!!$       end do
!!$!       w(istate,jstate,kstate,n) = wAdjRef
!!$       groundLevel = 1
!!$       sps = 1
!!$       nn=97
!!$
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       dw(:,:,:,:) = dwtemp(:,:,:,:)
!!$       w(:,:,:,:) = wtemp(:,:,:,:)
!!$       p(:,:,:) = ptemp(:,:,:)
!!$       alphaAdj = ExtraAdjRef
!!$
!!$       call adjustInflowAngleAdj(alphaAdj,betaAdj,velDirFreestreamAdj,&
!!$            liftDirectionAdj,dragDirectionAdj,liftIndex)
!!$       
!!$       liftDirection = liftDirectionAdj
!!$       veldirfreestream = veldirfreestreamAdj
!!$       dragdirection = dragdirectionAdj
!!$
!!$!******************************************
!!$!Now Repeat for Beta
!!$!******************************************
!!$!
!!$!     Compute d(dw)/d(w) using central finite differences
!!$!_______________________________________________________
!!$       
!!$       deltaw = 1.d-5
!!$       !print *, "deltaw=", deltaw
!!$       wFD2(:,:,:,:) = 0. 
!!$       call cpu_time(time(3))
!!$       
!!$       groundLevel = 1
!!$       sps = 1
!!$       nn=97
!!$       
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       ! if (istate==icell .or. jstate==jcell .or. kstate==kcell) then
!!$       !Remember current values
!!$       dwtemp(:,:,:,:) = dw(:,:,:,:)
!!$       wtemp(:,:,:,:) = w(:,:,:,:)
!!$       ptemp(:,:,:) = p(:,:,:)
!!$       ExtraAdjRef = betaadj
!!$       
!!$       betaadj = ExtraAdjRef + deltaw
!!$       
!!$       !call checkInputParam
!!$       call adjustInflowAngleAdj(alphaAdj,betaAdj,velDirFreestreamAdj,&
!!$            liftDirectionAdj,dragDirectionAdj,liftIndex)
!!$       
!!$       liftDirection = liftDirectionAdj
!!$       veldirfreestream = veldirfreestreamAdj
!!$       dragdirection = dragdirectionAdj
!!$       call referenceState
!!$       
!!$       call setFlowInfinityState
!!$              do mm=1,nTimeIntervalsSpectral
!!$          
!!$          ! Compute the time, which corresponds to this spectral solution.
!!$          ! For steady and unsteady mode this is simply the restart time;
!!$          ! for the spectral mode the periodic time must be taken into
!!$          ! account, which can be different for every section.
!!$          
!!$          t = timeUnsteadyRestart
!!$          
!!$          if(equationMode == timeSpectral) then
!!$             do nn=1,nSections
!!$                t(nn) = t(nn) + (mm-1)*sections(nn)%timePeriod &
!!$                     /         real(nTimeIntervalsSpectral,realType)
!!$             enddo
!!$          endif
!!$          
!!$          call gridVelocitiesFineLevel(.false., t, mm)
!!$          call gridVelocitiesCoarseLevels(mm)
!!$          call normalVelocitiesAllLevels(mm)
!!$          
!!$          call slipVelocitiesFineLevel(.false., t, mm)
!!$          call slipVelocitiesCoarseLevels(mm)
!!$          
!!$       enddo
!!$!       wAdjtemp = w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:)
!!$!       !                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$!       call computePressureAdj(wAdjtemp, pAdjtemp)
!!$!       p(istate, jstate, kstate) = pAdjtemp(0,0,0)
!!$       
!!$       !print *, "Calling Residual ="
!!$       call initres(1_intType, nwf)
!!$       
!!$       call applyAllBC(secondHalo)
!!$       
!!$       ! Exchange the solution. Either whalo1 or whalo2
!!$       ! must be called.
!!$       if( secondHalo ) then
!!$          call whalo2(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       else
!!$          call whalo1(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       endif
!!$       !if( secondHalo ) then
!!$       !   !  write(*,*)'2ndHalo..........'
!!$       !   call whalo2(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !else
!!$       !   ! write(*,*)'1stHalo..........'
!!$       !   call whalo1(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !endif
!!$       
!!$       call residual
!!$       !print *, "Called Residual =", dw(istate,jstate,kstate,n)
!!$       groundLevel = 1
!!$       sps = 1
!!$       nn=97
!!$
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       dwp(:,:,:,:) = dw(:,:,:,:)
!!$       dw(:,:,:,:) = dwtemp(:,:,:,:)
!!$       
!!$       betaadj = ExtraAdjRef - deltaw
!!$       
!!$       !call checkInputParam
!!$       call adjustInflowAngleAdj(alphaAdj,betaAdj,velDirFreestreamAdj,&
!!$            liftDirectionAdj,dragDirectionAdj,liftIndex)
!!$       
!!$       liftDirection = liftDirectionAdj
!!$       veldirfreestream = veldirfreestreamAdj
!!$       dragdirection = dragdirectionAdj
!!$       
!!$       call referenceState
!!$       
!!$       call setFlowInfinityState
!!$             do mm=1,nTimeIntervalsSpectral
!!$          
!!$          ! Compute the time, which corresponds to this spectral solution.
!!$          ! For steady and unsteady mode this is simply the restart time;
!!$          ! for the spectral mode the periodic time must be taken into
!!$          ! account, which can be different for every section.
!!$          
!!$          t = timeUnsteadyRestart
!!$          
!!$          if(equationMode == timeSpectral) then
!!$             do nn=1,nSections
!!$                t(nn) = t(nn) + (mm-1)*sections(nn)%timePeriod &
!!$                     /         real(nTimeIntervalsSpectral,realType)
!!$             enddo
!!$          endif
!!$          
!!$          call gridVelocitiesFineLevel(.false., t, mm)
!!$          call gridVelocitiesCoarseLevels(mm)
!!$          call normalVelocitiesAllLevels(mm)
!!$          
!!$          call slipVelocitiesFineLevel(.false., t, mm)
!!$          call slipVelocitiesCoarseLevels(mm)
!!$          
!!$       enddo 
!!$!       wAdjtemp = w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:)
!!$!       !                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$!       call computePressureAdj(wAdjtemp, pAdjtemp)
!!$       !                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$!       p(istate, jstate, kstate) = pAdjtemp(0,0,0)
!!$       
!!$       !print *, "Calling Residual ="
!!$       call initres(1_intType, nwf)
!!$       
!!$       call applyAllBC(secondHalo)
!!$       
!!$       ! Exchange the solution. Either whalo1 or whalo2
!!$       ! must be called.
!!$       if( secondHalo ) then
!!$          call whalo2(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       else
!!$          call whalo1(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       endif
!!$       !if( secondHalo ) then
!!$       !   !write(*,*)'2ndHalo..........'
!!$       !   call whalo2(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !else
!!$       !   !write(*,*)'1stHalo..........'
!!$       !   call whalo1(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !endif
!!$       
!!$       call residual
!!$       !print *, "Called Residual 2 =", dw(istate,jstate,kstate,n)
!!$       groundLevel = 1
!!$       sps = 1
!!$       nn=97
!!$
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       dwm(:,:,:,:) = dw(:,:,:,:)
!!$       dw(:,:,:,:) = dwtemp(:,:,:,:)
!!$       
!!$       !                             end if
!!$       
!!$       
!!$       !                   ival = 0
!!$       do kCell = 2, kl
!!$          do jCell = 2, jl
!!$             do iCell = 2, il
!!$                do iRes = 1, nw   
!!$                   ! Loop over location of output (R) cell of residual
!!$       
!!$                   wFD2(iCell-1, jCell-1, kCell-1, iRes) = (dwp(iCell,jCell,kCell,iRes)-dwm(iCell,jCell,kCell,iRes))/(2.0*deltaw)
!!$                   
!!$                   idxres   = globalCell(iCell,jCell,kCell)*nw+ires
!!$                   !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
!!$                   if( idxres>=0) then
!!$                      if (wFD2(iCell-1, jCell-1, kCell-1, iRes).ne. zero) then
!!$                         !dRdExtraFD(idxres,nDesignSSA,nn,sps) = wFD2(iCell-1, jCell-1, kCell-1, iRes)
!!$                         dRdExtraFD(idxres,nDesignSSA,1,sps) = wFD2(iCell-1, jCell-1, kCell-1, iRes)
!!$
!                         call MatSetValues(dRda, 1, idxres-1, 1, nDesignMach-1,   &
!                              dRdExtraFD(idxres,nDesignMach,nn,sps), INSERT_VALUES, PETScIerr)
!                         if( PETScIerr/=0 ) &
!                              print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
!!$                      end if
!!$                   endif
!!$                end do
!!$
!!$                
!!$             end do
!!$          end do
!!$       end do
!!$!       w(istate,jstate,kstate,n) = wAdjRef
!!$       groundLevel = 1
!!$       sps = 1
!!$       nn=97
!!$
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       dw(:,:,:,:) = dwtemp(:,:,:,:)
!!$       w(:,:,:,:) = wtemp(:,:,:,:)
!!$       p(:,:,:) = ptemp(:,:,:)
!!$       betaadj = ExtraAdjRef
!!$       call adjustInflowAngleAdj(alphaAdj,betaAdj,velDirFreestreamAdj,&
!!$            liftDirectionAdj,dragDirectionAdj,liftIndex)
!!$       
!!$       liftDirection = liftDirectionAdj
!!$       veldirfreestream = veldirfreestreamAdj
!!$       dragdirection = dragdirectionAdj
!!$!*************************
!!$!Now repeat for RotRate
!!$!*************************
!!$
!!$!
!!$!     Compute d(dw)/d(w) using central finite differences
!!$!_______________________________________________________
!!$       
!!$       deltaw = 1.d-5*timeref
!!$
!!$       !print *, "deltaw=", deltaw
!!$       wFD2(:,:,:,:) = 0. 
!!$       call cpu_time(time(3))
!!$       
!!$       groundLevel = 1
!!$       sps = 1
!!$       nn=97
!!$      
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       ! if (istate==icell .or. jstate==jcell .or. kstate==kcell) then
!!$       !Remember current values
!!$       dwtemp(:,:,:,:) = dw(:,:,:,:)
!!$       wtemp(:,:,:,:) = w(:,:,:,:)
!!$       ptemp(:,:,:) = p(:,:,:)
!!$       ExtraAdjRef = RotRateAdj(1)
!!$       
!!$       RotRateAdj(1) = ExtraAdjRef + deltaw
!!$
!!$       do nn=1,nDom
!!$          
!!$          ! Set the pointers for this block.
!!$          
!!$          call setPointers(nn, groundLevel, sps)
!!$       
!!$          cgnsDoms(nbkglobal)%rotRate(1) = RotRateAdj(1)/timeRef
!!$       enddo 
!!$       nn=97
!!$       !call checkInputParam
!!$       
!!$       call referenceState
!!$       
!!$       call setFlowInfinityState
!!$
!!$       do mm=1,nTimeIntervalsSpectral
!!$
!!$         ! Compute the time, which corresponds to this spectral solution.
!!$         ! For steady and unsteady mode this is simply the restart time;
!!$         ! for the spectral mode the periodic time must be taken into
!!$         ! account, which can be different for every section.
!!$
!!$         t = timeUnsteadyRestart
!!$
!!$         if(equationMode == timeSpectral) then
!!$           do nn=1,nSections
!!$             t(nn) = t(nn) + (mm-1)*sections(nn)%timePeriod &
!!$                   /         real(nTimeIntervalsSpectral,realType)
!!$           enddo
!!$         endif
!!$
!!$         call gridVelocitiesFineLevel(.false., t, mm)
!!$         call gridVelocitiesCoarseLevels(mm)
!!$         call normalVelocitiesAllLevels(mm)
!!$
!!$         call slipVelocitiesFineLevel(.false., t, mm)
!!$         call slipVelocitiesCoarseLevels(mm)
!!$
!!$       enddo
!!$       
!!$       
!!$!       wAdjtemp = w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:)
!!$!       !                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$!       call computePressureAdj(wAdjtemp, pAdjtemp)
!!$!       p(istate, jstate, kstate) = pAdjtemp(0,0,0)
!!$       
!!$       !print *, "Calling Residual ="
!!$       call initres(1_intType, nwf)
!!$       
!!$       call applyAllBC(secondHalo)
!!$       
!!$       ! Exchange the solution. Either whalo1 or whalo2
!!$       ! must be called.
!!$       if( secondHalo ) then
!!$          call whalo2(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       else
!!$          call whalo1(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       endif
!!$       !if( secondHalo ) then
!!$       !   !  write(*,*)'2ndHalo..........'
!!$       !   call whalo2(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !else
!!$       !   ! write(*,*)'1stHalo..........'
!!$       !   call whalo1(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !endif
!!$       
!!$       call residual
!!$       !print *, "Called Residual =", dw(istate,jstate,kstate,n)
!!$       groundLevel = 1
!!$       sps = 1
!!$       nn=97
!!$
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       dwp(:,:,:,:) = dw(:,:,:,:)
!!$       dw(:,:,:,:) = dwtemp(:,:,:,:)
!!$       
!!$       RotRateAdj(1) = ExtraAdjRef - deltaw
!!$       
!!$       do nn=1,nDom
!!$          
!!$          ! Set the pointers for this block.
!!$          
!!$          call setPointers(nn, groundLevel, sps)
!!$       
!!$          cgnsDoms(nbkglobal)%rotRate(1) = RotRateAdj(1)/timeRef
!!$       enddo 
!!$       nn=97
!!$       !call checkInputParam
!!$       
!!$       call referenceState
!!$       
!!$       call setFlowInfinityState
!!$       
!!$       do mm=1,nTimeIntervalsSpectral
!!$
!!$         ! Compute the time, which corresponds to this spectral solution.
!!$         ! For steady and unsteady mode this is simply the restart time;
!!$         ! for the spectral mode the periodic time must be taken into
!!$         ! account, which can be different for every section.
!!$
!!$         t = timeUnsteadyRestart
!!$
!!$         if(equationMode == timeSpectral) then
!!$           do nn=1,nSections
!!$             t(nn) = t(nn) + (mm-1)*sections(nn)%timePeriod &
!!$                   /         real(nTimeIntervalsSpectral,realType)
!!$           enddo
!!$         endif
!!$
!!$         call gridVelocitiesFineLevel(.false., t, mm)
!!$         call gridVelocitiesCoarseLevels(mm)
!!$         call normalVelocitiesAllLevels(mm)
!!$
!!$         call slipVelocitiesFineLevel(.false., t, mm)
!!$         call slipVelocitiesCoarseLevels(mm)
!!$
!!$      enddo
!!$
!!$!       wAdjtemp = w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:)
!!$!       !                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$!       call computePressureAdj(wAdjtemp, pAdjtemp)
!!$       !                   call computePressureAdj(w(istate-2:istate+2, jstate-2:jstate+2, kstate-2:kstate+2,:), pAdjtemp)
!!$!       p(istate, jstate, kstate) = pAdjtemp(0,0,0)
!!$      
!!$       !print *, "Calling Residual ="
!!$       call initres(1_intType, nwf)
!!$       
!!$       call applyAllBC(secondHalo)
!!$       
!!$       ! Exchange the solution. Either whalo1 or whalo2
!!$       ! must be called.
!!$       if( secondHalo ) then
!!$          call whalo2(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       else
!!$          call whalo1(currentLevel, 1_intType, nMGVar, .true., &
!!$               .true., .true.)
!!$       endif
!!$       !if( secondHalo ) then
!!$       !   !write(*,*)'2ndHalo..........'
!!$       !   call whalo2(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !else
!!$       !   !write(*,*)'1stHalo..........'
!!$       !   call whalo1(currentLevel, 1_intType, nVarInt, .true., &
!!$       !        .true., .true.)
!!$       !endif
!!$       
!!$       call residual
!!$       !print *, "Called Residual 2 =", dw(istate,jstate,kstate,n)
!!$       groundLevel = 1
!!$       sps = 1
!!$       nn=97
!!$
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       dwm(:,:,:,:) = dw(:,:,:,:)
!!$       dw(:,:,:,:) = dwtemp(:,:,:,:)
!!$       
!!$
!!$       do kCell = 2, kl
!!$          do jCell = 2, jl
!!$             do iCell = 2, il
!!$                do iRes = 1, nw         
!!$                   ! Loop over location of output (R) cell of residual
!!$                   
!!$                   wFD2(iCell-1, jCell-1, kCell-1, iRes) = (dwp(iCell,jCell,kCell,iRes)-dwm(iCell,jCell,kCell,iRes))/(2.0*deltaw)
!!$                   
!!$                   idxres   = globalCell(iCell,jCell,kCell)*nw+ires
!!$                   !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
!!$                   if( idxres>=0) then
!!$                      if (wFD2(iCell-1, jCell-1, kCell-1, iRes).ne. zero) then
!!$                         !dRdExtraFD(idxres,nDesignRotX,nn,sps) = wFD2(iCell-1, jCell-1, kCell-1, iRes)
!!$                         dRdExtraFD(idxres,nDesignRotX,1,sps) = wFD2(iCell-1, jCell-1, kCell-1, iRes)
!!$
!                         call MatSetValues(dRda, 1, idxres-1, 1, nDesignMach-1,   &
!                              dRdExtraFD(idxres,nDesignMach,nn,sps), INSERT_VALUES, PETScIerr)
!                         if( PETScIerr/=0 ) &
!                              print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
!!$                      end if
!!$                   endif
!!$                end do
!!$
!!$                
!!$             end do
!!$          end do
!!$       end do
!!$!       w(istate,jstate,kstate,n) = wAdjRef
!!$       groundLevel = 1
!!$       sps = 1
!!$       nn=97
!!$
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$       dw(:,:,:,:) = dwtemp(:,:,:,:)
!!$       w(:,:,:,:) = wtemp(:,:,:,:)
!!$       p(:,:,:) = ptemp(:,:,:)
!!$       rotRateAdj(1) = ExtraAdjRef
!!$
!!$       do nn=1,nDom
!!$          
!!$          ! Set the pointers for this block.
!!$          
!!$          call setPointers(nn, groundLevel, sps)
!!$       
!!$          cgnsDoms(nbkglobal)%rotRate(1) = RotRateAdj(1)/timeRef
!!$       enddo
!!$       nn=97
!!$!_______________________________________________________
!!$!
!!$!     Compute the errors in dR/dw
!!$!_______________________________________________________
!!$
!!$      !print *, "Computing the the error..."
!!$       groundLevel = 1
!!$       sps = 1
!!$       nn=97
!!$       
!!$       call setPointersAdj(nn,groundlevel,sps)
!!$    do m=1,nw
!!$      do iCell=2,il
!!$        do jCell=2,jl
!!$          do kCell=2,kl
!!$             idxres   = globalCell(iCell,jCell,kCell)*nw+m
!!$             if( idxres>=0) then
!!$                dRdExtraErr(idxres, ndesignMach, 1, 1) = dRdExtraAdj(idxres, nDesignMach, 1, 1) - dRdExtraFD(idxres, nDesignMach, 1, 1)
!!$                dRdExtraErr(idxres, ndesignAOA, 1, 1) = dRdExtraAdj(idxres, nDesignAOA, 1, 1) - dRdExtraFD(idxres, nDesignAOA, 1, 1)
!!$                dRdExtraErr(idxres, ndesignSSA, 1, 1) = dRdExtraAdj(idxres, nDesignSSA, 1, 1) - dRdExtraFD(idxres, nDesignSSA, 1, 1)
!!$                dRdExtraErr(idxres, ndesignRotX, 1, 1) = dRdExtraAdj(idxres, nDesignRotX, 1, 1) - dRdExtraFD(idxres, nDesignRotX, 1, 1)
!!$                !dRdExtraErr(idxres, ndesignRotY, 1, 1) = dRdExtraAdj(idxres, nDesignMach, 1, 1) - dRdExtraFD(idxres, nDesignMach, 1, 1)
!!$                !                           if(dRdwFD(idxres,idxstate,1,1).ne.0) then
!!$                if(abs(dRdExtraFD(idxres,nDesignMach,1,1))>1e-10) then
!!$                   !if((ii.ne.zero .and.jj.ne.zero).and.(ii.ne.zero .and.kk.ne.zero).and.(kk.ne.zero .and.jj.ne.zero))then
!!$                   !print *,'test',(ii.ne.zero .and.jj.ne.zero).and.(ii.ne.zero .and.kk.ne.zero).and.(kk.ne.zero .and.jj.ne.zero),(ii.ne.zero .and.jj.ne.zero),(ii.ne.zero .and.kk.ne.zero),(kk.ne.zero .and.jj.ne.zero)
!!$                   if (abs(drdextraerr(idxres,nDesignMach,1,1))>1.e-9)then
!!$                      write(*,*) 'mach',iCell,jCell,kCell, &
!!$                           dRdExtraAdj(idxres,nDesignMach,1,1), dRdextraFD(idxres,nDesignMach,1,1), dRdextraErr(idxres,nDesignMach,1,1)
!!$                   endif
!!$                end if
!!$                if(abs(dRdExtraFD(idxres,nDesignAOA,1,1))>1e-12) then
!!$                   if (abs(drdextraerr(idxres,nDesignAoA,1,1)/dRdExtraFD(idxres,nDesignAoA,1,1))>1.e-2)then
!!$                      write(*,*) 'AOA',m,iCell,jCell,kCell, &
!!$                           dRdExtraAdj(idxres,nDesignAOA,1,1), dRdextraFD(idxres,nDesignAOA,1,1), dRdextraErr(idxres,nDesignAOA,1,1)/dRdExtraFD(idxres,nDesignAoA,1,1)
!!$                   endif
!!$                end if
!!$                if(abs(dRdExtraFD(idxres,nDesignSSA,1,1))>1e-12) then
!!$                   if (abs(drdextraerr(idxres,nDesignSSA,1,1)/dRdExtraFD(idxres,nDesignSSA,1,1))>1.e-2)then
!!$                      write(*,*) 'SSA',m,iCell,jCell,kCell, &
!!$                           dRdExtraAdj(idxres,nDesignSSA,1,1), dRdextraFD(idxres,nDesignSSA,1,1), dRdextraErr(idxres,nDesignSSA,1,1)/dRdExtraFD(idxres,nDesignSSA,1,1)
!!$                   endif
!!$                end if
!!$                if(abs(dRdExtraFD(idxres,nDesignRotX,1,1))>1e-12) then
!!$                   if (abs(drdextraerr(idxres,nDesignRotx,1,1)/dRdExtraFD(idxres,nDesignRotX,1,1))>1.e-2)then
!!$                      write(*,*) 'RotX',m,iCell,jCell,kCell, &
!!$                           dRdExtraAdj(idxres,nDesignRotx,1,1), dRdextraFD(idxres,nDesignRotx,1,1), dRdextraErr(idxres,nDesignRotx,1,1)/dRdExtraFD(idxres,nDesignRotX,1,1)!, dRdExtraAdj(idxres,nDesignRotx,1,1)/dRdextraFD(idxres,nDesignRotx,1,1),timeref
!!$                   endif
!!$                end if
!!$             end if
!!$          enddo
!!$       end do
!!$    end do
!!$ end do

!      dRdwErrRel_1(:, :, :, :) = abs(dRdwErr_q(:, :, :, :)/dRdw_q(:,:,:,:))
!      print *
!      print *, "max, min of dRdExtraErr =", maxval(dRdExtraErr(:,:,1,1)), minval(dRdExtraErr(:,:,1,1))
!      print *, "max, min of dRdwErrRel =", maxval(dRdwErrRel_q(:,:,:,1)), minval(dRdwErrRel_q(:,:,:,1))
!      print *

           !deallocate memory for FD
      deallocatedomains: do nn = 1,ndom
         print *,'domain',nn
         groundLevel = 1
         sps = 1
         call setPointersAdj(nn,1,sps)
         deallocate(flowDoms(nn,level,sps)%dwp)
         deallocate(flowDoms(nn,level,sps)%dwm)
         deallocate(flowDoms(nn,level,sps)%dwtemp)
      end do deallocatedomains


       close(unitM)
       close(unitAoA)
       close(unitSSA)
       close(unitRotx)
       close(unitRoty)
       close(unitRotz)
      ! Output format.

99    format(a,1x,i6)
!111   format(4I4, 3ES22.14)



    end subroutine verifydRdExtra
    
