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
      subroutine verifydRdExtraFDAD(level)
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
      real(kind=realType) :: deltaw, ExtraAdjRef,dextra
      real(kind=realType) :: timeAdj, timeFD, timeResAdj,test

      integer(kind=intType), dimension(nDom) :: maxglobalcell
      integer(kind=intType) :: idx, ii, jj, kk, idxres, m, l,idxmgb,liftindex
      integer(kind=intType) :: sps, nTime, max_nTime, nHalo, nn, discr
      !real(kind=realType), allocatable, dimension(:,:,:,:) :: dRdExtraAdj,dRdExtraFD,dRdExtraErr

      integer :: ierr
      logical :: fineGrid, correctForK, exchangeTurb,secondHalo



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
101   format("FDMachfile",a,".out")
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
103   format("FDAoAfile",a,".out")
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
105   format("FDSSAfile",a,".out")
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
107   format("FDRotxfile",a,".out")
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
109   format("FDRotyfile",a,".out")
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
111   format("FDRotzfile",a,".out")
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
      !Mach Derivatives
      deltaw = 1.e-6
      !print *,'Entering Domain loop'
      domainLoopAD: do nn=1,nDom!97,97!1,1!nDom
         
         ! Loop over the number of time instances for this block.

         spectralLoop: do sps=1,nTimeIntervalsSpectral
  !          print *,'Setting Pointers'
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
                     !print *,'rotcenteradj',nn,icell,jcell,kcell,rotcenteradj
                     dwAdjb(:) = 0.
                     dwAdj(:)  = 0.
                     wAdjb(:,:,:,:)  = 0.  !dR(m)/dw
                     machadjb = 0
                     ExtraAdjRef = machadj

                     machadj = ExtraAdjRef+deltaw
                     
                     
                     call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
                          &  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
                          &  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
                          &  icell, jcell, kcell, nn, sps, correctfork, secondhalo, prefadj, &
                          &  rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj, &
                          &  rotrateadjb, rotcenteradj, murefadj, timerefadj, pinfcorradj, &
                          &  liftindex)
                     dwAdjP = dwadj

                     machadj = ExtraAdjRef-deltaw

                     call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
                          &  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
                          &  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
                          &  icell, jcell, kcell, nn, sps, correctfork, secondhalo, prefadj, &
                          &  rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj, &
                          &  rotrateadjb, rotcenteradj, murefadj, timerefadj, pinfcorradj, &
                          &  liftindex)

                     dwAdjM = dwadj
                     
                     mLoop: do m = 1, nw           ! Loop over output cell residuals (R)

                        ! Store the block Jacobians (by rows).
                        idxres   = globalCell(iCell,jCell,kCell)*nw+m
                        !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
                        if( idxres>=0) then
                           
                           dextra = (dwAdjP(m)- dwAdjM(m))/(2.0*deltaw)
                           !print *,'dwAdjp',dwAdjP(m),dwAdjM(m),dwAdjP(m)- dwAdjM(m),dextra
                           write(unitM,10) dextra,nn,icell,jcell,kcell,m,idxres
10                         format(1x,'Mach ',f18.10,6I8)

                        endif
                                       
                     end do mLoop
                     machadj = ExtraAdjRef
                    
                  end do !iCell
               end do !jCell
            end do! kCell
         end do spectralLoop
      end do domainLoopAD

      !AoA Derivatives
      deltaw = 1e-6
      !print *,'Entering Domain loop'
      domainLoopAD1: do nn=1,nDom!97,97!1,1!nDom
         
         ! Loop over the number of time instances for this block.

         spectralLoop1: do sps=1,nTimeIntervalsSpectral
  !          print *,'Setting Pointers'
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
                     !print *,'rotcenteradj',nn,icell,jcell,kcell,rotcenteradj
                     dwAdjb(:) = 0.
                     dwAdj(:)  = 0.
                     wAdjb(:,:,:,:)  = 0.  !dR(m)/dw
                     machadjb = 0
                     ExtraAdjRef = alphaadj

                     alphaadj = ExtraAdjRef+deltaw
                     
                     
                     call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
                          &  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
                          &  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
                          &  icell, jcell, kcell, nn, sps, correctfork, secondhalo, prefadj, &
                          &  rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj, &
                          &  rotrateadjb, rotcenteradj, murefadj, timerefadj, pinfcorradj, &
                          &  liftindex)

                     dwAdjP = dwadj

                     alphaadj = ExtraAdjRef-deltaw

                     call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
                          &  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
                          &  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
                          &  icell, jcell, kcell, nn, sps, correctfork, secondhalo, prefadj, &
                          &  rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj, &
                          &  rotrateadjb, rotcenteradj, murefadj, timerefadj, pinfcorradj, &
                          &  liftindex)

                     dwAdjM = dwadj
                     
                     mLoop1: do m = 1, nw           ! Loop over output cell residuals (R)

                        ! Store the block Jacobians (by rows).
                        idxres   = globalCell(iCell,jCell,kCell)*nw+m
                        !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
                        if( idxres>=0) then
                           dextra = (dwAdjP(m)- dwAdjM(m))/(2.0*deltaw)!2*deltaw

                           write(unitAoA,11) dextra,nn,icell,jcell,kcell,m,idxres
11                         format(1x,'AoA ',f18.10,6I8)

                        endif
                                       
                     end do mLoop1
                     alphaadj = ExtraAdjRef
                    
                  end do !iCell
               end do !jCell
            end do! kCell
         end do spectralLoop1
      end do domainLoopAD1


      !SSA Derivatives
      deltaw = 1e-6
      !print *,'Entering Domain loop'
      domainLoopAD2: do nn=1,nDom!97,97!1,1!nDom
         
         ! Loop over the number of time instances for this block.

         spectralLoop2: do sps=1,nTimeIntervalsSpectral
  !          print *,'Setting Pointers'
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
                     !print *,'rotcenteradj',nn,icell,jcell,kcell,rotcenteradj
                     dwAdjb(:) = 0.
                     dwAdj(:)  = 0.
                     wAdjb(:,:,:,:)  = 0.  !dR(m)/dw
                     machadjb = 0
                     ExtraAdjRef = betaadj

                     betaadj = ExtraAdjRef+deltaw
                     
                     
                     call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
                          &  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
                          &  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
                          &  icell, jcell, kcell, nn, sps, correctfork, secondhalo, prefadj, &
                          &  rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj, &
                          &  rotrateadjb, rotcenteradj, murefadj, timerefadj, pinfcorradj, &
                          &  liftindex)

                     dwAdjP = dwadj

                     betaadj = ExtraAdjRef-deltaw

                     call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
                          &  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
                          &  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
                          &  icell, jcell, kcell, nn, sps, correctfork, secondhalo, prefadj, &
                          &  rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj, &
                          &  rotrateadjb, rotcenteradj, murefadj, timerefadj, pinfcorradj, &
                          &  liftindex)

                     dwAdjM = dwadj
                     
                     mLoop2: do m = 1, nw           ! Loop over output cell residuals (R)

                        ! Store the block Jacobians (by rows).
                        idxres   = globalCell(iCell,jCell,kCell)*nw+m
                        !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
                        if( idxres>=0) then
                           dextra = (dwAdjP(m)- dwAdjM(m))/(2.0*deltaw)!2*deltaw

                           write(unitSSA,12) dextra,nn,icell,jcell,kcell,m,idxres
12                         format(1x,'SSA ',f18.10,6I8)



                        endif
                                       
                     end do mLoop2
                     betaadj = ExtraAdjRef
                    
                  end do !iCell
               end do !jCell
            end do! kCell
         end do spectralLoop2
      end do domainLoopAD2

      
      !rotx Derivatives
      deltaw = 1e-6
      !print *,'Entering Domain loop'
      domainLoopAD3: do nn=1,nDom!97,97!1,1!nDom
         
         ! Loop over the number of time instances for this block.

         spectralLoop3: do sps=1,nTimeIntervalsSpectral
  !          print *,'Setting Pointers'
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
                     !print *,'rotcenteradj',nn,icell,jcell,kcell,rotcenteradj
                     dwAdjb(:) = 0.
                     dwAdj(:)  = 0.
                     wAdjb(:,:,:,:)  = 0.  !dR(m)/dw
                     machadjb = 0
                     ExtraAdjRef = rotrateadj(1)

                     rotrateadj(1) = ExtraAdjRef+deltaw
                     
                     call computeRAdjoint(wAdj,xAdj,xBlockCornerAdj,dwAdj,alphaAdj,&
                          betaAdj,MachAdj, &
                          MachCoefAdj,machGridAdj,iCell, jCell,  kCell, &
                          nn,sps, correctForK,secondHalo,prefAdj,&
                          rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                          rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
                          murefAdj, timerefAdj,pInfCorrAdj,liftIndex)

                     dwAdjP = dwadj

                     rotrateadj(1) = ExtraAdjRef-deltaw

                     call computeRAdjoint(wAdj,xAdj,xBlockCornerAdj,dwAdj,alphaAdj,&
                          betaAdj,MachAdj, &
                          MachCoefAdj,machGridAdj,iCell, jCell,  kCell, &
                          nn,sps, correctForK,secondHalo,prefAdj,&
                          rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                          rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
                          murefAdj, timerefAdj,pInfCorrAdj,liftIndex)

                     dwAdjM = dwadj
                     
                     mLoop3: do m = 1, nw           ! Loop over output cell residuals (R)

                        ! Store the block Jacobians (by rows).
                        idxres   = globalCell(iCell,jCell,kCell)*nw+m
                        !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
                        if( idxres>=0) then
                           dextra = (dwAdjP(m)- dwAdjM(m))/(2.0*deltaw)!2*deltaw
                           write(unitRotx,13) dextra,nn,icell,jcell,kcell,m,idxres
13                         format(1x,'Rotx ',f18.10,6I8)

                        endif
                                       
                     end do mLoop3
!!$                     mLoop3: do m = 1, 3           ! Loop over output cell residuals (R)
!!$
!!$                        ! Store the block Jacobians (by rows).
!!$                        idxres   = globalCell(iCell,jCell,kCell)*3+m
!!$                        !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
!!$                        if( idxres>=0) then
!!$                           dextra = (dwAdjP(m)- dwAdjM(m))/(2.0*deltaw)!2*deltaw
!!$                           !if( m==2 )then
!!$                           !   dextra = nn
!!$                           !endif
!!$                           write(unitRotx,13) dextra,nn,icell,jcell,kcell,m,idxres
!!$13                         format(1x,'S ',f18.10,6I8)
!!$
!!$                        endif
!!$                                       
!!$                     end do mLoop3
                     rotrateadj(1) = ExtraAdjRef
                    
                  end do !iCell
               end do !jCell
            end do! kCell
         end do spectralLoop3
      end do domainLoopAD3
      
      !roty Derivatives
      deltaw = 1e-6
      !print *,'Entering Domain loop'
      domainLoopAD4: do nn=1,nDom!97,97!1,1!nDom
         
         ! Loop over the number of time instances for this block.

         spectralLoop4: do sps=1,nTimeIntervalsSpectral
  !          print *,'Setting Pointers'
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
                     !print *,'rotcenteradj',nn,icell,jcell,kcell,rotcenteradj
                     dwAdjb(:) = 0.
                     dwAdj(:)  = 0.
                     wAdjb(:,:,:,:)  = 0.  !dR(m)/dw
                     machadjb = 0
                     ExtraAdjRef = rotrateadj(2)

                     rotrateadj(2)= ExtraAdjRef+deltaw
                     
                     
                     call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
                          &  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
                          &  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
                          &  icell, jcell, kcell, nn, sps, correctfork, secondhalo, prefadj, &
                          &  rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj, &
                          &  rotrateadjb, rotcenteradj, murefadj, timerefadj, pinfcorradj, &
                          &  liftindex)

                     dwAdjP = dwadj

                     rotrateadj(2)= ExtraAdjRef-deltaw

                     call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
                          &  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
                          &  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
                          &  icell, jcell, kcell, nn, sps, correctfork, secondhalo, prefadj, &
                          &  rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj, &
                          &  rotrateadjb, rotcenteradj, murefadj, timerefadj, pinfcorradj, &
                          &  liftindex)

                     dwAdjM = dwadj
                     
                     mLoop4: do m = 1, nw           ! Loop over output cell residuals (R)

                        ! Store the block Jacobians (by rows).
                        idxres   = globalCell(iCell,jCell,kCell)*nw+m
                        !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
                        if( idxres>=0) then
                           dextra = (dwAdjP(m)- dwAdjM(m))/(2.0*deltaw)!2*deltaw

                           write(unitRoty,14) dextra,nn,icell,jcell,kcell,m,idxres
14                         format(1x,'Roty ',f18.10,6I8)

                        endif
                                       
                     end do mLoop4
                     rotrateadj(2)= ExtraAdjRef
                    
                  end do !iCell
               end do !jCell
            end do! kCell
         end do spectralLoop4
      end do domainLoopAD4

      !rotz Derivatives
      deltaw = 1e-6
      !print *,'Entering Domain loop'
      domainLoopAD5: do nn=1,nDom!97,97!1,1!nDom
         
         ! Loop over the number of time instances for this block.

         spectralLoop5: do sps=1,nTimeIntervalsSpectral
  !          print *,'Setting Pointers'
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
                     !print *,'rotcenteradj',nn,icell,jcell,kcell,rotcenteradj
                     dwAdjb(:) = 0.
                     dwAdj(:)  = 0.
                     wAdjb(:,:,:,:)  = 0.  !dR(m)/dw
                     machadjb = 0
                     ExtraAdjRef = rotrateadj(3)

                     rotrateadj(3) = ExtraAdjRef+deltaw
                     
                     
                     call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
                          &  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
                          &  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
                          &  icell, jcell, kcell, nn, sps, correctfork, secondhalo, prefadj, &
                          &  rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj, &
                          &  rotrateadjb, rotcenteradj, murefadj, timerefadj, pinfcorradj, &
                          &  liftindex)

                     dwAdjP = dwadj

                     rotrateadj(3) = ExtraAdjRef-deltaw

                     call COMPUTERADJOINT_B(wadj, wadjb, xadj, xadjb, xblockcorneradj, &
                          &  xblockcorneradjb, dwadj, dwadjb, alphaadj, alphaadjb, betaadj, &
                          &  betaadjb, machadj, machadjb, machcoefadj, machgridadj, machgridadjb, &
                          &  icell, jcell, kcell, nn, sps, correctfork, secondhalo, prefadj, &
                          &  rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, pinfadj, rotrateadj, &
                          &  rotrateadjb, rotcenteradj, murefadj, timerefadj, pinfcorradj, &
                          &  liftindex)

                     dwAdjM = dwadj
                     
                     mLoop5: do m = 1, nw           ! Loop over output cell residuals (R)

                        ! Store the block Jacobians (by rows).
                        idxres   = globalCell(iCell,jCell,kCell)*nw+m
                        !print *,'globalindices',idxstate,idxres,shape(dRdwAdj)
                        if( idxres>=0) then
                           dextra = (dwAdjP(m)- dwAdjM(m))/(2.0*deltaw)!2*deltaw

                           write(unitRotz,15) dextra,nn,icell,jcell,kcell,m,idxres
15                         format(1x,'Rotz ',f18.10,6I8)

                        endif
                                       
                     end do mLoop5
                     rotrateadj(3) = ExtraAdjRef
                    
                  end do !iCell
               end do !jCell
            end do! kCell
         end do spectralLoop5
      end do domainLoopAD5
      !print *,'AD Completed'
      ! Get new time and compute the elapsed AD time.
!      stop
      call mpi_barrier(SUmb_comm_world, ierr)
      if(myID == 0) then
         call cpu_time(time(2))
         timeAdj = time(2)-time(1)
      endif



       close(unitM)
       close(unitAoA)
       close(unitSSA)
       close(unitRotx)
       close(unitRoty)
       close(unitRotz)
      ! Output format.

99    format(a,1x,i6)
!111   format(4I4, 3ES22.14)



     end subroutine verifydRdExtraFDAD
    
