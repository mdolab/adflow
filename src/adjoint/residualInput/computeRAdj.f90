!
!      ******************************************************************
!      *                                                                *
!      * File:          computeRAdj.f90                                 *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 02-01-2008                                      *
!      * Last modified: 04-23-2008                                      *
!      *                                                                *
!      ******************************************************************
!

subroutine computeRAdjoint(wAdj,xAdj,xBlockCornerAdj,dwAdj,alphaAdj,&
                          betaAdj,MachAdj, &
                          MachCoefAdj,machGridAdj,iCell, jCell,  kCell, &
                          nn,level,sps, correctForK,secondHalo,prefAdj,&
                          rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                          rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
                          murefAdj, timerefAdj,pInfCorrAdj,liftIndex)
  
!      Set Use Modules
  use blockPointers
  use flowVarRefState
  use inputTimeSpectral !nTimeIntervalsSpectral
  use inputPhysics
  use section !nsection
  use monitor
  implicit none


!      Set Passed in Variables

  integer(kind=intType), intent(in) :: iCell, jCell, kCell,nn,level,sps
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral), &
       intent(in) :: wAdj
  real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral), &
       intent(in) :: xAdj

  real(kind=realType), dimension(nw,nTimeIntervalsSpectral)                :: dwAdj
  real(kind=realType), dimension(3),intent(in) ::rotRateAdj,rotCenterAdj
  real(kind=realType), dimension(2,2,2,3,nTimeIntervalsSpectral) ::xBlockCornerAdj

  logical :: secondHalo, correctForK,useOldCoor=.false.

!      Set Local Variables

  !variables for test loops
  integer(kind=intType)::i,j,k,ii,jj,kk,liftIndex,nnn,sps2
  integer(kind=intType) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

  real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral) :: pAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,3,nTimeIntervalsSpectral) :: sAdj
  real(kind=realType), dimension(nBocos,-2:2,-2:2,3,nTimeIntervalsSpectral) :: normAdj
  real(kind=realType), dimension(nBocos,-2:2,-2:2,nTimeIntervalsSpectral) ::rFaceAdj
  real(kind=realType),dimension(nTimeIntervalsSpectral):: volAdj
!  real(kind=realType), dimension(-2:2,-2:2,-2:2,3) :: siAdj, sjAdj, skAdj
  real(kind=realType), dimension(-1:1,-1:1,-1:1,nTimeIntervalsSpectral) :: radIAdj,radJAdj,radKAdj
  real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral) :: siAdj, sjAdj, skAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral) ::sFaceIAdj,sFaceJAdj,sFaceKAdj

  real(kind=realType), dimension(3) :: velDirFreestreamAdj
  real(kind=realType), dimension(3) :: liftDirectionAdj
  real(kind=realType), dimension(3) :: dragDirectionAdj
  real(kind=realType) :: MachAdj,MachCoefAdj,uInfAdj,pInfCorrAdj,machGridAdj
  real(kind=realType), dimension(nw)::wInfAdj 
  REAL(KIND=REALTYPE) :: prefAdj, rhorefAdj
  REAL(KIND=REALTYPE) :: pinfdimAdj, rhoinfdimAdj
  REAL(KIND=REALTYPE) :: rhoinfAdj, pinfAdj
  REAL(KIND=REALTYPE) :: murefAdj, timerefAdj

  real(kind=realType) :: alphaAdj, betaAdj

  real(kind=realType), dimension(nSections) :: t


!!$!File Parameters remove for AD
!!$      integer :: unitxAD = 15,ierror
!!$      integer ::iii,iiii,jjj,jjjj,kkk,kkkk,nnnn,istart2,jstart2,kstart2,iend2,jend2,kend2,n
!!$      character(len = 16)::outfile
!!$      
!!$      outfile = "xAD.txt"
!!$      
!!$      open (UNIT=unitxAD,File=outfile,status='old',position='append',action='write',iostat=ierror)
!!$      if(ierror /= 0)                        &
!!$           call terminate("verifyResiduals", &
!!$           "Something wrong when &
!!$           &calling open")
! *************************************************************************
!      Begin Execution
! *************************************************************************
  !print *,'in computeRadj',wadj(:,:,:,irho)!
!      call the initialization routines to calculate the effect of Mach and alpha
       call adjustInflowAngleAdj(alphaAdj,betaAdj,velDirFreestreamAdj,&
            liftDirectionAdj,dragDirectionAdj,liftIndex)
  
       call checkInputParamAdj(velDirFreestreamAdj,liftDirectionAdj,&
            dragDirectionAdj, Machadj, MachCoefAdj)

       call referenceStateAdj(Machadj, MachCoefAdj,uInfAdj,prefAdj,&
            rhorefAdj, pinfdimAdj, rhoinfdimAdj, rhoinfAdj, pinfAdj,&
            murefAdj, timerefAdj)
       !call referenceStateAdj(velDirFreestreamAdj,liftDirectionAdj,&
       !     dragDirectionAdj, Machadj, MachCoefAdj,uInfAdj,prefAdj,&
       !     rhorefAdj, pinfdimAdj, rhoinfdimAdj, rhoinfAdj, pinfAdj,&
       !     murefAdj, timerefAdj)
       !(velDirFreestreamAdj,liftDirectionAdj,&
       !     dragDirectionAdj, Machadj, MachCoefAdj,uInfAdj)

       call setFlowInfinityStateAdj(velDirFreestreamAdj,liftDirectionAdj,&
            dragDirectionAdj, Machadj, MachCoefAdj,uInfAdj,wInfAdj,prefAdj,&
            rhorefAdj, pinfdimAdj, rhoinfdimAdj, rhoinfAdj, pinfAdj,&
            murefAdj, timerefAdj,pInfCorrAdj)

       do sps2 = 1,nTimeIntervalsSpectral
!      Call the metric routines to generate the areas, volumes and surface normals for the stencil.
       !   print *,'nn xhalo',nn,secondhalo
       call xhaloAdj(xAdj,xBlockCornerAdj,icell,jcell,kcell,nn,level,sps,sps2)
       !print *,'nn metric',nn,secondhalo
       call metricAdj(xAdj,siAdj,sjAdj,skAdj,volAdj,normAdj, &
            iCell,jCell,kCell,nn,level,sps,sps2)
       !print *,'siAdj',siAdj(:,0,0,1,:)
!call the gridVelocities function to get the cell center ,face center and boundary mesh velocities.

       ! Compute the time, which corresponds to this spectral solution.
       ! For steady and unsteady mode this is simply the restart time;
       ! for the spectral mode the periodic time must be taken into
       ! account, which can be different for every section.
       
       t = timeUnsteadyRestart
       
       if(equationMode == timeSpectral) then
          do nnn=1,nSections
             t(nnn) = t(nnn) + (sps2-1)*sections(nnn)%timePeriod &
                   /         real(nTimeIntervalsSpectral,realType)
          enddo
       endif
       !first two arguments needed for time spectral.just set to initial values for the current steady case...
       !print *,'grid velocities',il,jl,kl,nn,secondhalo
       call gridVelocitiesFineLevelAdj(.false., t, sps,xAdj,&
            siAdj, sjAdj, skAdj,rotCenterAdj, rotRateAdj,sAdj,sFaceIAdj,&
            sFaceJAdj,sFaceKAdj,machGridAdj,velDirFreestreamAdj,&
            iCell, jCell, kCell,nn,level,sps2)
       !print *,'gvsiAdj',siAdj(:,0,0,1,:)
!for debugging intermediate stages      
!       dwadj(1:3) = sAdj(0,0,0,:)
!       return

       !print *,'normalVelocities',il,jl,kl,nn,secondhalo

       call normalVelocitiesAllLevelsAdj(sps,iCell, jCell, kCell,sFaceIAdj,&
            sFaceJAdj,sFaceKAdj,siAdj, sjAdj, skAdj,rFaceAdj,nn,level,sps2)

     ! print *,'nvsiAdj',siAdj(:,0,0,1,:)
 
       !needed for uSlip in Viscous Calculations
       !call slipVelocitiesFineLevel(.false., t, mm)

!      Mimic the Residual calculation in the main code

       !Compute the Pressure in the stencil based on the current 
       !States
       
       !print *,'Calling computepressure',il,jl,kl,nn,secondhalo!,wadj(:,:,:,irho)!
       ! replace with Compute Pressure Adjoint!
       call computePressureAdj(wAdj, pAdj,nn,level,sps,sps2)
!!$       !print out pAdj
!!$       istart2 = -2
!!$       jstart2 = -2
!!$       kstart2 = -2
!!$       iend2 = 2
!!$       jend2 = 2
!!$       kend2 = 2 
!!$       if(icell==2) istart2=-1
!!$       if(jcell==2) jstart2=-1
!!$       if(kcell==2) kstart2=-1
!!$       if(icell==il) iend2=1
!!$       if(jcell==jl) jend2=1
!!$       if(kcell==kl) kend2=1
!!$       do iiii = istart2,iend2
!!$          do jjjj = jstart2,jend2
!!$             do kkkk = kstart2,kend2
!!$                !do n = 1,3!nw
!!$                   do n = 1,1!nw
!!$                   !do n = 1,nw 
!!$                   !do sps2 = 1,nTimeIntervalsSpectral
!!$                      i = icell+iiii
!!$                      j = jcell+jjjj
!!$                      k = kcell+kkkk
!!$                                 !write(unitxAD,11) i,j,k,n,nnn,sps,sps2,wAdj(iiii,jjjj,kkkk,n,sps2) 
!!$                      !write(unitxAD,11) i,j,k,n,nn,sps,sps2,pAdj(iiii,jjjj,kkkk,sps2)
!!$                      write(unitxAD,11) i,j,k,n,nn,sps,sps2,sFaceIAdj(iiii,jjjj,kkkk,sps2)
!!$                      !write(unitxAD,11) i,j,k,n,nn,sps,sps2,sAdj(iiii,jjjj,kkkk,n,sps2)
!!$
!!$11                    format(1x,'wadj',7I8,f20.14)
!!$                   !enddo
!!$                enddo
!!$             enddo
!!$          enddo
!!$       enddo
       
       ! Apply all boundary conditions to stencil.
       ! In case of a full mg mode, and a segegated turbulent solver,
       ! first call the turbulent boundary conditions, such that the
       ! turbulent kinetic energy is properly initialized in the halo's.

!###! Ignore Viscous for now
!###!       if(turbSegregated .and. (.not. corrections)) &
!###!         call applyAllTurbBCAdj(secondHalo)

       ! Apply all boundary conditions of the mean flow.
       !print *,'applying bcs',nn,secondhalo,sps2
!******************************************
       call applyAllBCAdj(wInfAdj,pInfCorrAdj,wAdj, pAdj,sAdj, &
            siAdj, sjAdj, skAdj, volAdj, normAdj, &
            rFaceAdj,iCell, jCell, kCell,secondHalo,nn,level,sps,sps2)
       !print *,'bcsiAdj',siAdj(:,0,0,1,:)
!#!#Shouldn't need this section for derivatives...
!#!$       ! In case this routine is called in full mg mode call the mean
!#!$       ! flow boundary conditions again such that the normal momentum
!#!$       ! boundary condition is treated correctly.
!#!$
!#!$       if(.not. corrections) call applyAllBCAdj(wAdj, pAdj, &
!#!$                              siAdj, sjAdj, skAdj, volAdj, normAdj, &
!#!$                              iCell, jCell, kCell,secondHalo)

       !Leave out State exchanges for now. If there are discrepancies 
       !Later, this may be a source...
!#!$       ! Exchange the solution. Either whalo1 or whalo2
!#!$       ! must be called.
!#!$
!#!$       if( secondHalo ) then
!#!$         call whalo2(currentLevel, 1_intType, nVarInt, .true., &
!#!$                     .true., .true.)
!#!$       else
!#!$         call whalo1(currentLevel, 1_intType, nVarInt, .true., &
!#!$                     .true., .true.)
!#!$       endif

!Again this should not be required, so leave out for now...
       ! For full multigrid mode the bleeds must be determined, the
       ! boundary conditions must be applied one more time and the
       ! solution must be exchanged again.

!#!$       if(.not. corrections) then
!#!$         call BCDataMassBleedOutflowAdj(.true., .true.)
!#!$         call applyAllBCAdj(secondHalo)
!#!$
!#!$       !Leave out State exchanges for now. If there are discrepancies 
!#!$       !Later, this may be a source...

!#!$!         if( secondHalo ) then
!#!$!           call whalo2(currentLevel, 1_intType, nVarInt, .true., &
!#!$!                       .true., .true.)
!#!$!         else
!#!!$           call whalo1(currentLevel, 1_intType, nVarInt, .true., &
!#!!$                       .true., .true.)
!#!!$         endif
!#!$       endif
!#!$
!#!$
!#!$
!#!$       ! Reset the values of rkStage and currentLevel, such that
!#!$       ! they correspond to a new iteration.
!#!$
!#!$       rkStage = 0
!#!$       currentLevel = groundLevel
!#!$
!#!$       ! Compute the latest values of the skin friction velocity.
!#!$       ! The currently stored values are of the previous iteration.
!#!$
!#!$       call computeUtauAdj
!#!$
!#!$       ! Apply an iteration to the turbulent transport equations in
!#!$       ! case these must be solved segregatedly.
!#!$
!#!$       if( turbSegregated ) call turbSolveSegregatedAdj
!#!$
       ! Compute the time step.

       !call timeStepAdj(.false.)
       !print *,'nntimestep',nn
       call timeStepAdj(.true.,wAdj,pAdj,siAdj, sjAdj, skAdj,&
            sFaceIAdj,sFaceJAdj,sFaceKAdj,volAdj,radIAdj,radJAdj,radKAdj,&
            iCell, jCell, kCell,pInfCorrAdj,rhoInfAdj,nn,level,sps,sps2)
       !print *,'tstepsiAdj',siAdj(:,0,0,1,:)
    end do
!#!$
!#!$       ! Compute the residual of the new solution on the ground level.
!#!$
!#!$       if( turbCoupled ) then
!#!$         call initresAdj(nt1MG, nMGVar)
!#!$         call turbResidualAdj
!#!$       endif
!#!$
     !do sps2 = 1,nTimeIntervalsSpectral
     !  print *,'calculating initres',nn
       !call initresAdj(1_intType, nwf,sps,dwAdj)
    !dwAdj(:,sps) = 0.0
       call initresAdj(1, nwf,wAdj,volAdj,dwAdj,nn,level,sps)
       !print *,'dwadj',dwadj,icell,jcell,kcell
       
    
     !  print *,'calculating residuals',nn
       call residualAdj(wAdj,pAdj,siAdj,sjAdj,skAdj,volAdj,normAdj,&
                              sFaceIAdj,sFaceJAdj,sFaceKAdj,&
                              radIAdj,radJAdj,radKAdj,&
                              dwAdj, iCell, jCell, kCell,  &  
                              rotRateAdj,correctForK,nn,level,sps)
    !end do
      ! print *,'nn end',nn
!stop
!close (UNIT=unitxAD)
     end subroutine computeRAdjoint
