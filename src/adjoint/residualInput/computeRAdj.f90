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

subroutine computeRAdjoint(wAdj,xAdj,dwAdj,alphaAdj,betaAdj,MachAdj, &
                          MachCoefAdj,iCell, jCell,  kCell, &
                          nn,sps, correctForK,secondHalo,prefAdj,&
                          rhorefAdj, pinfdimAdj, rhoinfdimAdj,&
                          rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,&
                          murefAdj, timerefAdj,pInfCorrAdj,liftIndex)
  
!      Set Use Modules
  use blockPointers
  use flowVarRefState



!      Set Passed in Variables

  integer(kind=intType), intent(in) :: iCell, jCell, kCell,nn,sps
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw), &
       intent(in) :: wAdj
  real(kind=realType), dimension(-3:2,-3:2,-3:2,3), &
       intent(in) :: xAdj

  real(kind=realType), dimension(nw)                :: dwAdj
  real(kind=realType), dimension(3),intent(in) ::rotRateAdj,rotCenterAdj

  logical :: secondHalo, correctForK,useOldCoor=.false.

!      Set Local Variables

  !variables for test loops
  integer(kind=intType)::i,j,k,ii,jj,kk,liftIndex
  integer(kind=intType) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

  real(kind=realType), dimension(-2:2,-2:2,-2:2) :: pAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,3) :: sAdj
  real(kind=realType), dimension(nBocos,-2:2,-2:2,3) :: normAdj
  real(kind=realType), dimension(nBocos,-2:2,-2:2) ::rFaceAdj
  real(kind=realType):: volAdj
!  real(kind=realType), dimension(-2:2,-2:2,-2:2,3) :: siAdj, sjAdj, skAdj
  real(kind=realType), dimension(-3:2,-3:2,-3:2,3) :: siAdj, sjAdj, skAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2) ::sFaceIAdj,sFaceJAdj,sFaceKAdj

  real(kind=realType), dimension(3) :: velDirFreestreamAdj
  real(kind=realType), dimension(3) :: liftDirectionAdj
  real(kind=realType), dimension(3) :: dragDirectionAdj
  real(kind=realType) :: MachAdj,MachCoefAdj,uInfAdj,pInfCorrAdj
  real(kind=realType), dimension(nw)::wInfAdj 
  REAL(KIND=REALTYPE) :: prefAdj, rhorefAdj
  REAL(KIND=REALTYPE) :: pinfdimAdj, rhoinfdimAdj
  REAL(KIND=REALTYPE) :: rhoinfAdj, pinfAdj
  REAL(KIND=REALTYPE) :: murefAdj, timerefAdj

  real(kind=realType) :: alphaAdj, betaAdj


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


!      Call the metric routines to generate the areas, volumes and surface normals for the stencil.
       
       call metricAdj(xAdj,siAdj,sjAdj,skAdj,volAdj,normAdj, &
            iCell,jCell,kCell)
      
!call the gridVelocities function to get the cell center ,face center and boundary mesh velocities.

       !first two arguments needed for time spectral.just set to initial values for the current steady case...
       !print *,'grid velocities'
       call gridVelocitiesFineLevelAdj(.false., zero, sps,xAdj,&
            siAdj, sjAdj, skAdj,rotCenterAdj, rotRateAdj,sAdj,sFaceIAdj,&
            sFaceJAdj,sFaceKAdj, iCell, jCell, kCell)
       
       !print *,'normalVelocities'

       call normalVelocitiesAllLevelsAdj(sps,iCell, jCell, kCell,sFaceIAdj,&
            sFaceJAdj,sFaceKAdj,siAdj, sjAdj, skAdj,rFaceAdj)
 
       !needed for uSlip in Viscous Calculations
       !call slipVelocitiesFineLevel(.false., t, mm)

!      Mimic the Residual calculation in the main code

       !Compute the Pressure in the stencil based on the current 
       !States
       
       !print *,'Calling computepressure',wadj(:,:,:,irho)!
       ! replace with Compute Pressure Adjoint!
       call computePressureAdj(wAdj, pAdj)
       
      
       ! Apply all boundary conditions to stencil.
       ! In case of a full mg mode, and a segegated turbulent solver,
       ! first call the turbulent boundary conditions, such that the
       ! turbulent kinetic energy is properly initialized in the halo's.

!###! Ignore Viscous for now
!###!       if(turbSegregated .and. (.not. corrections)) &
!###!         call applyAllTurbBCAdj(secondHalo)

       ! Apply all boundary conditions of the mean flow.
      ! print *,'applying bcs'
!******************************************
       call applyAllBCAdj(wInfAdj,pInfCorrAdj,wAdj, pAdj,sAdj, &
            siAdj, sjAdj, skAdj, volAdj, normAdj, &
            rFaceAdj,iCell, jCell, kCell,secondHalo)

!!#Shouldn't need this section for derivatives...
!!$       ! In case this routine is called in full mg mode call the mean
!!$       ! flow boundary conditions again such that the normal momentum
!!$       ! boundary condition is treated correctly.
!!$
!!$       if(.not. corrections) call applyAllBCAdj(wAdj, pAdj, &
!!$                              siAdj, sjAdj, skAdj, volAdj, normAdj, &
!!$                              iCell, jCell, kCell,secondHalo)

       !Leave out State exchanges for now. If there are discrepancies 
       !Later, this may be a source...
!!$       ! Exchange the solution. Either whalo1 or whalo2
!!$       ! must be called.
!!$
!!$       if( secondHalo ) then
!!$         call whalo2(currentLevel, 1_intType, nVarInt, .true., &
!!$                     .true., .true.)
!!$       else
!!$         call whalo1(currentLevel, 1_intType, nVarInt, .true., &
!!$                     .true., .true.)
!!$       endif

!Again this should not be required, so leave out for now...
       ! For full multigrid mode the bleeds must be determined, the
       ! boundary conditions must be applied one more time and the
       ! solution must be exchanged again.

!!$       if(.not. corrections) then
!!$         call BCDataMassBleedOutflowAdj(.true., .true.)
!!$         call applyAllBCAdj(secondHalo)
!!$
!!$       !Leave out State exchanges for now. If there are discrepancies 
!!$       !Later, this may be a source...

!!$!         if( secondHalo ) then
!!$!           call whalo2(currentLevel, 1_intType, nVarInt, .true., &
!!$!                       .true., .true.)
!!$!         else
!!!$           call whalo1(currentLevel, 1_intType, nVarInt, .true., &
!!!$                       .true., .true.)
!!!$         endif
!!$       endif
!!$
!!$
!!$
!!$       ! Reset the values of rkStage and currentLevel, such that
!!$       ! they correspond to a new iteration.
!!$
!!$       rkStage = 0
!!$       currentLevel = groundLevel
!!$
!!$       ! Compute the latest values of the skin friction velocity.
!!$       ! The currently stored values are of the previous iteration.
!!$
!!$       call computeUtauAdj
!!$
!!$       ! Apply an iteration to the turbulent transport equations in
!!$       ! case these must be solved segregatedly.
!!$
!!$       if( turbSegregated ) call turbSolveSegregatedAdj
!!$
!!$       ! Compute the time step.
!!$
!!$       call timeStepAdj(.false.)
!!$
!!$       ! Compute the residual of the new solution on the ground level.
!!$
!!$       if( turbCoupled ) then
!!$         call initresAdj(nt1MG, nMGVar)
!!$         call turbResidualAdj
!!$       endif
!!$

       !print *,'calculating residuals'
       !call initresAdj(1_intType, nwf,sps,dwAdj)
       call initresAdj(1, nwf,sps,dwAdj)
       
       call residualAdj(wAdj,pAdj,siAdj,sjAdj,skAdj,volAdj,normAdj,&
                              sFaceIAdj,sFaceJAdj,sFaceKAdj,&
                              dwAdj, iCell, jCell, kCell,  &  
                              rotRateAdj,correctForK)

!stop

     end subroutine computeRAdjoint
