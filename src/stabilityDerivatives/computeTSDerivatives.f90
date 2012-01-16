  !
!     ******************************************************************
!     *                                                                *
!     * File:          computeTSdetrivatives.f90                       *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 09-15-2009                                      *
!     * Last modified: 02-14-2011                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine computeTSDerivatives(coef0,dcdalpha,dcdalphadot,dcdq,dcdqdot)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Computes an estimate for the stability derivatives from the    *
  !     * coefficient values of the Time Spectral Solution.              *
  !     *                                                                *
  !     ******************************************************************
  !
  use blockPointers
  use iteration
  use precision
  use inputMotion
  use inputPhysics
  use inputTimeSpectral
  use inputTSStabDeriv
  use flowvarrefstate !Timeref
  use monitor
  use section
  use communication !myid
  use costFunctions
  implicit none

  !
  !     Subroutine arguments.
  !

  !
  !     Local variables.
  !
!  real(kind=realType),dimension(nTimeIntervalsSpectral)  :: CL, CD,&
 !      CFx, CFy,CFz, CMx, CMy, CMz
  real(kind=realType), dimension(3) :: cfp, cfv, cmp, cmv
  real(kind=realType) ::  yplusMax
  integer :: ierr

  real(kind=realType),dimension(nTimeIntervalsSpectral)  :: dPhix,&
       dPhiy,dPhiz

  real(kind=realType),dimension(nTimeIntervalsSpectral)  :: intervalAlpha,intervalAlphadot
  real(kind=realType),dimension(nTimeIntervalsSpectral)  :: intervalMach,intervalMachdot

!  real(kind=realType),dimension(nTimeIntervalsSpectral)  :: resCL,&
!       resCD, resCFx, resCFy,resCFz, resCMx, resCMy, resCMz

  real(kind=realType),dimension(nTimeIntervalsSpectral)  :: dPhixdot,&
       dPhiydot,dPhizdot

  real(kind=realType),dimension(8)::dcdp,dcdpdot,dcdq,dcdqdot,dcdr,dcdrdot
  !real(kind=realType)::dcldp,dcddp,dcfxdp,dcfydp,dcfzdp,dcmxdp,dcmydp,dcmzdp
  !real(kind=realType)::dcldpdot,dcddpdot,dcfxdpdot,dcfydpdot,dcfzdpdot,dcmxdpdot,dcmydpdot,dcmzdpdot
  !real(kind=realType)::dcldq,dcddq,dcfxdq,dcfydq,dcfzdq,dcmxdq,dcmydq,dcmzdq
  !real(kind=realType)::dcldqdot,dcddqdot,dcfxdqdot,dcfydqdot,dcfzdqdot,dcmxdqdot,dcmydqdot,dcmzdqdot
  !real(kind=realType)::dcldr,dcddr,dcfxdr,dcfydr,dcfzdr,dcmxdr,dcmydr,dcmzdr
  !real(kind=realType)::dcldrdot,dcddrdot,dcfxdrdot,dcfydrdot,dcfzdrdot,dcmxdrdot,dcmydrdot,dcmzdrdot
  !real(kind=realType)::dcldp,dcldpdot,dcddp,dcddpdot,dcmzdp,dcmzdpdot         
  !real(kind=realType)::dcldq,dcldqdot,dcddq,dcddqdot,dcmzdq,dcmzdqdot
  !real(kind=realType)::dcldr,dcldrdot,dcddr,dcddrdot,dcmzdr,dcmzdrdot
  real(kind=realType),dimension(8)::dcdalpha,dcdalphadot,dcdbeta,dcdbetadot,dcdMach,dcdMachdot
  !real(kind=realType)::dcldalpha,dcddalpha,dcfxdalpha,dcfydalpha,dcfzdalpha,dcmxdalpha,dcmydalpha,dcmzdalpha
  !real(kind=realType)::dcldalphadot,dcddalphadot,dcfxdalphadot,dcfydalphadot,dcfzdalphadot,dcmxdalphadot,dcmydalphadot,dcmzdalphadot
  !real(kind=realType)::dcldalpha,dcldalphadot,dcddalpha,dcddalphadot,dcmzdalpha,dcmzdalphadot
  !real(kind=realType)::dcldbeta,dcddbeta,dcfxdbeta,dcfydbeta,dcfzdbeta,dcmxdbeta,dcmydbeta,dcmzdbeta
  !real(kind=realType)::dcldbetadot,dcddbetadot,dcfxdbetadot,dcfydbetadot,dcfzdbetadot,dcmxdbetadot,dcmydbetadot,dcmzdbetadot
  !
  !real(kind=realType)::dcldMach,dcddMach,dcfxdMach,dcfydMach,dcfzdMach,dcmxdMach,dcmydMach,dcmzdMach
  !real(kind=realType)::dcldMachdot,dcddMachdot,dcfxdMachdot,dcfydMachdot,dcfzdMachdot,dcmxdMachdot,dcmydMachdot,dcmzdMachdot
  !real(kind=realType)::dcldMach,dcldMachdot,dcddMach,dcddMachdot,dcmzdMach,dcmzdMachdot
  real(kind=realType),dimension(8)::Coef0,Coef0dot
  !real(kind=realType)::cl0,cd0,cfx0,cfy0,cfz0,cmx0,cmy0,cmz0
  !real(kind=realType)::cl0dot,cd0dot,cfx0dot,cfy0dot,cfz0dot,cmx0dot,cmy0dot,cmz0dot
  !real(kind=realType)::cl0,cl0dot,cd0,cd0dot,cmz0,cmz0dot
  real(kind=realType),dimension(nTimeIntervalsSpectral,8)::BaseCoef,ResBaseCoef, BaseCoefLocal
  real(kind=realType),dimension(nTimeIntervalsSpectral):: storeAlpha
 
  real(kind=realType), dimension(nSections) :: t
  integer(kind=inttype)::level,sps,i,nn

  !function definitions
  real(kind=realType)::derivativeRigidRotAngle,&
       secondDerivativeRigidRotAngle
  !real(kind=realType)::t
  !real(kind=realType):: TSAlpha,TSAlphadot
  !real(kind=realType):: TSMach,TSMachdot
  real(kind=realType), dimension(nCostFunction)::globalCFVals

  ! updated lift and drag directions
  real(kind=realType), dimension(3)::velDir,liftDirTS,DragDirTS

  !Rotation variables
  real(kind=realType), dimension(3)   :: rotationPoint
  real(kind=realType), dimension(3,3) :: rotationMatrix,&
       derivRotationMatrix

  !lift direction variables
  integer(kind=intType) :: liftIndex
  real(kind=realType) :: alpha,beta,alphaTS,alphaIncrement,&
       betaTS,betaIncrement

  !Time Variables
  real(kind=realType) :: tNew, tOld

  !speed of sound: for normalization of q derivatives
  real(kind=realType)::a
  !Function Definitions
       
  real(kind=realType) :: TSAlpha,TSBeta,TSMach,TSAlphadot,TSMachdot
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !

  if(myID==0) print *,'in compute TS deriv...',nTimeintervalsSpectral
  !Compute and store the aero coef. Values for each TS level


  do sps =1,nTimeIntervalsSpectral
     
!!$     level = 1
!!$     call computeAeroCoef(globalCFVals,sps)
!!$
!!$     BaseCoef(sps,1) = globalCFVals(costFuncLiftCoef)
!!$     BaseCoef(sps,2) = globalCFVals(costFuncDragCoef)
!!$     BaseCoef(sps,3) = globalCFVals(costFuncForceXCoef)
!!$     BaseCoef(sps,4) = globalCFVals(costFuncForceYCoef)
!!$     BaseCoef(sps,5) = globalCFVals(costFuncForceZCoef)
!!$     BaseCoef(sps,6) = globalCFVals(costFuncMomXCoef)
!!$     BaseCoef(sps,7) = globalCFVals(costFuncMomYCoef)
!!$     BaseCoef(sps,8) = globalCFVals(costFuncMomZCoef)

     !update the lift vector and drag vector to account for changing 
     !angles of attack
     
     ! compute the time of this interval
     t = timeUnsteadyRestart
     
     if(equationMode == timeSpectral) then
        do nn=1,nSections
           t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                /         real(nTimeIntervalsSpectral,realType)
        enddo
     endif
     ! Determine the time values of the old and new time level.
     ! It is assumed that the rigid body rotation of the mesh is only
     ! used when only 1 section is present.
     tNew = timeUnsteady + timeUnsteadyRestart
     tOld = tNew - t(1)
    
    
     if(TSpMode.or. TSqMode .or.TSrMode) then
        ! Compute the rotation matrix of the rigid body rotation as
        ! well as the rotation point; the latter may vary in time due
        ! to rigid body translation.
                
        call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)
        
        liftDirTS(1) = rotationMatrix(1,1)*liftDirection(1) &
                  + rotationMatrix(1,2)*liftDirection(2) &
                  + rotationMatrix(1,3)*liftDirection(3)
        liftDirTS(2) = rotationMatrix(2,1)*liftDirection(1) &
                  + rotationMatrix(2,2)*liftDirection(2) &
                  + rotationMatrix(2,3)*liftDirection(3)
        liftDirTS(3) = rotationMatrix(3,1)*liftDirection(1) &
                  + rotationMatrix(3,2)*liftDirection(2) &
                  + rotationMatrix(3,3)*liftDirection(3)
        dragDirTS(1) = rotationMatrix(1,1)*dragDirection(1) &
                  + rotationMatrix(1,2)*dragDirection(2) &
                  + rotationMatrix(1,3)*dragDirection(3)
        dragDirTS(2) = rotationMatrix(2,1)*dragDirection(1) &
                  + rotationMatrix(2,2)*dragDirection(2) &
                  + rotationMatrix(2,3)*dragDirection(3)
        dragDirTS(3) = rotationMatrix(3,1)*dragDirection(1) &
                  + rotationMatrix(3,2)*dragDirection(2) &
                  + rotationMatrix(3,3)*dragDirection(3)
     elseif(tsAlphaMode)then
        ! get the baseline alpha and determine the liftIndex
        call getDirAngle(velDirFreestream,liftDirection,liftIndex,alpha,beta)
             !Determine the alpha for this time instance
        alphaIncrement = TSAlpha(degreePolAlpha,   coefPolAlpha,       &
             degreeFourAlpha,  omegaFourAlpha,     &
             cosCoefFourAlpha, sinCoefFourAlpha, t(1))
        
        alphaTS = alpha+alphaIncrement
        !Determine the grid velocity for this alpha
        call adjustInflowAngleAdj(alphaTS,beta,velDir,liftDirTS,dragDirTS,&
                  liftIndex)
        !do I need to update the lift direction and drag direction as well? yes!!!
        !if (myID==0) print *,'LiftDirTS',liftDirTS
        storeAlpha(sps) = alphaTS
        
     elseif(tsBetaMode)then
        ! get the baseline alpha and determine the liftIndex
        call getDirAngle(velDirFreestream,liftDirection,liftIndex,alpha,beta)
        
        !Determine the alpha for this time instance
        betaIncrement = TSBeta(degreePolBeta,   coefPolBeta,       &
             degreeFourBeta,  omegaFourBeta,     &
             cosCoefFourBeta, sinCoefFourBeta, t(1))
        
        betaTS = beta+betaIncrement
        !Determine the grid velocity for this alpha
        call adjustInflowAngleAdj(alpha,betaTS,velDir,liftDirTS,dragDirTS,&
             liftIndex)
       
     elseif(TSMachMode)then
        !determine the mach number at this time interval
        IntervalMach = TSMach(degreePolMach,   coefPolMach,       &
             degreeFourMach,  omegaFourMach,     &
             cosCoefFourMach, sinCoefFourMach, t(1))
        !no update needed
        
     elseif(TSAltitudeMode)then
        call terminate('gridVelocityFineLevel','altitude motion not yet implemented...')
     else
        call terminate('gridVelocityFineLevel','Not a recognized Stability Motion')
     end if
  
     ! Initialize the local monitoring variables to zero.
     !print *,'zeroing',sps,groundlevel
     groundlevel = 1
     BaseCoeflocal(sps,:) = zero
     BaseCoef(sps,:) = zero
     ! Loop over the blocks.
     
     domains: do nn=1,nDom
        !print *,'in domain',nn
        ! Set the pointers for this block.
        
        call setPointers(nn, groundLevel, sps)
        
        ! Compute the forces and moments for this block.
        
        call forcesAndMoments(cfp, cfv, cmp, cmv, yplusMax)
        BaseCoeflocal(sps,1) = BaseCoeflocal(sps,1)                       &
                          + (cfp(1) + cfv(1))*liftDirTS(1) &
                          + (cfp(2) + cfv(2))*liftDirTS(2) &
                          + (cfp(3) + cfv(3))*liftDirTS(3)

        BaseCoeflocal(sps,2) = BaseCoeflocal(sps,2)                        &
                          + (cfp(1) + cfv(1))*dragDirTS(1) &
                          + (cfp(2) + cfv(2))*dragDirTS(2) &
                          + (cfp(3) + cfv(3))*dragDirTS(3)
        
        BaseCoeflocal(sps,3) = BaseCoeflocal(sps,3) + cfp(1) + cfv(1)

        BaseCoeflocal(sps,4) = BaseCoeflocal(sps,4) + cfp(2) + cfv(2)

        BaseCoeflocal(sps,5) = BaseCoeflocal(sps,5) + cfp(3) + cfv(3)

        BaseCoeflocal(sps,6) = BaseCoeflocal(sps,6) + cmp(1) + cmv(1)

        BaseCoeflocal(sps,7) = BaseCoeflocal(sps,7) + cmp(2) + cmv(2)

        BaseCoeflocal(sps,8) = BaseCoeflocal(sps,8) + cmp(3) + cmv(3)
     enddo domains

     ! Determine the global sum of the coefficient
     ! variables. 
     
     call mpi_allreduce(BaseCoefLocal(sps,:), BaseCoef(sps,:), 8, sumb_real, &
                           mpi_sum, SUmb_comm_world, ierr)

  end do

  if (myID==0)then
     print *,'Lift coefficients: ', BaseCoef(:,1)
     print *,'Drag coefficients: ', BaseCoef(:,2)
     print *,'Alphas: ',storeAlpha
  end if

  if (TSpMode)then
     !P is roll
     do sps =1,nTimeIntervalsSpectral
        !compute the time of this interval
        t = timeUnsteadyRestart

        if(equationMode == timeSpectral) then
           do nn=1,nSections
              t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                   /         real(nTimeIntervalsSpectral,realType)
           enddo
        endif
        !print *,'t',t
        ! Compute the time derivative of the rotation angles around the
        ! x-axis. i.e. compute p

        dphiX(sps) = derivativeRigidRotAngle(degreePolXRot,   &
             coefPolXRot,     &
             degreeFourXRot,  &
             omegaFourXRot,   &
             cosCoefFourXRot, &
             sinCoefFourXRot, t)
        !if(myID==0)print *,'dphix',dphix

        !add in pdot computation here!
        dphiXdot(sps) = secondDerivativeRigidRotAngle(degreePolXRot,   &
             coefPolXRot,     &
             degreeFourXRot,  &
             omegaFourXRot,   &
             cosCoefFourXRot, &
             sinCoefFourXRot, t)
        !if(myID==0)print *,'dphixdot',dphixdot
     enddo
     !if(myID==0)print *,'dphix',dphix
     !if(myID==0)print *,'dphixdot',dphixdot
     !now compute dCl/dp
!!$     if (myid==0) then
!!$        print *,'CL',CL
!!$        print *,'dphix',dphix
!!$     endif

     do i =1,8
        call computeLeastSquaresRegression(BaseCoef(:,i),dphix,nTimeIntervalsSpectral,dcdp(i),coef0(i))
     end do
!     call computeLeastSquaresRegression(cl,dphix,nTimeIntervalsSpectral,dcldp,cl0)

!     !now compute dCd/dp

!     call computeLeastSquaresRegression(cd,dphix,nTimeIntervalsSpectral,dcddp,cd0)

 !    !now compute dCmz/dp

 !    call computeLeastSquaresRegression(cmz,dphix,nTimeIntervalsSpectral,dcmzdp,cmz0)
     if(myID==0)then
        print *,'CL estimates:'
        print *,'Clp = : ',dcdp(1),' cl0 = : ',Coef0(1)
        print *,'CD estimates:'
        print *,'CDp = : ',dcdp(2),' cd0 = : ',Coef0(2)
        print *,'CMz estimates:'
        print *,'CMzp = : ',dcdp(8),' cmz0 = : ',Coef0(8)
     endif

     ! now subtract off estimated cl,cmz and use remainder to compute 
     ! clpdot and cmzpdot.
     do i = 1,8
        do sps = 1,nTimeIntervalsSpectral
           ResBaseCoef(sps,i) = BaseCoef(sps,i)-(dcdp(i)*dphix(sps)+Coef0(i))
        enddo
        !ResCL(i) = Cl(i)-(dcldp*dphix(i)+cl0)
        !ResCD(i) = Cd(i)-(dcddp*dphix(i)+cd0)
        !ResCMz(i) = Cmz(i)-(dcmzdp*dphix(i)+cmz0)
     enddo
     if(myID==0) print *,'ResCL',ResBaseCoef(:,1)!resCL
     if(myID==0) print *,'ResCd',ResBaseCoef(:,2)!resCD
     if(myID==0) print *,'resCMZ',ResBaseCoef(:,8)!resCMZ

     !now compute dCl/dpdot
     do i = 1,8
        call computeLeastSquaresRegression(ResBaseCoef(:,i),dphixdot,nTimeIntervalsSpectral,dcdpdot(i),Coef0dot(i))
     enddo
     !call computeLeastSquaresRegression(Rescl,dphixdot,nTimeIntervalsSpectral,dcldpdot,cl0dot)

     !now compute dCd/dpdot

     !call computeLeastSquaresRegression(Rescd,dphixdot,nTimeIntervalsSpectral,dcddpdot,cl0dot)

     !now compute dCmz/dpdot

     !call computeLeastSquaresRegression(Rescmz,dphixdot,nTimeIntervalsSpectral,dcmzdpdot,cmz0dot)

     if(myID==0)then
        print *,'CL estimates:'
        print *,'Clpdot = : ',dcdpdot(1),' cl0dot = : ',coef0dot(1)
        print *,'Cd estimates:'
        print *,'Cdpdot = : ',dcdpdot(2),' cd0dot = : ',coef0dot(2)
        print *,'CMz estimates:'
        print *,'CMzpdot = : ',dcdpdot(8),' cmz0dot = : ',coef0dot(8)
     end if

  elseif(TSqMode)then
     !q is pitch
     do sps =1,nTimeIntervalsSpectral
        !compute the time of this interval
        t = timeUnsteadyRestart

        if(equationMode == timeSpectral) then
           do nn=1,nSections
              t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                   /         real(nTimeIntervalsSpectral,realType)
           enddo
        endif
        !if (myid==0)  print *,'t',t

        ! Compute the time derivative of the rotation angles around the
        ! z-axis. i.e. compute q

        dphiZ(sps) = derivativeRigidRotAngle(degreePolZRot,   &
             coefPolZRot,     &
             degreeFourZRot,  &
             omegaFourZRot,   &
             cosCoefFourZRot, &
             sinCoefFourZRot, t)
        !if(myID==0)print *,'dphiz',dphiz

        ! add in q_dot computation
        dphiZdot(sps) = secondDerivativeRigidRotAngle(degreePolZRot,   &
             coefPolZRot,     &
             degreeFourZRot,  &
             omegaFourZRot,   &
             cosCoefFourZRot, &
             sinCoefFourZRot, t)
        !if(myID==0)print *,'dphizdot',dphizdot
     end do
     !print *,'MYID', myID
!!$     if(myID==0)print *,'dphiz',dphiz
!!$     if(myID==0)print *,'dphizdot',dphizdot
!!$     if(myID==0)print *,'BaseCoef',BaseCoef

     do i =1,8
        call computeLeastSquaresRegression(BaseCoef(:,i),dphiz,nTimeIntervalsSpectral,dcdq(i),coef0(i))
     end do

     !if(myID==0)print *,'dcdq',dcdq
     ! now subtract off estimated cl,cmz and use remainder to compute 
     ! clqdot and cmzqdot.
     do i = 1,8
        do sps = 1,nTimeIntervalsSpectral
           ResBaseCoef(sps,i) = BaseCoef(sps,i)-(dcdq(i)*dphiz(sps)+Coef0(i))
        enddo
        !ResCL(i) = Cl(i)-(dcldp*dphix(i)+cl0)
        !ResCD(i) = Cd(i)-(dcddp*dphix(i)+cd0)
        !ResCMz(i) = Cmz(i)-(dcmzdp*dphix(i)+cmz0)
     enddo
     
!!$     if(myID==0)then
!!$        a  = sqrt(gammaInf*pInfDim/rhoInfDim)
!!$        print *,'before normalization',timeRef,(machGrid*a),lengthRef
!!$        print *,'CL estimates:'
!!$        print *,'Clq = : ',dcdq(1),' cl0 = : ',coef0(1)
!!$        print *,'CD estimates:'
!!$        print *,'Cdq = : ',dcdq(2),' cd0 = : ',coef0(2)
!!$        print *,'CMz estimates:'
!!$        print *,'CMzq = : ',dcdq(8),' cmz0 = : ',coef0(8)
!!$     endif
     !now normalize the results...
     
     a  = sqrt(gammaInf*pInfDim/rhoInfDim)
     dcdq = dcdq*timeRef*2*(machGrid*a)/lengthRef
!!$     if(myID==0)print *,'normalization',dcdq,timeRef,(machGrid*a),lengthRef
!!$     if(myID==0)then
!!$        !a  = sqrt(gammaInf*pInfDim/rhoInfDim)
!!$        print *,'normalization',timeRef,(machGrid*a),lengthRef
!!$        print *,'CL estimates:'
!!$        print *,'Clq = : ',dcdq(1),' cl0 = : ',coef0(1)
!!$        print *,'CD estimates:'
!!$        print *,'Cdq = : ',dcdq(2),' cd0 = : ',coef0(2)
!!$        print *,'CMz estimates:'
!!$        print *,'CMzq = : ',dcdq(8),' cmz0 = : ',coef0(8)
!!$     endif

    
!!$     if(myID==0) print *,'ResCL',ResBaseCoef(:,1)!resCL
!!$     if(myID==0) print *,'ResCd',ResBaseCoef(:,2)!resCD
!!$     if(myID==0) print *,'resCMZ',ResBaseCoef(:,8)!resCMZ

     !now compute dCl/dpdot
     do i = 1,8
        call computeLeastSquaresRegression(ResBaseCoef(:,i),dphizdot,nTimeIntervalsSpectral,dcdqdot(i),Coef0dot(i))
     enddo

!!$
!!$     if(myID==0)then
!!$        print *,'CL estimates:'
!!$        print *,'Clqdot = : ',dcdqdot(1),' cl0dot = : ',coef0dot(1)
!!$        print *,'Cd estimates:'
!!$        print *,'Cdqdot = : ',dcdqdot(2),' cd0dot = : ',coef0dot(2)
!!$        print *,'CMz estimates:'
!!$        print *,'CMzqdot = : ',dcdqdot(8),' cmz0dot = : ',coef0dot(8)
!!$     end if

  elseif(TSrMode)then
     !r is yaw
     do sps =1,nTimeIntervalsSpectral
        !compute the time of this interval
        t = timeUnsteadyRestart

        if(equationMode == timeSpectral) then
           do nn=1,nSections
              t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                   /         real(nTimeIntervalsSpectral,realType)
           enddo
        endif
        !print *,'t',t

        ! Compute the time derivative of the rotation angles around the
        ! y-axis and z-axis. i.e. compute r

        dphiY(sps) = derivativeRigidRotAngle(degreePolYRot,   &
             coefPolYRot,     &
             degreeFourYRot,  &
             omegaFourYRot,   &
             cosCoefFourYRot, &
             sinCoefFourYRot, t)

        !if(myID==0)print *,'dphiy',dphiy

        ! add in r_dot computation
        dphiYdot(sps) = secondDerivativeRigidRotAngle(degreePolYRot,   &
             coefPolYRot,     &
             degreeFourYRot,  &
             omegaFourYRot,   &
             cosCoefFourYRot, &
             sinCoefFourYRot, t)

        !if(myID==0)print *,'dphiydot',dphiydot

     enddo

     !if(myID==0)print *,'dphiy',dphiy
     !if(myID==0)print *,'dphiydot',dphiydot

!!$     if (myid==0) then
!!$        print *,'CL',CL
!!$        print *,'dphiy',dphiy
!!$     endif

     do i =1,8
        call computeLeastSquaresRegression(BaseCoef(:,i),dphiy,nTimeIntervalsSpectral,dcdr(i),coef0(i))
     end do
!!$     !now compute dCl/dr
!!$
!!$     call computeLeastSquaresRegression(cl,dphiy,nTimeIntervalsSpectral,dcldr,cl0)
!!$
!!$     !now compute dCD/dr
!!$
!!$     call computeLeastSquaresRegression(cd,dphiy,nTimeIntervalsSpectral,dcddr,cd0)
!!$
!!$     !now compute dCmz/dr
!!$
!!$     call computeLeastSquaresRegression(cmz,dphiy,nTimeIntervalsSpectral,dcmzdr,cmz0)
     if(myID==0)then
        print *,'CL estimates:'
        print *,'Clr = : ',dcdr(1),' cl0 = : ',coef0(1)
        print *,'Cd estimates:'
        print *,'Cdr = : ',dcdr(2),' cd0 = : ',coef0(2)
        print *,'CMz estimates:'
        print *,'CMzr = : ',dcdr(8),' cmz0 = : ',coef0(8)
     endif

     ! now subtract off estimated cl,cmz and use remainder to compute 
     ! clrdot and cmzrdot.
     do i = 1,8
        do sps = 1,nTimeIntervalsSpectral
           ResBaseCoef(sps,i) = BaseCoef(sps,i)-(dcdr(i)*dphiy(sps)+Coef0(i))
        enddo
        !ResCL(i) = Cl(i)-(dcldp*dphix(i)+cl0)
        !ResCD(i) = Cd(i)-(dcddp*dphix(i)+cd0)
        !ResCMz(i) = Cmz(i)-(dcmzdp*dphix(i)+cmz0)
     enddo
     if(myID==0) print *,'ResCL',ResBaseCoef(:,1)!resCL
     if(myID==0) print *,'ResCd',ResBaseCoef(:,2)!resCD
     if(myID==0) print *,'resCMZ',ResBaseCoef(:,8)!resCMZ

     !now compute dCi/drdot
     do i = 1,8
        call computeLeastSquaresRegression(ResBaseCoef(:,i),dphiydot,nTimeIntervalsSpectral,dcdrdot(i),Coef0dot(i))
     enddo
!!$     do i = 1,nTimeIntervalsSpectral
!!$        ResCL(i) = Cl(i)-(dcldr*dphiy(i)+cl0)
!!$        ResCd(i) = Cd(i)-(dcddr*dphiy(i)+cd0)
!!$        ResCMz(i) = Cmz(i)-(dcmzdr*dphiy(i)+cmz0)
!!$     enddo
!!$     if(myID==0) print *,'ResCL',resCL
!!$     if(myID==0) print *,'ResCd',resCd
!!$     if(myID==0) print *,'resCMZ',resCMZ


  elseif(TSAlphaMode)then

     !compute the alphas and alphadots
     do sps =1,nTimeIntervalsSpectral

        !compute the time of this interval
        t = timeUnsteadyRestart

        if(equationMode == timeSpectral) then
           do nn=1,nSections
              t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                   /         real(nTimeIntervalsSpectral,realType)
           enddo
        endif
        !print *,'t',t
        intervalAlpha(sps) = TSAlpha(degreePolAlpha,   coefPolAlpha,       &
             degreeFourAlpha,  omegaFourAlpha,     &
             cosCoefFourAlpha, sinCoefFourAlpha, t)

        intervalAlphadot(sps) = TSAlphadot(degreePolAlpha,   coefPolAlpha,       &
             degreeFourAlpha,  omegaFourAlpha,     &
             cosCoefFourAlpha, sinCoefFourAlpha, t)
     end do
     !if(myID==0)print *,'intervalAlpha',intervalAlpha
     do i =1,8
        call computeLeastSquaresRegression(BaseCoef(:,i),intervalAlpha,nTimeIntervalsSpectral,dcdAlpha(i),coef0(i))
     end do



!!$     !Derivatives are per radian!
!!$     if(myID==0)then
!!$        print *,'CL estimates:'
!!$        print *,'ClAlpha = : ',dcdalpha(1),' cl0 = : ',coef0(1)
!!$        print *,'Normal Force Coef: CFy'
!!$        print *,'CFyAlpha = : ',dcdalpha(4),' cFy0 = : ',coef0(4)
!!$        print *,'Cd estimates:'
!!$        print *,'CdAlpha = : ',dcdalpha(2),' cd0 = : ',coef0(2)
!!$        print *,'CMz estimates:'
!!$        print *,'CMzAlpha = : ',dcdalpha(8),' cmz0 = : ',coef0(8)
!!$     endif

     ! now subtract off estimated cl,cmz and use remainder to compute 
     ! clalphadot and cmzalphadot.
     do i = 1,8
        do sps = 1,nTimeIntervalsSpectral
           ResBaseCoef(sps,i) = BaseCoef(sps,i)-(dcdalpha(i)*intervalAlpha(sps)+Coef0(i))
        enddo
     enddo
!!$     if(myID==0) print *,'ResCL',ResBaseCoef(:,1)!resCL
!!$     if(myID==0) print *,'ResCd',ResBaseCoef(:,2)!resCD
!!$     if(myID==0) print *,'resCMZ',ResBaseCoef(:,8)!resCMZ

     !now compute dCi/dalphadot
     do i = 1,8
        call computeLeastSquaresRegression(ResBaseCoef(:,i),intervalAlphadot,nTimeIntervalsSpectral,dcdalphadot(i),Coef0dot(i))
     enddo

     a  = sqrt(gammaInf*pInfDim/rhoInfDim)
     dcdalphadot = dcdalphadot*2*(machGrid*a)/lengthRef
!!$     if(myID==0)then
!!$        print *,'normalization',timeRef,(machGrid*a),lengthRef
!!$        print *,'CL estimates:'
!!$        print *,'ClAlphadot = : ',dcdalphadot(1),' cl0dot = : ',coef0dot(1)
!!$        print *,'Normal Force Coef: CFy'
!!$        print *,'CFyAlphaDot = : ',dcdalphadot(4),' cFy0 = : ',coef0dot(4)
!!$        print *,'Cd estimates:'
!!$        print *,'CdAlphadot = : ',dcdalphadot(2),' cl0dot = : ',coef0dot(2)
!!$        print *,'CMz estimates:'
!!$        print *,'CMzAlphadot = : ',dcdalphadot(8),' cmz0dot = : ',coef0dot(8)
!!$     end if

  elseif(TSBetaMode)then
     print *,'Beta mode not yet implemented'
  elseif(TSMachMode)then

     !compute the Machs and Machdots
     do sps =1,nTimeIntervalsSpectral

        !compute the time of this interval
        t = timeUnsteadyRestart

        if(equationMode == timeSpectral) then
           do nn=1,nSections
              t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                   /         real(nTimeIntervalsSpectral,realType)
           enddo
        endif
        !if(myID==0)print *,'t',t
        intervalMach(sps) = machgrid+TSMach(degreePolMach,   coefPolMach,       &
             degreeFourMach,  omegaFourMach,     &
             cosCoefFourMach, sinCoefFourMach, t)

        intervalMachdot(sps) = TSMachdot(degreePolMach,   coefPolMach,       &
             degreeFourMach,  omegaFourMach,     &
             cosCoefFourMach, sinCoefFourMach, t)
     end do
     !if(myID==0) print *,'Mach', intervalMach,machgrid
     !now compute dCi/dMach
     do i =1,8
        call computeLeastSquaresRegression(BaseCoef(:,i),intervalMach,nTimeIntervalsSpectral,dcdMach(i),coef0(i))
     end do
!!$     call computeLeastSquaresRegression(cl,intervalMach,nTimeIntervalsSpectral,dcldMach,cl0)
!!$
!!$     !now compute dCd/dmach
!!$
!!$     call computeLeastSquaresRegression(cd,intervalMach,nTimeIntervalsSpectral,dcddMach,cd0)
!!$
!!$     !now compute dCmz/dmach
!!$
!!$     call computeLeastSquaresRegression(cmz,intervalMach,nTimeIntervalsSpectral,dcmzdMach,cmz0)
     if(myID==0)then
        print *,'CL estimates:'
        print *,'ClMach = : ',dcdMach(1),' cl0 = : ',coef0(1)
        print *,'Cd estimates:'
        print *,'CdMach = : ',dcdMach(2),' cl0 = : ',coef0(2)
        print *,'CMz estimates:'
        print *,'CMzMach = : ',dcdMach(8),' cmz0 = : ',coef0(8)
     endif

     ! now subtract off estimated cl,cmz and use remainder to compute 
     ! clMachdot and cmzMachdot.
     do i = 1,8
        do sps = 1,nTimeIntervalsSpectral
           ResBaseCoef(sps,i) = BaseCoef(sps,i)-(dcdMach(i)*intervalMach(sps)+Coef0(i))
        enddo
     enddo
     if(myID==0) print *,'ResCL',ResBaseCoef(:,1)!resCL
     if(myID==0) print *,'ResCd',ResBaseCoef(:,2)!resCD
     if(myID==0) print *,'resCMZ',ResBaseCoef(:,8)!resCMZ

     !now compute dCi/dMachdot
     do i = 1,8
        call computeLeastSquaresRegression(ResBaseCoef(:,i),intervalMachdot,nTimeIntervalsSpectral,dcdmachdot(i),Coef0dot(i))
     enddo
!!$     do i = 1,nTimeIntervalsSpectral
!!$        ResCL(i) = Cl(i)-(dcldalpha*IntervalMach(i)+cl0)
!!$        ResCd(i) = Cd(i)-(dcddalpha*IntervalMach(i)+cd0)
!!$        ResCMz(i) = Cmz(i)-(dcmzdalpha*IntervalMach(i)+cmz0)
!!$     enddo
!!$     if(myID==0) print *,'ResCL',resCL
!!$     if(myID==0) print *,'ResCd',resCd
!!$     if(myID==0) print *,'resCMZ',resCMZ
!!$     !now compute dCl/dqdot
!!$
!!$     call computeLeastSquaresRegression(Rescl,intervalMachdot,nTimeIntervalsSpectral,dcldMachdot,cl0dot)
!!$
!!$     !now compute dCd/dqdot
!!$
!!$     call computeLeastSquaresRegression(Rescd,intervalMachdot,nTimeIntervalsSpectral,dcddMachdot,cd0dot)
!!$
!!$     !now compute dCmz/dqdot
!!$
!!$     call computeLeastSquaresRegression(Rescmz,intervalAlphadot,nTimeIntervalsSpectral,dcmzdalphadot,cmz0dot)

     if(myID==0)then
        print *,'CL estimates:'
        print *,'ClMachdot = : ',dcdMachdot(1),' cl0dot = : ',coef0dot(1)
        print *,'Cd estimates:'
        print *,'CdMachdot = : ',dcdMachdot(2),' cd0dot = : ',coef0dot(2)
        print *,'CMz estimates:'
        print *,'CMzMachdot = : ',dcdMachdot(8),' cmz0dot = : ',coef0dot(8)
     end if

  else
     call terminate('computeTSDerivatives','Not a valid stability motion')
  endif
end subroutine computeTSDerivatives

