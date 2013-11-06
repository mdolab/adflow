!
!     ******************************************************************
!     *                                                                *
!     * File:          computeTSDerivatives.f90                        *
!     * Author:        C.A.(Sandy) Mader, G. Kenway                    *
!     * Starting date: 11-25-2009                                      *
!     * Last modified: 11-26-2009                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine computeTSDerivatives(force, moment, liftIndex, coef0, dcdalpha, &
     dcdalphadot, dcdq, dcdqdot)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Computes the stability derivatives based on the time spectral  *
  !     * solution of a given mesh. Takes in the force coefficients at   *
  !     * all time instantces and computes the agregate parameters       *
  !     *                                                                *
  !     ******************************************************************
  !
  use communication     
  use inputPhysics      
  use inputTimeSpectral 
  use inputTSStabDeriv
  use flowvarrefstate   
  use monitor           
  use section           
  use inputMotion

  implicit none

  !
  !     Subroutine arguments.
  !
  real(kind=realType), dimension(3, nTimeIntervalsSpectral) :: force, moment
  real(kind=realType), dimension(8):: dcdq, dcdqdot
  real(kind=realType), dimension(8):: dcdalpha,dcdalphadot
  real(kind=realType), dimension(8):: Coef0
  integer(kind=intType) :: liftIndex

  ! Working Variables
  real(kind=realType), dimension(nTimeIntervalsSpectral, 8) :: baseCoef
  real(kind=realType), dimension(8) ::coef0dot
  real(kind=realType), dimension(nTimeIntervalsSpectral,8)::ResBaseCoef
  real(kind=realType), dimension(nTimeIntervalsSpectral)  :: intervalAlpha,intervalAlphadot
  real(kind=realType), dimension(nTimeIntervalsSpectral)  :: intervalMach,intervalMachdot
  real(kind=realType), dimension(nSections) :: t
  real(kind=realType) :: alpha, beta
  integer(kind=intType):: i,sps,nn
  !speed of sound: for normalization of q derivatives
  real(kind=realType)::a
  real(kind=realType) :: scaleDim, fact, factMoment
  ! Functions
  real(kind=realType),dimension(nTimeIntervalsSpectral)  :: dPhix, dPhiy, dphiz
  real(kind=realType),dimension(nTimeIntervalsSpectral)  :: dPhixdot, dPhiydot, dphizdot
  real(kind=realType)::derivativeRigidRotAngle, secondDerivativeRigidRotAngle
  real(kind=realType):: TSAlpha,TSAlphadot


  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !

  scaleDim = pRef/pInf

  fact = two/(gammaInf*pInf*MachCoef**2 &
       *surfaceRef*LRef**2*scaleDim)
  factMoment = fact/(lengthRef*LRef)


  if (TSqMode)then

     print *,'TS Q Mode code needs to be updated in computeTSDerivatives!'
     stop

     ! !q is pitch
     ! do sps =1,nTimeIntervalsSpectral
     !    !compute the time of this intervavc
     !    t = timeUnsteadyRestart

     !    if(equationMode == timeSpectral) then
     !       do nn=1,nSections
     !          t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
     !               /         (nTimeIntervalsSpectral*1.0)
     !       enddo
     !    endif

     !    ! Compute the time derivative of the rotation angles around the
     !    ! z-axis. i.e. compute q

     !    dphiZ(sps) = derivativeRigidRotAngle(degreePolZRot,   &
     !         coefPolZRot,     &
     !         degreeFourZRot,  &
     !         omegaFourZRot,   &
     !         cosCoefFourZRot, &
     !         sinCoefFourZRot, t)

     !    ! add in q_dot computation
     !    dphiZdot(sps) = secondDerivativeRigidRotAngle(degreePolZRot,   &
     !         coefPolZRot,     &
     !         degreeFourZRot,  &
     !         omegaFourZRot,   &
     !         cosCoefFourZRot, &
     !         sinCoefFourZRot, t)
     ! end do

     ! !now compute dCl/dq
     ! do i =1,8
     !    call computeLeastSquaresRegression(BaseCoef(:,i),dphiz,nTimeIntervalsSpectral,dcdq(i),coef0(i))
     ! end do

     ! ! now subtract off estimated cl,cmz and use remainder to compute 
     ! ! clqdot and cmzqdot.
     ! do i = 1,8
     !    do sps = 1,nTimeIntervalsSpectral
     !       ResBaseCoef(sps,i) = BaseCoef(sps,i)-(dcdq(i)*dphiz(sps)+Coef0(i))
     !    enddo
     ! enddo

     ! !now normalize the results...
     ! a  = sqrt(gammaInf*pInfDim/rhoInfDim)
     ! dcdq = dcdq*timeRef*2*(machGrid*a)/lengthRef

     ! !now compute dCl/dpdot
     ! do i = 1,8
     !    call computeLeastSquaresRegression(ResBaseCoef(:,i),dphizdot,nTimeIntervalsSpectral,dcdqdot(i),Coef0dot(i))
     ! enddo

  elseif(TSAlphaMode)then

     do sps=1,nTimeIntervalsSpectral

        !compute the time of this interval
        t = timeUnsteadyRestart

        if(equationMode == timeSpectral) then
           do nn=1,nSections
              t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                   /         (nTimeIntervalsSpectral*1.0)
           enddo
        endif

        intervalAlpha(sps) = TSAlpha(degreePolAlpha,   coefPolAlpha,       &
             degreeFourAlpha,  omegaFourAlpha,     &
             cosCoefFourAlpha, sinCoefFourAlpha, t)

        intervalAlphadot(sps) = TSAlphadot(degreePolAlpha,   coefPolAlpha,       &
             degreeFourAlpha,  omegaFourAlpha,     &
             cosCoefFourAlpha, sinCoefFourAlpha, t)

        call getDirAngle(velDirFreestream,liftDirection,liftIndex,alpha+intervalAlpha(sps), beta)

        BaseCoef(sps,1) = fact*(&
             force(1, sps)*liftDirection(1) + &
             force(2, sps)*liftDirection(2) + &
             force(3, sps)*liftDIrection(3))
        BaseCoef(sps,2) = fact*(&
             force(1, sps)*dragDirection(1) + &
             force(2, sps)*dragDirection(2) + &
             force(3, sps)*dragDIrection(3))
        BaseCoef(sps,3) = force(1, sps)*fact
        BaseCoef(sps,4) = force(2, sps)*fact
        BaseCoef(sps,5) = force(3, sps)*fact
        BaseCoef(sps,6) = moment(1, sps)*factMoment
        BaseCoef(sps,7) = moment(2, sps)*factMoment
        BaseCoef(sps,8) = moment(3, sps)*factMoment
     end do

     !now compute dCl/dalpha
     do i =1,8
        call computeLeastSquaresRegression(BaseCoef(:,i),intervalAlpha,nTimeIntervalsSpectral,dcdAlpha(i),coef0(i))
     end do
     
     ! now subtract off estimated cl,cmz and use remainder to compute 
     ! clalphadot and cmzalphadot.
     do i = 1,8
        do sps = 1,nTimeIntervalsSpectral
           ResBaseCoef(sps,i) = BaseCoef(sps,i)-(dcdalpha(i)*intervalAlpha(sps)+Coef0(i))
        enddo
     enddo
     
     !now compute dCi/dalphadot
     do i = 1,8
        call computeLeastSquaresRegression(ResBaseCoef(:,i),intervalAlphadot,nTimeIntervalsSpectral,dcdalphadot(i),Coef0dot(i))
     enddo
     
     a  = sqrt(gammaInf*pInfDim/rhoInfDim)
     dcdalphadot = dcdalphadot*2*(machGrid*a)/lengthRef

  else
     call terminate('computeTSDerivatives','Not a valid stability motion')
  endif

end subroutine computeTSDerivatives


