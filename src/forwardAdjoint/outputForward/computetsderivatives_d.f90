   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of computetsderivatives in forward (tangent) mode (with options i4 dr8 r8):
   !   variations   of useful results: dcdalphadot coef0 dcdalpha
   !   with respect to varying inputs: machgrid lengthref machcoef
   !                dragdirection liftdirection gammainf pinf rhoinfdim
   !                pinfdim pref moment force
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
   SUBROUTINE COMPUTETSDERIVATIVES_D(force, forced, moment, momentd, &
   & liftindex, coef0, coef0d, dcdalpha, dcdalphad, dcdalphadot, &
   & dcdalphadotd, dcdq, dcdqdot)
   !
   !     ******************************************************************
   !     *                                                                *
   !     * Computes the stability derivatives based on the time spectral  *
   !     * solution of a given mesh. Takes in the force coefficients at   *
   !     * all time instantces and computes the agregate parameters       *
   !     *                                                                *
   !     ******************************************************************
   !
   USE COMMUNICATION
   USE INPUTPHYSICS
   USE INPUTTIMESPECTRAL
   USE INPUTTSSTABDERIV
   USE FLOWVARREFSTATE
   USE MONITOR
   USE SECTION
   USE INPUTMOTION
   IMPLICIT NONE
   !
   !     Subroutine arguments.
   !
   REAL(kind=realtype), DIMENSION(3, ntimeintervalsspectral) :: force, &
   & moment
   REAL(kind=realtype), DIMENSION(3, ntimeintervalsspectral) :: forced, &
   & momentd
   REAL(kind=realtype), DIMENSION(8) :: dcdq, dcdqdot
   REAL(kind=realtype), DIMENSION(8) :: dcdalpha, dcdalphadot
   REAL(kind=realtype), DIMENSION(8) :: dcdalphad, dcdalphadotd
   REAL(kind=realtype), DIMENSION(8) :: coef0
   REAL(kind=realtype), DIMENSION(8) :: coef0d
   INTEGER(kind=inttype) :: liftindex
   ! Working Variables
   REAL(kind=realtype), DIMENSION(ntimeintervalsspectral, 8) :: basecoef
   REAL(kind=realtype), DIMENSION(ntimeintervalsspectral, 8) :: basecoefd
   REAL(kind=realtype), DIMENSION(8) :: coef0dot
   REAL(kind=realtype), DIMENSION(8) :: coef0dotd
   REAL(kind=realtype), DIMENSION(ntimeintervalsspectral, 8) :: &
   & resbasecoef
   REAL(kind=realtype), DIMENSION(ntimeintervalsspectral, 8) :: &
   & resbasecoefd
   REAL(kind=realtype), DIMENSION(ntimeintervalsspectral) :: &
   & intervalalpha, intervalalphadot
   REAL(kind=realtype), DIMENSION(ntimeintervalsspectral) :: intervalmach&
   & , intervalmachdot
   REAL(kind=realtype), DIMENSION(nsections) :: t
   REAL(kind=realtype) :: alpha, beta
   INTEGER(kind=inttype) :: i, sps, nn
   !speed of sound: for normalization of q derivatives
   REAL(kind=realtype) :: a
   REAL(kind=realtype) :: ad
   REAL(kind=realtype) :: scaledim, fact, factmoment
   REAL(kind=realtype) :: scaledimd, factd, factmomentd
   ! Functions
   REAL(kind=realtype), DIMENSION(ntimeintervalsspectral) :: dphix, dphiy&
   & , dphiz
   REAL(kind=realtype), DIMENSION(ntimeintervalsspectral) :: dphixdot, &
   & dphiydot, dphizdot
   REAL(kind=realtype) :: derivativerigidrotangle, &
   & secondderivativerigidrotangle
   REAL(kind=realtype) :: DERIVATIVERIGIDROTANGLE_D
   REAL(kind=realtype) :: TSALPHA, TSALPHADOT
   INTRINSIC SQRT
   REAL(kind=realtype) :: arg1
   REAL(kind=realtype) :: arg1d
   !
   !     ******************************************************************
   !     *                                                                *
   !     * Begin execution.                                               *
   !     *                                                                *
   !     ******************************************************************
   !
   scaledimd = (prefd*pinf-pref*pinfd)/pinf**2
   scaledim = pref/pinf
   factd = -(two*surfaceref*lref**2*((gammainfd*pinf+gammainf*pinfd)*&
   &   machcoef**2*scaledim+gammainf*pinf*(2*machcoef*machcoefd*scaledim+&
   &   machcoef**2*scaledimd))/(gammainf*pinf*machcoef**2*surfaceref*lref**&
   &   2*scaledim)**2)
   fact = two/(gammainf*pinf*machcoef**2*surfaceref*lref**2*scaledim)
   factmomentd = (factd*lengthref*lref-fact*lref*lengthrefd)/(lengthref*&
   &   lref)**2
   factmoment = fact/(lengthref*lref)
   CALL GETDIRANGLE(veldirfreestream, liftdirection, liftindex, alpha&
   &               , beta)
   IF (tsqmode) THEN
   PRINT*, &
   &   'TS Q Mode code needs to be updated in computeTSDerivatives!'
   STOP
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
   ELSE IF (tsalphamode) THEN
   basecoefd = 0.0_8
   DO sps=1,ntimeintervalsspectral
   !compute the time of this interval
   t = timeunsteadyrestart
   IF (equationmode .EQ. timespectral) THEN
   DO nn=1,nsections
   t(nn) = t(nn) + (sps-1)*sections(nn)%timeperiod/(&
   &           ntimeintervalsspectral*1.0)
   END DO
   END IF
   intervalalpha(sps) = TSALPHA(degreepolalpha, coefpolalpha, &
   &       degreefouralpha, omegafouralpha, coscoeffouralpha, &
   &       sincoeffouralpha, t)
   intervalalphadot(sps) = TSALPHADOT(degreepolalpha, coefpolalpha&
   &       , degreefouralpha, omegafouralpha, coscoeffouralpha, &
   &       sincoeffouralpha, t)
   CALL GETDIRANGLE(veldirfreestream, liftdirection, liftindex, &
   &                   alpha + intervalalpha(sps), beta)
   basecoefd(sps, 1) = factd*(force(1, sps)*liftdirection(1)+force(2&
   &       , sps)*liftdirection(2)+force(3, sps)*liftdirection(3)) + fact*(&
   &       forced(1, sps)*liftdirection(1)+force(1, sps)*liftdirectiond(1)+&
   &       forced(2, sps)*liftdirection(2)+force(2, sps)*liftdirectiond(2)+&
   &       forced(3, sps)*liftdirection(3)+force(3, sps)*liftdirectiond(3))
   basecoef(sps, 1) = fact*(force(1, sps)*liftdirection(1)+force(2, &
   &       sps)*liftdirection(2)+force(3, sps)*liftdirection(3))
   basecoefd(sps, 2) = factd*(force(1, sps)*dragdirection(1)+force(2&
   &       , sps)*dragdirection(2)+force(3, sps)*dragdirection(3)) + fact*(&
   &       forced(1, sps)*dragdirection(1)+force(1, sps)*dragdirectiond(1)+&
   &       forced(2, sps)*dragdirection(2)+force(2, sps)*dragdirectiond(2)+&
   &       forced(3, sps)*dragdirection(3)+force(3, sps)*dragdirectiond(3))
   basecoef(sps, 2) = fact*(force(1, sps)*dragdirection(1)+force(2, &
   &       sps)*dragdirection(2)+force(3, sps)*dragdirection(3))
   basecoefd(sps, 3) = forced(1, sps)*fact + force(1, sps)*factd
   basecoef(sps, 3) = force(1, sps)*fact
   basecoefd(sps, 4) = forced(2, sps)*fact + force(2, sps)*factd
   basecoef(sps, 4) = force(2, sps)*fact
   basecoefd(sps, 5) = forced(3, sps)*fact + force(3, sps)*factd
   basecoef(sps, 5) = force(3, sps)*fact
   basecoefd(sps, 6) = momentd(1, sps)*factmoment + moment(1, sps)*&
   &       factmomentd
   basecoef(sps, 6) = moment(1, sps)*factmoment
   basecoefd(sps, 7) = momentd(2, sps)*factmoment + moment(2, sps)*&
   &       factmomentd
   basecoef(sps, 7) = moment(2, sps)*factmoment
   basecoefd(sps, 8) = momentd(3, sps)*factmoment + moment(3, sps)*&
   &       factmomentd
   basecoef(sps, 8) = moment(3, sps)*factmoment
   END DO
   coef0d = 0.0_8
   dcdalphad = 0.0_8
   !now compute dCl/dalpha
   DO i=1,8
   CALL COMPUTELEASTSQUARESREGRESSION_D(basecoef(:, i), basecoefd(:, &
   &                                    i), intervalalpha, &
   &                                    ntimeintervalsspectral, dcdalpha(i)&
   &                                    , dcdalphad(i), coef0(i), coef0d(i)&
   &                                   )
   END DO
   resbasecoefd = 0.0_8
   ! now subtract off estimated cl,cmz and use remainder to compute 
   ! clalphadot and cmzalphadot.
   DO i=1,8
   DO sps=1,ntimeintervalsspectral
   resbasecoefd(sps, i) = basecoefd(sps, i) - intervalalpha(sps)*&
   &         dcdalphad(i) - coef0d(i)
   resbasecoef(sps, i) = basecoef(sps, i) - (dcdalpha(i)*&
   &         intervalalpha(sps)+coef0(i))
   END DO
   END DO
   dcdalphadotd = 0.0_8
   !now compute dCi/dalphadot
   DO i=1,8
   CALL COMPUTELEASTSQUARESREGRESSION_D(resbasecoef(:, i), &
   &                                    resbasecoefd(:, i), &
   &                                    intervalalphadot, &
   &                                    ntimeintervalsspectral, dcdalphadot&
   &                                    (i), dcdalphadotd(i), coef0dot(i), &
   &                                    coef0dotd(i))
   END DO
   arg1d = ((gammainfd*pinfdim+gammainf*pinfdimd)*rhoinfdim-gammainf*&
   &     pinfdim*rhoinfdimd)/rhoinfdim**2
   arg1 = gammainf*pinfdim/rhoinfdim
   IF (arg1 .EQ. 0.0_8) THEN
   ad = 0.0_8
   ELSE
   ad = arg1d/(2.0*SQRT(arg1))
   END IF
   a = SQRT(arg1)
   dcdalphadotd = (2*((dcdalphadotd*machgrid+dcdalphadot*machgridd)*a+&
   &     dcdalphadot*machgrid*ad)*lengthref-dcdalphadot*2*machgrid*a*&
   &     lengthrefd)/lengthref**2
   dcdalphadot = dcdalphadot*2*(machgrid*a)/lengthref
   ELSE
   CALL TERMINATE('computeTSDerivatives', &
   &               'Not a valid stability motion')
   dcdalphadotd = 0.0_8
   coef0d = 0.0_8
   dcdalphad = 0.0_8
   END IF
   END SUBROUTINE COMPUTETSDERIVATIVES_D
