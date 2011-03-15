!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade - Version 2.2 (r1239) - Wed 28 Jun 2006 04:59:55 PM CEST
!  
!  Differentiation of computetsstabilityderivadj in reverse (adjoint) mode:
!   gradient, with respect to input variables: dcdalphadot basecoef
!                coef0 dcdq dcdqdot dcdalpha
!   of linear combination of output variables: dcdalphadot coef0
!                dcdq dcdqdot dcdalpha
!
!     ******************************************************************
!     *                                                                *
!     * File:          computeTSStabilityDerivAdj.f90                  *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 11-25-2009                                      *
!     * Last modified: 11-26-2009                                      *
!     *                                                                *
!     ******************************************************************
!
SUBROUTINE COMPUTETSSTABILITYDERIVADJ_B(basecoef, basecoefb, coef0, &
&  coef0b, dcdalpha, dcdalphab, dcdalphadot, dcdalphadotb, dcdq, dcdqb, &
&  dcdqdot, dcdqdotb)
  USE communication
  USE inputmotion
  USE inputphysics
  USE inputtimespectral
  USE inputtsstabderiv
  USE monitor
  USE section
  IMPLICIT NONE
  REAL(KIND=REALTYPE) :: basecoef(ntimeintervalsspectral, 8), basecoefb(&
&  ntimeintervalsspectral, 8)
  REAL(KIND=REALTYPE) :: coef0(8), coef0b(8)
  REAL(KIND=REALTYPE) :: dcdalpha(8), dcdalphab(8), dcdalphadot(8), &
&  dcdalphadotb(8)
  REAL(KIND=REALTYPE) :: dcdq(8), dcdqb(8), dcdqdot(8), dcdqdotb(8)
  REAL(KIND=REALTYPE) :: a
  INTEGER :: branch
  REAL(KIND=REALTYPE) :: coef0dot(8), coef0dotb(8)
  REAL(KIND=REALTYPE) :: dcdbeta(8), dcdbetadot(8), dcdmach(8), dcdmachb&
&  (8), dcdmachdot(8)
  REAL(KIND=REALTYPE) :: dcdp(8), dcdpb(8), dcdpdot(8), dcdr(8), dcdrb(8&
&  ), dcdrdot(8)
  REAL(KIND=REALTYPE) :: DERIVATIVERIGIDROTANGLE, res, res0, &
&  SECONDDERIVATIVERIGIDROTANGLE
  REAL(KIND=REALTYPE) :: dphix(ntimeintervalsspectral), dphiy(&
&  ntimeintervalsspectral), dphiz(ntimeintervalsspectral)
  REAL(KIND=REALTYPE) :: dphixdot(ntimeintervalsspectral), dphiydot(&
&  ntimeintervalsspectral), dphizdot(ntimeintervalsspectral)
  REAL :: gammainf
  REAL(KIND=REALTYPE) :: intervalalpha(ntimeintervalsspectral), &
&  intervalalphadot(ntimeintervalsspectral)
  REAL(KIND=REALTYPE) :: intervalmach(ntimeintervalsspectral), &
&  intervalmachdot(ntimeintervalsspectral)
  INTEGER(KIND=INTTYPE) :: i, nn, sps
  REAL :: pinfdim
  REAL(KIND=REALTYPE) :: resbasecoef(ntimeintervalsspectral, 8), &
&  resbasecoefb(ntimeintervalsspectral, 8)
  REAL :: rhoinfdim
  REAL(KIND=REALTYPE) :: t(nsections)
  REAL :: timeref
  REAL(KIND=REALTYPE) :: TSALPHA, TSALPHADOT
  REAL(KIND=REALTYPE) :: res1, result1, TSMACH, TSMACHDOT
  INTRINSIC SQRT
  EXTERNAL TSALPHA, TSMACHDOT, DERIVATIVERIGIDROTANGLE, TSMACH, &
&      TSALPHADOT, SECONDDERIVATIVERIGIDROTANGLE, TERMINATE
!Given
  IF (tspmode) THEN
!P is roll
    DO sps=1,ntimeintervalsspectral
!compute the time of this interval
      t = timeunsteadyrestart
      IF (equationmode .EQ. timespectral) THEN
        DO nn=1,nsections
!t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
!     /         real(nTimeIntervalsSpectral,realType)
          t(nn) = t(nn) + (sps-1)*sections(nn)%timeperiod/(&
&            ntimeintervalsspectral*1.0)
        END DO
        CALL PUSHINTEGER4(1)
      ELSE
        CALL PUSHINTEGER4(0)
      END IF
! Compute the time derivative of the rotation angles around the
! x-axis. i.e. compute p
      dphix(sps) = DERIVATIVERIGIDROTANGLE(degreepolxrot, coefpolxrot, &
&        degreefourxrot, omegafourxrot, coscoeffourxrot, sincoeffourxrot&
&        , t)
!if(myID==0)print *,'dphix',dphix
!add in pdot computation here!
      res = SECONDDERIVATIVERIGIDROTANGLE(degreepolxrot, coefpolxrot, &
&        degreefourxrot, omegafourxrot, coscoeffourxrot, sincoeffourxrot&
&        , t)
    END DO
!now compute dCl/dp
    DO i=1,8
      CALL COMPUTELEASTSQUARESREGRESSION(basecoef(:, i), dphix, &
&                                   ntimeintervalsspectral, dcdp(i), &
&                                   coef0(i))
    END DO
    basecoefb(1:ntimeintervalsspectral, 1:8) = 0.0
    DO i=8,1,-1
      dcdpb(:) = 0.0
      CALL COMPUTELEASTSQUARESREGRESSION_B(basecoef(:, i), basecoefb(:, &
&                                     i), dphix, ntimeintervalsspectral&
&                                     , dcdp(i), dcdpb(i), coef0(i), &
&                                     coef0b(i))
      coef0b(i) = 0.0
    END DO
    DO sps=ntimeintervalsspectral,1,-1
      CALL POPINTEGER4(branch)
      IF (.NOT.branch .LT. 1) nn = 0
    END DO
  ELSE IF (tsqmode) THEN
!!$         if(myID==0)then
!!$            print *,'CL estimates:'
!!$            print *,'Clpdot = : ',dcldpdot,' cl0dot = : ',cl0dot
!!$            print *,'Cd estimates:'
!!$            print *,'Cdpdot = : ',dcddpdot,' cd0dot = : ',cd0dot
!!$            print *,'CMz estimates:'
!!$            print *,'CMzpdot = : ',dcmzdpdot,' cmz0dot = : ',cmz0dot
!!$         end if
!q is pitch
    DO sps=1,ntimeintervalsspectral
!compute the time of this interval
      t = timeunsteadyrestart
      IF (equationmode .EQ. timespectral) THEN
        DO nn=1,nsections
!t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
!     /         real(nTimeIntervalsSpectral,realType)
          t(nn) = t(nn) + (sps-1)*sections(nn)%timeperiod/(&
&            ntimeintervalsspectral*1.0)
        END DO
        CALL PUSHINTEGER4(1)
      ELSE
        CALL PUSHINTEGER4(0)
      END IF
! Compute the time derivative of the rotation angles around the
! z-axis. i.e. compute q
      dphiz(sps) = DERIVATIVERIGIDROTANGLE(degreepolzrot, coefpolzrot, &
&        degreefourzrot, omegafourzrot, coscoeffourzrot, sincoeffourzrot&
&        , t)
!if(myID==0)print *,'dphiz',dphiz
! add in q_dot computation
      dphizdot(sps) = SECONDDERIVATIVERIGIDROTANGLE(degreepolzrot, &
&        coefpolzrot, degreefourzrot, omegafourzrot, coscoeffourzrot, &
&        sincoeffourzrot, t)
    END DO
!if(myID==0)print *,'dphizdot',dphizdot
!!$         if(myID==0)print *,'dphiz',dphiz
!!$         if(myID==0)print *,'dphizdot',dphizdot
!now compute dCl/dq
    DO i=1,8
      CALL COMPUTELEASTSQUARESREGRESSION(basecoef(:, i), dphiz, &
&                                   ntimeintervalsspectral, dcdq(i), &
&                                   coef0(i))
    END DO
!ResCL(i) = Cl(i)-(dcldp*dphix(i)+cl0)
!ResCD(i) = Cd(i)-(dcddp*dphix(i)+cd0)
!ResCMz(i) = Cmz(i)-(dcmzdp*dphix(i)+cmz0)
!now normalize the results...
    a = SQRT(gammainf*pinfdim/rhoinfdim)
!now compute dCl/dpdot
    DO i=1,8
      CALL COMPUTELEASTSQUARESREGRESSION(resbasecoef(:, i), dphizdot, &
&                                   ntimeintervalsspectral, dcdqdot(i), &
&                                   coef0dot(i))
    END DO
    resbasecoefb(1:ntimeintervalsspectral, 1:8) = 0.0
    DO i=8,1,-1
      coef0dotb(:) = 0.0
      CALL COMPUTELEASTSQUARESREGRESSION_B(resbasecoef(:, i), &
&                                     resbasecoefb(:, i), dphizdot, &
&                                     ntimeintervalsspectral, dcdqdot(i)&
&                                     , dcdqdotb(i), coef0dot(i), &
&                                     coef0dotb(i))
      dcdqdotb(i) = 0.0
    END DO
    dcdqb = timeref*2*machgrid*a*dcdqb/lengthref
    basecoefb(1:ntimeintervalsspectral, 1:8) = 0.0
    DO i=8,1,-1
      DO sps=ntimeintervalsspectral,1,-1
        basecoefb(sps, i) = basecoefb(sps, i) + resbasecoefb(sps, i)
        dcdqb(i) = dcdqb(i) - dphiz(sps)*resbasecoefb(sps, i)
        coef0b(i) = coef0b(i) - resbasecoefb(sps, i)
        resbasecoefb(sps, i) = 0.0
      END DO
    END DO
    DO i=8,1,-1
      CALL COMPUTELEASTSQUARESREGRESSION_B(basecoef(:, i), basecoefb(:, &
&                                     i), dphiz, ntimeintervalsspectral&
&                                     , dcdq(i), dcdqb(i), coef0(i), &
&                                     coef0b(i))
      coef0b(i) = 0.0
      dcdqb(i) = 0.0
    END DO
    DO sps=ntimeintervalsspectral,1,-1
      CALL POPINTEGER4(branch)
      IF (.NOT.branch .LT. 1) nn = 0
    END DO
  ELSE IF (tsrmode) THEN
!!$         if(myID==0)then
!!$            print *,'CL estimates:'
!!$            print *,'Clqdot = : ',dcldqdot,' cl0dot = : ',cl0dot
!!$            print *,'Cd estimates:'
!!$            print *,'Cdqdot = : ',dcddqdot,' cd0dot = : ',cd0dot
!!$            print *,'CMz estimates:'
!!$            print *,'CMzqdot = : ',dcmzdqdot,' cmz0dot = : ',cmz0dot
!!$         end if
!r is yaw
    DO sps=1,ntimeintervalsspectral
!compute the time of this interval
      t = timeunsteadyrestart
      IF (equationmode .EQ. timespectral) THEN
        DO nn=1,nsections
!t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
!     /         real(nTimeIntervalsSpectral,realType)
          t(nn) = t(nn) + (sps-1)*sections(nn)%timeperiod/(&
&            ntimeintervalsspectral*1.0)
        END DO
        CALL PUSHINTEGER4(1)
      ELSE
        CALL PUSHINTEGER4(0)
      END IF
! Compute the time derivative of the rotation angles around the
! y-axis and z-axis. i.e. compute r
      dphiy(sps) = DERIVATIVERIGIDROTANGLE(degreepolyrot, coefpolyrot, &
&        degreefouryrot, omegafouryrot, coscoeffouryrot, sincoeffouryrot&
&        , t)
!if(myID==0)print *,'dphiy',dphiy
! add in r_dot computation
      res0 = SECONDDERIVATIVERIGIDROTANGLE(degreepolyrot, coefpolyrot, &
&        degreefouryrot, omegafouryrot, coscoeffouryrot, sincoeffouryrot&
&        , t)
    END DO
    DO i=1,8
      CALL COMPUTELEASTSQUARESREGRESSION(basecoef(:, i), dphiy, &
&                                   ntimeintervalsspectral, dcdr(i), &
&                                   coef0(i))
    END DO
    basecoefb(1:ntimeintervalsspectral, 1:8) = 0.0
    DO i=8,1,-1
      dcdrb(:) = 0.0
      CALL COMPUTELEASTSQUARESREGRESSION_B(basecoef(:, i), basecoefb(:, &
&                                     i), dphiy, ntimeintervalsspectral&
&                                     , dcdr(i), dcdrb(i), coef0(i), &
&                                     coef0b(i))
      coef0b(i) = 0.0
    END DO
    DO sps=ntimeintervalsspectral,1,-1
      CALL POPINTEGER4(branch)
      IF (.NOT.branch .LT. 1) nn = 0
    END DO
  ELSE IF (tsalphamode) THEN
!compute the alphas and alphadots
    DO sps=1,ntimeintervalsspectral
!compute the time of this interval
      t = timeunsteadyrestart
      IF (equationmode .EQ. timespectral) THEN
        DO nn=1,nsections
!t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
!     /         real(nTimeIntervalsSpectral,realType)
          t(nn) = t(nn) + (sps-1)*sections(nn)%timeperiod/(&
&            ntimeintervalsspectral*1.0)
        END DO
        CALL PUSHINTEGER4(1)
      ELSE
        CALL PUSHINTEGER4(0)
      END IF
!print *,'t',t
      intervalalpha(sps) = TSALPHA(degreepolalpha, coefpolalpha, &
&        degreefouralpha, omegafouralpha, coscoeffouralpha, &
&        sincoeffouralpha, t)
      intervalalphadot(sps) = TSALPHADOT(degreepolalpha, coefpolalpha, &
&        degreefouralpha, omegafouralpha, coscoeffouralpha, &
&        sincoeffouralpha, t)
    END DO
!print *,'Alpha', intervalAlpha
!now compute dCl/dalpha
    DO i=1,8
      CALL COMPUTELEASTSQUARESREGRESSION(basecoef(:, i), intervalalpha, &
&                                   ntimeintervalsspectral, dcdalpha(i)&
&                                   , coef0(i))
    END DO
!now compute dCi/dalphadot
    DO i=1,8
      CALL COMPUTELEASTSQUARESREGRESSION(resbasecoef(:, i), &
&                                   intervalalphadot, &
&                                   ntimeintervalsspectral, dcdalphadot(&
&                                   i), coef0dot(i))
    END DO
    resbasecoefb(1:ntimeintervalsspectral, 1:8) = 0.0
    DO i=8,1,-1
      coef0dotb(:) = 0.0
      CALL COMPUTELEASTSQUARESREGRESSION_B(resbasecoef(:, i), &
&                                     resbasecoefb(:, i), &
&                                     intervalalphadot, &
&                                     ntimeintervalsspectral, &
&                                     dcdalphadot(i), dcdalphadotb(i), &
&                                     coef0dot(i), coef0dotb(i))
      dcdalphadotb(i) = 0.0
    END DO
    basecoefb(1:ntimeintervalsspectral, 1:8) = 0.0
    DO i=8,1,-1
      DO sps=ntimeintervalsspectral,1,-1
        basecoefb(sps, i) = basecoefb(sps, i) + resbasecoefb(sps, i)
        dcdalphab(i) = dcdalphab(i) - intervalalpha(sps)*resbasecoefb(&
&          sps, i)
        coef0b(i) = coef0b(i) - resbasecoefb(sps, i)
        resbasecoefb(sps, i) = 0.0
      END DO
    END DO
    DO i=8,1,-1
      CALL COMPUTELEASTSQUARESREGRESSION_B(basecoef(:, i), basecoefb(:, &
&                                     i), intervalalpha, &
&                                     ntimeintervalsspectral, dcdalpha(i&
&                                     ), dcdalphab(i), coef0(i), coef0b(&
&                                     i))
      coef0b(i) = 0.0
      dcdalphab(i) = 0.0
    END DO
    DO sps=ntimeintervalsspectral,1,-1
      CALL POPINTEGER4(branch)
      IF (.NOT.branch .LT. 1) nn = 0
    END DO
  ELSE IF (tsbetamode) THEN
    basecoefb(1:ntimeintervalsspectral, 1:8) = 0.0
  ELSE IF (tsmachmode) THEN
!compute the alphas and alphadots
    DO sps=1,ntimeintervalsspectral
!compute the time of this interval
      t = timeunsteadyrestart
      IF (equationmode .EQ. timespectral) THEN
        DO nn=1,nsections
!t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
!     /         real(nTimeIntervalsSpectral,realType)
          t(nn) = t(nn) + (sps-1)*sections(nn)%timeperiod/(&
&            ntimeintervalsspectral*1.0)
        END DO
        CALL PUSHINTEGER4(1)
      ELSE
        CALL PUSHINTEGER4(0)
      END IF
!if(myID==0)print *,'t',t
      result1 = TSMACH(degreepolmach, coefpolmach, degreefourmach, &
&        omegafourmach, coscoeffourmach, sincoeffourmach, t)
      intervalmach(sps) = machgrid + result1
      res1 = TSMACHDOT(degreepolmach, coefpolmach, degreefourmach, &
&        omegafourmach, coscoeffourmach, sincoeffourmach, t)
    END DO
!if(myID==0) print *,'Mach', intervalMach,machgrid
!now compute dCl/dalpha
    DO i=1,8
      CALL COMPUTELEASTSQUARESREGRESSION(basecoef(:, i), intervalmach, &
&                                   ntimeintervalsspectral, dcdmach(i), &
&                                   coef0(i))
    END DO
    basecoefb(1:ntimeintervalsspectral, 1:8) = 0.0
    DO i=8,1,-1
      dcdmachb(:) = 0.0
      CALL COMPUTELEASTSQUARESREGRESSION_B(basecoef(:, i), basecoefb(:, &
&                                     i), intervalmach, &
&                                     ntimeintervalsspectral, dcdmach(i)&
&                                     , dcdmachb(i), coef0(i), coef0b(i)&
&                                    )
      coef0b(i) = 0.0
    END DO
    DO sps=ntimeintervalsspectral,1,-1
      CALL POPINTEGER4(branch)
      IF (.NOT.branch .LT. 1) nn = 0
    END DO
  ELSE
    basecoefb(1:ntimeintervalsspectral, 1:8) = 0.0
  END IF
  dcdalphadotb(1:8) = 0.0
  coef0b(1:8) = 0.0
  dcdqb(1:8) = 0.0
  dcdqdotb(1:8) = 0.0
  dcdalphab(1:8) = 0.0
END SUBROUTINE COMPUTETSSTABILITYDERIVADJ_B
