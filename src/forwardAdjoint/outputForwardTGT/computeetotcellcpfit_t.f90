   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.7 (r4786) - 21 Feb 2013 15:53
   !
   !  Differentiation of computeetotcellcpfit in forward (tangent) mode (with options debugTangent i4 dr8 r8):
   !   variations   of useful results: *gamma *w
   !   with respect to varying inputs: *p *gamma *w tref rgas scale
   !   Plus diff mem management of: p:in gamma:in w:in
   !      ==================================================================
   SUBROUTINE COMPUTEETOTCELLCPFIT_T(i, j, k, scale, scaled, correctfork)
   USE FLOWVARREFSTATE
   USE BLOCKPOINTERS_D
   USE CPCURVEFITS
   USE DIFFSIZES
   !  Hint: ISIZE4OFDrfw should be the size of dimension 4 of array *w
   !  Hint: ISIZE3OFDrfw should be the size of dimension 3 of array *w
   !  Hint: ISIZE2OFDrfw should be the size of dimension 2 of array *w
   !  Hint: ISIZE1OFDrfw should be the size of dimension 1 of array *w
   !  Hint: ISIZE3OFDrfp should be the size of dimension 3 of array *p
   !  Hint: ISIZE2OFDrfp should be the size of dimension 2 of array *p
   !  Hint: ISIZE1OFDrfp should be the size of dimension 1 of array *p
   !  Hint: ISIZE3OFDrfgamma should be the size of dimension 3 of array *gamma
   !  Hint: ISIZE2OFDrfgamma should be the size of dimension 2 of array *gamma
   !  Hint: ISIZE1OFDrfgamma should be the size of dimension 1 of array *gamma
   IMPLICIT NONE
   !
   !      ******************************************************************
   !      *                                                                *
   !      * ComputeEtotCellCpfit will compute the total energy for the     *
   !      * given cell of the block given by the current pointers with the *
   !      * cp temperature curve fit model.                                *
   !      *                                                                *
   !      ******************************************************************
   !
   !
   !      Local parameter.
   !
   REAL(kind=realtype), PARAMETER :: twothird=two*third
   !
   !      Subroutine arguments.
   !
   INTEGER(kind=inttype), INTENT(IN) :: i, j, k
   REAL(kind=realtype), INTENT(IN) :: scale
   REAL(kind=realtype), INTENT(IN) :: scaled
   LOGICAL, INTENT(IN) :: correctfork
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: nn, mm, ii, start
   REAL(kind=realtype) :: pp, t, t2, cv, eint
   REAL(kind=realtype) :: ppd, td, t2d, cvd, eintd
   EXTERNAL DEBUG_TGT_HERE
   LOGICAL :: DEBUG_TGT_HERE
   INTRINSIC LOG
   IF (.TRUE. .AND. DEBUG_TGT_HERE('entry', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('p', p, pd, ISIZE1OFDrfp*ISIZE2OFDrfp*&
   &                        ISIZE3OFDrfp)
   CALL DEBUG_TGT_REAL8ARRAY('w', w, wd, ISIZE1OFDrfw*ISIZE2OFDrfw*&
   &                        ISIZE3OFDrfw*ISIZE4OFDrfw)
   CALL DEBUG_TGT_REAL8('tref', tref, trefd)
   CALL DEBUG_TGT_REAL8('rgas', rgas, rgasd)
   CALL DEBUG_TGT_REAL8('scale', scale, scaled)
   CALL DEBUG_TGT_DISPLAY('entry')
   END IF
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Compute the dimensional temperature.
   ppd = pd(i, j, k)
   pp = p(i, j, k)
   IF (correctfork) THEN
   ppd = ppd - twothird*(wd(i, j, k, irho)*w(i, j, k, itu1)+w(i, j, k, &
   &      irho)*wd(i, j, k, itu1))
   pp = pp - twothird*w(i, j, k, irho)*w(i, j, k, itu1)
   END IF
   td = ((trefd*pp+tref*ppd)*rgas*w(i, j, k, irho)-tref*pp*(rgasd*w(i, j&
   &    , k, irho)+rgas*wd(i, j, k, irho)))/(rgas*w(i, j, k, irho))**2
   t = tref*pp/(rgas*w(i, j, k, irho))
   ! Determine the case we are having here.
   IF (t .LE. cptrange(0)) THEN
   ! Temperature is less than the smallest temperature
   ! in the curve fits. Use extrapolation using
   ! constant cv.
   eintd = scaled*(cpeint(0)+cv0*(t-cptrange(0))) + scale*cv0*td
   eint = scale*(cpeint(0)+cv0*(t-cptrange(0)))
   gammad(i, j, k) = 0.0_8
   gamma(i, j, k) = (cv0+one)/cv0
   ELSE IF (t .GE. cptrange(cpnparts)) THEN
   ! Temperature is larger than the largest temperature
   ! in the curve fits. Use extrapolation using
   ! constant cv.
   eintd = scaled*(cpeint(cpnparts)+cvn*(t-cptrange(cpnparts))) + scale&
   &      *cvn*td
   eint = scale*(cpeint(cpnparts)+cvn*(t-cptrange(cpnparts)))
   gammad(i, j, k) = 0.0_8
   gamma(i, j, k) = (cvn+one)/cvn
   ELSE
   ! Temperature is in the curve fit range.
   ! First find the valid range.
   ii = cpnparts
   start = 1
   interval:DO 
   ! Next guess for the interval.
   nn = start + ii/2
   ! Determine the situation we are having here.
   IF (t .GT. cptrange(nn)) THEN
   ! Temperature is larger than the upper boundary of
   ! the current interval. Update the lower boundary.
   start = nn + 1
   ii = ii - 1
   ELSE
   IF (.TRUE. .AND. DEBUG_TGT_HERE('middle', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('gamma', gamma, gammad, &
   &                              ISIZE1OFDrfgamma*ISIZE2OFDrfgamma*&
   &                              ISIZE3OFDrfgamma)
   CALL DEBUG_TGT_REAL8ARRAY('w', w, wd, ISIZE1OFDrfw*&
   &                              ISIZE2OFDrfw*ISIZE3OFDrfw*ISIZE4OFDrfw)
   CALL DEBUG_TGT_REAL8('scale', scale, scaled)
   CALL DEBUG_TGT_REAL8('t', t, td)
   CALL DEBUG_TGT_DISPLAY('middle')
   END IF
   IF (t .GE. cptrange(nn-1)) GOTO 100
   END IF
   ! This is the correct range. Exit the do-loop.
   ! Modify ii for the next branch to search.
   ii = ii/2
   END DO interval
   ! Nn contains the correct curve fit interval.
   ! Integrate cv to compute eint.
   100 eintd = -td
   eint = cptempfit(nn)%eint0 - t
   cv = -one
   cvd = 0.0_8
   DO ii=1,cptempfit(nn)%nterm
   IF (t .GT. 0.0_8 .OR. (t .LT. 0.0_8 .AND. cptempfit(nn)%exponents(&
   &          ii) .EQ. INT(cptempfit(nn)%exponents(ii)))) THEN
   t2d = cptempfit(nn)%exponents(ii)*t**(cptempfit(nn)%exponents(ii&
   &          )-1)*td
   ELSE IF (t .EQ. 0.0_8 .AND. cptempfit(nn)%exponents(ii) .EQ. 1.0) &
   &      THEN
   t2d = td
   ELSE
   t2d = 0.0_8
   END IF
   t2 = t**cptempfit(nn)%exponents(ii)
   cvd = cvd + cptempfit(nn)%constants(ii)*t2d
   cv = cv + cptempfit(nn)%constants(ii)*t2
   IF (cptempfit(nn)%exponents(ii) .EQ. -1) THEN
   eintd = eintd + cptempfit(nn)%constants(ii)*td/t
   eint = eint + cptempfit(nn)%constants(ii)*LOG(t)
   ELSE
   mm = cptempfit(nn)%exponents(ii) + 1
   t2d = td*t2 + t*t2d
   t2 = t*t2
   eintd = eintd + cptempfit(nn)%constants(ii)*t2d/mm
   eint = eint + cptempfit(nn)%constants(ii)*t2/mm
   END IF
   END DO
   eintd = scaled*eint + scale*eintd
   eint = scale*eint
   gammad(i, j, k) = (cvd*cv-(cv+one)*cvd)/cv**2
   gamma(i, j, k) = (cv+one)/cv
   END IF
   ! Compute the total energy per unit volume.
   wd(i, j, k, irhoe) = wd(i, j, k, irho)*(eint+half*(w(i, j, k, ivx)**2+&
   &    w(i, j, k, ivy)**2+w(i, j, k, ivz)**2)) + w(i, j, k, irho)*(eintd+&
   &    half*(2*w(i, j, k, ivx)*wd(i, j, k, ivx)+2*w(i, j, k, ivy)*wd(i, j, &
   &    k, ivy)+2*w(i, j, k, ivz)*wd(i, j, k, ivz)))
   w(i, j, k, irhoe) = w(i, j, k, irho)*(eint+half*(w(i, j, k, ivx)**2+w(&
   &    i, j, k, ivy)**2+w(i, j, k, ivz)**2))
   IF (correctfork) THEN
   wd(i, j, k, irhoe) = wd(i, j, k, irhoe) + wd(i, j, k, irho)*w(i, j, &
   &      k, itu1) + w(i, j, k, irho)*wd(i, j, k, itu1)
   w(i, j, k, irhoe) = w(i, j, k, irhoe) + w(i, j, k, irho)*w(i, j, k, &
   &      itu1)
   END IF
   IF (.TRUE. .AND. DEBUG_TGT_HERE('exit', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('gamma', gamma, gammad, ISIZE1OFDrfgamma*&
   &                        ISIZE2OFDrfgamma*ISIZE3OFDrfgamma)
   CALL DEBUG_TGT_REAL8ARRAY('w', w, wd, ISIZE1OFDrfw*ISIZE2OFDrfw*&
   &                        ISIZE3OFDrfw*ISIZE4OFDrfw)
   CALL DEBUG_TGT_DISPLAY('exit')
   END IF
   END SUBROUTINE COMPUTEETOTCELLCPFIT_T
