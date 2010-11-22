!
!     ******************************************************************
!     *                                                                *
!     * File:          computeTSdetrivatives.f90                       *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 09-15-2009                                      *
!     * Last modified: 11-26-2009                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine computeTSDerivatives(cl0,cd0,cmz0,dcldalpha,dcddalpha,&
           dcmzdalpha,dcldalphadot,dcddalphadot,dcmzdalphadot,dcldq,&
           dcddq,dcmzdq,dcldqdot,dcddqdot,dcmzdqdot)
!
!     ******************************************************************
!     *                                                                *
!     * Computes an estimate for the stability derivatives from the    *
!     * coefficient values of the Time Spectral Solution.              *
!     *                                                                *
!     ******************************************************************
!
      use precision
      use inputMotion
      use inputPhysics
      use inputTimeSpectral
      use inputTSStabDeriv
      use flowvarrefstate !Timeref
      use monitor
      use section
      use communication !myid
      implicit none

!
!     Subroutine arguments.
!

!
!     Local variables.
!
         real(kind=realType),dimension(nTimeIntervalsSpectral)  :: CL, CD,&
              CFx, CFy,CFz, CMx, CMy, CMz

         real(kind=realType),dimension(nTimeIntervalsSpectral)  :: dPhix,&
              dPhiy,dPhiz
         
         real(kind=realType),dimension(nTimeIntervalsSpectral)  :: intervalAlpha,intervalAlphadot
         real(kind=realType),dimension(nTimeIntervalsSpectral)  :: intervalMach,intervalMachdot

         real(kind=realType),dimension(nTimeIntervalsSpectral)  :: resCL,&
              resCD, resCFx, resCFy,resCFz, resCMx, resCMy, resCMz

         real(kind=realType),dimension(nTimeIntervalsSpectral)  :: dPhixdot,&
              dPhiydot,dPhizdot

         real(kind=realType)::dcldp,dcldpdot,dcddp,dcddpdot,dcmzdp,dcmzdpdot         
         real(kind=realType)::dcldq,dcldqdot,dcddq,dcddqdot,dcmzdq,dcmzdqdot
         real(kind=realType)::dcldr,dcldrdot,dcddr,dcddrdot,dcmzdr,dcmzdrdot
         real(kind=realType)::dcldalpha,dcldalphadot,dcddalpha,dcddalphadot,dcmzdalpha,dcmzdalphadot
         real(kind=realType)::dcldMach,dcldMachdot,dcddMach,dcddMachdot,dcmzdMach,dcmzdMachdot
         real(kind=realType)::cl0,cl0dot,cd0,cd0dot,cmz0,cmz0dot

         real(kind=realType), dimension(nSections) :: t
         integer(kind=inttype)::level,sps,i,nn

         !function definitions
         real(kind=realType)::derivativeRigidRotAngle,&
              secondDerivativeRigidRotAngle
         !real(kind=realType)::t
         real(kind=realType):: TSAlpha,TSAlphadot
         real(kind=realType):: TSMach,TSMachdot

         !speed of sound: for normalization of q derivatives
         real(kind=realType)::a

!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!

      if(myID==0) print *,'in compute TS deriv...'
         !Compute and store the aero coef. Values for each TS level
         do sps =1,nTimeIntervalsSpectral
            !if(myID==0) print *,'sps',sps
            level = 1
            call computeAeroCoef(CL(sps),CD(sps),CFx(sps),CFy(sps),CFz(sps),&
                 CMx(sps),CMy(sps),CMz(sps),level,sps)
            !if(myID==0) print *,'CL',cl
           
         end do
         
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
         if (myid==0) then
            print *,'CL',CL
            print *,'dphix',dphix
         endif
         call computeLeastSquaresRegression(cl,dphix,nTimeIntervalsSpectral,dcldp,cl0)

         !now compute dCd/dp
         
         call computeLeastSquaresRegression(cd,dphix,nTimeIntervalsSpectral,dcddp,cd0)
         
         !now compute dCmz/dp
         
         call computeLeastSquaresRegression(cmz,dphix,nTimeIntervalsSpectral,dcmzdp,cmz0)
         if(myID==0)then
            print *,'CL estimates:'
            print *,'Clp = : ',dcldp,' cl0 = : ',cl0
            print *,'CD estimates:'
            print *,'CDp = : ',dcddp,' cd0 = : ',cd0
            print *,'CMz estimates:'
            print *,'CMzp = : ',dcmzdp,' cmz0 = : ',cmz0
         endif
         
         ! now subtract off estimated cl,cmz and use remainder to compute 
         ! clpdot and cmzpdot.
         
         do i = 1,nTimeIntervalsSpectral
            ResCL(i) = Cl(i)-(dcldp*dphix(i)+cl0)
            ResCD(i) = Cd(i)-(dcddp*dphix(i)+cd0)
            ResCMz(i) = Cmz(i)-(dcmzdp*dphix(i)+cmz0)
         enddo
         if(myID==0) print *,'ResCL',resCL
         if(myID==0) print *,'ResCd',resCD
         if(myID==0) print *,'resCMZ',resCMZ
         
         !now compute dCl/dpdot
         
         call computeLeastSquaresRegression(Rescl,dphixdot,nTimeIntervalsSpectral,dcldpdot,cl0dot)

         !now compute dCd/dpdot
         
         call computeLeastSquaresRegression(Rescd,dphixdot,nTimeIntervalsSpectral,dcddpdot,cl0dot)
         
         !now compute dCmz/dpdot
         
         call computeLeastSquaresRegression(Rescmz,dphixdot,nTimeIntervalsSpectral,dcmzdpdot,cmz0dot)
         
         if(myID==0)then
            print *,'CL estimates:'
            print *,'Clpdot = : ',dcldpdot,' cl0dot = : ',cl0dot
            print *,'Cd estimates:'
            print *,'Cdpdot = : ',dcddpdot,' cd0dot = : ',cd0dot
            print *,'CMz estimates:'
            print *,'CMzpdot = : ',dcmzdpdot,' cmz0dot = : ',cmz0dot
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
            if (myid==0)  print *,'t',t
            
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
         !if(myID==0)print *,'dphiz',dphiz
         !if(myID==0)print *,'dphizdot',dphizdot
         
         if (myid==0) then
            print *,'CL',CL
            print *,'dphiz',dphiz
            print *,'dphizdot',dphizdot
         endif
         !now compute dCl/dq
         
         call computeLeastSquaresRegression(cl,dphiz,nTimeIntervalsSpectral,dcldq,cl0)

         !now compute dCd/dq
         
         call computeLeastSquaresRegression(cd,dphiz,nTimeIntervalsSpectral,dcddq,cd0)
         
         !now compute dCmz/dq
         
         call computeLeastSquaresRegression(cmz,dphiz,nTimeIntervalsSpectral,dcmzdq,cmz0)
         if(myID==0)then
            a  = sqrt(gammaInf*pInfDim/rhoInfDim)
            print *,'normalization',timeRef,(machGrid*a),lengthRef
            print *,'CL estimates:'
            print *,'Clq = : ',dcldq*timeRef*2*(machGrid*a)/lengthRef,' cl0 = : ',cl0
            print *,'CD estimates:'
            print *,'Cdq = : ',dcddq*timeRef*2*(machGrid*a)/lengthRef,' cd0 = : ',cd0
            print *,'CMz estimates:'
            print *,'CMzq = : ',dcmzdq*timeRef*2*(machGrid*a)/lengthRef,' cmz0 = : ',cmz0
         endif
         
         ! now subtract off estimated cl,cmz and use remainder to compute 
         ! clqdot and cmzqdot.
         
         do i = 1,nTimeIntervalsSpectral
            ResCL(i) = Cl(i)-(dcldq*dphiz(i)+cl0)
            ResCd(i) = Cd(i)-(dcddq*dphiz(i)+cd0)
            ResCMz(i) = Cmz(i)-(dcmzdq*dphiz(i)+cmz0)
         enddo
         if(myID==0) print *,'ResCL',resCL
         if(myID==0) print *,'ResCd',resCd
         if(myID==0) print *,'resCMZ',resCMZ
         
         !now compute dCl/dqdot
            
         call computeLeastSquaresRegression(Rescl,dphizdot,nTimeIntervalsSpectral,dcldqdot,cl0dot)

      
         !now compute dCd/dqdot
            
         call computeLeastSquaresRegression(Rescd,dphizdot,nTimeIntervalsSpectral,dcddqdot,cd0dot)
         
         !now compute dCmz/dqdot
         
         call computeLeastSquaresRegression(Rescmz,dphizdot,nTimeIntervalsSpectral,dcmzdqdot,cmz0dot)
         
         if(myID==0)then
            print *,'CL estimates:'
            print *,'Clqdot = : ',dcldqdot,' cl0dot = : ',cl0dot
            print *,'Cd estimates:'
            print *,'Cdqdot = : ',dcddqdot,' cd0dot = : ',cd0dot
            print *,'CMz estimates:'
            print *,'CMzqdot = : ',dcmzdqdot,' cmz0dot = : ',cmz0dot
         end if
         
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
         
         if (myid==0) then
            print *,'CL',CL
            print *,'dphiy',dphiy
         endif
         !now compute dCl/dr
         
         call computeLeastSquaresRegression(cl,dphiy,nTimeIntervalsSpectral,dcldr,cl0)
         
         !now compute dCD/dr
         
         call computeLeastSquaresRegression(cd,dphiy,nTimeIntervalsSpectral,dcddr,cd0)
         
         !now compute dCmz/dr
         
         call computeLeastSquaresRegression(cmz,dphiy,nTimeIntervalsSpectral,dcmzdr,cmz0)
         if(myID==0)then
            print *,'CL estimates:'
            print *,'Clr = : ',dcldr,' cl0 = : ',cl0
            print *,'Cd estimates:'
            print *,'Cdr = : ',dcddr,' cd0 = : ',cd0
            print *,'CMz estimates:'
            print *,'CMzr = : ',dcmzdr,' cmz0 = : ',cmz0
         endif
         
         ! now subtract off estimated cl,cmz and use remainder to compute 
         ! clrdot and cmzrdot.
         
         do i = 1,nTimeIntervalsSpectral
            ResCL(i) = Cl(i)-(dcldr*dphiy(i)+cl0)
            ResCd(i) = Cd(i)-(dcddr*dphiy(i)+cd0)
            ResCMz(i) = Cmz(i)-(dcmzdr*dphiy(i)+cmz0)
         enddo
         if(myID==0) print *,'ResCL',resCL
         if(myID==0) print *,'ResCd',resCd
         if(myID==0) print *,'resCMZ',resCMZ
         
         
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
         !print *,'Alpha', intervalAlpha
         !now compute dCl/dalpha
         if (myid==0) then
            print *,'CL',CL
            print *,'intervalAlpha',intervalAlpha
         endif
         call computeLeastSquaresRegression(cl,intervalAlpha,nTimeIntervalsSpectral,dcldalpha,cl0)
         
         !now compute dCd/dalpha
         
         call computeLeastSquaresRegression(cd,intervalAlpha,nTimeIntervalsSpectral,dcddalpha,cd0)
         
         !now compute dCmz/dq
         
         call computeLeastSquaresRegression(cmz,intervalAlpha,nTimeIntervalsSpectral,dcmzdalpha,cmz0)
         if(myID==0)then
            print *,'CL estimates:'
            print *,'ClAlpha = : ',dcldalpha,' cl0 = : ',cl0
            print *,'Cd estimates:'
            print *,'CdAlpha = : ',dcddalpha,' cd0 = : ',cd0
            print *,'CMz estimates:'
            print *,'CMzAlpha = : ',dcmzdalpha,' cmz0 = : ',cmz0
         endif
         
         ! now subtract off estimated cl,cmz and use remainder to compute 
         ! clqdot and cmzqdot.
         
         do i = 1,nTimeIntervalsSpectral
            ResCL(i) = Cl(i)-(dcldalpha*IntervalAlpha(i)+cl0)
            ResCd(i) = Cd(i)-(dcddalpha*IntervalAlpha(i)+cd0)
            ResCMz(i) = Cmz(i)-(dcmzdalpha*IntervalAlpha(i)+cmz0)
         enddo
         if(myID==0) print *,'ResCL',resCL
         if(myID==0) print *,'ResCd',resCd
         if(myID==0) print *,'resCMZ',resCMZ

         if(myID==0) print *,'intervalAlphaDot',intervalAlphadot
         !now compute dCl/dqdot
         
         call computeLeastSquaresRegression(Rescl,intervalAlphadot,nTimeIntervalsSpectral,dcldalphadot,cl0dot)

         !now compute dCd/dqdot
         
         call computeLeastSquaresRegression(Rescd,intervalAlphadot,nTimeIntervalsSpectral,dcddalphadot,cd0dot)
         
         !now compute dCmz/dqdot
         
         call computeLeastSquaresRegression(Rescmz,intervalAlphadot,nTimeIntervalsSpectral,dcmzdalphadot,cmz0dot)
         
         if(myID==0)then
            print *,'CL estimates:'
            print *,'ClAlphadot = : ',dcldalphadot,' cl0dot = : ',cl0dot
            print *,'Cd estimates:'
            print *,'CdAlphadot = : ',dcddalphadot,' cl0dot = : ',cd0dot
            print *,'CMz estimates:'
            print *,'CMzAlphadot = : ',dcmzdalphadot,' cmz0dot = : ',cmz0dot
         end if
         
      elseif(TSMachMode)then
         
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
            !if(myID==0)print *,'t',t
            intervalMach(sps) = machgrid+TSMach(degreePolMach,   coefPolMach,       &
                 degreeFourMach,  omegaFourMach,     &
                 cosCoefFourMach, sinCoefFourMach, t)
            
            intervalMachdot(sps) = TSMachdot(degreePolMach,   coefPolMach,       &
                 degreeFourMach,  omegaFourMach,     &
                 cosCoefFourMach, sinCoefFourMach, t)
         end do
         !if(myID==0) print *,'Mach', intervalMach,machgrid
         !now compute dCl/dalpha
         
         call computeLeastSquaresRegression(cl,intervalMach,nTimeIntervalsSpectral,dcldMach,cl0)

         !now compute dCd/dalpha
         
         call computeLeastSquaresRegression(cd,intervalMach,nTimeIntervalsSpectral,dcddMach,cd0)
         
         !now compute dCmz/dq
         
         call computeLeastSquaresRegression(cmz,intervalMach,nTimeIntervalsSpectral,dcmzdMach,cmz0)
         if(myID==0)then
            print *,'CL estimates:'
            print *,'ClMach = : ',dcldMach,' cl0 = : ',cl0
            print *,'Cd estimates:'
            print *,'CdMach = : ',dcddMach,' cl0 = : ',cd0
            print *,'CMz estimates:'
            print *,'CMzMach = : ',dcmzdMach,' cmz0 = : ',cmz0
         endif
         
         ! now subtract off estimated cl,cmz and use remainder to compute 
         ! clqdot and cmzqdot.
         
         do i = 1,nTimeIntervalsSpectral
            ResCL(i) = Cl(i)-(dcldalpha*IntervalMach(i)+cl0)
            ResCd(i) = Cd(i)-(dcddalpha*IntervalMach(i)+cd0)
            ResCMz(i) = Cmz(i)-(dcmzdalpha*IntervalMach(i)+cmz0)
         enddo
         if(myID==0) print *,'ResCL',resCL
         if(myID==0) print *,'ResCd',resCd
         if(myID==0) print *,'resCMZ',resCMZ
         !now compute dCl/dqdot
         
         call computeLeastSquaresRegression(Rescl,intervalMachdot,nTimeIntervalsSpectral,dcldMachdot,cl0dot)

         !now compute dCd/dqdot
         
         call computeLeastSquaresRegression(Rescd,intervalMachdot,nTimeIntervalsSpectral,dcddMachdot,cd0dot)
         
         !now compute dCmz/dqdot
         
         call computeLeastSquaresRegression(Rescmz,intervalAlphadot,nTimeIntervalsSpectral,dcmzdalphadot,cmz0dot)
         
         if(myID==0)then
            print *,'CL estimates:'
            print *,'ClMachdot = : ',dcldMachdot,' cl0dot = : ',cl0dot
            print *,'Cd estimates:'
            print *,'CdMachdot = : ',dcddMachdot,' cd0dot = : ',cd0dot
            print *,'CMz estimates:'
            print *,'CMzMachdot = : ',dcmzdMachdot,' cmz0dot = : ',cmz0dot
         end if
         
      else
         call terminate('computeTSDerivatives','Not a valid stability motion')
      endif
    end subroutine computeTSDerivatives
    
