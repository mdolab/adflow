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
      subroutine computeTSStabilityDerivAdj(cFxAdj,cFyAdj,cFzAdj,cMxAdj,&
                       cMyAdj,cMzAdj,CLAdj,CDAdj,&
                       cl0,cmz0,dcldalpha,dcmzdalpha)
 
!
!     ******************************************************************
!     *                                                                *
!     * Computes the stability derivatives based on the time spectral  *
!     * solution of a given mesh. Takes in the force coefficients at   *
!     * all time instantces and computes the agregate parameters       *
!     *                                                                *
!     ******************************************************************
!
!      use blockPointers        ! ie,je,ke
      use communication        ! procHalo(currentLevel)%nProcSend, myID
      use inputPhysics         ! equations
      use inputTimeSpectral    ! nTimeIntervalsSpectral!nTimeInstancesMax
      use inputTSStabDeriv
      use monitor              !timeunsteadyrestart,timeunsteady
      use section              !nsections, sections%
      use inputMotion          ! degreePol*
!      use bcTypes              ! EulerWall, ...
!      use flowvarrefstate      !nw

      implicit none

!
!     Subroutine arguments.
!
      
      real(kind=realType), dimension(nTimeIntervalsSpectral)::       &
                             ClAdj,CdAdj,CfxAdj,CfyAdj,CfzAdj,   &
                             CmxAdj,CmyAdj,CmzAdj
     
      real(kind=realType)::dcldp,dcldpdot,dcmzdp,dcmzdpdot         
      real(kind=realType)::dcldq,dcldqdot,dcmzdq,dcmzdqdot
      real(kind=realType)::dcldr,dcldrdot,dcmzdr,dcmzdrdot
      real(kind=realType)::dcldalpha,dcldalphadot,dcmzdalpha,dcmzdalphadot
      real(kind=realType)::dcldMach,dcldMachdot,dcmzdMach,dcmzdMachdot
      real(kind=realType)::cl0,cl0dot,cmz0,cmz0dot
      
      
!
!     Local variables.
!
     

      !Stability derivative variables
      real(kind=realType),dimension(nTimeIntervalsSpectral)  :: dPhix,&
           dPhiy,dPhiz
      
      real(kind=realType),dimension(nTimeIntervalsSpectral)  :: intervalAlpha,intervalAlphadot
      real(kind=realType),dimension(nTimeIntervalsSpectral)  :: intervalMach,intervalMachdot
      
      real(kind=realType),dimension(nTimeIntervalsSpectral)  :: resCL,&
           resCD, resCFx, resCFy,resCFz, resCMx, resCMy, resCMz
      
      real(kind=realType),dimension(nTimeIntervalsSpectral)  :: dPhixdot,&
           dPhiydot,dPhizdot
       
      
      
      real(kind=realType), dimension(nSections) :: t

      integer(kind=intType):: i,sps,nn

      !function definitions
      real(kind=realType)::derivativeRigidRotAngle,&
           secondDerivativeRigidRotAngle
      
      real(kind=realType):: TSAlpha,TSAlphadot
      real(kind=realType):: TSMach,TSMachdot
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!

      if(myID==0) print *,'in compute TS deriv Adj...'
      !Given
 
      if (TSpMode)then
         !P is roll
         do sps =1,nTimeIntervalsSpectral
            !compute the time of this interval
            t = timeUnsteadyRestart
            
            if(equationMode == timeSpectral) then
               do nn=1,nSections
                  !t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                  !     /         real(nTimeIntervalsSpectral,realType)
                  t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                       /         (nTimeIntervalsSpectral*1.0)
               enddo
            endif
            print *,'t',t
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
         if(myID==0)print *,'dphix',dphix
         if(myID==0)print *,'dphixdot',dphixdot
         !now compute dCl/dp
         
         call computeLeastSquaresRegression(cladj,dphix,nTimeIntervalsSpectral,dcldp,cl0)
         
         !now compute dCmz/dp
         
         call computeLeastSquaresRegression(cmzadj,dphix,nTimeIntervalsSpectral,dcmzdp,cmz0)
         if(myID==0)then
            print *,'CL estimates:'
            print *,'Clp = : ',dcldp,' cl0 = : ',cl0
            print *,'CMz estimates:'
            print *,'CMzp = : ',dcmzdp,' cmz0 = : ',cmz0
         endif
         
         ! now subtract off estimated cl,cmz and use remainder to compute 
         ! clpdot and cmzpdot.
         
         do i = 1,nTimeIntervalsSpectral
            ResCL(i) = ClAdj(i)-(dcldp*dphix(i)+cl0)
            ResCMz(i) = CmzAdj(i)-(dcmzdp*dphix(i)+cmz0)
         enddo
         if(myID==0) print *,'ResCL',resCL
         if(myID==0) print *,'resCMZ',resCMZ
         
         !now compute dCl/dpdot
         
         call computeLeastSquaresRegression(Rescl,dphixdot,nTimeIntervalsSpectral,dcldpdot,cl0dot)
         
         !now compute dCmz/dpdot
         
         call computeLeastSquaresRegression(Rescmz,dphixdot,nTimeIntervalsSpectral,dcmzdpdot,cmz0dot)
         
         if(myID==0)then
            print *,'CL estimates:'
            print *,'Clpdot = : ',dcldpdot,' cl0dot = : ',cl0dot
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
                  !t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                  !     /         real(nTimeIntervalsSpectral,realType)
                  t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                       /         (nTimeIntervalsSpectral*1.0)
               enddo
            endif
            print *,'t',t
            
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
         if(myID==0)print *,'dphiz',dphiz
         if(myID==0)print *,'dphizdot',dphizdot
         
         !now compute dCl/dq
         
         call computeLeastSquaresRegression(clAdj,dphiz,nTimeIntervalsSpectral,dcldq,cl0)
         
         !now compute dCmz/dq
         
         call computeLeastSquaresRegression(cmzAdj,dphiz,nTimeIntervalsSpectral,dcmzdq,cmz0)
         if(myID==0)then
            print *,'CL estimates:'
            print *,'Clq = : ',dcldq,' cl0 = : ',cl0
            print *,'CMz estimates:'
            print *,'CMzq = : ',dcmzdq,' cmz0 = : ',cmz0
         endif
         
         ! now subtract off estimated cl,cmz and use remainder to compute 
         ! clqdot and cmzqdot.
         
         do i = 1,nTimeIntervalsSpectral
            ResCL(i) = ClAdj(i)-(dcldq*dphiz(i)+cl0)
            ResCMz(i) = CmzAdj(i)-(dcmzdq*dphiz(i)+cmz0)
         enddo
         if(myID==0) print *,'ResCL',resCL
         if(myID==0) print *,'resCMZ',resCMZ
         
         !now compute dCl/dqdot
         
         call computeLeastSquaresRegression(Rescl,dphizdot,nTimeIntervalsSpectral,dcldqdot,cl0dot)
         
         !now compute dCmz/dqdot
         
         call computeLeastSquaresRegression(Rescmz,dphizdot,nTimeIntervalsSpectral,dcmzdqdot,cmz0dot)
         
         if(myID==0)then
            print *,'CL estimates:'
            print *,'Clqdot = : ',dcldqdot,' cl0dot = : ',cl0dot
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
                  !t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                  !     /         real(nTimeIntervalsSpectral,realType)
                  t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                       /         (nTimeIntervalsSpectral*1.0)
               enddo
            endif
            print *,'t',t
            
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
         
         if(myID==0)print *,'dphiy',dphiy
         if(myID==0)print *,'dphiydot',dphiydot
         
         !now compute dCl/dr
         
         call computeLeastSquaresRegression(clAdj,dphiy,nTimeIntervalsSpectral,dcldr,cl0)
         
         !now compute dCmz/dr
         
         call computeLeastSquaresRegression(cmzAdj,dphiy,nTimeIntervalsSpectral,dcmzdr,cmz0)
         if(myID==0)then
            print *,'CL estimates:'
            print *,'Clr = : ',dcldr,' cl0 = : ',cl0
            print *,'CMz estimates:'
            print *,'CMzr = : ',dcmzdr,' cmz0 = : ',cmz0
         endif
         
         ! now subtract off estimated cl,cmz and use remainder to compute 
         ! clrdot and cmzrdot.
         
         do i = 1,nTimeIntervalsSpectral
            ResCL(i) = ClAdj(i)-(dcldr*dphiy(i)+cl0)
            ResCMz(i) = CmzAdj(i)-(dcmzdr*dphiy(i)+cmz0)
         enddo
         if(myID==0) print *,'ResCL',resCL
         if(myID==0) print *,'resCMZ',resCMZ
         
         
      elseif(TSAlphaMode)then
         
         !compute the alphas and alphadots
         do sps =1,nTimeIntervalsSpectral
            
            !compute the time of this interval
            t = timeUnsteadyRestart
            
            if(equationMode == timeSpectral) then
               do nn=1,nSections
                  !t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                  !     /         real(nTimeIntervalsSpectral,realType)
                  t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                       /         (nTimeIntervalsSpectral*1.0)
               enddo
            endif
            print *,'t',t
            intervalAlpha(sps) = TSAlpha(degreePolAlpha,   coefPolAlpha,       &
                 degreeFourAlpha,  omegaFourAlpha,     &
                 cosCoefFourAlpha, sinCoefFourAlpha, t)
            
            intervalAlphadot(sps) = TSAlphadot(degreePolAlpha,   coefPolAlpha,       &
                 degreeFourAlpha,  omegaFourAlpha,     &
                 cosCoefFourAlpha, sinCoefFourAlpha, t)
         end do
         print *,'Alpha', intervalAlpha
         !now compute dCl/dalpha
         
         call computeLeastSquaresRegression(clAdj,intervalAlpha,nTimeIntervalsSpectral,dcldalpha,cl0)
         
         !now compute dCmz/dq
         
         call computeLeastSquaresRegression(cmzAdj,intervalAlpha,nTimeIntervalsSpectral,dcmzdalpha,cmz0)
         if(myID==0)then
            print *,'CL estimates:'
            print *,'ClAlpha = : ',dcldalpha,' cl0 = : ',cl0
            print *,'CMz estimates:'
            print *,'CMzAlpha = : ',dcmzdalpha,' cmz0 = : ',cmz0
         endif
         
         ! now subtract off estimated cl,cmz and use remainder to compute 
         ! clqdot and cmzqdot.
         
         do i = 1,nTimeIntervalsSpectral
            ResCL(i) = ClAdj(i)-(dcldalpha*IntervalAlpha(i)+cl0)
            ResCMz(i) = CmzAdj(i)-(dcmzdalpha*IntervalAlpha(i)+cmz0)
         enddo
         if(myID==0) print *,'ResCL',resCL
         if(myID==0) print *,'resCMZ',resCMZ
         !now compute dCl/dqdot
         
         call computeLeastSquaresRegression(Rescl,intervalAlphadot,nTimeIntervalsSpectral,dcldalphadot,cl0dot)
         
         !now compute dCmz/dqdot
         
         call computeLeastSquaresRegression(Rescmz,intervalAlphadot,nTimeIntervalsSpectral,dcmzdalphadot,cmz0dot)
         
         if(myID==0)then
            print *,'CL estimates:'
            print *,'ClAlphadot = : ',dcldalphadot,' cl0dot = : ',cl0dot
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
                  !t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                  !     /         real(nTimeIntervalsSpectral,realType)
                  t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                       /         (nTimeIntervalsSpectral*1.0)
               enddo
            endif
            if(myID==0)print *,'t',t
            intervalMach(sps) = machgrid+TSMach(degreePolMach,   coefPolMach,       &
                 degreeFourMach,  omegaFourMach,     &
                 cosCoefFourMach, sinCoefFourMach, t)
            
            intervalMachdot(sps) = TSMachdot(degreePolMach,   coefPolMach,       &
                 degreeFourMach,  omegaFourMach,     &
                 cosCoefFourMach, sinCoefFourMach, t)
         end do
         if(myID==0) print *,'Mach', intervalMach,machgrid
         !now compute dCl/dalpha
         
         call computeLeastSquaresRegression(clAdj,intervalMach,nTimeIntervalsSpectral,dcldMach,cl0)
         
         !now compute dCmz/dq
         
         call computeLeastSquaresRegression(cmzAdj,intervalMach,nTimeIntervalsSpectral,dcmzdMach,cmz0)
         if(myID==0)then
            print *,'CL estimates:'
            print *,'ClMach = : ',dcldMach,' cl0 = : ',cl0
            print *,'CMz estimates:'
            print *,'CMzMach = : ',dcmzdMach,' cmz0 = : ',cmz0
         endif
         
         ! now subtract off estimated cl,cmz and use remainder to compute 
         ! clqdot and cmzqdot.
         
         do i = 1,nTimeIntervalsSpectral
            ResCL(i) = ClAdj(i)-(dcldalpha*IntervalMach(i)+cl0)
            ResCMz(i) = CmzAdj(i)-(dcmzdalpha*IntervalMach(i)+cmz0)
         enddo
         if(myID==0) print *,'ResCL',resCL
         if(myID==0) print *,'resCMZ',resCMZ
         !now compute dCl/dqdot
         
         call computeLeastSquaresRegression(Rescl,intervalMachdot,nTimeIntervalsSpectral,dcldMachdot,cl0dot)
         
         !now compute dCmz/dqdot
      
         call computeLeastSquaresRegression(Rescmz,intervalAlphadot,nTimeIntervalsSpectral,dcmzdalphadot,cmz0dot)
         
         if(myID==0)then
            print *,'CL estimates:'
            print *,'ClMachdot = : ',dcldMachdot,' cl0dot = : ',cl0dot
            print *,'CMz estimates:'
            print *,'CMzMachdot = : ',dcmzdMachdot,' cmz0dot = : ',cmz0dot
         end if
      
      else
         call terminate('computeTSDerivatives','Not a valid stability motion')
      endif
      
    end subroutine computeTSStabilityDerivAdj
