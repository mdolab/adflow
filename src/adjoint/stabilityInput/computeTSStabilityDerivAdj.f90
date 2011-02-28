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
      subroutine computeTSStabilityDerivAdj(BaseCoef,coef0,dcdalpha,&
           dcdalphadot,dcdq,dcdqdot)
!!$        computeTSStabilityDerivAdj(cFxAdj,cFyAdj,cFzAdj,cMxAdj,&
!!$                       cMyAdj,cMzAdj,CLAdj,CDAdj,&
!!$                       cl0,cd0,cmz0,dcldalpha,dcddalpha,dcmzdalpha,&
!!$                       dcmzdalphadot,dcmzdq)
 
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
      use flowvarrefstate      !timeref
      use monitor              !timeunsteadyrestart,timeunsteady
      use section              !nsections, sections%
      use inputMotion          ! degreePol*
!      use bcTypes              ! EulerWall, ...
!      use flowvarrefstate      !nw

      implicit none

!
!     Subroutine arguments.
!
      real(kind=realType),dimension(nTimeIntervalsSpectral,8)::BaseCoef
!!$      real(kind=realType), dimension(nTimeIntervalsSpectral)::       &
!!$                             ClAdj,CdAdj,CfxAdj,CfyAdj,CfzAdj,   &
!!$                             CmxAdj,CmyAdj,CmzAdj
     
!!$      real(kind=realType)::dcldp,dcldpdot,dcddp,dcddpdot,dcmzdp,dcmzdpdot
!!$      real(kind=realType)::dcldq,dcldqdot,dcddq,dcddqdot,dcmzdq,dcmzdqdot
!!$      real(kind=realType)::dcldr,dcldrdot,dcddr,dcddrdot,dcmzdr,dcmzdrdot
!!$      real(kind=realType)::dcldalpha,dcldalphadot,dcddalpha,dcddalphadot,dcmzdalpha,dcmzdalphadot
!!$      real(kind=realType)::dcldMach,dcldMachdot,dcddMach,dcddMachdot,dcmzdMach,dcmzdMachdot
!!$      real(kind=realType)::cl0,cl0dot,cd0,cd0dot,cmz0,cmz0dot
      real(kind=realType),dimension(8)::dcdp,dcdpdot,dcdq,dcdqdot,dcdr,dcdrdot
      real(kind=realType),dimension(8)::dcdalpha,dcdalphadot,dcdbeta,dcdbetadot,dcdMach,dcdMachdot
      real(kind=realType),dimension(8)::Coef0,Coef0dot
      real(kind=realType),dimension(nTimeIntervalsSpectral,8)::ResBaseCoef
!
!     Local variables.
!
     

      !Stability derivative variables
      real(kind=realType),dimension(nTimeIntervalsSpectral)  :: dPhix,&
           dPhiy,dPhiz
      
      real(kind=realType),dimension(nTimeIntervalsSpectral)  :: intervalAlpha,intervalAlphadot
      real(kind=realType),dimension(nTimeIntervalsSpectral)  :: intervalMach,intervalMachdot
      
!!$      real(kind=realType),dimension(nTimeIntervalsSpectral)  :: resCL,&
!!$           resCD, resCFx, resCFy,resCFz, resCMx, resCMy, resCMz
      
      real(kind=realType),dimension(nTimeIntervalsSpectral)  :: dPhixdot,&
           dPhiydot,dPhizdot
       
      
      
      real(kind=realType), dimension(nSections) :: t

      integer(kind=intType):: i,sps,nn

      !function definitions
      real(kind=realType)::derivativeRigidRotAngle,&
           secondDerivativeRigidRotAngle
      
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
         do i =1,8
            call computeLeastSquaresRegression(BaseCoef(:,i),dphix,nTimeIntervalsSpectral,dcdp(i),coef0(i))
         end do
        
        
         
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
!!$
!!$         do i = 1,nTimeIntervalsSpectral
!!$            ResCL(i) = ClAdj(i)-(dcldp*dphix(i)+cl0)
!!$            ResCd(i) = CdAdj(i)-(dcddp*dphix(i)+cd0)
!!$            ResCMz(i) = CmzAdj(i)-(dcmzdp*dphix(i)+cmz0)
!!$         enddo
!!$         if(myID==0) print *,'ResCL',resCL
!!$         if(myID==0) print *,'ResCd',resCd
!!$         if(myID==0) print *,'resCMZ',resCMZ
         
          !now compute dCl/dpdot
         do i = 1,8
            call computeLeastSquaresRegression(ResBaseCoef(:,i),dphixdot,nTimeIntervalsSpectral,dcdpdot(i),Coef0dot(i))
         enddo
         
!!$         if(myID==0)then
!!$            print *,'CL estimates:'
!!$            print *,'Clpdot = : ',dcldpdot,' cl0dot = : ',cl0dot
!!$            print *,'Cd estimates:'
!!$            print *,'Cdpdot = : ',dcddpdot,' cd0dot = : ',cd0dot
!!$            print *,'CMz estimates:'
!!$            print *,'CMzpdot = : ',dcmzdpdot,' cmz0dot = : ',cmz0dot
!!$         end if
         
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
!!$         if(myID==0)print *,'dphiz',dphiz
!!$         if(myID==0)print *,'dphizdot',dphizdot
         
         !now compute dCl/dq
         do i =1,8
            call computeLeastSquaresRegression(BaseCoef(:,i),dphiz,nTimeIntervalsSpectral,dcdq(i),coef0(i))
         end do
       
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
     
         !now normalize the results...
         dcdq = dcdq*timeRef*2*(machGrid*a)/lengthRef
    
         
         !now compute dCl/dpdot
         do i = 1,8
            call computeLeastSquaresRegression(ResBaseCoef(:,i),dphizdot,nTimeIntervalsSpectral,dcdqdot(i),Coef0dot(i))
         enddo
!!$         if(myID==0)then
!!$            print *,'CL estimates:'
!!$            print *,'Clqdot = : ',dcldqdot,' cl0dot = : ',cl0dot
!!$            print *,'Cd estimates:'
!!$            print *,'Cdqdot = : ',dcddqdot,' cd0dot = : ',cd0dot
!!$            print *,'CMz estimates:'
!!$            print *,'CMzqdot = : ',dcmzdqdot,' cmz0dot = : ',cmz0dot
!!$         end if
         
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
         do i =1,8
            call computeLeastSquaresRegression(BaseCoef(:,i),dphiy,nTimeIntervalsSpectral,dcdr(i),coef0(i))
         end do
        
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
   
      elseif(TSBetaMode)then
         print *,'Beta mode not yet implemented'   
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
         do i =1,8
            call computeLeastSquaresRegression(BaseCoef(:,i),intervalMach,nTimeIntervalsSpectral,dcdMach(i),coef0(i))
         end do

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
                
      else
         call terminate('computeTSDerivatives','Not a valid stability motion')
      endif
      
    end subroutine computeTSStabilityDerivAdj
