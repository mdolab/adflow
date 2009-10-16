!
!      ******************************************************************
!      *                                                                *
!      * File:          computeTSdetrivatives.f90                       *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 09-15-2009                                      *
!      * Last modified: 09-15-2009                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine computeTSDerivatives
!
!      ******************************************************************
!      *                                                                *
!      * Computes an estimate for the stability derivatives from the    *
!      * coefficient values of the Time Spectral Solution.              *
!      * for now just the q derivative...                               *
!      ******************************************************************
!
         use precision
         use inputMotion
         use inputPhysics
         use inputTimeSpectral
         use monitor
         use section
         implicit none

         !Subroutine arguments



         !local Variables
         real(kind=realType),dimension(nTimeIntervalsSpectral)  :: CL, CD,&
              CFx, CFy,CFz, CMx, CMy, CMz

         real(kind=realType),dimension(nTimeIntervalsSpectral)  :: dPhix,&
              dPhiy,dPhiz

         real(kind=realType),dimension(nTimeIntervalsSpectral)  :: resCL,&
              resCD, resCFx, resCFy,resCFz, resCMx, resCMy, resCMz

         real(kind=realType),dimension(nTimeIntervalsSpectral)  :: dPhixdot,&
              dPhiydot,dPhizdot
         
         real(kind=realType)::dcldq,dcldqdot,dcmzdq,dcmzdqdot
         real(kind=realType)::cl0,cl0dot,cmz0,cmz0dot

         real(kind=realType), dimension(nSections) :: t
         integer(kind=inttype)::level,sps,i,nn
         real(kind=realType)::derivativeRigidRotAngle,&
              secondDerivativeRigidRotAngle
         !real(kind=realType)::t
         !begin execution

         print *,'in compute TS deriv...'
         !Compute and store the aero coef. Values for each TS level
         do sps =1,nTimeIntervalsSpectral
            print *,'sps',sps
            level = 1
            call computeAeroCoef(CL(sps),CD(sps),CFx(sps),CFy(sps),CFz(sps),&
                 CMx(sps),CMy(sps),CMz(sps),level,sps)
            print *,'CL',cl
            t = timeUnsteadyRestart
            
            if(equationMode == timeSpectral) then
               do nn=1,nSections
                  t(nn) = t(nn) + (sps-1)*sections(nn)%timePeriod &
                       /         real(nTimeIntervalsSpectral,realType)
               enddo
            endif
            ! Compute the time derivative of the rotation angles around the
            ! x-axis, y-axis and z-axis. i.e. compute p,q,r
            
            dphiX(sps) = derivativeRigidRotAngle(degreePolXRot,   &
                 coefPolXRot,     &
                 degreeFourXRot,  &
                 omegaFourXRot,   &
                 cosCoefFourXRot, &
                 sinCoefFourXRot, t)
            print *,'dphix',dphix
            
            dphiY(sps) = derivativeRigidRotAngle(degreePolYRot,   &
                 coefPolYRot,     &
                 degreeFourYRot,  &
                 omegaFourYRot,   &
                 cosCoefFourYRot, &
                 sinCoefFourYRot, t)
            
            print *,'dphiy',dphiy
            dphiZ(sps) = derivativeRigidRotAngle(degreePolZRot,   &
                 coefPolZRot,     &
                 degreeFourZRot,  &
                 omegaFourZRot,   &
                 cosCoefFourZRot, &
                 sinCoefFourZRot, t)
            print *,'dphiz',dphiz
            
            !add in qdot computation here!
            dphiXdot(sps) = secondDerivativeRigidRotAngle(degreePolXRot,   &
                 coefPolXRot,     &
                 degreeFourXRot,  &
                 omegaFourXRot,   &
                 cosCoefFourXRot, &
                 sinCoefFourXRot, t)
            print *,'dphixdot',dphixdot
            
            dphiYdot(sps) = secondDerivativeRigidRotAngle(degreePolYRot,   &
                 coefPolYRot,     &
                 degreeFourYRot,  &
                 omegaFourYRot,   &
                 cosCoefFourYRot, &
                 sinCoefFourYRot, t)
            
            print *,'dphiydot',dphiydot
            dphiZdot(sps) = secondDerivativeRigidRotAngle(degreePolZRot,   &
                 coefPolZRot,     &
                 degreeFourZRot,  &
                 omegaFourZRot,   &
                 cosCoefFourZRot, &
                 sinCoefFourZRot, t)
            print *,'dphizdot',dphizdot
         enddo


         !now compute dCl/dq
         
         call computeLeastSquaresRegression(cl,dphiz,nTimeIntervalsSpectral,dcldq,cl0)

         !now compute dCmz/dq
         
         call computeLeastSquaresRegression(cmz,dphiz,nTimeIntervalsSpectral,dcmzdq,cmz0)

         print *,'CL estimates:'
         print *,'Clq = : ',dcldq,' cl0 = : ',cl0
         print *,'CMz estimates:'
         print *,'CMzq = : ',dcmzdq,' cmz0 = : ',cmz0

         ! now subtract off estimated cl,cmz and use remainder to compute 
         ! clqdot and cmzqdot.

         do i = 1,nTimeIntervalsSpectral
            ResCL(i) = Cl(i)-(dcldq*dphiz(i)+cl0)
            ResCMz(i) = Cmz(i)-(dcmzdq*dphiz(i)+cmz0)
         enddo

         !now compute dCl/dqdot
         
         call computeLeastSquaresRegression(Rescl,dphizdot,nTimeIntervalsSpectral,dcldqdot,cl0dot)

         !now compute dCmz/dqdot
         
         call computeLeastSquaresRegression(Rescmz,dphizdot,nTimeIntervalsSpectral,dcmzdqdot,cmz0dot)
         
         print *,'CL estimates:'
         print *,'Clqdot = : ',dcldqdot,' cl0dot = : ',cl0dot
         print *,'CMz estimates:'
         print *,'CMzqdot = : ',dcmzdqdot,' cmz0dot = : ',cmz0dot


       end subroutine computeTSDerivatives
