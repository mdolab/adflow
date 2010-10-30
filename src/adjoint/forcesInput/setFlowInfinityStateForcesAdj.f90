!
!      ******************************************************************
!      *                                                                *
!      * File:          setFlowInfinityStateAdj.f90                     *
!      * Author:        Edwin van der Weide, Georgi Kalitzin            *
!      *                C.A.(Sandy) Mader
!      * Starting date: 02-21-2003                                      *
!      * Last modified: 05-14-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setFlowInfinityStateForcesAdj(velDirFreestreamAdj,liftDirectionAdj,&
            dragDirectionAdj, Machadj, MachCoefAdj,uInfAdj,wInfAdj,prefAdj,&
            rhorefAdj, pinfdimAdj, rhoinfdimAdj, rhoinfAdj, pinfAdj,&
            murefAdj, timerefAdj,pInfCorrAdj)
!
!      ******************************************************************
!      *                                                                *
!      * setFlowInfinityState sets the free stream state vector of      *
!      * the flow variables. If nothing is specified for each of the    *
!      * farfield boundaries, these values will be taken to define the  *
!      * free stream.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use flowVarRefState
!!$       use inputPhysics
!!$       use paramTurb
       implicit none

!
!      subroutine Variables
!
       
       real(kind=realType), dimension(3), intent(inout) :: velDirFreestreamAdj
       real(kind=realType), dimension(3), intent(inout) :: liftDirectionAdj
       real(kind=realType), dimension(3), intent(inout) :: dragDirectionAdj
       real(kind=realType), intent(inout) :: MachAdj,MachCoefAdj,uInfAdj,pInfCorrAdj
       REAL(KIND=REALTYPE) :: prefAdj, rhorefAdj
       REAL(KIND=REALTYPE) :: pinfdimAdj, rhoinfdimAdj
       REAL(KIND=REALTYPE) :: rhoinfAdj, pinfAdj
       REAL(KIND=REALTYPE) :: murefAdj, timerefAdj
       real(kind=realType), dimension(nw),intent(out)::wInfAdj

!
!      Local variables
!
       integer :: ierr

       real(kind=realType) :: nuInf, ktmp, uInf2, vtmp, wtmp
!
!      Function definition
!
       real(kind=realType) :: saNuKnownEddyRatio
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Compute the velocity squared based on MachCoef;
       ! needed for the initialization of the turbulent energy,
       ! especially for moving geometries.

       uInf2 = MachCoefAdj*MachCoefAdj*gammaInf*pInf/rhoInfAdj

!!$       ! Allocate the memory for wInf.
!!$
!!$       allocate(wInf(nw), stat=ierr)
!!$       if(ierr /= 0)                             &
!!$         call terminate("setFlowReferenceState", &
!!$                        "Memory allocation failure for wInf")

       ! Set the reference value of the flow variables, except the total
       ! energy. This will be computed at the end of this routine.

       wInfAdj(irho) = rhoInfAdj
       wInfAdj(ivx)  = uInfAdj*velDirFreestreamAdj(1)
       wInfAdj(ivy)  = uInfAdj*velDirFreestreamAdj(2)
       wInfAdj(ivz)  = uInfAdj*velDirFreestreamAdj(3)

!!$       ! Set the turbulent variables if transport variables are
!!$       ! to be solved.
!!$
!!$       if(equations == RANSEquations) then
!!$
!!$         nuInf  = muInf/rhoInf
!!$
!!$         select case(turbModel)
!!$
!!$           case (spalartAllmaras, spalartAllmarasEdwards)
!!$
!!$             wInf(itu1) = saNuKnownEddyRatio(eddyVisInfRatio, nuInf)
!!$         !   wInf(itu1) = 1.341946*nuInf   ! eddyVis = 0.009*lamVis
!!$
!!$           !=============================================================
!!$
!!$           case (komegaWilcox, komegaModified, menterSST)
!!$
!!$             wInf(itu1) = 1.5_realType*uInf2*turbIntensityInf**2
!!$             wInf(itu2) = wInf(itu1)/(eddyVisInfRatio*nuInf)
!!$
!!$           !=============================================================
!!$
!!$           case (ktau)
!!$
!!$             wInf(itu1) = 1.5_realType*uInf2*turbIntensityInf**2
!!$             wInf(itu2) = eddyVisInfRatio*nuInf/wInf(itu1)
!!$
!!$           !=============================================================
!!$
!!$           case (v2f)
!!$
!!$             wInf(itu1) = 1.5_realType*uInf2*turbIntensityInf**2
!!$             wInf(itu2) = 0.09_realType*wInf(itu1)**2 &
!!$                        / (eddyVisInfRatio*nuInf)
!!$             wInf(itu3) = 0.666666_realType*wInf(itu1)
!!$             wInf(itu4) = 0.0_realType
!!$
!!$         end select
!!$
!!$       endif

       ! Set the value of pInfCorr. In case a k-equation is present
       ! add 2/3 times rho*k.

       pInfCorrAdj = pInf
       if( kPresent ) pInfCorrAdj = pInf + two*third*rhoInfAdj*wInfAdj(itu1)

       ! Compute the free stream total energy.

       ktmp = 0.0!zero
 !      if( kPresent ) ktmp = wInfAdj(itu1)

       vtmp = 0.0
       wtmp = 0.0
!       call etotArrayAdj(rhoInf, uInfAdj, zero, zero, pInfCorr, ktmp, &
!                       wInfAdj(irhoE), kPresent, 1_intType)
       call etotArrayForcesAdj(rhoInfAdj, uInfAdj, vtmp, wtmp, pInfCorrAdj, ktmp, &
                   wInfAdj(irhoE), kPresent, 1)

     end subroutine setFlowInfinityStateForcesAdj
