!
!      ******************************************************************
!      *                                                                *
!      * File:          setFlowInfinityState.f90                        *
!      * Author:        Edwin van der Weide, Georgi Kalitzin            *
!      * Starting date: 02-21-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setFlowInfinityState
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
       use inputPhysics
       use paramTurb
       implicit none
!
!      Local variables
!
       integer :: ierr

       real(kind=realType) :: nuInf, ktmp, uInf2
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

       uInf2 = MachCoef*MachCoef*gammaInf*pInf/rhoInf

       ! Allocate the memory for wInf.

       allocate(wInf(nw), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("setFlowReferenceState", &
                        "Memory allocation failure for wInf")

       ! Set the reference value of the flow variables, except the total
       ! energy. This will be computed at the end of this routine.

       wInf(irho) = rhoInf
       wInf(ivx)  = uInf*velDirFreestream(1)
       wInf(ivy)  = uInf*velDirFreestream(2)
       wInf(ivz)  = uInf*velDirFreestream(3)

       ! Set the turbulent variables if transport variables are
       ! to be solved.

       if(equations == RANSEquations) then

         nuInf  = muInf/rhoInf

         select case(turbModel)

           case (spalartAllmaras, spalartAllmarasEdwards)

             wInf(itu1) = saNuKnownEddyRatio(eddyVisInfRatio, nuInf)
         !   wInf(itu1) = 1.341946*nuInf   ! eddyVis = 0.009*lamVis

           !=============================================================

           case (komegaWilcox, komegaModified, menterSST)

             wInf(itu1) = 1.5_realType*uInf2*turbIntensityInf**2
             wInf(itu2) = wInf(itu1)/(eddyVisInfRatio*nuInf)

           !=============================================================

           case (ktau)

             wInf(itu1) = 1.5_realType*uInf2*turbIntensityInf**2
             wInf(itu2) = eddyVisInfRatio*nuInf/wInf(itu1)

           !=============================================================

           case (v2f)

             wInf(itu1) = 1.5_realType*uInf2*turbIntensityInf**2
             wInf(itu2) = 0.09_realType*wInf(itu1)**2 &
                        / (eddyVisInfRatio*nuInf)
             wInf(itu3) = 0.666666_realType*wInf(itu1)
             wInf(itu4) = 0.0_realType

         end select

       endif

       ! Set the value of pInfCorr. In case a k-equation is present
       ! add 2/3 times rho*k.

       pInfCorr = pInf
       if( kPresent ) pInfCorr = pInf + two*third*rhoInf*wInf(itu1)

       ! Compute the free stream total energy.

       ktmp = zero
       if( kPresent ) ktmp = wInf(itu1)
       call etotArray(rhoInf, uInf, zero, zero, pInfCorr, ktmp, &
                       wInf(irhoE), kPresent, 1_intType)

       end subroutine setFlowInfinityState
