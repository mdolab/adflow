!
!      ******************************************************************
!      *                                                                *
!      * File:          saNuFromEddyRatio.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-22-2003                                      *
!      * Last modified: 04-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       function saNuKnownEddyRatio(eddyRatio, nuLam)
!
!      ******************************************************************
!      *                                                                *
!      * saNuKnownEddyRatio computes the Spalart-Allmaras transport     *
!      * variable nu for the given eddy viscosity ratio.                *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use paramTurb
       implicit none
!
!      Function type.
!
       real(kind=realType) :: saNuKnownEddyRatio
!
!      Function arguments.
!
       real(kind=realType), intent(in) :: eddyRatio, nuLam
!
!      Local variables.
!
       real(kind=realType) :: cv13, chi, chi2, chi3, chi4, f, df, dchi
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Take care of the exceptional cases.

       if(eddyRatio <= zero) then
         saNuKnownEddyRatio = zero
         return
       endif

       ! Set the value of cv1^3, which is the constant appearing in the
       ! sa function fv1 to compute the eddy viscosity

       cv13 = rsaCv1**3

       ! Determine the value of chi, which is given by the quartic
       ! polynomial chi^4 - ratio*(chi^3 + cv1^3) = 0.
       ! First determine the start value, depending on the eddyRatio.

       if(eddyRatio < 1.e-4_realType) then
         chi = 0.5_realType
       else if(eddyRatio < 1.0_realType) then
         chi = 5.0_realType
       else if(eddyRatio < 10.0_realType) then
         chi = 10.0_realType
       else
         chi = eddyRatio
       endif

       ! The actual newton algorithm.

       do
         ! Compute the function value and the derivative.

         chi2 = chi*chi
         chi3 = chi*chi2
         chi4 = chi*chi3

         f  = chi4 - eddyRatio*(chi3 + cv13)
         df = four*chi3 - three*eddyRatio*chi2

         ! Compute the negative update and the new value of chi.

         dchi = f/df
         chi  = chi - dchi

         ! Condition to exit the loop.

         if(abs(dchi/chi) <= thresholdReal) exit
       enddo

       ! Chi is the ratio of the spalart allmaras transport variable and
       ! the laminar viscosity. So multiply chi with the laminar viscosity
       ! to obtain the correct value.

       saNuKnownEddyRatio = nuLam*chi

       end function saNuKnownEddyRatio
