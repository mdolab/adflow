!
!     ******************************************************************
!     *                                                                *
!     * File:          stabilityDerivativeDriver.f90                   *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 12-02-2009                                      *
!     * Last modified: 02-14-2011                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine stabilityDerivativeDriver
!
!     ******************************************************************
!     *                                                                *
!     * Runs the Time spectral stability derivative routines from the  *
!     * main program file                                              *
!     *                                                                *
!     ******************************************************************
!
      use precision
      implicit none
!
!     Local variables.
! 
      real(kind=realType),dimension(8)::dcdalpha,dcdalphadot,dcdbeta,&
           dcdbetadot,dcdMach,dcdMachdot
      real(kind=realType),dimension(8)::dcdp,dcdpdot,dcdq,dcdqdot,dcdr,dcdrdot
      real(kind=realType),dimension(8)::Coef0,Coef0dot

      !call computeTSDerivatives(coef0,dcdalpha,dcdalphadot,dcdq,dcdqdot)
 
   end subroutine stabilityDerivativeDriver
