!
!     ******************************************************************
!     *                                                                *
!     * File:          stabilityDerivativeDriver.f90                   *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 12-02-2009                                      *
!     * Last modified: 12-02-2009                                      *
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
!     Subroutine arguments.
! 


!
!     Local variables.
!
      real(kind=realType)::dcldp,dcldpdot,dcmzdp,dcmzdpdot         
      real(kind=realType)::dcldq,dcldqdot,dcmzdq,dcmzdqdot
      real(kind=realType)::dcldr,dcldrdot,dcmzdr,dcmzdrdot
      real(kind=realType)::dcldalpha,dcldalphadot,dcddalpha,dcmzdalpha,dcmzdalphadot
      real(kind=realType)::dcldMach,dcldMachdot,dcmzdMach,dcmzdMachdot
      real(kind=realType)::cl0,cl0dot,cD0,cmz0,cmz0dot
      

!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!

      call computeTSDerivatives(cl0,cd0,cmz0,dcldalpha,dcddalpha,dcmzdalpha,&
           dcmzdalphadot,dcmzdq)


    end subroutine stabilityDerivativeDriver
