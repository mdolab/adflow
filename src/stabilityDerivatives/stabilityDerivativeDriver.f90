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
!     Subroutine arguments.
! 


!
!     Local variables.
! 
      real(kind=realType),dimension(8)::dcdalpha,dcdalphadot,dcdbeta,&
           dcdbetadot,dcdMach,dcdMachdot
      real(kind=realType),dimension(8)::dcdp,dcdpdot,dcdq,dcdqdot,dcdr,dcdrdot
      real(kind=realType),dimension(8)::Coef0,Coef0dot
!!$      real(kind=realType)::dcldp,dcldpdot,dcmzdp,dcmzdpdot         
!!$      real(kind=realType)::dcldq,dcldqdot,dcddq,dcddqdot,dcmzdq,dcmzdqdot
!!$      real(kind=realType)::dcldr,dcldrdot,dcmzdr,dcmzdrdot
!!$      real(kind=realType)::dcldalpha,dcldalphadot,dcddalpha,dcddalphadot,dcmzdalpha,dcmzdalphadot
!!$      real(kind=realType)::dcldMach,dcldMachdot,dcmzdMach,dcmzdMachdot
!!$      real(kind=realType)::cl0,cl0dot,cD0,cmz0,cmz0dot
      

!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
      call computeTSDerivatives(coef0,dcdalpha,dcdalphadot,dcdq,dcdqdot)
      !(cl0,cd0,cmz0,dcldalpha,dcddalpha,&
      !     dcmzdalpha,dcldalphadot,dcddalphadot,dcmzdalphadot,dcldq,&
      !     dcddq,dcmzdq,dcldqdot,dcddqdot,dcmzdqdot)

    end subroutine stabilityDerivativeDriver
