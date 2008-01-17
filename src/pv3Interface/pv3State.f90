!
!      ******************************************************************
!      *                                                                *
!      * File:          pv3state.f90                                    *
!      * Author:        Juan J. Alonso                                  *
!      * Starting date: 05-15-2004                                      *
!      * Last modified: 03-24-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module PV3state
!
!      ******************************************************************
!      *                                                                *
!      * This module contains the value of the current programmed       *
!      * cutting plane being requested by PV3.  This is needed so that  *
!      * subroutine pvxyprime() in pv3Interface/pv3Routines.f90 can     *
!      * know the current cutting surface, as PV3 does not pass that    *
!      * information to it, yet it ensures that pvzprime(), where this  *
!      * value is set, is called prior to pvxyprime()                   *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
!
       integer(kind=intPV3Type) :: cutSurface
!
       end module PV3state
