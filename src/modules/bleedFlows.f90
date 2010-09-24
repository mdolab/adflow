!
!      ******************************************************************
!      *                                                                *
!      * File:          bleedFlows.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-11-2005                                      *
!      * Last modified: 09-12-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       module bleedFlows
!
!      ******************************************************************
!      *                                                                *
!      * Module which contains the derived data types as well as the    *
!      * corresponding arrays to store the information needed to model  *
!      * the bleed flows. Both inflow bleeds and outflow bleeds are     *
!      * possible.                                                      *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type for the flow bleeds.   *
!      *                                                                *
!      ******************************************************************
!
       type bleedflowType

         ! familyID:    Corresponding family ID.
         ! massFlux:    Prescribed mass flux, 
         ! curMassFlux: Current mass flux. When converged this should
         !              be equal to massFlux.

         integer(kind=intType) :: familyID
         real(kind=realType)   :: massFlux, curMassFlux

       end type bleedflowType
!
!      ******************************************************************
!      *                                                                *
!      * Variables stored in this module.                               *
!      *                                                                *
!      ******************************************************************
!
       ! nInflowBleeds:    Number of inflow bleeds present.
       ! nOutflowBleeds:   Number of outflow bleeds present.
       ! inflowBleeds(:):  Array with the information for the inflow
       !                   bleeds.
       ! outflowBleeds(:): Array with the information for the outflow
       !                   bleeds.

       integer(kind=intType) :: nInflowBleeds
       integer(kind=intType) :: nOutflowBleeds

       type(bleedflowType), dimension(:), allocatable :: inflowBleeds
       type(bleedflowType), dimension(:), allocatable :: outflowBleeds

       end module bleedFlows
