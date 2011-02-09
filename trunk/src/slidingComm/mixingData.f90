!
!      ******************************************************************
!      *                                                                *
!      * File:          mixingData.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-28-2005                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module mixingData
!
!      ******************************************************************
!      *                                                                *
!      * Local module which defines/stores variables used to determine  *
!      * communication/interpolation pattern of a mixing plane.         *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * Definition of the data type used to construct the              *
!      * interpolation interval.                                        *
!      *                                                                *
!      ******************************************************************
!
       type mixingIntervalType

         ! rMin:        The minimum of the two radial coordinates which
         !              define the edge.
         ! rMax:        Idem for the maximum.
         ! regularEdge: Whether or not the edge is a regular one.
         !              Regular means that the edge is in between
         !              opposite edges of the corresponding
         !              quadrilateral face.

         real(kind=realType) :: rMin, rMax

         logical :: regularEdge

       end type mixingIntervalType

       ! Interface for the extension of the operators <= and <.
       ! These are needed for the sorting of mixingIntervalType.
       ! Note that the = operator does not need to be defined, because
       ! mixingIntervalType only contains primitive types.

       interface operator(<=)
         module procedure lessEqualMixingInterval
       end interface

       interface operator(<)
         module procedure lessMixingInterval
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Variables stored in this module.                               *
!      *                                                                *
!      ******************************************************************
!
       ! radialInterface: Whether this is a radial or an axial
       !                  interface.

       logical :: radialInterface

       ! nMixingPoints:     Number of points used to define the
       !                    interpolation intervals.
       ! mixingPoints(nn):  The corresponding points; nn equals
       !                    nMixingPoints.
       ! mixingCells(nn-1): The coordinates of the centers of the
       !                    interpolation intervals.

       integer(kind=intType) :: nMixingPoints

       real(kind=realType), dimension(:), allocatable :: mixingPoints
       real(kind=realType), dimension(:), allocatable :: mixingCells

       !=================================================================

       contains

         !===============================================================

         logical function lessEqualMixingInterval(g1,g2)
!
!        ****************************************************************
!        *                                                              *
!        * LessEqualMixingInterval returns .true. if g1 <= g2 and       *
!        * .false. otherwise. The comparison is firstly based on rMin,  *
!        * followed by rMax and finally on regularEdge.                 *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(mixingIntervalType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Compare rMin.

         if(g1%rMin < g2%rMin) then
           lessEqualMixingInterval = .true.
           return
         else if(g1%rMin > g2%rMin) then
           lessEqualMixingInterval = .false.
           return
         endif

         ! rMin's are equal. Compare rMax.

         if(g1%rMax < g2%rMax) then
           lessEqualMixingInterval = .true.
           return
         else if(g1%rMax > g2%rMax) then
           lessEqualMixingInterval = .false.
           return
         endif

         ! Finally compare the logicals regularEdge.

         if(g1%regularEdge .and. (.not. g2%regularEdge)) then
           lessEqualMixingInterval = .false.
           return
         endif

         ! In all remaining cases g1 <= g2.
         ! Set lessEqualMixingInterval to .true.

         lessEqualMixingInterval = .true.

         end function lessEqualMixingInterval

         !===============================================================

         logical function lessMixingInterval(g1,g2)
!
!        ****************************************************************
!        *                                                              *
!        * LessMixingInterval returns .true. If g1 < g2 and .false.     *
!        * otherwise. The comparison is firstly based on rMin, followed *
!        * by rMax and finally on regularEdge.                          *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(mixingIntervalType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Compare rMin.

         if(g1%rMin < g2%rMin) then
           lessMixingInterval = .true.
           return
         else if(g1%rMin > g2%rMin) then
           lessMixingInterval = .false.
           return
         endif

         ! rMin's are equal. Compare rMax.

         if(g1%rMax < g2%rMax) then
           lessMixingInterval = .true.
           return
         else if(g1%rMax > g2%rMax) then
           lessMixingInterval = .false.
           return
         endif

         ! Finally compare the logicals regularEdge.

         if((.not. g1%regularEdge) .and. g2%regularEdge) then
           lessMixingInterval = .true.
           return
         endif

         ! In all remaining cases g1 >= g2.
         ! Set lessMixingInterval to .False.

         lessMixingInterval = .true.

         end function lessMixingInterval

       end module mixingData
