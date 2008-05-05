!
!      ******************************************************************
!      *                                                                *
!      * File:          section.f90                                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-09-2003                                      *
!      * Last modified: 06-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module section
!
!      ******************************************************************
!      *                                                                *
!      * This module contains the definition of the derived data type   *
!      * sectionType, which stores the information of a section of the  *
!      * grid. It also contains the array to store the info of all the  *
!      * sections.                                                      *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type sectionType.          *
!      *                                                                *
!      ******************************************************************
!
       type sectionType

         ! periodic:       Whether or not the section is periodic.
         ! rotating:       Whether or not the section is rotating.
         ! nSlices:        Number of periodic slices to obtain
         !                 a full rotation.
         ! timePeriod:     The physical time of one period.
         ! rotCenter(3):   Cartesian coordinates of the rotation center.
         ! rotMatrix(3,3): Rotation matrix of the periodic
         !                 transformation.
         ! translation(3): Translation vector of the periodic
         !                 transformation.
         ! rotAxis(3):     The rotation axis.
         ! rotRate(3):     The rotation rate of the section.

         logical :: periodic, rotating

         integer(kind=intType) :: nSlices

         real(kind=realType) :: timePeriod

         real(kind=realType), dimension(3)   :: rotCenter, translation
         real(kind=realType), dimension(3)   :: rotAxis, rotRate
         real(kind=realType), dimension(3,3) :: rotMatrix

       end type sectionType

       ! nSections:           Number of different sections in the grid.
       ! sections(nSections): The info of the corresponding sections.

       integer(kind=intType) :: nSections

       type(sectionType), dimension(:), allocatable :: sections

       end module section
