!
!      ******************************************************************
!      *                                                                *
!      * File:          commMixing.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-27-2005                                      *
!      * Last modified: 03-21-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module commMixing
!
!      ******************************************************************
!      *                                                                *
!      * This module contains the derived data type to store the        *
!      * interpolation information for mixing plane halo cells.         *
!      * A three-dimensional array of this derived type to store the    *
!      * info for all grid levels and sliding mesh interfaces is also   *
!      * defined in this module.                                        *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the data type comm_mixing_type, which stores *
!      * the interpolation data of both sides of a sliding interface.   *
!      *                                                                *
!      ******************************************************************
!
       type commMixingType

         ! nInter: number of interpolation points in the list

         integer(kind=intType) :: nInter
!
!        ****************************************************************
!        *                                                              *
!        * The donor data.                                              *
!        *                                                              *
!        ****************************************************************
!
         ! nDonor: number of local donor cells.

         integer(kind=intType) :: nDonor

         ! blockDonor(nDonor): The block ID of the donor cells

         integer(kind=intType), dimension(:), pointer :: blockDonor

         ! indD(nDonor,3,2): The i,j,k indices of the donor cells.
         !                   There is a third index with dimension 2
         !                   to store the 1st and thesecond level
         !                   donor cells.

         integer(kind=intType), dimension(:,:,:), pointer :: indD

         ! rotMatDonor(nDonor,3,3): The rotation matrix, which
         !                          transforms the cartesian velocity
         !                          components of the donors into the
         !                          componenents of the local
         !                          cylindrical frame.

         real(kind=realType), dimension(:,:,:), pointer :: rotMatDonor

         ! nIntervalsDonor(0:nDonor): Number of interpolation intervals
         !                            to which every donor contributes.
         !                            This array is in cumulative
         !                            storage format.

         integer(kind=intType), dimension(:), pointer :: nIntervalsDonor

         ! indListDonor(nn): The point in the list of interpolated
         !                   data to which the donor contributes.
         !                   The size of this array is given by nn,
         !                   which equals nIntervalsDonor(nDonor).

         integer(kind=intType), dimension(:), pointer :: indListDonor

         ! weightDonor(nn): the interpolation weight of the donors
         !                  to the corresponding point in the list.

         real(kind=realType), dimension(:), pointer :: weightDonor
!
!        ****************************************************************
!        *                                                              *
!        * The halo data.                                               *
!        *                                                              *
!        ****************************************************************
!
         ! nHalo: Number of local halo cells.

         integer(kind=intType) :: nHalo

         ! indListHalo(nHalo,2): The two points in the list of
         !                       interpolated data from which the halo
         !                       value is interpolated.

         integer(kind=intType), dimension(:,:), pointer :: indListHalo

         ! blockHalo(nHalo): The block ID of the halo cells.

         integer(kind=intType), dimension(:), pointer :: blockHalo

         ! indH(nHalo,3,2): The i,j,k indices of the halo cells. There
         !                  is a third index with dimension 2 to store
         !                  the 1st and 2nd level halo cells.

         integer(kind=intType), dimension(:,:,:), pointer :: indH

         ! weightHalo(nHalo,2): The weights of the two points in the
         !                      list of interpolated data from which
         !                      the halo value is interpolated.

         real(kind=realType), dimension(:,:), pointer :: weightHalo

         ! rotMatHalo(nHalo,3,3): The rotation matrix, which
         !                        transforms the cylindrical velocity
         !                        components of the halo cells into the
         !                        componenents of the cartesian frame.

         real(kind=realType), dimension(:,:,:), pointer :: rotMatHalo

       end type commMixingType
!
!      ******************************************************************
!      *                                                                *
!      * Variables stored in this module.                               *
!      *                                                                *
!      ******************************************************************
!
       ! commPatternMixing(nLevel,nInterface_groups,2): The 3D array to
       !                                                store all the
       !                                                mixing plane
       !                                                interpolations.

       type(commMixingType), allocatable, dimension(:,:,:) :: &
                                                      commPatternMixing

       end module commMixing
