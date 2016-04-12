!
!      ******************************************************************
!      *                                                                *
!      * File:          extraOutput.f90                                 *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 03-25-2003                                      *
!      * Last modified: 07-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module extraOutput
!
!      ******************************************************************
!      *                                                                *
!      * This module contains the logicals which define the variables   *
!      * to be written to the solution file. Both the surface variables *
!      * to be written as well as the extra volume variables are stored *
!      * in this module.                                                *
!      *                                                                *
!      ******************************************************************
!
        use constants
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * The logical variables, which define the surface variables to   *
!      * be written.                                                    *
!      *                                                                *
!      ******************************************************************
!
       logical :: surfWriteRho,   surfWriteP,        surfWriteTemp
       logical :: surfWriteVx,    surfWriteVy,       surfWriteVz
       logical :: surfWriteRVx,   surfWriteRVy,      surfWriteRVz
       logical :: surfWriteCp,    surfWritePtotLoss, surfWriteMach
       logical :: surfWriteRMach
       logical :: surfWriteCf,    surfWriteCh,       surfWriteYPlus
       logical :: surfWriteCfx,   surfWriteCfy,      surfWriteCfz
       logical :: surfWriteBlank, surfWriteLift,     surfWriteSepSensor
       logical :: surfWriteCavitation, surfWriteDrag
       logical :: surfWriteGC

!
!      ******************************************************************
!      *                                                                *
!      * The logical variables, which define the extra volume variables *
!      * to be written.                                                 *
!      *                                                                *
!      ******************************************************************
!
       logical :: volWriteMx,           volWriteMy,       volWriteMz
       logical :: volWriteRVx,          volWriteRVy,      volWriteRVz
       logical :: volWriteRhoE,         volWriteTemp,     volWriteCp
       logical :: volWriteMach,         volWriteMachTurb, volWriteEddyVis
       logical :: volWriteRMach
       logical :: volWriteRatioEddyVis, volWriteDist,     volWriteVortx
       logical :: volWritevorty,        volWritevortz,    volWriteVort
       logical :: volWritePtotLoss,     volWriteResRho,   volWriteresMom
       logical :: volWriteResRhoE,      volWriteResTurb,  volWriteBlank
       logical :: volWriteShock,        volWriteFilteredShock
       logical :: volWriteGC
!
!      ******************************************************************
!      *                                                                *
!      * The logical variables, which define the isosurface variables   *
!      * to be written.                                                 *
!      *                                                                *
!      ******************************************************************
!
       logical :: isoWriteRho,          isoWriteVx,       isoWriteVy
       logical :: isoWriteVz,           isoWriteP,        isoWriteTurb
       logical :: isoWriteMx,           isoWriteMy,       isoWriteMz
       logical :: isoWriteRVx,          isoWriteRVy,      isoWriteRVz
       logical :: isoWriteRhoE,         isoWriteTemp,     isoWriteCp
       logical :: isoWriteMach,         isoWriteMachTurb, isoWriteEddyVis
       logical :: isoWriteRMach
       logical :: isoWriteRatioEddyVis, isoWriteDist,     isoWriteVortx
       logical :: isoWritevorty,        isoWritevortz,    isoWriteVort
       logical :: isoWritePtotLoss,     isoWriteResRho,   isoWriteresMom
       logical :: isoWriteResRhoE,      isoWriteResTurb,  isoWriteBlank
       logical :: isoWriteShock,        isoWriteFilteredShock
!
!      ******************************************************************
!      *                                                                *
!      * Extra variables defining the type and number of iso surfaces   *
!      * to be written.                                                 *
!      *                                                                *
!      ******************************************************************
!
       integer(kind=intType) :: nIsoSurface = 0
       real(kind=realType), dimension(:), allocatable :: isoValues
       character(len=maxCGNSNameLen), dimension(:), allocatable :: isoSurfaceNames


       end module extraOutput
