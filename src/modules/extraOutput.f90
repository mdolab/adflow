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
       logical :: surfWriteCp,    surfWritePtotLoss, surfWriteMach
       logical :: surfWriteCf,    surfWriteCh,       surfWriteYPlus
       logical :: surfWriteCfx,   surfWriteCfy,      surfWriteCfz
       logical :: surfWriteBlank
!
!      ******************************************************************
!      *                                                                *
!      * The logical variables, which define the extra volume variables *
!      * to be written.                                                 *
!      *                                                                *
!      ******************************************************************
!
       logical :: volWriteMx,           volWriteMy,       volWriteMz
       logical :: volWriteRhoE,         volWriteTemp,     volWriteCp
       logical :: volWriteMach,         volWriteMachTurb, volWriteEddyVis
       logical :: volWriteRatioEddyVis, volWriteDist,     volWriteVortx
       logical :: volWritevorty,        volWritevortz,    volWriteVort
       logical :: volWritePtotLoss,     volWriteResRho,   volWriteresMom
       logical :: volWriteResRhoE,      volWriteResTurb,  volWriteBlank

       end module extraOutput
