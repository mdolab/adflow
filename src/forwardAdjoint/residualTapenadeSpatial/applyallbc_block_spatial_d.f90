!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
!
!  Differentiation of applyallbc_block in forward (tangent) mode:
!   variations   of useful results: *p *gamma *w
!   with respect to varying inputs: *x *si *sj *sk
!
!      ******************************************************************
!      *                                                                *
!      * File:          applyAllBC.f90                                  *
!      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
!      * Starting date: 03-07-2003                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
SUBROUTINE APPLYALLBC_BLOCK_SPATIAL_D(secondhalo)
  USE ITERATION_SPATIAL_D
  USE INPUTDISCRETIZATION_SPATIAL_D
  USE INPUTTIMESPECTRAL_SPATIAL_D
  USE BLOCKPOINTERS_SPATIAL_D
  USE FLOWVARREFSTATE_SPATIAL_D
  IMPLICIT NONE
!   ! Domain-interface boundary conditions,
!   ! when coupled with other solvers.
!   call bcDomainInterface(secondHalo, correctForK)
!   ! Supersonic inflow boundary conditions.
!   call bcSupersonicInflow(secondHalo, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * applyAllBC applies all boundary conditions for the all         *
!      * blocks on the grid level currentLevel.                         *
!      *                                                                *
!      ******************************************************************
!
!
!      Subroutine arguments.
!
  LOGICAL, INTENT(IN) :: secondhalo
!
!      Local variables.
!
  INTEGER(kind=inttype) :: nn, sps
  LOGICAL :: correctfork
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
! Determine whether or not the total energy must be corrected
! for the presence of the turbulent kinetic energy.
  IF (kpresent) THEN
    IF (currentlevel .LE. groundlevel .OR. turbcoupled) THEN
      correctfork = .true.
    ELSE
      correctfork = .false.
    END IF
  ELSE
    correctfork = .false.
  END IF
! Apply all the boundary conditions. The order is important.
! The symmetry boundary conditions.
  CALL BCSYMM_CD(secondhalo)
  CALL BCSYMMPOLAR_SPATIAL_D(secondhalo)
!call bcEulerWall(secondHalo, correctForK)
! The viscous wall boundary conditions.
!call bcNSWallAdiabatic( secondHalo, correctForK)
!call bcNSWallIsothermal(secondHalo, correctForK)
! The farfield is a special case, because the treatment
! differs when preconditioning is used. Make that distinction
! and call the appropriate routine.
  SELECT CASE  (precond) 
  CASE (noprecond) 
    CALL BCFARFIELD_SPATIAL_D(secondhalo, correctfork)
  CASE (turkel) 
    CALL TERMINATE_CD('applyAllBC', &
&                'Farfield boundary conditions for Turkel ', &
&                'preconditioner not implemented')
    pd = 0.0
    gammad = 0.0
  CASE (choimerkle) 
    CALL TERMINATE_CD('applyAllBC', &
&                'Farfield boundary conditions for Choi and ', &
&                'Merkle preconditioner not implemented')
    pd = 0.0
    gammad = 0.0
  CASE DEFAULT
    pd = 0.0
    gammad = 0.0
  END SELECT
!   ! Subsonic outflow and bleed outflow boundaries.
!   call bcSubsonicOutflow(secondHalo, correctForK)
!   ! Subsonic inflow boundary.
!   call bcSubsonicInflow(secondHalo, correctForK)
!   ! Bleed inflow regions.
!   call bcBleedInflow( secondHalo, correctForK)
!   ! Engine boundary conditions. Not implemented yet.
!   call bcMdot(secondHalo, correctForK)
!   call bcThrust(secondHalo, correctForK)
!   ! Extrapolation boundary conditions; this also includes
!   ! the supersonic outflow boundary conditions. The difference
!   ! between the two is that the extrap boundary conditions
!   ! correspond to singular lines and supersonic outflow
!   ! boundaries to physical boundaries. The treatment however
!   ! is identical.
!   call bcExtrap(secondHalo, correctForK)
! Inviscid wall boundary conditions.
  CALL BCEULERWALL_SPATIAL_D(secondhalo, correctfork)
END SUBROUTINE APPLYALLBC_BLOCK_SPATIAL_D
