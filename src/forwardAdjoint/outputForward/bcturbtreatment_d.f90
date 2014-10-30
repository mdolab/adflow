   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of bcturbtreatment in forward (tangent) mode (with options i4 dr8 r8):
   !   variations   of useful results: *bvtj1 *bvtj2 *bmtk1 *bmtk2
   !                *bvtk1 *bvtk2 *bmti1 *bmti2 *bvti1 *bvti2 *bmtj1
   !                *bmtj2
   !   with respect to varying inputs: *bvtj1 *bvtj2 *bmtk1 *w *bmtk2
   !                *rlv *bvtk1 *bvtk2 *bmti1 *bmti2 *bvti1 *bvti2
   !                *bmtj1 *bmtj2 winf
   !   Plus diff mem management of: bvtj1:in bvtj2:in bmtk1:in w:in
   !                bmtk2:in rlv:in bvtk1:in bvtk2:in bmti1:in bmti2:in
   !                bvti1:in bvti2:in bmtj1:in bmtj2:in bcdata:in
   !                *bcdata.norm:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          bcTurbTreatment.f90                             *
   !      * Author:        Georgi Kalitzin, Edwin van der Weide            *
   !      *                Seonghyeon Hahn                                 *
   !      * Starting date: 06-13-2003                                      *
   !      * Last modified: 08-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE BCTURBTREATMENT_D()
   !
   !      ******************************************************************
   !      *                                                                *
   !      * bcTurbTreatment sets the arrays bmti1, bvti1, etc, such that   *
   !      * the physical boundary conditions are treated correctly.        *
   !      * It is assumed that the variables in blockPointers already      *
   !      * point to the correct block.                                    *
   !      *                                                                *
   !      * The turbulent variable in the halo is computed as follows:     *
   !      * wHalo = -bmt*wInternal + bvt for every block facer. As it is   *
   !      * possible to have a coupling in the boundary conditions bmt     *
   !      * actually are matrices. If there is no coupling between the     *
   !      * boundary conditions of the turbulence equations bmt is a       *
   !      * diagonal matrix.                                               *
   !      *                                                                *
   !      ******************************************************************
   !
   USE BCTYPES
   USE BLOCKPOINTERS_D
   USE FLOWVARREFSTATE
   IMPLICIT NONE
   !
   !      Local variable.
   !
   INTEGER(kind=inttype) :: nn, i, j, k, l, m
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Initialize the arrays for the boundary condition treatment
   ! to zero, such that internal block boundaries are solved
   ! correctly (i.e. explicitly).
   DO k=1,ke
   DO j=1,je
   DO l=nt1,nt2
   DO m=nt1,nt2
   bmti1d(j, k, l, m) = 0.0_8
   bmti1(j, k, l, m) = zero
   bmti2d(j, k, l, m) = 0.0_8
   bmti2(j, k, l, m) = zero
   END DO
   bvti1d(j, k, l) = 0.0_8
   bvti1(j, k, l) = zero
   bvti2d(j, k, l) = 0.0_8
   bvti2(j, k, l) = zero
   END DO
   END DO
   END DO
   DO k=1,ke
   DO i=1,ie
   DO l=nt1,nt2
   DO m=nt1,nt2
   bmtj1d(i, k, l, m) = 0.0_8
   bmtj1(i, k, l, m) = zero
   bmtj2d(i, k, l, m) = 0.0_8
   bmtj2(i, k, l, m) = zero
   END DO
   bvtj1d(i, k, l) = 0.0_8
   bvtj1(i, k, l) = zero
   bvtj2d(i, k, l) = 0.0_8
   bvtj2(i, k, l) = zero
   END DO
   END DO
   END DO
   DO j=1,je
   DO i=1,ie
   DO l=nt1,nt2
   DO m=nt1,nt2
   bmtk1d(i, j, l, m) = 0.0_8
   bmtk1(i, j, l, m) = zero
   bmtk2d(i, j, l, m) = 0.0_8
   bmtk2(i, j, l, m) = zero
   END DO
   bvtk1d(i, j, l) = 0.0_8
   bvtk1(i, j, l) = zero
   bvtk2d(i, j, l) = 0.0_8
   bvtk2(i, j, l) = zero
   END DO
   END DO
   END DO
   ! Loop over the boundary condition subfaces of this block.
   bocos:DO nn=1,nbocos
   ! Determine the kind of boundary condition for this subface.
   SELECT CASE  (bctype(nn)) 
   CASE (nswalladiabatic, nswallisothermal) 
   ! Viscous wall. There is no difference between an adiabatic
   ! and an isothermal wall for the turbulent equations.
   ! Set the implicit treatment of the wall boundary conditions.
   CALL BCTURBWALL_D(nn)
   CASE (symm, symmpolar, eulerwall) 
   !=============================================================
   !=============================================================
   ! Symmetry, polar symmetry or inviscid wall. Treatment of
   ! the turbulent equations is identical.
   CALL BCTURBSYMM_D(nn)
   CASE (farfield) 
   !=============================================================
   ! Farfield. The kind of boundary condition to be applied,
   ! inflow or outflow, depends on the local conditions.
   CALL BCTURBFARFIELD_D(nn)
   CASE DEFAULT
   !=============================================================
   !=============================================================
   CALL TERMINATE('bcTurbTreatment', 'Unknown boundary condition')
   END SELECT
   END DO bocos
   END SUBROUTINE BCTURBTREATMENT_D
