   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of bcturbfarfield in reverse (adjoint) mode (with options i4 dr8 r8 noISIZE):
   !   gradient     of useful results: winf *bvtj1 *bvtj2 *bvtk1 *bvtk2
   !                *bvti1 *bvti2
   !   with respect to varying inputs: winf *bvtj1 *bvtj2 *bvtk1 *bvtk2
   !                *bvti1 *bvti2
   !   Plus diff mem management of: bvtj1:in bvtj2:in bvtk1:in bvtk2:in
   !                bvti1:in bvti2:in bcdata:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          bcTurbFarfield.f90                              *
   !      * Author:        Georgi Kalitzin, Edwin van der Weide            *
   !      * Starting date: 06-15-2003                                      *
   !      * Last modified: 06-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE BCTURBFARFIELD_B(nn)
   !
   !      ******************************************************************
   !      *                                                                *
   !      * bcTurbFarfield applies the implicit treatment of the           *
   !      * farfield boundary condition to subface nn. As the farfield     *
   !      * boundary condition is independent of the turbulence model,     *
   !      * this routine is valid for all models. It is assumed that the   *
   !      * pointers in blockPointers are already set to the correct       *
   !      * block on the correct grid level.                               *
   !      *                                                                *
   !      ******************************************************************
   !
   USE BLOCKPOINTERS
   USE BCTYPES
   USE CONSTANTS
   USE FLOWVARREFSTATE
   IMPLICIT NONE
   !
   !      Subroutine arguments.
   !
   INTEGER(kind=inttype), INTENT(IN) :: nn
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i, j, l
   REAL(kind=realtype) :: nnx, nny, nnz, dot
   INTEGER :: branch
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Loop over the faces of the subfaces and set the values of
   ! bmt and bvt for an implicit treatment.
   DO j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
   DO i=bcdata(nn)%icbeg,bcdata(nn)%icend
   ! Determine the dot product between the outward pointing
   ! normal and the free stream velocity direction and add the
   ! possible grid velocity.
   dot = bcdata(nn)%norm(i, j, 1)*winf(ivx) + bcdata(nn)%norm(i, j, 2&
   &       )*winf(ivy) + bcdata(nn)%norm(i, j, 3)*winf(ivz) - bcdata(nn)%&
   &       rface(i, j)
   ! Determine whether we are dealing with an inflow or
   ! outflow boundary here.
   IF (dot .GT. zero) THEN
   CALL PUSHCONTROL1B(1)
   ELSE
   ! Inflow. Turbulent variables are prescribed.
   DO l=nt1,nt2
   SELECT CASE  (bcfaceid(nn)) 
   CASE (imin) 
   CALL PUSHCONTROL3B(5)
   CASE (imax) 
   CALL PUSHCONTROL3B(4)
   CASE (jmin) 
   CALL PUSHCONTROL3B(3)
   CASE (jmax) 
   CALL PUSHCONTROL3B(2)
   CASE (kmin) 
   CALL PUSHCONTROL3B(1)
   CASE (kmax) 
   CALL PUSHCONTROL3B(0)
   CASE DEFAULT
   CALL PUSHCONTROL3B(6)
   END SELECT
   END DO
   CALL PUSHCONTROL1B(0)
   END IF
   END DO
   END DO
   DO j=bcdata(nn)%jcend,bcdata(nn)%jcbeg,-1
   DO i=bcdata(nn)%icend,bcdata(nn)%icbeg,-1
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   DO l=nt2,nt1,-1
   CALL POPCONTROL3B(branch)
   IF (branch .LT. 3) THEN
   IF (branch .EQ. 0) THEN
   winfd(l) = winfd(l) + bvtk2d(i, j, l)
   bvtk2d(i, j, l) = 0.0_8
   ELSE IF (branch .EQ. 1) THEN
   winfd(l) = winfd(l) + bvtk1d(i, j, l)
   bvtk1d(i, j, l) = 0.0_8
   ELSE
   winfd(l) = winfd(l) + bvtj2d(i, j, l)
   bvtj2d(i, j, l) = 0.0_8
   END IF
   ELSE IF (branch .LT. 5) THEN
   IF (branch .EQ. 3) THEN
   winfd(l) = winfd(l) + bvtj1d(i, j, l)
   bvtj1d(i, j, l) = 0.0_8
   ELSE
   winfd(l) = winfd(l) + bvti2d(i, j, l)
   bvti2d(i, j, l) = 0.0_8
   END IF
   ELSE IF (branch .EQ. 5) THEN
   winfd(l) = winfd(l) + bvti1d(i, j, l)
   bvti1d(i, j, l) = 0.0_8
   END IF
   END DO
   END IF
   END DO
   END DO
   END SUBROUTINE BCTURBFARFIELD_B
