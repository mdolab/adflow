   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of bcturbfarfield in forward (tangent) mode (with options i4 dr8 r8):
   !   variations   of useful results: *bvtj1 *bvtj2 *bmtk1 *bmtk2
   !                *bvtk1 *bvtk2 *bmti1 *bmti2 *bvti1 *bvti2 *bmtj1
   !                *bmtj2
   !   with respect to varying inputs: *bvtj1 *bvtj2 *bmtk1 *bmtk2
   !                *bvtk1 *bvtk2 *bmti1 *bmti2 *bvti1 *bvti2 *bmtj1
   !                *bmtj2
   !   Plus diff mem management of: bvtj1:in bvtj2:in bmtk1:in bmtk2:in
   !                bvtk1:in bvtk2:in bmti1:in bmti2:in bvti1:in bvti2:in
   !                bmtj1:in bmtj2:in bcdata:in *bcdata.norm:in *bcdata.rface:in
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
   SUBROUTINE BCTURBFARFIELD_D(nn)
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
   USE BLOCKPOINTERS_D
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
   REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmt
   REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmtd
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: bvt
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: bvtd
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Set the pointers for bmt and bvt, depending on the block face
   ! on which the subface is located.
   SELECT CASE  (bcfaceid(nn)) 
   CASE (imin) 
   bmtd => bmti1d
   bmt => bmti1
   bvtd => bvti1d
   bvt => bvti1
   CASE (imax) 
   bmtd => bmti2d
   bmt => bmti2
   bvtd => bvti2d
   bvt => bvti2
   CASE (jmin) 
   bmtd => bmtj1d
   bmt => bmtj1
   bvtd => bvtj1d
   bvt => bvtj1
   CASE (jmax) 
   bmtd => bmtj2d
   bmt => bmtj2
   bvtd => bvtj2d
   bvt => bvtj2
   CASE (kmin) 
   bmtd => bmtk1d
   bmt => bmtk1
   bvtd => bvtk1d
   bvt => bvtk1
   CASE (kmax) 
   bmtd => bmtk2d
   bmt => bmtk2
   bvtd => bvtk2d
   bvt => bvtk2
   END SELECT
   ! Loop over the faces of the subfaces and set the values of
   ! bmt and bvt for an implicit treatment.
   DO j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
   DO i=bcdata(nn)%icbeg,bcdata(nn)%icend
   ! Store the three components of the unit normal a bit easier.
   nnx = bcdata(nn)%norm(i, j, 1)
   nny = bcdata(nn)%norm(i, j, 2)
   nnz = bcdata(nn)%norm(i, j, 3)
   ! Determine the dot product between the outward pointing
   ! normal and the free stream velocity direction and add the
   ! possible grid velocity.
   dot = nnx*winf(ivx) + nny*winf(ivy) + nnz*winf(ivz) - bcdata(nn)%&
   &       rface(i, j)
   ! Determine whether we are dealing with an inflow or
   ! outflow boundary here.
   IF (dot .GT. zero) THEN
   ! Outflow. Simply extrapolation or zero Neumann BC
   ! of the turbulent variables.
   DO l=nt1,nt2
   bmtd(i, j, l, l) = 0.0_8
   bmt(i, j, l, l) = -one
   END DO
   ELSE
   ! Inflow. Turbulent variables are prescribed.
   DO l=nt1,nt2
   bvtd(i, j, l) = 0.0_8
   bvt(i, j, l) = winf(l)
   END DO
   END IF
   END DO
   END DO
   END SUBROUTINE BCTURBFARFIELD_D
