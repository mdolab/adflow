!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
!
!
!      ******************************************************************
!      *                                                                *
!      * File:          bcTurbSymm.F90                                  *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 06-11-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
SUBROUTINE BCTURBSYMM_CD(nn)
  USE BCTYPES_SPATIAL_D
  USE BLOCKPOINTERS_SPATIAL_D
  USE FLOWVARREFSTATE_SPATIAL_D
  IMPLICIT NONE
!
!      ******************************************************************
!      *                                                                *
!      * bcTurbSymm applies the implicit treatment of the symmetry      *
!      * boundary condition (or inviscid wall) to subface nn. As the    *
!      * symmetry boundary condition is independent of the turbulence   *
!      * model, this routine is valid for all models. It is assumed     *
!      * that the pointers in blockPointers are already set to the      *
!      * correct block on the correct grid level.                       *
!      *                                                                *
!      ******************************************************************
!
!
!      Subroutine arguments.
!
  INTEGER(kind=inttype), INTENT(IN) :: nn
!
!      Local variables.
!
  INTEGER(kind=inttype) :: i, j, l
  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmt
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
! Set the pointer for bmt, depending on the block face on which
! the subface is located.
  SELECT CASE  (bcfaceid(nn)) 
  CASE (imin) 
    bmt => bmti1
  CASE (imax) 
    bmt => bmti2
  CASE (jmin) 
    bmt => bmtj1
  CASE (jmax) 
    bmt => bmtj2
  CASE (kmin) 
    bmt => bmtk1
  CASE (kmax) 
    bmt => bmtk2
  END SELECT
! Loop over the faces of the subfaces and set the values of bmt
! for an implicit treatment. For a symmetry face this means
! that the halo value is set to the internal value.
  DO j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
    DO i=bcdata(nn)%icbeg,bcdata(nn)%icend
      DO l=nt1,nt2
        bmt(i, j, l, l) = -one
      END DO
    END DO
  END DO
END SUBROUTINE BCTURBSYMM_CD
