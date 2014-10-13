   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of turbbcnswall in reverse (adjoint) mode (with options i4 dr8 r8 noISIZE):
   !   gradient     of useful results: *rev *bvtj1 *bvtj2 *bmtk1 *w
   !                *bmtk2 *rlv *bvtk1 *bvtk2 *bmti1 *bmti2 *bvti1
   !                *bvti2 *bmtj1 *bmtj2
   !   with respect to varying inputs: *rev *bvtj1 *bvtj2 *bmtk1 *w
   !                *bmtk2 *rlv *bvtk1 *bvtk2 *bmti1 *bmti2 *bvti1
   !                *bvti2 *bmtj1 *bmtj2
   !   Plus diff mem management of: rev:in bvtj1:in bvtj2:in bmtk1:in
   !                w:in bmtk2:in rlv:in bvtk1:in bvtk2:in bmti1:in
   !                bmti2:in bvti1:in bvti2:in bmtj1:in bmtj2:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          turbBCNSWall.f90                                *
   !      * Author:        Edwin van der Weide                             *
   !      * Starting date: 05-30-2003                                      *
   !      * Last modified: 06-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE TURBBCNSWALL_B(secondhalo)
   !
   !      ******************************************************************
   !      *                                                                *
   !      * turbBCNSWall applies the viscous wall boundary conditions      *
   !      * of the turbulent transport equations to a block. It is assumed *
   !      * that the pointers in blockPointers are already set to the      *
   !      * correct block on the correct grid level.                       *
   !      *                                                                *
   !      ******************************************************************
   !
   USE BCTYPES
   USE BLOCKPOINTERS_B
   USE FLOWVARREFSTATE
   IMPLICIT NONE
   !
   !      Subroutine argument.
   !
   LOGICAL, INTENT(IN) :: secondhalo
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: nn, i, j, l, m
   REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmt
   REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmtb
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: bvt
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: bvtb
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww0, ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww0b, ww1b, ww2b
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev0, rev1, rev2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev0b, rev1b, rev2b
   REAL(kind=realtype) :: tmp
   REAL(kind=realtype) :: tmp0
   REAL(kind=realtype) :: tmp1
   REAL(kind=realtype) :: tmp2
   INTEGER :: ad_from
   INTEGER :: ad_to
   INTEGER :: ad_from0
   INTEGER :: ad_to0
   INTEGER :: ad_from1
   INTEGER :: ad_to1
   INTEGER :: ad_from2
   INTEGER :: ad_to2
   INTEGER :: ad_from3
   INTEGER :: ad_to3
   INTEGER :: ad_from4
   INTEGER :: ad_to4
   INTEGER :: ad_from5
   INTEGER :: ad_to5
   INTEGER :: ad_from6
   INTEGER :: ad_to6
   INTEGER :: branch
   INTERFACE 
   SUBROUTINE PUSHPOINTER4(pp)
   REAL, POINTER :: pp
   END SUBROUTINE PUSHPOINTER4
   SUBROUTINE LOOKPOINTER4(pp)
   REAL, POINTER :: pp
   END SUBROUTINE LOOKPOINTER4
   SUBROUTINE POPPOINTER4(pp)
   REAL, POINTER :: pp
   END SUBROUTINE POPPOINTER4
   END INTERFACE
      REAL(kind=realtype) :: tmpb
   REAL(kind=realtype) :: tmpb2
   REAL(kind=realtype) :: tmpb1
   REAL(kind=realtype) :: tmpb0
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Loop over the viscous subfaces of this block.
   bocos:DO nn=1,nviscbocos
   ! Set the corresponding arrays.
   CALL PUSHREAL8ARRAY(bmtj2, SIZE(bmtj2, 1)*SIZE(bmtj2, 2)*SIZE(bmtj2&
   &                 , 3)*SIZE(bmtj2, 4))
   CALL PUSHREAL8ARRAY(bmtj1, SIZE(bmtj1, 1)*SIZE(bmtj1, 2)*SIZE(bmtj1&
   &                 , 3)*SIZE(bmtj1, 4))
   CALL PUSHREAL8ARRAY(bmti2, SIZE(bmti2, 1)*SIZE(bmti2, 2)*SIZE(bmti2&
   &                 , 3)*SIZE(bmti2, 4))
   CALL PUSHREAL8ARRAY(bmti1, SIZE(bmti1, 1)*SIZE(bmti1, 2)*SIZE(bmti1&
   &                 , 3)*SIZE(bmti1, 4))
   CALL PUSHREAL8ARRAY(bmtk2, SIZE(bmtk2, 1)*SIZE(bmtk2, 2)*SIZE(bmtk2&
   &                 , 3)*SIZE(bmtk2, 4))
   CALL PUSHREAL8ARRAY(bmtk1, SIZE(bmtk1, 1)*SIZE(bmtk1, 2)*SIZE(bmtk1&
   &                 , 3)*SIZE(bmtk1, 4))
   CALL BCTURBWALL(nn)
   ! Determine the block face on which this subface is located
   ! and set some pointers accordingly.
   CALL PUSHCONTROL4B(12)
   CALL PUSHPOINTER4(bmtb)
   bmtb => bmti1b
   CALL PUSHPOINTER4(bmt)
   bmt => bmti1
   CALL PUSHPOINTER4(bvtb)
   bvtb => bvti1b
   CALL PUSHPOINTER4(bvt)
   bvt => bvti1
   CALL PUSHPOINTER4(ww0b)
   ww0b => wb(0, 1:, 1:, :)
   CALL PUSHPOINTER4(ww1b)
   ww1b => wb(1, 1:, 1:, :)
   CALL PUSHPOINTER4(ww1)
   ww1 => w(1, 1:, 1:, :)
   CALL PUSHPOINTER4(ww2b)
   ww2b => wb(2, 1:, 1:, :)
   CALL PUSHPOINTER4(ww2)
   ww2 => w(2, 1:, 1:, :)
   IF (eddymodel) THEN
   CALL PUSHPOINTER4(rev0b)
   rev0b => revb(0, 1:, 1:)
   CALL PUSHPOINTER4(rev1b)
   rev1b => revb(1, 1:, 1:)
   CALL PUSHPOINTER4(rev1)
   rev1 => rev(1, 1:, 1:)
   CALL PUSHPOINTER4(rev2b)
   rev2b => revb(2, 1:, 1:)
   CALL PUSHPOINTER4(rev2)
   rev2 => rev(2, 1:, 1:)
   CALL PUSHCONTROL4B(11)
   ELSE
   CALL PUSHCONTROL4B(10)
   END IF
   ad_from0 = bcdata(nn)%jcbeg
   ! Loop over the faces and set the state in
   ! the turbulent halo cells.
   DO j=ad_from0,bcdata(nn)%jcend
   ad_from = bcdata(nn)%icbeg
   DO i=ad_from,bcdata(nn)%icend
   DO l=nt1,nt2
   CALL PUSHREAL8(ww1(i, j, l))
   ww1(i, j, l) = bvt(i, j, l)
   DO m=nt1,nt2
   tmp = ww1(i, j, l) - bmt(i, j, l, m)*ww2(i, j, m)
   CALL PUSHREAL8(ww1(i, j, l))
   ww1(i, j, l) = tmp
   END DO
   END DO
   END DO
   CALL PUSHINTEGER4(i - 1)
   CALL PUSHINTEGER4(ad_from)
   END DO
   CALL PUSHINTEGER4(j - 1)
   CALL PUSHINTEGER4(ad_from0)
   ! Use constant extrapolation if the state in the second halo
   ! must be computed.
   IF (secondhalo) THEN
   ad_from2 = bcdata(nn)%jcbeg
   DO j=ad_from2,bcdata(nn)%jcend
   ad_from1 = bcdata(nn)%icbeg
   DO i=ad_from1,bcdata(nn)%icend
   DO l=nt1,nt2
   tmp0 = ww1(i, j, l)
   CALL PUSHREAL8(ww0(i, j, l))
   ww0(i, j, l) = tmp0
   END DO
   END DO
   CALL PUSHINTEGER4(i - 1)
   CALL PUSHINTEGER4(ad_from1)
   END DO
   CALL PUSHINTEGER4(j - 1)
   CALL PUSHINTEGER4(ad_from2)
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHCONTROL1B(1)
   END IF
   ! Set the eddy viscosity for an eddy model.
   IF (eddymodel) THEN
   ad_from4 = bcdata(nn)%jcbeg
   DO j=ad_from4,bcdata(nn)%jcend
   ad_from3 = bcdata(nn)%icbeg
   DO i=ad_from3,bcdata(nn)%icend
   tmp1 = -rev2(i, j)
   rev1(i, j) = tmp1
   END DO
   CALL PUSHINTEGER4(i - 1)
   CALL PUSHINTEGER4(ad_from3)
   END DO
   CALL PUSHINTEGER4(j - 1)
   CALL PUSHINTEGER4(ad_from4)
   IF (secondhalo) THEN
   ad_from6 = bcdata(nn)%jcbeg
   DO j=ad_from6,bcdata(nn)%jcend
   ad_from5 = bcdata(nn)%icbeg
   DO i=ad_from5,bcdata(nn)%icend
   tmp2 = rev1(i, j)
   rev0(i, j) = tmp2
   END DO
   CALL PUSHINTEGER4(i - 1)
   CALL PUSHINTEGER4(ad_from5)
   END DO
   CALL PUSHINTEGER4(j - 1)
   CALL PUSHINTEGER4(ad_from6)
   CALL PUSHCONTROL2B(2)
   ELSE
   CALL PUSHCONTROL2B(1)
   END IF
   ELSE
   CALL PUSHCONTROL2B(0)
   END IF
   END DO bocos
   DO nn=nviscbocos,1,-1
   CALL POPCONTROL2B(branch)
   IF (branch .NE. 0) THEN
   IF (branch .NE. 1) THEN
   CALL POPINTEGER4(ad_from6)
   CALL POPINTEGER4(ad_to6)
   DO j=ad_to6,ad_from6,-1
   CALL POPINTEGER4(ad_from5)
   CALL POPINTEGER4(ad_to5)
   DO i=ad_to5,ad_from5,-1
   tmpb2 = rev0b(i, j)
   rev0b(i, j) = 0.0_8
   rev1b(i, j) = rev1b(i, j) + tmpb2
   END DO
   END DO
   END IF
   CALL POPINTEGER4(ad_from4)
   CALL POPINTEGER4(ad_to4)
   DO j=ad_to4,ad_from4,-1
   CALL POPINTEGER4(ad_from3)
   CALL POPINTEGER4(ad_to3)
   DO i=ad_to3,ad_from3,-1
   tmpb1 = rev1b(i, j)
   rev1b(i, j) = 0.0_8
   rev2b(i, j) = rev2b(i, j) - tmpb1
   END DO
   END DO
   END IF
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPINTEGER4(ad_from2)
   CALL POPINTEGER4(ad_to2)
   DO j=ad_to2,ad_from2,-1
   CALL POPINTEGER4(ad_from1)
   CALL POPINTEGER4(ad_to1)
   DO i=ad_to1,ad_from1,-1
   DO l=nt2,nt1,-1
   CALL POPREAL8(ww0(i, j, l))
   tmpb0 = ww0b(i, j, l)
   ww0b(i, j, l) = 0.0_8
   ww1b(i, j, l) = ww1b(i, j, l) + tmpb0
   END DO
   END DO
   END DO
   END IF
   CALL POPINTEGER4(ad_from0)
   CALL POPINTEGER4(ad_to0)
   DO j=ad_to0,ad_from0,-1
   CALL POPINTEGER4(ad_from)
   CALL POPINTEGER4(ad_to)
   DO i=ad_to,ad_from,-1
   DO l=nt2,nt1,-1
   DO m=nt2,nt1,-1
   CALL POPREAL8(ww1(i, j, l))
   tmpb = ww1b(i, j, l)
   ww1b(i, j, l) = tmpb
   bmtb(i, j, l, m) = bmtb(i, j, l, m) - ww2(i, j, m)*tmpb
   ww2b(i, j, m) = ww2b(i, j, m) - bmt(i, j, l, m)*tmpb
   END DO
   CALL POPREAL8(ww1(i, j, l))
   bvtb(i, j, l) = bvtb(i, j, l) + ww1b(i, j, l)
   ww1b(i, j, l) = 0.0_8
   END DO
   END DO
   END DO
   CALL POPCONTROL4B(branch)
   IF (branch .LT. 6) THEN
   IF (branch .LT. 3) THEN
   IF (branch .NE. 0) THEN
   IF (branch .EQ. 1) THEN
   CALL POPPOINTER4(rev2)
   CALL POPPOINTER4(rev2b)
   CALL POPPOINTER4(rev1)
   CALL POPPOINTER4(rev1b)
   CALL POPPOINTER4(rev0b)
   ELSE
   GOTO 100
   END IF
   END IF
   CALL POPPOINTER4(ww2)
   CALL POPPOINTER4(ww2b)
   CALL POPPOINTER4(ww1)
   CALL POPPOINTER4(ww1b)
   CALL POPPOINTER4(ww0b)
   CALL POPPOINTER4(bvt)
   bvtk2b => bvtb
   CALL POPPOINTER4(bvtb)
   CALL POPPOINTER4(bmt)
   bmtk2b => bmtb
   CALL POPPOINTER4(bmtb)
   GOTO 120
   ELSE IF (branch .EQ. 3) THEN
   CALL POPPOINTER4(rev2)
   CALL POPPOINTER4(rev2b)
   CALL POPPOINTER4(rev1)
   CALL POPPOINTER4(rev1b)
   CALL POPPOINTER4(rev0b)
   ELSE
   IF (branch .NE. 4) THEN
   CALL POPPOINTER4(rev2)
   CALL POPPOINTER4(rev2b)
   CALL POPPOINTER4(rev1)
   CALL POPPOINTER4(rev1b)
   CALL POPPOINTER4(rev0b)
   END IF
   CALL POPPOINTER4(ww2)
   CALL POPPOINTER4(ww2b)
   CALL POPPOINTER4(ww1)
   CALL POPPOINTER4(ww1b)
   CALL POPPOINTER4(ww0b)
   CALL POPPOINTER4(bvt)
   bvtj2b => bvtb
   CALL POPPOINTER4(bvtb)
   CALL POPPOINTER4(bmt)
   bmtj2b => bmtb
   CALL POPPOINTER4(bmtb)
   GOTO 120
   END IF
   100  CALL POPPOINTER4(ww2)
   CALL POPPOINTER4(ww2b)
   CALL POPPOINTER4(ww1)
   CALL POPPOINTER4(ww1b)
   CALL POPPOINTER4(ww0b)
   CALL POPPOINTER4(bvt)
   bvtk1b => bvtb
   CALL POPPOINTER4(bvtb)
   CALL POPPOINTER4(bmt)
   bmtk1b => bmtb
   CALL POPPOINTER4(bmtb)
   ELSE
   IF (branch .LT. 9) THEN
   IF (branch .NE. 6) THEN
   IF (branch .EQ. 7) THEN
   CALL POPPOINTER4(rev2)
   CALL POPPOINTER4(rev2b)
   CALL POPPOINTER4(rev1)
   CALL POPPOINTER4(rev1b)
   CALL POPPOINTER4(rev0b)
   ELSE
   GOTO 110
   END IF
   END IF
   CALL POPPOINTER4(ww2)
   CALL POPPOINTER4(ww2b)
   CALL POPPOINTER4(ww1)
   CALL POPPOINTER4(ww1b)
   CALL POPPOINTER4(ww0b)
   CALL POPPOINTER4(bvt)
   bvtj1b => bvtb
   CALL POPPOINTER4(bvtb)
   CALL POPPOINTER4(bmt)
   bmtj1b => bmtb
   CALL POPPOINTER4(bmtb)
   ELSE
   IF (branch .LT. 11) THEN
   IF (branch .EQ. 9) THEN
   CALL POPPOINTER4(rev2)
   CALL POPPOINTER4(rev2b)
   CALL POPPOINTER4(rev1)
   CALL POPPOINTER4(rev1b)
   CALL POPPOINTER4(rev0b)
   GOTO 110
   END IF
   ELSE IF (branch .EQ. 11) THEN
   CALL POPPOINTER4(rev2)
   CALL POPPOINTER4(rev2b)
   CALL POPPOINTER4(rev1)
   CALL POPPOINTER4(rev1b)
   CALL POPPOINTER4(rev0b)
   ELSE
   GOTO 120
   END IF
   CALL POPPOINTER4(ww2)
   CALL POPPOINTER4(ww2b)
   CALL POPPOINTER4(ww1)
   CALL POPPOINTER4(ww1b)
   CALL POPPOINTER4(ww0b)
   CALL POPPOINTER4(bvt)
   bvti1b => bvtb
   CALL POPPOINTER4(bvtb)
   CALL POPPOINTER4(bmt)
   bmti1b => bmtb
   CALL POPPOINTER4(bmtb)
   END IF
   GOTO 120
   110  CALL POPPOINTER4(ww2)
   CALL POPPOINTER4(ww2b)
   CALL POPPOINTER4(ww1)
   CALL POPPOINTER4(ww1b)
   CALL POPPOINTER4(ww0b)
   CALL POPPOINTER4(bvt)
   bvti2b => bvtb
   CALL POPPOINTER4(bvtb)
   CALL POPPOINTER4(bmt)
   bmti2b => bmtb
   CALL POPPOINTER4(bmtb)
   END IF
   120 CALL POPREAL8ARRAY(bmtk1, SIZE(bmtk1, 1)*SIZE(bmtk1, 2)*SIZE(bmtk1&
   &                 , 3)*SIZE(bmtk1, 4))
   CALL POPREAL8ARRAY(bmtk2, SIZE(bmtk2, 1)*SIZE(bmtk2, 2)*SIZE(bmtk2, &
   &                3)*SIZE(bmtk2, 4))
   CALL POPREAL8ARRAY(bmti1, SIZE(bmti1, 1)*SIZE(bmti1, 2)*SIZE(bmti1, &
   &                3)*SIZE(bmti1, 4))
   CALL POPREAL8ARRAY(bmti2, SIZE(bmti2, 1)*SIZE(bmti2, 2)*SIZE(bmti2, &
   &                3)*SIZE(bmti2, 4))
   CALL POPREAL8ARRAY(bmtj1, SIZE(bmtj1, 1)*SIZE(bmtj1, 2)*SIZE(bmtj1, &
   &                3)*SIZE(bmtj1, 4))
   CALL POPREAL8ARRAY(bmtj2, SIZE(bmtj2, 1)*SIZE(bmtj2, 2)*SIZE(bmtj2, &
   &                3)*SIZE(bmtj2, 4))
   CALL BCTURBWALL_B(nn)
   END DO
   END SUBROUTINE TURBBCNSWALL_B
