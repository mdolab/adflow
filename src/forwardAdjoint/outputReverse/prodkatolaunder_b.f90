   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of prodkatolaunder in reverse (adjoint) mode (with options i4 dr8 r8 noISIZE):
   !   gradient     of useful results: *dw *w
   !   with respect to varying inputs: *dw *w
   !   Plus diff mem management of: dw:in w:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          prodKatoLaunder.f90                             *
   !      * Author:        Georgi Kalitzin, Edwin van der Weide            *
   !      * Starting date: 08-01-2003                                      *
   !      * Last modified: 06-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE PRODKATOLAUNDER_B()
   !
   !      ******************************************************************
   !      *                                                                *
   !      * prodKatoLaunder computes the turbulent production term using   *
   !      * the Kato-Launder formulation.                                  *
   !      *                                                                *
   !      ******************************************************************
   !
   USE BLOCKPOINTERS_B
   USE FLOWVARREFSTATE
   USE SECTION
   USE TURBMOD
   IMPLICIT NONE
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i, j, k
   REAL(kind=realtype) :: ux, uy, uz, vx, vy, vz, wx, wy, wz
   REAL(kind=realtype) :: uxb, uyb, uzb, vxb, vyb, vzb, wxb, wyb, wzb
   REAL(kind=realtype) :: qxx, qyy, qzz, qxy, qxz, qyz, sijsij
   REAL(kind=realtype) :: qxxb, qyyb, qzzb, qxyb, qxzb, qyzb, sijsijb
   REAL(kind=realtype) :: oxy, oxz, oyz, oijoij
   REAL(kind=realtype) :: oxyb, oxzb, oyzb, oijoijb
   REAL(kind=realtype) :: fact, omegax, omegay, omegaz
   INTRINSIC SQRT
   REAL(kind=realtype) :: tempb7
   REAL(kind=realtype) :: tempb6
   REAL(kind=realtype) :: tempb5
   REAL(kind=realtype) :: tempb4
   REAL(kind=realtype) :: tempb3
   REAL(kind=realtype) :: tempb2
   REAL(kind=realtype) :: tempb1
   REAL(kind=realtype) :: tempb0
   REAL(kind=realtype) :: tempb
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Determine the non-dimensional wheel speed of this block.
   ! The vorticity term, which appears in Kato-Launder is of course
   ! not frame invariant. To approximate frame invariance the wheel
   ! speed should be substracted from oxy, oxz and oyz, which results
   ! in the vorticity in the rotating frame. However some people
   ! claim that the absolute vorticity should be used to obtain the
   ! best results. In that omega should be set to zero.
   omegax = timeref*sections(sectionid)%rotrate(1)
   omegay = timeref*sections(sectionid)%rotrate(2)
   omegaz = timeref*sections(sectionid)%rotrate(3)
   ! Loop over the cell centers of the given block. It may be more
   ! efficient to loop over the faces and to scatter the gradient,
   ! but in that case the gradients for u, v and w must be stored.
   ! In the current approach no extra memory is needed.
   DO k=2,kl
   DO j=2,jl
   DO i=2,il
   ! Compute the gradient of u in the cell center. Use is made
   ! of the fact that the surrounding normals sum up to zero,
   ! such that the cell i,j,k does not give a contribution.
   ! The gradient is scaled by a factor 2*vol.
   ux = w(i+1, j, k, ivx)*si(i, j, k, 1) - w(i-1, j, k, ivx)*si(i-1&
   &         , j, k, 1) + w(i, j+1, k, ivx)*sj(i, j, k, 1) - w(i, j-1, k, &
   &         ivx)*sj(i, j-1, k, 1) + w(i, j, k+1, ivx)*sk(i, j, k, 1) - w(i&
   &         , j, k-1, ivx)*sk(i, j, k-1, 1)
   uy = w(i+1, j, k, ivx)*si(i, j, k, 2) - w(i-1, j, k, ivx)*si(i-1&
   &         , j, k, 2) + w(i, j+1, k, ivx)*sj(i, j, k, 2) - w(i, j-1, k, &
   &         ivx)*sj(i, j-1, k, 2) + w(i, j, k+1, ivx)*sk(i, j, k, 2) - w(i&
   &         , j, k-1, ivx)*sk(i, j, k-1, 2)
   uz = w(i+1, j, k, ivx)*si(i, j, k, 3) - w(i-1, j, k, ivx)*si(i-1&
   &         , j, k, 3) + w(i, j+1, k, ivx)*sj(i, j, k, 3) - w(i, j-1, k, &
   &         ivx)*sj(i, j-1, k, 3) + w(i, j, k+1, ivx)*sk(i, j, k, 3) - w(i&
   &         , j, k-1, ivx)*sk(i, j, k-1, 3)
   ! Idem for the gradient of v.
   vx = w(i+1, j, k, ivy)*si(i, j, k, 1) - w(i-1, j, k, ivy)*si(i-1&
   &         , j, k, 1) + w(i, j+1, k, ivy)*sj(i, j, k, 1) - w(i, j-1, k, &
   &         ivy)*sj(i, j-1, k, 1) + w(i, j, k+1, ivy)*sk(i, j, k, 1) - w(i&
   &         , j, k-1, ivy)*sk(i, j, k-1, 1)
   vy = w(i+1, j, k, ivy)*si(i, j, k, 2) - w(i-1, j, k, ivy)*si(i-1&
   &         , j, k, 2) + w(i, j+1, k, ivy)*sj(i, j, k, 2) - w(i, j-1, k, &
   &         ivy)*sj(i, j-1, k, 2) + w(i, j, k+1, ivy)*sk(i, j, k, 2) - w(i&
   &         , j, k-1, ivy)*sk(i, j, k-1, 2)
   vz = w(i+1, j, k, ivy)*si(i, j, k, 3) - w(i-1, j, k, ivy)*si(i-1&
   &         , j, k, 3) + w(i, j+1, k, ivy)*sj(i, j, k, 3) - w(i, j-1, k, &
   &         ivy)*sj(i, j-1, k, 3) + w(i, j, k+1, ivy)*sk(i, j, k, 3) - w(i&
   &         , j, k-1, ivy)*sk(i, j, k-1, 3)
   ! And for the gradient of w.
   wx = w(i+1, j, k, ivz)*si(i, j, k, 1) - w(i-1, j, k, ivz)*si(i-1&
   &         , j, k, 1) + w(i, j+1, k, ivz)*sj(i, j, k, 1) - w(i, j-1, k, &
   &         ivz)*sj(i, j-1, k, 1) + w(i, j, k+1, ivz)*sk(i, j, k, 1) - w(i&
   &         , j, k-1, ivz)*sk(i, j, k-1, 1)
   wy = w(i+1, j, k, ivz)*si(i, j, k, 2) - w(i-1, j, k, ivz)*si(i-1&
   &         , j, k, 2) + w(i, j+1, k, ivz)*sj(i, j, k, 2) - w(i, j-1, k, &
   &         ivz)*sj(i, j-1, k, 2) + w(i, j, k+1, ivz)*sk(i, j, k, 2) - w(i&
   &         , j, k-1, ivz)*sk(i, j, k-1, 2)
   wz = w(i+1, j, k, ivz)*si(i, j, k, 3) - w(i-1, j, k, ivz)*si(i-1&
   &         , j, k, 3) + w(i, j+1, k, ivz)*sj(i, j, k, 3) - w(i, j-1, k, &
   &         ivz)*sj(i, j-1, k, 3) + w(i, j, k+1, ivz)*sk(i, j, k, 3) - w(i&
   &         , j, k-1, ivz)*sk(i, j, k-1, 3)
   ! Compute the strain and vorticity terms. The multiplication
   ! is present to obtain the correct gradients. Note that
   ! the wheel speed is substracted from the vorticity terms.
   fact = half/vol(i, j, k)
   CALL PUSHREAL8(qxx)
   qxx = fact*ux
   CALL PUSHREAL8(qyy)
   qyy = fact*vy
   CALL PUSHREAL8(qzz)
   qzz = fact*wz
   CALL PUSHREAL8(qxy)
   qxy = fact*half*(uy+vx)
   CALL PUSHREAL8(qxz)
   qxz = fact*half*(uz+wx)
   CALL PUSHREAL8(qyz)
   qyz = fact*half*(vz+wy)
   CALL PUSHREAL8(oxy)
   oxy = fact*half*(vx-uy) - omegaz
   CALL PUSHREAL8(oxz)
   oxz = fact*half*(uz-wx) - omegay
   CALL PUSHREAL8(oyz)
   oyz = fact*half*(wy-vz) - omegax
   ! Compute the summation of the strain and vorticity tensors.
   CALL PUSHREAL8(sijsij)
   sijsij = two*(qxy**2+qxz**2+qyz**2) + qxx**2 + qyy**2 + qzz**2
   ! Compute the production term.
   END DO
   END DO
   END DO
   DO k=kl,2,-1
   DO j=jl,2,-1
   DO i=il,2,-1
   oijoij = two*(oxy**2+oxz**2+oyz**2)
   IF (sijsij*oijoij .EQ. 0.0_8) THEN
   tempb = 0.0
   ELSE
   tempb = two*dwb(i, j, k, iprod)/(2.0*SQRT(sijsij*oijoij))
   END IF
   sijsijb = oijoij*tempb
   oijoijb = sijsij*tempb
   dwb(i, j, k, iprod) = 0.0_8
   tempb0 = two*oijoijb
   oxyb = 2*oxy*tempb0
   oxzb = 2*oxz*tempb0
   oyzb = 2*oyz*tempb0
   CALL POPREAL8(sijsij)
   tempb1 = two*sijsijb
   qxyb = 2*qxy*tempb1
   qxzb = 2*qxz*tempb1
   qyzb = 2*qyz*tempb1
   qxxb = 2*qxx*sijsijb
   qyyb = 2*qyy*sijsijb
   qzzb = 2*qzz*sijsijb
   fact = half/vol(i, j, k)
   CALL POPREAL8(oyz)
   tempb2 = fact*half*oyzb
   CALL POPREAL8(oxz)
   tempb4 = fact*half*oxzb
   CALL POPREAL8(oxy)
   tempb6 = fact*half*oxyb
   CALL POPREAL8(qyz)
   tempb3 = fact*half*qyzb
   wyb = tempb3 + tempb2
   vzb = tempb3 - tempb2
   CALL POPREAL8(qxz)
   tempb5 = fact*half*qxzb
   uzb = tempb5 + tempb4
   wxb = tempb5 - tempb4
   CALL POPREAL8(qxy)
   tempb7 = fact*half*qxyb
   vxb = tempb7 + tempb6
   uyb = tempb7 - tempb6
   CALL POPREAL8(qzz)
   wzb = fact*qzzb
   CALL POPREAL8(qyy)
   vyb = fact*qyyb
   CALL POPREAL8(qxx)
   uxb = fact*qxxb
   wb(i+1, j, k, ivz) = wb(i+1, j, k, ivz) + si(i, j, k, 3)*wzb
   wb(i-1, j, k, ivz) = wb(i-1, j, k, ivz) - si(i-1, j, k, 3)*wzb
   wb(i, j+1, k, ivz) = wb(i, j+1, k, ivz) + sj(i, j, k, 3)*wzb
   wb(i, j, k+1, ivz) = wb(i, j, k+1, ivz) + sk(i, j, k, 3)*wzb
   wb(i, j-1, k, ivz) = wb(i, j-1, k, ivz) - sj(i, j-1, k, 3)*wzb
   wb(i, j, k-1, ivz) = wb(i, j, k-1, ivz) - sk(i, j, k-1, 3)*wzb
   wb(i+1, j, k, ivz) = wb(i+1, j, k, ivz) + si(i, j, k, 2)*wyb
   wb(i-1, j, k, ivz) = wb(i-1, j, k, ivz) - si(i-1, j, k, 2)*wyb
   wb(i, j+1, k, ivz) = wb(i, j+1, k, ivz) + sj(i, j, k, 2)*wyb
   wb(i, j, k+1, ivz) = wb(i, j, k+1, ivz) + sk(i, j, k, 2)*wyb
   wb(i, j-1, k, ivz) = wb(i, j-1, k, ivz) - sj(i, j-1, k, 2)*wyb
   wb(i, j, k-1, ivz) = wb(i, j, k-1, ivz) - sk(i, j, k-1, 2)*wyb
   wb(i+1, j, k, ivz) = wb(i+1, j, k, ivz) + si(i, j, k, 1)*wxb
   wb(i-1, j, k, ivz) = wb(i-1, j, k, ivz) - si(i-1, j, k, 1)*wxb
   wb(i, j+1, k, ivz) = wb(i, j+1, k, ivz) + sj(i, j, k, 1)*wxb
   wb(i, j, k+1, ivz) = wb(i, j, k+1, ivz) + sk(i, j, k, 1)*wxb
   wb(i, j-1, k, ivz) = wb(i, j-1, k, ivz) - sj(i, j-1, k, 1)*wxb
   wb(i, j, k-1, ivz) = wb(i, j, k-1, ivz) - sk(i, j, k-1, 1)*wxb
   wb(i+1, j, k, ivy) = wb(i+1, j, k, ivy) + si(i, j, k, 3)*vzb
   wb(i-1, j, k, ivy) = wb(i-1, j, k, ivy) - si(i-1, j, k, 3)*vzb
   wb(i, j+1, k, ivy) = wb(i, j+1, k, ivy) + sj(i, j, k, 3)*vzb
   wb(i, j, k+1, ivy) = wb(i, j, k+1, ivy) + sk(i, j, k, 3)*vzb
   wb(i, j-1, k, ivy) = wb(i, j-1, k, ivy) - sj(i, j-1, k, 3)*vzb
   wb(i, j, k-1, ivy) = wb(i, j, k-1, ivy) - sk(i, j, k-1, 3)*vzb
   wb(i+1, j, k, ivy) = wb(i+1, j, k, ivy) + si(i, j, k, 2)*vyb
   wb(i-1, j, k, ivy) = wb(i-1, j, k, ivy) - si(i-1, j, k, 2)*vyb
   wb(i, j+1, k, ivy) = wb(i, j+1, k, ivy) + sj(i, j, k, 2)*vyb
   wb(i, j, k+1, ivy) = wb(i, j, k+1, ivy) + sk(i, j, k, 2)*vyb
   wb(i, j-1, k, ivy) = wb(i, j-1, k, ivy) - sj(i, j-1, k, 2)*vyb
   wb(i, j, k-1, ivy) = wb(i, j, k-1, ivy) - sk(i, j, k-1, 2)*vyb
   wb(i+1, j, k, ivy) = wb(i+1, j, k, ivy) + si(i, j, k, 1)*vxb
   wb(i-1, j, k, ivy) = wb(i-1, j, k, ivy) - si(i-1, j, k, 1)*vxb
   wb(i, j+1, k, ivy) = wb(i, j+1, k, ivy) + sj(i, j, k, 1)*vxb
   wb(i, j, k+1, ivy) = wb(i, j, k+1, ivy) + sk(i, j, k, 1)*vxb
   wb(i, j-1, k, ivy) = wb(i, j-1, k, ivy) - sj(i, j-1, k, 1)*vxb
   wb(i, j, k-1, ivy) = wb(i, j, k-1, ivy) - sk(i, j, k-1, 1)*vxb
   wb(i+1, j, k, ivx) = wb(i+1, j, k, ivx) + si(i, j, k, 3)*uzb
   wb(i-1, j, k, ivx) = wb(i-1, j, k, ivx) - si(i-1, j, k, 3)*uzb
   wb(i, j+1, k, ivx) = wb(i, j+1, k, ivx) + sj(i, j, k, 3)*uzb
   wb(i, j, k+1, ivx) = wb(i, j, k+1, ivx) + sk(i, j, k, 3)*uzb
   wb(i, j-1, k, ivx) = wb(i, j-1, k, ivx) - sj(i, j-1, k, 3)*uzb
   wb(i, j, k-1, ivx) = wb(i, j, k-1, ivx) - sk(i, j, k-1, 3)*uzb
   wb(i+1, j, k, ivx) = wb(i+1, j, k, ivx) + si(i, j, k, 2)*uyb
   wb(i-1, j, k, ivx) = wb(i-1, j, k, ivx) - si(i-1, j, k, 2)*uyb
   wb(i, j+1, k, ivx) = wb(i, j+1, k, ivx) + sj(i, j, k, 2)*uyb
   wb(i, j, k+1, ivx) = wb(i, j, k+1, ivx) + sk(i, j, k, 2)*uyb
   wb(i, j-1, k, ivx) = wb(i, j-1, k, ivx) - sj(i, j-1, k, 2)*uyb
   wb(i, j, k-1, ivx) = wb(i, j, k-1, ivx) - sk(i, j, k-1, 2)*uyb
   wb(i+1, j, k, ivx) = wb(i+1, j, k, ivx) + si(i, j, k, 1)*uxb
   wb(i-1, j, k, ivx) = wb(i-1, j, k, ivx) - si(i-1, j, k, 1)*uxb
   wb(i, j+1, k, ivx) = wb(i, j+1, k, ivx) + sj(i, j, k, 1)*uxb
   wb(i, j, k+1, ivx) = wb(i, j, k+1, ivx) + sk(i, j, k, 1)*uxb
   wb(i, j-1, k, ivx) = wb(i, j-1, k, ivx) - sj(i, j-1, k, 1)*uxb
   wb(i, j, k-1, ivx) = wb(i, j, k-1, ivx) - sk(i, j, k-1, 1)*uxb
   END DO
   END DO
   END DO
   END SUBROUTINE PRODKATOLAUNDER_B
