!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of volume_block in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: *vol
!   with respect to varying inputs: *x
!   plus diff mem management of: x:in vol:in
subroutine volume_block_d()
! this is copy of metric.f90. it was necessary to copy this file
! since there is debugging stuff in the original that is not
! necessary for ad.
  use bctypes
  use blockpointers
  use cgnsgrid
  use communication
  use inputtimespectral
  implicit none
!
!      local parameter.
!
  real(kind=realtype), parameter :: thresvolume=1.e-2_realtype
  real(kind=realtype), parameter :: halocellratio=1e-10_realtype
!
!      local variables.
!
  integer(kind=inttype) :: i, j, k, n, m, l, ii
  integer(kind=inttype) :: mm
  real(kind=realtype) :: fact, mult
  real(kind=realtype) :: xp, yp, zp, vp1, vp2, vp3, vp4, vp5, vp6
  real(kind=realtype) :: xpd, ypd, zpd, vp1d, vp2d, vp3d, vp4d, vp5d, &
& vp6d
  real(kind=realtype) :: xxp, yyp, zzp
  real(kind=realtype), dimension(3) :: v1, v2
  intrinsic abs
!
!      ******************************************************************
!      *                                                                *
!      * begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
! compute the volumes. the hexahedron is split into 6 pyramids
! whose volumes are computed. the volume is positive for a
! right handed block.
! initialize the volumes to zero. the reasons is that the second
! level halo's must be initialized to zero and for convenience
! all the volumes are set to zero.
  vold = 0.0_8
  vol = zero
  vold = 0.0_8
  vp1d = 0.0_8
  vp2d = 0.0_8
  vp3d = 0.0_8
  vp4d = 0.0_8
  vp5d = 0.0_8
  vp6d = 0.0_8
  do k=1,ke
    n = k - 1
    do j=1,je
      m = j - 1
      do i=1,ie
        l = i - 1
! compute the coordinates of the center of gravity.
        xpd = eighth*(xd(i, j, k, 1)+xd(i, m, k, 1)+xd(i, m, n, 1)+xd(i&
&         , j, n, 1)+xd(l, j, k, 1)+xd(l, m, k, 1)+xd(l, m, n, 1)+xd(l, &
&         j, n, 1))
        xp = eighth*(x(i, j, k, 1)+x(i, m, k, 1)+x(i, m, n, 1)+x(i, j, n&
&         , 1)+x(l, j, k, 1)+x(l, m, k, 1)+x(l, m, n, 1)+x(l, j, n, 1))
        ypd = eighth*(xd(i, j, k, 2)+xd(i, m, k, 2)+xd(i, m, n, 2)+xd(i&
&         , j, n, 2)+xd(l, j, k, 2)+xd(l, m, k, 2)+xd(l, m, n, 2)+xd(l, &
&         j, n, 2))
        yp = eighth*(x(i, j, k, 2)+x(i, m, k, 2)+x(i, m, n, 2)+x(i, j, n&
&         , 2)+x(l, j, k, 2)+x(l, m, k, 2)+x(l, m, n, 2)+x(l, j, n, 2))
        zpd = eighth*(xd(i, j, k, 3)+xd(i, m, k, 3)+xd(i, m, n, 3)+xd(i&
&         , j, n, 3)+xd(l, j, k, 3)+xd(l, m, k, 3)+xd(l, m, n, 3)+xd(l, &
&         j, n, 3))
        zp = eighth*(x(i, j, k, 3)+x(i, m, k, 3)+x(i, m, n, 3)+x(i, j, n&
&         , 3)+x(l, j, k, 3)+x(l, m, k, 3)+x(l, m, n, 3)+x(l, j, n, 3))
! compute the volumes of the 6 sub pyramids. the
! arguments of volpym must be such that for a (regular)
! right handed hexahedron all volumes are positive.
        call volpym_d(x(i, j, k, 1), xd(i, j, k, 1), x(i, j, k, 2), xd(i&
&               , j, k, 2), x(i, j, k, 3), xd(i, j, k, 3), x(i, j, n, 1)&
&               , xd(i, j, n, 1), x(i, j, n, 2), xd(i, j, n, 2), x(i, j&
&               , n, 3), xd(i, j, n, 3), x(i, m, n, 1), xd(i, m, n, 1), &
&               x(i, m, n, 2), xd(i, m, n, 2), x(i, m, n, 3), xd(i, m, n&
&               , 3), x(i, m, k, 1), xd(i, m, k, 1), x(i, m, k, 2), xd(i&
&               , m, k, 2), x(i, m, k, 3), xd(i, m, k, 3), vp1, vp1d)
        call volpym_d(x(l, j, k, 1), xd(l, j, k, 1), x(l, j, k, 2), xd(l&
&               , j, k, 2), x(l, j, k, 3), xd(l, j, k, 3), x(l, m, k, 1)&
&               , xd(l, m, k, 1), x(l, m, k, 2), xd(l, m, k, 2), x(l, m&
&               , k, 3), xd(l, m, k, 3), x(l, m, n, 1), xd(l, m, n, 1), &
&               x(l, m, n, 2), xd(l, m, n, 2), x(l, m, n, 3), xd(l, m, n&
&               , 3), x(l, j, n, 1), xd(l, j, n, 1), x(l, j, n, 2), xd(l&
&               , j, n, 2), x(l, j, n, 3), xd(l, j, n, 3), vp2, vp2d)
        call volpym_d(x(i, j, k, 1), xd(i, j, k, 1), x(i, j, k, 2), xd(i&
&               , j, k, 2), x(i, j, k, 3), xd(i, j, k, 3), x(l, j, k, 1)&
&               , xd(l, j, k, 1), x(l, j, k, 2), xd(l, j, k, 2), x(l, j&
&               , k, 3), xd(l, j, k, 3), x(l, j, n, 1), xd(l, j, n, 1), &
&               x(l, j, n, 2), xd(l, j, n, 2), x(l, j, n, 3), xd(l, j, n&
&               , 3), x(i, j, n, 1), xd(i, j, n, 1), x(i, j, n, 2), xd(i&
&               , j, n, 2), x(i, j, n, 3), xd(i, j, n, 3), vp3, vp3d)
        call volpym_d(x(i, m, k, 1), xd(i, m, k, 1), x(i, m, k, 2), xd(i&
&               , m, k, 2), x(i, m, k, 3), xd(i, m, k, 3), x(i, m, n, 1)&
&               , xd(i, m, n, 1), x(i, m, n, 2), xd(i, m, n, 2), x(i, m&
&               , n, 3), xd(i, m, n, 3), x(l, m, n, 1), xd(l, m, n, 1), &
&               x(l, m, n, 2), xd(l, m, n, 2), x(l, m, n, 3), xd(l, m, n&
&               , 3), x(l, m, k, 1), xd(l, m, k, 1), x(l, m, k, 2), xd(l&
&               , m, k, 2), x(l, m, k, 3), xd(l, m, k, 3), vp4, vp4d)
        call volpym_d(x(i, j, k, 1), xd(i, j, k, 1), x(i, j, k, 2), xd(i&
&               , j, k, 2), x(i, j, k, 3), xd(i, j, k, 3), x(i, m, k, 1)&
&               , xd(i, m, k, 1), x(i, m, k, 2), xd(i, m, k, 2), x(i, m&
&               , k, 3), xd(i, m, k, 3), x(l, m, k, 1), xd(l, m, k, 1), &
&               x(l, m, k, 2), xd(l, m, k, 2), x(l, m, k, 3), xd(l, m, k&
&               , 3), x(l, j, k, 1), xd(l, j, k, 1), x(l, j, k, 2), xd(l&
&               , j, k, 2), x(l, j, k, 3), xd(l, j, k, 3), vp5, vp5d)
        call volpym_d(x(i, j, n, 1), xd(i, j, n, 1), x(i, j, n, 2), xd(i&
&               , j, n, 2), x(i, j, n, 3), xd(i, j, n, 3), x(l, j, n, 1)&
&               , xd(l, j, n, 1), x(l, j, n, 2), xd(l, j, n, 2), x(l, j&
&               , n, 3), xd(l, j, n, 3), x(l, m, n, 1), xd(l, m, n, 1), &
&               x(l, m, n, 2), xd(l, m, n, 2), x(l, m, n, 3), xd(l, m, n&
&               , 3), x(i, m, n, 1), xd(i, m, n, 1), x(i, m, n, 2), xd(i&
&               , m, n, 2), x(i, m, n, 3), xd(i, m, n, 3), vp6, vp6d)
! set the volume to 1/6 of the sum of the volumes of the
! pyramid. remember that volpym computes 6 times the
! volume.
        vold(i, j, k) = sixth*(vp1d+vp2d+vp3d+vp4d+vp5d+vp6d)
        vol(i, j, k) = sixth*(vp1+vp2+vp3+vp4+vp5+vp6)
        if (vol(i, j, k) .ge. 0.) then
          vol(i, j, k) = vol(i, j, k)
        else
          vold(i, j, k) = -vold(i, j, k)
          vol(i, j, k) = -vol(i, j, k)
        end if
      end do
    end do
  end do
! some additional safety stuff for halo volumes.
  do k=2,kl
    do j=2,jl
      if (vol(1, j, k)/vol(2, j, k) .lt. halocellratio) then
        vold(1, j, k) = vold(2, j, k)
        vol(1, j, k) = vol(2, j, k)
      end if
      if (vol(ie, j, k)/vol(il, j, k) .lt. halocellratio) then
        vold(ie, j, k) = vold(il, j, k)
        vol(ie, j, k) = vol(il, j, k)
      end if
    end do
  end do
  do k=2,kl
    do i=1,ie
      if (vol(i, 1, k)/vol(i, 2, k) .lt. halocellratio) then
        vold(i, 1, k) = vold(i, 2, k)
        vol(i, 1, k) = vol(i, 2, k)
      end if
      if (vol(i, je, k)/vol(i, jl, k) .lt. halocellratio) then
        vold(i, je, k) = vold(i, jl, k)
        vol(i, je, k) = vol(i, jl, k)
      end if
    end do
  end do
  do j=1,je
    do i=1,ie
      if (vol(i, j, 1)/vol(i, j, 2) .lt. halocellratio) then
        vold(i, j, 1) = vold(i, j, 2)
        vol(i, j, 1) = vol(i, j, 2)
      end if
      if (vol(i, j, ke)/vol(i, j, kl) .lt. halocellratio) then
        vold(i, j, ke) = vold(i, j, kl)
        vol(i, j, ke) = vol(i, j, kl)
      end if
    end do
  end do

contains
!  differentiation of volpym in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: volume
!   with respect to varying inputs: xp yp zp xa xb xc xd ya yb
!                yc yd za zb zc zd
  subroutine volpym_d(xa, xad, ya, yad, za, zad, xb, xbd, yb, ybd, zb, &
&   zbd, xc, xcd, yc, ycd, zc, zcd, xd, xdd, yd, ydd, zd, zdd, volume, &
&   volumed)
!
!        ****************************************************************
!        *                                                              *
!        * volpym computes 6 times the volume of a pyramid. node p,     *
!        * whose coordinates are set in the subroutine metric itself,   *
!        * is the top node and a-b-c-d is the quadrilateral surface.    *
!        * it is assumed that the cross product vca * vdb points in     *
!        * the direction of the top node. here vca is the diagonal      *
!        * running from node c to node a and vdb the diagonal from      *
!        * node d to node b.                                            *
!        *                                                              *
!        ****************************************************************
!
    use precision
    implicit none
!
!        function type.
!
    real(kind=realtype) :: volume
    real(kind=realtype) :: volumed
!
!        function arguments.
!
    real(kind=realtype), intent(in) :: xa, ya, za, xb, yb, zb
    real(kind=realtype), intent(in) :: xad, yad, zad, xbd, ybd, zbd
    real(kind=realtype), intent(in) :: xc, yc, zc, xd, yd, zd
    real(kind=realtype), intent(in) :: xcd, ycd, zcd, xdd, ydd, zdd
!
!        ****************************************************************
!        *                                                              *
!        * begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
    volumed = (xpd-fourth*(xad+xbd+xcd+xdd))*((ya-yc)*(zb-zd)-(za-zc)*(&
&     yb-yd)) + (xp-fourth*(xa+xb+xc+xd))*((yad-ycd)*(zb-zd)+(ya-yc)*(&
&     zbd-zdd)-(zad-zcd)*(yb-yd)-(za-zc)*(ybd-ydd)) + (ypd-fourth*(yad+&
&     ybd+ycd+ydd))*((za-zc)*(xb-xd)-(xa-xc)*(zb-zd)) + (yp-fourth*(ya+&
&     yb+yc+yd))*((zad-zcd)*(xb-xd)+(za-zc)*(xbd-xdd)-(xad-xcd)*(zb-zd)-&
&     (xa-xc)*(zbd-zdd)) + (zpd-fourth*(zad+zbd+zcd+zdd))*((xa-xc)*(yb-&
&     yd)-(ya-yc)*(xb-xd)) + (zp-fourth*(za+zb+zc+zd))*((xad-xcd)*(yb-yd&
&     )+(xa-xc)*(ybd-ydd)-(yad-ycd)*(xb-xd)-(ya-yc)*(xbd-xdd))
    volume = (xp-fourth*(xa+xb+xc+xd))*((ya-yc)*(zb-zd)-(za-zc)*(yb-yd))&
&     + (yp-fourth*(ya+yb+yc+yd))*((za-zc)*(xb-xd)-(xa-xc)*(zb-zd)) + (&
&     zp-fourth*(za+zb+zc+zd))*((xa-xc)*(yb-yd)-(ya-yc)*(xb-xd))
  end subroutine volpym_d
  subroutine volpym(xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd, &
&   volume)
!
!        ****************************************************************
!        *                                                              *
!        * volpym computes 6 times the volume of a pyramid. node p,     *
!        * whose coordinates are set in the subroutine metric itself,   *
!        * is the top node and a-b-c-d is the quadrilateral surface.    *
!        * it is assumed that the cross product vca * vdb points in     *
!        * the direction of the top node. here vca is the diagonal      *
!        * running from node c to node a and vdb the diagonal from      *
!        * node d to node b.                                            *
!        *                                                              *
!        ****************************************************************
!
    use precision
    implicit none
!
!        function type.
!
    real(kind=realtype) :: volume
!
!        function arguments.
!
    real(kind=realtype), intent(in) :: xa, ya, za, xb, yb, zb
    real(kind=realtype), intent(in) :: xc, yc, zc, xd, yd, zd
!
!        ****************************************************************
!        *                                                              *
!        * begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
    volume = (xp-fourth*(xa+xb+xc+xd))*((ya-yc)*(zb-zd)-(za-zc)*(yb-yd))&
&     + (yp-fourth*(ya+yb+yc+yd))*((za-zc)*(xb-xd)-(xa-xc)*(zb-zd)) + (&
&     zp-fourth*(za+zb+zc+zd))*((xa-xc)*(yb-yd)-(ya-yc)*(xb-xd))
  end subroutine volpym
end subroutine volume_block_d
