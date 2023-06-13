!        generated by tapenade     (inria, ecuador team)
!  tapenade 3.16 (develop) -  2 may 2023 13:20
!
module oversetutilities_d
  implicit none

contains
!  differentiation of fractoweights in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: weights
!   with respect to varying inputs: frac
!   rw status of diff variables: weights:out frac:in
! --------------------------------------------------
!           tapenade routine below this point
! --------------------------------------------------
  subroutine fractoweights_d(frac, fracd, weights, weightsd)
    use constants
    implicit none
    real(kind=realtype), dimension(3), intent(in) :: frac
    real(kind=realtype), dimension(3), intent(in) :: fracd
    real(kind=realtype), dimension(8), intent(out) :: weights
    real(kind=realtype), dimension(8), intent(out) :: weightsd
    real(kind=realtype) :: temp
    weightsd = 0.0_8
    temp = (one-frac(1))*(one-frac(2))
    weightsd(1) = (one-frac(3))*(-((one-frac(2))*fracd(1))-(one-frac(1))&
&     *fracd(2)) - temp*fracd(3)
    weights(1) = temp*(one-frac(3))
    temp = frac(1)*(one-frac(2))
    weightsd(2) = (one-frac(3))*((one-frac(2))*fracd(1)-frac(1)*fracd(2)&
&     ) - temp*fracd(3)
    weights(2) = temp*(one-frac(3))
    temp = (one-frac(1))*frac(2)
    weightsd(3) = (one-frac(3))*((one-frac(1))*fracd(2)-frac(2)*fracd(1)&
&     ) - temp*fracd(3)
    weights(3) = temp*(one-frac(3))
    weightsd(4) = (one-frac(3))*(frac(2)*fracd(1)+frac(1)*fracd(2)) - &
&     frac(1)*frac(2)*fracd(3)
    weights(4) = frac(1)*frac(2)*(one-frac(3))
    temp = (one-frac(1))*frac(3)
    weightsd(5) = (one-frac(2))*((one-frac(1))*fracd(3)-frac(3)*fracd(1)&
&     ) - temp*fracd(2)
    weights(5) = temp*(one-frac(2))
    weightsd(6) = (one-frac(2))*(frac(3)*fracd(1)+frac(1)*fracd(3)) - &
&     frac(1)*frac(3)*fracd(2)
    weights(6) = frac(1)*(one-frac(2))*frac(3)
    weightsd(7) = (one-frac(1))*(frac(3)*fracd(2)+frac(2)*fracd(3)) - &
&     frac(2)*frac(3)*fracd(1)
    weights(7) = (one-frac(1))*frac(2)*frac(3)
    weightsd(8) = frac(3)*(frac(2)*fracd(1)+frac(1)*fracd(2)) + frac(1)*&
&     frac(2)*fracd(3)
    weights(8) = frac(1)*frac(2)*frac(3)
  end subroutine fractoweights_d

! --------------------------------------------------
!           tapenade routine below this point
! --------------------------------------------------
  subroutine fractoweights(frac, weights)
    use constants
    implicit none
    real(kind=realtype), dimension(3), intent(in) :: frac
    real(kind=realtype), dimension(8), intent(out) :: weights
    weights(1) = (one-frac(1))*(one-frac(2))*(one-frac(3))
    weights(2) = frac(1)*(one-frac(2))*(one-frac(3))
    weights(3) = (one-frac(1))*frac(2)*(one-frac(3))
    weights(4) = frac(1)*frac(2)*(one-frac(3))
    weights(5) = (one-frac(1))*(one-frac(2))*frac(3)
    weights(6) = frac(1)*(one-frac(2))*frac(3)
    weights(7) = (one-frac(1))*frac(2)*frac(3)
    weights(8) = frac(1)*frac(2)*frac(3)
  end subroutine fractoweights

  subroutine fractoweights2(frac, weights)
    use constants
    implicit none
    real(kind=realtype), dimension(3), intent(in) :: frac
    real(kind=realtype), dimension(8), intent(out) :: weights
    weights(1) = (one-frac(1))*(one-frac(2))*(one-frac(3))
    weights(2) = frac(1)*(one-frac(2))*(one-frac(3))
    weights(3) = frac(1)*frac(2)*(one-frac(3))
    weights(4) = (one-frac(1))*frac(2)*(one-frac(3))
    weights(5) = (one-frac(1))*(one-frac(2))*frac(3)
    weights(6) = frac(1)*(one-frac(2))*frac(3)
    weights(7) = frac(1)*frac(2)*frac(3)
    weights(8) = (one-frac(1))*frac(2)*frac(3)
  end subroutine fractoweights2

!  differentiation of newtonupdate in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: frac
!   with respect to varying inputs: xcen blk frac
!   rw status of diff variables: xcen:in blk:in frac:in-out
  subroutine newtonupdate_d(xcen, xcend, blk, blkd, frac0, frac, fracd)
! this routine performs the newton update to recompute the new
! "frac" (u,v,w) for the point xcen. the actual search is performed
! on the the dual cell formed by the cell centers of the 3x3x3 block
! of primal nodes. this routine is ad'd with tapenade in both
! forward and reverse.
    use constants
    implicit none
! input
    real(kind=realtype), dimension(3), intent(in) :: xcen
    real(kind=realtype), dimension(3), intent(in) :: xcend
    real(kind=realtype), dimension(3, 3, 3, 3), intent(in) :: blk
    real(kind=realtype), dimension(3, 3, 3, 3), intent(in) :: blkd
    real(kind=realtype), dimension(3), intent(in) :: frac0
! output
    real(kind=realtype), dimension(3), intent(out) :: frac
    real(kind=realtype), dimension(3), intent(out) :: fracd
! working
    real(kind=realtype), dimension(3, 8) :: xn
    real(kind=realtype), dimension(3, 8) :: xnd
    real(kind=realtype) :: u, v, w, uv, uw, vw, wvu, du, dv, dw
    real(kind=realtype) :: ud, vd, wd, uvd, uwd, vwd, wvud, dud, dvd, &
&   dwd
    real(kind=realtype) :: a11, a12, a13, a21, a22, a23, a31, a32, a33, &
&   val
    real(kind=realtype) :: a11d, a12d, a13d, a21d, a22d, a23d, a31d, &
&   a32d, a33d, vald
    real(kind=realtype) :: f(3), x(3)
    real(kind=realtype) :: fd(3), xd(3)
    integer(kind=inttype), dimension(8), parameter :: indices=(/1, 2, 4&
&     , 3, 5, 6, 8, 7/)
    integer(kind=inttype) :: i, j, k, ii, ll
    real(kind=realtype), parameter :: adteps=1.e-25_realtype
    real(kind=realtype), parameter :: thresconv=1.e-10_realtype
    intrinsic sign
    intrinsic abs
    intrinsic max
    intrinsic sqrt
    real(kind=realtype) :: x1
    real(kind=realtype) :: x1d
    real(kind=realtype) :: max1
    real(kind=realtype) :: max1d
    real(kind=realtype) :: arg1
    real(kind=realtype) :: temp
    real(kind=realtype) :: temp0
    real(kind=realtype) :: temp1
    real(kind=realtype) :: temp2
! compute the cell center locations for the 8 nodes describing the
! dual cell. note that this must be counter-clockwise ordering.
    ii = 0
    xnd = 0.0_8
    do k=1,2
      do j=1,2
        do i=1,2
          ii = ii + 1
          xnd(:, indices(ii)) = eighth*(blkd(i, j, k, :)+blkd(i+1, j, k&
&           , :)+blkd(i, j+1, k, :)+blkd(i+1, j+1, k, :)+blkd(i, j, k+1&
&           , :)+blkd(i+1, j, k+1, :)+blkd(i, j+1, k+1, :)+blkd(i+1, j+1&
&           , k+1, :))
          xn(:, indices(ii)) = eighth*(blk(i, j, k, :)+blk(i+1, j, k, :)&
&           +blk(i, j+1, k, :)+blk(i+1, j+1, k, :)+blk(i, j, k+1, :)+blk&
&           (i+1, j, k+1, :)+blk(i, j+1, k+1, :)+blk(i+1, j+1, k+1, :))
        end do
      end do
    end do
! compute the coordinates relative to node 1.
    do i=2,8
      xnd(:, i) = xnd(:, i) - xnd(:, 1)
      xn(:, i) = xn(:, i) - xn(:, 1)
    end do
! compute the location of our seach point relative to the first node.
    xd = xcend - xnd(:, 1)
    x = xcen - xn(:, 1)
! modify the coordinates of node 3, 6, 8 and 7 such that
! they correspond to the weights of the u*v, u*w, v*w and
! u*v*w term in the transformation respectively.
    xnd(1, 7) = xnd(1, 7) + xnd(1, 2) + xnd(1, 4) + xnd(1, 5) - xnd(1, 3&
&     ) - xnd(1, 6) - xnd(1, 8)
    xn(1, 7) = xn(1, 7) + xn(1, 2) + xn(1, 4) + xn(1, 5) - xn(1, 3) - xn&
&     (1, 6) - xn(1, 8)
    xnd(2, 7) = xnd(2, 7) + xnd(2, 2) + xnd(2, 4) + xnd(2, 5) - xnd(2, 3&
&     ) - xnd(2, 6) - xnd(2, 8)
    xn(2, 7) = xn(2, 7) + xn(2, 2) + xn(2, 4) + xn(2, 5) - xn(2, 3) - xn&
&     (2, 6) - xn(2, 8)
    xnd(3, 7) = xnd(3, 7) + xnd(3, 2) + xnd(3, 4) + xnd(3, 5) - xnd(3, 3&
&     ) - xnd(3, 6) - xnd(3, 8)
    xn(3, 7) = xn(3, 7) + xn(3, 2) + xn(3, 4) + xn(3, 5) - xn(3, 3) - xn&
&     (3, 6) - xn(3, 8)
    xnd(1, 3) = xnd(1, 3) - xnd(1, 2) - xnd(1, 4)
    xn(1, 3) = xn(1, 3) - xn(1, 2) - xn(1, 4)
    xnd(2, 3) = xnd(2, 3) - xnd(2, 2) - xnd(2, 4)
    xn(2, 3) = xn(2, 3) - xn(2, 2) - xn(2, 4)
    xnd(3, 3) = xnd(3, 3) - xnd(3, 2) - xnd(3, 4)
    xn(3, 3) = xn(3, 3) - xn(3, 2) - xn(3, 4)
    xnd(1, 6) = xnd(1, 6) - xnd(1, 2) - xnd(1, 5)
    xn(1, 6) = xn(1, 6) - xn(1, 2) - xn(1, 5)
    xnd(2, 6) = xnd(2, 6) - xnd(2, 2) - xnd(2, 5)
    xn(2, 6) = xn(2, 6) - xn(2, 2) - xn(2, 5)
    xnd(3, 6) = xnd(3, 6) - xnd(3, 2) - xnd(3, 5)
    xn(3, 6) = xn(3, 6) - xn(3, 2) - xn(3, 5)
    xnd(1, 8) = xnd(1, 8) - xnd(1, 4) - xnd(1, 5)
    xn(1, 8) = xn(1, 8) - xn(1, 4) - xn(1, 5)
    xnd(2, 8) = xnd(2, 8) - xnd(2, 4) - xnd(2, 5)
    xn(2, 8) = xn(2, 8) - xn(2, 4) - xn(2, 5)
    xnd(3, 8) = xnd(3, 8) - xnd(3, 4) - xnd(3, 5)
    xn(3, 8) = xn(3, 8) - xn(3, 4) - xn(3, 5)
! set the starting values of u, v and w based on our previous values
    u = frac0(1)
    v = frac0(2)
    w = frac0(3)
    fd = 0.0_8
    ud = 0.0_8
    vd = 0.0_8
    wd = 0.0_8
! the newton algorithm to determine the parametric
! weights u, v and w for the given coordinate.
newtonhexa:do ll=1,15
! compute the rhs.
      uvd = v*ud + u*vd
      uv = u*v
      uwd = w*ud + u*wd
      uw = u*w
      vwd = w*vd + v*wd
      vw = v*w
      wvud = w*(v*ud+u*vd) + u*v*wd
      wvu = u*v*w
      fd(1) = u*xnd(1, 2) + xn(1, 2)*ud + v*xnd(1, 4) + xn(1, 4)*vd + w*&
&       xnd(1, 5) + xn(1, 5)*wd + uv*xnd(1, 3) + xn(1, 3)*uvd + uw*xnd(1&
&       , 6) + xn(1, 6)*uwd + vw*xnd(1, 8) + xn(1, 8)*vwd + wvu*xnd(1, 7&
&       ) + xn(1, 7)*wvud - xd(1)
      f(1) = xn(1, 2)*u + xn(1, 4)*v + xn(1, 5)*w + xn(1, 3)*uv + xn(1, &
&       6)*uw + xn(1, 8)*vw + xn(1, 7)*wvu - x(1)
      fd(2) = u*xnd(2, 2) + xn(2, 2)*ud + v*xnd(2, 4) + xn(2, 4)*vd + w*&
&       xnd(2, 5) + xn(2, 5)*wd + uv*xnd(2, 3) + xn(2, 3)*uvd + uw*xnd(2&
&       , 6) + xn(2, 6)*uwd + vw*xnd(2, 8) + xn(2, 8)*vwd + wvu*xnd(2, 7&
&       ) + xn(2, 7)*wvud - xd(2)
      f(2) = xn(2, 2)*u + xn(2, 4)*v + xn(2, 5)*w + xn(2, 3)*uv + xn(2, &
&       6)*uw + xn(2, 8)*vw + xn(2, 7)*wvu - x(2)
      fd(3) = u*xnd(3, 2) + xn(3, 2)*ud + v*xnd(3, 4) + xn(3, 4)*vd + w*&
&       xnd(3, 5) + xn(3, 5)*wd + uv*xnd(3, 3) + xn(3, 3)*uvd + uw*xnd(3&
&       , 6) + xn(3, 6)*uwd + vw*xnd(3, 8) + xn(3, 8)*vwd + wvu*xnd(3, 7&
&       ) + xn(3, 7)*wvud - xd(3)
      f(3) = xn(3, 2)*u + xn(3, 4)*v + xn(3, 5)*w + xn(3, 3)*uv + xn(3, &
&       6)*uw + xn(3, 8)*vw + xn(3, 7)*wvu - x(3)
! compute the jacobian.
      a11d = xnd(1, 2) + v*xnd(1, 3) + xn(1, 3)*vd + w*xnd(1, 6) + xn(1&
&       , 6)*wd + vw*xnd(1, 7) + xn(1, 7)*vwd
      a11 = xn(1, 2) + xn(1, 3)*v + xn(1, 6)*w + xn(1, 7)*vw
      a12d = xnd(1, 4) + u*xnd(1, 3) + xn(1, 3)*ud + w*xnd(1, 8) + xn(1&
&       , 8)*wd + uw*xnd(1, 7) + xn(1, 7)*uwd
      a12 = xn(1, 4) + xn(1, 3)*u + xn(1, 8)*w + xn(1, 7)*uw
      a13d = xnd(1, 5) + u*xnd(1, 6) + xn(1, 6)*ud + v*xnd(1, 8) + xn(1&
&       , 8)*vd + uv*xnd(1, 7) + xn(1, 7)*uvd
      a13 = xn(1, 5) + xn(1, 6)*u + xn(1, 8)*v + xn(1, 7)*uv
      a21d = xnd(2, 2) + v*xnd(2, 3) + xn(2, 3)*vd + w*xnd(2, 6) + xn(2&
&       , 6)*wd + vw*xnd(2, 7) + xn(2, 7)*vwd
      a21 = xn(2, 2) + xn(2, 3)*v + xn(2, 6)*w + xn(2, 7)*vw
      a22d = xnd(2, 4) + u*xnd(2, 3) + xn(2, 3)*ud + w*xnd(2, 8) + xn(2&
&       , 8)*wd + uw*xnd(2, 7) + xn(2, 7)*uwd
      a22 = xn(2, 4) + xn(2, 3)*u + xn(2, 8)*w + xn(2, 7)*uw
      a23d = xnd(2, 5) + u*xnd(2, 6) + xn(2, 6)*ud + v*xnd(2, 8) + xn(2&
&       , 8)*vd + uv*xnd(2, 7) + xn(2, 7)*uvd
      a23 = xn(2, 5) + xn(2, 6)*u + xn(2, 8)*v + xn(2, 7)*uv
      a31d = xnd(3, 2) + v*xnd(3, 3) + xn(3, 3)*vd + w*xnd(3, 6) + xn(3&
&       , 6)*wd + vw*xnd(3, 7) + xn(3, 7)*vwd
      a31 = xn(3, 2) + xn(3, 3)*v + xn(3, 6)*w + xn(3, 7)*vw
      a32d = xnd(3, 4) + u*xnd(3, 3) + xn(3, 3)*ud + w*xnd(3, 8) + xn(3&
&       , 8)*wd + uw*xnd(3, 7) + xn(3, 7)*uwd
      a32 = xn(3, 4) + xn(3, 3)*u + xn(3, 8)*w + xn(3, 7)*uw
      a33d = xnd(3, 5) + u*xnd(3, 6) + xn(3, 6)*ud + v*xnd(3, 8) + xn(3&
&       , 8)*vd + uv*xnd(3, 7) + xn(3, 7)*uvd
      a33 = xn(3, 5) + xn(3, 6)*u + xn(3, 8)*v + xn(3, 7)*uv
! compute the determinant. make sure that it is not zero
! and invert the value. the cut off is needed to be able
! to handle exceptional cases for degenerate elements.
      temp = a22*a33 - a32*a23
      temp0 = a13*a32 - a12*a33
      temp1 = a12*a23 - a13*a22
      vald = temp*a11d + a11*(a33*a22d+a22*a33d-a23*a32d-a32*a23d) + &
&       temp0*a21d + a21*(a32*a13d+a13*a32d-a33*a12d-a12*a33d) + temp1*&
&       a31d + a31*(a23*a12d+a12*a23d-a22*a13d-a13*a22d)
      val = a11*temp + a21*temp0 + a31*temp1
      if (val .ge. 0.) then
        x1d = vald
        x1 = val
      else
        x1d = -vald
        x1 = -val
      end if
      if (x1 .lt. adteps) then
        max1 = adteps
        max1d = 0.0_8
      else
        max1d = x1d
        max1 = x1
      end if
      temp1 = sign(one, val)/max1
      vald = -(temp1*max1d/max1)
      val = temp1
! compute the new values of u, v and w.
      temp1 = a12*a23 - a13*a22
      temp0 = a13*a32 - a12*a33
      temp = a22*a33 - a23*a32
      temp2 = temp*f(1) + temp0*f(2) + temp1*f(3)
      dud = temp2*vald + val*(f(1)*(a33*a22d+a22*a33d-a32*a23d-a23*a32d)&
&       +temp*fd(1)+f(2)*(a32*a13d+a13*a32d-a33*a12d-a12*a33d)+temp0*fd(&
&       2)+f(3)*(a23*a12d+a12*a23d-a22*a13d-a13*a22d)+temp1*fd(3))
      du = val*temp2
      temp2 = a13*a21 - a11*a23
      temp1 = a11*a33 - a13*a31
      temp0 = a23*a31 - a21*a33
      temp = temp0*f(1) + temp1*f(2) + temp2*f(3)
      dvd = temp*vald + val*(f(1)*(a31*a23d+a23*a31d-a33*a21d-a21*a33d)+&
&       temp0*fd(1)+f(2)*(a33*a11d+a11*a33d-a31*a13d-a13*a31d)+temp1*fd(&
&       2)+f(3)*(a21*a13d+a13*a21d-a23*a11d-a11*a23d)+temp2*fd(3))
      dv = val*temp
      temp2 = a11*a22 - a12*a21
      temp1 = a12*a31 - a11*a32
      temp0 = a21*a32 - a22*a31
      temp = temp0*f(1) + temp1*f(2) + temp2*f(3)
      dwd = temp*vald + val*(f(1)*(a32*a21d+a21*a32d-a31*a22d-a22*a31d)+&
&       temp0*fd(1)+f(2)*(a31*a12d+a12*a31d-a32*a11d-a11*a32d)+temp1*fd(&
&       2)+f(3)*(a22*a11d+a11*a22d-a21*a12d-a12*a21d)+temp2*fd(3))
      dw = val*temp
      ud = ud - dud
      u = u - du
      vd = vd - dvd
      v = v - dv
      wd = wd - dwd
      w = w - dw
! exit the loop if the update of the parametric
! weights is below the threshold
      arg1 = du*du + dv*dv + dw*dw
      val = sqrt(arg1)
      if (val .le. thresconv) exit
    end do newtonhexa
! we would *like* that all solutions fall inside the hexa, but we
! can't be picky here since we are not changing the donors. so
! whatever the u,v,w is we have to accept. even if it is greater than
! 1 or less than zero, it shouldn't be by much.
    fracd(1) = ud
    frac(1) = u
    fracd(2) = vd
    frac(2) = v
    fracd(3) = wd
    frac(3) = w
  end subroutine newtonupdate_d

  subroutine newtonupdate(xcen, blk, frac0, frac)
! this routine performs the newton update to recompute the new
! "frac" (u,v,w) for the point xcen. the actual search is performed
! on the the dual cell formed by the cell centers of the 3x3x3 block
! of primal nodes. this routine is ad'd with tapenade in both
! forward and reverse.
    use constants
    implicit none
! input
    real(kind=realtype), dimension(3), intent(in) :: xcen
    real(kind=realtype), dimension(3, 3, 3, 3), intent(in) :: blk
    real(kind=realtype), dimension(3), intent(in) :: frac0
! output
    real(kind=realtype), dimension(3), intent(out) :: frac
! working
    real(kind=realtype), dimension(3, 8) :: xn
    real(kind=realtype) :: u, v, w, uv, uw, vw, wvu, du, dv, dw
    real(kind=realtype) :: a11, a12, a13, a21, a22, a23, a31, a32, a33, &
&   val
    real(kind=realtype) :: f(3), x(3)
    integer(kind=inttype), dimension(8), parameter :: indices=(/1, 2, 4&
&     , 3, 5, 6, 8, 7/)
    integer(kind=inttype) :: i, j, k, ii, ll
    real(kind=realtype), parameter :: adteps=1.e-25_realtype
    real(kind=realtype), parameter :: thresconv=1.e-10_realtype
    intrinsic sign
    intrinsic abs
    intrinsic max
    intrinsic sqrt
    real(kind=realtype) :: x1
    real(kind=realtype) :: max1
    real(kind=realtype) :: arg1
! compute the cell center locations for the 8 nodes describing the
! dual cell. note that this must be counter-clockwise ordering.
    ii = 0
    do k=1,2
      do j=1,2
        do i=1,2
          ii = ii + 1
          xn(:, indices(ii)) = eighth*(blk(i, j, k, :)+blk(i+1, j, k, :)&
&           +blk(i, j+1, k, :)+blk(i+1, j+1, k, :)+blk(i, j, k+1, :)+blk&
&           (i+1, j, k+1, :)+blk(i, j+1, k+1, :)+blk(i+1, j+1, k+1, :))
        end do
      end do
    end do
! compute the coordinates relative to node 1.
    do i=2,8
      xn(:, i) = xn(:, i) - xn(:, 1)
    end do
! compute the location of our seach point relative to the first node.
    x = xcen - xn(:, 1)
! modify the coordinates of node 3, 6, 8 and 7 such that
! they correspond to the weights of the u*v, u*w, v*w and
! u*v*w term in the transformation respectively.
    xn(1, 7) = xn(1, 7) + xn(1, 2) + xn(1, 4) + xn(1, 5) - xn(1, 3) - xn&
&     (1, 6) - xn(1, 8)
    xn(2, 7) = xn(2, 7) + xn(2, 2) + xn(2, 4) + xn(2, 5) - xn(2, 3) - xn&
&     (2, 6) - xn(2, 8)
    xn(3, 7) = xn(3, 7) + xn(3, 2) + xn(3, 4) + xn(3, 5) - xn(3, 3) - xn&
&     (3, 6) - xn(3, 8)
    xn(1, 3) = xn(1, 3) - xn(1, 2) - xn(1, 4)
    xn(2, 3) = xn(2, 3) - xn(2, 2) - xn(2, 4)
    xn(3, 3) = xn(3, 3) - xn(3, 2) - xn(3, 4)
    xn(1, 6) = xn(1, 6) - xn(1, 2) - xn(1, 5)
    xn(2, 6) = xn(2, 6) - xn(2, 2) - xn(2, 5)
    xn(3, 6) = xn(3, 6) - xn(3, 2) - xn(3, 5)
    xn(1, 8) = xn(1, 8) - xn(1, 4) - xn(1, 5)
    xn(2, 8) = xn(2, 8) - xn(2, 4) - xn(2, 5)
    xn(3, 8) = xn(3, 8) - xn(3, 4) - xn(3, 5)
! set the starting values of u, v and w based on our previous values
    u = frac0(1)
    v = frac0(2)
    w = frac0(3)
! the newton algorithm to determine the parametric
! weights u, v and w for the given coordinate.
newtonhexa:do ll=1,15
! compute the rhs.
      uv = u*v
      uw = u*w
      vw = v*w
      wvu = u*v*w
      f(1) = xn(1, 2)*u + xn(1, 4)*v + xn(1, 5)*w + xn(1, 3)*uv + xn(1, &
&       6)*uw + xn(1, 8)*vw + xn(1, 7)*wvu - x(1)
      f(2) = xn(2, 2)*u + xn(2, 4)*v + xn(2, 5)*w + xn(2, 3)*uv + xn(2, &
&       6)*uw + xn(2, 8)*vw + xn(2, 7)*wvu - x(2)
      f(3) = xn(3, 2)*u + xn(3, 4)*v + xn(3, 5)*w + xn(3, 3)*uv + xn(3, &
&       6)*uw + xn(3, 8)*vw + xn(3, 7)*wvu - x(3)
! compute the jacobian.
      a11 = xn(1, 2) + xn(1, 3)*v + xn(1, 6)*w + xn(1, 7)*vw
      a12 = xn(1, 4) + xn(1, 3)*u + xn(1, 8)*w + xn(1, 7)*uw
      a13 = xn(1, 5) + xn(1, 6)*u + xn(1, 8)*v + xn(1, 7)*uv
      a21 = xn(2, 2) + xn(2, 3)*v + xn(2, 6)*w + xn(2, 7)*vw
      a22 = xn(2, 4) + xn(2, 3)*u + xn(2, 8)*w + xn(2, 7)*uw
      a23 = xn(2, 5) + xn(2, 6)*u + xn(2, 8)*v + xn(2, 7)*uv
      a31 = xn(3, 2) + xn(3, 3)*v + xn(3, 6)*w + xn(3, 7)*vw
      a32 = xn(3, 4) + xn(3, 3)*u + xn(3, 8)*w + xn(3, 7)*uw
      a33 = xn(3, 5) + xn(3, 6)*u + xn(3, 8)*v + xn(3, 7)*uv
! compute the determinant. make sure that it is not zero
! and invert the value. the cut off is needed to be able
! to handle exceptional cases for degenerate elements.
      val = a11*(a22*a33-a32*a23) + a21*(a13*a32-a12*a33) + a31*(a12*a23&
&       -a13*a22)
      if (val .ge. 0.) then
        x1 = val
      else
        x1 = -val
      end if
      if (x1 .lt. adteps) then
        max1 = adteps
      else
        max1 = x1
      end if
      val = sign(one, val)/max1
! compute the new values of u, v and w.
      du = val*((a22*a33-a23*a32)*f(1)+(a13*a32-a12*a33)*f(2)+(a12*a23-&
&       a13*a22)*f(3))
      dv = val*((a23*a31-a21*a33)*f(1)+(a11*a33-a13*a31)*f(2)+(a13*a21-&
&       a11*a23)*f(3))
      dw = val*((a21*a32-a22*a31)*f(1)+(a12*a31-a11*a32)*f(2)+(a11*a22-&
&       a12*a21)*f(3))
      u = u - du
      v = v - dv
      w = w - dw
! exit the loop if the update of the parametric
! weights is below the threshold
      arg1 = du*du + dv*dv + dw*dw
      val = sqrt(arg1)
      if (val .le. thresconv) goto 100
    end do newtonhexa
! we would *like* that all solutions fall inside the hexa, but we
! can't be picky here since we are not changing the donors. so
! whatever the u,v,w is we have to accept. even if it is greater than
! 1 or less than zero, it shouldn't be by much.
 100 frac(1) = u
    frac(2) = v
    frac(3) = w
  end subroutine newtonupdate

end module oversetutilities_d

