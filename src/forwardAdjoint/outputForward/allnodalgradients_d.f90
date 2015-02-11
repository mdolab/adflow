!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of allnodalgradients in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: *wx *wy *wz *qx *qy *qz *ux
!                *uy *uz *vx *vy *vz
!   with respect to varying inputs: *aa *w *vol *si *sj *sk
!   plus diff mem management of: aa:in wx:in wy:in wz:in w:in qx:in
!                qy:in qz:in ux:in vol:in uy:in uz:in si:in sj:in
!                sk:in vx:in vy:in vz:in
subroutine allnodalgradients_d()
!
!        ****************************************************************
!        *                                                              *
!        * nodalgradients computes the nodal velocity gradients and     *
!        * minus the gradient of the speed of sound squared. the minus  *
!        * sign is present, because this is the definition of the heat  *
!        * flux. these gradients are computed for all nodes.            * 
!        *                                                              *
!        ****************************************************************
!
  use blockpointers
  implicit none
!        local variables.
  integer(kind=inttype) :: i, j, k
  integer(kind=inttype) :: k1, kk
  integer(kind=inttype) :: istart, iend, isize, ii
  integer(kind=inttype) :: jstart, jend, jsize
  integer(kind=inttype) :: kstart, kend, ksize
  real(kind=realtype) :: oneoverv, ubar, vbar, wbar, a2
  real(kind=realtype) :: oneovervd, ubard, vbard, wbard, a2d
  real(kind=realtype) :: sx, sx1, sy, sy1, sz, sz1
  real(kind=realtype) :: sxd, syd, szd
!
!        ****************************************************************
!        *                                                              *
!        * begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
! zero all nodeal gradients:
  uxd = 0.0_8
  ux = zero
  uyd = 0.0_8
  uy = zero
  uzd = 0.0_8
  uz = zero
  vxd = 0.0_8
  vx = zero
  vyd = 0.0_8
  vy = zero
  vzd = 0.0_8
  vz = zero
  wxd = 0.0_8
  wx = zero
  wyd = 0.0_8
  wy = zero
  wzd = 0.0_8
  wz = zero
  qxd = 0.0_8
  qx = zero
  qyd = 0.0_8
  qy = zero
  qzd = 0.0_8
  qz = zero
  wxd = 0.0_8
  wyd = 0.0_8
  wzd = 0.0_8
  qxd = 0.0_8
  qyd = 0.0_8
  qzd = 0.0_8
  uxd = 0.0_8
  uyd = 0.0_8
  uzd = 0.0_8
  vxd = 0.0_8
  vyd = 0.0_8
  vzd = 0.0_8
! first part. contribution in the k-direction.
! the contribution is scattered to both the left and right node
! in k-direction.
  do k=1,ke
    do j=1,jl
      do i=1,il
! compute 8 times the average normal for this part of
! the control volume. the factor 8 is taken care of later
! on when the division by the volume takes place.
        sxd = skd(i, j, k-1, 1) + skd(i+1, j, k-1, 1) + skd(i, j+1, k-1&
&         , 1) + skd(i+1, j+1, k-1, 1) + skd(i, j, k, 1) + skd(i+1, j, k&
&         , 1) + skd(i, j+1, k, 1) + skd(i+1, j+1, k, 1)
        sx = sk(i, j, k-1, 1) + sk(i+1, j, k-1, 1) + sk(i, j+1, k-1, 1) &
&         + sk(i+1, j+1, k-1, 1) + sk(i, j, k, 1) + sk(i+1, j, k, 1) + &
&         sk(i, j+1, k, 1) + sk(i+1, j+1, k, 1)
        syd = skd(i, j, k-1, 2) + skd(i+1, j, k-1, 2) + skd(i, j+1, k-1&
&         , 2) + skd(i+1, j+1, k-1, 2) + skd(i, j, k, 2) + skd(i+1, j, k&
&         , 2) + skd(i, j+1, k, 2) + skd(i+1, j+1, k, 2)
        sy = sk(i, j, k-1, 2) + sk(i+1, j, k-1, 2) + sk(i, j+1, k-1, 2) &
&         + sk(i+1, j+1, k-1, 2) + sk(i, j, k, 2) + sk(i+1, j, k, 2) + &
&         sk(i, j+1, k, 2) + sk(i+1, j+1, k, 2)
        szd = skd(i, j, k-1, 3) + skd(i+1, j, k-1, 3) + skd(i, j+1, k-1&
&         , 3) + skd(i+1, j+1, k-1, 3) + skd(i, j, k, 3) + skd(i+1, j, k&
&         , 3) + skd(i, j+1, k, 3) + skd(i+1, j+1, k, 3)
        sz = sk(i, j, k-1, 3) + sk(i+1, j, k-1, 3) + sk(i, j+1, k-1, 3) &
&         + sk(i+1, j+1, k-1, 3) + sk(i, j, k, 3) + sk(i+1, j, k, 3) + &
&         sk(i, j+1, k, 3) + sk(i+1, j+1, k, 3)
! compute the average velocities and speed of sound squared
! for this integration point. node that these variables are
! stored in w(ivx), w(ivy), w(ivz) and p.
        ubard = fourth*(wd(i, j, k, ivx)+wd(i+1, j, k, ivx)+wd(i, j+1, k&
&         , ivx)+wd(i+1, j+1, k, ivx))
        ubar = fourth*(w(i, j, k, ivx)+w(i+1, j, k, ivx)+w(i, j+1, k, &
&         ivx)+w(i+1, j+1, k, ivx))
        vbard = fourth*(wd(i, j, k, ivy)+wd(i+1, j, k, ivy)+wd(i, j+1, k&
&         , ivy)+wd(i+1, j+1, k, ivy))
        vbar = fourth*(w(i, j, k, ivy)+w(i+1, j, k, ivy)+w(i, j+1, k, &
&         ivy)+w(i+1, j+1, k, ivy))
        wbard = fourth*(wd(i, j, k, ivz)+wd(i+1, j, k, ivz)+wd(i, j+1, k&
&         , ivz)+wd(i+1, j+1, k, ivz))
        wbar = fourth*(w(i, j, k, ivz)+w(i+1, j, k, ivz)+w(i, j+1, k, &
&         ivz)+w(i+1, j+1, k, ivz))
        a2d = fourth*(aad(i, j, k)+aad(i+1, j, k)+aad(i, j+1, k)+aad(i+1&
&         , j+1, k))
        a2 = fourth*(aa(i, j, k)+aa(i+1, j, k)+aa(i, j+1, k)+aa(i+1, j+1&
&         , k))
! add the contributions to the surface integral to the node
! j-1 and substract it from the node j. for the heat flux it
! is reversed, because the negative of the gradient of the
! speed of sound must be computed.
        if (k .gt. 1) then
          uxd(i, j, k-1) = uxd(i, j, k-1) + ubard*sx + ubar*sxd
          ux(i, j, k-1) = ux(i, j, k-1) + ubar*sx
          uyd(i, j, k-1) = uyd(i, j, k-1) + ubard*sy + ubar*syd
          uy(i, j, k-1) = uy(i, j, k-1) + ubar*sy
          uzd(i, j, k-1) = uzd(i, j, k-1) + ubard*sz + ubar*szd
          uz(i, j, k-1) = uz(i, j, k-1) + ubar*sz
          vxd(i, j, k-1) = vxd(i, j, k-1) + vbard*sx + vbar*sxd
          vx(i, j, k-1) = vx(i, j, k-1) + vbar*sx
          vyd(i, j, k-1) = vyd(i, j, k-1) + vbard*sy + vbar*syd
          vy(i, j, k-1) = vy(i, j, k-1) + vbar*sy
          vzd(i, j, k-1) = vzd(i, j, k-1) + vbard*sz + vbar*szd
          vz(i, j, k-1) = vz(i, j, k-1) + vbar*sz
          wxd(i, j, k-1) = wxd(i, j, k-1) + wbard*sx + wbar*sxd
          wx(i, j, k-1) = wx(i, j, k-1) + wbar*sx
          wyd(i, j, k-1) = wyd(i, j, k-1) + wbard*sy + wbar*syd
          wy(i, j, k-1) = wy(i, j, k-1) + wbar*sy
          wzd(i, j, k-1) = wzd(i, j, k-1) + wbard*sz + wbar*szd
          wz(i, j, k-1) = wz(i, j, k-1) + wbar*sz
          qxd(i, j, k-1) = qxd(i, j, k-1) - a2d*sx - a2*sxd
          qx(i, j, k-1) = qx(i, j, k-1) - a2*sx
          qyd(i, j, k-1) = qyd(i, j, k-1) - a2d*sy - a2*syd
          qy(i, j, k-1) = qy(i, j, k-1) - a2*sy
          qzd(i, j, k-1) = qzd(i, j, k-1) - a2d*sz - a2*szd
          qz(i, j, k-1) = qz(i, j, k-1) - a2*sz
        end if
        if (k .lt. ke) then
          uxd(i, j, k) = uxd(i, j, k) - ubard*sx - ubar*sxd
          ux(i, j, k) = ux(i, j, k) - ubar*sx
          uyd(i, j, k) = uyd(i, j, k) - ubard*sy - ubar*syd
          uy(i, j, k) = uy(i, j, k) - ubar*sy
          uzd(i, j, k) = uzd(i, j, k) - ubard*sz - ubar*szd
          uz(i, j, k) = uz(i, j, k) - ubar*sz
          vxd(i, j, k) = vxd(i, j, k) - vbard*sx - vbar*sxd
          vx(i, j, k) = vx(i, j, k) - vbar*sx
          vyd(i, j, k) = vyd(i, j, k) - vbard*sy - vbar*syd
          vy(i, j, k) = vy(i, j, k) - vbar*sy
          vzd(i, j, k) = vzd(i, j, k) - vbard*sz - vbar*szd
          vz(i, j, k) = vz(i, j, k) - vbar*sz
          wxd(i, j, k) = wxd(i, j, k) - wbard*sx - wbar*sxd
          wx(i, j, k) = wx(i, j, k) - wbar*sx
          wyd(i, j, k) = wyd(i, j, k) - wbard*sy - wbar*syd
          wy(i, j, k) = wy(i, j, k) - wbar*sy
          wzd(i, j, k) = wzd(i, j, k) - wbard*sz - wbar*szd
          wz(i, j, k) = wz(i, j, k) - wbar*sz
          qxd(i, j, k) = qxd(i, j, k) + a2d*sx + a2*sxd
          qx(i, j, k) = qx(i, j, k) + a2*sx
          qyd(i, j, k) = qyd(i, j, k) + a2d*sy + a2*syd
          qy(i, j, k) = qy(i, j, k) + a2*sy
          qzd(i, j, k) = qzd(i, j, k) + a2d*sz + a2*szd
          qz(i, j, k) = qz(i, j, k) + a2*sz
        end if
      end do
    end do
  end do
! second part. contribution in the j-direction.
! the contribution is scattered to both the left and right node
! in j-direction.
  do k=1,kl
    do j=1,je
      do i=1,il
! compute 8 times the average normal for this part of
! the control volume. the factor 8 is taken care of later
! on when the division by the volume takes place.
        sxd = sjd(i, j-1, k, 1) + sjd(i+1, j-1, k, 1) + sjd(i, j-1, k+1&
&         , 1) + sjd(i+1, j-1, k+1, 1) + sjd(i, j, k, 1) + sjd(i+1, j, k&
&         , 1) + sjd(i, j, k+1, 1) + sjd(i+1, j, k+1, 1)
        sx = sj(i, j-1, k, 1) + sj(i+1, j-1, k, 1) + sj(i, j-1, k+1, 1) &
&         + sj(i+1, j-1, k+1, 1) + sj(i, j, k, 1) + sj(i+1, j, k, 1) + &
&         sj(i, j, k+1, 1) + sj(i+1, j, k+1, 1)
        syd = sjd(i, j-1, k, 2) + sjd(i+1, j-1, k, 2) + sjd(i, j-1, k+1&
&         , 2) + sjd(i+1, j-1, k+1, 2) + sjd(i, j, k, 2) + sjd(i+1, j, k&
&         , 2) + sjd(i, j, k+1, 2) + sjd(i+1, j, k+1, 2)
        sy = sj(i, j-1, k, 2) + sj(i+1, j-1, k, 2) + sj(i, j-1, k+1, 2) &
&         + sj(i+1, j-1, k+1, 2) + sj(i, j, k, 2) + sj(i+1, j, k, 2) + &
&         sj(i, j, k+1, 2) + sj(i+1, j, k+1, 2)
        szd = sjd(i, j-1, k, 3) + sjd(i+1, j-1, k, 3) + sjd(i, j-1, k+1&
&         , 3) + sjd(i+1, j-1, k+1, 3) + sjd(i, j, k, 3) + sjd(i+1, j, k&
&         , 3) + sjd(i, j, k+1, 3) + sjd(i+1, j, k+1, 3)
        sz = sj(i, j-1, k, 3) + sj(i+1, j-1, k, 3) + sj(i, j-1, k+1, 3) &
&         + sj(i+1, j-1, k+1, 3) + sj(i, j, k, 3) + sj(i+1, j, k, 3) + &
&         sj(i, j, k+1, 3) + sj(i+1, j, k+1, 3)
! compute the average velocities and speed of sound squared
! for this integration point. node that these variables are
! stored in w(ivx), w(ivy), w(ivz) and p.
        ubard = fourth*(wd(i, j, k, ivx)+wd(i+1, j, k, ivx)+wd(i, j, k+1&
&         , ivx)+wd(i+1, j, k+1, ivx))
        ubar = fourth*(w(i, j, k, ivx)+w(i+1, j, k, ivx)+w(i, j, k+1, &
&         ivx)+w(i+1, j, k+1, ivx))
        vbard = fourth*(wd(i, j, k, ivy)+wd(i+1, j, k, ivy)+wd(i, j, k+1&
&         , ivy)+wd(i+1, j, k+1, ivy))
        vbar = fourth*(w(i, j, k, ivy)+w(i+1, j, k, ivy)+w(i, j, k+1, &
&         ivy)+w(i+1, j, k+1, ivy))
        wbard = fourth*(wd(i, j, k, ivz)+wd(i+1, j, k, ivz)+wd(i, j, k+1&
&         , ivz)+wd(i+1, j, k+1, ivz))
        wbar = fourth*(w(i, j, k, ivz)+w(i+1, j, k, ivz)+w(i, j, k+1, &
&         ivz)+w(i+1, j, k+1, ivz))
        a2d = fourth*(aad(i, j, k)+aad(i+1, j, k)+aad(i, j, k+1)+aad(i+1&
&         , j, k+1))
        a2 = fourth*(aa(i, j, k)+aa(i+1, j, k)+aa(i, j, k+1)+aa(i+1, j, &
&         k+1))
! add the contributions to the surface integral to the node
! j-1 and substract it from the node j. for the heat flux it
! is reversed, because the negative of the gradient of the
! speed of sound must be computed.
        if (j .gt. 1) then
          uxd(i, j-1, k) = uxd(i, j-1, k) + ubard*sx + ubar*sxd
          ux(i, j-1, k) = ux(i, j-1, k) + ubar*sx
          uyd(i, j-1, k) = uyd(i, j-1, k) + ubard*sy + ubar*syd
          uy(i, j-1, k) = uy(i, j-1, k) + ubar*sy
          uzd(i, j-1, k) = uzd(i, j-1, k) + ubard*sz + ubar*szd
          uz(i, j-1, k) = uz(i, j-1, k) + ubar*sz
          vxd(i, j-1, k) = vxd(i, j-1, k) + vbard*sx + vbar*sxd
          vx(i, j-1, k) = vx(i, j-1, k) + vbar*sx
          vyd(i, j-1, k) = vyd(i, j-1, k) + vbard*sy + vbar*syd
          vy(i, j-1, k) = vy(i, j-1, k) + vbar*sy
          vzd(i, j-1, k) = vzd(i, j-1, k) + vbard*sz + vbar*szd
          vz(i, j-1, k) = vz(i, j-1, k) + vbar*sz
          wxd(i, j-1, k) = wxd(i, j-1, k) + wbard*sx + wbar*sxd
          wx(i, j-1, k) = wx(i, j-1, k) + wbar*sx
          wyd(i, j-1, k) = wyd(i, j-1, k) + wbard*sy + wbar*syd
          wy(i, j-1, k) = wy(i, j-1, k) + wbar*sy
          wzd(i, j-1, k) = wzd(i, j-1, k) + wbard*sz + wbar*szd
          wz(i, j-1, k) = wz(i, j-1, k) + wbar*sz
          qxd(i, j-1, k) = qxd(i, j-1, k) - a2d*sx - a2*sxd
          qx(i, j-1, k) = qx(i, j-1, k) - a2*sx
          qyd(i, j-1, k) = qyd(i, j-1, k) - a2d*sy - a2*syd
          qy(i, j-1, k) = qy(i, j-1, k) - a2*sy
          qzd(i, j-1, k) = qzd(i, j-1, k) - a2d*sz - a2*szd
          qz(i, j-1, k) = qz(i, j-1, k) - a2*sz
        end if
        if (j .lt. je) then
          uxd(i, j, k) = uxd(i, j, k) - ubard*sx - ubar*sxd
          ux(i, j, k) = ux(i, j, k) - ubar*sx
          uyd(i, j, k) = uyd(i, j, k) - ubard*sy - ubar*syd
          uy(i, j, k) = uy(i, j, k) - ubar*sy
          uzd(i, j, k) = uzd(i, j, k) - ubard*sz - ubar*szd
          uz(i, j, k) = uz(i, j, k) - ubar*sz
          vxd(i, j, k) = vxd(i, j, k) - vbard*sx - vbar*sxd
          vx(i, j, k) = vx(i, j, k) - vbar*sx
          vyd(i, j, k) = vyd(i, j, k) - vbard*sy - vbar*syd
          vy(i, j, k) = vy(i, j, k) - vbar*sy
          vzd(i, j, k) = vzd(i, j, k) - vbard*sz - vbar*szd
          vz(i, j, k) = vz(i, j, k) - vbar*sz
          wxd(i, j, k) = wxd(i, j, k) - wbard*sx - wbar*sxd
          wx(i, j, k) = wx(i, j, k) - wbar*sx
          wyd(i, j, k) = wyd(i, j, k) - wbard*sy - wbar*syd
          wy(i, j, k) = wy(i, j, k) - wbar*sy
          wzd(i, j, k) = wzd(i, j, k) - wbard*sz - wbar*szd
          wz(i, j, k) = wz(i, j, k) - wbar*sz
          qxd(i, j, k) = qxd(i, j, k) + a2d*sx + a2*sxd
          qx(i, j, k) = qx(i, j, k) + a2*sx
          qyd(i, j, k) = qyd(i, j, k) + a2d*sy + a2*syd
          qy(i, j, k) = qy(i, j, k) + a2*sy
          qzd(i, j, k) = qzd(i, j, k) + a2d*sz + a2*szd
          qz(i, j, k) = qz(i, j, k) + a2*sz
        end if
      end do
    end do
  end do
! third part. contribution in the i-direction.
! the contribution is scattered to both the left and right node
! in i-direction.
  do k=1,kl
    do j=1,jl
      do i=1,ie
! compute 8 times the average normal for this part of
! the control volume. the factor 8 is taken care of later
! on when the division by the volume takes place.
        sxd = sid(i-1, j, k, 1) + sid(i-1, j+1, k, 1) + sid(i-1, j, k+1&
&         , 1) + sid(i-1, j+1, k+1, 1) + sid(i, j, k, 1) + sid(i, j+1, k&
&         , 1) + sid(i, j, k+1, 1) + sid(i, j+1, k+1, 1)
        sx = si(i-1, j, k, 1) + si(i-1, j+1, k, 1) + si(i-1, j, k+1, 1) &
&         + si(i-1, j+1, k+1, 1) + si(i, j, k, 1) + si(i, j+1, k, 1) + &
&         si(i, j, k+1, 1) + si(i, j+1, k+1, 1)
        syd = sid(i-1, j, k, 2) + sid(i-1, j+1, k, 2) + sid(i-1, j, k+1&
&         , 2) + sid(i-1, j+1, k+1, 2) + sid(i, j, k, 2) + sid(i, j+1, k&
&         , 2) + sid(i, j, k+1, 2) + sid(i, j+1, k+1, 2)
        sy = si(i-1, j, k, 2) + si(i-1, j+1, k, 2) + si(i-1, j, k+1, 2) &
&         + si(i-1, j+1, k+1, 2) + si(i, j, k, 2) + si(i, j+1, k, 2) + &
&         si(i, j, k+1, 2) + si(i, j+1, k+1, 2)
        szd = sid(i-1, j, k, 3) + sid(i-1, j+1, k, 3) + sid(i-1, j, k+1&
&         , 3) + sid(i-1, j+1, k+1, 3) + sid(i, j, k, 3) + sid(i, j+1, k&
&         , 3) + sid(i, j, k+1, 3) + sid(i, j+1, k+1, 3)
        sz = si(i-1, j, k, 3) + si(i-1, j+1, k, 3) + si(i-1, j, k+1, 3) &
&         + si(i-1, j+1, k+1, 3) + si(i, j, k, 3) + si(i, j+1, k, 3) + &
&         si(i, j, k+1, 3) + si(i, j+1, k+1, 3)
! compute the average velocities and speed of sound squared
! for this integration point. node that these variables are
! stored in w(ivx), w(ivy), w(ivz) and p.
        ubard = fourth*(wd(i, j, k, ivx)+wd(i, j+1, k, ivx)+wd(i, j, k+1&
&         , ivx)+wd(i, j+1, k+1, ivx))
        ubar = fourth*(w(i, j, k, ivx)+w(i, j+1, k, ivx)+w(i, j, k+1, &
&         ivx)+w(i, j+1, k+1, ivx))
        vbard = fourth*(wd(i, j, k, ivy)+wd(i, j+1, k, ivy)+wd(i, j, k+1&
&         , ivy)+wd(i, j+1, k+1, ivy))
        vbar = fourth*(w(i, j, k, ivy)+w(i, j+1, k, ivy)+w(i, j, k+1, &
&         ivy)+w(i, j+1, k+1, ivy))
        wbard = fourth*(wd(i, j, k, ivz)+wd(i, j+1, k, ivz)+wd(i, j, k+1&
&         , ivz)+wd(i, j+1, k+1, ivz))
        wbar = fourth*(w(i, j, k, ivz)+w(i, j+1, k, ivz)+w(i, j, k+1, &
&         ivz)+w(i, j+1, k+1, ivz))
        a2d = fourth*(aad(i, j, k)+aad(i, j+1, k)+aad(i, j, k+1)+aad(i, &
&         j+1, k+1))
        a2 = fourth*(aa(i, j, k)+aa(i, j+1, k)+aa(i, j, k+1)+aa(i, j+1, &
&         k+1))
! add the contributions to the surface integral to the node
! j-1 and substract it from the node j. for the heat flux it
! is reversed, because the negative of the gradient of the
! speed of sound must be computed.
        if (i .gt. 1) then
          uxd(i-1, j, k) = uxd(i-1, j, k) + ubard*sx + ubar*sxd
          ux(i-1, j, k) = ux(i-1, j, k) + ubar*sx
          uyd(i-1, j, k) = uyd(i-1, j, k) + ubard*sy + ubar*syd
          uy(i-1, j, k) = uy(i-1, j, k) + ubar*sy
          uzd(i-1, j, k) = uzd(i-1, j, k) + ubard*sz + ubar*szd
          uz(i-1, j, k) = uz(i-1, j, k) + ubar*sz
          vxd(i-1, j, k) = vxd(i-1, j, k) + vbard*sx + vbar*sxd
          vx(i-1, j, k) = vx(i-1, j, k) + vbar*sx
          vyd(i-1, j, k) = vyd(i-1, j, k) + vbard*sy + vbar*syd
          vy(i-1, j, k) = vy(i-1, j, k) + vbar*sy
          vzd(i-1, j, k) = vzd(i-1, j, k) + vbard*sz + vbar*szd
          vz(i-1, j, k) = vz(i-1, j, k) + vbar*sz
          wxd(i-1, j, k) = wxd(i-1, j, k) + wbard*sx + wbar*sxd
          wx(i-1, j, k) = wx(i-1, j, k) + wbar*sx
          wyd(i-1, j, k) = wyd(i-1, j, k) + wbard*sy + wbar*syd
          wy(i-1, j, k) = wy(i-1, j, k) + wbar*sy
          wzd(i-1, j, k) = wzd(i-1, j, k) + wbard*sz + wbar*szd
          wz(i-1, j, k) = wz(i-1, j, k) + wbar*sz
          qxd(i-1, j, k) = qxd(i-1, j, k) - a2d*sx - a2*sxd
          qx(i-1, j, k) = qx(i-1, j, k) - a2*sx
          qyd(i-1, j, k) = qyd(i-1, j, k) - a2d*sy - a2*syd
          qy(i-1, j, k) = qy(i-1, j, k) - a2*sy
          qzd(i-1, j, k) = qzd(i-1, j, k) - a2d*sz - a2*szd
          qz(i-1, j, k) = qz(i-1, j, k) - a2*sz
        end if
        if (i .lt. ie) then
          uxd(i, j, k) = uxd(i, j, k) - ubard*sx - ubar*sxd
          ux(i, j, k) = ux(i, j, k) - ubar*sx
          uyd(i, j, k) = uyd(i, j, k) - ubard*sy - ubar*syd
          uy(i, j, k) = uy(i, j, k) - ubar*sy
          uzd(i, j, k) = uzd(i, j, k) - ubard*sz - ubar*szd
          uz(i, j, k) = uz(i, j, k) - ubar*sz
          vxd(i, j, k) = vxd(i, j, k) - vbard*sx - vbar*sxd
          vx(i, j, k) = vx(i, j, k) - vbar*sx
          vyd(i, j, k) = vyd(i, j, k) - vbard*sy - vbar*syd
          vy(i, j, k) = vy(i, j, k) - vbar*sy
          vzd(i, j, k) = vzd(i, j, k) - vbard*sz - vbar*szd
          vz(i, j, k) = vz(i, j, k) - vbar*sz
          wxd(i, j, k) = wxd(i, j, k) - wbard*sx - wbar*sxd
          wx(i, j, k) = wx(i, j, k) - wbar*sx
          wyd(i, j, k) = wyd(i, j, k) - wbard*sy - wbar*syd
          wy(i, j, k) = wy(i, j, k) - wbar*sy
          wzd(i, j, k) = wzd(i, j, k) - wbard*sz - wbar*szd
          wz(i, j, k) = wz(i, j, k) - wbar*sz
          qxd(i, j, k) = qxd(i, j, k) + a2d*sx + a2*sxd
          qx(i, j, k) = qx(i, j, k) + a2*sx
          qyd(i, j, k) = qyd(i, j, k) + a2d*sy + a2*syd
          qy(i, j, k) = qy(i, j, k) + a2*sy
          qzd(i, j, k) = qzd(i, j, k) + a2d*sz + a2*szd
          qz(i, j, k) = qz(i, j, k) + a2*sz
        end if
      end do
    end do
  end do
! divide by 8 times the volume to obtain the correct gradients.
  do k=1,kl
    do j=1,jl
      do i=1,il
! compute the inverse of 8 times the volume for this node.
        oneovervd = -(one*(vold(i, j, k)+vold(i, j, k+1)+vold(i+1, j, k)&
&         +vold(i+1, j, k+1)+vold(i, j+1, k)+vold(i, j+1, k+1)+vold(i+1&
&         , j+1, k)+vold(i+1, j+1, k+1))/(vol(i, j, k)+vol(i, j, k+1)+&
&         vol(i+1, j, k)+vol(i+1, j, k+1)+vol(i, j+1, k)+vol(i, j+1, k+1&
&         )+vol(i+1, j+1, k)+vol(i+1, j+1, k+1))**2)
        oneoverv = one/(vol(i, j, k)+vol(i, j, k+1)+vol(i+1, j, k)+vol(i&
&         +1, j, k+1)+vol(i, j+1, k)+vol(i, j+1, k+1)+vol(i+1, j+1, k)+&
&         vol(i+1, j+1, k+1))
! compute the correct velocity gradients and "unit" heat
! fluxes. the velocity gradients are stored in ux, etc.
        uxd(i, j, k) = uxd(i, j, k)*oneoverv + ux(i, j, k)*oneovervd
        ux(i, j, k) = ux(i, j, k)*oneoverv
        uyd(i, j, k) = uyd(i, j, k)*oneoverv + uy(i, j, k)*oneovervd
        uy(i, j, k) = uy(i, j, k)*oneoverv
        uzd(i, j, k) = uzd(i, j, k)*oneoverv + uz(i, j, k)*oneovervd
        uz(i, j, k) = uz(i, j, k)*oneoverv
        vxd(i, j, k) = vxd(i, j, k)*oneoverv + vx(i, j, k)*oneovervd
        vx(i, j, k) = vx(i, j, k)*oneoverv
        vyd(i, j, k) = vyd(i, j, k)*oneoverv + vy(i, j, k)*oneovervd
        vy(i, j, k) = vy(i, j, k)*oneoverv
        vzd(i, j, k) = vzd(i, j, k)*oneoverv + vz(i, j, k)*oneovervd
        vz(i, j, k) = vz(i, j, k)*oneoverv
        wxd(i, j, k) = wxd(i, j, k)*oneoverv + wx(i, j, k)*oneovervd
        wx(i, j, k) = wx(i, j, k)*oneoverv
        wyd(i, j, k) = wyd(i, j, k)*oneoverv + wy(i, j, k)*oneovervd
        wy(i, j, k) = wy(i, j, k)*oneoverv
        wzd(i, j, k) = wzd(i, j, k)*oneoverv + wz(i, j, k)*oneovervd
        wz(i, j, k) = wz(i, j, k)*oneoverv
        qxd(i, j, k) = qxd(i, j, k)*oneoverv + qx(i, j, k)*oneovervd
        qx(i, j, k) = qx(i, j, k)*oneoverv
        qyd(i, j, k) = qyd(i, j, k)*oneoverv + qy(i, j, k)*oneovervd
        qy(i, j, k) = qy(i, j, k)*oneoverv
        qzd(i, j, k) = qzd(i, j, k)*oneoverv + qz(i, j, k)*oneovervd
        qz(i, j, k) = qz(i, j, k)*oneoverv
      end do
    end do
  end do
end subroutine allnodalgradients_d
