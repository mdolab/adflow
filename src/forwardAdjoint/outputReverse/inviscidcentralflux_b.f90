!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of inviscidcentralflux in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *p *dw *w *vol *si *sj *sk
!   with respect to varying inputs: timeref *p *dw *w *vol *si
!                *sj *sk
!   plus diff mem management of: p:in dw:in w:in vol:in si:in sj:in
!                sk:in
!
!      ******************************************************************
!      *                                                                *
!      * file:          inviscidcentralflux.f90                         *
!      * author:        edwin van der weide                             *
!      * starting date: 03-24-2003                                      *
!      * last modified: 10-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine inviscidcentralflux_b()
!
!      ******************************************************************
!      *                                                                *
!      * inviscidcentralflux computes the euler fluxes using a central  *
!      * discretization for a given block. therefore it is assumed that *
!      * the pointers in block pointer already point to the correct     *
!      * block on the correct multigrid level.                          *
!      *                                                                *
!      ******************************************************************
!
  use blockpointers
  use cgnsgrid
  use constants
  use flowvarrefstate
  use inputphysics
  implicit none
!
!      local variables.
!
  integer(kind=inttype) :: i, j, k, ind, ii
  real(kind=realtype) :: qsp, qsm, rqsp, rqsm, porvel, porflux
  real(kind=realtype) :: qspd, qsmd, rqspd, rqsmd
  real(kind=realtype) :: pa, fs, sface, vnp, vnm
  real(kind=realtype) :: pad, fsd, vnpd, vnmd
  real(kind=realtype) :: wwx, wwy, wwz, rvol
  real(kind=realtype) :: wwxd, wwyd, wwzd, rvold
  intrinsic mod
  integer :: branch
  real(kind=realtype) :: temp3
  real(kind=realtype) :: temp2
  real(kind=realtype) :: temp1
  real(kind=realtype) :: temp0
  real(kind=realtype) :: tempd
  real(kind=realtype) :: tempd4
  real(kind=realtype) :: tempd3
  real(kind=realtype) :: tempd2
  real(kind=realtype) :: tempd1
  real(kind=realtype) :: tempd0
  real(kind=realtype) :: temp
  real(kind=realtype) :: temp4
  call pushinteger4(i)
  call pushinteger4(j)
  call pushreal8(vnm)
  call pushreal8(vnp)
  call pushreal8(porflux)
  call pushreal8(porvel)
  call pushreal8(qsm)
  call pushreal8(qsp)
  call pushinteger4(i)
  call pushinteger4(j)
  call pushreal8(vnm)
  call pushreal8(vnp)
  call pushreal8(porflux)
  call pushreal8(porvel)
  call pushreal8(qsm)
  call pushreal8(qsp)
  call pushinteger4(i)
  call pushinteger4(j)
  timerefd = 0.0_8
  if (blockismoving .and. equationmode .eq. steady) then
! compute the three nondimensional angular velocities.
    wwx = timeref*cgnsdoms(nbkglobal)%rotrate(1)
    wwy = timeref*cgnsdoms(nbkglobal)%rotrate(2)
    wwz = timeref*cgnsdoms(nbkglobal)%rotrate(3)
    wwxd = 0.0_8
    wwyd = 0.0_8
    wwzd = 0.0_8
    do ii=0,nx*ny*nz-1
      i = mod(ii, nx) + 2
      j = mod(ii/nx, ny) + 2
      k = ii/(nx*ny) + 2
      rvol = w(i, j, k, irho)*vol(i, j, k)
      temp4 = w(i, j, k, ivx)
      temp3 = w(i, j, k, ivy)
      tempd2 = rvol*dwd(i, j, k, imz)
      rvold = (wwx*temp3-wwy*temp4)*dwd(i, j, k, imz)
      wd(i, j, k, ivy) = wd(i, j, k, ivy) + wwx*tempd2
      wd(i, j, k, ivx) = wd(i, j, k, ivx) - wwy*tempd2
      temp2 = w(i, j, k, ivz)
      temp1 = w(i, j, k, ivx)
      tempd3 = rvol*dwd(i, j, k, imy)
      wwxd = wwxd + temp3*tempd2 - temp2*tempd3
      rvold = rvold + (wwz*temp1-wwx*temp2)*dwd(i, j, k, imy)
      wd(i, j, k, ivx) = wd(i, j, k, ivx) + wwz*tempd3
      wd(i, j, k, ivz) = wd(i, j, k, ivz) - wwx*tempd3
      temp0 = w(i, j, k, ivy)
      temp = w(i, j, k, ivz)
      tempd4 = rvol*dwd(i, j, k, imx)
      wwyd = wwyd + temp*tempd4 - temp4*tempd2
      wwzd = wwzd + temp1*tempd3 - temp0*tempd4
      rvold = rvold + (wwy*temp-wwz*temp0)*dwd(i, j, k, imx)
      wd(i, j, k, ivz) = wd(i, j, k, ivz) + wwy*tempd4
      wd(i, j, k, ivy) = wd(i, j, k, ivy) - wwz*tempd4
      wd(i, j, k, irho) = wd(i, j, k, irho) + vol(i, j, k)*rvold
      vold(i, j, k) = vold(i, j, k) + w(i, j, k, irho)*rvold
    end do
    timerefd = cgnsdoms(nbkglobal)%rotrate(2)*wwyd + cgnsdoms(nbkglobal)&
&     %rotrate(1)*wwxd + cgnsdoms(nbkglobal)%rotrate(3)*wwzd
  else
    timerefd = 0.0_8
  end if
  call popinteger4(j)
  call popinteger4(i)
  sface = zero
  do ii=0,nx*ny*kl-1
    i = mod(ii, nx) + 2
    j = mod(ii/nx, ny) + 2
    k = ii/(nx*ny) + 1
! set the dot product of the grid velocity and the
! normal in k-direction for a moving face.
    if (addgridvelocities) sface = sfacek(i, j, k)
! compute the normal velocities of the left and right state.
    vnp = w(i, j, k+1, ivx)*sk(i, j, k, 1) + w(i, j, k+1, ivy)*sk(i, j, &
&     k, 2) + w(i, j, k+1, ivz)*sk(i, j, k, 3)
    vnm = w(i, j, k, ivx)*sk(i, j, k, 1) + w(i, j, k, ivy)*sk(i, j, k, 2&
&     ) + w(i, j, k, ivz)*sk(i, j, k, 3)
! set the values of the porosities for this face.
! porvel defines the porosity w.r.t. velocity;
! porflux defines the porosity w.r.t. the entire flux.
! the latter is only zero for a discontinuous block
! block boundary that must be treated conservatively.
! the default value of porflux is 0.5, such that the
! correct central flux is scattered to both cells.
! in case of a boundflux the normal velocity is set
! to sface.
    porvel = one
    porflux = half
    if (pork(i, j, k) .eq. noflux) porflux = zero
    if (pork(i, j, k) .eq. boundflux) then
      porvel = zero
      vnp = sface
      vnm = sface
      call pushcontrol1b(0)
    else
      call pushcontrol1b(1)
    end if
! incorporate porflux in porvel.
    porvel = porvel*porflux
! compute the normal velocities for the face as well as the
! mass fluxes.
    qsp = (vnp-sface)*porvel
    qsm = (vnm-sface)*porvel
    rqsp = qsp*w(i, j, k+1, irho)
    rqsm = qsm*w(i, j, k, irho)
! compute the sum of the pressure multiplied by porflux.
! for the default value of porflux, 0.5, this leads to
! the average pressure.
    pa = porflux*(p(i, j, k+1)+p(i, j, k))
! compute the fluxes and scatter them to the cells
! i,j,k and i,j,k+1. store the density flux in the
! mass flow of the appropriate sliding mesh interface.
    fsd = dwd(i, j, k, irhoe) - dwd(i, j, k+1, irhoe)
    tempd1 = porflux*fsd
    qspd = w(i, j, k+1, irhoe)*fsd
    wd(i, j, k+1, irhoe) = wd(i, j, k+1, irhoe) + qsp*fsd
    qsmd = w(i, j, k, irhoe)*fsd
    wd(i, j, k, irhoe) = wd(i, j, k, irhoe) + qsm*fsd
    pd(i, j, k+1) = pd(i, j, k+1) + vnp*tempd1
    pd(i, j, k) = pd(i, j, k) + vnm*tempd1
    fsd = dwd(i, j, k, imz) - dwd(i, j, k+1, imz)
    rqspd = w(i, j, k+1, ivz)*fsd
    wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + rqsp*fsd
    rqsmd = w(i, j, k, ivz)*fsd
    wd(i, j, k, ivz) = wd(i, j, k, ivz) + rqsm*fsd
    pad = sk(i, j, k, 3)*fsd
    skd(i, j, k, 3) = skd(i, j, k, 3) + pa*fsd
    fsd = dwd(i, j, k, imy) - dwd(i, j, k+1, imy)
    rqspd = rqspd + w(i, j, k+1, ivy)*fsd
    wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + rqsp*fsd
    rqsmd = rqsmd + w(i, j, k, ivy)*fsd
    wd(i, j, k, ivy) = wd(i, j, k, ivy) + rqsm*fsd
    pad = pad + sk(i, j, k, 2)*fsd
    skd(i, j, k, 2) = skd(i, j, k, 2) + pa*fsd
    fsd = dwd(i, j, k, imx) - dwd(i, j, k+1, imx)
    rqspd = rqspd + w(i, j, k+1, ivx)*fsd
    wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + rqsp*fsd
    rqsmd = rqsmd + w(i, j, k, ivx)*fsd
    wd(i, j, k, ivx) = wd(i, j, k, ivx) + rqsm*fsd
    pad = pad + sk(i, j, k, 1)*fsd
    skd(i, j, k, 1) = skd(i, j, k, 1) + pa*fsd
    fsd = dwd(i, j, k, irho) - dwd(i, j, k+1, irho)
    rqspd = rqspd + fsd
    rqsmd = rqsmd + fsd
    pd(i, j, k+1) = pd(i, j, k+1) + porflux*pad
    pd(i, j, k) = pd(i, j, k) + porflux*pad
    qsmd = qsmd + w(i, j, k, irho)*rqsmd
    vnmd = porvel*qsmd + p(i, j, k)*tempd1
    wd(i, j, k, irho) = wd(i, j, k, irho) + qsm*rqsmd
    qspd = qspd + w(i, j, k+1, irho)*rqspd
    vnpd = porvel*qspd + p(i, j, k+1)*tempd1
    wd(i, j, k+1, irho) = wd(i, j, k+1, irho) + qsp*rqspd
    call popcontrol1b(branch)
    if (branch .eq. 0) then
      vnmd = 0.0_8
      vnpd = 0.0_8
    end if
    wd(i, j, k, ivx) = wd(i, j, k, ivx) + sk(i, j, k, 1)*vnmd
    skd(i, j, k, 1) = skd(i, j, k, 1) + w(i, j, k, ivx)*vnmd
    wd(i, j, k, ivy) = wd(i, j, k, ivy) + sk(i, j, k, 2)*vnmd
    skd(i, j, k, 2) = skd(i, j, k, 2) + w(i, j, k, ivy)*vnmd
    wd(i, j, k, ivz) = wd(i, j, k, ivz) + sk(i, j, k, 3)*vnmd
    skd(i, j, k, 3) = skd(i, j, k, 3) + w(i, j, k, ivz)*vnmd
    wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + sk(i, j, k, 1)*vnpd
    skd(i, j, k, 1) = skd(i, j, k, 1) + w(i, j, k+1, ivx)*vnpd
    wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + sk(i, j, k, 2)*vnpd
    skd(i, j, k, 2) = skd(i, j, k, 2) + w(i, j, k+1, ivy)*vnpd
    wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + sk(i, j, k, 3)*vnpd
    skd(i, j, k, 3) = skd(i, j, k, 3) + w(i, j, k+1, ivz)*vnpd
  end do
  call popreal8(qsp)
  call popreal8(qsm)
  call popreal8(porvel)
  call popreal8(porflux)
  call popreal8(vnp)
  call popreal8(vnm)
  call popinteger4(j)
  call popinteger4(i)
  sface = zero
  do ii=0,nx*jl*nz-1
    i = mod(ii, nx) + 2
    j = mod(ii/nx, jl) + 1
    k = ii/(nx*jl) + 2
! set the dot product of the grid velocity and the
! normal in j-direction for a moving face.
    if (addgridvelocities) sface = sfacej(i, j, k)
! compute the normal velocities of the left and right state.
    vnp = w(i, j+1, k, ivx)*sj(i, j, k, 1) + w(i, j+1, k, ivy)*sj(i, j, &
&     k, 2) + w(i, j+1, k, ivz)*sj(i, j, k, 3)
    vnm = w(i, j, k, ivx)*sj(i, j, k, 1) + w(i, j, k, ivy)*sj(i, j, k, 2&
&     ) + w(i, j, k, ivz)*sj(i, j, k, 3)
! set the values of the porosities for this face.
! porvel defines the porosity w.r.t. velocity;
! porflux defines the porosity w.r.t. the entire flux.
! the latter is only zero for a discontinuous block
! boundary that must be treated conservatively.
! the default value of porflux is 0.5, such that the
! correct central flux is scattered to both cells.
! in case of a boundflux the normal velocity is set
! to sface.
    porvel = one
    porflux = half
    if (porj(i, j, k) .eq. noflux) porflux = zero
    if (porj(i, j, k) .eq. boundflux) then
      porvel = zero
      vnp = sface
      vnm = sface
      call pushcontrol1b(0)
    else
      call pushcontrol1b(1)
    end if
! incorporate porflux in porvel.
    porvel = porvel*porflux
! compute the normal velocities for the face as well as the
! mass fluxes.
    qsp = (vnp-sface)*porvel
    qsm = (vnm-sface)*porvel
    rqsp = qsp*w(i, j+1, k, irho)
    rqsm = qsm*w(i, j, k, irho)
! compute the sum of the pressure multiplied by porflux.
! for the default value of porflux, 0.5, this leads to
! the average pressure.
    pa = porflux*(p(i, j+1, k)+p(i, j, k))
! compute the fluxes and scatter them to the cells
! i,j,k and i,j+1,k. store the density flux in the
! mass flow of the appropriate sliding mesh interface.
    fsd = dwd(i, j, k, irhoe) - dwd(i, j+1, k, irhoe)
    tempd0 = porflux*fsd
    qspd = w(i, j+1, k, irhoe)*fsd
    wd(i, j+1, k, irhoe) = wd(i, j+1, k, irhoe) + qsp*fsd
    qsmd = w(i, j, k, irhoe)*fsd
    wd(i, j, k, irhoe) = wd(i, j, k, irhoe) + qsm*fsd
    pd(i, j+1, k) = pd(i, j+1, k) + vnp*tempd0
    pd(i, j, k) = pd(i, j, k) + vnm*tempd0
    fsd = dwd(i, j, k, imz) - dwd(i, j+1, k, imz)
    rqspd = w(i, j+1, k, ivz)*fsd
    wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + rqsp*fsd
    rqsmd = w(i, j, k, ivz)*fsd
    wd(i, j, k, ivz) = wd(i, j, k, ivz) + rqsm*fsd
    pad = sj(i, j, k, 3)*fsd
    sjd(i, j, k, 3) = sjd(i, j, k, 3) + pa*fsd
    fsd = dwd(i, j, k, imy) - dwd(i, j+1, k, imy)
    rqspd = rqspd + w(i, j+1, k, ivy)*fsd
    wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + rqsp*fsd
    rqsmd = rqsmd + w(i, j, k, ivy)*fsd
    wd(i, j, k, ivy) = wd(i, j, k, ivy) + rqsm*fsd
    pad = pad + sj(i, j, k, 2)*fsd
    sjd(i, j, k, 2) = sjd(i, j, k, 2) + pa*fsd
    fsd = dwd(i, j, k, imx) - dwd(i, j+1, k, imx)
    rqspd = rqspd + w(i, j+1, k, ivx)*fsd
    wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + rqsp*fsd
    rqsmd = rqsmd + w(i, j, k, ivx)*fsd
    wd(i, j, k, ivx) = wd(i, j, k, ivx) + rqsm*fsd
    pad = pad + sj(i, j, k, 1)*fsd
    sjd(i, j, k, 1) = sjd(i, j, k, 1) + pa*fsd
    fsd = dwd(i, j, k, irho) - dwd(i, j+1, k, irho)
    rqspd = rqspd + fsd
    rqsmd = rqsmd + fsd
    pd(i, j+1, k) = pd(i, j+1, k) + porflux*pad
    pd(i, j, k) = pd(i, j, k) + porflux*pad
    qsmd = qsmd + w(i, j, k, irho)*rqsmd
    vnmd = porvel*qsmd + p(i, j, k)*tempd0
    wd(i, j, k, irho) = wd(i, j, k, irho) + qsm*rqsmd
    qspd = qspd + w(i, j+1, k, irho)*rqspd
    vnpd = porvel*qspd + p(i, j+1, k)*tempd0
    wd(i, j+1, k, irho) = wd(i, j+1, k, irho) + qsp*rqspd
    call popcontrol1b(branch)
    if (branch .eq. 0) then
      vnmd = 0.0_8
      vnpd = 0.0_8
    end if
    wd(i, j, k, ivx) = wd(i, j, k, ivx) + sj(i, j, k, 1)*vnmd
    sjd(i, j, k, 1) = sjd(i, j, k, 1) + w(i, j, k, ivx)*vnmd
    wd(i, j, k, ivy) = wd(i, j, k, ivy) + sj(i, j, k, 2)*vnmd
    sjd(i, j, k, 2) = sjd(i, j, k, 2) + w(i, j, k, ivy)*vnmd
    wd(i, j, k, ivz) = wd(i, j, k, ivz) + sj(i, j, k, 3)*vnmd
    sjd(i, j, k, 3) = sjd(i, j, k, 3) + w(i, j, k, ivz)*vnmd
    wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + sj(i, j, k, 1)*vnpd
    sjd(i, j, k, 1) = sjd(i, j, k, 1) + w(i, j+1, k, ivx)*vnpd
    wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + sj(i, j, k, 2)*vnpd
    sjd(i, j, k, 2) = sjd(i, j, k, 2) + w(i, j+1, k, ivy)*vnpd
    wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + sj(i, j, k, 3)*vnpd
    sjd(i, j, k, 3) = sjd(i, j, k, 3) + w(i, j+1, k, ivz)*vnpd
  end do
  call popreal8(qsp)
  call popreal8(qsm)
  call popreal8(porvel)
  call popreal8(porflux)
  call popreal8(vnp)
  call popreal8(vnm)
  call popinteger4(j)
  call popinteger4(i)
!
!      ******************************************************************
!      *                                                                *
!      * begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
! initialize sface to zero. this value will be used if the
! block is not moving.
  sface = zero
  do ii=0,il*ny*nz-1
    i = mod(ii, il) + 1
    j = mod(ii/il, ny) + 2
    k = ii/(il*ny) + 2
! set the dot product of the grid velocity and the
! normal in i-direction for a moving face.
    if (addgridvelocities) sface = sfacei(i, j, k)
! compute the normal velocities of the left and right state.
    vnp = w(i+1, j, k, ivx)*si(i, j, k, 1) + w(i+1, j, k, ivy)*si(i, j, &
&     k, 2) + w(i+1, j, k, ivz)*si(i, j, k, 3)
    vnm = w(i, j, k, ivx)*si(i, j, k, 1) + w(i, j, k, ivy)*si(i, j, k, 2&
&     ) + w(i, j, k, ivz)*si(i, j, k, 3)
! set the values of the porosities for this face.
! porvel defines the porosity w.r.t. velocity;
! porflux defines the porosity w.r.t. the entire flux.
! the latter is only zero for a discontinuous block
! boundary that must be treated conservatively.
! the default value of porflux is 0.5, such that the
! correct central flux is scattered to both cells.
! in case of a boundflux the normal velocity is set
! to sface.
    porvel = one
    porflux = half
    if (pori(i, j, k) .eq. noflux) porflux = zero
    if (pori(i, j, k) .eq. boundflux) then
      porvel = zero
      vnp = sface
      vnm = sface
      call pushcontrol1b(0)
    else
      call pushcontrol1b(1)
    end if
! incorporate porflux in porvel.
    porvel = porvel*porflux
! compute the normal velocities relative to the grid for
! the face as well as the mass fluxes.
    qsp = (vnp-sface)*porvel
    qsm = (vnm-sface)*porvel
    rqsp = qsp*w(i+1, j, k, irho)
    rqsm = qsm*w(i, j, k, irho)
! compute the sum of the pressure multiplied by porflux.
! for the default value of porflux, 0.5, this leads to
! the average pressure.
    pa = porflux*(p(i+1, j, k)+p(i, j, k))
! compute the fluxes and scatter them to the cells
! i,j,k and i+1,j,k. store the density flux in the
! mass flow of the appropriate sliding mesh interface.
    fsd = dwd(i, j, k, irhoe) - dwd(i+1, j, k, irhoe)
    tempd = porflux*fsd
    qspd = w(i+1, j, k, irhoe)*fsd
    wd(i+1, j, k, irhoe) = wd(i+1, j, k, irhoe) + qsp*fsd
    qsmd = w(i, j, k, irhoe)*fsd
    wd(i, j, k, irhoe) = wd(i, j, k, irhoe) + qsm*fsd
    pd(i+1, j, k) = pd(i+1, j, k) + vnp*tempd
    pd(i, j, k) = pd(i, j, k) + vnm*tempd
    fsd = dwd(i, j, k, imz) - dwd(i+1, j, k, imz)
    rqspd = w(i+1, j, k, ivz)*fsd
    wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + rqsp*fsd
    rqsmd = w(i, j, k, ivz)*fsd
    wd(i, j, k, ivz) = wd(i, j, k, ivz) + rqsm*fsd
    pad = si(i, j, k, 3)*fsd
    sid(i, j, k, 3) = sid(i, j, k, 3) + pa*fsd
    fsd = dwd(i, j, k, imy) - dwd(i+1, j, k, imy)
    rqspd = rqspd + w(i+1, j, k, ivy)*fsd
    wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + rqsp*fsd
    rqsmd = rqsmd + w(i, j, k, ivy)*fsd
    wd(i, j, k, ivy) = wd(i, j, k, ivy) + rqsm*fsd
    pad = pad + si(i, j, k, 2)*fsd
    sid(i, j, k, 2) = sid(i, j, k, 2) + pa*fsd
    fsd = dwd(i, j, k, imx) - dwd(i+1, j, k, imx)
    rqspd = rqspd + w(i+1, j, k, ivx)*fsd
    wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + rqsp*fsd
    rqsmd = rqsmd + w(i, j, k, ivx)*fsd
    wd(i, j, k, ivx) = wd(i, j, k, ivx) + rqsm*fsd
    pad = pad + si(i, j, k, 1)*fsd
    sid(i, j, k, 1) = sid(i, j, k, 1) + pa*fsd
    fsd = dwd(i, j, k, irho) - dwd(i+1, j, k, irho)
    rqspd = rqspd + fsd
    rqsmd = rqsmd + fsd
    pd(i+1, j, k) = pd(i+1, j, k) + porflux*pad
    pd(i, j, k) = pd(i, j, k) + porflux*pad
    qsmd = qsmd + w(i, j, k, irho)*rqsmd
    vnmd = porvel*qsmd + p(i, j, k)*tempd
    wd(i, j, k, irho) = wd(i, j, k, irho) + qsm*rqsmd
    qspd = qspd + w(i+1, j, k, irho)*rqspd
    vnpd = porvel*qspd + p(i+1, j, k)*tempd
    wd(i+1, j, k, irho) = wd(i+1, j, k, irho) + qsp*rqspd
    call popcontrol1b(branch)
    if (branch .eq. 0) then
      vnmd = 0.0_8
      vnpd = 0.0_8
    end if
    wd(i, j, k, ivx) = wd(i, j, k, ivx) + si(i, j, k, 1)*vnmd
    sid(i, j, k, 1) = sid(i, j, k, 1) + w(i, j, k, ivx)*vnmd
    wd(i, j, k, ivy) = wd(i, j, k, ivy) + si(i, j, k, 2)*vnmd
    sid(i, j, k, 2) = sid(i, j, k, 2) + w(i, j, k, ivy)*vnmd
    wd(i, j, k, ivz) = wd(i, j, k, ivz) + si(i, j, k, 3)*vnmd
    sid(i, j, k, 3) = sid(i, j, k, 3) + w(i, j, k, ivz)*vnmd
    wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + si(i, j, k, 1)*vnpd
    sid(i, j, k, 1) = sid(i, j, k, 1) + w(i+1, j, k, ivx)*vnpd
    wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + si(i, j, k, 2)*vnpd
    sid(i, j, k, 2) = sid(i, j, k, 2) + w(i+1, j, k, ivy)*vnpd
    wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + si(i, j, k, 3)*vnpd
    sid(i, j, k, 3) = sid(i, j, k, 3) + w(i+1, j, k, ivz)*vnpd
  end do
end subroutine inviscidcentralflux_b
