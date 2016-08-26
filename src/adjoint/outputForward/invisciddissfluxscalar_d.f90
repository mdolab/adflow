!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of invisciddissfluxscalar in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: *fw
!   with respect to varying inputs: gammainf rhoinf pinfcorr *p
!                *w *radi *radj *radk
!   plus diff mem management of: p:in w:in fw:in radi:in radj:in
!                radk:in
subroutine invisciddissfluxscalar_d()
!
!       invisciddissfluxscalar computes the scalar artificial          
!       dissipation, see aiaa paper 81-1259, for a given block.        
!       therefore it is assumed that the pointers in  blockpointers    
!       already point to the correct block.                            
!
  use constants
  use blockpointers, only : nx, ny, nz, il, jl, kl, ie, je, ke, ib, jb&
& , kb, w, wd, p, pd, pori, porj, pork, fw, fwd, radi, radid, radj, &
& radjd, radk, radkd, gamma
  use flowvarrefstate, only : gammainf, gammainfd, pinfcorr, pinfcorrd&
& , rhoinf, rhoinfd
  use inputdiscretization, only : vis2, vis4
  use inputphysics, only : equations
  use iteration, only : rfil
  use utils_d, only : mydim, mydim_d
  implicit none
!
!      local parameter.
!
  real(kind=realtype), parameter :: dssmax=0.25_realtype
!
!      local variables.
!
  integer(kind=inttype) :: i, j, k, ind, ii
  real(kind=realtype) :: sslim, rhoi
  real(kind=realtype) :: sslimd
  real(kind=realtype) :: sfil, fis2, fis4
  real(kind=realtype) :: ppor, rrad, dis2, dis4
  real(kind=realtype) :: rradd, dis2d, dis4d
  real(kind=realtype) :: ddw1, ddw2, ddw3, ddw4, ddw5, fs
  real(kind=realtype) :: ddw1d, ddw2d, ddw3d, ddw4d, ddw5d, fsd
  real(kind=realtype), dimension(ie, je, ke, 3) :: dss
  real(kind=realtype), dimension(ie, je, ke, 3) :: dssd
  real(kind=realtype), dimension(0:ib, 0:jb, 0:kb) :: ss
  real(kind=realtype), dimension(0:ib, 0:jb, 0:kb) :: ssd
  intrinsic abs
  intrinsic max
  intrinsic min
  real(kind=realtype) :: pwr1
  real(kind=realtype) :: pwr1d
  real(kind=realtype) :: min3
  real(kind=realtype) :: min2
  real(kind=realtype) :: min1
  real(kind=realtype) :: min1d
  real(kind=realtype) :: x3
  real(kind=realtype) :: x2
  real(kind=realtype) :: x2d
  real(kind=realtype) :: x1
  real(kind=realtype) :: y3d
  real(kind=realtype) :: x1d
  real(kind=realtype) :: min3d
  real(kind=realtype) :: y2d
  real(kind=realtype) :: abs0
  real(kind=realtype) :: min2d
  real(kind=realtype) :: y3
  real(kind=realtype) :: y2
  real(kind=realtype) :: x3d
  real(kind=realtype) :: y1
  real(kind=realtype) :: y1d
  if (rfil .ge. 0.) then
    abs0 = rfil
  else
    abs0 = -rfil
  end if
! check if rfil == 0. if so, the dissipative flux needs not to
! be computed.
  if (abs0 .lt. thresholdreal) then
    fwd = 0.0_8
    return
  else
! determine the variables used to compute the switch.
! for the inviscid case this is the pressure; for the viscous
! case it is the entropy.
    select case  (equations) 
    case (eulerequations) 
! inviscid case. pressure switch is based on the pressure.
! also set the value of sslim. to be fully consistent this
! must have the dimension of pressure and it is therefore
! set to a fraction of the free stream value.
      sslimd = 0.001_realtype*pinfcorrd
      sslim = 0.001_realtype*pinfcorr
! copy the pressure in ss. only need the entries used in the
! discretization, i.e. not including the corner halo's, but we'll
! just copy all anyway. 
      ssd = pd
      ss = p
      dssd = 0.0_8
    case (nsequations, ransequations) 
!===============================================================
! viscous case. pressure switch is based on the entropy.
! also set the value of sslim. to be fully consistent this
! must have the dimension of entropy and it is therefore
! set to a fraction of the free stream value.
      if (rhoinf .gt. 0.0_8) then
        pwr1d = rhoinf**gammainf*(log(rhoinf)*gammainfd+gammainf*rhoinfd&
&         /rhoinf)
      else if (rhoinf .eq. 0.0_8) then
        if (gammainf .eq. 1.0) then
          pwr1d = rhoinfd
        else
          pwr1d = 0.0_8
        end if
      else if (gammainf .eq. int(gammainf)) then
        pwr1d = gammainf*rhoinf**(gammainf-1)*rhoinfd
      else
        pwr1d = 0.0_8
      end if
      pwr1 = rhoinf**gammainf
      sslimd = (0.001_realtype*pinfcorrd*pwr1-0.001_realtype*pinfcorr*&
&       pwr1d)/pwr1**2
      sslim = 0.001_realtype*pinfcorr/pwr1
      ssd = 0.0_8
! store the entropy in ss. see above. 
      do k=0,kb
        do j=0,jb
          do i=0,ib
            if (w(i, j, k, irho) .gt. 0.0_8 .or. (w(i, j, k, irho) .lt. &
&               0.0_8 .and. gamma(i, j, k) .eq. int(gamma(i, j, k)))) &
&           then
              pwr1d = gamma(i, j, k)*w(i, j, k, irho)**(gamma(i, j, k)-1&
&               )*wd(i, j, k, irho)
            else if (w(i, j, k, irho) .eq. 0.0_8 .and. gamma(i, j, k) &
&               .eq. 1.0) then
              pwr1d = wd(i, j, k, irho)
            else
              pwr1d = 0.0_8
            end if
            pwr1 = w(i, j, k, irho)**gamma(i, j, k)
            ssd(i, j, k) = (pd(i, j, k)*pwr1-p(i, j, k)*pwr1d)/pwr1**2
            ss(i, j, k) = p(i, j, k)/pwr1
          end do
        end do
      end do
      dssd = 0.0_8
    case default
      sslimd = 0.0_8
      ssd = 0.0_8
      dssd = 0.0_8
    end select
! compute the pressure sensor for each cell, in each direction:
    do k=1,ke
      do j=1,je
        do i=1,ie
          x1d = ((ssd(i+1, j, k)-two*ssd(i, j, k)+ssd(i-1, j, k))*(ss(i+&
&           1, j, k)+two*ss(i, j, k)+ss(i-1, j, k)+sslim)-(ss(i+1, j, k)&
&           -two*ss(i, j, k)+ss(i-1, j, k))*(ssd(i+1, j, k)+two*ssd(i, j&
&           , k)+ssd(i-1, j, k)+sslimd))/(ss(i+1, j, k)+two*ss(i, j, k)+&
&           ss(i-1, j, k)+sslim)**2
          x1 = (ss(i+1, j, k)-two*ss(i, j, k)+ss(i-1, j, k))/(ss(i+1, j&
&           , k)+two*ss(i, j, k)+ss(i-1, j, k)+sslim)
          if (x1 .ge. 0.) then
            dssd(i, j, k, 1) = x1d
            dss(i, j, k, 1) = x1
          else
            dssd(i, j, k, 1) = -x1d
            dss(i, j, k, 1) = -x1
          end if
          x2d = ((ssd(i, j+1, k)-two*ssd(i, j, k)+ssd(i, j-1, k))*(ss(i&
&           , j+1, k)+two*ss(i, j, k)+ss(i, j-1, k)+sslim)-(ss(i, j+1, k&
&           )-two*ss(i, j, k)+ss(i, j-1, k))*(ssd(i, j+1, k)+two*ssd(i, &
&           j, k)+ssd(i, j-1, k)+sslimd))/(ss(i, j+1, k)+two*ss(i, j, k)&
&           +ss(i, j-1, k)+sslim)**2
          x2 = (ss(i, j+1, k)-two*ss(i, j, k)+ss(i, j-1, k))/(ss(i, j+1&
&           , k)+two*ss(i, j, k)+ss(i, j-1, k)+sslim)
          if (x2 .ge. 0.) then
            dssd(i, j, k, 2) = x2d
            dss(i, j, k, 2) = x2
          else
            dssd(i, j, k, 2) = -x2d
            dss(i, j, k, 2) = -x2
          end if
          x3d = ((ssd(i, j, k+1)-two*ssd(i, j, k)+ssd(i, j, k-1))*(ss(i&
&           , j, k+1)+two*ss(i, j, k)+ss(i, j, k-1)+sslim)-(ss(i, j, k+1&
&           )-two*ss(i, j, k)+ss(i, j, k-1))*(ssd(i, j, k+1)+two*ssd(i, &
&           j, k)+ssd(i, j, k-1)+sslimd))/(ss(i, j, k+1)+two*ss(i, j, k)&
&           +ss(i, j, k-1)+sslim)**2
          x3 = (ss(i, j, k+1)-two*ss(i, j, k)+ss(i, j, k-1))/(ss(i, j, k&
&           +1)+two*ss(i, j, k)+ss(i, j, k-1)+sslim)
          if (x3 .ge. 0.) then
            dssd(i, j, k, 3) = x3d
            dss(i, j, k, 3) = x3
          else
            dssd(i, j, k, 3) = -x3d
            dss(i, j, k, 3) = -x3
          end if
        end do
      end do
    end do
! set a couple of constants for the scheme.
    fis2 = rfil*vis2
    fis4 = rfil*vis4
    sfil = one - rfil
! initialize the dissipative residual to a certain times,
! possibly zero, the previously stored value. owned cells
! only, because the halo values do not matter.
    fwd = 0.0_8
    fw = sfil*fw
    fwd = 0.0_8
!
!       dissipative fluxes in the i-direction.                         
!
    do k=2,kl
      do j=2,jl
        do i=1,il
! compute the dissipation coefficients for this face.
          ppor = zero
          if (pori(i, j, k) .eq. normalflux) ppor = half
          rradd = ppor*(radid(i, j, k)+radid(i+1, j, k))
          rrad = ppor*(radi(i, j, k)+radi(i+1, j, k))
          if (dss(i, j, k, 1) .lt. dss(i+1, j, k, 1)) then
            y1d = dssd(i+1, j, k, 1)
            y1 = dss(i+1, j, k, 1)
          else
            y1d = dssd(i, j, k, 1)
            y1 = dss(i, j, k, 1)
          end if
          if (dssmax .gt. y1) then
            min1d = y1d
            min1 = y1
          else
            min1 = dssmax
            min1d = 0.0_8
          end if
          dis2d = fis2*(rradd*min1+rrad*min1d)
          dis2 = fis2*rrad*min1
          dis4d = mydim_d(fis4*rrad, fis4*rradd, dis2, dis2d, dis4)
! compute and scatter the dissipative flux.
! density. store it in the mass flow of the
! appropriate sliding mesh interface.
          ddw1d = wd(i+1, j, k, irho) - wd(i, j, k, irho)
          ddw1 = w(i+1, j, k, irho) - w(i, j, k, irho)
          fsd = dis2d*ddw1 + dis2*ddw1d - dis4d*(w(i+2, j, k, irho)-w(i-&
&           1, j, k, irho)-three*ddw1) - dis4*(wd(i+2, j, k, irho)-wd(i-&
&           1, j, k, irho)-three*ddw1d)
          fs = dis2*ddw1 - dis4*(w(i+2, j, k, irho)-w(i-1, j, k, irho)-&
&           three*ddw1)
          fwd(i+1, j, k, irho) = fwd(i+1, j, k, irho) + fsd
          fw(i+1, j, k, irho) = fw(i+1, j, k, irho) + fs
          fwd(i, j, k, irho) = fwd(i, j, k, irho) - fsd
          fw(i, j, k, irho) = fw(i, j, k, irho) - fs
! x-momentum.
          ddw2d = wd(i+1, j, k, ivx)*w(i+1, j, k, irho) + w(i+1, j, k, &
&           ivx)*wd(i+1, j, k, irho) - wd(i, j, k, ivx)*w(i, j, k, irho)&
&           - w(i, j, k, ivx)*wd(i, j, k, irho)
          ddw2 = w(i+1, j, k, ivx)*w(i+1, j, k, irho) - w(i, j, k, ivx)*&
&           w(i, j, k, irho)
          fsd = dis2d*ddw2 + dis2*ddw2d - dis4d*(w(i+2, j, k, ivx)*w(i+2&
&           , j, k, irho)-w(i-1, j, k, ivx)*w(i-1, j, k, irho)-three*&
&           ddw2) - dis4*(wd(i+2, j, k, ivx)*w(i+2, j, k, irho)+w(i+2, j&
&           , k, ivx)*wd(i+2, j, k, irho)-wd(i-1, j, k, ivx)*w(i-1, j, k&
&           , irho)-w(i-1, j, k, ivx)*wd(i-1, j, k, irho)-three*ddw2d)
          fs = dis2*ddw2 - dis4*(w(i+2, j, k, ivx)*w(i+2, j, k, irho)-w(&
&           i-1, j, k, ivx)*w(i-1, j, k, irho)-three*ddw2)
          fwd(i+1, j, k, imx) = fwd(i+1, j, k, imx) + fsd
          fw(i+1, j, k, imx) = fw(i+1, j, k, imx) + fs
          fwd(i, j, k, imx) = fwd(i, j, k, imx) - fsd
          fw(i, j, k, imx) = fw(i, j, k, imx) - fs
! y-momentum.
          ddw3d = wd(i+1, j, k, ivy)*w(i+1, j, k, irho) + w(i+1, j, k, &
&           ivy)*wd(i+1, j, k, irho) - wd(i, j, k, ivy)*w(i, j, k, irho)&
&           - w(i, j, k, ivy)*wd(i, j, k, irho)
          ddw3 = w(i+1, j, k, ivy)*w(i+1, j, k, irho) - w(i, j, k, ivy)*&
&           w(i, j, k, irho)
          fsd = dis2d*ddw3 + dis2*ddw3d - dis4d*(w(i+2, j, k, ivy)*w(i+2&
&           , j, k, irho)-w(i-1, j, k, ivy)*w(i-1, j, k, irho)-three*&
&           ddw3) - dis4*(wd(i+2, j, k, ivy)*w(i+2, j, k, irho)+w(i+2, j&
&           , k, ivy)*wd(i+2, j, k, irho)-wd(i-1, j, k, ivy)*w(i-1, j, k&
&           , irho)-w(i-1, j, k, ivy)*wd(i-1, j, k, irho)-three*ddw3d)
          fs = dis2*ddw3 - dis4*(w(i+2, j, k, ivy)*w(i+2, j, k, irho)-w(&
&           i-1, j, k, ivy)*w(i-1, j, k, irho)-three*ddw3)
          fwd(i+1, j, k, imy) = fwd(i+1, j, k, imy) + fsd
          fw(i+1, j, k, imy) = fw(i+1, j, k, imy) + fs
          fwd(i, j, k, imy) = fwd(i, j, k, imy) - fsd
          fw(i, j, k, imy) = fw(i, j, k, imy) - fs
! z-momentum.
          ddw4d = wd(i+1, j, k, ivz)*w(i+1, j, k, irho) + w(i+1, j, k, &
&           ivz)*wd(i+1, j, k, irho) - wd(i, j, k, ivz)*w(i, j, k, irho)&
&           - w(i, j, k, ivz)*wd(i, j, k, irho)
          ddw4 = w(i+1, j, k, ivz)*w(i+1, j, k, irho) - w(i, j, k, ivz)*&
&           w(i, j, k, irho)
          fsd = dis2d*ddw4 + dis2*ddw4d - dis4d*(w(i+2, j, k, ivz)*w(i+2&
&           , j, k, irho)-w(i-1, j, k, ivz)*w(i-1, j, k, irho)-three*&
&           ddw4) - dis4*(wd(i+2, j, k, ivz)*w(i+2, j, k, irho)+w(i+2, j&
&           , k, ivz)*wd(i+2, j, k, irho)-wd(i-1, j, k, ivz)*w(i-1, j, k&
&           , irho)-w(i-1, j, k, ivz)*wd(i-1, j, k, irho)-three*ddw4d)
          fs = dis2*ddw4 - dis4*(w(i+2, j, k, ivz)*w(i+2, j, k, irho)-w(&
&           i-1, j, k, ivz)*w(i-1, j, k, irho)-three*ddw4)
          fwd(i+1, j, k, imz) = fwd(i+1, j, k, imz) + fsd
          fw(i+1, j, k, imz) = fw(i+1, j, k, imz) + fs
          fwd(i, j, k, imz) = fwd(i, j, k, imz) - fsd
          fw(i, j, k, imz) = fw(i, j, k, imz) - fs
! energy.
          ddw5d = wd(i+1, j, k, irhoe) + pd(i+1, j, k) - wd(i, j, k, &
&           irhoe) - pd(i, j, k)
          ddw5 = w(i+1, j, k, irhoe) + p(i+1, j, k) - (w(i, j, k, irhoe)&
&           +p(i, j, k))
          fsd = dis2d*ddw5 + dis2*ddw5d - dis4d*(w(i+2, j, k, irhoe)+p(i&
&           +2, j, k)-(w(i-1, j, k, irhoe)+p(i-1, j, k))-three*ddw5) - &
&           dis4*(wd(i+2, j, k, irhoe)+pd(i+2, j, k)-wd(i-1, j, k, irhoe&
&           )-pd(i-1, j, k)-three*ddw5d)
          fs = dis2*ddw5 - dis4*(w(i+2, j, k, irhoe)+p(i+2, j, k)-(w(i-1&
&           , j, k, irhoe)+p(i-1, j, k))-three*ddw5)
          fwd(i+1, j, k, irhoe) = fwd(i+1, j, k, irhoe) + fsd
          fw(i+1, j, k, irhoe) = fw(i+1, j, k, irhoe) + fs
          fwd(i, j, k, irhoe) = fwd(i, j, k, irhoe) - fsd
          fw(i, j, k, irhoe) = fw(i, j, k, irhoe) - fs
        end do
      end do
    end do
!
!       dissipative fluxes in the j-direction.                         
!
    do k=2,kl
      do j=1,jl
        do i=2,il
! compute the dissipation coefficients for this face.
          ppor = zero
          if (porj(i, j, k) .eq. normalflux) ppor = half
          rradd = ppor*(radjd(i, j, k)+radjd(i, j+1, k))
          rrad = ppor*(radj(i, j, k)+radj(i, j+1, k))
          if (dss(i, j, k, 2) .lt. dss(i, j+1, k, 2)) then
            y2d = dssd(i, j+1, k, 2)
            y2 = dss(i, j+1, k, 2)
          else
            y2d = dssd(i, j, k, 2)
            y2 = dss(i, j, k, 2)
          end if
          if (dssmax .gt. y2) then
            min2d = y2d
            min2 = y2
          else
            min2 = dssmax
            min2d = 0.0_8
          end if
          dis2d = fis2*(rradd*min2+rrad*min2d)
          dis2 = fis2*rrad*min2
          dis4d = mydim_d(fis4*rrad, fis4*rradd, dis2, dis2d, dis4)
! compute and scatter the dissipative flux.
! density. store it in the mass flow of the
! appropriate sliding mesh interface.
          ddw1d = wd(i, j+1, k, irho) - wd(i, j, k, irho)
          ddw1 = w(i, j+1, k, irho) - w(i, j, k, irho)
          fsd = dis2d*ddw1 + dis2*ddw1d - dis4d*(w(i, j+2, k, irho)-w(i&
&           , j-1, k, irho)-three*ddw1) - dis4*(wd(i, j+2, k, irho)-wd(i&
&           , j-1, k, irho)-three*ddw1d)
          fs = dis2*ddw1 - dis4*(w(i, j+2, k, irho)-w(i, j-1, k, irho)-&
&           three*ddw1)
          fwd(i, j+1, k, irho) = fwd(i, j+1, k, irho) + fsd
          fw(i, j+1, k, irho) = fw(i, j+1, k, irho) + fs
          fwd(i, j, k, irho) = fwd(i, j, k, irho) - fsd
          fw(i, j, k, irho) = fw(i, j, k, irho) - fs
! x-momentum.
          ddw2d = wd(i, j+1, k, ivx)*w(i, j+1, k, irho) + w(i, j+1, k, &
&           ivx)*wd(i, j+1, k, irho) - wd(i, j, k, ivx)*w(i, j, k, irho)&
&           - w(i, j, k, ivx)*wd(i, j, k, irho)
          ddw2 = w(i, j+1, k, ivx)*w(i, j+1, k, irho) - w(i, j, k, ivx)*&
&           w(i, j, k, irho)
          fsd = dis2d*ddw2 + dis2*ddw2d - dis4d*(w(i, j+2, k, ivx)*w(i, &
&           j+2, k, irho)-w(i, j-1, k, ivx)*w(i, j-1, k, irho)-three*&
&           ddw2) - dis4*(wd(i, j+2, k, ivx)*w(i, j+2, k, irho)+w(i, j+2&
&           , k, ivx)*wd(i, j+2, k, irho)-wd(i, j-1, k, ivx)*w(i, j-1, k&
&           , irho)-w(i, j-1, k, ivx)*wd(i, j-1, k, irho)-three*ddw2d)
          fs = dis2*ddw2 - dis4*(w(i, j+2, k, ivx)*w(i, j+2, k, irho)-w(&
&           i, j-1, k, ivx)*w(i, j-1, k, irho)-three*ddw2)
          fwd(i, j+1, k, imx) = fwd(i, j+1, k, imx) + fsd
          fw(i, j+1, k, imx) = fw(i, j+1, k, imx) + fs
          fwd(i, j, k, imx) = fwd(i, j, k, imx) - fsd
          fw(i, j, k, imx) = fw(i, j, k, imx) - fs
! y-momentum.
          ddw3d = wd(i, j+1, k, ivy)*w(i, j+1, k, irho) + w(i, j+1, k, &
&           ivy)*wd(i, j+1, k, irho) - wd(i, j, k, ivy)*w(i, j, k, irho)&
&           - w(i, j, k, ivy)*wd(i, j, k, irho)
          ddw3 = w(i, j+1, k, ivy)*w(i, j+1, k, irho) - w(i, j, k, ivy)*&
&           w(i, j, k, irho)
          fsd = dis2d*ddw3 + dis2*ddw3d - dis4d*(w(i, j+2, k, ivy)*w(i, &
&           j+2, k, irho)-w(i, j-1, k, ivy)*w(i, j-1, k, irho)-three*&
&           ddw3) - dis4*(wd(i, j+2, k, ivy)*w(i, j+2, k, irho)+w(i, j+2&
&           , k, ivy)*wd(i, j+2, k, irho)-wd(i, j-1, k, ivy)*w(i, j-1, k&
&           , irho)-w(i, j-1, k, ivy)*wd(i, j-1, k, irho)-three*ddw3d)
          fs = dis2*ddw3 - dis4*(w(i, j+2, k, ivy)*w(i, j+2, k, irho)-w(&
&           i, j-1, k, ivy)*w(i, j-1, k, irho)-three*ddw3)
          fwd(i, j+1, k, imy) = fwd(i, j+1, k, imy) + fsd
          fw(i, j+1, k, imy) = fw(i, j+1, k, imy) + fs
          fwd(i, j, k, imy) = fwd(i, j, k, imy) - fsd
          fw(i, j, k, imy) = fw(i, j, k, imy) - fs
! z-momentum.
          ddw4d = wd(i, j+1, k, ivz)*w(i, j+1, k, irho) + w(i, j+1, k, &
&           ivz)*wd(i, j+1, k, irho) - wd(i, j, k, ivz)*w(i, j, k, irho)&
&           - w(i, j, k, ivz)*wd(i, j, k, irho)
          ddw4 = w(i, j+1, k, ivz)*w(i, j+1, k, irho) - w(i, j, k, ivz)*&
&           w(i, j, k, irho)
          fsd = dis2d*ddw4 + dis2*ddw4d - dis4d*(w(i, j+2, k, ivz)*w(i, &
&           j+2, k, irho)-w(i, j-1, k, ivz)*w(i, j-1, k, irho)-three*&
&           ddw4) - dis4*(wd(i, j+2, k, ivz)*w(i, j+2, k, irho)+w(i, j+2&
&           , k, ivz)*wd(i, j+2, k, irho)-wd(i, j-1, k, ivz)*w(i, j-1, k&
&           , irho)-w(i, j-1, k, ivz)*wd(i, j-1, k, irho)-three*ddw4d)
          fs = dis2*ddw4 - dis4*(w(i, j+2, k, ivz)*w(i, j+2, k, irho)-w(&
&           i, j-1, k, ivz)*w(i, j-1, k, irho)-three*ddw4)
          fwd(i, j+1, k, imz) = fwd(i, j+1, k, imz) + fsd
          fw(i, j+1, k, imz) = fw(i, j+1, k, imz) + fs
          fwd(i, j, k, imz) = fwd(i, j, k, imz) - fsd
          fw(i, j, k, imz) = fw(i, j, k, imz) - fs
! energy.
          ddw5d = wd(i, j+1, k, irhoe) + pd(i, j+1, k) - wd(i, j, k, &
&           irhoe) - pd(i, j, k)
          ddw5 = w(i, j+1, k, irhoe) + p(i, j+1, k) - (w(i, j, k, irhoe)&
&           +p(i, j, k))
          fsd = dis2d*ddw5 + dis2*ddw5d - dis4d*(w(i, j+2, k, irhoe)+p(i&
&           , j+2, k)-(w(i, j-1, k, irhoe)+p(i, j-1, k))-three*ddw5) - &
&           dis4*(wd(i, j+2, k, irhoe)+pd(i, j+2, k)-wd(i, j-1, k, irhoe&
&           )-pd(i, j-1, k)-three*ddw5d)
          fs = dis2*ddw5 - dis4*(w(i, j+2, k, irhoe)+p(i, j+2, k)-(w(i, &
&           j-1, k, irhoe)+p(i, j-1, k))-three*ddw5)
          fwd(i, j+1, k, irhoe) = fwd(i, j+1, k, irhoe) + fsd
          fw(i, j+1, k, irhoe) = fw(i, j+1, k, irhoe) + fs
          fwd(i, j, k, irhoe) = fwd(i, j, k, irhoe) - fsd
          fw(i, j, k, irhoe) = fw(i, j, k, irhoe) - fs
        end do
      end do
    end do
!
!       dissipative fluxes in the k-direction.                         
!
    do k=1,kl
      do j=2,jl
        do i=2,il
! compute the dissipation coefficients for this face.
          ppor = zero
          if (pork(i, j, k) .eq. normalflux) ppor = half
          rradd = ppor*(radkd(i, j, k)+radkd(i, j, k+1))
          rrad = ppor*(radk(i, j, k)+radk(i, j, k+1))
          if (dss(i, j, k, 3) .lt. dss(i, j, k+1, 3)) then
            y3d = dssd(i, j, k+1, 3)
            y3 = dss(i, j, k+1, 3)
          else
            y3d = dssd(i, j, k, 3)
            y3 = dss(i, j, k, 3)
          end if
          if (dssmax .gt. y3) then
            min3d = y3d
            min3 = y3
          else
            min3 = dssmax
            min3d = 0.0_8
          end if
          dis2d = fis2*(rradd*min3+rrad*min3d)
          dis2 = fis2*rrad*min3
          dis4d = mydim_d(fis4*rrad, fis4*rradd, dis2, dis2d, dis4)
! compute and scatter the dissipative flux.
! density. store it in the mass flow of the
! appropriate sliding mesh interface.
          ddw1d = wd(i, j, k+1, irho) - wd(i, j, k, irho)
          ddw1 = w(i, j, k+1, irho) - w(i, j, k, irho)
          fsd = dis2d*ddw1 + dis2*ddw1d - dis4d*(w(i, j, k+2, irho)-w(i&
&           , j, k-1, irho)-three*ddw1) - dis4*(wd(i, j, k+2, irho)-wd(i&
&           , j, k-1, irho)-three*ddw1d)
          fs = dis2*ddw1 - dis4*(w(i, j, k+2, irho)-w(i, j, k-1, irho)-&
&           three*ddw1)
          fwd(i, j, k+1, irho) = fwd(i, j, k+1, irho) + fsd
          fw(i, j, k+1, irho) = fw(i, j, k+1, irho) + fs
          fwd(i, j, k, irho) = fwd(i, j, k, irho) - fsd
          fw(i, j, k, irho) = fw(i, j, k, irho) - fs
! x-momentum.
          ddw2d = wd(i, j, k+1, ivx)*w(i, j, k+1, irho) + w(i, j, k+1, &
&           ivx)*wd(i, j, k+1, irho) - wd(i, j, k, ivx)*w(i, j, k, irho)&
&           - w(i, j, k, ivx)*wd(i, j, k, irho)
          ddw2 = w(i, j, k+1, ivx)*w(i, j, k+1, irho) - w(i, j, k, ivx)*&
&           w(i, j, k, irho)
          fsd = dis2d*ddw2 + dis2*ddw2d - dis4d*(w(i, j, k+2, ivx)*w(i, &
&           j, k+2, irho)-w(i, j, k-1, ivx)*w(i, j, k-1, irho)-three*&
&           ddw2) - dis4*(wd(i, j, k+2, ivx)*w(i, j, k+2, irho)+w(i, j, &
&           k+2, ivx)*wd(i, j, k+2, irho)-wd(i, j, k-1, ivx)*w(i, j, k-1&
&           , irho)-w(i, j, k-1, ivx)*wd(i, j, k-1, irho)-three*ddw2d)
          fs = dis2*ddw2 - dis4*(w(i, j, k+2, ivx)*w(i, j, k+2, irho)-w(&
&           i, j, k-1, ivx)*w(i, j, k-1, irho)-three*ddw2)
          fwd(i, j, k+1, imx) = fwd(i, j, k+1, imx) + fsd
          fw(i, j, k+1, imx) = fw(i, j, k+1, imx) + fs
          fwd(i, j, k, imx) = fwd(i, j, k, imx) - fsd
          fw(i, j, k, imx) = fw(i, j, k, imx) - fs
! y-momentum.
          ddw3d = wd(i, j, k+1, ivy)*w(i, j, k+1, irho) + w(i, j, k+1, &
&           ivy)*wd(i, j, k+1, irho) - wd(i, j, k, ivy)*w(i, j, k, irho)&
&           - w(i, j, k, ivy)*wd(i, j, k, irho)
          ddw3 = w(i, j, k+1, ivy)*w(i, j, k+1, irho) - w(i, j, k, ivy)*&
&           w(i, j, k, irho)
          fsd = dis2d*ddw3 + dis2*ddw3d - dis4d*(w(i, j, k+2, ivy)*w(i, &
&           j, k+2, irho)-w(i, j, k-1, ivy)*w(i, j, k-1, irho)-three*&
&           ddw3) - dis4*(wd(i, j, k+2, ivy)*w(i, j, k+2, irho)+w(i, j, &
&           k+2, ivy)*wd(i, j, k+2, irho)-wd(i, j, k-1, ivy)*w(i, j, k-1&
&           , irho)-w(i, j, k-1, ivy)*wd(i, j, k-1, irho)-three*ddw3d)
          fs = dis2*ddw3 - dis4*(w(i, j, k+2, ivy)*w(i, j, k+2, irho)-w(&
&           i, j, k-1, ivy)*w(i, j, k-1, irho)-three*ddw3)
          fwd(i, j, k+1, imy) = fwd(i, j, k+1, imy) + fsd
          fw(i, j, k+1, imy) = fw(i, j, k+1, imy) + fs
          fwd(i, j, k, imy) = fwd(i, j, k, imy) - fsd
          fw(i, j, k, imy) = fw(i, j, k, imy) - fs
! z-momentum.
          ddw4d = wd(i, j, k+1, ivz)*w(i, j, k+1, irho) + w(i, j, k+1, &
&           ivz)*wd(i, j, k+1, irho) - wd(i, j, k, ivz)*w(i, j, k, irho)&
&           - w(i, j, k, ivz)*wd(i, j, k, irho)
          ddw4 = w(i, j, k+1, ivz)*w(i, j, k+1, irho) - w(i, j, k, ivz)*&
&           w(i, j, k, irho)
          fsd = dis2d*ddw4 + dis2*ddw4d - dis4d*(w(i, j, k+2, ivz)*w(i, &
&           j, k+2, irho)-w(i, j, k-1, ivz)*w(i, j, k-1, irho)-three*&
&           ddw4) - dis4*(wd(i, j, k+2, ivz)*w(i, j, k+2, irho)+w(i, j, &
&           k+2, ivz)*wd(i, j, k+2, irho)-wd(i, j, k-1, ivz)*w(i, j, k-1&
&           , irho)-w(i, j, k-1, ivz)*wd(i, j, k-1, irho)-three*ddw4d)
          fs = dis2*ddw4 - dis4*(w(i, j, k+2, ivz)*w(i, j, k+2, irho)-w(&
&           i, j, k-1, ivz)*w(i, j, k-1, irho)-three*ddw4)
          fwd(i, j, k+1, imz) = fwd(i, j, k+1, imz) + fsd
          fw(i, j, k+1, imz) = fw(i, j, k+1, imz) + fs
          fwd(i, j, k, imz) = fwd(i, j, k, imz) - fsd
          fw(i, j, k, imz) = fw(i, j, k, imz) - fs
! energy.
          ddw5d = wd(i, j, k+1, irhoe) + pd(i, j, k+1) - wd(i, j, k, &
&           irhoe) - pd(i, j, k)
          ddw5 = w(i, j, k+1, irhoe) + p(i, j, k+1) - (w(i, j, k, irhoe)&
&           +p(i, j, k))
          fsd = dis2d*ddw5 + dis2*ddw5d - dis4d*(w(i, j, k+2, irhoe)+p(i&
&           , j, k+2)-(w(i, j, k-1, irhoe)+p(i, j, k-1))-three*ddw5) - &
&           dis4*(wd(i, j, k+2, irhoe)+pd(i, j, k+2)-wd(i, j, k-1, irhoe&
&           )-pd(i, j, k-1)-three*ddw5d)
          fs = dis2*ddw5 - dis4*(w(i, j, k+2, irhoe)+p(i, j, k+2)-(w(i, &
&           j, k-1, irhoe)+p(i, j, k-1))-three*ddw5)
          fwd(i, j, k+1, irhoe) = fwd(i, j, k+1, irhoe) + fsd
          fw(i, j, k+1, irhoe) = fw(i, j, k+1, irhoe) + fs
          fwd(i, j, k, irhoe) = fwd(i, j, k, irhoe) - fsd
          fw(i, j, k, irhoe) = fw(i, j, k, irhoe) - fs
        end do
      end do
    end do
  end if
end subroutine invisciddissfluxscalar_d
