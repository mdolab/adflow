   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
   !
   !  Differentiation of utauwf in forward (tangent) mode:
   !   variations   of useful results: *fw
   !   with respect to varying inputs: *w *fw *(*viscsubface.tau)
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          utauWF.f90                                      *
   !      * Author:        Georgi Kalitzin, Edwin van der Weide            *
   !      * Starting date: 10-01-2003                                      *
   !      * Last modified: 06-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE UTAUWF_SPATIAL_D(rfilv)
   USE BCTYPES
   USE INPUTPHYSICS
   USE BLOCKPOINTERS_D
   USE FLOWVARREFSTATE
   IMPLICIT NONE
   !
   !      ******************************************************************
   !      *                                                                *
   !      * utauWF substitutes the wall shear stress with values from a    *
   !      * look-up table, if desired.                                     *
   !      *                                                                *
   !      ******************************************************************
   !
   !
   !      Subroutine argument.
   !
   REAL(kind=realtype), INTENT(IN) :: rfilv
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i, j, nn
   REAL(kind=realtype) :: fact
   REAL(kind=realtype) :: tauxx, tauyy, tauzz
   REAL(kind=realtype) :: tauxxd, tauyyd, tauzzd
   REAL(kind=realtype) :: tauxy, tauxz, tauyz
   REAL(kind=realtype) :: tauxyd, tauxzd, tauyzd
   REAL(kind=realtype) :: rbar, ubar, vbar, wbar, vx, vy, vz
   REAL(kind=realtype) :: rbard, ubard, vbard, wbard, vxd, vyd, vzd
   REAL(kind=realtype) :: fmx, fmy, fmz, frhoe
   REAL(kind=realtype) :: fmxd, fmyd, fmzd, frhoed
   REAL(kind=realtype) :: veln, velnx, velny, velnz, tx, ty, tz
   REAL(kind=realtype) :: velnd, velnxd, velnyd, velnzd, txd, tyd, tzd
   REAL(kind=realtype) :: veltx, velty, veltz, veltmag
   REAL(kind=realtype) :: veltxd, veltyd, veltzd, veltmagd
   REAL(kind=realtype) :: txnx, txny, txnz, tynx, tyny, tynz
   REAL(kind=realtype) :: txnxd, txnyd, txnzd, tynxd, tynyd, tynzd
   REAL(kind=realtype) :: tznx, tzny, tznz
   REAL(kind=realtype) :: tznxd, tznyd, tznzd
   REAL(kind=realtype) :: tautn, tauwall, utau, re
   REAL(kind=realtype) :: tautnd, tauwalld, utaud, red
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1d, ww2d
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ss, rres
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: rresd
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: norm
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rrlv2, dd2wall2
   !
   !      Function definition.
   !
   REAL(kind=realtype) :: CURVEUPRE
   REAL(kind=realtype) :: CURVEUPRE_SPATIAL_D
   REAL(kind=realtype) :: arg1
   REAL(kind=realtype) :: arg1d
   INTRINSIC MAX
   REAL(kind=realtype) :: x1
   REAL(kind=realtype) :: max1d
   REAL(kind=realtype) :: x1d
   INTRINSIC SQRT
   REAL(kind=realtype) :: max1
   REAL(kind=realtype) :: y1
   REAL(kind=realtype) :: y1d
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Return immediately if no wall functions must be used.
   IF (.NOT.wallfunctions) THEN
   RETURN
   ELSE
   ! Loop over the viscous subfaces of this block.
   viscsubfaces:DO nn=1,nviscbocos
   ! Set a bunch of variables depending on the face id to make
   ! a generic treatment possible.
   SELECT CASE  (bcfaceid(nn)) 
   CASE (imin) 
   fact = -one
   ss => si(1, :, :, :)
   rresd => fwd(2, 1:, 1:, :)
   rres => fw(2, 1:, 1:, :)
   ww2d => wd(2, 1:, 1:, :)
   ww2 => w(2, 1:, 1:, :)
   ww1d => wd(1, 1:, 1:, :)
   ww1 => w(1, 1:, 1:, :)
   dd2wall2 => d2wall(2, :, :)
   rrlv2 => rlv(2, 1:, 1:)
   CASE (imax) 
   !===========================================================
   fact = one
   ss => si(il, :, :, :)
   rresd => fwd(il, 1:, 1:, :)
   rres => fw(il, 1:, 1:, :)
   ww2d => wd(il, 1:, 1:, :)
   ww2 => w(il, 1:, 1:, :)
   ww1d => wd(ie, 1:, 1:, :)
   ww1 => w(ie, 1:, 1:, :)
   dd2wall2 => d2wall(il, :, :)
   rrlv2 => rlv(il, 1:, 1:)
   CASE (jmin) 
   !===========================================================
   fact = -one
   ss => sj(:, 1, :, :)
   rresd => fwd(1:, 2, 1:, :)
   rres => fw(1:, 2, 1:, :)
   ww2d => wd(1:, 2, 1:, :)
   ww2 => w(1:, 2, 1:, :)
   ww1d => wd(1:, 1, 1:, :)
   ww1 => w(1:, 1, 1:, :)
   dd2wall2 => d2wall(:, 2, :)
   rrlv2 => rlv(1:, 2, 1:)
   CASE (jmax) 
   !===========================================================
   fact = one
   ss => sj(:, jl, :, :)
   rresd => fwd(1:, jl, 1:, :)
   rres => fw(1:, jl, 1:, :)
   ww2d => wd(1:, jl, 1:, :)
   ww2 => w(1:, jl, 1:, :)
   ww1d => wd(1:, je, 1:, :)
   ww1 => w(1:, je, 1:, :)
   dd2wall2 => d2wall(:, jl, :)
   rrlv2 => rlv(1:, jl, 1:)
   CASE (kmin) 
   !===========================================================
   fact = -one
   ss => sk(:, :, 1, :)
   rresd => fwd(1:, 1:, 2, :)
   rres => fw(1:, 1:, 2, :)
   ww2d => wd(1:, 1:, 2, :)
   ww2 => w(1:, 1:, 2, :)
   ww1d => wd(1:, 1:, 1, :)
   ww1 => w(1:, 1:, 1, :)
   dd2wall2 => d2wall(:, :, 2)
   rrlv2 => rlv(1:, 1:, 2)
   CASE (kmax) 
   !===========================================================
   fact = one
   ss => sk(:, :, kl, :)
   rresd => fwd(1:, 1:, kl, :)
   rres => fw(1:, 1:, kl, :)
   ww2d => wd(1:, 1:, kl, :)
   ww2 => w(1:, 1:, kl, :)
   ww1d => wd(1:, 1:, ke, :)
   ww1 => w(1:, 1:, ke, :)
   dd2wall2 => d2wall(:, :, kl)
   rrlv2 => rlv(1:, 1:, kl)
   END SELECT
   ! Set the pointer for the unit outward normals.
   norm => bcdata(nn)%norm
   ! Loop over the quadrilateral faces of the subface. Note
   ! that the nodal range of BCData must be used and not the
   ! cell range, because the latter may include the halo's in i
   ! and j-direction. The offset +1 is there, because inBeg and
   ! jnBeg refer to nodal ranges and not to cell ranges.
   DO j=bcdata(nn)%jnbeg+1,bcdata(nn)%jnend
   DO i=bcdata(nn)%inbeg+1,bcdata(nn)%inend
   ! Store the viscous stress tensor a bit easier.
   tauxxd = viscsubfaced(nn)%tau(i, j, 1)
   tauxx = viscsubface(nn)%tau(i, j, 1)
   tauyyd = viscsubfaced(nn)%tau(i, j, 2)
   tauyy = viscsubface(nn)%tau(i, j, 2)
   tauzzd = viscsubfaced(nn)%tau(i, j, 3)
   tauzz = viscsubface(nn)%tau(i, j, 3)
   tauxyd = viscsubfaced(nn)%tau(i, j, 4)
   tauxy = viscsubface(nn)%tau(i, j, 4)
   tauxzd = viscsubfaced(nn)%tau(i, j, 5)
   tauxz = viscsubface(nn)%tau(i, j, 5)
   tauyzd = viscsubfaced(nn)%tau(i, j, 6)
   tauyz = viscsubface(nn)%tau(i, j, 6)
   ! Compute the velocities at the wall face; these are only
   ! non-zero for moving a block. Also compute the density,
   ! which is needed to compute the wall shear stress via
   ! wall functions.
   rbard = half*(ww2d(i, j, irho)+ww1d(i, j, irho))
   rbar = half*(ww2(i, j, irho)+ww1(i, j, irho))
   ubard = half*(ww2d(i, j, ivx)+ww1d(i, j, ivx))
   ubar = half*(ww2(i, j, ivx)+ww1(i, j, ivx))
   vbard = half*(ww2d(i, j, ivy)+ww1d(i, j, ivy))
   vbar = half*(ww2(i, j, ivy)+ww1(i, j, ivy))
   wbard = half*(ww2d(i, j, ivz)+ww1d(i, j, ivz))
   wbar = half*(ww2(i, j, ivz)+ww1(i, j, ivz))
   ! Compute the velocity difference between the internal cell
   ! and the wall.
   vxd = ww2d(i, j, ivx) - ubard
   vx = ww2(i, j, ivx) - ubar
   vyd = ww2d(i, j, ivy) - vbard
   vy = ww2(i, j, ivy) - vbar
   vzd = ww2d(i, j, ivz) - wbard
   vz = ww2(i, j, ivz) - wbar
   ! Compute the normal velocity of the internal cell.
   velnd = norm(i, j, 1)*vxd + norm(i, j, 2)*vyd + norm(i, j, 3)*&
   &            vzd
   veln = vx*norm(i, j, 1) + vy*norm(i, j, 2) + vz*norm(i, j, 3)
   velnxd = norm(i, j, 1)*velnd
   velnx = veln*norm(i, j, 1)
   velnyd = norm(i, j, 2)*velnd
   velny = veln*norm(i, j, 2)
   velnzd = norm(i, j, 3)*velnd
   velnz = veln*norm(i, j, 3)
   ! Compute the tangential velocity, its magnitude and its
   ! unit vector of the internal cell.
   veltxd = vxd - velnxd
   veltx = vx - velnx
   veltyd = vyd - velnyd
   velty = vy - velny
   veltzd = vzd - velnzd
   veltz = vz - velnz
   arg1d = 2*veltx*veltxd + 2*velty*veltyd + 2*veltz*veltzd
   arg1 = veltx**2 + velty**2 + veltz**2
   IF (arg1 .EQ. 0.0) THEN
   y1d = 0.0
   ELSE
   y1d = arg1d/(2.0*SQRT(arg1))
   END IF
   y1 = SQRT(arg1)
   IF (eps .LT. y1) THEN
   veltmagd = y1d
   veltmag = y1
   ELSE
   veltmag = eps
   veltmagd = 0.0
   END IF
   txd = (veltxd*veltmag-veltx*veltmagd)/veltmag**2
   tx = veltx/veltmag
   tyd = (veltyd*veltmag-velty*veltmagd)/veltmag**2
   ty = velty/veltmag
   tzd = (veltzd*veltmag-veltz*veltmagd)/veltmag**2
   tz = veltz/veltmag
   ! Compute some coefficients needed for the transformation
   ! between the cartesian frame and the frame defined by the
   ! tangential direction (tx,ty,tz) and the normal direction.
   ! The minus sign is present, because for this transformation
   ! the normal direction should be inward pointing and norm
   ! is outward pointing.
   txnxd = -(norm(i, j, 1)*txd)
   txnx = -(tx*norm(i, j, 1))
   txnyd = -(norm(i, j, 2)*txd)
   txny = -(tx*norm(i, j, 2))
   txnzd = -(norm(i, j, 3)*txd)
   txnz = -(tx*norm(i, j, 3))
   tynxd = -(norm(i, j, 1)*tyd)
   tynx = -(ty*norm(i, j, 1))
   tynyd = -(norm(i, j, 2)*tyd)
   tyny = -(ty*norm(i, j, 2))
   tynzd = -(norm(i, j, 3)*tyd)
   tynz = -(ty*norm(i, j, 3))
   tznxd = -(norm(i, j, 1)*tzd)
   tznx = -(tz*norm(i, j, 1))
   tznyd = -(norm(i, j, 2)*tzd)
   tzny = -(tz*norm(i, j, 2))
   tznzd = -(norm(i, j, 3)*tzd)
   tznz = -(tz*norm(i, j, 3))
   ! Compute the tn component of the wall shear stress
   ! tensor. Normally this is the only nonzero shear
   ! stress component in the t-n frame.
   tautnd = tauxxd*txnx + tauxx*txnxd + tauyyd*tyny + tauyy*tynyd&
   &            + tauzzd*tznz + tauzz*tznzd + tauxyd*(txny+tynx) + tauxy*(&
   &            txnyd+tynxd) + tauxzd*(txnz+tznx) + tauxz*(txnzd+tznxd) + &
   &            tauyzd*(tynz+tzny) + tauyz*(tynzd+tznyd)
   tautn = tauxx*txnx + tauyy*tyny + tauzz*tznz + tauxy*(txny+&
   &            tynx) + tauxz*(txnz+tznx) + tauyz*(tynz+tzny)
   ! Compute the Reynolds number using the velocity, density,
   ! laminar viscosity and wall distance. Note that an offset
   ! of -1 must be used in dd2Wall2, because the original array
   ! d2Wall starts at 2.
   red = dd2wall2(i-1, j-1)*(ww2d(i, j, irho)*veltmag+ww2(i, j, &
   &            irho)*veltmagd)/rrlv2(i, j)
   re = ww2(i, j, irho)*veltmag*dd2wall2(i-1, j-1)/rrlv2(i, j)
   x1d = CURVEUPRE_SPATIAL_D(re, red, x1)
   IF (x1 .LT. eps) THEN
   max1 = eps
   max1d = 0.0
   ELSE
   max1d = x1d
   max1 = x1
   END IF
   ! Determine the friction velocity from the table and
   ! compute the wall shear stress from it.
   utaud = (veltmagd*max1-veltmag*max1d)/max1**2
   utau = veltmag/max1
   tauwalld = (rbard*utau+rbar*utaud)*utau + rbar*utau*utaud
   tauwall = rbar*utau*utau
   ! Compute the correction to the wall shear stress tautn and
   ! transform this correction back to the cartesian frame.
   ! Take rFilv into account, such that the correction to the
   ! stress tensor is computed correctly.
   tautnd = rfilv*tauwalld - tautnd
   tautn = rfilv*tauwall - tautn
   tauxxd = two*(tautnd*txnx+tautn*txnxd)
   tauxx = two*tautn*txnx
   tauyyd = two*(tautnd*tyny+tautn*tynyd)
   tauyy = two*tautn*tyny
   tauzzd = two*(tautnd*tznz+tautn*tznzd)
   tauzz = two*tautn*tznz
   tauxyd = tautnd*(txny+tynx) + tautn*(txnyd+tynxd)
   tauxy = tautn*(txny+tynx)
   tauxzd = tautnd*(txnz+tznx) + tautn*(txnzd+tznxd)
   tauxz = tautn*(txnz+tznx)
   tauyzd = tautnd*(tynz+tzny) + tautn*(tynzd+tznyd)
   tauyz = tautn*(tynz+tzny)
   ! Compute the correction to the viscous flux at the wall.
   fmxd = ss(i, j, 1)*tauxxd + ss(i, j, 2)*tauxyd + ss(i, j, 3)*&
   &            tauxzd
   fmx = tauxx*ss(i, j, 1) + tauxy*ss(i, j, 2) + tauxz*ss(i, j, 3&
   &            )
   fmyd = ss(i, j, 1)*tauxyd + ss(i, j, 2)*tauyyd + ss(i, j, 3)*&
   &            tauyzd
   fmy = tauxy*ss(i, j, 1) + tauyy*ss(i, j, 2) + tauyz*ss(i, j, 3&
   &            )
   fmzd = ss(i, j, 1)*tauxzd + ss(i, j, 2)*tauyzd + ss(i, j, 3)*&
   &            tauzzd
   fmz = tauxz*ss(i, j, 1) + tauyz*ss(i, j, 2) + tauzz*ss(i, j, 3&
   &            )
   frhoed = ss(i, j, 1)*(ubard*tauxx+ubar*tauxxd+vbard*tauxy+vbar&
   &            *tauxyd+wbard*tauxz+wbar*tauxzd) + ss(i, j, 2)*(ubard*tauxy+&
   &            ubar*tauxyd+vbard*tauyy+vbar*tauyyd+wbard*tauyz+wbar*tauyzd)&
   &            + ss(i, j, 3)*(ubard*tauxz+ubar*tauxzd+vbard*tauyz+vbar*&
   &            tauyzd+wbard*tauzz+wbar*tauzzd)
   frhoe = (ubar*tauxx+vbar*tauxy+wbar*tauxz)*ss(i, j, 1) + (ubar&
   &            *tauxy+vbar*tauyy+wbar*tauyz)*ss(i, j, 2) + (ubar*tauxz+vbar&
   &            *tauyz+wbar*tauzz)*ss(i, j, 3)
   ! Add them to the residual. Note that now the factor rFilv
   ! is already taken into account via tau. Fact is present to
   ! take inward/outward pointing normals into account
   rresd(i, j, imx) = rresd(i, j, imx) - fact*fmxd
   rres(i, j, imx) = rres(i, j, imx) - fact*fmx
   rresd(i, j, imy) = rresd(i, j, imy) - fact*fmyd
   rres(i, j, imy) = rres(i, j, imy) - fact*fmy
   rresd(i, j, imz) = rresd(i, j, imz) - fact*fmzd
   rres(i, j, imz) = rres(i, j, imz) - fact*fmz
   rresd(i, j, irhoe) = rresd(i, j, irhoe) - fact*frhoed
   rres(i, j, irhoe) = rres(i, j, irhoe) - fact*frhoe
   ! Store the friction velocity for later use.
   viscsubface(nn)%utau(i, j) = utau
   ! Also add the correction to the wall stress tensor.
   viscsubfaced(nn)%tau(i, j, 1) = viscsubfaced(nn)%tau(i, j, 1) &
   &            + tauxxd
   viscsubface(nn)%tau(i, j, 1) = viscsubface(nn)%tau(i, j, 1) + &
   &            tauxx
   viscsubfaced(nn)%tau(i, j, 2) = viscsubfaced(nn)%tau(i, j, 2) &
   &            + tauyyd
   viscsubface(nn)%tau(i, j, 2) = viscsubface(nn)%tau(i, j, 2) + &
   &            tauyy
   viscsubfaced(nn)%tau(i, j, 3) = viscsubfaced(nn)%tau(i, j, 3) &
   &            + tauzzd
   viscsubface(nn)%tau(i, j, 3) = viscsubface(nn)%tau(i, j, 3) + &
   &            tauzz
   viscsubfaced(nn)%tau(i, j, 4) = viscsubfaced(nn)%tau(i, j, 4) &
   &            + tauxyd
   viscsubface(nn)%tau(i, j, 4) = viscsubface(nn)%tau(i, j, 4) + &
   &            tauxy
   viscsubfaced(nn)%tau(i, j, 5) = viscsubfaced(nn)%tau(i, j, 5) &
   &            + tauxzd
   viscsubface(nn)%tau(i, j, 5) = viscsubface(nn)%tau(i, j, 5) + &
   &            tauxz
   viscsubfaced(nn)%tau(i, j, 6) = viscsubfaced(nn)%tau(i, j, 6) &
   &            + tauyzd
   viscsubface(nn)%tau(i, j, 6) = viscsubface(nn)%tau(i, j, 6) + &
   &            tauyz
   END DO
   END DO
   END DO viscsubfaces
   END IF
   END SUBROUTINE UTAUWF_SPATIAL_D
