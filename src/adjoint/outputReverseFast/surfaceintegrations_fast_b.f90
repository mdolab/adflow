!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
module surfaceintegrations_fast_b
  implicit none
! ----------------------------------------------------------------------
!                                                                      |
!                    no tapenade routine below this line               |
!                                                                      |
! ----------------------------------------------------------------------

contains
  subroutine getcostfunctions(globalvals, funcvalues)
    use constants
    use inputtimespectral, only : ntimeintervalsspectral
    use flowvarrefstate, only : pref, rhoref, tref, lref, gammainf
    use inputphysics, only : liftdirection, dragdirection, surfaceref,&
&   machcoef, lengthref
    use inputtsstabderiv, only : tsstability
    use utils_fast_b, only : computetsderivatives
    implicit none
! input/output
    real(kind=realtype), dimension(:, :), intent(in) :: globalvals
    real(kind=realtype), dimension(:), intent(out) :: funcvalues
! working
    real(kind=realtype) :: fact, factmoment, ovrnts
    real(kind=realtype), dimension(3, ntimeintervalsspectral) :: force, &
&   moment, cforce, cmoment
    real(kind=realtype) :: mavgptot, mavgttot, mavgps, mflow, mflow2, &
&   mavgmn, sigmamn, sigmaptot
    integer(kind=inttype) :: sps
    real(kind=realtype), dimension(8) :: dcdq, dcdqdot
    real(kind=realtype), dimension(8) :: dcdalpha, dcdalphadot
    real(kind=realtype), dimension(8) :: coef0
    intrinsic sqrt
! factor used for time-averaged quantities.
    ovrnts = one/ntimeintervalsspectral
! sum pressure and viscous contributions
    force = globalvals(ifp:ifp+2, :) + globalvals(ifv:ifv+2, :) + &
&     globalvals(iflowfm:iflowfm+2, :)
    moment = globalvals(imp:imp+2, :) + globalvals(imv:imv+2, :) + &
&     globalvals(iflowmm:iflowmm+2, :)
    fact = two/(gammainf*machcoef*machcoef*surfaceref*lref*lref*pref)
    cforce = fact*force
! moment factor has an extra lengthref
    fact = fact/(lengthref*lref)
    cmoment = fact*moment
! zero values since we are summing.
    funcvalues = zero
! here we finally assign the final function values
    do sps=1,ntimeintervalsspectral
      funcvalues(costfuncforcex) = funcvalues(costfuncforcex) + ovrnts*&
&       force(1, sps)
      funcvalues(costfuncforcey) = funcvalues(costfuncforcey) + ovrnts*&
&       force(2, sps)
      funcvalues(costfuncforcez) = funcvalues(costfuncforcez) + ovrnts*&
&       force(3, sps)
      funcvalues(costfuncforcexcoef) = funcvalues(costfuncforcexcoef) + &
&       ovrnts*cforce(1, sps)
      funcvalues(costfuncforceycoef) = funcvalues(costfuncforceycoef) + &
&       ovrnts*cforce(2, sps)
      funcvalues(costfuncforcezcoef) = funcvalues(costfuncforcezcoef) + &
&       ovrnts*cforce(3, sps)
      funcvalues(costfuncmomx) = funcvalues(costfuncmomx) + ovrnts*&
&       moment(1, sps)
      funcvalues(costfuncmomy) = funcvalues(costfuncmomy) + ovrnts*&
&       moment(2, sps)
      funcvalues(costfuncmomz) = funcvalues(costfuncmomz) + ovrnts*&
&       moment(3, sps)
      funcvalues(costfuncmomxcoef) = funcvalues(costfuncmomxcoef) + &
&       ovrnts*cmoment(1, sps)
      funcvalues(costfuncmomycoef) = funcvalues(costfuncmomycoef) + &
&       ovrnts*cmoment(2, sps)
      funcvalues(costfuncmomzcoef) = funcvalues(costfuncmomzcoef) + &
&       ovrnts*cmoment(3, sps)
      funcvalues(costfuncsepsensor) = funcvalues(costfuncsepsensor) + &
&       ovrnts*globalvals(isepsensor, sps)
      funcvalues(costfunccavitation) = funcvalues(costfunccavitation) + &
&       ovrnts*globalvals(icavitation, sps)
      funcvalues(costfuncsepsensoravgx) = funcvalues(&
&       costfuncsepsensoravgx) + ovrnts*globalvals(isepavg, sps)
      funcvalues(costfuncsepsensoravgy) = funcvalues(&
&       costfuncsepsensoravgy) + ovrnts*globalvals(isepavg+1, sps)
      funcvalues(costfuncsepsensoravgz) = funcvalues(&
&       costfuncsepsensoravgz) + ovrnts*globalvals(isepavg+2, sps)
      funcvalues(costfuncpk) = funcvalues(costfuncpk) + ovrnts*&
&       globalvals(ipk, sps)
! mass flow like objective
      mflow = globalvals(imassflow, sps)
      if (mflow .ne. zero) then
        mavgptot = globalvals(imassptot, sps)/mflow
        mavgttot = globalvals(imassttot, sps)/mflow
        mavgps = globalvals(imassps, sps)/mflow
        mavgmn = globalvals(imassmn, sps)/mflow
        mflow2 = globalvals(imassflow, sps)*sqrt(pref/rhoref)
        sigmamn = sqrt(globalvals(isigmamn, sps)/mflow)
        sigmaptot = sqrt(globalvals(isigmaptot, sps)/mflow)
      else
        mavgptot = zero
        mavgttot = zero
        mavgps = zero
        mavgmn = zero
        sigmamn = zero
        sigmaptot = zero
        mflow2 = zero
      end if
      funcvalues(costfuncmdot) = funcvalues(costfuncmdot) + ovrnts*mflow
      funcvalues(costfuncmavgptot) = funcvalues(costfuncmavgptot) + &
&       ovrnts*mavgptot
      funcvalues(costfuncmavgttot) = funcvalues(costfuncmavgttot) + &
&       ovrnts*mavgttot
      funcvalues(costfuncmavgps) = funcvalues(costfuncmavgps) + ovrnts*&
&       mavgps
      funcvalues(costfuncmavgmn) = funcvalues(costfuncmavgmn) + ovrnts*&
&       mavgmn
      funcvalues(costfuncsigmamn) = funcvalues(costfuncsigmamn) + ovrnts&
&       *sigmamn
      funcvalues(costfuncsigmaptot) = funcvalues(costfuncsigmaptot) + &
&       ovrnts*sigmaptot
      funcvalues(costfuncpk) = funcvalues(costfuncpk) + ovrnts*&
&       globalvals(ipk, sps)
    end do
! bending moment calc - also broken. 
! call computerootbendingmoment(cforce, cmoment, liftindex, bendingmoment)
! funcvalues(costfuncbendingcoef) = funcvalues(costfuncbendingcoef) + ovrnts*bendingmoment
! lift and drag (coefficients): dot product with the lift/drag direction.
    funcvalues(costfunclift) = funcvalues(costfuncforcex)*liftdirection(&
&     1) + funcvalues(costfuncforcey)*liftdirection(2) + funcvalues(&
&     costfuncforcez)*liftdirection(3)
    funcvalues(costfuncdrag) = funcvalues(costfuncforcex)*dragdirection(&
&     1) + funcvalues(costfuncforcey)*dragdirection(2) + funcvalues(&
&     costfuncforcez)*dragdirection(3)
    funcvalues(costfuncliftcoef) = funcvalues(costfuncforcexcoef)*&
&     liftdirection(1) + funcvalues(costfuncforceycoef)*liftdirection(2)&
&     + funcvalues(costfuncforcezcoef)*liftdirection(3)
    funcvalues(costfuncdragcoef) = funcvalues(costfuncforcexcoef)*&
&     dragdirection(1) + funcvalues(costfuncforceycoef)*dragdirection(2)&
&     + funcvalues(costfuncforcezcoef)*dragdirection(3)
! -------------------- time spectral objectives ------------------
    if (tsstability) then
      print*, &
&     'error: tsstabilityderivatives are *broken*. they need to be ', &
&     'completely verifed from scratch'
      stop
    end if
  end subroutine getcostfunctions
  subroutine wallintegrationface(localvalues, mm)
!
!       wallintegrations computes the contribution of the block
!       given by the pointers in blockpointers to the force and
!       moment of the geometry. a distinction is made
!       between the inviscid and viscous parts. in case the maximum
!       yplus value must be monitored (only possible for rans), this
!       value is also computed. the separation sensor and the cavita-
!       tion sensor is also computed
!       here.
!
    use constants
    use communication
    use blockpointers
    use flowvarrefstate
    use inputcostfunctions
    use inputphysics, only : machcoef, pointref, veldirfreestream, &
&   equations
    use bcpointers_fast_b
    implicit none
! input/output variables
    real(kind=realtype), dimension(nlocalvalues), intent(inout) :: &
&   localvalues
    integer(kind=inttype) :: mm
! local variables.
    real(kind=realtype), dimension(3) :: fp, fv, mp, mv
    real(kind=realtype) :: yplusmax, sepsensor, sepsensoravg(3), &
&   cavitation
    integer(kind=inttype) :: i, j, ii, blk
    real(kind=realtype) :: pm1, fx, fy, fz, fn, sigma
    real(kind=realtype) :: xc, yc, zc, qf(3)
    real(kind=realtype) :: fact, rho, mul, yplus, dwall
    real(kind=realtype) :: v(3), sensor, sensor1, cp, tmp, plocal
    real(kind=realtype) :: tauxx, tauyy, tauzz
    real(kind=realtype) :: tauxy, tauxz, tauyz
    real(kind=realtype), dimension(3) :: refpoint
    real(kind=realtype) :: mx, my, mz, cellarea
    intrinsic mod
    intrinsic max
    intrinsic sqrt
    intrinsic exp
    select case  (bcfaceid(mm)) 
    case (imin, jmin, kmin) 
      fact = -one
    case (imax, jmax, kmax) 
      fact = one
    end select
! determine the reference point for the moment computation in
! meters.
    refpoint(1) = lref*pointref(1)
    refpoint(2) = lref*pointref(2)
    refpoint(3) = lref*pointref(3)
! initialize the force and moment coefficients to 0 as well as
! yplusmax.
    fp = zero
    fv = zero
    mp = zero
    mv = zero
    yplusmax = zero
    sepsensor = zero
    cavitation = zero
    sepsensoravg = zero
!
!         integrate the inviscid contribution over the solid walls,
!         either inviscid or viscous. the integration is done with
!         cp. for closed contours this is equal to the integration
!         of p; for open contours this is not the case anymore.
!         question is whether a force for an open contour is
!         meaningful anyway.
!
! loop over the quadrilateral faces of the subface. note that
! the nodal range of bcdata must be used and not the cell
! range, because the latter may include the halo's in i and
! j-direction. the offset +1 is there, because inbeg and jnbeg
! refer to nodal ranges and not to cell ranges. the loop
! (without the ad stuff) would look like:
!
! do j=(bcdata(mm)%jnbeg+1),bcdata(mm)%jnend
!    do i=(bcdata(mm)%inbeg+1),bcdata(mm)%inend
    do ii=0,(bcdata(mm)%jnend-bcdata(mm)%jnbeg)*(bcdata(mm)%inend-bcdata&
&       (mm)%inbeg)-1
      i = mod(ii, bcdata(mm)%inend - bcdata(mm)%inbeg) + bcdata(mm)%&
&       inbeg + 1
      j = ii/(bcdata(mm)%inend-bcdata(mm)%inbeg) + bcdata(mm)%jnbeg + 1
! compute the average pressure minus 1 and the coordinates
! of the centroid of the face relative from from the
! moment reference point. due to the usage of pointers for
! the coordinates, whose original array starts at 0, an
! offset of 1 must be used. the pressure is multipled by
! fact to account for the possibility of an inward or
! outward pointing normal.
      pm1 = fact*(half*(pp2(i, j)+pp1(i, j))-pinf)*pref
      xc = fourth*(xx(i, j, 1)+xx(i+1, j, 1)+xx(i, j+1, 1)+xx(i+1, j+1, &
&       1)) - refpoint(1)
      yc = fourth*(xx(i, j, 2)+xx(i+1, j, 2)+xx(i, j+1, 2)+xx(i+1, j+1, &
&       2)) - refpoint(2)
      zc = fourth*(xx(i, j, 3)+xx(i+1, j, 3)+xx(i, j+1, 3)+xx(i+1, j+1, &
&       3)) - refpoint(3)
      if (bcdata(mm)%iblank(i, j) .lt. 0) then
        blk = 0
      else
        blk = bcdata(mm)%iblank(i, j)
      end if
      fx = pm1*ssi(i, j, 1)
      fy = pm1*ssi(i, j, 2)
      fz = pm1*ssi(i, j, 3)
! iblank forces
      fx = fx*blk
      fy = fy*blk
      fz = fz*blk
! update the inviscid force and moment coefficients.
      fp(1) = fp(1) + fx
      fp(2) = fp(2) + fy
      fp(3) = fp(3) + fz
      mx = yc*fz - zc*fy
      my = zc*fx - xc*fz
      mz = xc*fy - yc*fx
      mp(1) = mp(1) + mx
      mp(2) = mp(2) + my
      mp(3) = mp(3) + mz
! save the face-based forces and area
      bcdata(mm)%fp(i, j, 1) = fx
      bcdata(mm)%fp(i, j, 2) = fy
      bcdata(mm)%fp(i, j, 3) = fz
      cellarea = sqrt(ssi(i, j, 1)**2 + ssi(i, j, 2)**2 + ssi(i, j, 3)**&
&       2)
      bcdata(mm)%area(i, j) = cellarea
! get normalized surface velocity:
      v(1) = ww2(i, j, ivx)
      v(2) = ww2(i, j, ivy)
      v(3) = ww2(i, j, ivz)
      v = v/(sqrt(v(1)**2+v(2)**2+v(3)**2)+1e-16)
! dot product with free stream
      sensor = -(v(1)*veldirfreestream(1)+v(2)*veldirfreestream(2)+v(3)*&
&       veldirfreestream(3))
!now run through a smooth heaviside function:
      sensor = one/(one+exp(-(2*sepsensorsharpness*(sensor-&
&       sepsensoroffset))))
! and integrate over the area of this cell and save:
      sensor = sensor*cellarea
      sepsensor = sepsensor + sensor
! also accumulate into the sepsensoravg
      xc = fourth*(xx(i, j, 1)+xx(i+1, j, 1)+xx(i, j+1, 1)+xx(i+1, j+1, &
&       1))
      yc = fourth*(xx(i, j, 2)+xx(i+1, j, 2)+xx(i, j+1, 2)+xx(i+1, j+1, &
&       2))
      zc = fourth*(xx(i, j, 3)+xx(i+1, j, 3)+xx(i, j+1, 3)+xx(i+1, j+1, &
&       3))
      sepsensoravg(1) = sepsensoravg(1) + sensor*xc
      sepsensoravg(2) = sepsensoravg(2) + sensor*yc
      sepsensoravg(3) = sepsensoravg(3) + sensor*zc
      plocal = pp2(i, j)
      tmp = two/(gammainf*machcoef*machcoef)
      cp = tmp*(plocal-pinf)
      sigma = 1.4
      sensor1 = -cp - sigma
      sensor1 = one/(one+exp(-(2*10*sensor1)))
      sensor1 = sensor1*cellarea
      cavitation = cavitation + sensor1
    end do
!
! integration of the viscous forces.
! only for viscous boundaries.
!
    if (bctype(mm) .eq. nswalladiabatic .or. bctype(mm) .eq. &
&       nswallisothermal) then
! initialize dwall for the laminar case and set the pointer
! for the unit normals.
      dwall = zero
! loop over the quadrilateral faces of the subface and
! compute the viscous contribution to the force and
! moment and update the maximum value of y+.
      do ii=0,(bcdata(mm)%jnend-bcdata(mm)%jnbeg)*(bcdata(mm)%inend-&
&         bcdata(mm)%inbeg)-1
        i = mod(ii, bcdata(mm)%inend - bcdata(mm)%inbeg) + bcdata(mm)%&
&         inbeg + 1
        j = ii/(bcdata(mm)%inend-bcdata(mm)%inbeg) + bcdata(mm)%jnbeg + &
&         1
        if (bcdata(mm)%iblank(i, j) .lt. 0) then
          blk = 0
        else
          blk = bcdata(mm)%iblank(i, j)
        end if
        tauxx = viscsubface(mm)%tau(i, j, 1)
        tauyy = viscsubface(mm)%tau(i, j, 2)
        tauzz = viscsubface(mm)%tau(i, j, 3)
        tauxy = viscsubface(mm)%tau(i, j, 4)
        tauxz = viscsubface(mm)%tau(i, j, 5)
        tauyz = viscsubface(mm)%tau(i, j, 6)
! compute the viscous force on the face. a minus sign
! is now present, due to the definition of this force.
        fx = -(fact*(tauxx*ssi(i, j, 1)+tauxy*ssi(i, j, 2)+tauxz*ssi(i, &
&         j, 3))*pref)
        fy = -(fact*(tauxy*ssi(i, j, 1)+tauyy*ssi(i, j, 2)+tauyz*ssi(i, &
&         j, 3))*pref)
        fz = -(fact*(tauxz*ssi(i, j, 1)+tauyz*ssi(i, j, 2)+tauzz*ssi(i, &
&         j, 3))*pref)
! iblank forces after saving for zipper mesh
        tauxx = tauxx*blk
        tauyy = tauyy*blk
        tauzz = tauzz*blk
        tauxy = tauxy*blk
        tauxz = tauxz*blk
        tauyz = tauyz*blk
        fx = fx*blk
        fy = fy*blk
        fz = fz*blk
! compute the coordinates of the centroid of the face
! relative from the moment reference point. due to the
! usage of pointers for xx and offset of 1 is present,
! because x originally starts at 0.
        xc = fourth*(xx(i, j, 1)+xx(i+1, j, 1)+xx(i, j+1, 1)+xx(i+1, j+1&
&         , 1)) - refpoint(1)
        yc = fourth*(xx(i, j, 2)+xx(i+1, j, 2)+xx(i, j+1, 2)+xx(i+1, j+1&
&         , 2)) - refpoint(2)
        zc = fourth*(xx(i, j, 3)+xx(i+1, j, 3)+xx(i, j+1, 3)+xx(i+1, j+1&
&         , 3)) - refpoint(3)
! update the viscous force and moment coefficients.
        fv(1) = fv(1) + fx
        fv(2) = fv(2) + fy
        fv(3) = fv(3) + fz
        mx = yc*fz - zc*fy
        my = zc*fx - xc*fz
        mz = xc*fy - yc*fx
        mv(1) = mv(1) + mx
        mv(2) = mv(2) + my
        mv(3) = mv(3) + mz
! save the face based forces for the slice operations
        bcdata(mm)%fv(i, j, 1) = fx
        bcdata(mm)%fv(i, j, 2) = fy
        bcdata(mm)%fv(i, j, 3) = fz
! compute the tangential component of the stress tensor,
! which is needed to monitor y+. the result is stored
! in fx, fy, fz, although it is not really a force.
! as later on only the magnitude of the tangential
! component is important, there is no need to take the
! sign into account (it should be a minus sign).
        fx = tauxx*bcdata(mm)%norm(i, j, 1) + tauxy*bcdata(mm)%norm(i, j&
&         , 2) + tauxz*bcdata(mm)%norm(i, j, 3)
        fy = tauxy*bcdata(mm)%norm(i, j, 1) + tauyy*bcdata(mm)%norm(i, j&
&         , 2) + tauyz*bcdata(mm)%norm(i, j, 3)
        fz = tauxz*bcdata(mm)%norm(i, j, 1) + tauyz*bcdata(mm)%norm(i, j&
&         , 2) + tauzz*bcdata(mm)%norm(i, j, 3)
        fn = fx*bcdata(mm)%norm(i, j, 1) + fy*bcdata(mm)%norm(i, j, 2) +&
&         fz*bcdata(mm)%norm(i, j, 3)
        fx = fx - fn*bcdata(mm)%norm(i, j, 1)
        fy = fy - fn*bcdata(mm)%norm(i, j, 2)
        fz = fz - fn*bcdata(mm)%norm(i, j, 3)
      end do
    else
! compute the local value of y+. due to the usage
! of pointers there is on offset of -1 in dd2wall..
! if we had no viscous force, set the viscous component to zero
      bcdata(mm)%fv = zero
    end if
! increment the local values array with the values we computed here.
    localvalues(ifp:ifp+2) = localvalues(ifp:ifp+2) + fp
    localvalues(ifv:ifv+2) = localvalues(ifv:ifv+2) + fv
    localvalues(imp:imp+2) = localvalues(imp:imp+2) + mp
    localvalues(imv:imv+2) = localvalues(imv:imv+2) + mv
    localvalues(isepsensor) = localvalues(isepsensor) + sepsensor
    localvalues(icavitation) = localvalues(icavitation) + cavitation
    localvalues(isepavg:isepavg+2) = localvalues(isepavg:isepavg+2) + &
&     sepsensoravg
  end subroutine wallintegrationface
  subroutine flowintegrationface(isinflow, localvalues, mm, withgathered&
&   , funcvalues)
    use constants
    use blockpointers, only : bctype, bcfaceid, bcdata, &
&   addgridvelocities
    use flowvarrefstate, only : pref, pinf, rhoref, timeref, lref, &
&   tref, rgas, uref, uinf
    use inputphysics, only : pointref, flowtype, veldirfreestream, &
&   alpha, beta, liftindex
    use flowutils_fast_b, only : computeptot, computettot, getdirvector
    use bcpointers_fast_b, only : ssi, sface, ww1, ww2, pp1, pp2, xx, gamma1,&
&   gamma2
    use utils_fast_b, only : mynorm2
    implicit none
! input/output variables
    logical, intent(in) :: isinflow
    logical, intent(in) :: withgathered
    real(kind=realtype), dimension(nlocalvalues), intent(inout) :: &
&   localvalues
    integer(kind=inttype), intent(in) :: mm
    real(kind=realtype), dimension(:), optional, intent(in) :: &
&   funcvalues
! local variables
    real(kind=realtype) :: massflowrate, mass_ptot, mass_ttot, mass_ps, &
&   mass_mn
    real(kind=realtype) :: mredim, pk, sigma_mn, sigma_ptot
    integer(kind=inttype) :: i, j, ii, blk
    real(kind=realtype) :: internalflowfact, inflowfact, fact, xc, yc, &
&   zc, cellarea, mx, my, mz
    real(kind=realtype) :: sf, vmag, vnm, vxm, vym, vzm, fx, fy, fz, u, &
&   v, w
    real(kind=realtype) :: pm, ptot, ttot, rhom, gammam, a2
    real(kind=realtype), dimension(3) :: fp, mp, fmom, mmom, refpoint, &
&   vcoordref, vfreestreamref, sfacefreestreamref
    real(kind=realtype) :: mnm, massflowratelocal
    intrinsic sqrt
    intrinsic mod
    intrinsic max
    refpoint(1) = lref*pointref(1)
    refpoint(2) = lref*pointref(2)
    refpoint(3) = lref*pointref(3)
! note that these are *opposite* of force integrations. the reason
! is that we want positive mass flow into the domain and negative
! mass flow out of the domain. since the low faces have ssi
! vectors pointining into the domain, this is correct. the high
! end faces need to flip this. 
    select case  (bcfaceid(mm)) 
    case (imin, jmin, kmin) 
      fact = one
    case (imax, jmax, kmax) 
      fact = -one
    end select
! the sign of momentum forces are flipped for internal flows
    internalflowfact = one
    if (flowtype .eq. internalflow) internalflowfact = -one
    inflowfact = one
    if (isinflow) inflowfact = -one
! loop over the quadrilateral faces of the subface. note that
! the nodal range of bcdata must be used and not the cell
! range, because the latter may include the halo's in i and
! j-direction. the offset +1 is there, because inbeg and jnbeg
! refer to nodal ranges and not to cell ranges. the loop
! (without the ad stuff) would look like:
!
! do j=(bcdata(mm)%jnbeg+1),bcdata(mm)%jnend
!    do i=(bcdata(mm)%inbeg+1),bcdata(mm)%inend
    fp = zero
    mp = zero
    fmom = zero
    mmom = zero
    pk = zero
    mredim = sqrt(pref*rhoref)
    massflowrate = zero
    mass_ptot = zero
    mass_ttot = zero
    mass_ps = zero
    mass_mn = zero
    sigma_mn = zero
    sigma_ptot = zero
    do ii=0,(bcdata(mm)%jnend-bcdata(mm)%jnbeg)*(bcdata(mm)%inend-bcdata&
&       (mm)%inbeg)-1
      i = mod(ii, bcdata(mm)%inend - bcdata(mm)%inbeg) + bcdata(mm)%&
&       inbeg + 1
      j = ii/(bcdata(mm)%inend-bcdata(mm)%inbeg) + bcdata(mm)%jnbeg + 1
      if (addgridvelocities) then
        sf = sface(i, j)
      else
        sf = zero
      end if
      if (bcdata(mm)%iblank(i, j) .lt. 0) then
        blk = 0
      else
        blk = bcdata(mm)%iblank(i, j)
      end if
      vxm = half*(ww1(i, j, ivx)+ww2(i, j, ivx))
      vym = half*(ww1(i, j, ivy)+ww2(i, j, ivy))
      vzm = half*(ww1(i, j, ivz)+ww2(i, j, ivz))
      rhom = half*(ww1(i, j, irho)+ww2(i, j, irho))
      pm = half*(pp1(i, j)+pp2(i, j))
      gammam = half*(gamma1(i, j)+gamma2(i, j))
      vnm = vxm*ssi(i, j, 1) + vym*ssi(i, j, 2) + vzm*ssi(i, j, 3) - sf
      vmag = sqrt(vxm**2 + vym**2 + vzm**2) - sf
! a = sqrt(gamma*p/rho); sqrt(v**2/a**2)
      mnm = vmag/sqrt(gammam*pm/rhom)
      call computeptot(rhom, vxm, vym, vzm, pm, ptot)
      call computettot(rhom, vxm, vym, vzm, pm, ttot)
      massflowratelocal = rhom*vnm*blk*fact*mredim
      if (withgathered) then
        sigma_mn = sigma_mn + massflowratelocal*(mnm-funcvalues(&
&         costfuncmavgmn))**2
        sigma_ptot = sigma_ptot + massflowratelocal*(ptot-funcvalues(&
&         costfuncmavgptot))**2
      else
        massflowrate = massflowrate + massflowratelocal
        pk = pk + (pm-pinf+half*rhom*(vmag**2-uinf**2))*vnm*pref*uref*&
&         fact
! re-dimentionalize quantities
        pm = pm*pref
        mass_ptot = mass_ptot + ptot*massflowratelocal*pref
        mass_ttot = mass_ttot + ttot*massflowratelocal*tref
        mass_ps = mass_ps + pm*massflowratelocal
        mass_mn = mass_mn + mnm*massflowratelocal
        xc = fourth*(xx(i, j, 1)+xx(i+1, j, 1)+xx(i, j+1, 1)+xx(i+1, j+1&
&         , 1)) - refpoint(1)
        yc = fourth*(xx(i, j, 2)+xx(i+1, j, 2)+xx(i, j+1, 2)+xx(i+1, j+1&
&         , 2)) - refpoint(2)
        zc = fourth*(xx(i, j, 3)+xx(i+1, j, 3)+xx(i, j+1, 3)+xx(i+1, j+1&
&         , 3)) - refpoint(3)
! pressure forces. note that these need a *negative* and to subtract 
! the reference pressure sign to be consistent with the force 
! computation on the walls. 
        pm = -((pm-pinf*pref)*fact*blk)
        fx = pm*ssi(i, j, 1)
        fy = pm*ssi(i, j, 2)
        fz = pm*ssi(i, j, 3)
! update the pressure force and moment coefficients.
        fp(1) = fp(1) + fx
        fp(2) = fp(2) + fy
        fp(3) = fp(3) + fz
        mx = yc*fz - zc*fy
        my = zc*fx - xc*fz
        mz = xc*fy - yc*fx
        mp(1) = mp(1) + mx
        mp(2) = mp(2) + my
        mp(3) = mp(3) + mz
! momentum forces are a little tricky.  we negate because 
! have to re-apply fact to massflowratelocal to undoo it, because 
! we need the signed behavior of ssi to get the momentum forces correct. 
! also, the sign is flipped between inflow and outflow types 
        cellarea = mynorm2(ssi(i, j, :))
        massflowratelocal = massflowratelocal*fact/timeref*blk/cellarea*&
&         internalflowfact*inflowfact
        fx = massflowratelocal*ssi(i, j, 1)*vxm
        fy = massflowratelocal*ssi(i, j, 2)*vym
        fz = massflowratelocal*ssi(i, j, 3)*vzm
        fmom(1) = fmom(1) + fx
        fmom(2) = fmom(2) + fy
        fmom(3) = fmom(3) + fz
        mx = yc*fz - zc*fy
        my = zc*fx - xc*fz
        mz = xc*fy - yc*fx
        mmom(1) = mmom(1) + mx
        mmom(2) = mmom(2) + my
        mmom(3) = mmom(3) + mz
! ! computes the normalized vector maped into the freestream direction, so we multiply by the magnitude after
! vcoordref(1) = vxm
! vcoordref(2) = vym
! vcoordref(3) = vzm
! call getdirvector(vcoordref, -alpha, -beta, vfreestreamref, liftindex)
! vfreestreamref = vfreestreamref * vmag
! !project the face normal into the freestream velocity and scale by the face
! call getdirvector(ssi(i,j,:), -alpha, -beta, sfacefreestreamref, liftindex)
! sfacefreestreamref = sfacefreestreamref * sf
! ! compute the pertubations of the flow from the free-stream velocity
! u = vfreestreamref(1) - sfacefreestreamref(1) - uinf
! v = vfreestreamref(2) - sfacefreestreamref(2)
! w = vfreestreamref(3) - sfacefreestreamref(3)
! !edota = edota + half*(rhom)
      end if
    end do
    if (withgathered) then
      localvalues(isigmamn) = localvalues(isigmamn) + sigma_mn
      localvalues(isigmaptot) = localvalues(isigmaptot) + sigma_ptot
    else
! increment the local values array with what we computed here
      localvalues(imassflow) = localvalues(imassflow) + massflowrate
      localvalues(imassptot) = localvalues(imassptot) + mass_ptot
      localvalues(imassttot) = localvalues(imassttot) + mass_ttot
      localvalues(imassps) = localvalues(imassps) + mass_ps
      localvalues(imassmn) = localvalues(imassmn) + mass_mn
      localvalues(ipk) = localvalues(ipk) + pk
      localvalues(ifp:ifp+2) = localvalues(ifp:ifp+2) + fp
      localvalues(iflowfm:iflowfm+2) = localvalues(iflowfm:iflowfm+2) + &
&       fmom
      localvalues(iflowmp:iflowmp+2) = localvalues(iflowmp:iflowmp+2) + &
&       mp
      localvalues(iflowmm:iflowmm+2) = localvalues(iflowmm:iflowmm+2) + &
&       mmom
    end if
  end subroutine flowintegrationface
end module surfaceintegrations_fast_b
