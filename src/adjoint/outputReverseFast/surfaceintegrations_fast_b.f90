!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
module surfaceintegrations_fast_b
  use constants
  use communication, only : commtype, internalcommtype
  implicit none
! ----------------------------------------------------------------------
!                                                                      |
!                    no tapenade routine below this line               |
!                                                                      |
! ----------------------------------------------------------------------
! data required on each proc:
! ndonor: the number of donor points the proc will provide
! frac (3, ndonor) : the uvw coordinates of the interpolation point
! donorinfo(4, ndonor) : donor information. 1 is the local block id and 2-4 is the 
!    starting i,j,k indices for the interpolation. 
! procsizes(0:nproc-1) : the number of donors on each proc
! procdisps(0:nproc) : cumulative form of procsizes
! inv(nconn) : array allocated only on root processor used to
! reorder the nodes or elements back to the original order. 
  type usersurfcommtype
      integer(kind=inttype) :: ndonor
      real(kind=realtype), dimension(:, :), allocatable :: frac
      integer(kind=inttype), dimension(:, :), allocatable :: donorinfo
      integer(kind=inttype), dimension(:), allocatable :: procsizes, &
&     procdisps
      integer(kind=inttype), dimension(:), allocatable :: inv
      logical, dimension(:), allocatable :: valid
  end type usersurfcommtype
! two separate commes: one for the nodes (based on the primal
! mesh) and one for the variables (based on the dual mesh)
  type userintsurf
      character(len=maxstringlen) :: famname
      integer(kind=inttype) :: famid
      real(kind=realtype), dimension(:, :), allocatable :: pts
      integer(kind=inttype), dimension(:, :), allocatable :: conn
      type(usersurfcommtype) :: nodecomm, facecomm
  end type userintsurf
  integer(kind=inttype), parameter :: nuserintsurfsmax=25
  type(userintsurf), dimension(nuserintsurfsmax), target :: userintsurfs
  integer(kind=inttype) :: nuserintsurfs=0

contains
  subroutine integratesurfaces(localvalues, famlist)
! this is a shell routine that calls the specific surface
! integration routines. currently we have have the forceandmoment
! routine as well as the flow properties routine. this routine
! takes care of setting pointers, while the actual computational
! routine just acts on a specific fast pointed to by pointers. 
    use constants
    use blockpointers, only : nbocos, bcdata, bctype, sk, sj, si, x, &
&   rlv, sfacei, sfacej, sfacek, gamma, rev, p, viscsubface
    use utils_fast_b, only : setbcpointers, iswalltype
    use sorting, only : bsearchintegers
    use costfunctions, only : nlocalvalues
! tapenade needs to see these modules that the callees use.
    use bcpointers_fast_b
    use flowvarrefstate
    use inputphysics
    implicit none
! input/output variables
    real(kind=realtype), dimension(nlocalvalues), intent(inout) :: &
&   localvalues
    integer(kind=inttype), dimension(:), intent(in) :: famlist
! working variables
    integer(kind=inttype) :: mm
! loop over all possible boundary conditions
bocos:do mm=1,nbocos
! determine if this boundary condition is to be incldued in the
! currently active group
      if (bsearchintegers(bcdata(mm)%famid, famlist) .gt. 0) then
! set a bunch of pointers depending on the face id to make
! a generic treatment possible. 
        call setbcpointers(mm, .true.)
        if (iswalltype(bctype(mm))) call wallintegrationface(localvalues&
&                                                      , mm)
        if (bctype(mm) .eq. subsonicinflow .or. bctype(mm) .eq. &
&           supersonicinflow) then
          call flowintegrationface(.true., localvalues, mm)
        else if (bctype(mm) .eq. subsonicoutflow .or. bctype(mm) .eq. &
&           supersonicoutflow) then
          call flowintegrationface(.false., localvalues, mm)
        end if
      end if
    end do bocos
  end subroutine integratesurfaces
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
    use costfunctions
    use communication
    use blockpointers
    use flowvarrefstate
    use inputphysics, only : machcoef, pointref, veldirfreestream, &
&   equations
    use costfunctions, only : nlocalvalues, ifp, ifv, imp, imv, &
&   isepsensor, isepavg, icavitation, sepsensorsharpness, &
&   sepsensoroffset, iyplus
    use sorting, only : bsearchintegers
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
  subroutine flowintegrationface(isinflow, localvalues, mm)
    use constants
    use costfunctions
    use blockpointers, only : bctype, bcfaceid, bcdata, &
&   addgridvelocities
    use costfunctions, only : nlocalvalues, imassflow, imassptot, &
&   imassttot, imassps
    use sorting, only : bsearchintegers
    use flowvarrefstate, only : pref, pinf, rhoref, timeref, lref, &
&   tref, rgas
    use inputphysics, only : pointref, flowtype
    use flowutils_fast_b, only : computeptot, computettot
    use bcpointers_fast_b, only : ssi, sface, ww1, ww2, pp1, pp2, xx
    use utils_fast_b, only : mynorm2
    implicit none
! input/output variables
    logical, intent(in) :: isinflow
    real(kind=realtype), dimension(nlocalvalues), intent(inout) :: &
&   localvalues
    integer(kind=inttype), intent(in) :: mm
! local variables
    real(kind=realtype) :: massflowrate, mass_ptot, mass_ttot, mass_ps
    integer(kind=inttype) :: i, j, ii, blk
    real(kind=realtype) :: internalflowfact, fact, xc, yc, zc, cellarea&
&   , mx, my, mz
    real(kind=realtype) :: sf, vnm, vxm, vym, vzm, mredim, fx, fy, fz
    real(kind=realtype) :: pm, ptot, ttot, rhom, massflowratelocal
    real(kind=realtype), dimension(3) :: fp, mp, fmom, mmom, refpoint
    intrinsic sqrt
    intrinsic mod
    intrinsic max
    massflowrate = zero
    mass_ptot = zero
    mass_ttot = zero
    mass_ps = zero
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
! loop over the quadrilateral faces of the subface. note that
! the nodal range of bcdata must be used and not the cell
! range, because the latter may include the halo's in i and
! j-direction. the offset +1 is there, because inbeg and jnbeg
! refer to nodal ranges and not to cell ranges. the loop
! (without the ad stuff) would look like:
!
! do j=(bcdata(mm)%jnbeg+1),bcdata(mm)%jnend
!    do i=(bcdata(mm)%inbeg+1),bcdata(mm)%inend
    mredim = sqrt(pref*rhoref)
    fp = zero
    mp = zero
    fmom = zero
    mmom = zero
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
      vnm = vxm*ssi(i, j, 1) + vym*ssi(i, j, 2) + vzm*ssi(i, j, 3) - sf
      call computeptot(rhom, vxm, vym, vzm, pm, ptot)
      call computettot(rhom, vxm, vym, vzm, pm, ttot)
      pm = pm*pref
      massflowratelocal = rhom*vnm*mredim*blk*fact
      massflowrate = massflowrate + massflowratelocal
      mass_ptot = mass_ptot + ptot*massflowratelocal*pref
      mass_ttot = mass_ttot + ttot*massflowratelocal*tref
      mass_ps = mass_ps + pm*massflowratelocal
      xc = fourth*(xx(i, j, 1)+xx(i+1, j, 1)+xx(i, j+1, 1)+xx(i+1, j+1, &
&       1)) - refpoint(1)
      yc = fourth*(xx(i, j, 2)+xx(i+1, j, 2)+xx(i, j+1, 2)+xx(i+1, j+1, &
&       2)) - refpoint(2)
      zc = fourth*(xx(i, j, 3)+xx(i+1, j, 3)+xx(i, j+1, 3)+xx(i+1, j+1, &
&       3)) - refpoint(3)
! pressure forces. note that these need a *negative* and to subtract 
! the reference pressure sign to be consistent with the force 
! computation on the walls. 
      pm = -((pm-pinf*pref)*fact*blk)
! print *, "foo", mm, pm, pinf*pref
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
      select case  (bctype(mm)) 
      case (subsonicinflow, supersonicinflow) 
        massflowratelocal = -(massflowratelocal*fact/timeref*blk/&
&         cellarea*internalflowfact)
      case (subsonicoutflow, supersonicoutflow) 
        massflowratelocal = massflowratelocal*fact/timeref*blk/cellarea*&
&         internalflowfact
      end select
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
    end do
! increment the local values array with what we computed here
    localvalues(imassflow) = localvalues(imassflow) + massflowrate
    localvalues(imassptot) = localvalues(imassptot) + mass_ptot
    localvalues(imassttot) = localvalues(imassttot) + mass_ttot
    localvalues(imassps) = localvalues(imassps) + mass_ps
    localvalues(iflowfp:iflowfp+2) = localvalues(iflowfp:iflowfp+2) + fp
    localvalues(iflowfm:iflowfm+2) = localvalues(iflowfm:iflowfm+2) + &
&     fmom
    localvalues(iflowmp:iflowmp+2) = localvalues(iflowmp:iflowmp+2) + mp
    localvalues(iflowmm:iflowmm+2) = localvalues(iflowmm:iflowmm+2) + &
&     mmom
  end subroutine flowintegrationface
  subroutine flowintegrationzipper(isinflow, zipper, vars, localvalues, &
&   famlist, sps)
! integrate over the trianges for the inflow/outflow conditions. 
    use constants
    use costfunctions, only : nlocalvalues, imassflow, imassptot, &
&   imassttot, imassps, iflowmm, iflowmp, iflowfm, iflowfp
    use blockpointers, only : bctype
    use sorting, only : bsearchintegers
    use flowvarrefstate, only : pref, pinf, rhoref, pref, timeref, &
&   lref, tref, rgas
    use inputphysics, only : pointref, flowtype
    use flowutils_fast_b, only : computeptot, computettot
    use overset, only : zippermeshes, zippermesh
    use surfacefamilies, only : familyexchange, bcfamexchange
    use utils_fast_b, only : mynorm2, cross_prod
    implicit none
! input/output variables
    logical, intent(in) :: isinflow
    type(zippermesh), intent(in) :: zipper
    real(kind=realtype), dimension(:, :), intent(in) :: vars
    real(kind=realtype), dimension(nlocalvalues), intent(inout) :: &
&   localvalues
    integer(kind=inttype), dimension(:), intent(in) :: famlist
    integer(kind=inttype), intent(in) :: sps
! working variables
    integer(kind=inttype) :: i, j
    real(kind=realtype) :: sf, vnm, vxm, vym, vzm, mredim, fx, fy, fz
    real(kind=realtype), dimension(3) :: fp, mp, fmom, mmom, refpoint, &
&   ss, x1, x2, x3, norm
    real(kind=realtype) :: pm, ptot, ttot, rhom, massflowratelocal
    real(kind=realtype) :: massflowrate, mass_ptot, mass_ttot, mass_ps
    real(kind=realtype) :: internalflowfact, inflowfact, xc, yc, zc, &
&   cellarea, mx, my, mz
    real(kind=realtype), dimension(:), pointer :: localptr
    intrinsic sqrt
    intrinsic size
    real(kind=realtype), dimension(3) :: arg1
    real(kind=realtype), dimension(3) :: arg2
    real(kind=realtype) :: result1
    massflowrate = zero
    mass_ptot = zero
    mass_ttot = zero
    mass_ps = zero
    refpoint(1) = lref*pointref(1)
    refpoint(2) = lref*pointref(2)
    refpoint(3) = lref*pointref(3)
    mredim = sqrt(pref*rhoref)
    fp = zero
    mp = zero
    fmom = zero
    mmom = zero
    internalflowfact = one
    if (flowtype .eq. internalflow) internalflowfact = -one
    inflowfact = one
    if (isinflow) inflowfact = -one
    do i=1,size(zipper%conn, 2)
      if (bsearchintegers(zipper%fam(i), famlist) .gt. 0) then
! compute the averaged values for this trianlge
        vxm = zero
        vym = zero
        vzm = zero
        rhom = zero
        pm = zero
        sf = zero
        do j=1,3
          rhom = rhom + vars(zipper%conn(j, i), irho)
          vxm = vxm + vars(zipper%conn(j, i), ivx)
          vym = vym + vars(zipper%conn(j, i), ivy)
          vzm = vzm + vars(zipper%conn(j, i), ivz)
          pm = pm + vars(zipper%conn(j, i), irhoe)
          sf = sf + vars(zipper%conn(j, i), 6)
        end do
! divide by 3 due to the summation above:
        rhom = third*rhom
        vxm = third*vxm
        vym = third*vym
        vzm = third*vzm
        pm = third*pm
        sf = third*sf
! get the nodes of triangle.
        x1 = vars(zipper%conn(1, i), 7:9)
        x2 = vars(zipper%conn(2, i), 7:9)
        x3 = vars(zipper%conn(3, i), 7:9)
        arg1(:) = x2 - x1
        arg2(:) = x3 - x1
        call cross_prod(arg1(:), arg2(:), norm)
        ss = -(half*norm)
        call computeptot(rhom, vxm, vym, vzm, pm, ptot)
        call computettot(rhom, vxm, vym, vzm, pm, ttot)
        pm = pm*pref
        vnm = vxm*ss(1) + vym*ss(2) + vzm*ss(3) - sf
        massflowratelocal = rhom*vnm*mredim
        massflowrate = massflowrate + massflowratelocal
        mass_ptot = mass_ptot + ptot*massflowratelocal*pref
        mass_ttot = mass_ttot + ttot*massflowratelocal*tref
        mass_ps = mass_ps + pm*massflowratelocal
! compute the average cell center. 
        xc = zero
        yc = zero
        zc = zero
        do j=1,3
          xc = xc + vars(zipper%conn(1, i), 7)
          yc = yc + vars(zipper%conn(2, i), 8)
          zc = zc + vars(zipper%conn(3, i), 9)
        end do
! finish average for cell center
        xc = third*xc
        yc = third*yc
        zc = third*zc
        xc = xc - refpoint(1)
        yc = yc - refpoint(2)
        zc = zc - refpoint(3)
        pm = -(pm-pinf*pref)
        fx = pm*ss(1)
        fy = pm*ss(2)
        fz = pm*ss(3)
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
! momentum forces
! get unit normal vector. 
        result1 = mynorm2(ss)
        ss = ss/result1
        massflowratelocal = massflowratelocal/timeref*internalflowfact*&
&         inflowfact
        fx = massflowratelocal*ss(1)*vxm/timeref
        fy = massflowratelocal*ss(2)*vym/timeref
        fz = massflowratelocal*ss(3)*vzm/timeref
! note: momentum forces have opposite sign to pressure forces
        fmom(1) = fmom(1) - fx
        fmom(2) = fmom(2) - fy
        fmom(3) = fmom(3) - fz
        mx = yc*fz - zc*fy
        my = zc*fx - xc*fz
        mz = xc*fy - yc*fx
        mmom(1) = mmom(1) + mx
        mmom(2) = mmom(2) + my
        mmom(3) = mmom(3) + mz
      end if
    end do
! increment the local values array with what we computed here
    localvalues(imassflow) = localvalues(imassflow) + massflowrate
    localvalues(imassptot) = localvalues(imassptot) + mass_ptot
    localvalues(imassttot) = localvalues(imassttot) + mass_ttot
    localvalues(imassps) = localvalues(imassps) + mass_ps
    localvalues(iflowfp:iflowfp+2) = localvalues(iflowfp:iflowfp+2) + fp
    localvalues(iflowfm:iflowfm+2) = localvalues(iflowfm:iflowfm+2) + &
&     fmom
    localvalues(iflowmp:iflowmp+2) = localvalues(iflowmp:iflowmp+2) + mp
    localvalues(iflowmm:iflowmm+2) = localvalues(iflowmm:iflowmm+2) + &
&     mmom
  end subroutine flowintegrationzipper
  subroutine wallintegrationzipper(zipper, vars, localvalues, famlist, &
&   sps)
    use constants
    use costfunctions
    use sorting, only : bsearchintegers
    use flowvarrefstate, only : lref
    use inputphysics, only : pointref
    use overset, only : zippermeshes, zippermesh
    use utils_fast_b, only : mynorm2, cross_prod
    implicit none
! input/output
    type(zippermesh), intent(in) :: zipper
    real(kind=realtype), dimension(:, :), intent(in) :: vars
    real(kind=realtype), intent(inout) :: localvalues(nlocalvalues)
    integer(kind=inttype), dimension(:), intent(in) :: famlist
    integer(kind=inttype), intent(in) :: sps
! working
    real(kind=realtype), dimension(3) :: fp, fv, mp, mv
    integer(kind=inttype) :: i, j
    real(kind=realtype), dimension(3) :: ss, norm, refpoint
    real(kind=realtype), dimension(3) :: p1, p2, p3, v1, v2, v3, x1, x2&
&   , x3
    real(kind=realtype) :: fact, triarea, fx, fy, fz, mx, my, mz, xc, yc&
&   , zc
    intrinsic size
    real(kind=realtype), dimension(3) :: arg1
    real(kind=realtype), dimension(3) :: arg2
    real(kind=realtype) :: result1
! determine the reference point for the moment computation in
! meters.
    refpoint(1) = lref*pointref(1)
    refpoint(2) = lref*pointref(2)
    refpoint(3) = lref*pointref(3)
    fp = zero
    fv = zero
    mp = zero
    mv = zero
    do i=1,size(zipper%conn, 2)
      if (bsearchintegers(zipper%fam(i), famlist) .gt. 0) then
! get the nodes of triangle. the *3 is becuase of the
! blanket third above. 
        x1 = vars(zipper%conn(1, i), 7:9)
        x2 = vars(zipper%conn(2, i), 7:9)
        x3 = vars(zipper%conn(3, i), 7:9)
        arg1(:) = x2 - x1
        arg2(:) = x3 - x1
        call cross_prod(arg1(:), arg2(:), norm)
        ss = half*norm
! the third here is to account for the summation of p1, p2
! and p3
        result1 = mynorm2(ss)
        triarea = result1*third
! compute the average cell center. 
        xc = third*(x1(1)+x2(1)+x3(1))
        xc = third*(x1(2)+x2(2)+x3(2))
        xc = third*(x1(3)+x2(3)+x3(3))
        xc = xc - refpoint(1)
        yc = yc - refpoint(2)
        zc = zc - refpoint(3)
! update the pressure force and moment coefficients.
        p1 = vars(zipper%conn(1, i), 1:3)
        p2 = vars(zipper%conn(2, i), 1:3)
        p3 = vars(zipper%conn(3, i), 1:3)
        fx = (p1(1)+p2(1)+p3(1))*triarea
        fy = (p1(2)+p2(2)+p3(2))*triarea
        fz = (p1(3)+p2(3)+p3(3))*triarea
        fp(1) = fp(1) + fx
        fp(2) = fp(2) + fy
        fp(3) = fp(3) + fz
        mx = yc*fz - zc*fy
        my = zc*fx - xc*fz
        mz = xc*fy - yc*fx
        mp(1) = mp(1) + mx
        mp(2) = mp(2) + my
        mp(3) = mp(3) + mz
! update the viscous force and moment coefficients
        v1 = vars(zipper%conn(1, i), 4:6)
        v2 = vars(zipper%conn(2, i), 4:6)
        v3 = vars(zipper%conn(3, i), 4:6)
        fx = (v1(1)+v2(1)+v3(1))*triarea
        fy = (v1(2)+v2(2)+v3(2))*triarea
        fz = (v1(3)+v2(3)+v3(3))*triarea
! note: momentum forces have opposite sign to pressure forces
        fv(1) = fv(1) + fx
        fv(2) = fv(2) + fy
        fv(3) = fv(3) + fz
        mx = yc*fz - zc*fy
        my = zc*fx - xc*fz
        mz = xc*fy - yc*fx
        mv(1) = mv(1) + mx
        mv(2) = mv(2) + my
        mv(3) = mv(3) + mz
      end if
    end do
! increment into the local vector
    localvalues(ifp:ifp+2) = localvalues(ifp:ifp+2) + fp
    localvalues(ifv:ifv+2) = localvalues(ifv:ifv+2) + fv
    localvalues(imp:imp+2) = localvalues(imp:imp+2) + mp
    localvalues(imv:imv+2) = localvalues(imv:imv+2) + mv
  end subroutine wallintegrationzipper
end module surfaceintegrations_fast_b
