!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of slipvelocitiesfinelevel_block in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: *(*bcdata.uslip)
!   with respect to varying inputs: gammainf pinf timeref rhoinf
!                *x veldirfreestream machgrid
!   plus diff mem management of: x:in bcdata:in *bcdata.uslip:in
!
!      ******************************************************************
!      *                                                                *
!      * file:          slipvelocities.f90                              *
!      * author:        edwin van der weide                             *
!      * starting date: 02-12-2004                                      *
!      * last modified: 06-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine slipvelocitiesfinelevel_block_d(useoldcoor, t, sps)
!
!      ******************************************************************
!      *                                                                *
!      * slipvelocitiesfinelevel computes the slip velocities for       *
!      * viscous subfaces on all viscous boundaries on groundlevel for  *
!      * the given spectral solution. if useoldcoor is .true. the       *
!      * velocities are determined using the unsteady time integrator;  *
!      * otherwise the analytic form is used.                           *
!      *                                                                *
!      ******************************************************************
!
  use bctypes
  use inputtimespectral
  use blockpointers
  use cgnsgrid
  use flowvarrefstate
  use inputmotion
  use inputunsteady
  use iteration
  use inputphysics
  use inputtsstabderiv
  use monitor
  use communication
  use diffsizes
!  hint: isize1ofdrfbcdata should be the size of dimension 1 of array *bcdata
  implicit none
!
!      subroutine arguments.
!
  integer(kind=inttype), intent(in) :: sps
  logical, intent(in) :: useoldcoor
  real(kind=realtype), dimension(*), intent(in) :: t
!
!      local variables.
!
  integer(kind=inttype) :: nn, mm, i, j, level
  real(kind=realtype) :: oneover4dt
  real(kind=realtype) :: oneover4dtd
  real(kind=realtype) :: velxgrid, velygrid, velzgrid, ainf
  real(kind=realtype) :: velxgridd, velygridd, velzgridd, ainfd
  real(kind=realtype) :: velxgrid0, velygrid0, velzgrid0
  real(kind=realtype) :: velxgrid0d, velygrid0d, velzgrid0d
  real(kind=realtype), dimension(3) :: xc, xxc
  real(kind=realtype), dimension(3) :: xcd, xxcd
  real(kind=realtype), dimension(3) :: rotcenter, rotrate
  real(kind=realtype), dimension(3) :: rotrated
  real(kind=realtype), dimension(3) :: rotationpoint
  real(kind=realtype), dimension(3, 3) :: rotationmatrix, &
& derivrotationmatrix
  real(kind=realtype), dimension(3, 3) :: derivrotationmatrixd
  real(kind=realtype) :: tnew, told
  real(kind=realtype), dimension(:, :, :), pointer :: uslip
  real(kind=realtype), dimension(:, :, :), pointer :: uslipd
  real(kind=realtype), dimension(:, :, :), pointer :: xface
  real(kind=realtype), dimension(:, :, :), pointer :: xfaced
  real(kind=realtype), dimension(:, :, :, :), pointer :: xfaceold
  integer(kind=inttype) :: liftindex
  real(kind=realtype) :: alpha, beta, intervalmach, alphats, &
& alphaincrement, betats, betaincrement
  real(kind=realtype) :: alphad, betad, alphatsd, betatsd
  real(kind=realtype), dimension(3) :: veldir
  real(kind=realtype), dimension(3) :: veldird
  real(kind=realtype), dimension(3) :: refdirection
!function definitions
  real(kind=realtype) :: tsalpha, tsbeta, tsmach
  intrinsic sqrt
  real(kind=realtype) :: arg1
  real(kind=realtype) :: arg1d
  integer :: ii1
!
!      ******************************************************************
!      *                                                                *
!      * begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
! determine the situation we are having here.
  if (useoldcoor) then
! the velocities must be determined via a finite difference
! formula using the coordinates of the old levels.
! set the coefficients for the time integrator and store the
! inverse of the physical nondimensional time step, divided
! by 4, a bit easier.
    call setcoeftimeintegrator()
    oneover4dtd = fourth*timerefd/deltat
    oneover4dt = fourth*timeref/deltat
    do ii1=1,isize1ofdrfbcdata
      bcdatad(ii1)%uslip = 0.0_8
    end do
    xcd = 0.0_8
! loop over the number of viscous subfaces.
bocoloop1:do mm=1,nviscbocos
! set the pointer for uslip to make the code more
! readable.
      uslipd => bcdatad(mm)%uslip
      uslip => bcdata(mm)%uslip
! determine the grid face on which the subface is located
! and set some variables accordingly.
      select case  (bcfaceid(mm)) 
      case (imin) 
        xfaced => xd(1, :, :, :)
        xface => x(1, :, :, :)
        xfaceold => xold(:, 1, :, :, :)
      case (imax) 
        xfaced => xd(il, :, :, :)
        xface => x(il, :, :, :)
        xfaceold => xold(:, il, :, :, :)
      case (jmin) 
        xfaced => xd(:, 1, :, :)
        xface => x(:, 1, :, :)
        xfaceold => xold(:, :, 1, :, :)
      case (jmax) 
        xfaced => xd(:, jl, :, :)
        xface => x(:, jl, :, :)
        xfaceold => xold(:, :, jl, :, :)
      case (kmin) 
        xfaced => xd(:, :, 1, :)
        xface => x(:, :, 1, :)
        xfaceold => xold(:, :, :, 1, :)
      case (kmax) 
        xfaced => xd(:, :, kl, :)
        xface => x(:, :, kl, :)
        xfaceold => xold(:, :, :, kl, :)
      end select
! some boundary faces have a different rotation speed than
! the corresponding block. this happens e.g. in the tip gap
! region of turbomachinary problems where the casing does
! not rotate. as the coordinate difference corresponds to
! the rotation rate of the block, a correction must be
! computed. therefore compute the difference in rotation
! rate and store the rotation center a bit easier. note that
! the rotation center of subface is taken, because if there
! is a difference in rotation rate this info for the subface
! must always be specified.
      j = nbkglobal
      i = cgnssubface(mm)
      rotcenter = cgnsdoms(j)%bocoinfo(i)%rotcenter
      rotrated = (cgnsdoms(j)%bocoinfo(i)%rotrate-cgnsdoms(j)%rotrate)*&
&       timerefd
      rotrate = timeref*(cgnsdoms(j)%bocoinfo(i)%rotrate-cgnsdoms(j)%&
&       rotrate)
! loop over the quadrilateral faces of the viscous subface.
! note that due to the usage of the pointers xface and
! xfaceold an offset of +1 must be used in the coordinate
! arrays, because x and xold originally start at 0 for the
! i, j and k indices.
      do j=bcdata(mm)%jcbeg,bcdata(mm)%jcend
        do i=bcdata(mm)%icbeg,bcdata(mm)%icend
! determine the coordinates of the centroid of the
! face, multiplied by 4.
          xcd(1) = xfaced(i+1, j+1, 1) + xfaced(i+1, j, 1) + xfaced(i, j&
&           +1, 1) + xfaced(i, j, 1)
          xc(1) = xface(i+1, j+1, 1) + xface(i+1, j, 1) + xface(i, j+1, &
&           1) + xface(i, j, 1)
          xcd(2) = xfaced(i+1, j+1, 2) + xfaced(i+1, j, 2) + xfaced(i, j&
&           +1, 2) + xfaced(i, j, 2)
          xc(2) = xface(i+1, j+1, 2) + xface(i+1, j, 2) + xface(i, j+1, &
&           2) + xface(i, j, 2)
          xcd(3) = xfaced(i+1, j+1, 3) + xfaced(i+1, j, 3) + xfaced(i, j&
&           +1, 3) + xfaced(i, j, 3)
          xc(3) = xface(i+1, j+1, 3) + xface(i+1, j, 3) + xface(i, j+1, &
&           3) + xface(i, j, 3)
! multiply the sum of the 4 vertex coordinates with
! coeftime(0) to obtain the contribution for the
! current time level. the division by 4*deltat will
! take place later. this is both more efficient and
! more accurate for extremely small time steps.
          uslipd(i, j, 1) = coeftime(0)*xcd(1)
          uslip(i, j, 1) = coeftime(0)*xc(1)
          uslipd(i, j, 2) = coeftime(0)*xcd(2)
          uslip(i, j, 2) = coeftime(0)*xc(2)
          uslipd(i, j, 3) = coeftime(0)*xcd(3)
          uslip(i, j, 3) = coeftime(0)*xc(3)
! loop over the older time levels and take their
! contribution into account.
          do level=1,noldlevels
            uslip(i, j, 1) = uslip(i, j, 1) + coeftime(level)*(xfaceold(&
&             level, i+1, j+1, 1)+xfaceold(level, i+1, j, 1)+xfaceold(&
&             level, i, j+1, 1)+xfaceold(level, i, j, 1))
            uslip(i, j, 2) = uslip(i, j, 2) + coeftime(level)*(xfaceold(&
&             level, i+1, j+1, 2)+xfaceold(level, i+1, j, 2)+xfaceold(&
&             level, i, j+1, 2)+xfaceold(level, i, j, 2))
            uslip(i, j, 3) = uslip(i, j, 3) + coeftime(level)*(xfaceold(&
&             level, i+1, j+1, 3)+xfaceold(level, i+1, j, 3)+xfaceold(&
&             level, i, j+1, 3)+xfaceold(level, i, j, 3))
          end do
! divide by 4 times the time step to obtain the
! correct velocity.
          uslipd(i, j, 1) = uslipd(i, j, 1)*oneover4dt + uslip(i, j, 1)*&
&           oneover4dtd
          uslip(i, j, 1) = uslip(i, j, 1)*oneover4dt
          uslipd(i, j, 2) = uslipd(i, j, 2)*oneover4dt + uslip(i, j, 2)*&
&           oneover4dtd
          uslip(i, j, 2) = uslip(i, j, 2)*oneover4dt
          uslipd(i, j, 3) = uslipd(i, j, 3)*oneover4dt + uslip(i, j, 3)*&
&           oneover4dtd
          uslip(i, j, 3) = uslip(i, j, 3)*oneover4dt
! determine the correction due to the difference
! in rotation rate between the block and subface.
! first determine the coordinates relative to the
! rotation center. remember that 4 times this value
! is currently stored in xc.
          xcd(1) = fourth*xcd(1)
          xc(1) = fourth*xc(1) - rotcenter(1)
          xcd(2) = fourth*xcd(2)
          xc(2) = fourth*xc(2) - rotcenter(2)
          xcd(3) = fourth*xcd(3)
          xc(3) = fourth*xc(3) - rotcenter(3)
! compute the velocity, which is the cross product
! of rotrate and xc and add it to uslip.
          uslipd(i, j, 1) = uslipd(i, j, 1) + rotrated(2)*xc(3) + &
&           rotrate(2)*xcd(3) - rotrated(3)*xc(2) - rotrate(3)*xcd(2)
          uslip(i, j, 1) = uslip(i, j, 1) + rotrate(2)*xc(3) - rotrate(3&
&           )*xc(2)
          uslipd(i, j, 2) = uslipd(i, j, 2) + rotrated(3)*xc(1) + &
&           rotrate(3)*xcd(1) - rotrated(1)*xc(3) - rotrate(1)*xcd(3)
          uslip(i, j, 2) = uslip(i, j, 2) + rotrate(3)*xc(1) - rotrate(1&
&           )*xc(3)
          uslipd(i, j, 3) = uslipd(i, j, 3) + rotrated(1)*xc(2) + &
&           rotrate(1)*xcd(2) - rotrated(2)*xc(1) - rotrate(2)*xcd(1)
          uslip(i, j, 3) = uslip(i, j, 3) + rotrate(1)*xc(2) - rotrate(2&
&           )*xc(1)
        end do
      end do
    end do bocoloop1
  else
! the velocities must be determined analytically.
! compute the mesh velocity from the given mesh mach number.
!  ainf = sqrt(gammainf*pinf/rhoinf)
!  velxgrid = ainf*machgrid(1)
!  velygrid = ainf*machgrid(2)
!  velzgrid = ainf*machgrid(3)
    arg1d = ((gammainfd*pinf+gammainf*pinfd)*rhoinf-gammainf*pinf*&
&     rhoinfd)/rhoinf**2
    arg1 = gammainf*pinf/rhoinf
    if (arg1 .eq. 0.0_8) then
      ainfd = 0.0_8
    else
      ainfd = arg1d/(2.0*sqrt(arg1))
    end if
    ainf = sqrt(arg1)
    velxgrid0d = -((ainfd*machgrid+ainf*machgridd)*veldirfreestream(1)) &
&     - ainf*machgrid*veldirfreestreamd(1)
    velxgrid0 = ainf*machgrid*(-veldirfreestream(1))
    velygrid0d = -((ainfd*machgrid+ainf*machgridd)*veldirfreestream(2)) &
&     - ainf*machgrid*veldirfreestreamd(2)
    velygrid0 = ainf*machgrid*(-veldirfreestream(2))
    velzgrid0d = -((ainfd*machgrid+ainf*machgridd)*veldirfreestream(3)) &
&     - ainf*machgrid*veldirfreestreamd(3)
    velzgrid0 = ainf*machgrid*(-veldirfreestream(3))
! compute the derivative of the rotation matrix and the rotation
! point; needed for velocity due to the rigid body rotation of
! the entire grid. it is assumed that the rigid body motion of
! the grid is only specified if there is only 1 section present.
    call derivativerotmatrixrigid_d(derivrotationmatrix, &
&                             derivrotationmatrixd, rotationpoint, t(1))
!compute the rotation matrix to update the velocities for the time
!spectral stability derivative case...
    if (tsstability) then
! determine the time values of the old and new time level.
! it is assumed that the rigid body rotation of the mesh is only
! used when only 1 section is present.
      tnew = timeunsteady + timeunsteadyrestart
      told = tnew - t(1)
      if ((tspmode .or. tsqmode) .or. tsrmode) then
! compute the rotation matrix of the rigid body rotation as
! well as the rotation point; the latter may vary in time due
! to rigid body translation.
        call rotmatrixrigidbody(tnew, told, rotationmatrix, &
&                            rotationpoint)
        velxgrid0d = rotationmatrix(1, 1)*velxgrid0d + rotationmatrix(1&
&         , 2)*velygrid0d + rotationmatrix(1, 3)*velzgrid0d
        velxgrid0 = rotationmatrix(1, 1)*velxgrid0 + rotationmatrix(1, 2&
&         )*velygrid0 + rotationmatrix(1, 3)*velzgrid0
        velygrid0d = rotationmatrix(2, 1)*velxgrid0d + rotationmatrix(2&
&         , 2)*velygrid0d + rotationmatrix(2, 3)*velzgrid0d
        velygrid0 = rotationmatrix(2, 1)*velxgrid0 + rotationmatrix(2, 2&
&         )*velygrid0 + rotationmatrix(2, 3)*velzgrid0
        velzgrid0d = rotationmatrix(3, 1)*velxgrid0d + rotationmatrix(3&
&         , 2)*velygrid0d + rotationmatrix(3, 3)*velzgrid0d
        velzgrid0 = rotationmatrix(3, 1)*velxgrid0 + rotationmatrix(3, 2&
&         )*velygrid0 + rotationmatrix(3, 3)*velzgrid0
        do ii1=1,isize1ofdrfbcdata
          bcdatad(ii1)%uslip = 0.0_8
        end do
        xcd = 0.0_8
        xxcd = 0.0_8
      else if (tsalphamode) then
! get the baseline alpha and determine the liftindex
        call getdirangle_d(veldirfreestream, veldirfreestreamd, &
&                    liftdirection, liftindex, alpha, alphad, beta, &
&                    betad)
!determine the alpha for this time instance
        alphaincrement = tsalpha(degreepolalpha, coefpolalpha, &
&         degreefouralpha, omegafouralpha, coscoeffouralpha, &
&         sincoeffouralpha, t(1))
        alphatsd = alphad
        alphats = alpha + alphaincrement
!determine the grid velocity for this alpha
        refdirection(:) = zero
        refdirection(1) = one
        call getdirvector_d(refdirection, alphats, alphatsd, beta, betad&
&                     , veldir, veldird, liftindex)
!do i need to update the lift direction and drag direction as well?
!set the effictive grid velocity for this time interval
        velxgrid0d = -((ainfd*machgrid+ainf*machgridd)*veldir(1)) - ainf&
&         *machgrid*veldird(1)
        velxgrid0 = ainf*machgrid*(-veldir(1))
        velygrid0d = -((ainfd*machgrid+ainf*machgridd)*veldir(2)) - ainf&
&         *machgrid*veldird(2)
        velygrid0 = ainf*machgrid*(-veldir(2))
        velzgrid0d = -((ainfd*machgrid+ainf*machgridd)*veldir(3)) - ainf&
&         *machgrid*veldird(3)
        velzgrid0 = ainf*machgrid*(-veldir(3))
        do ii1=1,isize1ofdrfbcdata
          bcdatad(ii1)%uslip = 0.0_8
        end do
        xcd = 0.0_8
        xxcd = 0.0_8
      else if (tsbetamode) then
! get the baseline alpha and determine the liftindex
        call getdirangle_d(veldirfreestream, veldirfreestreamd, &
&                    liftdirection, liftindex, alpha, alphad, beta, &
&                    betad)
!determine the alpha for this time instance
        betaincrement = tsbeta(degreepolbeta, coefpolbeta, &
&         degreefourbeta, omegafourbeta, coscoeffourbeta, &
&         sincoeffourbeta, t(1))
        betatsd = betad
        betats = beta + betaincrement
!determine the grid velocity for this alpha
        refdirection(:) = zero
        refdirection(1) = one
        call getdirvector_d(refdirection, alpha, alphad, betats, betatsd&
&                     , veldir, veldird, liftindex)
!do i need to update the lift direction and drag direction as well?
!set the effictive grid velocity for this time interval
        velxgrid0d = -((ainfd*machgrid+ainf*machgridd)*veldir(1)) - ainf&
&         *machgrid*veldird(1)
        velxgrid0 = ainf*machgrid*(-veldir(1))
        velygrid0d = -((ainfd*machgrid+ainf*machgridd)*veldir(2)) - ainf&
&         *machgrid*veldird(2)
        velygrid0 = ainf*machgrid*(-veldir(2))
        velzgrid0d = -((ainfd*machgrid+ainf*machgridd)*veldir(3)) - ainf&
&         *machgrid*veldird(3)
        velzgrid0 = ainf*machgrid*(-veldir(3))
        do ii1=1,isize1ofdrfbcdata
          bcdatad(ii1)%uslip = 0.0_8
        end do
        xcd = 0.0_8
        xxcd = 0.0_8
      else if (tsmachmode) then
!determine the mach number at this time interval
        intervalmach = tsmach(degreepolmach, coefpolmach, &
&         degreefourmach, omegafourmach, coscoeffourmach, &
&         sincoeffourmach, t(1))
!set the effective grid velocity
        velxgrid0d = -((ainfd*(intervalmach+machgrid)+ainf*machgridd)*&
&         veldirfreestream(1)) - ainf*(intervalmach+machgrid)*&
&         veldirfreestreamd(1)
        velxgrid0 = ainf*(intervalmach+machgrid)*(-veldirfreestream(1))
        velygrid0d = -((ainfd*(intervalmach+machgrid)+ainf*machgridd)*&
&         veldirfreestream(2)) - ainf*(intervalmach+machgrid)*&
&         veldirfreestreamd(2)
        velygrid0 = ainf*(intervalmach+machgrid)*(-veldirfreestream(2))
        velzgrid0d = -((ainfd*(intervalmach+machgrid)+ainf*machgridd)*&
&         veldirfreestream(3)) - ainf*(intervalmach+machgrid)*&
&         veldirfreestreamd(3)
        velzgrid0 = ainf*(intervalmach+machgrid)*(-veldirfreestream(3))
        do ii1=1,isize1ofdrfbcdata
          bcdatad(ii1)%uslip = 0.0_8
        end do
        xcd = 0.0_8
        xxcd = 0.0_8
      else if (tsaltitudemode) then
        call returnFail('gridvelocityfinelevel', &
&                   'altitude motion not yet implemented...')
        do ii1=1,isize1ofdrfbcdata
          bcdatad(ii1)%uslip = 0.0_8
        end do
        xcd = 0.0_8
        xxcd = 0.0_8
      else
        call returnFail('gridvelocityfinelevel', &
&                   'not a recognized stability motion')
        do ii1=1,isize1ofdrfbcdata
          bcdatad(ii1)%uslip = 0.0_8
        end do
        xcd = 0.0_8
        xxcd = 0.0_8
      end if
    else
      do ii1=1,isize1ofdrfbcdata
        bcdatad(ii1)%uslip = 0.0_8
      end do
      xcd = 0.0_8
      xxcd = 0.0_8
    end if
! loop over the number of viscous subfaces.
bocoloop2:do mm=1,nviscbocos
! determine the grid face on which the subface is located
! and set some variables accordingly.
      select case  (bcfaceid(mm)) 
      case (imin) 
        xfaced => xd(1, :, :, :)
        xface => x(1, :, :, :)
      case (imax) 
        xfaced => xd(il, :, :, :)
        xface => x(il, :, :, :)
      case (jmin) 
        xfaced => xd(:, 1, :, :)
        xface => x(:, 1, :, :)
      case (jmax) 
        xfaced => xd(:, jl, :, :)
        xface => x(:, jl, :, :)
      case (kmin) 
        xfaced => xd(:, :, 1, :)
        xface => x(:, :, 1, :)
      case (kmax) 
        xfaced => xd(:, :, kl, :)
        xface => x(:, :, kl, :)
      end select
! store the rotation center and the rotation rate
! for this subface.
      j = nbkglobal
      i = cgnssubface(mm)
      rotcenter = cgnsdoms(j)%bocoinfo(i)%rotcenter
      rotrated = cgnsdoms(j)%bocoinfo(i)%rotrate*timerefd
      rotrate = timeref*cgnsdoms(j)%bocoinfo(i)%rotrate
! usewindaxis should go back here!
      velxgridd = velxgrid0d
      velxgrid = velxgrid0
      velygridd = velygrid0d
      velygrid = velygrid0
      velzgridd = velzgrid0d
      velzgrid = velzgrid0
! loop over the quadrilateral faces of the viscous
! subface.
      do j=bcdata(mm)%jcbeg,bcdata(mm)%jcend
        do i=bcdata(mm)%icbeg,bcdata(mm)%icend
! compute the coordinates of the centroid of the face.
! normally this would be an average of i-1 and i, but
! due to the usage of the pointer xface and the fact
! that x starts at index 0 this is shifted 1 index.
          xcd(1) = fourth*(xfaced(i+1, j+1, 1)+xfaced(i+1, j, 1)+xfaced(&
&           i, j+1, 1)+xfaced(i, j, 1))
          xc(1) = fourth*(xface(i+1, j+1, 1)+xface(i+1, j, 1)+xface(i, j&
&           +1, 1)+xface(i, j, 1))
          xcd(2) = fourth*(xfaced(i+1, j+1, 2)+xfaced(i+1, j, 2)+xfaced(&
&           i, j+1, 2)+xfaced(i, j, 2))
          xc(2) = fourth*(xface(i+1, j+1, 2)+xface(i+1, j, 2)+xface(i, j&
&           +1, 2)+xface(i, j, 2))
          xcd(3) = fourth*(xfaced(i+1, j+1, 3)+xfaced(i+1, j, 3)+xfaced(&
&           i, j+1, 3)+xfaced(i, j, 3))
          xc(3) = fourth*(xface(i+1, j+1, 3)+xface(i+1, j, 3)+xface(i, j&
&           +1, 3)+xface(i, j, 3))
! determine the coordinates relative to the center
! of rotation.
          xxcd(1) = xcd(1)
          xxc(1) = xc(1) - rotcenter(1)
          xxcd(2) = xcd(2)
          xxc(2) = xc(2) - rotcenter(2)
          xxcd(3) = xcd(3)
          xxc(3) = xc(3) - rotcenter(3)
! compute the velocity, which is the cross product
! of rotrate and xc.
          bcdatad(mm)%uslip(i, j, 1) = rotrated(2)*xxc(3) + rotrate(2)*&
&           xxcd(3) - rotrated(3)*xxc(2) - rotrate(3)*xxcd(2)
          bcdata(mm)%uslip(i, j, 1) = rotrate(2)*xxc(3) - rotrate(3)*xxc&
&           (2)
          bcdatad(mm)%uslip(i, j, 2) = rotrated(3)*xxc(1) + rotrate(3)*&
&           xxcd(1) - rotrated(1)*xxc(3) - rotrate(1)*xxcd(3)
          bcdata(mm)%uslip(i, j, 2) = rotrate(3)*xxc(1) - rotrate(1)*xxc&
&           (3)
          bcdatad(mm)%uslip(i, j, 3) = rotrated(1)*xxc(2) + rotrate(1)*&
&           xxcd(2) - rotrated(2)*xxc(1) - rotrate(2)*xxcd(1)
          bcdata(mm)%uslip(i, j, 3) = rotrate(1)*xxc(2) - rotrate(2)*xxc&
&           (1)
! determine the coordinates relative to the
! rigid body rotation point.
          xxcd(1) = xcd(1)
          xxc(1) = xc(1) - rotationpoint(1)
          xxcd(2) = xcd(2)
          xxc(2) = xc(2) - rotationpoint(2)
          xxcd(3) = xcd(3)
          xxc(3) = xc(3) - rotationpoint(3)
! determine the total velocity of the cell center.
! this is a combination of rotation speed of this
! block and the entire rigid body rotation.
          bcdatad(mm)%uslip(i, j, 1) = bcdatad(mm)%uslip(i, j, 1) + &
&           velxgridd + derivrotationmatrixd(1, 1)*xxc(1) + &
&           derivrotationmatrix(1, 1)*xxcd(1) + derivrotationmatrixd(1, &
&           2)*xxc(2) + derivrotationmatrix(1, 2)*xxcd(2) + &
&           derivrotationmatrixd(1, 3)*xxc(3) + derivrotationmatrix(1, 3&
&           )*xxcd(3)
          bcdata(mm)%uslip(i, j, 1) = bcdata(mm)%uslip(i, j, 1) + &
&           velxgrid + derivrotationmatrix(1, 1)*xxc(1) + &
&           derivrotationmatrix(1, 2)*xxc(2) + derivrotationmatrix(1, 3)&
&           *xxc(3)
          bcdatad(mm)%uslip(i, j, 2) = bcdatad(mm)%uslip(i, j, 2) + &
&           velygridd + derivrotationmatrixd(2, 1)*xxc(1) + &
&           derivrotationmatrix(2, 1)*xxcd(1) + derivrotationmatrixd(2, &
&           2)*xxc(2) + derivrotationmatrix(2, 2)*xxcd(2) + &
&           derivrotationmatrixd(2, 3)*xxc(3) + derivrotationmatrix(2, 3&
&           )*xxcd(3)
          bcdata(mm)%uslip(i, j, 2) = bcdata(mm)%uslip(i, j, 2) + &
&           velygrid + derivrotationmatrix(2, 1)*xxc(1) + &
&           derivrotationmatrix(2, 2)*xxc(2) + derivrotationmatrix(2, 3)&
&           *xxc(3)
          bcdatad(mm)%uslip(i, j, 3) = bcdatad(mm)%uslip(i, j, 3) + &
&           velzgridd + derivrotationmatrixd(3, 1)*xxc(1) + &
&           derivrotationmatrix(3, 1)*xxcd(1) + derivrotationmatrixd(3, &
&           2)*xxc(2) + derivrotationmatrix(3, 2)*xxcd(2) + &
&           derivrotationmatrixd(3, 3)*xxc(3) + derivrotationmatrix(3, 3&
&           )*xxcd(3)
          bcdata(mm)%uslip(i, j, 3) = bcdata(mm)%uslip(i, j, 3) + &
&           velzgrid + derivrotationmatrix(3, 1)*xxc(1) + &
&           derivrotationmatrix(3, 2)*xxc(2) + derivrotationmatrix(3, 3)&
&           *xxc(3)
        end do
      end do
    end do bocoloop2
  end if
end subroutine slipvelocitiesfinelevel_block_d
