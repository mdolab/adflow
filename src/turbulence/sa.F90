! This module contains the source code related to the SA turbulence
! model. It is slightly more modularized than the original which makes
! performing reverse mode AD simplier.

module sa

  use constants
  real(kind=realType) :: cv13, kar2Inv, cw36, cb3Inv
  real(kind=realType), dimension(:,:,:), allocatable :: qq
  real(kind=realType), dimension(:,:,:), pointer :: ddw, ww, ddvt
  real(kind=realType), dimension(:,:),   pointer :: rrlv
  real(kind=realType), dimension(:,:),   pointer :: dd2Wall

contains
#ifndef USE_TAPENADE
  subroutine sa_block(resOnly)
    !
    !       sa solves the transport equation for the Spalart-Allmaras
    !       turbulence model in a decoupled manner using a diagonal
    !       dominant ADI-scheme. Note that the scratch and boundary
    !       matrix values are not strictly, but tapande would like to
    !       see them becuase it must save them.
    !
    use constants
    use blockPointers, only : nDom, il, jl, kl, scratch, bmtj1, bmtj2, &
         bmti1, bmti2, bmtk1, bmtk2
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use iteration, only : currentLevel
    use inputPhysics, only : turbProd
    use paramTurb
    use turbutils
    use turbBCRoutines
    implicit none
    !
    !      Subroutine argument.
    !
    logical, intent(in) :: resOnly
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, sps


    ! Set the arrays for the boundary condition treatment.
    call bcTurbTreatment

    ! Alloc central jacobian memory
    allocate(qq(2:il,2:jl,2:kl))

    ! Source Terms
    call saSource

    ! Advection Term
    nn = itu1 - 1
    call turbAdvection(1_intType, 1_intType, nn, qq)

    ! Unsteady Term
    call unsteadyTurbTerm(1_intType, 1_intType, nn, qq)

    ! Viscous Terms
    call saViscous

    ! Perform the residual scaling
    call saResScale

    ! We need to do an acutal solve. Solve and update the eddy
    ! viscosity and the boundary conditions

    if(.not. resOnly ) then

       ! Do solve
       call saSolve

       ! Compute the corresponding eddy viscosity.

       call saEddyViscosity(2, il, 2, jl, 2, kl)

       ! Set the halo values for the turbulent variables.
       ! We are on the finest mesh, so the second layer of halo
       ! cells must be computed as well.

       call applyAllTurbBCThisBlock(.true.)
    endif

    deallocate(qq)
  end subroutine sa_block
#endif

  subroutine saSource
    !
    !  Source terms.
    !  Determine the source term and its derivative w.r.t. nuTilde
    !  for all internal cells of the block.
    !  Remember that the SA field variable nuTilde = w(i,j,k,itu1)

    use blockPointers
    use constants
    use paramTurb
    use section
    use inputPhysics
    use inputDiscretization, only : approxSA
    use flowVarRefState
    implicit none

    ! Local parameters
    real(kind=realType), parameter :: f23 = two*third

    ! Local variables.
    integer(kind=intType) :: i, j, k, nn, ii
    real(kind=realType) :: fv1, fv2, ft2
    real(kind=realType) :: ss, sst, nu, dist2Inv, chi, chi2, chi3
    real(kind=realType) :: rr, gg, gg6, termFw, fwSa, term1, term2
    real(kind=realType) :: dfv1, dfv2, dft2, drr, dgg, dfw
    real(kind=realType) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
    real(kind=realType) :: div2, fact, sxx, syy, szz, sxy, sxz, syz
    real(kind=realType) :: vortx, vorty, vortz
    real(kind=realType) :: omegax, omegay, omegaz
    real(kind=realType) :: strainMag2, strainProd, vortProd
    real(kind=realType), parameter :: xminn = 1.e-10_realType

    ! Set model constants
    cv13    = rsaCv1**3
    kar2Inv = one/(rsaK**2)
    cw36    = rsaCw3**6
    cb3Inv  = one/rsaCb3

    ! Determine the non-dimensional wheel speed of this block.

    omegax = timeRef*sections(sectionID)%rotRate(1)
    omegay = timeRef*sections(sectionID)%rotRate(2)
    omegaz = timeRef*sections(sectionID)%rotRate(3)

    ! Create switches to production term depending on the variable that
    ! should be used
    if (turbProd .eq. katoLaunder) then
       print *,'katoLaunder production term not supported for SA'
       stop
    end if

#ifdef TAPENADE_REVERSE
    !$AD II-LOOP
    do ii=0,nx*ny*nz-1
       i = mod(ii, nx) + 2
       j = mod(ii/nx, ny) + 2
       k = ii/(nx*ny) + 2
#else
       do k=2, kl
          do j=2, jl
             do i=2, il
#endif
                ! Compute the gradient of u in the cell center. Use is made
                ! of the fact that the surrounding normals sum up to zero,
                ! such that the cell i,j,k does not give a contribution.
                ! The gradient is scaled by the factor 2*vol.

                uux = w(i+1,j,k,ivx)*si(i,j,k,1) - w(i-1,j,k,ivx)*si(i-1,j,k,1) &
                     + w(i,j+1,k,ivx)*sj(i,j,k,1) - w(i,j-1,k,ivx)*sj(i,j-1,k,1) &
                     + w(i,j,k+1,ivx)*sk(i,j,k,1) - w(i,j,k-1,ivx)*sk(i,j,k-1,1)
                uuy = w(i+1,j,k,ivx)*si(i,j,k,2) - w(i-1,j,k,ivx)*si(i-1,j,k,2) &
                     + w(i,j+1,k,ivx)*sj(i,j,k,2) - w(i,j-1,k,ivx)*sj(i,j-1,k,2) &
                     + w(i,j,k+1,ivx)*sk(i,j,k,2) - w(i,j,k-1,ivx)*sk(i,j,k-1,2)
                uuz = w(i+1,j,k,ivx)*si(i,j,k,3) - w(i-1,j,k,ivx)*si(i-1,j,k,3) &
                     + w(i,j+1,k,ivx)*sj(i,j,k,3) - w(i,j-1,k,ivx)*sj(i,j-1,k,3) &
                     + w(i,j,k+1,ivx)*sk(i,j,k,3) - w(i,j,k-1,ivx)*sk(i,j,k-1,3)

                ! Idem for the gradient of v.

                vvx = w(i+1,j,k,ivy)*si(i,j,k,1) - w(i-1,j,k,ivy)*si(i-1,j,k,1) &
                     + w(i,j+1,k,ivy)*sj(i,j,k,1) - w(i,j-1,k,ivy)*sj(i,j-1,k,1) &
                     + w(i,j,k+1,ivy)*sk(i,j,k,1) - w(i,j,k-1,ivy)*sk(i,j,k-1,1)
                vvy = w(i+1,j,k,ivy)*si(i,j,k,2) - w(i-1,j,k,ivy)*si(i-1,j,k,2) &
                     + w(i,j+1,k,ivy)*sj(i,j,k,2) - w(i,j-1,k,ivy)*sj(i,j-1,k,2) &
                     + w(i,j,k+1,ivy)*sk(i,j,k,2) - w(i,j,k-1,ivy)*sk(i,j,k-1,2)
                vvz = w(i+1,j,k,ivy)*si(i,j,k,3) - w(i-1,j,k,ivy)*si(i-1,j,k,3) &
                     + w(i,j+1,k,ivy)*sj(i,j,k,3) - w(i,j-1,k,ivy)*sj(i,j-1,k,3) &
                     + w(i,j,k+1,ivy)*sk(i,j,k,3) - w(i,j,k-1,ivy)*sk(i,j,k-1,3)

                ! And for the gradient of w.

                wwx = w(i+1,j,k,ivz)*si(i,j,k,1) - w(i-1,j,k,ivz)*si(i-1,j,k,1) &
                     + w(i,j+1,k,ivz)*sj(i,j,k,1) - w(i,j-1,k,ivz)*sj(i,j-1,k,1) &
                     + w(i,j,k+1,ivz)*sk(i,j,k,1) - w(i,j,k-1,ivz)*sk(i,j,k-1,1)
                wwy = w(i+1,j,k,ivz)*si(i,j,k,2) - w(i-1,j,k,ivz)*si(i-1,j,k,2) &
                     + w(i,j+1,k,ivz)*sj(i,j,k,2) - w(i,j-1,k,ivz)*sj(i,j-1,k,2) &
                     + w(i,j,k+1,ivz)*sk(i,j,k,2) - w(i,j,k-1,ivz)*sk(i,j,k-1,2)
                wwz = w(i+1,j,k,ivz)*si(i,j,k,3) - w(i-1,j,k,ivz)*si(i-1,j,k,3) &
                     + w(i,j+1,k,ivz)*sj(i,j,k,3) - w(i,j-1,k,ivz)*sj(i,j-1,k,3) &
                     + w(i,j,k+1,ivz)*sk(i,j,k,3) - w(i,j,k-1,ivz)*sk(i,j,k-1,3)

                ! Compute the components of the stress tensor.
                ! The combination of the current scaling of the velocity
                ! gradients (2*vol) and the definition of the stress tensor,
                ! leads to the factor 1/(4*vol).

                fact = fourth/vol(i,j,k)

                if (turbProd .eq. strain) then

                   sxx = two*fact*uux
                   syy = two*fact*vvy
                   szz = two*fact*wwz

                   sxy = fact*(uuy + vvx)
                   sxz = fact*(uuz + wwx)
                   syz = fact*(vvz + wwy)

                   ! Compute 2/3 * divergence of velocity squared

                   div2 = f23*(sxx+syy+szz)**2

                   ! Compute strain production term

                   strainMag2 = two*(sxy**2 + sxz**2 + syz**2) &
                        +           sxx**2 + syy**2 + szz**2

                   strainProd = two*strainMag2 - div2

                   ss = sqrt(strainProd)

                else if (turbProd .eq. vorticity) then

                   ! Compute the three components of the vorticity vector.
                   ! Substract the part coming from the rotating frame.

                   vortx = two*fact*(wwy - vvz) - two*omegax
                   vorty = two*fact*(uuz - wwx) - two*omegay
                   vortz = two*fact*(vvx - uuy) - two*omegaz

                   ! Compute the vorticity production term

                   vortProd = vortx**2 + vorty**2 + vortz**2

                   ! First take the square root of the production term to
                   ! obtain the correct production term for spalart-allmaras.
                   ! We do this to avoid if statements.

                   ss = sqrt(vortProd)

                end if

                ! Compute the laminar kinematic viscosity, the inverse of
                ! wall distance squared, the ratio chi (ratio of nuTilde
                ! and nu) and the functions fv1 and fv2. The latter corrects
                ! the production term near a viscous wall.

                nu       = rlv(i,j,k)/w(i,j,k,irho)
                dist2Inv = one/(d2Wall(i,j,k)**2)
                chi      = w(i,j,k,itu1)/nu
                chi2     = chi*chi
                chi3     = chi*chi2
                fv1      = chi3/(chi3+cv13)
                fv2      = one - chi/(one + chi*fv1)

                ! The function ft2, which is designed to keep a laminar
                ! solution laminar. When running in fully turbulent mode
                ! this function should be set to 0.0.

                if (useft2SA) then
                   ft2 = rsaCt3*exp(-rsaCt4*chi2)
                else
                   ft2 = zero
                end if

                ! Correct the production term to account for the influence
                ! of the wall.

                sst = ss + w(i,j,k,itu1)*fv2*kar2Inv*dist2Inv

                ! Add rotation term (useRotationSA defined in inputParams.F90)

                if (useRotationSA) then
                   sst = sst + rsaCrot*min(zero,sqrt(two*strainMag2))
                end if

                ! Make sure that this term remains positive
                ! (the function fv2 is negative between chi = 1 and 18.4,
                ! which can cause sst to go negative, which is undesirable).

                sst = max(sst,xminn)

                ! Compute the function fw. The argument rr is cut off at 10
                ! to avoid numerical problems. This is ok, because the
                ! asymptotical value of fw is then already reached.

                rr     = w(i,j,k,itu1)*kar2Inv*dist2Inv/sst
                rr     = min(rr,10.0_realType)
                gg     = rr + rsaCw2*(rr**6 - rr)
                gg6    = gg**6
                termFw = ((one + cw36)/(gg6 + cw36))**sixth
                fwSa   = gg*termFw

                ! Compute the source term; some terms are saved for the
                ! linearization. The source term is stored in dvt.

                if (approxSA) then
                  term1 = zero
                else
                  term1 = rsaCb1*(one-ft2)*ss
                end if
                term2 = dist2Inv*(kar2Inv*rsaCb1*((one-ft2)*fv2 + ft2) &
                     -           rsaCw1*fwSa)

                scratch(i,j,k,idvt) = (term1 + term2*w(i,j,k,itu1))*w(i,j,k,itu1)

#ifndef USE_TAPENADE
                ! Compute some derivatives w.r.t. nuTilde. These will occur
                ! in the left hand side, i.e. the matrix for the implicit
                ! treatment.

                dfv1 = three*chi2*cv13/((chi3+cv13)**2)
                dfv2 = (chi2*dfv1 - one)/(nu*((one + chi*fv1)**2))
                dft2 = -two*rsaCt4*chi*ft2/nu

                drr = (one - rr*(fv2 + w(i,j,k,itu1)*dfv2)) &
                     * kar2Inv*dist2Inv/sst
                dgg = (one - rsaCw2 + six*rsaCw2*(rr**5))*drr
                dfw = (cw36/(gg6 + cw36))*termFw*dgg

                ! Compute the source term jacobian. Note that the part
                ! containing term1 is treated explicitly. The reason is that
                ! implicit treatment of this part leads to a decrease of the
                ! diagonal dominance of the jacobian and it thus decreases
                ! the stability. You may want to play around and try to
                ! take this term into account in the jacobian.
                ! Note that -dsource/dnu is stored.
                qq(i,j,k) = -two*term2*w(i,j,k,itu1)                      &
                     -  dist2Inv*w(i,j,k,itu1)*w(i,j,k,itu1)         &
                     * (rsaCb1*kar2Inv*(dfv2-ft2*dfv2-fv2*dft2+dft2) &
                     -  rsaCw1*dfw)

                ! A couple of terms in qq may lead to a negative
                ! contribution. Clip qq to zero, if the total is negative.

                qq(i,j,k) = max(qq(i,j,k), zero)
#endif
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif
  end subroutine saSource

  subroutine saViscous
    !
    !  Viscous term.
    !  Determine the viscous contribution to the residual
    !  for all internal cells of the block.

    use blockPointers
    use paramTurb
    implicit none
    ! Local variables.
    integer(kind=intType) :: i, j, k, nn, ii
    real(kind=realType) :: nu
    real(kind=realType) :: fv1, fv2, ft2
    real(kind=realType) :: voli, volmi, volpi, xm, ym, zm, xp, yp, zp
    real(kind=realType) :: xa, ya, za, ttm, ttp, cnud, cam, cap
    real(kind=realType) :: nutm, nutp, num, nup, cdm, cdp
    real(kind=realType) :: c1m, c1p, c10, b1, c1, d1, qs

    ! Set model constants
    cv13    = rsaCv1**3
    kar2Inv = one/(rsaK**2)
    cw36    = rsaCw3**6
    cb3Inv  = one/rsaCb3

    !
    !       Viscous terms in k-direction.
    !
#ifdef TAPENADE_REVERSE
    !$AD II-LOOP
    do ii=0,nx*ny*nz-1
       i = mod(ii, nx) + 2
       j = mod(ii/nx, ny) + 2
       k = ii/(nx*ny) + 2
#else
       do k=2, kl
          do j=2, jl
             do i=2, il
#endif
                ! Compute the metrics in zeta-direction, i.e. along the
                ! line k = constant.

                voli  = one/vol(i,j,k)
                volmi = two/(vol(i,j,k) + vol(i,j,k-1))
                volpi = two/(vol(i,j,k) + vol(i,j,k+1))

                xm = sk(i,j,k-1,1)*volmi
                ym = sk(i,j,k-1,2)*volmi
                zm = sk(i,j,k-1,3)*volmi
                xp = sk(i,j,k,  1)*volpi
                yp = sk(i,j,k,  2)*volpi
                zp = sk(i,j,k,  3)*volpi

                xa  = half*(sk(i,j,k,1) + sk(i,j,k-1,1))*voli
                ya  = half*(sk(i,j,k,2) + sk(i,j,k-1,2))*voli
                za  = half*(sk(i,j,k,3) + sk(i,j,k-1,3))*voli
                ttm = xm*xa + ym*ya + zm*za
                ttp = xp*xa + yp*ya + zp*za

                ! ttm and ttp ~ 1/deltaX^2

                ! Computation of the viscous terms in zeta-direction; note
                ! that cross-derivatives are neglected, i.e. the mesh is
                ! assumed to be orthogonal.
                ! Furthermore, the grad(nu)**2 has been rewritten as
                ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
                ! The second derivative in zeta-direction is constructed as
                ! the central difference of the first order derivatives, i.e.
                ! d^2/dzeta^2 = d/dzeta (d/dzeta k+1/2 - d/dzeta k-1/2).
                ! In this way the metric can be taken into account.

                ! Compute the diffusion coefficients multiplying the nodes
                ! k+1, k and k-1 in the second derivative. Make sure that
                ! these coefficients are nonnegative.

                cnud = -rsaCb2*w(i,j,k,itu1)*cb3Inv
                cam  =  ttm*cnud
                cap  =  ttp*cnud

                ! Compute nuTilde at the faces

                nutm = half*(w(i,j,k-1,itu1) + w(i,j,k,itu1))
                nutp = half*(w(i,j,k+1,itu1) + w(i,j,k,itu1))

                ! Compute nu at the faces

                nu   = rlv(i,j,k)/w(i,j,k,irho)
                num  = half*(rlv(i,j,k-1)/w(i,j,k-1,irho) + nu)
                nup  = half*(rlv(i,j,k+1)/w(i,j,k+1,irho) + nu)

                cdm  = (num + (one + rsaCb2)*nutm)*ttm*cb3Inv
                cdp  = (nup + (one + rsaCb2)*nutp)*ttp*cb3Inv

                c1m = max(cdm+cam, zero)
                c1p = max(cdp+cap, zero)
                c10 = c1m + c1p

                ! Update the residual for this cell and store the possible
                ! coefficients for the matrix in b1, c1 and d1.

                scratch(i,j,k,idvt) = scratch(i,j,k,idvt)      + c1m*w(i,j,k-1,itu1) &
                     - c10*w(i,j,k,itu1) + c1p*w(i,j,k+1,itu1)
#ifndef USE_TAPENADE
                b1 = -c1m
                c1 =  c10
                d1 = -c1p

                ! Update the central jacobian. For nonboundary cells this
                ! is simply c1. For boundary cells this is slightly more
                ! complicated, because the boundary conditions are treated
                ! implicitly and the off-diagonal terms b1 and d1 must be
                ! taken into account.
                ! The boundary conditions are only treated implicitly if
                ! the diagonal dominance of the matrix is increased.

                if(k == 2) then
                   qq(i,j,k) = qq(i,j,k) + c1 &
                        - b1*max(bmtk1(i,j,itu1,itu1),zero)
                else if(k == kl) then
                   qq(i,j,k) = qq(i,j,k) + c1 &
                        - d1*max(bmtk2(i,j,itu1,itu1),zero)
                else
                   qq(i,j,k) = qq(i,j,k) + c1
                endif
#endif
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif
    !
    !       Viscous terms in j-direction.
    !
#ifdef TAPENADE_REVERSE
    !$AD II-LOOP
    do ii=0,nx*ny*nz-1
       i = mod(ii, nx) + 2
       j = mod(ii/nx, ny) + 2
       k = ii/(nx*ny) + 2
#else
       do k=2, kl
          do j=2, jl
             do i=2, il
#endif
                ! Compute the metrics in eta-direction, i.e. along the
                ! line j = constant.

                voli  = one/vol(i,j,k)
                volmi = two/(vol(i,j,k) + vol(i,j-1,k))
                volpi = two/(vol(i,j,k) + vol(i,j+1,k))

                xm = sj(i,j-1,k,1)*volmi
                ym = sj(i,j-1,k,2)*volmi
                zm = sj(i,j-1,k,3)*volmi
                xp = sj(i,j,  k,1)*volpi
                yp = sj(i,j,  k,2)*volpi
                zp = sj(i,j,  k,3)*volpi

                xa  = half*(sj(i,j,k,1) + sj(i,j-1,k,1))*voli
                ya  = half*(sj(i,j,k,2) + sj(i,j-1,k,2))*voli
                za  = half*(sj(i,j,k,3) + sj(i,j-1,k,3))*voli
                ttm = xm*xa + ym*ya + zm*za
                ttp = xp*xa + yp*ya + zp*za

                ! Computation of the viscous terms in eta-direction; note
                ! that cross-derivatives are neglected, i.e. the mesh is
                ! assumed to be orthogonal.
                ! Furthermore, the grad(nu)**2 has been rewritten as
                ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
                ! The second derivative in eta-direction is constructed as
                ! the central difference of the first order derivatives, i.e.
                ! d^2/deta^2 = d/deta (d/deta j+1/2 - d/deta j-1/2).
                ! In this way the metric can be taken into account.

                ! Compute the diffusion coefficients multiplying the nodes
                ! j+1, j and j-1 in the second derivative. Make sure that
                ! these coefficients are nonnegative.

                cnud = -rsaCb2*w(i,j,k,itu1)*cb3Inv
                cam  =  ttm*cnud
                cap  =  ttp*cnud

                nutm = half*(w(i,j-1,k,itu1) + w(i,j,k,itu1))
                nutp = half*(w(i,j+1,k,itu1) + w(i,j,k,itu1))
                nu   = rlv(i,j,k)/w(i,j,k,irho)
                num  = half*(rlv(i,j-1,k)/w(i,j-1,k,irho) + nu)
                nup  = half*(rlv(i,j+1,k)/w(i,j+1,k,irho) + nu)
                cdm  = (num + (one + rsaCb2)*nutm)*ttm*cb3Inv
                cdp  = (nup + (one + rsaCb2)*nutp)*ttp*cb3Inv

                c1m = max(cdm+cam, zero)
                c1p = max(cdp+cap, zero)
                c10 = c1m + c1p

                ! Update the residual for this cell and store the possible
                ! coefficients for the matrix in b1, c1 and d1.

                scratch(i,j,k,idvt) = scratch(i,j,k,idvt)      + c1m*w(i,j-1,k,itu1) &
                     - c10*w(i,j,k,itu1) + c1p*w(i,j+1,k,itu1)
#ifndef USE_TAPENADE
                b1 = -c1m
                c1 =  c10
                d1 = -c1p

                ! Update the central jacobian. For nonboundary cells this
                ! is simply c1. For boundary cells this is slightly more
                ! complicated, because the boundary conditions are treated
                ! implicitly and the off-diagonal terms b1 and d1 must be
                ! taken into account.
                ! The boundary conditions are only treated implicitly if
                ! the diagonal dominance of the matrix is increased.

                if(j == 2) then
                   qq(i,j,k) = qq(i,j,k) + c1 &
                        - b1*max(bmtj1(i,k,itu1,itu1),zero)
                else if(j == jl) then
                   qq(i,j,k) = qq(i,j,k) + c1 &
                        - d1*max(bmtj2(i,k,itu1,itu1),zero)
                else
                   qq(i,j,k) = qq(i,j,k) + c1
                endif
#endif
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif
    !
    !       Viscous terms in i-direction.
    !
#ifdef TAPENADE_REVERSE
    !$AD II-LOOP
    do ii=0,nx*ny*nz-1
       i = mod(ii, nx) + 2
       j = mod(ii/nx, ny) + 2
       k = ii/(nx*ny) + 2
#else
       do k=2, kl
          do j=2, jl
             do i=2, il
#endif
                ! Compute the metrics in xi-direction, i.e. along the
                ! line i = constant.

                voli  = one/vol(i,j,k)
                volmi = two/(vol(i,j,k) + vol(i-1,j,k))
                volpi = two/(vol(i,j,k) + vol(i+1,j,k))

                xm = si(i-1,j,k,1)*volmi
                ym = si(i-1,j,k,2)*volmi
                zm = si(i-1,j,k,3)*volmi
                xp = si(i,  j,k,1)*volpi
                yp = si(i,  j,k,2)*volpi
                zp = si(i,  j,k,3)*volpi

                xa  = half*(si(i,j,k,1) + si(i-1,j,k,1))*voli
                ya  = half*(si(i,j,k,2) + si(i-1,j,k,2))*voli
                za  = half*(si(i,j,k,3) + si(i-1,j,k,3))*voli
                ttm = xm*xa + ym*ya + zm*za
                ttp = xp*xa + yp*ya + zp*za

                ! Computation of the viscous terms in xi-direction; note
                ! that cross-derivatives are neglected, i.e. the mesh is
                ! assumed to be orthogonal.
                ! Furthermore, the grad(nu)**2 has been rewritten as
                ! div(nu grad(nu)) - nu div(grad nu) to enhance stability.
                ! The second derivative in xi-direction is constructed as
                ! the central difference of the first order derivatives, i.e.
                ! d^2/dxi^2 = d/dxi (d/dxi i+1/2 - d/dxi i-1/2).
                ! In this way the metric can be taken into account.

                ! Compute the diffusion coefficients multiplying the nodes
                ! i+1, i and i-1 in the second derivative. Make sure that
                ! these coefficients are nonnegative.

                cnud = -rsaCb2*w(i,j,k,itu1)*cb3Inv
                cam  =  ttm*cnud
                cap  =  ttp*cnud

                nutm = half*(w(i-1,j,k,itu1) + w(i,j,k,itu1))
                nutp = half*(w(i+1,j,k,itu1) + w(i,j,k,itu1))
                nu   = rlv(i,j,k)/w(i,j,k,irho)
                num  = half*(rlv(i-1,j,k)/w(i-1,j,k,irho) + nu)
                nup  = half*(rlv(i+1,j,k)/w(i+1,j,k,irho) + nu)
                cdm  = (num + (one + rsaCb2)*nutm)*ttm*cb3Inv
                cdp  = (nup + (one + rsaCb2)*nutp)*ttp*cb3Inv

                c1m = max(cdm+cam, zero)
                c1p = max(cdp+cap, zero)
                c10 = c1m + c1p

                ! Update the residual for this cell and store the possible
                ! coefficients for the matrix in b1, c1 and d1.

                scratch(i,j,k,idvt) = scratch(i,j,k,idvt)      + c1m*w(i-1,j,k,itu1) &
                     - c10*w(i,j,k,itu1) + c1p*w(i+1,j,k,itu1)
#ifndef USE_TAPENADE
                b1 = -c1m
                c1 =  c10
                d1 = -c1p

                ! Update the central jacobian. For nonboundary cells this
                ! is simply c1. For boundary cells this is slightly more
                ! complicated, because the boundary conditions are treated
                ! implicitly and the off-diagonal terms b1 and d1 must be
                ! taken into account.
                ! The boundary conditions are only treated implicitly if
                ! the diagonal dominance of the matrix is increased.

                if(i == 2) then
                   qq(i,j,k) = qq(i,j,k) + c1 &
                        - b1*max(bmti1(j,k,itu1,itu1),zero)
                else if(i == il) then
                   qq(i,j,k) = qq(i,j,k) + c1 &
                        - d1*max(bmti2(j,k,itu1,itu1),zero)
                else
                   qq(i,j,k) = qq(i,j,k) + c1
                endif
#endif
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif
  end subroutine saViscous

  subroutine saResScale

    !
    !  Multiply the residual by the volume and store this in dw; this
    ! * is done for monitoring reasons only. The multiplication with the
    ! * volume is present to be consistent with the flow residuals; also
    !  the negative value is taken, again to be consistent with the
    ! * flow equations. Also multiply by iblank so that no updates occur
    !  in holes or the overset boundary.
    use blockPointers
    implicit none

    ! Local variables
    integer(kind=intType) :: i,j,k,ii
    real(kind=realType) :: rblank

#ifdef TAPENADE_REVERSE
    !$AD II-LOOP
    do ii=0,nx*ny*nz-1
       i = mod(ii, nx) + 2
       j = mod(ii/nx, ny) + 2
       k = ii/(nx*ny) + 2
#else
       do k=2, kl
          do j=2, jl
             do i=2, il
#endif
                rblank = max(real(iblank(i,j,k), realType), zero)
                dw(i,j,k,itu1) = -volRef(i,j,k)*scratch(i,j,k,idvt)*rblank
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif
  end subroutine saResScale

#ifndef USE_TAPENADE
  subroutine saSolve
    !
    !  saSolve solves the turbulent transport equation for the
    !  original Spalart-Allmaras model in a decoupled manner using
    !  a diagonal dominant ADI-scheme.
    use blockPointers
    use inputIteration
    use inputPhysics
    use paramTurb
    use turbutils
    use turbCurveFits, only : curveTupYp
    implicit none

    integer(kind=intType) :: i, j, k, nn, ii
    real(kind=realType), dimension(2:max(kl,il,jl)) :: bb, cc, dd, ff
    real(kind=realType) :: voli, volmi, volpi, xm, ym, zm, xp, yp, zp
    real(kind=realType) :: xa, ya, za, ttm, ttp, cnud, cam, cap
    real(kind=realType) :: nutm, nutp, num, nup, cdm, cdp
    real(kind=realType) :: c1m, c1p, c10, b1, c1, d1, qs
    real(kind=realType) :: uu, um, up, factor, f, tu1p(1), nu, rblank

    logical, dimension(2:jl,2:kl), target :: flagI2, flagIl
    logical, dimension(2:il,2:kl), target :: flagJ2, flagJl
    logical, dimension(2:il,2:jl), target :: flagK2, flagKl
    logical, dimension(:,:), pointer :: flag

    ! Initialize the wall function flags to .false.

    flagI2 = .false.
    flagIl = .false.
    flagJ2 = .false.
    flagJl = .false.
    flagK2 = .false.
    flagKl = .false.

    ! Modify the rhs of the 1st internal cell, if wall functions
    ! are used; their value is determined by the table.

    testWallFunctions: if( wallFunctions ) then

       bocos: do nn=1,nViscBocos

          ! Determine the block face on which the subface is located
          ! and set some variables. As flag points to the entire array
          ! flagI2, etc., its starting indices are the starting indices
          ! of its target and not 1.

          select case (BCFaceID(nn))
          case (iMin)
             flag    => flagI2
             ddw     => dw(2,1:,1:,1:); ddvt => scratch(2,1:,1:,idvt:)
             ww      => w(2,1:,1:,1:);  rrlv => rlv(2,1:,1:)
             dd2Wall => d2Wall(2,:,:)

          case (iMax)
             flag    => flagIl
             ddw     => dw(il,1:,1:,1:); ddvt => scratch(il,1:,1:,idvt:)
             ww      => w(il,1:,1:,1:);  rrlv => rlv(il,1:,1:)
             dd2Wall => d2Wall(il,:,:)

          case (jMin)
             flag    => flagJ2
             ddw     => dw(1:,2,1:,1:); ddvt => scratch(1:,2,1:,idvt:)
             ww      => w(1:,2,1:,1:);  rrlv => rlv(1:,2,1:)
             dd2Wall => d2Wall(:,2,:)

          case (jMax)
             flag    => flagJl
             ddw     => dw(1:,jl,1:,1:); ddvt => scratch(1:,jl,1:,idvt:)
             ww      => w(1:,jl,1:,1:);  rrlv => rlv(1:,jl,1:)
             dd2Wall => d2Wall(:,jl,:)

          case (kMin)
             flag    => flagK2
             ddw     => dw(1:,1:,2,1:); ddvt => scratch(1:,1:,2,idvt:)
             ww      => w(1:,1:,2,1:);  rrlv => rlv(1:,1:,2)
             dd2Wall => d2Wall(:,:,2)

          case (kMax)
             flag    => flagKl
             ddw     => dw(1:,1:,kl,:); ddvt => scratch(1:,1:,kl,idvt:)
             ww      => w(1:,1:,kl,1:); rrlv => rlv(1:,1:,kl)
             dd2Wall => d2Wall(:,:,kl)

          end select

          ! Loop over the owned faces of this subface. Therefore the
          ! nodal range of BCData must be used. The offset of +1 is
          ! present, because the starting index of the cell range is
          ! 1 larger than the starting index of the nodal range.

          do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
             do i=(BCData(nn)%inBeg+1),BCData(nn)%inEnd

                ! Set ddw to zero.

                ddw(i,j,itu1) = zero

                ! Enforce nu tilde in the 1st internal cell from the
                ! wall function table. There is an offset of -1 in the
                ! wall distance. Note that the offset compared to the
                ! current value must be stored, because dvt contains
                ! the update. Also note that the curve fits contain the
                ! non-dimensional value.

                yp = ww(i,j,irho)*dd2Wall(i-1,j-1) &
                     * viscSubface(nn)%utau(i,j)/rrlv(i,j)

                call curveTupYp(tu1p, yp, itu1, itu1)
                ddvt(i,j,1) = tu1p(1)*rrlv(i,j)/ww(i,j,irho) - ww(i,j,itu1)

                ! Set the wall flag to .true.

                flag(i,j) = .true.

             enddo
          enddo

       enddo bocos
    endif testWallFunctions

    ! For implicit relaxation take the local time step into account,
    ! where dt is the inverse of the central jacobian times the cfl
    ! number. The following system is solved:
    ! (I/dt + cc + bb + dd)*dw = rhs, in which I/dt = cc/cfl. As in
    ! the rest of the algorithm only the modified central jacobian is
    ! used, stored it now.

    ! Compute the factor multiplying the central jacobian, which
    ! is 1 + 1/cfl (implicit relaxation only).

    factor = one
    if(turbRelax == turbRelaxImplicit) &
         factor = one + (one-alfaTurb)/alfaTurb

    do k=2,kl
       do j=2,jl
          do i=2,il

             qq(i,j,k) = factor*qq(i,j,k)

             ! Set qq to 1 if the value is determined by the
             ! wall function table.

             if((i == 2  .and. flagI2(j,k)) .or. &
                  (i == il .and. flagIl(j,k)) .or. &
                  (j == 2  .and. flagJ2(i,k)) .or. &
                  (j == jl .and. flagJl(i,k)) .or. &
                  (k == 2  .and. flagK2(i,j)) .or. &
                  (k == kl .and. flagKl(i,j))) qq(i,j,k) = one

          enddo
       enddo
    enddo

    ! Initialize the grid velocity to zero. This value will be used
    ! if the block is not moving.

    qs = zero
    !
    !       dd-ADI step in j-direction. There is no particular reason to
    !       start in j-direction, it just happened to be so. As we solve
    !       in j-direction, the j-loop is the innermost loop.
    !
    do k=2,kl
       do i=2,il
          do j=2,jl

             ! More or less the same code is executed here as above when
             ! the residual was built. However, now the off-diagonal
             ! terms for the dd-ADI must be built and stored. This could
             ! have been done earlier, but then all the coefficients had
             ! to be stored. To save memory, they are recomputed.
             ! Consequently, see the j-loop to build the residual for
             ! the comments.

             voli  = one/vol(i,j,k)
             volmi = two/(vol(i,j,k) + vol(i,j-1,k))
             volpi = two/(vol(i,j,k) + vol(i,j+1,k))

             xm = sj(i,j-1,k,1)*volmi
             ym = sj(i,j-1,k,2)*volmi
             zm = sj(i,j-1,k,3)*volmi
             xp = sj(i,j,  k,1)*volpi
             yp = sj(i,j,  k,2)*volpi
             zp = sj(i,j,  k,3)*volpi

             xa  = half*(sj(i,j,k,1) + sj(i,j-1,k,1))*voli
             ya  = half*(sj(i,j,k,2) + sj(i,j-1,k,2))*voli
             za  = half*(sj(i,j,k,3) + sj(i,j-1,k,3))*voli
             ttm = xm*xa + ym*ya + zm*za
             ttp = xp*xa + yp*ya + zp*za

             cnud = -rsaCb2*w(i,j,k,itu1)*cb3Inv
             cam  =  ttm*cnud
             cap  =  ttp*cnud

             ! Off-diagonal terms due to the diffusion terms
             ! in j-direction.

             nutm = half*(w(i,j-1,k,itu1) + w(i,j,k,itu1))
             nutp = half*(w(i,j+1,k,itu1) + w(i,j,k,itu1))
             nu   = rlv(i,j,k)/w(i,j,k,irho)
             num  = half*(rlv(i,j-1,k)/w(i,j-1,k,irho) + nu)
             nup  = half*(rlv(i,j+1,k)/w(i,j+1,k,irho) + nu)
             cdm  = (num + (one + rsaCb2)*nutm)*ttm*cb3Inv
             cdp  = (nup + (one + rsaCb2)*nutp)*ttp*cb3Inv

             c1m = max(cdm+cam, zero)
             c1p = max(cdp+cap, zero)

             bb(j) = -c1m
             dd(j) = -c1p

             ! Compute the grid velocity if present.
             ! It is taken as the average of j and j-1,

             if( addGridVelocities ) &
                  qs = half*(sFaceJ(i,j,k) + sFaceJ(i,j-1,k))*voli

             ! Off-diagonal terms due to the advection term in
             ! j-direction. First order approximation.

             uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs
             um = zero
             up = zero
             if(uu < zero) um = uu
             if(uu > zero) up = uu

             bb(j) = bb(j) - up
             dd(j) = dd(j) + um

             ! Store the central jacobian and rhs in cc and ff.
             ! Multiply the off-diagonal terms and rhs by the iblank
             ! value so the update determined for iblank = 0 is zero.

             rblank = max(real(iblank(i,j,k), realType), zero)

             cc(j) = qq(i,j,k)
             ff(j) = scratch(i,j,k,idvt)*rblank

             bb(j) = bb(j)*rblank
             dd(j) = dd(j)*rblank

             ! Set the off diagonal terms to zero if the wall is flagged.

             if((i == 2  .and. flagI2(j,k)) .or. &
                  (i == il .and. flagIl(j,k)) .or. &
                  (j == 2  .and. flagJ2(i,k)) .or. &
                  (j == jl .and. flagJl(i,k)) .or. &
                  (k == 2  .and. flagK2(i,j)) .or. &
                  (k == kl .and. flagKl(i,j))) then
                bb(j) = zero
                dd(j) = zero
             endif

          enddo

          ! Solve the tri-diagonal system in j-direction.
          ! First the backward sweep to eliMinate the upper diagonal dd.

          do j=ny,2,-1
             f     = dd(j)/cc(j+1)
             cc(j) = cc(j) - f*bb(j+1)
             ff(j) = ff(j) - f*ff(j+1)
          enddo

          ! The matrix is now in lower block bi-diagonal form.
          ! Perform a forward sweep to compute the solution.

          ff(2) = ff(2)/cc(2)
          do j=3,jl
             ff(j) = ff(j) - bb(j)*ff(j-1)
             ff(j) = ff(j)/cc(j)
          enddo

          ! Determine the new rhs for the next direction.

          do j=2,jl
             scratch(i,j,k,idvt) = ff(j)*qq(i,j,k)
          enddo

       enddo
    enddo
    !
    !       dd-ADI step in i-direction. As we solve in i-direction, the
    !       i-loop is the innermost loop.
    !
    do k=2,kl
       do j=2,jl
          do i=2,il

             ! More or less the same code is executed here as above when
             ! the residual was built. However, now the off-diagonal
             ! terms for the dd-ADI must be built and stored. This could
             ! have been done earlier, but then all the coefficients had
             ! to be stored. To save memory, they are recomputed.
             ! Consequently, see the i-loop to build the residual for
             ! the comments.

             voli  = one/vol(i,j,k)
             volmi = two/(vol(i,j,k) + vol(i-1,j,k))
             volpi = two/(vol(i,j,k) + vol(i+1,j,k))

             xm = si(i-1,j,k,1)*volmi
             ym = si(i-1,j,k,2)*volmi
             zm = si(i-1,j,k,3)*volmi
             xp = si(i,  j,k,1)*volpi
             yp = si(i,  j,k,2)*volpi
             zp = si(i,  j,k,3)*volpi

             xa  = half*(si(i,j,k,1) + si(i-1,j,k,1))*voli
             ya  = half*(si(i,j,k,2) + si(i-1,j,k,2))*voli
             za  = half*(si(i,j,k,3) + si(i-1,j,k,3))*voli
             ttm = xm*xa + ym*ya + zm*za
             ttp = xp*xa + yp*ya + zp*za

             cnud = -rsaCb2*w(i,j,k,itu1)*cb3Inv
             cam  =  ttm*cnud
             cap  =  ttp*cnud

             ! Off-diagonal terms due to the diffusion terms
             ! in i-direction.

             nutm = half*(w(i-1,j,k,itu1) + w(i,j,k,itu1))
             nutp = half*(w(i+1,j,k,itu1) + w(i,j,k,itu1))
             nu   = rlv(i,j,k)/w(i,j,k,irho)
             num  = half*(rlv(i-1,j,k)/w(i-1,j,k,irho) + nu)
             nup  = half*(rlv(i+1,j,k)/w(i+1,j,k,irho) + nu)
             cdm  = (num + (one + rsaCb2)*nutm)*ttm*cb3Inv
             cdp  = (nup + (one + rsaCb2)*nutp)*ttp*cb3Inv

             c1m = max(cdm+cam, zero)
             c1p = max(cdp+cap, zero)

             bb(i) = -c1m
             dd(i) = -c1p

             ! Compute the grid velocity if present.
             ! It is taken as the average of i and i-1,

             if( addGridVelocities ) &
                  qs = half*(sFaceI(i,j,k) + sFaceI(i-1,j,k))*voli

             ! Off-diagonal terms due to the advection term in
             ! i-direction. First order approximation.

             uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs
             um = zero
             up = zero
             if(uu < zero) um = uu
             if(uu > zero) up = uu

             bb(i) = bb(i) - up
             dd(i) = dd(i) + um

             ! Store the central jacobian and rhs in cc and ff.
             ! Multiply the off-diagonal terms and rhs by the iblank
             ! value so the update determined for iblank = 0 is zero.

             rblank = max(real(iblank(i,j,k), realType), zero)

             cc(i) = qq(i,j,k)
             ff(i) = scratch(i,j,k,idvt)*rblank

             bb(i) = bb(i)*rblank
             dd(i) = dd(i)*rblank

             ! Set the off diagonal terms to zero if the wall is flagged.

             if((i == 2  .and. flagI2(j,k)) .or. &
                  (i == il .and. flagIl(j,k)) .or. &
                  (j == 2  .and. flagJ2(i,k)) .or. &
                  (j == jl .and. flagJl(i,k)) .or. &
                  (k == 2  .and. flagK2(i,j)) .or. &
                  (k == kl .and. flagKl(i,j))) then
                bb(i) = zero
                dd(i) = zero
             endif

          enddo

          ! Solve the tri-diagonal system in i-direction.
          ! First the backward sweep to eliMinate the upper diagonal dd.

          do i=nx,2,-1
             f     = dd(i)/cc(i+1)
             cc(i) = cc(i) - f*bb(i+1)
             ff(i) = ff(i) - f*ff(i+1)
          enddo

          ! The matrix is now in lower block bi-diagonal form.
          ! Perform a forward sweep to compute the solution.

          ff(2) = ff(2)/cc(2)
          do i=3,il
             ff(i) = ff(i) - bb(i)*ff(i-1)
             ff(i) = ff(i)/cc(i)
          enddo

          ! Determine the new rhs for the next direction.

          do i=2,il
             scratch(i,j,k,idvt) = ff(i)*qq(i,j,k)
          enddo

       enddo
    enddo
    !
    !       dd-ADI step in k-direction. As we solve in k-direction, the
    !       k-loop is the innermost loop.
    !
    do j=2,jl
       do i=2,il
          do k=2,kl

             ! More or less the same code is executed here as above when
             ! the residual was built. However, now the off-diagonal
             ! terms for the dd-ADI must be built and stored. This could
             ! have been done earlier, but then all the coefficients had
             ! to be stored. To save memory, they are recomputed.
             ! Consequently, see the k-loop to build the residual for
             ! the comments.

             voli  = one/vol(i,j,k)
             volmi = two/(vol(i,j,k) + vol(i,j,k-1))
             volpi = two/(vol(i,j,k) + vol(i,j,k+1))

             xm = sk(i,j,k-1,1)*volmi
             ym = sk(i,j,k-1,2)*volmi
             zm = sk(i,j,k-1,3)*volmi
             xp = sk(i,j,k,  1)*volpi
             yp = sk(i,j,k,  2)*volpi
             zp = sk(i,j,k,  3)*volpi

             xa  = half*(sk(i,j,k,1) + sk(i,j,k-1,1))*voli
             ya  = half*(sk(i,j,k,2) + sk(i,j,k-1,2))*voli
             za  = half*(sk(i,j,k,3) + sk(i,j,k-1,3))*voli
             ttm = xm*xa + ym*ya + zm*za
             ttp = xp*xa + yp*ya + zp*za

             cnud = -rsaCb2*w(i,j,k,itu1)*cb3Inv
             cam  =  ttm*cnud
             cap  =  ttp*cnud

             ! Off-diagonal terms due to the diffusion terms
             ! in k-direction.

             nutm = half*(w(i,j,k-1,itu1) + w(i,j,k,itu1))
             nutp = half*(w(i,j,k+1,itu1) + w(i,j,k,itu1))
             nu   = rlv(i,j,k)/w(i,j,k,irho)
             num  = half*(rlv(i,j,k-1)/w(i,j,k-1,irho) + nu)
             nup  = half*(rlv(i,j,k+1)/w(i,j,k+1,irho) + nu)
             cdm  = (num + (one + rsaCb2)*nutm)*ttm*cb3Inv
             cdp  = (nup + (one + rsaCb2)*nutp)*ttp*cb3Inv

             c1m = max(cdm+cam, zero)
             c1p = max(cdp+cap, zero)

             bb(k) = -c1m
             dd(k) = -c1p

             ! Compute the grid velocity if present.
             ! It is taken as the average of k and k-1,

             if( addGridVelocities ) &
                  qs = half*(sFaceK(i,j,k) + sFaceK(i,j,k-1))*voli

             ! Off-diagonal terms due to the advection term in
             ! k-direction. First order approximation.

             uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs
             um = zero
             up = zero
             if(uu < zero) um = uu
             if(uu > zero) up = uu

             bb(k) = bb(k) - up
             dd(k) = dd(k) + um

             ! Store the central jacobian and rhs in cc and ff.
             ! Multiply the off-diagonal terms and rhs by the iblank
             ! value so the update determined for iblank = 0 is zero.

             rblank = max(real(iblank(i,j,k), realType), zero)

             cc(k) = qq(i,j,k)
             ff(k) = scratch(i,j,k,idvt)*rblank

             bb(k) = bb(k)*rblank
             dd(k) = dd(k)*rblank

             ! Set the off diagonal terms to zero if the wall is flagged.

             if((i == 2  .and. flagI2(j,k)) .or. &
                  (i == il .and. flagIl(j,k)) .or. &
                  (j == 2  .and. flagJ2(i,k)) .or. &
                  (j == jl .and. flagJl(i,k)) .or. &
                  (k == 2  .and. flagK2(i,j)) .or. &
                  (k == kl .and. flagKl(i,j))) then
                bb(k) = zero
                dd(k) = zero
             endif

          enddo

          ! Solve the tri-diagonal system in k-direction.
          ! First the backward sweep to eliMinate the upper diagonal dd.

          do k=nz,2,-1
             f     = dd(k)/cc(k+1)
             cc(k) = cc(k) - f*bb(k+1)
             ff(k) = ff(k) - f*ff(k+1)
          enddo

          ! The matrix is now in lower block bi-diagonal form.
          ! Perform a forward sweep to compute the solution.

          ff(2) = ff(2)/cc(2)
          do k=3,kl
             ff(k) = ff(k) - bb(k)*ff(k-1)
             ff(k) = ff(k)/cc(k)
          enddo

          ! Store the update in dvt.

          do k=2,kl
             scratch(i,j,k,idvt) = ff(k)
          enddo

       enddo
    enddo
    !
    !       Update the turbulent variables. For explicit relaxation the
    !       update must be relaxed; for implicit relaxation this has been
    !       done via the time step.
    !
    factor = one
    if(turbRelax == turbRelaxExplicit) factor = alfaTurb

    do k=2,kl
       do j=2,jl
          do i=2,il
             w(i,j,k,itu1) = w(i,j,k,itu1) + factor*scratch(i,j,k,idvt)
             w(i,j,k,itu1) = max(w(i,j,k,itu1), zero)
          enddo
       enddo
    enddo

  end subroutine saSolve
#endif
end module sa
