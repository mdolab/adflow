module turbUtils

contains

  subroutine prodKatoLaunder
    !
    !       prodKatoLaunder computes the turbulent production term using
    !       the Kato-Launder formulation.
    !
    use constants
    use blockPointers, only : nx, ny, nz, il, jl, kl, w, si, sj, sk, vol, sectionID, scratch
    use flowVarRefState, only : timeRef
    use section, only : sections
    use turbMod, only : prod
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ii

    real(kind=realType) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
    real(kind=realType) :: qxx, qyy, qzz, qxy, qxz, qyz, sijsij
    real(kind=realType) :: oxy, oxz, oyz, oijoij
    real(kind=realType) :: fact, omegax, omegay, omegaz

    ! Determine the non-dimensional wheel speed of this block.
    ! The vorticity term, which appears in Kato-Launder is of course
    ! not frame invariant. To approximate frame invariance the wheel
    ! speed should be substracted from oxy, oxz and oyz, which results
    ! in the vorticity in the rotating frame. However some people
    ! claim that the absolute vorticity should be used to obtain the
    ! best results. In that omega should be set to zero.

    omegax = timeRef*sections(sectionID)%rotRate(1)
    omegay = timeRef*sections(sectionID)%rotRate(2)
    omegaz = timeRef*sections(sectionID)%rotRate(3)

    ! Loop over the cell centers of the given block. It may be more
    ! efficient to loop over the faces and to scatter the gradient,
    ! but in that case the gradients for u, v and w must be stored.
    ! In the current approach no extra memory is needed.
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
                ! The gradient is scaled by a factor 2*vol.

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

                ! Compute the strain and vorticity terms. The multiplication
                ! is present to obtain the correct gradients. Note that
                ! the wheel speed is substracted from the vorticity terms.

                fact = half/vol(i,j,k)

                qxx = fact*uux
                qyy = fact*vvy
                qzz = fact*wwz

                qxy = fact*half*(uuy + vvx)
                qxz = fact*half*(uuz + wwx)
                qyz = fact*half*(vvz + wwy)

                oxy = fact*half*(vvx - uuy) - omegaz
                oxz = fact*half*(uuz - wwx) - omegay
                oyz = fact*half*(wwy - vvz) - omegax

                ! Compute the summation of the strain and vorticity tensors.

                sijsij = two*(qxy**2 + qxz**2 + qyz**2) &
                     +      qxx**2 + qyy**2 + qzz**2
                oijoij = two*(oxy**2 + oxz**2 + oyz**2)

                ! Compute the production term.

                scratch(i,j,k, iprod) = two*sqrt(sijsij*oijoij)
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif
  end subroutine prodKatoLaunder

  subroutine prodSmag2
    !
    !       prodSmag2 computes the term:
    !              2*sij*sij - 2/3 div(u)**2 with  sij=0.5*(duidxj+dujdxi)
    !       which is used for the turbulence equations.
    !       It is assumed that the pointer prod, stored in turbMod, is
    !       already set to the correct entry.
    !
    use constants
    use blockPointers, only : nx, ny, nz, il, jl, kl, w, si, sj, sk, vol, sectionID, scratch
    implicit none
    !
    !      Local parameter
    !
    real(kind=realType), parameter :: f23 = two*third
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ii
    real(kind=realType)   :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
    real(kind=realType)   :: div2, fact, sxx, syy, szz, sxy, sxz, syz

    ! Loop over the cell centers of the given block. It may be more
    ! efficient to loop over the faces and to scatter the gradient,
    ! but in that case the gradients for u, v and w must be stored.
    ! In the current approach no extra memory is needed.

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

                sxx = two*fact*uux
                syy = two*fact*vvy
                szz = two*fact*wwz

                sxy = fact*(uuy + vvx)
                sxz = fact*(uuz + wwx)
                syz = fact*(vvz + wwy)

                ! Compute 2/3 * divergence of velocity squared

                div2 = f23*(sxx+syy+szz)**2

                ! Store the square of strain as the production term.

                scratch(i,j,k, iprod) = two*(two*(sxy**2 + sxz**2 + syz**2) &
                     +           sxx**2 + syy**2 + szz**2) - div2
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif
  end subroutine prodSmag2

  subroutine prodWmag2
    !
    !       prodWmag2 computes the term:
    !          2*oij*oij  with oij=0.5*(duidxj - dujdxi).
    !       This is equal to the magnitude squared of the vorticity.
    !       It is assumed that the pointer vort, stored in turbMod, is
    !       already set to the correct entry.
    !
    use constants
    use blockPointers, only : nx, ny, nz, il, jl, kl, w, si, sj, sk, vol, sectionID, scratch
    use flowVarRefState, only : timeRef
    use section, only : sections
    implicit none
    !
    !      Local variables.
    !
    integer :: i, j, k, ii

    real(kind=realType) :: uuy, uuz, vvx, vvz, wwx, wwy
    real(kind=realType) :: fact, vortx, vorty, vortz
    real(kind=realType) :: omegax, omegay, omegaz

    ! Determine the non-dimensional wheel speed of this block.

    omegax = timeRef*sections(sectionID)%rotRate(1)
    omegay = timeRef*sections(sectionID)%rotRate(2)
    omegaz = timeRef*sections(sectionID)%rotRate(3)

    ! Loop over the cell centers of the given block. It may be more
    ! efficient to loop over the faces and to scatter the gradient,
    ! but in that case the gradients for u, v and w must be stored.
    ! In the current approach no extra memory is needed.
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

                ! Compute the necessary derivatives of u in the cell center.
                ! Use is made of the fact that the surrounding normals sum up
                ! to zero, such that the cell i,j,k does not give a
                ! contribution. The gradient is scaled by a factor 2*vol.

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

                ! Compute the three components of the vorticity vector.
                ! Substract the part coming from the rotating frame.

                fact = half/vol(i,j,k)

                vortx = fact*(wwy - vvz) - two*omegax
                vorty = fact*(uuz - wwx) - two*omegay
                vortz = fact*(vvx - uuy) - two*omegaz

                ! Compute the magnitude squared of the vorticity.

                scratch(i,j,k,ivort) = vortx**2 + vorty**2 + vortz**2
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif
  end subroutine prodWmag2
  function saNuKnownEddyRatio(eddyRatio, nuLam)
    !
    !       saNuKnownEddyRatio computes the Spalart-Allmaras transport
    !       variable nu for the given eddy viscosity ratio.
    !
    use constants
    use paramTurb
    implicit none
    !
    !      Function type.
    !
    real(kind=realType) :: saNuKnownEddyRatio
    !
    !      Function arguments.
    !
    real(kind=realType), intent(in) :: eddyRatio, nuLam
    !
    !      Local variables.
    !
    real(kind=realType) :: cv13, chi, chi2, chi3, chi4, f, df, dchi

    ! Take care of the exceptional cases.

    if(eddyRatio <= zero) then
       saNuKnownEddyRatio = zero
       return
    endif

    ! Set the value of cv1^3, which is the constant appearing in the
    ! sa function fv1 to compute the eddy viscosity

    cv13 = rsaCv1**3

    ! Determine the value of chi, which is given by the quartic
    ! polynomial chi^4 - ratio*(chi^3 + cv1^3) = 0.
    ! First determine the start value, depending on the eddyRatio.

    if(eddyRatio < 1.e-4_realType) then
       chi = 0.5_realType
    else if(eddyRatio < 1.0_realType) then
       chi = 5.0_realType
    else if(eddyRatio < 10.0_realType) then
       chi = 10.0_realType
    else
       chi = eddyRatio
    endif

    ! The actual newton algorithm.

    do
       ! Compute the function value and the derivative.

       chi2 = chi*chi
       chi3 = chi*chi2
       chi4 = chi*chi3

       f  = chi4 - eddyRatio*(chi3 + cv13)
       df = four*chi3 - three*eddyRatio*chi2

       ! Compute the negative update and the new value of chi.

       dchi = f/df
       chi  = chi - dchi

       ! Condition to exit the loop.

       if(abs(dchi/chi) <= thresholdReal) exit
    enddo

    ! Chi is the ratio of the spalart allmaras transport variable and
    ! the laminar viscosity. So multiply chi with the laminar viscosity
    ! to obtain the correct value.

    saNuKnownEddyRatio = nuLam*chi

  end function saNuKnownEddyRatio

  subroutine unsteadyTurbTerm(mAdv, nAdv, offset, qq)
    !
    !       unsteadyTurbTerm discretizes the time derivative of the
    !       turbulence transport equations and add it to the residual.
    !       As the time derivative is the same for all turbulence models,
    !       this generic routine can be used; both the discretization of
    !       the time derivative and its contribution to the central
    !       jacobian are computed by this routine.
    !       Only nAdv equations are treated, while the actual system has
    !       size mAdv. The reason is that some equations for some
    !       turbulence equations do not have a time derivative, e.g. the
    !       f equation in the v2-f model. The argument offset indicates
    !       the offset in the w vector where this subsystem starts. As a
    !       consequence it is assumed that the indices of the current
    !       subsystem are contiguous, e.g. if a 2*2 system is solved the
    !       Last index in w is offset+1 and offset+2 respectively.
    !
    use blockPointers
    use flowVarRefState
    use inputPhysics
    use inputTimeSpectral
    use inputUnsteady
    use iteration
    use section
    use turbMod

    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: mAdv, nAdv, offset

    real(kind=realType), dimension(2:il,2:jl,2:kl,mAdv,mAdv), &
         intent(inout) :: qq
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ii, jj, nn

    real(kind=realType) :: oneOverDt, tmp

    ! Determine the equation mode.

    select case (equationMode)

    case (steady)

       ! Steady computation. No time derivative present.

       return

       !===============================================================

    case (unsteady)

       ! The time deritvative term depends on the integration
       ! scheme used.

       select case (timeIntegrationScheme)

       case (BDF)

          ! Backward difference formula is used as time
          ! integration scheme.

          ! Store the inverse of the physical nonDimensional
          ! time step a bit easier.

          oneOverDt = timeRef/deltaT

          ! Loop over the number of turbulent transport equations.

          nAdvLoopUnsteady: do ii=1,nAdv

             ! Store the index of the current turbulent variable in jj.

             jj = ii + offset

             ! Loop over the owned cells of this block to compute the
             ! time derivative.

             do k=2,kl
                do j=2,jl
                   do i=2,il

                      ! Initialize tmp to the value of the current
                      ! level multiplied by the corresponding coefficient
                      ! in the time integration scheme.

                      tmp = coefTime(0)*w(i,j,k,jj)

                      ! Loop over the old time levels and add the
                      ! corresponding contribution to tmp.

                      do nn=1,noldLevels
                         tmp = tmp + coefTime(nn)*wold(nn,i,j,k,jj)
                      enddo

                      ! Update the residual. Note that in the turbulent
                      ! routines the residual is defined with an opposite
                      ! sign compared to the residual of the flow equations.
                      ! Therefore the time derivative must be substracted
                      ! from dvt.

                      scratch(i,j,k,idvt+ii-1) = scratch(i,j,k,idvt+ii-1) - oneOverDt*tmp

                      ! Update the central jacobian.

                      qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) &
                           + coefTime(0)*oneOverDt
                   enddo
                enddo
             enddo

          enddo nAdvLoopUnsteady

          !===========================================================

       case (explicitRK)

          ! Explicit time integration scheme. The time derivative
          ! is handled differently.

          return

       end select

       !===============================================================

    case (timeSpectral)

       ! Time spectral method.

       ! Loop over the number of turbulent transport equations.

       nAdvLoopSpectral: do ii=1,nAdv

          ! Store the index of the current turbulent variable in jj.

          jj = ii + offset

          ! The time derivative has been computed earlier in
          ! unsteadyTurbSpectral and stored in entry jj of scratch.
          ! Substract this value for all owned cells. It must be
          ! substracted, because in the turbulent routines the
          ! residual is defined with an opposite sign compared to
          ! the residual of the flow equations.
          ! Also add a term to the diagonal matrix, which corresponds
          ! to to the contribution of the highest frequency. This is
          ! equivalent to an explicit treatment of the time derivative
          ! and may need to be changed.

          tmp = nTimeIntervalsSpectral*pi*timeRef &
               / sections(sectionID)%timePeriod

          do k=2,kl
             do j=2,jl
                do i=2,il
                   scratch(i,j,k,idvt+ii-1)   = scratch(i,j,k,idvt+ii-1)   - dw(i,j,k,jj)
                   qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) + tmp
                enddo
             enddo
          enddo

       enddo nAdvLoopSpectral

    end select

  end subroutine unsteadyTurbTerm

  subroutine computeEddyViscosity(includeHalos)
    !
    !       computeEddyViscosity computes the eddy viscosity in the
    !       owned cell centers of the given block. It is assumed that the
    !       pointes already point to the correct block before entering
    !       this subroutine.
    !
    use constants
    use flowVarRefState
    use inputPhysics
    use iteration
    use blockPointers
    implicit none

    ! Input Parameter
    logical, intent(in) :: includeHalos

    !
    !      Local variables.
    !
    logical :: returnImmediately
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    ! Check if an immediate return can be made.

    if( eddyModel ) then
       if((currentLevel <= groundLevel)) then
          returnImmediately = .false.
       else
          returnImmediately = .true.
       endif
    else
       returnImmediately = .true.
    endif

    if( returnImmediately ) return

    ! Determine the turbulence model and call the appropriate
    ! routine to compute the eddy viscosity.
    if (includeHalos) then
       iBeg = 1
       iEnd = ie
       jBeg = 1
       jEnd = je
       kBeg = 1
       kEnd = ke
    else
       iBeg = 2
       iEnd = il
       jBeg = 2
       jEnd = jl
       kBeg = 2
       kEnd = kl
    end if

    select case (turbModel)

    case (spalartAllmaras, spalartAllmarasEdwards)
       call saEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
#ifndef USE_TAPENADE

    case (v2f)
       call vfEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
    case (komegaWilcox, komegaModified)
       call kwEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)

    case (menterSST)
       call SSTEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)

    case (ktau)
       call ktEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
#endif
    end select

  end subroutine computeEddyViscosity

  subroutine saEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
    !
    !       saEddyViscosity computes the eddy-viscosity according to the
    !       Spalart-Allmaras model for the block given in blockPointers.
    !       This routine for both the original version as well as the
    !       modified version according to Edwards.
    !
    use constants
    use blockPointers
    use constants
    use paramTurb
    implicit none
    ! Input variables
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ii, iSize, jSize, kSize
    real(kind=realType)   :: chi, chi3, fv1, rnuSA, cv13

    ! Store the cv1^3; cv1 is a constant of the Spalart-Allmaras model.

    cv13 = rsaCv1**3

    ! Loop over the cells of this block and compute the eddy viscosity.
    ! Do not include halo's.
#ifdef TAPENADE_REVERSE
    iSize = (iEnd-iBeg)+1
    jSize = (jEnd-jBeg)+1
    kSize = (kEnd-kBeg)+1

    !$AD II-LOOP
    do ii=0, iSize*jSize*kSize-1
       i = mod(ii, iSize) + iBeg
       j = mod(ii/iSize, jSize) + jBeg
       k = ii/((iSize*jSize)) + kBeg
#else
       do k=kBeg, kEnd
          do j=jBeg, jEnd
             do i=iBeg, iEnd
#endif
                rnuSA      = w(i,j,k,itu1)*w(i,j,k,irho)
                chi        = rnuSA/rlv(i,j,k)
                chi3       = chi**3
                fv1        = chi3/(chi3+cv13)
                rev(i,j,k) = fv1*rnuSA
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif
  end subroutine saEddyViscosity

  subroutine kwEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
    !
    !       kwEddyViscosity computes the eddy viscosity according to the
    !       k-omega models (both the original Wilcox as well as the
    !       modified version) for the block given in blockPointers.
    !
    use constants
    use blockPointers
    implicit none
    ! Input variables
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ii, iSize, jSize, kSize


    ! Loop over the cells of this block and compute the eddy viscosity.
    ! Do not include halo's.
#ifdef TAPENADE_REVERSE
    iSize = (iEnd-iBeg)+1
    jSize = (jEnd-jBeg)+1
    kSize = (kEnd-kBeg)+1

    !$AD II-LOOP
    do ii=0, iSize*jSize*kSize-1
       i = mod(ii, iSize) + iBeg
       j = mod(ii/iSize, jSize) + jBeg
       k = ii/((iSize*jSize)) + kBeg
#else
       do k=kBeg, kEnd
          do j=jBeg, jEnd
             do i=iBeg, iEnd
#endif
                rev(i,j,k) = abs(w(i,j,k,irho)*w(i,j,k,itu1)/w(i,j,k,itu2))
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif

  end subroutine kwEddyViscosity

  subroutine SSTEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
    !
    !       SSTEddyViscosity computes the eddy viscosity according to
    !       menter's SST variant of the k-omega turbulence model for the
    !       block given in blockPointers.
    !
    use constants
    use blockPointers
    use paramTurb
    use turbMod
    implicit none
    ! Input variables
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ii, iSize, jSize, kSize
    real(kind=realType) :: t1, t2, arg2, f2, vortMag

    ! Compute the vorticity squared in the cell centers. The reason
    ! for computing the vorticity squared is that a routine exists
    ! for it; for the actual eddy viscosity computation the vorticity
    ! itself is needed.

    call prodWmag2

    ! Loop over the cells of this block and compute the eddy viscosity.
    ! Do not include halo's.
#ifdef TAPENADE_REVERSE
    iSize = (iEnd-iBeg)+1
    jSize = (jEnd-jBeg)+1
    kSize = (kEnd-kBeg)+1

    !$AD II-LOOP
    do ii=0, iSize*jSize*kSize-1
       i = mod(ii, iSize) + iBeg
       j = mod(ii/iSize, jSize) + jBeg
       k = ii/((iSize*jSize)) + kBeg
#else
       do k=kBeg, kEnd
          do j=jBeg, jEnd
             do i=iBeg, iEnd
#endif
                ! Compute the value of the function f2, which occurs in the
                ! eddy-viscosity computation.

                t1 = two*sqrt(w(i,j,k,itu1)) &
                     / (0.09_realType*w(i,j,k,itu2)*d2Wall(i,j,k))
                t2 = 500.0_realType*rlv(i,j,k) &
                     / (w(i,j,k,irho)*w(i,j,k,itu2)*d2Wall(i,j,k)**2)

                arg2 = max(t1,t2)
                f2   = tanh(arg2**2)

                ! And compute the eddy viscosity.

                vortMag    = sqrt(scratch(i,j,k,iprod))
                rev(i,j,k) = w(i,j,k,irho)*rSSTA1*w(i,j,k,itu1) &
                     / max(rSSTA1*w(i,j,k,itu2), f2*vortMag)
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif

  end subroutine SSTEddyViscosity

  subroutine turbAdvection(mAdv, nAdv, offset, qq)
    !
    !       turbAdvection discretizes the advection part of the turbulent
    !       transport equations. As the advection part is the same for all
    !       models, this generic routine can be used. Both the
    !       discretization and the central jacobian are computed in this
    !       subroutine. The former can either be 1st or 2nd order
    !       accurate; the latter is always based on the 1st order upwind
    !       discretization. When the discretization must be second order
    !       accurate, the fully upwind (kappa = -1) scheme in combination
    !       with the minmod limiter is used.
    !       Only nAdv equations are treated, while the actual system has
    !       size mAdv. The reason is that some equations for some
    !       turbulence equations do not have an advection part, e.g. the
    !       f equation in the v2-f model. The argument offset indicates
    !       the offset in the w vector where this subsystem starts. As a
    !       consequence it is assumed that the indices of the current
    !       subsystem are contiguous, e.g. if a 2*2 system is solved the
    !       Last index in w is offset+1 and offset+2 respectively.
    !
    use constants
    use blockPointers, only : nx, ny, nz, il, jl, kl, vol, sfaceI, sfaceJ, sfaceK, &
         w, si, sj, sk, addGridVelocities, bmti1, bmti2, bmtj1, bmtj2, &
         bmtk1, bmtk2, scratch
    use inputDiscretization, only : orderTurb
    use iteration, only : groundLevel
    use turbMod, only : secondOrd
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nAdv, mAdv, offset

    real(kind=realType), dimension(2:il,2:jl,2:kl,mAdv,mAdv), &
         intent(inout) :: qq
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ii, jj, kk, iii

    real(kind=realType) :: qs, voli, xa, ya, za
    real(kind=realType) :: uu, dwt, dwtm1, dwtp1, dwti, dwtj, dwtk

    real(kind=realType), dimension(mAdv) :: impl

    ! Determine whether or not a second order discretization for the
    ! advective terms must be used.
    secondOrd = .false.
    if(groundLevel == 1_intType .and. &
         orderTurb == secondOrder) secondOrd = .true.

    ! Initialize the grid velocity to zero. This value will be used
    ! if the block is not moving.
    continue
    !$AD CHECKPOINT-START
    qs = zero
    !
    !       Upwind discretization of the convective term in k (zeta)
    !       direction. Either the 1st order upwind or the second order
    !       fully upwind interpolation scheme, kappa = -1, is used in
    !       combination with the minmod limiter.
    !       The possible grid velocity must be taken into account.
    !
#ifdef TAPENADE_REVERSE
    !$AD II-LOOP
    do iii=0,nx*ny*nz-1
       i = mod(iii, nx) + 2
       j = mod(iii/nx, ny) + 2
       k = iii/(nx*ny) + 2
#else
       do k=2, kl
          do j=2, jl
             do i=2, il
#endif
                ! Compute the grid velocity if present.
                ! It is taken as the average of k and k-1,

                voli = half/vol(i,j,k)
                if( addGridVelocities ) &
                     qs = (sFaceK(i,j,k) + sFaceK(i,j,k-1))*voli

                ! Compute the normal velocity, where the normal direction
                ! is taken as the average of faces k and k-1.

                xa = (sk(i,j,k,1) + sk(i,j,k-1,1))*voli
                ya = (sk(i,j,k,2) + sk(i,j,k-1,2))*voli
                za = (sk(i,j,k,3) + sk(i,j,k-1,3))*voli

                uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs
                ! This term has unit: velocity/length

                ! Determine the situation we are having here, i.e. positive
                ! or negative normal velocity.

                velKdir: if(uu > zero) then

                   ! Velocity has a component in positive k-direction.
                   ! Loop over the number of advection equations.
                   !$AD II-LOOP
                   do ii=1,nAdv

                      ! Set the value of jj such that it corresponds to the
                      ! turbulent entry in w.

                      jj = ii + offset

                      ! Check whether a first or a second order discretization
                      ! must be used.

                      if( secondOrd ) then

                         ! Second order; store the three differences for the
                         ! discretization of the derivative in k-direction.

                         dwtm1 = w(i,j,k-1,jj) - w(i,j,k-2,jj)
                         dwt   = w(i,j,k,  jj) - w(i,j,k-1,jj)
                         dwtp1 = w(i,j,k+1,jj) - w(i,j,k,  jj)

                         ! Construct the derivative in this cell center. This
                         ! is the first order upwind derivative with two
                         ! nonlinear corrections.

                         dwtk = dwt

                         if(dwt*dwtp1 > zero) then
                            if(abs(dwt) < abs(dwtp1)) then
                               dwtk = dwtk + half*dwt
                            else
                               dwtk = dwtk + half*dwtp1
                            endif
                         endif

                         if(dwt*dwtm1 > zero) then
                            if(abs(dwt) < abs(dwtm1)) then
                               dwtk = dwtk - half*dwt
                            else
                               dwtk = dwtk - half*dwtm1
                            endif
                         endif

                      else

                         ! 1st order upwind scheme.

                         dwtk = w(i,j,k,jj) - w(i,j,k-1,jj)

                      endif

                      ! Update the residual. The convective term must be
                      ! substracted, because it appears on the other side of
                      ! the equation as the source and viscous terms.
                      ! uu*dwtk = (V.dot.face_normal)*delta(nuTilde)/delta(x)

                      scratch(i,j,k,idvt+ii-1) = scratch(i,j,k,idvt+ii-1) - uu*dwtk
#ifndef USE_TAPENADE
                      ! Update the central jacobian. First the term which is
                      ! always present, i.e. uu.

                      qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) + uu

                      ! For boundary cells k == 2, the implicit treatment must
                      ! be taken into account. Note that the implicit part
                      ! is only based on the 1st order discretization.
                      ! To improve stability the diagonal term is only taken
                      ! into account when it improves stability, i.e. when
                      ! it is positive.

                      if(k == 2) then
                         do kk=1,mAdv
                            impl(kk) = bmtk1(i,j,jj,kk+offset)
                         enddo

                         impl(ii) = max(impl(ii),zero)

                         do kk=1,mAdv
                            qq(i,j,k,ii,kk) = qq(i,j,k,ii,kk) + uu*impl(kk)
                         enddo
                      endif
#endif

                   enddo

                else velKdir

                   ! Velocity has a component in negative k-direction.
                   ! Loop over the number of advection equations
                   !$AD II-LOOP
                   do ii=1,nAdv

                      ! Set the value of jj such that it corresponds to the
                      ! turbulent entry in w.

                      jj = ii + offset

                      ! Check whether a first or a second order discretization
                      ! must be used.

                      if( secondOrd ) then

                         ! Store the three differences for the discretization of
                         ! the derivative in k-direction.

                         dwtm1 = w(i,j,k,  jj) - w(i,j,k-1,jj)
                         dwt   = w(i,j,k+1,jj) - w(i,j,k,  jj)
                         dwtp1 = w(i,j,k+2,jj) - w(i,j,k+1,jj)

                         ! Construct the derivative in this cell center. This is
                         ! the first order upwind derivative with two nonlinear
                         ! corrections.

                         dwtk = dwt

                         if(dwt*dwtp1 > zero) then
                            if(abs(dwt) < abs(dwtp1)) then
                               dwtk = dwtk - half*dwt
                            else
                               dwtk = dwtk - half*dwtp1
                            endif
                         endif

                         if(dwt*dwtm1 > zero) then
                            if(abs(dwt) < abs(dwtm1)) then
                               dwtk = dwtk + half*dwt
                            else
                               dwtk = dwtk + half*dwtm1
                            endif
                         endif

                      else

                         ! 1st order upwind scheme.

                         dwtk = w(i,j,k+1,jj) - w(i,j,k,jj)

                      endif

                      ! Update the residual. The convective term must be
                      ! substracted, because it appears on the other side
                      ! of the equation as the source and viscous terms.

                      scratch(i,j,k,idvt+ii-1) = scratch(i,j,k,idvt+ii-1) - uu*dwtk

                      ! Update the central jacobian. First the term which is
                      ! always present, i.e. -uu.
#ifndef USE_TAPENADE
                      qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) - uu

                      ! For boundary cells k == kl, the implicit treatment must
                      ! be taken into account. Note that the implicit part
                      ! is only based on the 1st order discretization.
                      ! To improve stability the diagonal term is only taken
                      ! into account when it improves stability, i.e. when
                      ! it is positive.

                      if(k == kl) then
                         do kk=1,mAdv
                            impl(kk) = bmtk2(i,j,jj,kk+offset)
                         enddo

                         impl(ii) = max(impl(ii),zero)

                         do kk=1,mAdv
                            qq(i,j,k,ii,kk) = qq(i,j,k,ii,kk) - uu*impl(kk)
                         enddo
                      endif
#endif
                   enddo

                endif velKdir
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif
    continue
    !$AD CHECKPOINT-END
    !
    !       Upwind discretization of the convective term in j (eta)
    !       direction. Either the 1st order upwind or the second order
    !       fully upwind interpolation scheme, kappa = -1, is used in
    !       combination with the minmod limiter.
    !       The possible grid velocity must be taken into account.
    !
    continue
    !$AD CHECKPOINT-START
    qs = zero
#ifdef TAPENADE_REVERSE
    !$AD II-LOOP
    do iii=0,nx*ny*nz-1
       i = mod(iii, nx) + 2
       j = mod(iii/nx, ny) + 2
       k = iii/(nx*ny) + 2
#else
       do k=2, kl
          do j=2, jl
             do i=2, il
#endif

                ! Compute the grid velocity if present.
                ! It is taken as the average of j and j-1,

                voli = half/vol(i,j,k)
                if( addGridVelocities ) &
                     qs = (sFaceJ(i,j,k) + sFaceJ(i,j-1,k))*voli

                ! Compute the normal velocity, where the normal direction
                ! is taken as the average of faces j and j-1.

                xa = (sj(i,j,k,1) + sj(i,j-1,k,1))*voli
                ya = (sj(i,j,k,2) + sj(i,j-1,k,2))*voli
                za = (sj(i,j,k,3) + sj(i,j-1,k,3))*voli

                uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs

                ! Determine the situation we are having here, i.e. positive
                ! or negative normal velocity.

                velJdir: if(uu > zero) then

                   ! Velocity has a component in positive j-direction.
                   ! Loop over the number of advection equations.
                   !$AD II-LOOP
                   do ii=1,nAdv

                      ! Set the value of jj such that it corresponds to the
                      ! turbulent entry in w.

                      jj = ii + offset

                      ! Check whether a first or a second order discretization
                      ! must be used.

                      if( secondOrd ) then

                         ! Second order; store the three differences for the
                         ! discretization of the derivative in j-direction.

                         dwtm1 = w(i,j-1,k,jj) - w(i,j-2,k,jj)
                         dwt   = w(i,j,  k,jj) - w(i,j-1,k,jj)
                         dwtp1 = w(i,j+1,k,jj) - w(i,j,  k,jj)

                         ! Construct the derivative in this cell center. This is
                         ! the first order upwind derivative with two nonlinear
                         ! corrections.

                         dwtj = dwt

                         if(dwt*dwtp1 > zero) then
                            if(abs(dwt) < abs(dwtp1)) then
                               dwtj = dwtj + half*dwt
                            else
                               dwtj = dwtj + half*dwtp1
                            endif
                         endif

                         if(dwt*dwtm1 > zero) then
                            if(abs(dwt) < abs(dwtm1)) then
                               dwtj = dwtj - half*dwt
                            else
                               dwtj = dwtj - half*dwtm1
                            endif
                         endif

                      else

                         ! 1st order upwind scheme.

                         dwtj = w(i,j,k,jj) - w(i,j-1,k,jj)

                      endif

                      ! Update the residual. The convective term must be
                      ! substracted, because it appears on the other side of
                      ! the equation as the source and viscous terms.

                      scratch(i,j,k,idvt+ii-1) = scratch(i,j,k,idvt+ii-1) - uu*dwtj

                      ! Update the central jacobian. First the term which is
                      ! always present, i.e. uu.
#ifndef USE_TAPENADE
                      qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) + uu

                      ! For boundary cells j == 2, the implicit treatment must
                      ! be taken into account. Note that the implicit part
                      ! is only based on the 1st order discretization.
                      ! To improve stability the diagonal term is only taken
                      ! into account when it improves stability, i.e. when
                      ! it is positive.

                      if(j == 2) then
                         do kk=1,mAdv
                            impl(kk) = bmtj1(i,k,jj,kk+offset)
                         enddo

                         impl(ii) = max(impl(ii),zero)

                         do kk=1,mAdv
                            qq(i,j,k,ii,kk) = qq(i,j,k,ii,kk) + uu*impl(kk)
                         enddo
                      endif
#endif
                   enddo

                else velJdir

                   ! Velocity has a component in negative j-direction.
                   ! Loop over the number of advection equations.
                   !$AD II-LOOP
                   do ii=1,nAdv

                      ! Set the value of jj such that it corresponds to the
                      ! turbulent entry in w.

                      jj = ii + offset

                      ! Check whether a first or a second order discretization
                      ! must be used.

                      if( secondOrd ) then

                         ! Store the three differences for the discretization of
                         ! the derivative in j-direction.

                         dwtm1 = w(i,j,  k,jj) - w(i,j-1,k,jj)
                         dwt   = w(i,j+1,k,jj) - w(i,j,  k,jj)
                         dwtp1 = w(i,j+2,k,jj) - w(i,j+1,k,jj)

                         ! Construct the derivative in this cell center. This is
                         ! the first order upwind derivative with two nonlinear
                         ! corrections.

                         dwtj = dwt

                         if(dwt*dwtp1 > zero) then
                            if(abs(dwt) < abs(dwtp1)) then
                               dwtj = dwtj - half*dwt
                            else
                               dwtj = dwtj - half*dwtp1
                            endif
                         endif

                         if(dwt*dwtm1 > zero) then
                            if(abs(dwt) < abs(dwtm1)) then
                               dwtj = dwtj + half*dwt
                            else
                               dwtj = dwtj + half*dwtm1
                            endif
                         endif

                      else

                         ! 1st order upwind scheme.

                         dwtj = w(i,j+1,k,jj) - w(i,j,k,jj)

                      endif

                      ! Update the residual. The convective term must be
                      ! substracted, because it appears on the other side
                      ! of the equation as the source and viscous terms.

                      scratch(i,j,k,idvt+ii-1) = scratch(i,j,k,idvt+ii-1) - uu*dwtj

                      ! Update the central jacobian. First the term which is
                      ! always present, i.e. -uu.
#ifndef USE_TAPENADE
                      qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) - uu

                      ! For boundary cells j == jl, the implicit treatment must
                      ! be taken into account. Note that the implicit part
                      ! is only based on the 1st order discretization.
                      ! To improve stability the diagonal term is only taken
                      ! into account when it improves stability, i.e. when
                      ! it is positive.

                      if(j == jl) then
                         do kk=1,mAdv
                            impl(kk) = bmtj2(i,k,jj,kk+offset)
                         enddo

                         impl(ii) = max(impl(ii),zero)

                         do kk=1,mAdv
                            qq(i,j,k,ii,kk) = qq(i,j,k,ii,kk) - uu*impl(kk)
                         enddo
                      endif
#endif
                   enddo

                endif velJdir
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif
    continue
    !$AD CHECKPOINT-END
    !
    !       Upwind discretization of the convective term in i (xi)
    !       direction. Either the 1st order upwind or the second order
    !       fully upwind interpolation scheme, kappa = -1, is used in
    !       combination with the minmod limiter.
    !       The possible grid velocity must be taken into account.
    !
    continue
    !$AD CHECKPOINT-START
    qs = zero
#ifdef TAPENADE_REVERSE
    !$AD II-LOOP
    do iii=0,nx*ny*nz-1
       i = mod(iii, nx) + 2
       j = mod(iii/nx, ny) + 2
       k = iii/(nx*ny) + 2
#else
       do k=2, kl
          do j=2, jl
             do i=2, il
#endif


                ! Compute the grid velocity if present.
                ! It is taken as the average of i and i-1,

                voli = half/vol(i,j,k)
                if( addGridVelocities ) &
                     qs = (sFaceI(i,j,k) + sFaceI(i-1,j,k))*voli

                ! Compute the normal velocity, where the normal direction
                ! is taken as the average of faces i and i-1.

                xa = (si(i,j,k,1) + si(i-1,j,k,1))*voli
                ya = (si(i,j,k,2) + si(i-1,j,k,2))*voli
                za = (si(i,j,k,3) + si(i-1,j,k,3))*voli

                uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs

                ! Determine the situation we are having here, i.e. positive
                ! or negative normal velocity.

                velIdir: if(uu > zero) then

                   ! Velocity has a component in positive i-direction.
                   ! Loop over the number of advection equations.
                   !$AD II-LOOP
                   do ii=1,nAdv

                      ! Set the value of jj such that it corresponds to the
                      ! turbulent entry in w.

                      jj = ii + offset

                      ! Check whether a first or a second order discretization
                      ! must be used.

                      if( secondOrd ) then

                         ! Second order; store the three differences for the
                         ! discretization of the derivative in i-direction.

                         dwtm1 = w(i-1,j,k,jj) - w(i-2,j,k,jj)
                         dwt   = w(i,  j,k,jj) - w(i-1,j,k,jj)
                         dwtp1 = w(i+1,j,k,jj) - w(i,  j,k,jj)

                         ! Construct the derivative in this cell center. This is
                         ! the first order upwind derivative with two nonlinear
                         ! corrections.

                         dwti = dwt

                         if(dwt*dwtp1 > zero) then
                            if(abs(dwt) < abs(dwtp1)) then
                               dwti = dwti + half*dwt
                            else
                               dwti = dwti + half*dwtp1
                            endif
                         endif

                         if(dwt*dwtm1 > zero) then
                            if(abs(dwt) < abs(dwtm1)) then
                               dwti = dwti - half*dwt
                            else
                               dwti = dwti - half*dwtm1
                            endif
                         endif

                      else

                         ! 1st order upwind scheme.

                         dwti = w(i,j,k,jj) - w(i-1,j,k,jj)

                      endif

                      ! Update the residual. The convective term must be
                      ! substracted, because it appears on the other side of
                      ! the equation as the source and viscous terms.

                      scratch(i,j,k,idvt+ii-1) = scratch(i,j,k,idvt+ii-1) - uu*dwti

                      ! Update the central jacobian. First the term which is
                      ! always present, i.e. uu.
#ifndef USE_TAPENADE
                      qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) + uu

                      ! For boundary cells i == 2, the implicit treatment must
                      ! be taken into account. Note that the implicit part
                      ! is only based on the 1st order discretization.
                      ! To improve stability the diagonal term is only taken
                      ! into account when it improves stability, i.e. when
                      ! it is positive.

                      if(i == 2) then
                         do kk=1,mAdv
                            impl(kk) = bmti1(j,k,jj,kk+offset)
                         enddo

                         impl(ii) = max(impl(ii),zero)

                         do kk=1,mAdv
                            qq(i,j,k,ii,kk) = qq(i,j,k,ii,kk) + uu*impl(kk)
                         enddo
                      endif
#endif
                   enddo

                else velIdir

                   ! Velocity has a component in negative i-direction.
                   ! Loop over the number of advection equations.
                   !$AD II-LOOP
                   do ii=1,nAdv

                      ! Set the value of jj such that it corresponds to the
                      ! turbulent entry in w.

                      jj = ii + offset

                      ! Check whether a first or a second order discretization
                      ! must be used.

                      if( secondOrd ) then

                         ! Second order; store the three differences for the
                         ! discretization of the derivative in i-direction.

                         dwtm1 = w(i,  j,k,jj) - w(i-1,j,k,jj)
                         dwt   = w(i+1,j,k,jj) - w(i,  j,k,jj)
                         dwtp1 = w(i+2,j,k,jj) - w(i+1,j,k,jj)

                         ! Construct the derivative in this cell center. This is
                         ! the first order upwind derivative with two nonlinear
                         ! corrections.

                         dwti = dwt

                         if(dwt*dwtp1 > zero) then
                            if(abs(dwt) < abs(dwtp1)) then
                               dwti = dwti - half*dwt
                            else
                               dwti = dwti - half*dwtp1
                            endif
                         endif

                         if(dwt*dwtm1 > zero) then
                            if(abs(dwt) < abs(dwtm1)) then
                               dwti = dwti + half*dwt
                            else
                               dwti = dwti + half*dwtm1
                            endif
                         endif

                      else

                         ! 1st order upwind scheme.

                         dwti = w(i+1,j,k,jj) - w(i,j,k,jj)

                      endif

                      ! Update the residual. The convective term must be
                      ! substracted, because it appears on the other side
                      ! of the equation as the source and viscous terms.

                      scratch(i,j,k,idvt+ii-1) = scratch(i,j,k,idvt+ii-1) - uu*dwti

                      ! Update the central jacobian. First the term which is
                      ! always present, i.e. -uu.
#ifndef USE_TAPENADE
                      qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) - uu

                      ! For boundary cells i == il, the implicit treatment must
                      ! be taken into account. Note that the implicit part
                      ! is only based on the 1st order discretization.
                      ! To improve stability the diagonal term is only taken
                      ! into account when it improves stability, i.e. when
                      ! it is positive.

                      if(i == il) then
                         do kk=1,mAdv
                            impl(kk) = bmti2(j,k,jj,kk+offset)
                         enddo

                         impl(ii) = max(impl(ii),zero)

                         do kk=1,mAdv
                            qq(i,j,k,ii,kk) = qq(i,j,k,ii,kk) - uu*impl(kk)
                         enddo
                      endif
#endif
                   enddo

                endif velIdir
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif

    !$AD CHECKPOINT-END
    continue
  end subroutine turbAdvection


  ! ----------------------------------------------------------------------
  !                                                                      |
  !                    No Tapenade Routine below this line               |
  !                                                                      |
  ! ----------------------------------------------------------------------

#ifndef USE_TAPENADE

  subroutine tdia3(nb, ne, l, c, u, r)
    !
    !       tdia3 solves the tridiagonal linear system (l+c+u) v = r,
    !       where l is the lower, c the central and u the upper diagonal.
    !       Every entry in the matrix is a 2x2 block matrix, i.e.
    !                    x x  x 0  0 0  ........        = c(nb) u(nb)
    !                    x x  0 x  0 0  ........
    !                    x 0  x x  x 0  ........        = l(i) c(i) u(i)
    !                    0 x  x x  0 x  ........
    !                         ........  x 0  x x        = l(ne) c(ne)
    !                         ........  0 x  x x
    !               With c = x x     u,l = x 0
    !                        x x           0 x
    !
    use constants
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nb, ne

    real(kind=realType), dimension(2,nb:ne), intent(inout) :: l, u, r
    real(kind=realType), dimension(2,2,nb:ne), intent(inout) :: c
    !
    !      Local variables.
    !
    integer(kind=intType) :: n
    real(kind=realType)   :: deti, f11, f12, f21, f22, r1

    ! Perform the backward sweep to eliMinate the upper diagonal uu.
    ! f     = u(n)*c^-1(n+1),
    ! c'(n) = c(n) - f*l(n+1)
    ! r'(n) = r(n) - f*r(n+1)

    do n=ne-1,nb,-1

       deti =  one/(c(1,1,n+1)*c(2,2,n+1) - c(1,2,n+1)*c(2,1,n+1))
       f11  =  u(1,n)*c(2,2,n+1)*deti
       f12  = -u(1,n)*c(1,2,n+1)*deti
       f21  = -u(2,n)*c(2,1,n+1)*deti
       f22  =  u(2,n)*c(1,1,n+1)*deti

       c(1,1,n) = c(1,1,n) - f11*l(1,n+1)
       c(1,2,n) = c(1,2,n) - f12*l(2,n+1)
       c(2,1,n) = c(2,1,n) - f21*l(1,n+1)
       c(2,2,n) = c(2,2,n) - f22*l(2,n+1)

       r(1,n)  = r(1,n) - f11*r(1,n+1) - f12*r(2,n+1)
       r(2,n)  = r(2,n) - f21*r(1,n+1) - f22*r(2,n+1)

    enddo

    ! The matrix is now in low block bi-diagonal form and can thus
    ! be solved be a forward sweep. The solution is stored in r.

    deti    =  one/(c(1,1,nb)*c(2,2,nb) - c(1,2,nb)*c(2,1,nb))
    r1      =  r(1,nb)
    r(1,nb) =  deti*(c(2,2,nb)*r1 - c(1,2,nb)*r(2,nb))
    r(2,nb) = -deti*(c(2,1,nb)*r1 - c(1,1,nb)*r(2,nb))

    do n=nb+1,ne

       r(1,n) = r(1,n) - l(1,n)*r(1,n-1)
       r(2,n) = r(2,n) - l(2,n)*r(2,n-1)

       deti   =  one/(c(1,1,n)*c(2,2,n) - c(1,2,n)*c(2,1,n))
       r1     =  r(1,n)
       r(1,n) =  deti*(c(2,2,n)*r1 - c(1,2,n)*r(2,n))
       r(2,n) = -deti*(c(2,1,n)*r1 - c(1,1,n)*r(2,n))

    enddo

  end subroutine tdia3

  subroutine vfEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
    !
    !       vfEddyViscosity computes the eddy-viscosity according to the
    !       v2f turbulence model for the block given in blockPointers.
    !       This routine is for both the n=1 and n=6 version.
    !
    use constants
    use blockPointers
    use constants
    use paramTurb
    use turbMod
    use inputPhysics
    use turbCurveFits, only : curvetupyp
    implicit none
    ! Input variables
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, nn, ii, iSize, jSize, kSize
    real(kind=realType)   :: tke, tep, tkea, tepa, tepl, tv2, tv2a
    real(kind=realType)   :: yp, utau

    real(kind=realType), dimension(itu1:itu5) :: tup
    real(kind=realType), dimension(:,:,:), pointer :: ww
    real(kind=realType), dimension(:,:),   pointer :: rrlv, rrev
    real(kind=realType), dimension(:,:),   pointer :: dd2Wall


    ! Compute time and length scale

    call vfScale

    ! Loop over the cells of this block and compute the eddy viscosity.
    ! Do not include halo's.
#ifdef TAPENADE_REVERSE
    iSize = (iEnd-iBeg)+1
    jSize = (jEnd-jBeg)+1
    kSize = (kEnd-kBeg)+1

    !$AD II-LOOP
    do ii=0, iSize*jSize*kSize-1
       i = mod(ii, iSize) + iBeg
       j = mod(ii/iSize, jSize) + jBeg
       k = ii/((iSize*jSize)) + kBeg
#else
       do k=kBeg, kEnd
          do j=jBeg, jEnd
             do i=iBeg, iEnd
#endif

                tke   = w(i,j,k,itu1)
                tep   = w(i,j,k,itu2)
                tv2   = w(i,j,k,itu3)
                tkea  = abs(tke)
                tepa  = abs(tep)
                tv2a  = abs(tv2)
                tepl  = max(tepa,rvfLimitE)

                rev(i,j,k) = rvfCmu*w(i,j,k,irho)*tv2a/tepl*sct(i,j,k)
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif

    ! Modify the rhs of the 1st internal cell, if wall functions
    ! are used; their value is determined by the table.

    testWallFunctions: if( wallFunctions ) then

       bocos: do nn=1,nViscBocos

          ! Determine the block face on which the subface is located
          ! and set some variables. As flag points to the entire array
          ! flagI2, etc., its starting indices are the starting indices
          ! of its target and not 1.

          select case (bcFaceid(nn))
          case (iMin)
             ww      => w(2,1:,1:,1:);   rrlv => rlv(2,1:,1:)
             dd2Wall => d2Wall(2,:,:);   rrev => rev(2,1:,1:)

          case (iMax)
             ww      => w(il,1:,1:,1:);   rrlv => rlv(il,1:,1:)
             dd2Wall => d2Wall(il,:,:);   rrev => rev(il,1:,1:)

          case (jMin)
             ww      => w(1:,2,1:,1:);   rrlv => rlv(1:,2,1:)
             dd2Wall => d2Wall(:,2,:);   rrev => rev(1:,2,1:)

          case (jMax)
             ww      => w(1:,jl,1:,1:);   rrlv => rlv(1:,jl,1:)
             dd2Wall => d2Wall(:,jl,:);   rrev => rev(1:,jl,1:)

          case (kMin)
             ww      => w(1:,1:,2,1:);   rrlv => rlv(1:,1:,2)
             dd2Wall => d2Wall(:,:,2);   rrev => rev(1:,1:,2)

          case (kMax)
             ww      => w(1:,1:,kl,1:);   rrlv => rlv(1:,1:,kl)
             dd2Wall => d2Wall(:,:,kl);   rrev => rev(1:,1:,kl)

          end select

          ! Loop over the owned faces of this subface. Therefore the
          ! nodal range of bcData must be used. The offset of +1 is
          ! present, because the starting index of the cell range is
          ! 1 larger than the starting index of the nodal range.

          do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
             do i=(BCData(nn)%inBeg+1),BCData(nn)%inEnd

                ! Enforce k and epsilon in the 1st internal cell from
                ! the wall function table. There is an offset of -1 in
                ! the wall distance. Note that the offset compared to
                ! the current value must be stored. Also note that the
                ! curve fits contain the non-dimensional values.

                utau = viscSubface(nn)%utau(i,j)
                yp = ww(i,j,irho)*dd2Wall(i-1,j-1)*utau/rrlv(i,j)

                call curveTupYp(tup(itu5:itu5), yp, itu5, itu5)
                rrev(i,j) = tup(itu5)*rrlv(i,j)
             enddo
          enddo

       enddo bocos
    endif testWallFunctions

  end subroutine vfEddyViscosity
  subroutine vfScale
    !
    !       time and length scale definition for v2f turbulence model. The
    !       upper bound can be switched on by setting rvfB to .true.
    !       The strain squared is defined as: strain2 = 2 sij sij
    !
    use constants
    use blockPointers
    use inputPhysics
    use paramTurb
    use turbMod
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k

    real(kind=realType) :: sqrt3
    real(kind=realType) :: tkea, tepa, tv2a, supi, rn2
    real(kind=realType) :: rsct, rscl2, rnu, rstrain

    ! Some constants in the model.

    sqrt3 = sqrt(three)

    ! Set the pointer for dvt in dw, such that the code is more
    ! readable. Also set the pointers for the production term,
    ! vorticity, strain and the time and lenght scale of v2f.

    dvt     => scratch(1:,1:,1:,idvt:)
    prod    => scratch(1:,1:,1:,iprod)
    sct     => scratch(1:,1:,1:,isct)
    scl2    => scratch(1:,1:,1:,iscl2)
    vort    => prod
    strain2 => prod
    !
    !       Production term.
    !
    select case (turbProd)
    case (strain)
       call prodSmag2

    case (vorticity)
       call prodWmag2

    case (katoLaunder)
       call prodKatoLaunder

    end select
    !
    !       Compute the length and time scale for all internal cells.
    !
    if( rvfB ) then

       do k=2,kl
          do j=2,jl
             do i=2,il

                ! Compute the time and length scale with upper bound

                rstrain     = sqrt(strain2(i,j,k))
                rnu         = rlv(i,j,k)/w(i,j,k,irho)
                tkea        = abs(w(i,j,k,itu1))
                tepa        = abs(w(i,j,k,itu2))
                tv2a        = abs(w(i,j,k,itu3))
                supi        = tepa*tkea/max(sqrt3*tv2a*rvfCmu*rstrain,eps)
                rn2         = rvfCn**2*(rnu*tepa)**1.5_realType

                rsct        = max(tkea,six*sqrt(rnu*tepa))
                sct(i,j,k)  = min(rsct,0.6_realType*supi)
                rscl2       = tkea*min(tkea**2,supi**2)
                scl2(i,j,k) = rvfCl**2*max(rscl2,rn2)

             enddo
          enddo
       enddo

    else

       do k=2,kl
          do j=2,jl
             do i=2,il

                ! Compute the time and length scale without upper bound

                rnu         = rlv(i,j,k)/w(i,j,k,irho)
                tkea        = abs(w(i,j,k,itu1))
                tepa        = abs(w(i,j,k,itu2))
                rn2         = rvfCn**2*(rnu*tepa)**1.5_realType

                rsct        = max(tkea,six*sqrt(rnu*tepa))
                sct(i,j,k)  = rsct
                rscl2       = tkea**3
                scl2(i,j,k) = rvfCl**2*max(rscl2,rn2)
             enddo
          enddo
       enddo
    endif

  end subroutine vfScale

  subroutine ktEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd)
    !
    !       ktEddyViscosity computes the eddy viscosity according to the
    !       k-tau turbulence model for the block given in blockPointers.
    !
    use blockPointers
    use constants
    implicit none
    ! Input variables
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ii, iSize, jSize, kSize

    ! Loop over the cells of this block and compute the eddy viscosity.
    ! Do not include halo's.
#ifdef TAPENADE_REVERSE
    iSize = (iEnd-iBeg)+1
    jSize = (jEnd-jBeg)+1
    kSize = (kEnd-kBeg)+1

    !$AD II-LOOP
    do ii=0, iSize*jSize*kSize-1
       i = mod(ii, iSize) + iBeg
       j = mod(ii/iSize, jSize) + jBeg
       k = ii/((iSize*jSize)) + kBeg
#else
       do k=kBeg, kEnd
          do j=jBeg, jEnd
             do i=iBeg, iEnd
#endif
                rev(i,j,k) = abs(w(i,j,k,irho)*w(i,j,k,itu1)*w(i,j,k,itu2))
#ifdef TAPENADE_REVERSE
             end do
#else
          enddo
       enddo
    enddo
#endif
  end subroutine ktEddyViscosity


#ifndef USE_TAPENADE
  subroutine unsteadyTurbSpectral(ntu1, ntu2)
    use constants
    use blockPointers
    use inputPhysics
    use inputTimeSpectral
    use iteration
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: ntu1, ntu2
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, sps

    ! Return immediately if not the time spectral equations are to
    ! be solved.

    if(equationMode /= timeSpectral) return

    ! Loop over the number of spectral modes and local blocks.

    spectralLoop: do sps=1,nTimeIntervalsSpectral
       domains: do nn=1,nDom

          ! Set the pointers for this block.

          call setPointers(nn, currentLevel, sps)

          call unsteadyTurbSpectral_block(ntu1, ntu2, nn, sps)
       end do domains
    end do spectralLoop
  end subroutine unsteadyTurbSpectral
#endif

  subroutine unsteadyTurbSpectral_block(ntu1, ntu2, nn, sps)
    !
    !       unsteadyTurbSpectral determines the spectral time derivative
    !       for all owned cells. This routine is called before the actual
    !       solve routines, such that the treatment is identical for all
    !       spectral solutions. The results is stored in the corresponding
    !       entry in dw.
    !
    use constants
    use blockPointers
    use inputPhysics
    use inputTimeSpectral
    use iteration
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: ntu1, ntu2, nn, sps
    !
    !      Local variables.
    !
    integer(kind=intType) :: ii, mm, i, j, k

    ! Return immediately if not the time spectral equations are to
    ! be solved.

    if(equationMode /= timeSpectral) return

    ! Loop over the number of turbulent transport equations.

    nAdvLoop: do ii=ntu1, ntu2

       ! Initialize the time derivative to zero for the owned
       ! cell centers.

       do k=2,kl
          do j=2,jl
             do i=2,il
                dw(i,j,k,ii) = zero
             enddo
          enddo
       enddo

       ! Loop over the number of terms which contribute to the
       ! time derivative.

       do mm=1,nTimeIntervalsSpectral

          ! Add the contribution to the time derivative for
          ! all owned cells.

          do k=2,kl
             do j=2,jl
                do i=2,il
                   dw(i,j,k,ii) = dw(i,j,k,ii)              &
                        + dscalar(sectionID,sps,mm) &
                        * flowDoms(nn,currentLevel,mm)%w(i,j,k,ii)
                enddo
             enddo
          enddo

       enddo

    enddo nAdvLoop
  end subroutine unsteadyTurbSpectral_block

  subroutine initKOmega(pOffset)
    !
    !       initKOmega initializes the values of k and omega a bit more
    !       intelligently than just free-stream values.
    !       It is assumed that the pointers in blockPointers already
    !       point to the correct block on the correct level.
    !       The argument pOffset is present such that in case of restart
    !       a possible pointer offset is taken into account. For more
    !       details see the corresponding routines in initFlow.
    !
    use blockPointers
    use constants
    use flowVarRefState
    use paramTurb
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: pOffset
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ip, jp, kp
    real(kind=realType)   :: rhoi, tmp

    ! Loop over the owned cells of the block.

    do k=2,kl
       kp = k + pOffset
       do j=2,jl
          jp = j + pOffset
          do i=2,il
             ip = i + pOffset

             ! Compute a value of omega based on the wall distance.

             rhoi = one/w(ip,jp,kp,irho)
             tmp  = six*rhoi*rlv(i,j,k)/(rkwBeta1*(d2Wall(i,j,k)**2))

             ! Initialize omega using the value just computed; make sure
             ! that the free stream value is the lowest possible value.
             ! After that initialize k using the current value of omega
             ! and the eddy viscosity. Again the free-stream value is
             ! the lower bound.

             w(ip,jp,kp,itu2) = max(tmp,wInf(itu2))
             tmp              = rhoi*rev(i,j,k)*w(ip,jp,kp,itu2)
             w(ip,jp,kp,itu1) = max(tmp,wInf(itu1))

          enddo
       enddo
    enddo

  end subroutine initKOmega

  subroutine kwCDterm
    !
    !       kwCDterm computes the cross-diffusion term in the omega-eqn
    !       for the SST version as well as the modified k-omega turbulence
    !       model. It is assumed that the pointers in blockPointers and
    !       turbMod are already set.
    !
    use constants
    use blockPointers
    use turbMod
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k
    real(kind=realType)   :: kx, ky, kz, wwx, wwy, wwz
    real(kind=realType)   :: lnwip1, lnwim1, lnwjp1, lnwjm1
    real(kind=realType)   :: lnwkp1, lnwkm1

    ! Loop over the cell centers of the given block. It may be more
    ! efficient to loop over the faces and to scatter the gradient,
    ! but in that case the gradients for k and omega must be stored.
    ! In the current approach no extra memory is needed.

    do k=2,kl
       do j=2,jl
          do i=2,il

             ! Compute the gradient of k in the cell center. Use is made
             ! of the fact that the surrounding normals sum up to zero,
             ! such that the cell i,j,k does not give a contribution.
             ! The gradient is scaled by a factor 1/2vol.

             kx = w(i+1,j,k,itu1)*si(i,j,k,1) - w(i-1,j,k,itu1)*si(i-1,j,k,1) &
                  + w(i,j+1,k,itu1)*sj(i,j,k,1) - w(i,j-1,k,itu1)*sj(i,j-1,k,1) &
                  + w(i,j,k+1,itu1)*sk(i,j,k,1) - w(i,j,k-1,itu1)*sk(i,j,k-1,1)
             ky = w(i+1,j,k,itu1)*si(i,j,k,2) - w(i-1,j,k,itu1)*si(i-1,j,k,2) &
                  + w(i,j+1,k,itu1)*sj(i,j,k,2) - w(i,j-1,k,itu1)*sj(i,j-1,k,2) &
                  + w(i,j,k+1,itu1)*sk(i,j,k,2) - w(i,j,k-1,itu1)*sk(i,j,k-1,2)
             kz = w(i+1,j,k,itu1)*si(i,j,k,3) - w(i-1,j,k,itu1)*si(i-1,j,k,3) &
                  + w(i,j+1,k,itu1)*sj(i,j,k,3) - w(i,j-1,k,itu1)*sj(i,j-1,k,3) &
                  + w(i,j,k+1,itu1)*sk(i,j,k,3) - w(i,j,k-1,itu1)*sk(i,j,k-1,3)

             ! Compute the logarithm of omega in the points that
             ! contribute to the gradient in this cell.

             lnwip1 = log(abs(w(i+1,j,k,itu2)))
             lnwim1 = log(abs(w(i-1,j,k,itu2)))
             lnwjp1 = log(abs(w(i,j+1,k,itu2)))
             lnwjm1 = log(abs(w(i,j-1,k,itu2)))
             lnwkp1 = log(abs(w(i,j,k+1,itu2)))
             lnwkm1 = log(abs(w(i,j,k-1,itu2)))

             ! Compute the scaled gradient of ln omega.

             wwx = lnwip1*si(i,j,k,1) - lnwim1*si(i-1,j,k,1) &
                  + lnwjp1*sj(i,j,k,1) - lnwjm1*sj(i,j-1,k,1) &
                  + lnwkp1*sk(i,j,k,1) - lnwkm1*sk(i,j,k-1,1)
             wwy = lnwip1*si(i,j,k,2) - lnwim1*si(i-1,j,k,2) &
                  + lnwjp1*sj(i,j,k,2) - lnwjm1*sj(i,j-1,k,2) &
                  + lnwkp1*sk(i,j,k,2) - lnwkm1*sk(i,j,k-1,2)
             wwz = lnwip1*si(i,j,k,3) - lnwim1*si(i-1,j,k,3) &
                  + lnwjp1*sj(i,j,k,3) - lnwjm1*sj(i,j-1,k,3) &
                  + lnwkp1*sk(i,j,k,3) - lnwkm1*sk(i,j,k-1,3)

             ! Compute the dot product grad k grad ln omega.
             ! Multiply it by the correct scaling factor and store it.

             kwCD(i,j,k) = fourth*(kx*wwx + ky*wwy + kz*wwz)/(vol(i,j,k)**2)

          enddo
       enddo
    enddo

  end subroutine kwCDterm
#endif
end module turbUtils
