module blockette

  use constants
  ! This temporary module contains all cache-blocked code. It also
  ! contains the statically allocated variables on which the blocked
  ! code operates. 

  ! Dummy Block dimensions
  integer(kind=intType), parameter :: BS=8
  integer(kind=intType), parameter :: bbil=BS+1, bbjl=BS+1, bbkl=BS+1
  integer(kind=intType), parameter :: bbie=BS+2, bbje=BS+2, bbke=BS+2 
  integer(kind=intType), parameter :: bbib=BS+3, bbjb=BS+3, bbkb=BS+3

  ! Actual dimensions to execute
  integer(kind=intType) :: nx, ny, nz, il, jl, kl, ie, je, ke, ib, jb, kb


  ! Double halos
  real(kind=realType), dimension(0:bbib, 0:bbjb, 0:bbkb, 1:6) :: w
  real(kind=realType), dimension(0:bbib, 0:bbjb, 0:bbkb) :: P, gamma
  real(kind=realType), dimension(0:bbib, 0:bbjb, 0:bbkb) :: ss ! Entropy

  ! Single halos
  real(kind=realType), dimension(0:bbie, 0:bbje, 0:bbke, 3) :: x
  real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke):: rlv, rev, vol, aa
  real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke) :: radI, radJ, radK
  real(kind=realType),dimension(1:bbie, 1:bbje, 1:bbke, 3) :: dss ! Shock sensor

  ! No halos
  real(kind=realType), dimension(2:bbil, 2:bbjl, 2:bbkl) :: volRef, d2wall
  integer(kind=intType), dimension(2:bbil, 2:bbjl, 2:bbkl) :: iblank

  ! Face Porosities
  integer(kind=porType), dimension(1:bbil, 2:bbjl, 2:bbkl) :: porI
  integer(kind=porType), dimension(2:bbil, 1:bbjl, 2:bbkl) :: porJ
  integer(kind=porType), dimension(2:bbil, 2:bbjl, 1:bbkl) :: porK

  ! Single halos (only owned cells significant)
  real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke, 1:5) :: fw
  real(kind=realType), dimension(1:bbie, 1:bbje, 1:bbke, 1:6) :: dw

  ! Face projected areas
  real(kind=realType), dimension(0:bbie, 1:bbje, 1:bbke, 3) :: sI
  real(kind=realType), dimension(1:bbie, 0:bbje, 1:bbke, 3) :: sJ
  real(kind=realType), dimension(1:bbie, 1:bbje, 0:bbke, 3) :: sK

  ! Nodal gradients
  real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: ux, uy, uz
  real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: vx, vy, vz
  real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: wx, wy, wz
  real(kind=realType), dimension(1:bbil, 1:bbjl, 1:bbkl) :: qx, qy, qz 

  ! Make *all* of these variables tread-private
  !$OMP THREADPRIVATE(nx, ny, nz, il, jl, kl, ie, je, ke, ib, jb, kb)
  !$OMP THREADPRIVATE(w, p, gamma, ss, x, rlv, rev, vol, aa, radI, radJ, radK)
  !$OMP THREADPRIVATE(dss, volRef, d2wall, iblank, porI, porJ, porK, fw, dw)
  !$OMP THREADPRIVATE(sI, sJ, sK, ux, uy, uz, vx, vy, vz, wx, wy, wz, qx, qy, qz)
contains


  subroutine blocketteRes(dissApprox, viscApprox)

    ! Copy the values from blockPointers (assumed set) into the
    ! blockette

    use constants
    use blockPointers, only : bnx=>nx, bny=>ny, bnz=>nz, &
         bil=>il, bjl=>jl, bkl=>kl, &
         bie=>ie, bje=>je, bke=>ke, &
         bib=>ib, bjb=>jb, bkb=>kb, &
         bw=>w, bp=>p, bgamma=>gamma, &
         bx=>x, brlv=>rlv, brev=>rev, bvol=>vol, bVolRef=>volRef, bd2wall=>d2wall, &
         biblank=>iblank, bPorI=>porI, bPorJ=>porJ, bPorK=>porK, bdw=>dw, bfw=>fw, &
         bShockSensor=>shockSensor, &
         bsi=>si, bsj=>sj, bsk=>sk
    
    use block
    use surfaceFamilies, only : fullFamList
    use flowVarRefState, only : nwf, nw, viscous
    use communication, only : myid, adflow_comm_world
    implicit none

    ! Input
    logical :: dissApprox, viscApprox
    ! Misc
    integer(kind=intType) :: i, j, k, l, iSize, nn, nIter, ii,jj,kk, iiter, loopCount, ierr
    real(kind=realType) ::  timeA, timeB, tmp, tmp2
    integer(kind=intType) :: omp_get_thread_num

    ! Block loop over the owned cells
    !$OMP parallel do private(i,j,k,l) collapse(2) 
    do kk=2, bkl, BS
       do jj=2, bjl, BS
          do ii=2, bil, BS

             ! Determine the actual size this block will be and set
             ! the sizes in the blockette module for each of the
             ! subroutines.

             nx = min(ii+BS-1, bil) - ii + 1
             ny = min(jj+BS-1, bjl) - jj + 1
             nz = min(kk+BS-1, bkl) - kk + 1

             il = nx + 1; jl = ny + 1; kl = nz + 1
             ie = nx + 2; je = ny + 2; ke = nz + 2
             ib = nx + 3; jb = ny + 3; kb = nz + 3
             
             ! -------------------------------------
             !      Fill in all values
             ! -------------------------------------

             ! Double halos
             do k=0, kb
                do j=0, jb
                   do i=0, ib 
                      w(i,j,k,:) = bw(i+ii-2, j+jj-2, k+kk-2, :)
                      p(i,j,k) = bP(i+ii-2, j+jj-2, k+kk-2)
                      gamma(i,j,k) = bgamma(i+ii-2, j+jj-2, k+kk-2)
                      ss(i,j,k) = bShockSensor(i+ii-2, j+jj-2,k+kk-2)
                   end do
                end do
             end do

             ! Single halos
             do k=1, ke
                do j=1, je
                   do i=1, ie 
                      rlv(i,j,k) = brlv(i+ii-2, j+jj-2, k+kk-2)
                      rev(i,j,k) = brev(i+ii-2, j+jj-2, k+kk-2)
                      vol(i,j,k) = bvol(i+ii-2, j+jj-2, k+kk-2)
                   end do
                end do
             end do

             ! X
             do k=0, ke
                do j=0, je
                   do i=0, ie
                      x(i,j,k,:) = bx(i+ii-2, j+jj-2, k+kk-2, :)
                   end do
                end do
             end do

             ! No Halos (no change)
             do k=2, kl
                do j=2, jl
                   do i=2, il
                      iblank(i,j,k) = biblank(i+ii-2,j+jj-2,k+kk-2)
                      d2wall(i,j,k) = bd2wall(i+ii-2,j+jj-2,k+kk-2)
                      volRef(i,j,k) = bvolRef(i+ii-2,j+jj-2,k+kk-2)
                   end do
                end do
             end do

             ! Porosities (no change)
             do k=2, kl
                do j=2, jl
                   do i=1, il
                      porI(i,j,k) = bporI(i+ii-2,j+jj-2,k+kk-2)
                   end do
                end do
             end do

             do k=2, kl
                do j=1, jl
                   do i=2, il
                      PorJ(i,j,k) = bporJ(i+ii-2,j+jj-2,k+kk-2)
                   end do
                end do
             end do

             do k=1, kl
                do j=2, jl
                   do i=2, il
                      PorK(i,j,k) = bporK(i+ii-2,j+jj-2,k+kk-2)
                   end do
                end do
             end do

             ! fw would normally have a useful value.
             fw = zero

             ! Call the routines in order:
             call metrics
             call initres
             call spectralRadius
             call inviscidCentralFlux
             if (dissApprox) then 
                call inviscidDissFluxScalarApprox
             else
                call inviscidDissFluxScalar
             end if
              if (viscous) then 
                 call computeSpeedOfSoundSquared
                 if (viscApprox) then 
                    call viscousFluxApprox
                 else
                    call allNodalGradients
                    call viscousFlux
                 end if
              end if
             call sumDwandFw

             ! Now we can just set the part of dw we computed
             ! (owned cells only) and we're done!
             do l=1, nwf
                do k=2, kl
                   do j=2, jl
                      do i=2, il
                         bdw(i+ii-2,j+jj-2,k+kk-2,l) = dw(i,j,k,l)
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL DO 
  end subroutine blocketteRes


  subroutine metrics
    ! ---------------------------------------------
    !              Metric computation
    ! ---------------------------------------------

    use constants
    implicit none

    integer(kind=intType) :: i, j, k, l, m, n
    real(kind=realType), dimension(3) :: v1, v2
    real(kind=realType) :: fact

    ! Projected areas of cell faces in the i direction.
    fact = half
    do k=1,ke
       n = k -1
       do j=1,je
          m = j -1
          do i=0,ie

             ! Determine the two diagonal vectors of the face.

             v1(1) = x(i,j,n,1) - x(i,m,k,1)
             v1(2) = x(i,j,n,2) - x(i,m,k,2)
             v1(3) = x(i,j,n,3) - x(i,m,k,3)

             v2(1) = x(i,j,k,1) - x(i,m,n,1)
             v2(2) = x(i,j,k,2) - x(i,m,n,2)
             v2(3) = x(i,j,k,3) - x(i,m,n,3)

             ! The face normal, which is the cross product of the two
             ! diagonal vectors times fact; remember that fact is
             ! either -0.5 or 0.5.

             si(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
             si(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
             si(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

          enddo
       enddo
    enddo

    ! Projected areas of cell faces in the j direction.

    do k=1,ke
       n = k -1
       do j=0,je
          do i=1,ie
             l = i -1

             ! Determine the two diagonal vectors of the face.

             v1(1) = x(i,j,n,1) - x(l,j,k,1)
             v1(2) = x(i,j,n,2) - x(l,j,k,2)
             v1(3) = x(i,j,n,3) - x(l,j,k,3)

             v2(1) = x(l,j,n,1) - x(i,j,k,1)
             v2(2) = x(l,j,n,2) - x(i,j,k,2)
             v2(3) = x(l,j,n,3) - x(i,j,k,3)

             ! The face normal, which is the cross product of the two
             ! diagonal vectors times fact; remember that fact is
             ! either -0.5 or 0.5.

             sj(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
             sj(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
             sj(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

          enddo
       enddo
    enddo

    ! Projected areas of cell faces in the k direction.

    do k=0,ke
       do j=1,je
          m = j -1
          do i=1,ie
             l = i -1
 
             ! Determine the two diagonal vectors of the face.

             v1(1) = x(i,j,k,1) - x(l,m,k,1)
             v1(2) = x(i,j,k,2) - x(l,m,k,2)
             v1(3) = x(i,j,k,3) - x(l,m,k,3)

             v2(1) = x(l,j,k,1) - x(i,m,k,1)
             v2(2) = x(l,j,k,2) - x(i,m,k,2)
             v2(3) = x(l,j,k,3) - x(i,m,k,3)

             ! The face normal, which is the cross product of the two
             ! diagonal vectors times fact; remember that fact is
             ! either -0.5 or 0.5.

             sk(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
             sk(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
             sk(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

          enddo
       enddo
    enddo
  end subroutine metrics

  subroutine initRes
    ! ---------------------------------------------
    !                     Init Res
    ! ---------------------------------------------

    use constants
    implicit none


    ! Obviously this needs to be more complex for the actual code.
    dw = zero

  end subroutine initRes

  subroutine saSource
    ! ---------------------------------------------
    !                    SA Source Term
    ! ---------------------------------------------

    use constants 
    use paramTurb
    use blockPointers, only : sectionID
    use inputPhysics, only :useft2SA, useRotationSA, turbProd, equations
    use section, only : sections
    use sa, only : cv13, kar2Inv, cw36, cb3Inv
    use flowvarRefState, only : timeRef

    implicit none

    ! Variables for sa Souce
    real(kind=realType) :: fv1, fv2, ft2
    real(kind=realType) :: sst, nu, dist2Inv, chi, chi2, chi3
    real(kind=realType) :: rr, gg, gg6, termFw, fwSa, term1, term2
    real(kind=realType) :: dfv1, dfv2, dft2, drr, dgg, dfw, sqrtProd
    real(kind=realType) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
    real(kind=realType) :: div2, fact, sxx, syy, szz, sxy, sxz, syz
    real(kind=realType) :: vortx, vorty, vortz
    real(kind=realType) :: omegax, omegay, omegaz
    real(kind=realType) :: strainMag2, strainProd, vortProd
    real(kind=realType), parameter :: xminn = 1.e-10_realType
    real(kind=realType), parameter :: f23 = two*third
    integer(kind=intType) :: i, j, k

    ! Set model constants
    cv13    = rsaCv1**3
    kar2Inv = one/(rsaK**2)
    cw36    = rsaCw3**6
    cb3Inv  = one/rsaCb3

    ! Determine the non-dimensional wheel speed of this block.

    omegax = timeRef*sections(sectionID)%rotRate(1)
    omegay = timeRef*sections(sectionID)%rotRate(2)
    omegaz = timeRef*sections(sectionID)%rotRate(3)
    do k=2, kl
       do j=2, jl
          do i=2, il

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

                sqrtProd = sqrt(strainProd)

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

                sqrtProd = sqrt(vortProd)

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

             sst = sqrtProd + w(i,j,k,itu1)*fv2*kar2Inv*dist2Inv

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

             term1 = rsaCb1*(one-ft2)*sqrtProd
             term2 = dist2Inv*(kar2Inv*rsaCb1*((one-ft2)*fv2 + ft2) &
                  -           rsaCw1*fwSa)

             dw(i, j, k, itu1) = (term1 + term2*w(i,j,k,itu1))*w(i,j,k,itu1)

          enddo
       enddo
    enddo
  end subroutine saSource

  subroutine saViscous
    ! ---------------------------------------------
    !                    SA Viscous Term
    ! ---------------------------------------------

    use constants
    use sa, only : cv13, kar2Inv, cw36, cb3Inv
    use paramTurb
    implicit none

    ! Variables for sa Viscous
    real(kind=realType) :: voli, volmi, volpi, xm, ym, zm, xp, yp, zp
    real(kind=realType) :: xa, ya, za, ttm, ttp, cnud, cam, cap
    real(kind=realType) :: nutm, nutp, num, nup, cdm, cdp
    real(kind=realType) :: c1m, c1p, c10, b1, c1, d1, qs, nu
    integer(Kind=intType) :: i, j, k

    ! Set model constants
    cv13    = rsaCv1**3
    kar2Inv = one/(rsaK**2)
    cw36    = rsaCw3**6
    cb3Inv  = one/rsaCb3

    !
    !       Viscous terms in k-direction.
    !
    do k=2, kl
       do j=2, jl
          do i=2, il

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

             nutm = half*(w(i,j,k-1,itu1) + w(i,j,k,itu1))
             nutp = half*(w(i,j,k+1,itu1) + w(i,j,k,itu1))
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

             dw(i,j,k,itu1) = dw(i,j,k,itu1)      + c1m*w(i,j,k-1,itu1) &
                  - c10*w(i,j,k,itu1) + c1p*w(i,j,k+1,itu1)
          end do
       enddo
    enddo
    !
    !       Viscous terms in j-direction.
    !
    do k=2, kl
       do j=2, jl
          do i=2, il

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

             dw(i,j,k,itu1) = dw(i,j,k,itu1)      + c1m*w(i,j-1,k,itu1) &
                  - c10*w(i,j,k,itu1) + c1p*w(i,j+1,k,itu1)

          enddo
       enddo
    enddo
    !
    !       Viscous terms in i-direction.
    !
    do k=2, kl
       do j=2, jl
          do i=2, il

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

             dw(i,j,k,itu1) = dw(i,j,k,itu1)      + c1m*w(i-1,j,k,itu1) &
                  - c10*w(i,j,k,itu1) + c1p*w(i+1,j,k,itu1)
          enddo
       enddo
    enddo
  end subroutine saViscous

  subroutine saAdvection
    ! ---------------------------------------------
    !                SA Advection 
    ! ---------------------------------------------
    use constants
    use turbMod, only : secondOrd
    implicit none

    ! Variables for sa Advection
    real(kind=realType) :: uu, dwt, dwtm1, dwtp1, dwti, dwtj, dwtk, qs
    real(kind=realType) :: voli, xa, ya, za
    integer(kind=intType), parameter :: nAdv=1
    integer(kind=intType) :: offset, i, j, k, ii, jj


    offset=itu1-1
    do k=2, kl
       do j=2, jl
          do i=2, il

             ! Compute the grid velocity if present.
             ! It is taken as the average of k and k-1,

             voli = half/vol(i,j,k)

             ! if( addGridVelocities ) &
             !      qs = (sFaceK(i,j,k) + sFaceK(i,j,k-1))*voli

             ! Compute the normal velocity, where the normal direction
             ! is taken as the average of faces k and k-1.

             xa = (sk(i,j,k,1) + sk(i,j,k-1,1))*voli
             ya = (sk(i,j,k,2) + sk(i,j,k-1,2))*voli
             za = (sk(i,j,k,3) + sk(i,j,k-1,3))*voli

             uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs

             ! Determine the situation we are having here, i.e. positive
             ! or negative normal velocity.

             velKdir: if(uu > zero) then

                ! Velocity has a component in positive k-direction.
                ! Loop over the number of advection equations.

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

                   dw(i,j,k,itu1+ii-1) = dw(i,j,k,itu1+ii-1) - uu*dwtk
                enddo

             else velKdir

                ! Velocity has a component in negative k-direction.
                ! Loop over the number of advection equations
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

                   dw(i,j,k,itu1+ii-1) = dw(i,j,k,itu1+ii-1) - uu*dwtk
                end do
             endif velKdir
          enddo
       enddo
    enddo

    !
    !       Upwind discretization of the convective term in j (eta)        
    !       direction. Either the 1st order upwind or the second order     
    !       fully upwind interpolation scheme, kappa = -1, is used in      
    !       combination with the minmod limiter.                           
    !       The possible grid velocity must be taken into account.         
    !
    do k=2, kl
       do j=2, jl
          do i=2, il


             ! Compute the grid velocity if present.
             ! It is taken as the average of j and j-1,

             voli = half/vol(i,j,k)
             ! if( addGridVelocities ) &
             !      qs = (sFaceJ(i,j,k) + sFaceJ(i,j-1,k))*voli

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

                   dw(i,j,k,itu1+ii-1) = dw(i,j,k,itu1+ii-1) - uu*dwtj
                enddo

             else velJdir

                ! Velocity has a component in negative j-direction.
                ! Loop over the number of advection equations.
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

                   dw(i,j,k,itu1+ii-1) = dw(i,j,k,itu1+ii-1) - uu*dwtj
                enddo
             endif velJdir
          enddo
       enddo
    enddo
    !
    !       Upwind discretization of the convective term in i (xi)         
    !       direction. Either the 1st order upwind or the second order     
    !       fully upwind interpolation scheme, kappa = -1, is used in      
    !       combination with the minmod limiter.                           
    !       The possible grid velocity must be taken into account.         
    !
    qs = zero
    do k=2, kl
       do j=2, jl
          do i=2, il
             ! Compute the grid velocity if present.
             ! It is taken as the average of i and i-1,

             voli = half/vol(i,j,k)
             ! if( addGridVelocities ) &
             !         qs = (sFaceI(i,j,k) + sFaceI(i-1,j,k))*voli

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

                   dw(i,j,k,itu1+ii-1) = dw(i,j,k,itu1+ii-1) - uu*dwti
                enddo

             else velIdir

                ! Velocity has a component in negative i-direction.
                ! Loop over the number of advection equations.
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

                   dw(i,j,k,itu1+ii-1) = dw(i,j,k,itu1+ii-1) - uu*dwti

                   ! Update the central jacobian. First the term which is
                   ! always present, i.e. -uu.
                enddo

             endif velIdir
          enddo
       enddo
    enddo
  end subroutine saAdvection

  subroutine spectralRadius
    ! ---------------------------------------------
    !               Spectral Radius
    ! ---------------------------------------------

    use constants
    use flowvarRefState, only : pInfCorr, rhoInf, gammaInf
    use inputDiscretization, only : adis

    implicit none

    ! Variables for spectral Radius
    real(kind=realType) :: plim, rlim, clim2
    real(kind=realType) :: cc2, qsi, qsj, qsk, sx, sy, sz, rmu
    real(kind=realType) :: ri, rj, rk, rij, rjk, rki
    real(kind=realType) :: vsi, vsj, vsk, rfl, dpi, dpj, dpk
    real(kind=realType) :: sFace, tmp, uux, uuy, uuz
    logical :: doScaling
    integer(kind=intType) :: i, j, k

    ! Set the value of plim. To be fully consistent this must have
    ! the dimension of a pressure. Therefore a fraction of pInfCorr
    ! is used. Idem for rlim; compute clim2 as well.

    plim  = 0.001_realType*pInfCorr
    rlim  = 0.001_realType*rhoInf
    clim2 = 0.000001_realType*gammaInf*pInfCorr/rhoInf
    doScaling = .True.

    ! Initialize sFace to zero. This value will be used if the
    ! block is not moving.

    sFace = zero
    !
    !           Inviscid contribution, depending on the preconditioner.    
    !           Compute the cell centered values of the spectral radii.    
    !
    do k=1,ke
       do j=1,je
          do i=1,ie

             ! Compute the velocities and speed of sound squared.

             uux  = w(i,j,k,ivx)
             uuy  = w(i,j,k,ivy)
             uuz  = w(i,j,k,ivz)
             cc2 = gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho)
             cc2 = max(cc2,clim2)

             ! Set the dot product of the grid velocity and the
             ! normal in i-direction for a moving face. To avoid
             ! a number of multiplications by 0.5 simply the sum
             ! is taken.

             ! if( addGridVelocities ) &
             !      sFace = sFaceI(i-1,j,k) + sFaceI(i,j,k)

             ! Spectral radius in i-direction.

             sx = si(i-1,j,k,1) + si(i,j,k,1)
             sy = si(i-1,j,k,2) + si(i,j,k,2)
             sz = si(i-1,j,k,3) + si(i,j,k,3)

             qsi = uux*sx + uuy*sy + uuz*sz - sFace

             ri = half*(abs(qsi) &
                  +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

             ! The grid velocity in j-direction.

             ! if( addGridVelocities ) &
             !      sFace = sFaceJ(i,j-1,k) + sFaceJ(i,j,k)

             ! Spectral radius in j-direction.

             sx = sj(i,j-1,k,1) + sj(i,j,k,1)
             sy = sj(i,j-1,k,2) + sj(i,j,k,2)
             sz = sj(i,j-1,k,3) + sj(i,j,k,3)

             qsj = uux*sx + uuy*sy + uuz*sz - sFace

             rj = half*(abs(qsj) &
                  +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

             ! The grid velocity in k-direction.

             ! if( addGridVelocities ) &
             !      sFace = sFaceK(i,j,k-1) + sFaceK(i,j,k)

             ! Spectral radius in k-direction.

             sx = sk(i,j,k-1,1) + sk(i,j,k,1)
             sy = sk(i,j,k-1,2) + sk(i,j,k,2)
             sz = sk(i,j,k-1,3) + sk(i,j,k,3)

             qsk = uux*sx + uuy*sy + uuz*sz - sFace

             rk = half*(abs(qsk) &
                  +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

             !
             !           Adapt the spectral radii if directional scaling must be    
             !           applied.                                                   
             !

             ! Avoid division by zero by clipping radi, radJ and
             ! radK.

             ri = max(ri, eps)
             rj = max(rj, eps)
             rk = max(rk, eps)

             ! Compute the scaling in the three coordinate
             ! directions.

             rij = (ri/rj)**adis
             rjk = (rj/rk)**adis
             rki = (rk/ri)**adis

             ! Create the scaled versions of the aspect ratios.
             ! Note that the multiplication is done with radi, radJ
             ! and radK, such that the influence of the clipping
             ! is negligible.

             radi(i,j,k) = ri*(one + one/rij + rki)
             radJ(i,j,k) = rj*(one + one/rjk + rij)
             radK(i,j,k) = rk*(one + one/rki + rjk)
          end do
       enddo
    enddo
  end subroutine spectralRadius

  subroutine inviscidCentralFlux

    ! ---------------------------------------------
    !               Inviscid central flux
    ! ---------------------------------------------
    use constants
    implicit none

    ! Variables for inviscid central flux
    real(kind=realType) :: qsp, qsm, rqsp, rqsm, porVel, porFlux
    real(kind=realType) :: pa, vnp, vnm, fs, sFace
    integer(kind=intType) :: i, j, k

    ! Initialize sFace to zero. This value will be used if the
    ! block is not moving.
    sFace = zero

    do k=2, kl
       do j=2, jl
          do i=1, il

             ! Set the dot product of the grid velocity and the
             ! normal in i-direction for a moving face.

             !if( addGridVelocities ) sFace = sFaceI(i,j,k)

             ! Compute the normal velocities of the left and right state.

             vnp = w(i+1,j,k,ivx)*sI(i,j,k,1) &
                  + w(i+1,j,k,ivy)*sI(i,j,k,2) &
                  + w(i+1,j,k,ivz)*sI(i,j,k,3)
             vnm = w(i,  j,k,ivx)*sI(i,j,k,1) &
                  + w(i,  j,k,ivy)*sI(i,j,k,2) &
                  + w(i,  j,k,ivz)*sI(i,j,k,3)
             ! Set the values of the porosities for this face.
             ! porVel defines the porosity w.r.t. velocity;
             ! porFlux defines the porosity w.r.t. the entire flux.
             ! The latter is only zero for a discontinuous block
             ! boundary that must be treated conservatively.
             ! The default value of porFlux is 0.5, such that the
             ! correct central flux is scattered to both cells.
             ! In case of a boundFlux the normal velocity is set
             ! to sFace.

             porVel  = one
             porFlux = half
             if(porI(i,j,k) == noFlux)    porFlux = zero
             if(porI(i,j,k) == boundFlux) then
                porVel = zero
                vnp    = sFace
                vnm    = sFace
             endif

             ! Incorporate porFlux in porVel.

             porVel = porVel*porFlux

             ! Compute the normal velocities relative to the grid for
             ! the face as well as the mass fluxes.

             qsp = (vnp -sFace)*porVel
             qsm = (vnm -sFace)*porVel

             rqsp = qsp*w(i+1,j,k,irho)
             rqsm = qsm*w(i,  j,k,irho)

             ! Compute the sum of the pressure multiplied by porFlux.
             ! For the default value of porFlux, 0.5, this leads to
             ! the average pressure.

             pa = porFlux*(p(i+1,j,k) + p(i,j,k))

             ! Compute the fluxes and scatter them to the cells
             ! i,j,k and i+1,j,k. Store the density flux in the
             ! mass flow of the appropriate sliding mesh interface.

             fs = rqsp + rqsm
             dw(i+1,j,k,irho) = dw(i+1,j,k,irho) - fs
             dw(i,  j,k,irho) = dw(i,  j,k,irho) + fs

             fs = rqsp*w(i+1,j,k,ivx) + rqsm*w(i,j,k,ivx) &
                  + pa*sI(i,j,k,1)
             dw(i+1,j,k,imx) = dw(i+1,j,k,imx) - fs
             dw(i,  j,k,imx) = dw(i,  j,k,imx) + fs

             fs = rqsp*w(i+1,j,k,ivy) + rqsm*w(i,j,k,ivy) &
                  + pa*sI(i,j,k,2)
             dw(i+1,j,k,imy) = dw(i+1,j,k,imy) - fs
             dw(i,  j,k,imy) = dw(i,  j,k,imy) + fs

             fs = rqsp*w(i+1,j,k,ivz) + rqsm*w(i,j,k,ivz) &
                  + pa*sI(i,j,k,3)
             dw(i+1,j,k,imz) = dw(i+1,j,k,imz) - fs
             dw(i,  j,k,imz) = dw(i,  j,k,imz) + fs

             fs = qsp*w(i+1,j,k,irhoE) + qsm*w(i,j,k,irhoE) &
                  + porFlux*(vnp*p(i+1,j,k) + vnm*p(i,j,k))
             dw(i+1,j,k,irhoE) = dw(i+1,j,k,irhoE) - fs
             dw(i,  j,k,irhoE) = dw(i,  j,k,irhoE) + fs
          enddo
       enddo
    enddo

    sface = zero
    do k=2,kl
       do j=1,jl
          do i=2,il

             ! Set the dot product of the grid velocity and the
             ! normal in j-direction for a moving face.

             !if( addGridVelocities ) sFace = sFaceJ(i,j,k)

             ! Compute the normal velocities of the left and right state.

             vnp = w(i,j+1,k,ivx)*sJ(i,j,k,1) &
                  + w(i,j+1,k,ivy)*sJ(i,j,k,2) &
                  + w(i,j+1,k,ivz)*sJ(i,j,k,3)
             vnm = w(i,j,  k,ivx)*sJ(i,j,k,1) &
                  + w(i,j,  k,ivy)*sJ(i,j,k,2) &
                  + w(i,j,  k,ivz)*sJ(i,j,k,3)

             ! Set the values of the porosities for this face.
             ! porVel defines the porosity w.r.t. velocity;
             ! porFlux defines the porosity w.r.t. the entire flux.
             ! The latter is only zero for a discontinuous block
             ! boundary that must be treated conservatively.
             ! The default value of porFlux is 0.5, such that the
             ! correct central flux is scattered to both cells.
             ! In case of a boundFlux the normal velocity is set
             ! to sFace.

             porVel  = one
             porFlux = half
             if(porJ(i,j,k) == noFlux)    porFlux = zero
             if(porJ(i,j,k) == boundFlux) then
                porVel = zero
                vnp    = sFace
                vnm    = sFace
             endif

             ! Incorporate porFlux in porVel.

             porVel = porVel*porFlux

             ! Compute the normal velocities for the face as well as the
             ! mass fluxes.

             qsp = (vnp - sFace)*porVel
             qsm = (vnm - sFace)*porVel

             rqsp = qsp*w(i,j+1,k,irho)
             rqsm = qsm*w(i,j,  k,irho)

             ! Compute the sum of the pressure multiplied by porFlux.
             ! For the default value of porFlux, 0.5, this leads to
             ! the average pressure.

             pa = porFlux*(p(i,j+1,k) + p(i,j,k))

             ! Compute the fluxes and scatter them to the cells
             ! i,j,k and i,j+1,k. Store the density flux in the
             ! mass flow of the appropriate sliding mesh interface.

             fs = rqsp + rqsm
             dw(i,j+1,k,irho) = dw(i,j+1,k,irho) - fs
             dw(i,j,  k,irho) = dw(i,j,  k,irho) + fs

             fs = rqsp*w(i,j+1,k,ivx) + rqsm*w(i,j,k,ivx) &
                  + pa*sJ(i,j,k,1)
             dw(i,j+1,k,imx) = dw(i,j+1,k,imx) - fs
             dw(i,j,  k,imx) = dw(i,j,  k,imx) + fs

             fs = rqsp*w(i,j+1,k,ivy) + rqsm*w(i,j,k,ivy) &
                  + pa*sJ(i,j,k,2)
             dw(i,j+1,k,imy) = dw(i,j+1,k,imy) - fs
             dw(i,j,  k,imy) = dw(i,j,  k,imy) + fs

             fs = rqsp*w(i,j+1,k,ivz) + rqsm*w(i,j,k,ivz) &
                  + pa*sJ(i,j,k,3)
             dw(i,j+1,k,imz) = dw(i,j+1,k,imz) - fs
             dw(i,j,  k,imz) = dw(i,j,  k,imz) + fs

             fs = qsp*w(i,j+1,k,irhoE) + qsm*w(i,j,k,irhoE) &
                  + porFlux*(vnp*p(i,j+1,k) + vnm*p(i,j,k))
             dw(i,j+1,k,irhoE) = dw(i,j+1,k,irhoE) - fs
             dw(i,j,  k,irhoE) = dw(i,j,  k,irhoE) + fs
          enddo
       enddo
    enddo

    sface = zero

      do k=1,kl
          do j=2,jl
             do i=2,il

                ! Set the dot product of the grid velocity and the
                ! normal in k-direction for a moving face.
                
                !if( addGridVelocities ) sFace = sFaceK(i,j,k)

                ! Compute the normal velocities of the left and right state.

                vnp = w(i,j,k+1,ivx)*sK(i,j,k,1) &
                     + w(i,j,k+1,ivy)*sK(i,j,k,2) &
                     + w(i,j,k+1,ivz)*sK(i,j,k,3)
                vnm = w(i,j,k,  ivx)*sK(i,j,k,1) &
                     + w(i,j,k,  ivy)*sK(i,j,k,2) &
                     + w(i,j,k,  ivz)*sK(i,j,k,3)

                ! Set the values of the porosities for this face.
                ! porVel defines the porosity w.r.t. velocity;
                ! porFlux defines the porosity w.r.t. the entire flux.
                ! The latter is only zero for a discontinuous block
                ! block boundary that must be treated conservatively.
                ! The default value of porFlux is 0.5, such that the
                ! correct central flux is scattered to both cells.
                ! In case of a boundFlux the normal velocity is set
                ! to sFace.

                porVel  = one
                porFlux = half

                if(porK(i,j,k) == noFlux)    porFlux = zero
                if(porK(i,j,k) == boundFlux) then
                   porVel = zero
                   vnp    = sFace
                   vnm    = sFace
                endif

                ! Incorporate porFlux in porVel.

                porVel = porVel*porFlux

                ! Compute the normal velocities for the face as well as the
                ! mass fluxes.

                qsp = (vnp - sFace)*porVel
                qsm = (vnm - sFace)*porVel

                rqsp = qsp*w(i,j,k+1,irho)
                rqsm = qsm*w(i,j,k,  irho)

                ! Compute the sum of the pressure multiplied by porFlux.
                ! For the default value of porFlux, 0.5, this leads to
                ! the average pressure.

                pa = porFlux*(p(i,j,k+1) + p(i,j,k))

                ! Compute the fluxes and scatter them to the cells
                ! i,j,k and i,j,k+1. Store the density flux in the
                ! mass flow of the appropriate sliding mesh interface.

                fs = rqsp + rqsm
                dw(i,j,k+1,irho) = dw(i,j,k+1,irho) - fs
                dw(i,j,k,  irho) = dw(i,j,k,  irho) + fs

                fs = rqsp*w(i,j,k+1,ivx) + rqsm*w(i,j,k,ivx) &
                     + pa*sK(i,j,k,1)
                dw(i,j,k+1,imx) = dw(i,j,k+1,imx) - fs
                dw(i,j,k,  imx) = dw(i,j,k,  imx) + fs

                fs = rqsp*w(i,j,k+1,ivy) + rqsm*w(i,j,k,ivy) &
                     + pa*sK(i,j,k,2)
                dw(i,j,k+1,imy) = dw(i,j,k+1,imy) - fs
                dw(i,j,k,  imy) = dw(i,j,k,  imy) + fs

                fs = rqsp*w(i,j,k+1,ivz) + rqsm*w(i,j,k,ivz) &
                     + pa*sK(i,j,k,3)
                dw(i,j,k+1,imz) = dw(i,j,k+1,imz) - fs
                dw(i,j,k,  imz) = dw(i,j,k,  imz) + fs

                fs = qsp*w(i,j,k+1,irhoE) + qsm*w(i,j,k,irhoE) &
                     + porFlux*(vnp*p(i,j,k+1) + vnm*p(i,j,k))
                dw(i,j,k+1,irhoE) = dw(i,j,k+1,irhoE) - fs
                dw(i,j,k,  irhoE) = dw(i,j,k,  irhoE) + fs

          enddo
       enddo
    enddo
  end subroutine inviscidCentralFlux

  subroutine inviscidDissFluxScalar
    ! ---------------------------------------------
    !             Inviscid Diss Flux Scalar
    ! ---------------------------------------------

    use constants 
    use flowVarRefState, only : pInfCorr
    use inputDiscretization, only: vis2, vis4
    use inputPhysics, only : equations
    use iteration, only : rFil
    use flowVarRefState, only : gammaInf, pInfCorr, rhoInf
    implicit none

    ! Variables for inviscid diss flux scalar
    real(kind=realType), parameter :: dssMax = 0.25_realType
    real(kind=realType) :: sslim, rhoi
    real(kind=realType) :: sfil, fis2, fis4
    real(kind=realType) :: ppor, rrad, dis2, dis4, fs
    real(kind=realType) :: ddw1,ddw2,ddw3,ddw4,ddw5
    integer(kind=intType) :: i, j, k

    ! Determine the variables used to compute the switch.
    ! For the inviscid case this is the pressure; for the viscous
    ! case it is the entropy.

    select case (equations)
    case (EulerEquations)

       ! Inviscid case. Pressure switch is based on the pressure.
       ! Also set the value of sslim. To be fully consistent this
       ! must have the dimension of pressure and it is therefore
       ! set to a fraction of the free stream value.

       sslim = 0.001_realType*pInfCorr

       ! Copy the pressure in ss. Only need the entries used in the
       ! discretization, i.e. not including the corner halo's, but we'll
       ! just copy all anyway. 

       ss = P
       !===============================================================

    case (NSEquations, RANSEquations)

       ! Viscous case. Pressure switch is based on the entropy.
       ! Also set the value of sslim. To be fully consistent this
       ! must have the dimension of entropy and it is therefore
       ! set to a fraction of the free stream value.

       sslim = 0.001_realType*pInfCorr/(rhoInf**gammaInf)

       ! Store the entropy in ss. See above. 
       do k=0, kb
          do j=0, jb
             do i=0, ib
                ss(i,j,k) = p(i,j,k)/(w(i,j,k,irho)**gamma(i,j,k))
             end do
          end do
       end do
    end select

    ! Compute the pressure sensor for each cell, in each direction:
    do k=1,ke
       do j=1,je
          do i=1,ie
             dss(i,j,k,1) = abs((ss(i+1,j,k) - two*ss(i,j,k) + ss(i-1,j,k)) &
                  /     (ss(i+1,j,k) + two*ss(i,j,k) + ss(i-1,j,k) + sslim))

             dss(i,j,k,2) = abs((ss(i,j+1,k) - two*ss(i,j,k) + ss(i,j-1,k)) &
                  /     (ss(i,j+1,k) + two*ss(i,j,k) + ss(i,j-1,k) + sslim))

             dss(i,j,k,3) = abs((ss(i,j,k+1) - two*ss(i,j,k) + ss(i,j,k-1)) &
                  /     (ss(i,j,k+1) + two*ss(i,j,k) + ss(i,j,k-1) + sslim))
          end do
       end do
    end do

    ! Set a couple of constants for the scheme.

    fis2 = rFil*vis2
    fis4 = rFil*vis4
    sfil = one - rFil

    ! Initialize the dissipative residual to a certain times,
    ! possibly zero, the previously stored value. Owned cells
    ! only, because the halo values do not matter.

    fw = sfil*fw
    !
    !       Dissipative fluxes in the i-direction.                         
    !
    do k=2,kl
       do j=2,jl
          do i=1,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porI(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radI(i,j,k) + radI(i+1,j,k))

             dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,1), dss(i+1,j,k,1)))
             dis4 = dim(fis4*rrad, dis2)

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw1 = w(i+1,j,k,irho) - w(i,j,k,irho)
             fs  = dis2*ddw1 &
                  - dis4*(w(i+2,j,k,irho) - w(i-1,j,k,irho) - three*ddw1)

             fw(i+1,j,k,irho) = fw(i+1,j,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw2 = w(i+1,j,k,ivx)*w(i+1,j,k,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
             fs  = dis2*ddw2 &
                  - dis4*(w(i+2,j,k,ivx)*w(i+2,j,k,irho) - w(i-1,j,k,ivx)*w(i-1,j,k,irho) - three*ddw2)

             fw(i+1,j,k,imx) = fw(i+1,j,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw3 = w(i+1,j,k,ivy)*w(i+1,j,k,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
             fs  = dis2*ddw3 &
                  - dis4*(w(i+2,j,k,ivy)*w(i+2,j,k,irho) - w(i-1,j,k,ivy)*w(i-1,j,k,irho) - three*ddw3)

             fw(i+1,j,k,imy) = fw(i+1,j,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw4 = w(i+1,j,k,ivz)*w(i+1,j,k,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
             fs  = dis2*ddw4 &
                  - dis4*(w(i+2,j,k,ivz)*w(i+2,j,k,irho) - w(i-1,j,k,ivz)*w(i-1,j,k,irho) - three*ddw4)

             fw(i+1,j,k,imz) = fw(i+1,j,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             ddw5 = (w(i+1,j,k,irhoE) + P(i+1,j,K))- (w(i,j,k,irhoE) + P(i,j,k))
             fs  = dis2*ddw5 &
                  - dis4*((w(i+2,j,k,irhoE) + P(i+2,j,k)) - (w(i-1,j,k,irhoE) + P(i-1,j,k)) - three*ddw5)

             fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
          end do
       end do
    end do
    !
    !       Dissipative fluxes in the j-direction.                         
    !
    do k=2,kl
       do j=1,jl
          do i=2,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porJ(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radJ(i,j,k) + radJ(i,j+1,k))

             dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,2),dss(i,j+1,k,2)))
             dis4 = dim(fis4*rrad, dis2)

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw1 = w(i,j+1,k,irho) - w(i,j,k,irho)
             fs  = dis2*ddw1 &
                  - dis4*(w(i,j+2,k,irho) - w(i,j-1,k,irho) - three*ddw1)

             fw(i,j+1,k,irho) = fw(i,j+1,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw2 = w(i,j+1,k,ivx)*w(i,j+1,k,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
             fs  = dis2*ddw2 &
                  - dis4*(w(i,j+2,k,ivx)*w(i,j+2,k,irho) - w(i,j-1,k,ivx)*w(i,j-1,k,irho) - three*ddw2)

             fw(i,j+1,k,imx) = fw(i,j+1,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw3 = w(i,j+1,k,ivy)*w(i,j+1,k,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
             fs  = dis2*ddw3 &
                  - dis4*(w(i,j+2,k,ivy)*w(i,j+2,k,irho) - w(i,j-1,k,ivy)*w(i,j-1,k,irho) - three*ddw3)

             fw(i,j+1,k,imy) = fw(i,j+1,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw4 = w(i,j+1,k,ivz)*w(i,j+1,k,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
             fs  = dis2*ddw4 &
                  - dis4*(w(i,j+2,k,ivz)*w(i,j+2,k,irho) - w(i,j-1,k,ivz)*w(i,j-1,k,irho) - three*ddw4)

             fw(i,j+1,k,imz) = fw(i,j+1,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             ddw5 = (w(i,j+1,k,irhoE) + P(i,j+1,k)) - (w(i,j,k,irhoE) + P(i,j,k))
             fs  = dis2*ddw5 &
                  - dis4*((w(i,j+2,k,irhoE) + P(i,j+2,k)) - (w(i,j-1,k,irhoE) + P(i,j-1,k)) - three*ddw5)

             fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
          end do
       end do
    end do
    !
    !       Dissipative fluxes in the k-direction.                         
    !
    do k=1,kl
       do j=2,jl
          do i=2,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porK(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radK(i,j,k) + radK(i,j,k+1))

             dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,3), dss(i,j,k+1,3)))
             dis4 = dim(fis4*rrad, dis2)

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw1 = w(i,j,k+1,irho) - w(i,j,k,irho)
             fs  = dis2*ddw1 &
                  - dis4*(w(i,j,k+2,irho) - w(i,j,k-1,irho) - three*ddw1)

             fw(i,j,k+1,irho) = fw(i,j,k+1,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw2 = w(i,j,k+1,ivx)*w(i,j,k+1,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
             fs  = dis2*ddw2 &
                  - dis4*(w(i,j,k+2,ivx)*w(i,j,k+2,irho) - w(i,j,k-1,ivx)*w(i,j,k-1,irho) - three*ddw2)

             fw(i,j,k+1,imx) = fw(i,j,k+1,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw3 = w(i,j,k+1,ivy)*w(i,j,k+1,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
             fs  = dis2*ddw3 &
                  - dis4*(w(i,j,k+2,ivy)*w(i,j,k+2,irho) - w(i,j,k-1,ivy)*w(i,j,k-1,irho) - three*ddw3)

             fw(i,j,k+1,imy) = fw(i,j,k+1,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw4 = w(i,j,k+1,ivz)*w(i,j,k+1,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
             fs  = dis2*ddw4 &
                  - dis4*(w(i,j,k+2,ivz)*w(i,j,k+2,irho) - w(i,j,k-1,ivz)*w(i,j,k-1,irho) - three*ddw4)

             fw(i,j,k+1,imz) = fw(i,j,k+1,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             ddw5 = (w(i,j,k+1,irhoE) + P(i,j,k+1)) - (w(i,j,k,irhoE) + P(i,j,k))
             fs  = dis2*ddw5 &
                  - dis4*((w(i,j,k+2,irhoE) + P(i,j,k+2)) - (w(i,j,k-1,irhoE) + P(i,j,k-1)) - three*ddw5)

             fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
          end do
       end do
    end do
  end subroutine inviscidDissFluxScalar

  subroutine inviscidDissFluxScalarApprox
    ! ---------------------------------------------
    !             Inviscid Diss Flux Scalar
    ! ---------------------------------------------

    use constants 
    use flowVarRefState, only : pInfCorr
    use inputDiscretization, only: vis2, vis4, sigma
    use inputPhysics, only : equations
    use iteration, only : rFil
    use flowVarRefState, only : gammaInf, pInfCorr, rhoInf
    implicit none

    ! Variables for inviscid diss flux scalar
    real(kind=realType), parameter :: dssMax = 0.25_realType
    real(kind=realType) :: sslim, rhoi
    real(kind=realType) :: sfil, fis2, fis4
    real(kind=realType) :: ppor, rrad, dis2, dis4, fs
    real(kind=realType) :: ddw
    integer(kind=intType) :: i, j, k
    select case (equations)
    case (EulerEquations)

       ! Inviscid case. Pressure switch is based on the pressure.
       ! Also set the value of sslim. To be fully consistent this
       ! must have the dimension of pressure and it is therefore
       ! set to a fraction of the free stream value.

       sslim = 0.001_realType*pInfCorr

       !===============================================================

    case (NSEquations, RANSEquations)

       ! Viscous case. Pressure switch is based on the entropy.
       ! Also set the value of sslim. To be fully consistent this
       ! must have the dimension of entropy and it is therefore
       ! set to a fraction of the free stream value.

       sslim = 0.001_realType*pInfCorr/(rhoInf**gammaInf)

    end select

    ! Compute the pressure sensor for each cell, in each direction:
    do k=1,ke
       do j=1,je
          do i=1,ie
             dss(i,j,k,1) = abs((ss(i+1,j,k) - two*ss(i,j,k) + ss(i-1,j,k)) &
                  /     (ss(i+1,j,k) + two*ss(i,j,k) + ss(i-1,j,k) + sslim))

             dss(i,j,k,2) = abs((ss(i,j+1,k) - two*ss(i,j,k) + ss(i,j-1,k)) &
                  /     (ss(i,j+1,k) + two*ss(i,j,k) + ss(i,j-1,k) + sslim))

             dss(i,j,k,3) = abs((ss(i,j,k+1) - two*ss(i,j,k) + ss(i,j,k-1)) &
                  /     (ss(i,j,k+1) + two*ss(i,j,k) + ss(i,j,k-1) + sslim))
          end do
       end do
    end do

    ! Set a couple of constants for the scheme.
    fis2 = vis2
    fis4 = vis4
    !
    !       Dissipative fluxes in the i-direction.                         
    !
    do k=2,kl
       do j=2,jl
          do i=1,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porI(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radI(i,j,k) + radI(i+1,j,k))

             dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,1), dss(i+1,j,k,1))) + sigma*fis4*rrad

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = w(i+1,j,k,irho) - w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i+1,j,k,irho) = fw(i+1,j,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw = w(i+1,j,k,ivx)*w(i+1,j,k,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i+1,j,k,imx) = fw(i+1,j,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw = w(i+1,j,k,ivy)*w(i+1,j,k,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i+1,j,k,imy) = fw(i+1,j,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.
             ddw = w(i+1,j,k,ivz)*w(i+1,j,k,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i+1,j,k,imz) = fw(i+1,j,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.
             ddw = (w(i+1,j,k,irhoE) + P(i+1,j,K))- (w(i,j,k,irhoE) + P(i,j,k))
             fs  = dis2*ddw

             fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs

          end do
       end do
    end do
    !
    !       Dissipative fluxes in the j-direction.                         
    !
    do k=2,kl
       do j=1,jl
          do i=2,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porJ(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radJ(i,j,k) + radJ(i,j+1,k))

             dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,2),dss(i,j+1,k,2))) +sigma*fis4*rrad

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = w(i,j+1,k,irho) - w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j+1,k,irho) = fw(i,j+1,k,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw = w(i,j+1,k,ivx)*w(i,j+1,k,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j+1,k,imx) = fw(i,j+1,k,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw = w(i,j+1,k,ivy)*w(i,j+1,k,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j+1,k,imy) = fw(i,j+1,k,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw = w(i,j+1,k,ivz)*w(i,j+1,k,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j+1,k,imz) = fw(i,j+1,k,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.

             ddw = (w(i,j+1,k,irhoE) + P(i,j+1,k)) - (w(i,j,k,irhoE) + P(i,j,k))
             fs  = dis2*ddw

             fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
          end do
       end do
    end do
    !
    !       Dissipative fluxes in the k-direction.                         
    !
    do k=1,kl
       do j=2,jl
          do i=2,il

             ! Compute the dissipation coefficients for this face.

             ppor = zero
             if(porK(i,j,k) == normalFlux) ppor = half
             rrad = ppor*(radK(i,j,k) + radK(i,j,k+1))

             dis2 = fis2*rrad*min(dssMax, max(dss(i,j,k,3), dss(i,j,k+1,3))) + sigma*fis4*rrad

             ! Compute and scatter the dissipative flux.
             ! Density. Store it in the mass flow of the
             ! appropriate sliding mesh interface.

             ddw = w(i,j,k+1,irho) - w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j,k+1,irho) = fw(i,j,k+1,irho) + fs
             fw(i,j,k,irho)   = fw(i,j,k,irho)   - fs

             ! X-momentum.

             ddw = w(i,j,k+1,ivx)*w(i,j,k+1,irho) - w(i,j,k,ivx)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j,k+1,imx) = fw(i,j,k+1,imx) + fs
             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fs

             ! Y-momentum.

             ddw = w(i,j,k+1,ivy)*w(i,j,k+1,irho) - w(i,j,k,ivy)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j,k+1,imy) = fw(i,j,k+1,imy) + fs
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fs

             ! Z-momentum.

             ddw = w(i,j,k+1,ivz)*w(i,j,k+1,irho) - w(i,j,k,ivz)*w(i,j,k,irho)
             fs  = dis2*ddw

             fw(i,j,k+1,imz) = fw(i,j,k+1,imz) + fs
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fs

             ! Energy.
             ddw = (w(i,j,k+1,irhoE) + P(i,j,k+1)) - (w(i,j,k,irhoE) + P(i,j,k))
             fs  = dis2*ddw

             fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + fs
             fw(i,j,k,irhoE)   = fw(i,j,k,irhoE)   - fs
          end do
       end do
    end do
  end subroutine inviscidDissFluxScalarApprox

  subroutine computeSpeedOfSoundSquared
    ! ---------------------------------------------
    !           Compute the speed of sound squared
    ! ---------------------------------------------
    use constants
    use utils, only : getCorrectForK
    implicit none

    ! Variables for speed of sound
    logical :: correctForK
    real(kind=realType) :: pp
    real(kind=realType), parameter :: twoThird = two*third
    integer(kind=intType) :: i, j, k

    ! Determine if we need to correct for K
    correctForK = getCorrectForK()

    if (correctForK) then 
       do k=1,ke
          do j=1,je
             do i=1,ie
                pp = p(i,j,k) - twoThird*w(i,j,k,irho)*w(i,j,k,itu1)
                aa(i,j,k) = gamma(i,j,k)*pp/w(i,j,k,irho)
             enddo
          enddo
       enddo
    else
       do k=1,ke
          do j=1,je
             do i=1,ie
                aa(i,j,k) = gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho)
             enddo
          enddo
       enddo
    end if
  end subroutine computeSpeedOfSoundSquared

  subroutine allNodalGradients

    ! ---------------------------------------------
    !           Compute nodal gradients
    ! ---------------------------------------------
    use constants
    implicit none

    ! Variables for nodal gradients
    real(kind=realType) :: a2, oVol, uBar, vBar, wBar, sx, sy, sz
    integer(kind=intType) :: i, j, k

    ! Zero just the required part of the nodal gradients since the
    ! first value may be useful.
    ux(:, :, :) = zero
    uy(:, :, :) = zero 
    uz(:, :, :) = zero 

    vx(:, :, :) = zero
    vy(:, :, :) = zero 
    vz(:, :, :) = zero 

    wx(:, :, :) = zero
    wy(:, :, :) = zero 
    wz(:, :, :) = zero 

    qx(:, :, :) = zero
    qy(:, :, :) = zero 
    qz(:, :, :) = zero 

    ! First part. Contribution in the k-direction.
    ! The contribution is scattered to both the left and right node
    ! in k-direction.

    do k=1, ke
       do j=1, jl
          do i=1, il

             ! Compute 8 times the average normal for this part of
             ! the control volume. The factor 8 is taken care of later
             ! on when the division by the volume takes place.

             sx = sk(i,j,k-1,  1) + sk(i+1,j,k-1,  1) &
                  + sk(i,j+1,k-1,1) + sk(i+1,j+1,k-1,1) &
                  + sk(i,j,  k,  1) + sk(i+1,j,  k,  1) &
                  + sk(i,j+1,k  ,1) + sk(i+1,j+1,k  ,1)
             sy = sk(i,j,k-1,  2) + sk(i+1,j,k-1,  2) &
                  + sk(i,j+1,k-1,2) + sk(i+1,j+1,k-1,2) &
                  + sk(i,j,  k,  2) + sk(i+1,j,  k,  2) &
                  + sk(i,j+1,k  ,2) + sk(i+1,j+1,k  ,2)
             sz = sk(i,j,k-1,  3) + sk(i+1,j,k-1,  3) &
                  + sk(i,j+1,k-1,3) + sk(i+1,j+1,k-1,3) &
                  + sk(i,j,  k,  3) + sk(i+1,j,  k,  3) &
                  + sk(i,j+1,k  ,3) + sk(i+1,j+1,k  ,3)

             ! Compute the average velocities and speed of sound squared
             ! for this integration point. Node that these variables are
             ! stored in w(ivx), w(ivy), w(ivz) and p.

             ubar = fourth*(w(i,j,  k,ivx) + w(i+1,j,  k,ivx) &
                  +         w(i,j+1,k,ivx) + w(i+1,j+1,k,ivx))
             vbar = fourth*(w(i,j,  k,ivy) + w(i+1,j,  k,ivy) &
                  +         w(i,j+1,k,ivy) + w(i+1,j+1,k,ivy))
             wbar = fourth*(w(i,j,  k,ivz) + w(i+1,j,  k,ivz) &
                  +         w(i,j+1,k,ivz) + w(i+1,j+1,k,ivz))

             a2 = fourth*(aa(i,j,k) + aa(i+1,j,k) + aa(i,j+1,k) + aa(i+1,j+1,k))

             ! Add the contributions to the surface integral to the node
             ! j-1 and substract it from the node j. For the heat flux it
             ! is reversed, because the negative of the gradient of the
             ! speed of sound must be computed.

             if(k > 1) then
                ux(i,j,k-1) = ux(i,j,k-1) + ubar*sx
                uy(i,j,k-1) = uy(i,j,k-1) + ubar*sy
                uz(i,j,k-1) = uz(i,j,k-1) + ubar*sz

                vx(i,j,k-1) = vx(i,j,k-1) + vbar*sx
                vy(i,j,k-1) = vy(i,j,k-1) + vbar*sy
                vz(i,j,k-1) = vz(i,j,k-1) + vbar*sz

                wx(i,j,k-1) = wx(i,j,k-1) + wbar*sx
                wy(i,j,k-1) = wy(i,j,k-1) + wbar*sy
                wz(i,j,k-1) = wz(i,j,k-1) + wbar*sz

                qx(i,j,k-1) = qx(i,j,k-1) - a2*sx
                qy(i,j,k-1) = qy(i,j,k-1) - a2*sy
                qz(i,j,k-1) = qz(i,j,k-1) - a2*sz
             endif

             if(k < ke) then
                ux(i,j,k) = ux(i,j,k) - ubar*sx
                uy(i,j,k) = uy(i,j,k) - ubar*sy
                uz(i,j,k) = uz(i,j,k) - ubar*sz

                vx(i,j,k) = vx(i,j,k) - vbar*sx
                vy(i,j,k) = vy(i,j,k) - vbar*sy
                vz(i,j,k) = vz(i,j,k) - vbar*sz

                wx(i,j,k) = wx(i,j,k) - wbar*sx
                wy(i,j,k) = wy(i,j,k) - wbar*sy
                wz(i,j,k) = wz(i,j,k) - wbar*sz

                qx(i,j,k) = qx(i,j,k) + a2*sx
                qy(i,j,k) = qy(i,j,k) + a2*sy
                qz(i,j,k) = qz(i,j,k) + a2*sz
             endif
          end do
       enddo
    enddo


    ! Second part. Contribution in the j-direction.
    ! The contribution is scattered to both the left and right node
    ! in j-direction.

    do k=1, kl
       do j=1, je
          do i=1, il

             ! Compute 8 times the average normal for this part of
             ! the control volume. The factor 8 is taken care of later
             ! on when the division by the volume takes place.

             sx = sj(i,j-1,k,  1) + sj(i+1,j-1,k,  1) &
                  + sj(i,j-1,k+1,1) + sj(i+1,j-1,k+1,1) &
                  + sj(i,j,  k,  1) + sj(i+1,j,  k,  1) &
                  + sj(i,j,  k+1,1) + sj(i+1,j,  k+1,1)
             sy = sj(i,j-1,k,  2) + sj(i+1,j-1,k,  2) &
                  + sj(i,j-1,k+1,2) + sj(i+1,j-1,k+1,2) &
                  + sj(i,j,  k,  2) + sj(i+1,j,  k,  2) &
                  + sj(i,j,  k+1,2) + sj(i+1,j,  k+1,2)
             sz = sj(i,j-1,k,  3) + sj(i+1,j-1,k,  3) &
                  + sj(i,j-1,k+1,3) + sj(i+1,j-1,k+1,3) &
                  + sj(i,j,  k,  3) + sj(i+1,j,  k,  3) &
                  + sj(i,j,  k+1,3) + sj(i+1,j,  k+1,3)

             ! Compute the average velocities and speed of sound squared
             ! for this integration point. Node that these variables are
             ! stored in w(ivx), w(ivy), w(ivz) and p.

             ubar = fourth*(w(i,j,k,  ivx) + w(i+1,j,k,  ivx) &
                  +         w(i,j,k+1,ivx) + w(i+1,j,k+1,ivx))
             vbar = fourth*(w(i,j,k,  ivy) + w(i+1,j,k,  ivy) &
                  +         w(i,j,k+1,ivy) + w(i+1,j,k+1,ivy))
             wbar = fourth*(w(i,j,k,  ivz) + w(i+1,j,k,  ivz) &
                  +         w(i,j,k+1,ivz) + w(i+1,j,k+1,ivz))

             a2 = fourth*(aa(i,j,k) + aa(i+1,j,k) + aa(i,j,k+1) + aa(i+1,j,k+1))

             ! Add the contributions to the surface integral to the node
             ! j-1 and substract it from the node j. For the heat flux it
             ! is reversed, because the negative of the gradient of the
             ! speed of sound must be computed.

             if(j > 1) then
                ux(i,j-1,k) = ux(i,j-1,k) + ubar*sx
                uy(i,j-1,k) = uy(i,j-1,k) + ubar*sy
                uz(i,j-1,k) = uz(i,j-1,k) + ubar*sz

                vx(i,j-1,k) = vx(i,j-1,k) + vbar*sx
                vy(i,j-1,k) = vy(i,j-1,k) + vbar*sy
                vz(i,j-1,k) = vz(i,j-1,k) + vbar*sz

                wx(i,j-1,k) = wx(i,j-1,k) + wbar*sx
                wy(i,j-1,k) = wy(i,j-1,k) + wbar*sy
                wz(i,j-1,k) = wz(i,j-1,k) + wbar*sz

                qx(i,j-1,k) = qx(i,j-1,k) - a2*sx
                qy(i,j-1,k) = qy(i,j-1,k) - a2*sy
                qz(i,j-1,k) = qz(i,j-1,k) - a2*sz
             endif

             if(j < je) then
                ux(i,j,k) = ux(i,j,k) - ubar*sx
                uy(i,j,k) = uy(i,j,k) - ubar*sy
                uz(i,j,k) = uz(i,j,k) - ubar*sz

                vx(i,j,k) = vx(i,j,k) - vbar*sx
                vy(i,j,k) = vy(i,j,k) - vbar*sy
                vz(i,j,k) = vz(i,j,k) - vbar*sz

                wx(i,j,k) = wx(i,j,k) - wbar*sx
                wy(i,j,k) = wy(i,j,k) - wbar*sy
                wz(i,j,k) = wz(i,j,k) - wbar*sz

                qx(i,j,k) = qx(i,j,k) + a2*sx
                qy(i,j,k) = qy(i,j,k) + a2*sy
                qz(i,j,k) = qz(i,j,k) + a2*sz
             endif
          end do
       enddo
    enddo
    !
    ! Third part. Contribution in the i-direction.
    ! The contribution is scattered to both the left and right node
    ! in i-direction.
    !
    do k=1,kl
       do j=1,jl
          do i=1, ie

             ! Compute 8 times the average normal for this part of
             ! the control volume. The factor 8 is taken care of later
             ! on when the division by the volume takes place.

             sx = si(i-1,j,k,  1) + si(i-1,j+1,k,  1) &
                  + si(i-1,j,k+1,1) + si(i-1,j+1,k+1,1) &
                  + si(i,  j,k,  1) + si(i,  j+1,k,  1) &
                  + si(i,  j,k+1,1) + si(i,  j+1,k+1,1)
             sy = si(i-1,j,k,  2) + si(i-1,j+1,k,  2) &
                  + si(i-1,j,k+1,2) + si(i-1,j+1,k+1,2) &
                  + si(i,  j,k,  2) + si(i,  j+1,k,  2) &
                  + si(i,  j,k+1,2) + si(i,  j+1,k+1,2)
             sz = si(i-1,j,k,  3) + si(i-1,j+1,k,  3) &
                  + si(i-1,j,k+1,3) + si(i-1,j+1,k+1,3) &
                  + si(i,  j,k,  3) + si(i,  j+1,k,  3) &
                  + si(i,  j,k+1,3) + si(i,  j+1,k+1,3)

             ! Compute the average velocities and speed of sound squared
             ! for this integration point. Node that these variables are
             ! stored in w(ivx), w(ivy), w(ivz) and p.

             ubar = fourth*(w(i,j,k,  ivx) + w(i,j+1,k,  ivx) &
                  +         w(i,j,k+1,ivx) + w(i,j+1,k+1,ivx))
             vbar = fourth*(w(i,j,k,  ivy) + w(i,j+1,k,  ivy) &
                  +         w(i,j,k+1,ivy) + w(i,j+1,k+1,ivy))
             wbar = fourth*(w(i,j,k,  ivz) + w(i,j+1,k,  ivz) &
                  +         w(i,j,k+1,ivz) + w(i,j+1,k+1,ivz))

             a2 = fourth*(aa(i,j,k) + aa(i,j+1,k) + aa(i,j,k+1) + aa(i,j+1,k+1))

             ! Add the contributions to the surface integral to the node
             ! j-1 and substract it from the node j. For the heat flux it
             ! is reversed, because the negative of the gradient of the
             ! speed of sound must be computed.

             if(i > 1) then
                ux(i-1,j,k) = ux(i-1,j,k) + ubar*sx
                uy(i-1,j,k) = uy(i-1,j,k) + ubar*sy
                uz(i-1,j,k) = uz(i-1,j,k) + ubar*sz

                vx(i-1,j,k) = vx(i-1,j,k) + vbar*sx
                vy(i-1,j,k) = vy(i-1,j,k) + vbar*sy
                vz(i-1,j,k) = vz(i-1,j,k) + vbar*sz

                wx(i-1,j,k) = wx(i-1,j,k) + wbar*sx
                wy(i-1,j,k) = wy(i-1,j,k) + wbar*sy
                wz(i-1,j,k) = wz(i-1,j,k) + wbar*sz

                qx(i-1,j,k) = qx(i-1,j,k) - a2*sx
                qy(i-1,j,k) = qy(i-1,j,k) - a2*sy
                qz(i-1,j,k) = qz(i-1,j,k) - a2*sz
             endif

             if(i < ie) then
                ux(i,j,k) = ux(i,j,k) - ubar*sx
                uy(i,j,k) = uy(i,j,k) - ubar*sy
                uz(i,j,k) = uz(i,j,k) - ubar*sz

                vx(i,j,k) = vx(i,j,k) - vbar*sx
                vy(i,j,k) = vy(i,j,k) - vbar*sy
                vz(i,j,k) = vz(i,j,k) - vbar*sz

                wx(i,j,k) = wx(i,j,k) - wbar*sx
                wy(i,j,k) = wy(i,j,k) - wbar*sy
                wz(i,j,k) = wz(i,j,k) - wbar*sz

                qx(i,j,k) = qx(i,j,k) + a2*sx
                qy(i,j,k) = qy(i,j,k) + a2*sy
                qz(i,j,k) = qz(i,j,k) + a2*sz
             endif
          enddo
       enddo
    enddo

    ! Divide by 8 times the volume to obtain the correct gradients.

    do k=1,kl
       do j=1,jl
          do i=1,il

             ! Compute the inverse of 8 times the volume for this node.

             oVol = one/(vol(i,  j,  k) + vol(i,  j,  k+1) &
                  +      vol(i+1,j,  k) + vol(i+1,j,  k+1) &
                  +      vol(i,  j+1,k) + vol(i,  j+1,k+1) &
                  +      vol(i+1,j+1,k) + vol(i+1,j+1,k+1))

             ! Compute the correct velocity gradients and "unit" heat
             ! fluxes. The velocity gradients are stored in ux, etc.

             ux(i,j,k) = ux(i,j,k)*oVol
             uy(i,j,k) = uy(i,j,k)*oVol
             uz(i,j,k) = uz(i,j,k)*oVol

             vx(i,j,k) = vx(i,j,k)*oVol
             vy(i,j,k) = vy(i,j,k)*oVol
             vz(i,j,k) = vz(i,j,k)*oVol

             wx(i,j,k) = wx(i,j,k)*oVol
             wy(i,j,k) = wy(i,j,k)*oVol
             wz(i,j,k) = wz(i,j,k)*oVol

             qx(i,j,k) = qx(i,j,k)*oVol
             qy(i,j,k) = qy(i,j,k)*oVol
             qz(i,j,k) = qz(i,j,k)*oVol
          end do
       enddo
    enddo
  end subroutine allNodalGradients

  subroutine viscousFlux
    ! ---------------------------------------------
    !                   Viscous Flux
    ! ---------------------------------------------

    use constants
    use inputPhysics, only : useQCR, prandtl, prandtlturb
    use flowvarRefState, only : eddyModel
    use iteration, only : rFil
    implicit none

    ! Variables for viscous flux
    real(kind=realType) :: rFilv, por, mul, mue, mut, heatCoef
    real(kind=realType) :: gm1, factLamHeat, factTurbHeat
    real(kind=realType) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
    real(kind=realType) :: q_x, q_y, q_z
    real(kind=realType) :: corr, ssx, ssy, ssz, fracDiv, snrm
    real(kind=realType) :: tauxx, tauyy, tauzz
    real(kind=realType) :: tauxy, tauxz, tauyz
    real(kind=realType) :: tauxxS, tauyyS, tauzzS
    real(kind=realType) :: tauxyS, tauxzS, tauyzS
    real(kind=realType) :: ubar, vbar, wbar
    real(kind=realType) :: exx, eyy, ezz
    real(kind=realType) :: exy, exz, eyz
    real(kind=realType) :: Wxx, Wyy, Wzz
    real(kind=realType) :: Wxy, Wxz, Wyz, Wyx, Wzx, Wzy
    real(kind=realType) :: den, Ccr1
    real(kind=realType) :: fmx, fmy, fmz, frhoE, fact
    integer(kind=intType) :: i, j, k
    real(kind=realType), parameter :: xminn = 1.e-10_realType
    real(kind=realType), parameter :: twoThird = two*third

    ! Set QCR parameters
    Ccr1 = 0.3_realType
    rFilv = rFil
    
    ! The diagonals of the vorticity tensor components are always zero
    Wxx = zero
    Wyy = zero
    Wzz = zero
    !
    !         viscous fluxes in the k-direction.                           
    !
    mue = zero
    do k=1,kl
       do j=2,jl
          do i=2,il


                ! Set the value of the porosity. If not zero, it is set
                ! to average the eddy-viscosity and to take the factor
                ! rFilv into account.

                por = half*rFilv
                if(porK(i,j,k) == noFlux) por = zero

                ! Compute the laminar and (if present) the eddy viscosities
                ! multiplied by the porosity. Compute the factor in front of
                ! the gradients of the speed of sound squared for the heat
                ! flux.

                mul = por*(rlv(i,j,k) + rlv(i,j,k+1))
                mue = por*(rev(i,j,k) + rev(i,j,k+1))
                mut = mul + mue

                gm1          = half*(gamma(i,j,k) + gamma(i,j,k+1)) - one
                factLamHeat  = one/(prandtl*gm1)
                factTurbHeat = one/(prandtlTurb*gm1)

                heatCoef = mul*factLamHeat + mue*factTurbHeat

                ! Compute the gradients at the face by averaging the four
                ! nodal values.

                u_x = fourth*(ux(i-1,j-1,k) + ux(i,j-1,k) &
                     +         ux(i-1,j,  k) + ux(i,j,  k))
                u_y = fourth*(uy(i-1,j-1,k) + uy(i,j-1,k) &
                     +         uy(i-1,j,  k) + uy(i,j,  k))
                u_z = fourth*(uz(i-1,j-1,k) + uz(i,j-1,k) &
                     +         uz(i-1,j,  k) + uz(i,j,  k))

                v_x = fourth*(vx(i-1,j-1,k) + vx(i,j-1,k) &
                     +         vx(i-1,j,  k) + vx(i,j,  k))
                v_y = fourth*(vy(i-1,j-1,k) + vy(i,j-1,k) &
                     +         vy(i-1,j,  k) + vy(i,j,  k))
                v_z = fourth*(vz(i-1,j-1,k) + vz(i,j-1,k) &
                     +         vz(i-1,j,  k) + vz(i,j,  k))

                w_x = fourth*(wx(i-1,j-1,k) + wx(i,j-1,k) &
                     +         wx(i-1,j,  k) + wx(i,j,  k))
                w_y = fourth*(wy(i-1,j-1,k) + wy(i,j-1,k) &
                     +         wy(i-1,j,  k) + wy(i,j,  k))
                w_z = fourth*(wz(i-1,j-1,k) + wz(i,j-1,k) &
                     +         wz(i-1,j,  k) + wz(i,j,  k))

                q_x = fourth*(qx(i-1,j-1,k) + qx(i,j-1,k) &
                     +         qx(i-1,j,  k) + qx(i,j,  k))
                q_y = fourth*(qy(i-1,j-1,k) + qy(i,j-1,k) &
                     +         qy(i-1,j,  k) + qy(i,j,  k))
                q_z = fourth*(qz(i-1,j-1,k) + qz(i,j-1,k) &
                     +         qz(i-1,j,  k) + qz(i,j,  k))


                ! The gradients in the normal direction are corrected, such
                ! that no averaging takes places here.
                ! First determine the vector in the direction from the
                ! cell center k to cell center k+1.

                ssx = eighth*(x(i-1,j-1,k+1,1) - x(i-1,j-1,k-1,1) &
                     +         x(i-1,j,  k+1,1) - x(i-1,j,  k-1,1) &
                     +         x(i,  j-1,k+1,1) - x(i,  j-1,k-1,1) &
                     +         x(i,  j,  k+1,1) - x(i,  j,  k-1,1))
                ssy = eighth*(x(i-1,j-1,k+1,2) - x(i-1,j-1,k-1,2) &
                     +         x(i-1,j,  k+1,2) - x(i-1,j,  k-1,2) &
                     +         x(i,  j-1,k+1,2) - x(i,  j-1,k-1,2) &
                     +         x(i,  j,  k+1,2) - x(i,  j,  k-1,2))
                ssz = eighth*(x(i-1,j-1,k+1,3) - x(i-1,j-1,k-1,3) &
                     +         x(i-1,j,  k+1,3) - x(i-1,j,  k-1,3) &
                     +         x(i,  j-1,k+1,3) - x(i,  j-1,k-1,3) &
                     +         x(i,  j,  k+1,3) - x(i,  j,  k-1,3))

                ! Determine the length of this vector and create the
                ! unit normal.

                snrm  = one/sqrt(ssx*ssx + ssy*ssy + ssz*ssz)
                ssx = snrm*ssx
                ssy = snrm*ssy
                ssz = snrm*ssz

                ! Correct the gradients.

                corr = u_x*ssx + u_y*ssy + u_z*ssz        &
                     - (w(i,j,k+1,ivx) - w(i,j,k,ivx))*snrm
                u_x  = u_x - corr*ssx
                u_y  = u_y - corr*ssy
                u_z  = u_z - corr*ssz

                corr = v_x*ssx + v_y*ssy + v_z*ssz        &
                     - (w(i,j,k+1,ivy) - w(i,j,k,ivy))*snrm
                v_x  = v_x - corr*ssx
                v_y  = v_y - corr*ssy
                v_z  = v_z - corr*ssz

                corr = w_x*ssx + w_y*ssy + w_z*ssz        &
                     - (w(i,j,k+1,ivz) - w(i,j,k,ivz))*snrm
                w_x  = w_x - corr*ssx
                w_y  = w_y - corr*ssy
                w_z  = w_z - corr*ssz

                corr = q_x*ssx + q_y*ssy + q_z*ssz &
                     + (aa(i,j,k+1) - aa(i,j,k))*snrm
                q_x  = q_x - corr*ssx
                q_y  = q_y - corr*ssy
                q_z  = q_z - corr*ssz

                ! Compute the stress tensor and the heat flux vector.
                ! We remove the viscosity from the stress tensor (tau)
                ! to define tauS since we still need to separate between
                ! laminar and turbulent stress for QCR.
                ! Therefore, laminar tau = mue*tauS, turbulent
                ! tau = mue*tauS, and total tau = mut*tauS.

                fracDiv = twoThird*(u_x + v_y + w_z)

                tauxxS = two*u_x - fracDiv
                tauyyS = two*v_y - fracDiv
                tauzzS = two*w_z - fracDiv

                tauxyS = u_y + v_x
                tauxzS = u_z + w_x
                tauyzS = v_z + w_y

                q_x = heatCoef*q_x
                q_y = heatCoef*q_y
                q_z = heatCoef*q_z

                ! Add QCR corrections if necessary
                if (useQCR) then

                   ! In the QCR formulation, we add an extra term to the turbulent stress tensor:
                   !
                   ! tau_ij,QCR = tau_ij - e_ij
                   !
                   ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
                   !
                   ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
                   !
                   ! We are computing O_ik as follows:
                   !
                   ! O_ik = 2*W_ik/den
                   !
                   ! Remember that the tau_ij in e_ij should use only the eddy viscosity!

                   ! Compute denominator
                   den = sqrt(u_x*u_x + u_y*u_y + u_z*u_z + &
                        v_x*v_x + v_y*v_y + v_z*v_z + &
                        w_x*w_x + w_y*w_y + w_z*w_z)

                   ! Denominator should be limited to avoid division by zero in regions with
                   ! no gradients
                   den = max(den, xminn)

                   ! Compute factor that will multiply all tensor components.
                   ! Here we add the eddy viscosity that should multiply the stress tensor (tau)
                   ! components as well.
                   fact = mue*Ccr1/den

                   ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
                   ! The diagonals of the vorticity tensor components are always zero
                   Wxy = u_y - v_x
                   Wxz = u_z - w_x
                   Wyz = v_z - w_y
                   Wyx = -Wxy
                   Wzx = -Wxz
                   Wzy = -Wyz

                   ! Compute the extra terms of the Boussinesq relation
                   exx = fact*(Wxy*tauxyS + Wxz*tauxzS)*two
                   eyy = fact*(Wyx*tauxyS + Wyz*tauyzS)*two
                   ezz = fact*(Wzx*tauxzS + Wzy*tauyzS)*two

                   exy = fact*(Wxy*tauyyS + Wxz*tauyzS + &
                               Wyx*tauxxS + Wyz*tauxzS)
                   exz = fact*(Wxy*tauyzS + Wxz*tauzzS + &
                               Wzx*tauxxS + Wzy*tauxyS)
                   eyz = fact*(Wyx*tauxzS + Wyz*tauzzS + &
                               Wzx*tauxyS + Wzy*tauyyS)

                   ! Apply the total viscosity to the stress tensor and add extra terms
                   tauxx = mut*tauxxS - exx
                   tauyy = mut*tauyyS - eyy
                   tauzz = mut*tauzzS - ezz
                   tauxy = mut*tauxyS - exy
                   tauxz = mut*tauxzS - exz
                   tauyz = mut*tauyzS - eyz

                else

                   ! Just apply the total viscosity to the stress tensor
                   tauxx = mut*tauxxS
                   tauyy = mut*tauyyS
                   tauzz = mut*tauzzS
                   tauxy = mut*tauxyS
                   tauxz = mut*tauxzS
                   tauyz = mut*tauyzS

                end if

                ! Compute the average velocities for the face. Remember that
                ! the velocities are stored and not the momentum.

                ubar = half*(w(i,j,k,ivx) + w(i,j,k+1,ivx))
                vbar = half*(w(i,j,k,ivy) + w(i,j,k+1,ivy))
                wbar = half*(w(i,j,k,ivz) + w(i,j,k+1,ivz))

                ! Compute the viscous fluxes for this k-face.

                fmx   = tauxx*sk(i,j,k,1) + tauxy*sk(i,j,k,2) &
                     + tauxz*sk(i,j,k,3)
                fmy   = tauxy*sk(i,j,k,1) + tauyy*sk(i,j,k,2) &
                     + tauyz*sk(i,j,k,3)
                fmz   = tauxz*sk(i,j,k,1) + tauyz*sk(i,j,k,2) &
                     + tauzz*sk(i,j,k,3)
                frhoE =         (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sk(i,j,k,1)
                frhoE = frhoE + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sk(i,j,k,2)
                frhoE = frhoE + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sk(i,j,k,3)
                frhoE = frhoE -  q_x*sk(i,j,k,1) - q_y*sk(i,j,k,2) - q_z*sk(i,j,k,3)

                ! Update the residuals of cell k and k+1.

                fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
                fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
                fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

                fw(i,j,k+1,imx)   = fw(i,j,k+1,imx)   + fmx
                fw(i,j,k+1,imy)   = fw(i,j,k+1,imy)   + fmy
                fw(i,j,k+1,imz)   = fw(i,j,k+1,imz)   + fmz
                fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + frhoE

          end do
       enddo
    end do
    !
    !         Viscous fluxes in the j-direction.                           
    !
    do k=2,kl
       do j=1,jl
          do i=2,il


                ! Set the value of the porosity. If not zero, it is set
                ! to average the eddy-viscosity and to take the factor
                ! rFilv into account.

                por = half*rFilv
                if(porJ(i,j,k) == noFlux) por = zero

                ! Compute the laminar and (if present) the eddy viscosities
                ! multiplied by the porosity. Compute the factor in front of
                ! the gradients of the speed of sound squared for the heat
                ! flux.

                mul = por*(rlv(i,j,k) + rlv(i,j+1,k))
                mue = por*(rev(i,j,k) + rev(i,j+1,k))
                mut = mul + mue

                gm1          = half*(gamma(i,j,k) + gamma(i,j+1,k)) - one
                factLamHeat  = one/(prandtl*gm1)
                factTurbHeat = one/(prandtlTurb*gm1)

                heatCoef = mul*factLamHeat + mue*factTurbHeat

                ! Compute the gradients at the face by averaging the four
                ! nodal values.

                u_x = fourth*(ux(i-1,j,k-1) + ux(i,j,k-1) &
                     +        ux(i-1,j,k  ) + ux(i,j,k  ))
                u_y = fourth*(uy(i-1,j,k-1) + uy(i,j,k-1) &
                     +        uy(i-1,j,k  ) + uy(i,j,k  ))
                u_z = fourth*(uz(i-1,j,k-1) + uz(i,j,k-1) &
                     +        uz(i-1,j,k  ) + uz(i,j,k  ))

                v_x = fourth*(vx(i-1,j,k-1) + vx(i,j,k-1) &
                     +        vx(i-1,j,k  ) + vx(i,j,k  ))
                v_y = fourth*(vy(i-1,j,k-1) + vy(i,j,k-1) &
                     +        vy(i-1,j,k  ) + vy(i,j,k  ))
                v_z = fourth*(vz(i-1,j,k-1) + vz(i,j,k-1) &
                     +        vz(i-1,j,k  ) + vz(i,j,k  ))

                w_x = fourth*(wx(i-1,j,k-1) + wx(i,j,k-1) &
                     +        wx(i-1,j,k  ) + wx(i,j,k  ))
                w_y = fourth*(wy(i-1,j,k-1) + wy(i,j,k-1) &
                     +        wy(i-1,j,k  ) + wy(i,j,k  ))
                w_z = fourth*(wz(i-1,j,k-1) + wz(i,j,k-1) &
                     +        wz(i-1,j,k  ) + wz(i,j,k  ))

                q_x = fourth*(qx(i-1,j,k-1) + qx(i,j,k-1) &
                     +        qx(i-1,j,k  ) + qx(i,j,k  ))
                q_y = fourth*(qy(i-1,j,k-1) + qy(i,j,k-1) &
                     +        qy(i-1,j,k  ) + qy(i,j,k  ))
                q_z = fourth*(qz(i-1,j,k-1) + qz(i,j,k-1) &
                     +        qz(i-1,j,k  ) + qz(i,j,k  ))

                ! The gradients in the normal direction are corrected, such
                ! that no averaging takes places here.
                ! First determine the vector in the direction from the
                ! cell center j to cell center j+1.

                ssx = eighth*(x(i-1,j+1,k-1,1) - x(i-1,j-1,k-1,1) &
                     +         x(i-1,j+1,k,  1) - x(i-1,j-1,k,  1) &
                     +         x(i,  j+1,k-1,1) - x(i,  j-1,k-1,1) &
                     +         x(i,  j+1,k,  1) - x(i,  j-1,k,  1))
                ssy = eighth*(x(i-1,j+1,k-1,2) - x(i-1,j-1,k-1,2) &
                     +         x(i-1,j+1,k,  2) - x(i-1,j-1,k,  2) &
                     +         x(i,  j+1,k-1,2) - x(i,  j-1,k-1,2) &
                     +         x(i,  j+1,k,  2) - x(i,  j-1,k,  2))
                ssz = eighth*(x(i-1,j+1,k-1,3) - x(i-1,j-1,k-1,3) &
                     +         x(i-1,j+1,k,  3) - x(i-1,j-1,k,  3) &
                     +         x(i,  j+1,k-1,3) - x(i,  j-1,k-1,3) &
                     +         x(i,  j+1,k,  3) - x(i,  j-1,k,  3))

                ! Determine the length of this vector and create the
                ! unit normal.

                snrm  = one/sqrt(ssx*ssx + ssy*ssy + ssz*ssz)
                ssx = snrm*ssx
                ssy = snrm*ssy
                ssz = snrm*ssz

                ! Correct the gradients.

                corr = u_x*ssx + u_y*ssy + u_z*ssz        &
                     - (w(i,j+1,k,ivx) - w(i,j,k,ivx))*snrm
                u_x  = u_x - corr*ssx
                u_y  = u_y - corr*ssy
                u_z  = u_z - corr*ssz

                corr = v_x*ssx + v_y*ssy + v_z*ssz        &
                     - (w(i,j+1,k,ivy) - w(i,j,k,ivy))*snrm
                v_x  = v_x - corr*ssx
                v_y  = v_y - corr*ssy
                v_z  = v_z - corr*ssz

                corr = w_x*ssx + w_y*ssy + w_z*ssz        &
                     - (w(i,j+1,k,ivz) - w(i,j,k,ivz))*snrm
                w_x  = w_x - corr*ssx
                w_y  = w_y - corr*ssy
                w_z  = w_z - corr*ssz

                corr = q_x*ssx + q_y*ssy + q_z*ssz &
                     + (aa(i,j+1,k) - aa(i,j,k))*snrm
                q_x  = q_x - corr*ssx
                q_y  = q_y - corr*ssy
                q_z  = q_z - corr*ssz

                ! Compute the stress tensor and the heat flux vector.
                ! We remove the viscosity from the stress tensor (tau)
                ! to define tauS since we still need to separate between
                ! laminar and turbulent stress for QCR.
                ! Therefore, laminar tau = mue*tauS, turbulent
                ! tau = mue*tauS, and total tau = mut*tauS.

                fracDiv = twoThird*(u_x + v_y + w_z)

                tauxxS = two*u_x - fracDiv
                tauyyS = two*v_y - fracDiv
                tauzzS = two*w_z - fracDiv

                tauxyS = u_y + v_x
                tauxzS = u_z + w_x
                tauyzS = v_z + w_y

                q_x = heatCoef*q_x
                q_y = heatCoef*q_y
                q_z = heatCoef*q_z

                ! Add QCR corrections if necessary
                if (useQCR) then

                   ! In the QCR formulation, we add an extra term to the turbulent stress tensor:
                   !
                   ! tau_ij,QCR = tau_ij - e_ij
                   !
                   ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
                   !
                   ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
                   !
                   ! We are computing O_ik as follows:
                   !
                   ! O_ik = 2*W_ik/den
                   !
                   ! Remember that the tau_ij in e_ij should use only the eddy viscosity!

                   ! Compute denominator
                   den = sqrt(u_x*u_x + u_y*u_y + u_z*u_z + &
                        v_x*v_x + v_y*v_y + v_z*v_z + &
                        w_x*w_x + w_y*w_y + w_z*w_z)

                   ! Denominator should be limited to avoid division by zero in regions with
                   ! no gradients
                   den = max(den, xminn)

                   ! Compute factor that will multiply all tensor components.
                   ! Here we add the eddy viscosity that should multiply the stress tensor (tau)
                   ! components as well.
                   fact = mue*Ccr1/den

                   ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
                   ! The diagonals of the vorticity tensor components are always zero
                   Wxy = u_y - v_x
                   Wxz = u_z - w_x
                   Wyz = v_z - w_y
                   Wyx = -Wxy
                   Wzx = -Wxz
                   Wzy = -Wyz

                   ! Compute the extra terms of the Boussinesq relation
                   exx = fact*(Wxy*tauxyS + Wxz*tauxzS)*two
                   eyy = fact*(Wyx*tauxyS + Wyz*tauyzS)*two
                   ezz = fact*(Wzx*tauxzS + Wzy*tauyzS)*two

                   exy = fact*(Wxy*tauyyS + Wxz*tauyzS + &
                               Wyx*tauxxS + Wyz*tauxzS)
                   exz = fact*(Wxy*tauyzS + Wxz*tauzzS + &
                               Wzx*tauxxS + Wzy*tauxyS)
                   eyz = fact*(Wyx*tauxzS + Wyz*tauzzS + &
                               Wzx*tauxyS + Wzy*tauyyS)

                   ! Apply the total viscosity to the stress tensor and add extra terms
                   tauxx = mut*tauxxS - exx
                   tauyy = mut*tauyyS - eyy
                   tauzz = mut*tauzzS - ezz
                   tauxy = mut*tauxyS - exy
                   tauxz = mut*tauxzS - exz
                   tauyz = mut*tauyzS - eyz

                else

                   ! Just apply the total viscosity to the stress tensor
                   tauxx = mut*tauxxS
                   tauyy = mut*tauyyS
                   tauzz = mut*tauzzS
                   tauxy = mut*tauxyS
                   tauxz = mut*tauxzS
                   tauyz = mut*tauyzS

                end if

                ! Compute the average velocities for the face. Remember that
                ! the velocities are stored and not the momentum.

                ubar = half*(w(i,j,k,ivx) + w(i,j+1,k,ivx))
                vbar = half*(w(i,j,k,ivy) + w(i,j+1,k,ivy))
                wbar = half*(w(i,j,k,ivz) + w(i,j+1,k,ivz))

                ! Compute the viscous fluxes for this j-face.

                fmx   = tauxx*sj(i,j,k,1) + tauxy*sj(i,j,k,2) &
                     + tauxz*sj(i,j,k,3)
                fmy   = tauxy*sj(i,j,k,1) + tauyy*sj(i,j,k,2) &
                     + tauyz*sj(i,j,k,3)
                fmz   = tauxz*sj(i,j,k,1) + tauyz*sj(i,j,k,2) &
                     + tauzz*sj(i,j,k,3)
                frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sj(i,j,k,1) &
                     + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sj(i,j,k,2) &
                     + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sj(i,j,k,3) &
                     - q_x*sj(i,j,k,1) - q_y*sj(i,j,k,2) - q_z*sj(i,j,k,3)

                ! Update the residuals of cell j and j+1.

                fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
                fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
                fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

                fw(i,j+1,k,imx)   = fw(i,j+1,k,imx)   + fmx
                fw(i,j+1,k,imy)   = fw(i,j+1,k,imy)   + fmy
                fw(i,j+1,k,imz)   = fw(i,j+1,k,imz)   + fmz
                fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + frhoE

          enddo
       enddo
    enddo
    !
    !         Viscous fluxes in the i-direction.                           
    !
    do k=2, kl
       do j=2, jl
          do i=1, il
                ! Set the value of the porosity. If not zero, it is set
                ! to average the eddy-viscosity and to take the factor
                ! rFilv into account.

                por = half*rFilv
                if(porI(i,j,k) == noFlux) por = zero

                ! Compute the laminar and (if present) the eddy viscosities
                ! multiplied the porosity. Compute the factor in front of
                ! the gradients of the speed of sound squared for the heat
                ! flux.

                mul = por*(rlv(i,j,k) + rlv(i+1,j,k))
                mue = por*(rev(i,j,k) + rev(i+1,j,k))
                mut = mul + mue

                gm1          = half*(gamma(i,j,k) + gamma(i+1,j,k)) - one
                factLamHeat  = one/(prandtl*gm1)
                factTurbHeat = one/(prandtlTurb*gm1)

                heatCoef = mul*factLamHeat + mue*factTurbHeat

                ! Compute the gradients at the face by averaging the four
                ! nodal values.

                u_x = fourth*(ux(i,j-1,k-1) + ux(i,j,k-1) &
                     +        ux(i,j-1,k  ) + ux(i,j,k  ))
                u_y = fourth*(uy(i,j-1,k-1) + uy(i,j,k-1) &
                     +        uy(i,j-1,k  ) + uy(i,j,k  ))
                u_z = fourth*(uz(i,j-1,k-1) + uz(i,j,k-1) &
                     +        uz(i,j-1,k  ) + uz(i,j,k  ))

                v_x = fourth*(vx(i,j-1,k-1) + vx(i,j,k-1) &
                     +        vx(i,j-1,k  ) + vx(i,j,k  ))
                v_y = fourth*(vy(i,j-1,k-1) + vy(i,j,k-1) &
                     +        vy(i,j-1,k  ) + vy(i,j,k  ))
                v_z = fourth*(vz(i,j-1,k-1) + vz(i,j,k-1) &
                     +        vz(i,j-1,k  ) + vz(i,j,k  ))

                w_x = fourth*(wx(i,j-1,k-1) + wx(i,j,k-1) &
                     +        wx(i,j-1,k  ) + wx(i,j,k  ))
                w_y = fourth*(wy(i,j-1,k-1) + wy(i,j,k-1) &
                     +        wy(i,j-1,k  ) + wy(i,j,k  ))
                w_z = fourth*(wz(i,j-1,k-1) + wz(i,j,k-1) &
                     +        wz(i,j-1,k  ) + wz(i,j,k  ))

                q_x = fourth*(qx(i,j-1,k-1) + qx(i,j,k-1) &
                     +        qx(i,j-1,k  ) + qx(i,j,k  ))
                q_y = fourth*(qy(i,j-1,k-1) + qy(i,j,k-1) &
                     +        qy(i,j-1,k  ) + qy(i,j,k  ))
                q_z = fourth*(qz(i,j-1,k-1) + qz(i,j,k-1) &
                     +        qz(i,j-1,k  ) + qz(i,j,k  ))

                ! The gradients in the normal direction are corrected, such
                ! that no averaging takes places here.
                ! First determine the vector in the direction from the
                ! cell center i to cell center i+1.

                ssx = eighth*(x(i+1,j-1,k-1,1) - x(i-1,j-1,k-1,1) &
                     +         x(i+1,j-1,k,  1) - x(i-1,j-1,k,  1) &
                     +         x(i+1,j,  k-1,1) - x(i-1,j,  k-1,1) &
                     +         x(i+1,j,  k,  1) - x(i-1,j,  k,  1))
                ssy = eighth*(x(i+1,j-1,k-1,2) - x(i-1,j-1,k-1,2) &
                     +         x(i+1,j-1,k,  2) - x(i-1,j-1,k,  2) &
                     +         x(i+1,j,  k-1,2) - x(i-1,j,  k-1,2) &
                     +         x(i+1,j,  k,  2) - x(i-1,j,  k,  2))
                ssz = eighth*(x(i+1,j-1,k-1,3) - x(i-1,j-1,k-1,3) &
                     +         x(i+1,j-1,k,  3) - x(i-1,j-1,k,  3) &
                     +         x(i+1,j,  k-1,3) - x(i-1,j,  k-1,3) &
                     +         x(i+1,j,  k,  3) - x(i-1,j,  k,  3))

                ! Determine the length of this vector and create the
                ! unit normal.

                snrm  = one/sqrt(ssx*ssx + ssy*ssy + ssz*ssz)
                ssx = snrm*ssx
                ssy = snrm*ssy
                ssz = snrm*ssz

                ! Correct the gradients.

                corr = u_x*ssx + u_y*ssy + u_z*ssz        &
                     - (w(i+1,j,k,ivx) - w(i,j,k,ivx))*snrm
                u_x  = u_x - corr*ssx
                u_y  = u_y - corr*ssy
                u_z  = u_z - corr*ssz

                corr = v_x*ssx + v_y*ssy + v_z*ssz        &
                     - (w(i+1,j,k,ivy) - w(i,j,k,ivy))*snrm
                v_x  = v_x - corr*ssx
                v_y  = v_y - corr*ssy
                v_z  = v_z - corr*ssz

                corr = w_x*ssx + w_y*ssy + w_z*ssz        &
                     - (w(i+1,j,k,ivz) - w(i,j,k,ivz))*snrm
                w_x  = w_x - corr*ssx
                w_y  = w_y - corr*ssy
                w_z  = w_z - corr*ssz

                corr = q_x*ssx + q_y*ssy + q_z*ssz &
                     + (aa(i+1,j,k) - aa(i,j,k))*snrm
                q_x  = q_x - corr*ssx
                q_y  = q_y - corr*ssy
                q_z  = q_z - corr*ssz

                ! Compute the stress tensor and the heat flux vector.
                ! We remove the viscosity from the stress tensor (tau)
                ! to define tauS since we still need to separate between
                ! laminar and turbulent stress for QCR.
                ! Therefore, laminar tau = mue*tauS, turbulent
                ! tau = mue*tauS, and total tau = mut*tauS.

                fracDiv = twoThird*(u_x + v_y + w_z)

                tauxxS = two*u_x - fracDiv
                tauyyS = two*v_y - fracDiv
                tauzzS = two*w_z - fracDiv

                tauxyS = u_y + v_x
                tauxzS = u_z + w_x
                tauyzS = v_z + w_y

                q_x = heatCoef*q_x
                q_y = heatCoef*q_y
                q_z = heatCoef*q_z

                ! Add QCR corrections if necessary
                if (useQCR) then

                   ! In the QCR formulation, we add an extra term to the turbulent stress tensor:
                   !
                   ! tau_ij,QCR = tau_ij - e_ij
                   !
                   ! where, according to TMR website (http://turbmodels.larc.nasa.gov/spalart.html):
                   !
                   ! e_ij = Ccr1*(O_ik*tau_jk + O_jk*tau_ik)
                   !
                   ! We are computing O_ik as follows:
                   !
                   ! O_ik = 2*W_ik/den
                   !
                   ! Remember that the tau_ij in e_ij should use only the eddy viscosity!

                   ! Compute denominator
                   den = sqrt(u_x*u_x + u_y*u_y + u_z*u_z + &
                        v_x*v_x + v_y*v_y + v_z*v_z + &
                        w_x*w_x + w_y*w_y + w_z*w_z)

                   ! Denominator should be limited to avoid division by zero in regions with
                   ! no gradients
                   den = max(den, xminn)

                   ! Compute factor that will multiply all tensor components.
                   ! Here we add the eddy viscosity that should multiply the stress tensor (tau)
                   ! components as well.
                   fact = mue*Ccr1/den

                   ! Compute off-diagonal terms of vorticity tensor (we will ommit the 1/2)
                   ! The diagonals of the vorticity tensor components are always zero
                   Wxy = u_y - v_x
                   Wxz = u_z - w_x
                   Wyz = v_z - w_y
                   Wyx = -Wxy
                   Wzx = -Wxz
                   Wzy = -Wyz

                   ! Compute the extra terms of the Boussinesq relation
                   exx = fact*(Wxy*tauxyS + Wxz*tauxzS)*two
                   eyy = fact*(Wyx*tauxyS + Wyz*tauyzS)*two
                   ezz = fact*(Wzx*tauxzS + Wzy*tauyzS)*two

                   exy = fact*(Wxy*tauyyS + Wxz*tauyzS + &
                               Wyx*tauxxS + Wyz*tauxzS)
                   exz = fact*(Wxy*tauyzS + Wxz*tauzzS + &
                               Wzx*tauxxS + Wzy*tauxyS)
                   eyz = fact*(Wyx*tauxzS + Wyz*tauzzS + &
                               Wzx*tauxyS + Wzy*tauyyS)

                   ! Apply the total viscosity to the stress tensor and add extra terms
                   tauxx = mut*tauxxS - exx
                   tauyy = mut*tauyyS - eyy
                   tauzz = mut*tauzzS - ezz
                   tauxy = mut*tauxyS - exy
                   tauxz = mut*tauxzS - exz
                   tauyz = mut*tauyzS - eyz

                else

                   ! Just apply the total viscosity to the stress tensor
                   tauxx = mut*tauxxS
                   tauyy = mut*tauyyS
                   tauzz = mut*tauzzS
                   tauxy = mut*tauxyS
                   tauxz = mut*tauxzS
                   tauyz = mut*tauyzS

                end if

                ! Compute the average velocities for the face. Remember that
                ! the velocities are stored and not the momentum.

                ubar = half*(w(i,j,k,ivx) + w(i+1,j,k,ivx))
                vbar = half*(w(i,j,k,ivy) + w(i+1,j,k,ivy))
                wbar = half*(w(i,j,k,ivz) + w(i+1,j,k,ivz))

                ! Compute the viscous fluxes for this i-face.

                fmx   = tauxx*si(i,j,k,1) + tauxy*si(i,j,k,2) &
                     + tauxz*si(i,j,k,3)
                fmy   = tauxy*si(i,j,k,1) + tauyy*si(i,j,k,2) &
                     + tauyz*si(i,j,k,3)
                fmz   = tauxz*si(i,j,k,1) + tauyz*si(i,j,k,2) &
                     + tauzz*si(i,j,k,3)
                frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*si(i,j,k,1) &
                     + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*si(i,j,k,2) &
                     + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*si(i,j,k,3) &
                     - q_x*si(i,j,k,1) - q_y*si(i,j,k,2) - q_z*si(i,j,k,3)

                ! Update the residuals of cell i and i+1.

                fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
                fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
                fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
                fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

                fw(i+1,j,k,imx)   = fw(i+1,j,k,imx)   + fmx
                fw(i+1,j,k,imy)   = fw(i+1,j,k,imy)   + fmy
                fw(i+1,j,k,imz)   = fw(i+1,j,k,imz)   + fmz
                fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + frhoE

          enddo
       enddo
    enddo

  end subroutine viscousFlux

  subroutine viscousFluxApprox

    use constants
    use flowVarRefState
    use inputPhysics
    use iteration
    implicit none
    !
    !      Local parameter.
    !
    real(kind=realType), parameter :: twoThird = two*third
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k
    integer(kind=intType) :: ii, jj, kk

    real(kind=realType) :: rFilv, por, mul, mue, mut, heatCoef
    real(kind=realType) :: gm1, factLamHeat, factTurbHeat
    real(kind=realType) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
    real(kind=realType) :: q_x, q_y, q_z, ubar, vbar, wbar
    real(kind=realType) :: corr, ssx, ssy, ssz, ss, fracDiv
    real(kind=realType) :: tauxx, tauyy, tauzz
    real(kind=realType) :: tauxy, tauxz, tauyz
    real(kind=realType) :: fmx, fmy, fmz, frhoE
    real(kind=realType) :: dd
    logical :: correctForK

    mue = zero
    rFilv = rFil

    ! Viscous fluxes in the I-direction

    do k=2,kl
       do j=2,jl
          do i=1,il

             ! Compute the vector from the center of cell i to cell i+1
             ssx = eighth*(x(i+1,j-1,k-1,1) - x(i-1,j-1,k-1,1) &
                  +         x(i+1,j-1,k,  1) - x(i-1,j-1,k,  1) &
                  +         x(i+1,j,  k-1,1) - x(i-1,j,  k-1,1) &
                  +         x(i+1,j,  k,  1) - x(i-1,j,  k,  1))
             ssy = eighth*(x(i+1,j-1,k-1,2) - x(i-1,j-1,k-1,2) &
                  +         x(i+1,j-1,k,  2) - x(i-1,j-1,k,  2) &
                  +         x(i+1,j,  k-1,2) - x(i-1,j,  k-1,2) &
                  +         x(i+1,j,  k,  2) - x(i-1,j,  k,  2))
             ssz = eighth*(x(i+1,j-1,k-1,3) - x(i-1,j-1,k-1,3) &
                  +         x(i+1,j-1,k,  3) - x(i-1,j-1,k,  3) &
                  +         x(i+1,j,  k-1,3) - x(i-1,j,  k-1,3) &
                  +         x(i+1,j,  k,  3) - x(i-1,j,  k,  3))

             ! And determine one/ length of vector squared
             ss  = one/(ssx*ssx + ssy*ssy + ssz*ssz)
             ssx = ss*ssx
             ssy = ss*ssy
             ssz = ss*ssz

             ! Now compute each gradient
             dd = w(i+1, j, k, ivx)-w(i, j, k, ivx)
             u_x = dd*ssx
             u_y = dd*ssy
             u_z = dd*ssz

             dd = w(i+1, j, k, ivy)-w(i, j, k, ivy)
             v_x = dd*ssx
             v_y = dd*ssy
             v_z = dd*ssz

             dd = w(i+1, j, k, ivz)-w(i, j, k, ivz)
             w_x = dd*ssx
             w_y = dd*ssy
             w_z = dd*ssz

             dd = aa(i+1, j, k)-aa(i, j, k)
             q_x = -dd*ssx
             q_y = -dd*ssy
             q_z = -dd*ssz

             por = half*rFilv
             if(porI(i,j,k) == noFlux) por = zero

             ! Compute the laminar and (if present) the eddy viscosities
             ! multiplied by the porosity. Compute the factor in front of
             ! the gradients of the speed of sound squared for the heat
             ! flux.

             mul = por*(rlv(i,j,k) + rlv(i+1,j,k))
             mue = por*(rev(i,j,k) + rev(i+1,j,k))
             mut = mul + mue

             gm1          = half*(gamma(i,j,k) +gamma(i+1,j,k))- one
             factLamHeat  = one/(prandtl*gm1)
             factTurbHeat = one/(prandtlTurb*gm1)

             heatCoef = mul*factLamHeat + mue*factTurbHeat

             ! Compute the stress tensor and the heat flux vector.

             fracDiv = twoThird*(u_x + v_y + w_z)

             tauxx = mut*(two*u_x - fracDiv)
             tauyy = mut*(two*v_y - fracDiv)
             tauzz = mut*(two*w_z - fracDiv)

             tauxy = mut*(u_y + v_x)
             tauxz = mut*(u_z + w_x)
             tauyz = mut*(v_z + w_y)

             q_x = heatCoef*q_x
             q_y = heatCoef*q_y
             q_z = heatCoef*q_z

             ! Compute the average velocities for the face. Remember that
             ! the velocities are stored and not the momentum.

             ubar = half*(w(i,j,k,ivx) + w(i+1,j,k,ivx))
             vbar = half*(w(i,j,k,ivy) + w(i+1,j,k,ivy))
             wbar = half*(w(i,j,k,ivz) + w(i+1,j,k,ivz))

             ! Compute the viscous fluxes for this i-face.

             fmx   = tauxx*si(i,j,k,1) + tauxy*si(i,j,k,2) + tauxz*si(i,j,k,3)
             fmy   = tauxy*si(i,j,k,1) + tauyy*si(i,j,k,2) + tauyz*si(i,j,k,3)
             fmz   = tauxz*si(i,j,k,1) + tauyz*si(i,j,k,2) + tauzz*si(i,j,k,3)
             frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*si(i,j,k,1) &
                  + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*si(i,j,k,2) &
                  + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*si(i,j,k,3) &
                  - q_x*si(i,j,k,1) - q_y*si(i,j,k,2) - q_z*si(i,j,k,3)

             ! Update the residuals of cell i and i+1.
             fw(i+1,j,k,imx)   = fw(i+1,j,k,imx)   + fmx
             fw(i+1,j,k,imy)   = fw(i+1,j,k,imy)   + fmy
             fw(i+1,j,k,imz)   = fw(i+1,j,k,imz)   + fmz
             fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) + frhoE

             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
             fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE


          end do
       end do
    end do

    ! Viscous fluxes in the J-direction

    do k=2,kl
       do j=1,jl
          do i=2,il

             ! Compute the vector from the center of cell j to cell j+1
             ssx = eighth*(x(i-1,j+1,k-1,1) - x(i-1,j-1,k-1,1) &
                  +         x(i-1,j+1,k,  1) - x(i-1,j-1,k,  1) &
                  +         x(i,  j+1,k-1,1) - x(i,  j-1,k-1,1) &
                  +         x(i,  j+1,k,  1) - x(i,  j-1,k,  1))
             ssy = eighth*(x(i-1,j+1,k-1,2) - x(i-1,j-1,k-1,2) &
                  +         x(i-1,j+1,k,  2) - x(i-1,j-1,k,  2) &
                  +         x(i,  j+1,k-1,2) - x(i,  j-1,k-1,2) &
                  +         x(i,  j+1,k,  2) - x(i,  j-1,k,  2))
             ssz = eighth*(x(i-1,j+1,k-1,3) - x(i-1,j-1,k-1,3) &
                  +         x(i-1,j+1,k,  3) - x(i-1,j-1,k,  3) &
                  +         x(i,  j+1,k-1,3) - x(i,  j-1,k-1,3) &
                  +         x(i,  j+1,k,  3) - x(i,  j-1,k,  3))

             ! And determine one/ length of vector squared
             ss  = one/(ssx*ssx + ssy*ssy + ssz*ssz)
             ssx = ss*ssx
             ssy = ss*ssy
             ssz = ss*ssz

             ! Now compute each gradient
             dd = w(i, j+1, k, ivx)-w(i, j, k, ivx)
             u_x = dd*ssx
             u_y = dd*ssy
             u_z = dd*ssz

             dd = w(i, j+1, k, ivy)-w(i, j, k, ivy)
             v_x = dd*ssx
             v_y = dd*ssy
             v_z = dd*ssz

             dd = w(i, j+1, k, ivz)-w(i, j, k, ivz)
             w_x = dd*ssx
             w_y = dd*ssy
             w_z = dd*ssz

             dd = aa(i, j+1, k)-aa(i, j, k)
             q_x = -dd*ssx
             q_y = -dd*ssy
             q_z = -dd*ssz

             por = half*rFilv
             if(porJ(i,j,k) == noFlux) por = zero

             ! Compute the laminar and (if present) the eddy viscosities
             ! multiplied by the porosity. Compute the factor in front of
             ! the gradients of the speed of sound squared for the heat
             ! flux.

             mul = por*(rlv(i,j,k) + rlv(i,j+1,k))
             mue = por*(rev(i,j,k) + rev(i,j+1,k))
             mut = mul + mue

             gm1          = half*(gamma(i,j,k) + gamma(i,j+1,k)) - one
             factLamHeat  = one/(prandtl*gm1)
             factTurbHeat = one/(prandtlTurb*gm1)

             heatCoef = mul*factLamHeat + mue*factTurbHeat

             ! Compute the stress tensor and the heat flux vector.

             fracDiv = twoThird*(u_x + v_y + w_z)

             tauxx = mut*(two*u_x - fracDiv)
             tauyy = mut*(two*v_y - fracDiv)
             tauzz = mut*(two*w_z - fracDiv)

             tauxy = mut*(u_y + v_x)
             tauxz = mut*(u_z + w_x)
             tauyz = mut*(v_z + w_y)

             q_x = heatCoef*q_x
             q_y = heatCoef*q_y
             q_z = heatCoef*q_z

             ! Compute the average velocities for the face. Remember that
             ! the velocities are stored and not the momentum.

             ubar = half*(w(i,j,k,ivx) + w(i,j+1,k,ivx))
             vbar = half*(w(i,j,k,ivy) + w(i,j+1,k,ivy))
             wbar = half*(w(i,j,k,ivz) + w(i,j+1,k,ivz))

             ! Compute the viscous fluxes for this j-face.

             fmx   = tauxx*sj(i,j,k,1) + tauxy*sj(i,j,k,2) + tauxz*sj(i,j,k,3)
             fmy   = tauxy*sj(i,j,k,1) + tauyy*sj(i,j,k,2) + tauyz*sj(i,j,k,3)
             fmz   = tauxz*sj(i,j,k,1) + tauyz*sj(i,j,k,2) + tauzz*sj(i,j,k,3)
             frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sj(i,j,k,1) &
                  + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sj(i,j,k,2) &
                  + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sj(i,j,k,3) &
                  - q_x*sj(i,j,k,1) - q_y*sj(i,j,k,2) - q_z*sj(i,j,k,3)

             ! Update the residuals of cell j and j+1.

             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
             fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

             fw(i,j+1,k,imx)   = fw(i,j+1,k,imx)   + fmx
             fw(i,j+1,k,imy)   = fw(i,j+1,k,imy)   + fmy
             fw(i,j+1,k,imz)   = fw(i,j+1,k,imz)   + fmz
             fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) + frhoE

          end do
       end do
    end do

    ! Viscous fluxes in the K-direction

    do k=1,kl
       do j=2,jl
          do i=2,il

             ! Compute the vector from the center of cell k to cell k+1
             ssx = eighth*(x(i-1,j-1,k+1,1) - x(i-1,j-1,k-1,1) &
                  +         x(i-1,j,  k+1,1) - x(i-1,j,  k-1,1) &
                  +         x(i,  j-1,k+1,1) - x(i,  j-1,k-1,1) &
                  +         x(i,  j,  k+1,1) - x(i,  j,  k-1,1))
             ssy = eighth*(x(i-1,j-1,k+1,2) - x(i-1,j-1,k-1,2) &
                  +         x(i-1,j,  k+1,2) - x(i-1,j,  k-1,2) &
                  +         x(i,  j-1,k+1,2) - x(i,  j-1,k-1,2) &
                  +         x(i,  j,  k+1,2) - x(i,  j,  k-1,2))
             ssz = eighth*(x(i-1,j-1,k+1,3) - x(i-1,j-1,k-1,3) &
                  +         x(i-1,j,  k+1,3) - x(i-1,j,  k-1,3) &
                  +         x(i,  j-1,k+1,3) - x(i,  j-1,k-1,3) &
                  +         x(i,  j,  k+1,3) - x(i,  j,  k-1,3))
             ! And determine one/ length of vector squared
             ss  = one/(ssx*ssx + ssy*ssy + ssz*ssz)
             ssx = ss*ssx
             ssy = ss*ssy
             ssz = ss*ssz

             ! Now compute each gradient
             dd = w(i, j, k+1, ivx)-w(i, j, k, ivx)
             u_x = dd*ssx
             u_y = dd*ssy
             u_z = dd*ssz

             dd = w(i, j, k+1, ivy)-w(i, j, k, ivy)
             v_x = dd*ssx
             v_y = dd*ssy
             v_z = dd*ssz

             dd = w(i, j, k+1, ivz)-w(i, j, k, ivz)
             w_x = dd*ssx
             w_y = dd*ssy
             w_z = dd*ssz

             dd = aa(i, j, k+1)-aa(i, j, k)
             q_x = -dd*ssx
             q_y = -dd*ssy
             q_z = -dd*ssz

             por = half*rFilv
             if(porK(i,j,k) == noFlux) por = zero

             ! Compute the laminar and (if present) the eddy viscosities
             ! multiplied by the porosity. Compute the factor in front of
             ! the gradients of the speed of sound squared for the heat
             ! flux.

             mul = por*(rlv(i,j,k) + rlv(i,j,k+1))
             mue = por*(rev(i,j,k) + rev(i,j,k+1))
             mut = mul + mue

             gm1          = half*(gamma(i,j,k) + gamma(i,j,k+1)) - one
             factLamHeat  = one/(prandtl*gm1)
             factTurbHeat = one/(prandtlTurb*gm1)

             heatCoef = mul*factLamHeat + mue*factTurbHeat

             ! Compute the stress tensor and the heat flux vector.

             fracDiv = twoThird*(u_x + v_y + w_z)

             tauxx = mut*(two*u_x - fracDiv)
             tauyy = mut*(two*v_y - fracDiv)
             tauzz = mut*(two*w_z - fracDiv)

             tauxy = mut*(u_y + v_x)
             tauxz = mut*(u_z + w_x)
             tauyz = mut*(v_z + w_y)

             q_x = heatCoef*q_x
             q_y = heatCoef*q_y
             q_z = heatCoef*q_z

             ! Compute the average velocities for the face. Remember that
             ! the velocities are stored and not the momentum.

             ubar = half*(w(i,j,k,ivx) + w(i,j,k+1,ivx))
             vbar = half*(w(i,j,k,ivy) + w(i,j,k+1,ivy))
             wbar = half*(w(i,j,k,ivz) + w(i,j,k+1,ivz))

             ! Compute the viscous fluxes for this j-face.

             fmx   = tauxx*sk(i,j,k,1) + tauxy*sk(i,j,k,2) + tauxz*sk(i,j,k,3)
             fmy   = tauxy*sk(i,j,k,1) + tauyy*sk(i,j,k,2) + tauyz*sk(i,j,k,3)
             fmz   = tauxz*sk(i,j,k,1) + tauyz*sk(i,j,k,2) + tauzz*sk(i,j,k,3)
             frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*sk(i,j,k,1) &
                  + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*sk(i,j,k,2) &
                  + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*sk(i,j,k,3) &
                  - q_x*sk(i,j,k,1) - q_y*sk(i,j,k,2) - q_z*sk(i,j,k,3)

             ! Update the residuals of cell j and j+1.

             fw(i,j,k,imx)   = fw(i,j,k,imx)   - fmx
             fw(i,j,k,imy)   = fw(i,j,k,imy)   - fmy
             fw(i,j,k,imz)   = fw(i,j,k,imz)   - fmz
             fw(i,j,k,irhoE) = fw(i,j,k,irhoE) - frhoE

             fw(i,j,k+1,imx)   = fw(i,j,k+1,imx)   + fmx
             fw(i,j,k+1,imy)   = fw(i,j,k+1,imy)   + fmy
             fw(i,j,k+1,imz)   = fw(i,j,k+1,imz)   + fmz
             fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) + frhoE

          end do
       end do
    end do
  end subroutine viscousFluxApprox

  subroutine sumDwandFw

    ! ---------------------------------------------
    !                 Sum dw and fw/res scale
    ! ---------------------------------------------
    use constants
    use flowVarRefState, only : nw, nwf, nt1, nt2
    use inputIteration, only : turbResScale
    implicit none

    ! Variables for final summing
    integer(kind=intType) :: nTurb, i, j, k, l
    real(kind=realType) :: oVol, rBlank

    nTurb = nt2-nt1+1
    do l=1,nwf
       do k=2, kl
          do j=2, jl
             do i=2, il
                rblank = max(real(iblank(i,j,k), realType), zero)
                dw(i, j, k, l) = (dw(i, j, k, l) + fw(i, j, k, l))*rBlank
             end do
          end do
       end do
    end do
  end subroutine sumDwandFw
end module blockette
