!
!      ******************************************************************
!      *                                                                *
!      * File:          kwSolve.f90                                     *
!      * Author:        Georgi Kalitzin, Edwin van der Weide,           *
!      *                Steve Repsher (blanking)                        *
!      * Starting date: 06-27-2003                                      *
!      * Last modified: 07-05-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine kwSolve(resOnly)
!
!      ******************************************************************
!      *                                                                *
!      * kwSolve solves the k-omega transport equations of both         *
!      * the original and modified k-omega models                       *
!      * in a coupled manner using a diagonal dominant ADI-scheme.      *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use constants
       use flowVarRefState
       use inputIteration
       use inputPhysics
       use paramTurb
       use turbMod
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: resOnly
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, nn

       real(kind=realType) :: rkwGam1
       real(kind=realType) :: rhoi, ss, spk, sdk
       real(kind=realType) :: voli, volmi, volpi
       real(kind=realType) :: xm, ym, zm, xp, yp, zp, xa, ya, za
       real(kind=realType) :: ttm, ttp, mulm, mulp, muem, muep
       real(kind=realType) :: c1m, c1p, c10, c2m, c2p, c20
       real(kind=realType) :: b1, b2, c1, c2, d1, d2
       real(kind=realType) :: qs, uu, um, up, factor, utau, rblank

       real(kind=realType), dimension(itu1:itu2) :: tup

       real(kind=realType), dimension(2:il,2:jl,2:kl,2,2)  :: qq
       real(kind=realType), dimension(2,2:max(il,jl,kl))   :: bb, dd, ff
       real(kind=realType), dimension(2,2,2:max(il,jl,kl)) :: cc

       real(kind=realType), dimension(:,:,:), pointer :: ddw, ww, ddvt
       real(kind=realType), dimension(:,:),   pointer :: rrlv
       real(kind=realType), dimension(:,:),   pointer :: dd2Wall

       logical, dimension(2:jl,2:kl), target :: flagI2, flagIl
       logical, dimension(2:il,2:kl), target :: flagJ2, flagJl
       logical, dimension(2:il,2:jl), target :: flagK2, flagKl

       logical, dimension(:,:), pointer :: flag
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set model constants

       rkwGam1  = rkwBeta1/rkwBetas - rkwSigw1*rkwK*rkwK/sqrt(rkwBetas)
       sig1     = rkwSigk1
       sig2     = rkwSigw1

       ! Set the pointer for dvt in dw, such that the code is more
       ! readable. Also set the pointers for the production term,
       ! vorticity and the cross diffusion term.

       dvt  => dw(1:,1:,1:,idvt:)
       prod => dw(1:,1:,1:,iprod)
       vort => prod
       kwCD => dw(1:,1:,1:,icd)
!
!      ******************************************************************
!      *                                                                *
!      * Production term.                                               *
!      *                                                                *
!      ******************************************************************
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
!      ******************************************************************
!      *                                                                *
!      * Source terms.                                                  *
!      *                                                                *
!      * Determine the source term and its derivative w.r.t. k and      *
!      * omega for all internal cells of the block.                     *
!      *                                                                *
!      ******************************************************************
!
       do k=2,kl
         do j=2,jl
           do i=2,il

             ! Compute the source terms for both the k and the omega
             ! equation. Note that prod contains the unscaled production
             ! term. Furthermore the production term of k is limited to
             ! a certain times the destruction term.

             rhoi = one/w(i,j,k,irho)
             ss   = prod(i,j,k)
             spk  = rev(i,j,k)*ss*rhoi
             sdk  = rkwBetas*w(i,j,k,itu1)*w(i,j,k,itu2)
             spk  = min(spk, pklim*sdk)

             dvt(i,j,k,1) = spk - sdk
             dvt(i,j,k,2) = rkwGam1*ss - rkwBeta1*w(i,j,k,itu2)**2

             ! Compute the source term jacobian. Note that only the
             ! destruction terms are linearized to increase the diagonal
             ! dominance of the matrix. Furthermore minus the source
             ! term jacobian is stored.

             qq(i,j,k,1,1) = rkwBetas*w(i,j,k,itu2)
           ! qq(i,j,k,1,2) = rkwBetas*w(i,j,k,itu1)
             qq(i,j,k,1,2) = zero
             qq(i,j,k,2,1) = zero
             qq(i,j,k,2,2) = two*rkwBeta1*w(i,j,k,itu2)

           enddo
         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Compute the cross-diffusion term for the modified version of   *
!      * the k-omega model. It should cure the free-stream dependency   *
!      * of the original model.                                         *
!      *                                                                *
!      ******************************************************************
!
       if(turbModel == komegaModified) then

         ! Compute the cross diffusion term.

         call kwCDterm

         ! Multiply it with the correct constant and take it only into
         ! account if it is positive.

         do k=2,kl
           do j=2,jl
             do i=2,il
               dvt(i,j,k,2) = dvt(i,j,k,2) &
                            + rkwSigd1*max(zero, kwCD(i,j,k))
             enddo
           enddo
         enddo

       endif
!
!      ******************************************************************
!      *                                                                *
!      * Advection and unsteady terms.                                  *
!      *                                                                *
!      ******************************************************************
!
       nn = itu1 - 1
       call turbAdvection(2_intType, 2_intType, nn, qq)

       call unsteadyTurbTerm(2_intType, 2_intType, nn, qq)
!
!      ******************************************************************
!      *                                                                *
!      * Viscous terms in k-direction.                                  *
!      *                                                                *
!      ******************************************************************
!
       do k=2,kl
         do j=2,jl
           do i=2,il

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
             ! The second derivative in zeta-direction is constructed as
             ! the central difference of the first order derivatives, i.e.
             ! d^2/dzeta^2 = d/dzeta (d/dzeta k+1/2 - d/dzeta k-1/2).
             ! In this way the metric as well as the varying viscosity
             ! can be taken into account; the latter appears inside the
             ! d/dzeta derivative. The whole term is divided by rho to
             ! obtain the diffusion term for k and omega.

             rhoi = one/w(i,j,k,irho)
             mulm = half*(rlv(i,j,k-1) + rlv(i,j,k))
             mulp = half*(rlv(i,j,k+1) + rlv(i,j,k))
             muem = half*(rev(i,j,k-1) + rev(i,j,k))
             muep = half*(rev(i,j,k+1) + rev(i,j,k))

             c1m = ttm*(mulm + sig1*muem)*rhoi
             c1p = ttp*(mulp + sig1*muep)*rhoi
             c10 = c1m + c1p

             c2m = ttm*(mulm + sig2*muem)*rhoi
             c2p = ttp*(mulp + sig2*muep)*rhoi
             c20 = c2m + c2p

             ! Update the residual for this cell and store the possible
             ! coefficients for the matrix in b1, b2, c1, c2, d1 and d2.

             dvt(i,j,k,1) = dvt(i,j,k,1)      + c1m*w(i,j,k-1,itu1) &
                          - c10*w(i,j,k,itu1) + c1p*w(i,j,k+1,itu1)
             dvt(i,j,k,2) = dvt(i,j,k,2)      + c2m*w(i,j,k-1,itu2) &
                          - c20*w(i,j,k,itu2) + c2p*w(i,j,k+1,itu2)

             b1 = -c1m
             c1 =  c10
             d1 = -c1p

             b2 = -c2m
             c2 =  c20
             d2 = -c2p

             ! Update the central jacobian. For nonboundary cells this
             ! is simply c1 and c2. For boundary cells this is slightly
             ! more complicated, because the boundary conditions are
             ! treated implicitly and the off-diagonal terms b1, b2 and
             ! d1, d2 must be taken into account.
             ! The boundary conditions are only treated implicitly if
             ! the diagonal dominance of the matrix is increased.

             if(k == 2) then
               qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                             - b1*max(bmtk1(i,j,itu1,itu1),zero)
               qq(i,j,k,1,2) = qq(i,j,k,1,2) - b1*bmtk1(i,j,itu1,itu2)
               qq(i,j,k,2,1) = qq(i,j,k,2,1) - b2*bmtk1(i,j,itu2,itu1)
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                             - b2*max(bmtk1(i,j,itu2,itu2),zero)
             else if(k == kl) then
               qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                             - d1*max(bmtk2(i,j,itu1,itu1),zero)
               qq(i,j,k,1,2) = qq(i,j,k,1,2) - d1*bmtk2(i,j,itu1,itu2)
               qq(i,j,k,2,1) = qq(i,j,k,2,1) - d2*bmtk2(i,j,itu2,itu1)
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                             - d2*max(bmtk2(i,j,itu2,itu2),zero)
             else
               qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2
             endif

           enddo
         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Viscous terms in j-direction.                                  *
!      *                                                                *
!      ******************************************************************
!
       do k=2,kl
         do j=2,jl
           do i=2,il

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
             ! The second derivative in eta-direction is constructed as
             ! the central difference of the first order derivatives, i.e.
             ! d^2/deta^2 = d/deta (d/deta j+1/2 - d/deta j-1/2).
             ! In this way the metric as well as the varying viscosity
             ! can be taken into account; the latter appears inside the
             ! d/deta derivative. The whole term is divided by rho to
             ! obtain the diffusion term for k and omega.

             rhoi = one/w(i,j,k,irho)
             mulm = half*(rlv(i,j-1,k) + rlv(i,j,k))
             mulp = half*(rlv(i,j+1,k) + rlv(i,j,k))
             muem = half*(rev(i,j-1,k) + rev(i,j,k))
             muep = half*(rev(i,j+1,k) + rev(i,j,k))

             c1m = ttm*(mulm + sig1*muem)*rhoi
             c1p = ttp*(mulp + sig1*muep)*rhoi
             c10 = c1m + c1p

             c2m = ttm*(mulm + sig2*muem)*rhoi
             c2p = ttp*(mulp + sig2*muep)*rhoi
             c20 = c2m + c2p

             ! Update the residual for this cell and store the possible
             ! coefficients for the matrix in b1, b2, c1, c2, d1 and d2.

             dvt(i,j,k,1) = dvt(i,j,k,1)      + c1m*w(i,j-1,k,itu1) &
                          - c10*w(i,j,k,itu1) + c1p*w(i,j+1,k,itu1)
             dvt(i,j,k,2) = dvt(i,j,k,2)      + c2m*w(i,j-1,k,itu2) &
                          - c20*w(i,j,k,itu2) + c2p*w(i,j+1,k,itu2)

             b1 = -c1m
             c1 =  c10
             d1 = -c1p

             b2 = -c2m
             c2 =  c20
             d2 = -c2p

             ! Update the central jacobian. For nonboundary cells this
             ! is simply c1 and c2. For boundary cells this is slightly
             ! more complicated, because the boundary conditions are
             ! treated implicitly and the off-diagonal terms b1, b2 and
             ! d1, d2 must be taken into account.
             ! The boundary conditions are only treated implicitly if
             ! the diagonal dominance of the matrix is increased.

             if(j == 2) then
               qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                             - b1*max(bmtj1(i,k,itu1,itu1),zero)
               qq(i,j,k,1,2) = qq(i,j,k,1,2) - b1*bmtj1(i,k,itu1,itu2)
               qq(i,j,k,2,1) = qq(i,j,k,2,1) - b2*bmtj1(i,k,itu2,itu1)
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                             - b2*max(bmtj1(i,k,itu2,itu2),zero)
             else if(j == jl) then
               qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                             - d1*max(bmtj2(i,k,itu1,itu1),zero)
               qq(i,j,k,1,2) = qq(i,j,k,1,2) - d1*bmtj2(i,k,itu1,itu2)
               qq(i,j,k,2,1) = qq(i,j,k,2,1) - d2*bmtj2(i,k,itu2,itu1)
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                             - d2*max(bmtj2(i,k,itu2,itu2),zero)
             else
               qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2
             endif

           enddo
         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Viscous terms in i-direction.                                  *
!      *                                                                *
!      ******************************************************************
!
       do k=2,kl
         do j=2,jl
           do i=2,il

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
             ! The second derivative in xi-direction is constructed as
             ! the central difference of the first order derivatives, i.e.
             ! d^2/dxi^2 = d/dxi (d/dxi i+1/2 - d/dxi i-1/2).
             ! In this way the metric as well as the varying viscosity
             ! can be taken into account; the latter appears inside the
             ! d/dxi derivative. The whole term is divided by rho to
             ! obtain the diffusion term for k and omega.

             rhoi = one/w(i,j,k,irho)
             mulm = half*(rlv(i-1,j,k) + rlv(i,j,k))
             mulp = half*(rlv(i+1,j,k) + rlv(i,j,k))
             muem = half*(rev(i-1,j,k) + rev(i,j,k))
             muep = half*(rev(i+1,j,k) + rev(i,j,k))

             c1m = ttm*(mulm + sig1*muem)*rhoi
             c1p = ttp*(mulp + sig1*muep)*rhoi
             c10 = c1m + c1p

             c2m = ttm*(mulm + sig2*muem)*rhoi
             c2p = ttp*(mulp + sig2*muep)*rhoi
             c20 = c2m + c2p

             ! Update the residual for this cell and store the possible
             ! coefficients for the matrix in b1, b2, c1, c2, d1 and d2.

             dvt(i,j,k,1) = dvt(i,j,k,1)      + c1m*w(i-1,j,k,itu1) &
                          - c10*w(i,j,k,itu1) + c1p*w(i+1,j,k,itu1)
             dvt(i,j,k,2) = dvt(i,j,k,2)      + c2m*w(i-1,j,k,itu2) &
                          - c20*w(i,j,k,itu2) + c2p*w(i+1,j,k,itu2)

             b1 = -c1m
             c1 =  c10
             d1 = -c1p

             b2 = -c2m
             c2 =  c20
             d2 = -c2p

             ! Update the central jacobian. For nonboundary cells this
             ! is simply c1 and c2. For boundary cells this is slightly
             ! more complicated, because the boundary conditions are
             ! treated implicitly and the off-diagonal terms b1, b2 and
             ! d1, d2 must be taken into account.
             ! The boundary conditions are only treated implicitly if
             ! the diagonal dominance of the matrix is increased.

             if(i == 2) then
               qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                             - b1*max(bmti1(j,k,itu1,itu1),zero)
               qq(i,j,k,1,2) = qq(i,j,k,1,2) - b1*bmti1(j,k,itu1,itu2)
               qq(i,j,k,2,1) = qq(i,j,k,2,1) - b2*bmti1(j,k,itu2,itu1)
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                             - b2*max(bmti1(j,k,itu2,itu2),zero)
             else if(i == il) then
               qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                             - d1*max(bmti2(j,k,itu1,itu1),zero)
               qq(i,j,k,1,2) = qq(i,j,k,1,2) - d1*bmti2(j,k,itu1,itu2)
               qq(i,j,k,2,1) = qq(i,j,k,2,1) - d2*bmti2(j,k,itu2,itu1)
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                             - d2*max(bmti2(j,k,itu2,itu2),zero)
             else
               qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2
             endif

           enddo
         enddo
       enddo

       ! Multiply the residual by the volume and store this in dw; this
       ! is done for monitoring reasons only. The multiplication with the
       ! volume is present to be consistent with the flow residuals; also
       ! the negative value is taken, again to be consistent with the
       ! flow equations. Also multiply by iblank so that no updates occur
       ! in holes or the overset boundary.

       do k=2,kl
         do j=2,jl
           do i=2,il
             rblank = real(iblank(i,j,k), realType)
             dw(i,j,k,itu1) = -vol(i,j,k)*dvt(i,j,k,1)*rblank
             dw(i,j,k,itu2) = -vol(i,j,k)*dvt(i,j,k,2)*rblank
           enddo
         enddo
       enddo

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
               ddw     => dw(2,1:,1:,1:);  ddvt => dvt(2,1:,1:,1:)
               ww      => w(2,1:,1:,1:);   rrlv => rlv(2,1:,1:)
               dd2Wall => d2Wall(2,:,:)

             case (iMax)
               flag    => flagIl
               ddw     => dw(il,1:,1:,1:);  ddvt => dvt(il,1:,1:,1:)
               ww      => w(il,1:,1:,1:);   rrlv => rlv(il,1:,1:)
               dd2Wall => d2Wall(il,:,:)

             case (jMin)
               flag    => flagJ2
               ddw     => dw(1:,2,1:,1:);  ddvt => dvt(1:,2,1:,1:)
               ww      => w(1:,2,1:,1:);   rrlv => rlv(1:,2,1:)
               dd2Wall => d2Wall(:,2,:)

             case (jMax)
               flag    => flagJl
               ddw     => dw(1:,jl,1:,1:);  ddvt => dvt(1:,jl,1:,1:)
               ww      => w(1:,jl,1:,1:);   rrlv => rlv(1:,jl,1:)
               dd2Wall => d2Wall(:,jl,:)

             case (kMin)
               flag    => flagK2
               ddw     => dw(1:,1:,2,1:);  ddvt => dvt(1:,1:,2,1:)
               ww      => w(1:,1:,2,1:);   rrlv => rlv(1:,1:,2)
               dd2Wall => d2Wall(:,:,2)

             case (kMax)
               flag    => flagKl
               ddw     => dw(1:,1:,kl,1:);  ddvt => dvt(1:,1:,kl,1:)
               ww      => w(1:,1:,kl,1:);   rrlv => rlv(1:,1:,kl)
               dd2Wall => d2Wall(:,:,kl)

           end select

           ! Loop over the owned faces of this subface. Therefore the
           ! nodal range of BCData must be used. The offset of +1 is
           ! present, because the starting index of the cell range is
           ! 1 larger than the starting index of the nodal range.

           do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
             do i=(BCData(nn)%inBeg+1),BCData(nn)%inEnd

               ! Set ddw to zero for proper monitoring of the
               ! convergence.

               ddw(i,j,itu1) = zero
               ddw(i,j,itu2) = zero

               ! Enforce k and omega in the 1st internal cell from
               ! the wall function table. There is an offset of -1 in
               ! the wall distance. Note that the offset compared to
               ! the current value must be stored. Also note that the
               ! curve fits contain the non-dimensional values.

               utau = viscSubface(nn)%utau(i,j)
               yp = ww(i,j,irho)*dd2Wall(i-1,j-1)*utau/rrlv(i,j)

               call curveTupYp(tup, yp, itu1, itu2)

               tup(itu1) = tup(itu1)*utau**2
               tup(itu2) = tup(itu2)*utau**2/rrlv(i,j)*ww(i,j,irho)

               ddvt(i,j,1) = tup(itu1) - ww(i,j,itu1)
               ddvt(i,j,2) = tup(itu2) - ww(i,j,itu2)

               ! Set the wall flag to .true.

               flag(i,j) = .true.

             enddo
           enddo

         enddo bocos
       endif testWallFunctions

       ! Return if only the residual must be computed.

       if( resOnly ) return

       ! For implicit relaxation take the local time step into account,
       ! where dt is the inverse of the central jacobian times the cfl
       ! number. The following system is solved:
       ! (i/dt + cc + bb + dd)*dw = rhs, in which i/dt = cc/cfl. As in
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
             qq(i,j,k,1,1) = factor*qq(i,j,k,1,1)
             qq(i,j,k,1,2) = factor*qq(i,j,k,1,2)
             qq(i,j,k,2,1) = factor*qq(i,j,k,2,1)
             qq(i,j,k,2,2) = factor*qq(i,j,k,2,2)

             ! Set qq to 1 if the value is determined by the table.

             if((i == 2  .and. flagI2(j,k)) .or. &
                (i == il .and. flagIl(j,k)) .or. &
                (j == 2  .and. flagJ2(i,k)) .or. &
                (j == jl .and. flagJl(i,k)) .or. &
                (k == 2  .and. flagK2(i,j)) .or. &
                (k == kl .and. flagKl(i,j))) then
                    qq(i,j,k,1,1) = one
                    qq(i,j,k,1,2) = zero
                    qq(i,j,k,2,1) = zero
                    qq(i,j,k,2,2) = one
             endif
           enddo
         enddo
       enddo

       ! Initialize the grid velocity to zero. This value will be used
       ! if the block is not moving.

       qs = zero
!
!      ******************************************************************
!      *                                                                *
!      * dd-ADI step in j-direction. There is no particular reason to   *
!      * start in j-direction, it just happened to be so. As we solve   *
!      * in j-direction, the j-loop is the innermost loop.              *
!      *                                                                *
!      ******************************************************************
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

             ! Off-diagonal terms due to the diffusion terms
             ! in j-direction.

             rhoi = one/w(i,j,k,irho)
             mulm = half*(rlv(i,j-1,k) + rlv(i,j,k))
             mulp = half*(rlv(i,j+1,k) + rlv(i,j,k))
             muem = half*(rev(i,j-1,k) + rev(i,j,k))
             muep = half*(rev(i,j+1,k) + rev(i,j,k))

             c1m = ttm*(mulm + sig1*muem)*rhoi
             c1p = ttp*(mulp + sig1*muep)*rhoi

             c2m = ttm*(mulm + sig2*muem)*rhoi
             c2p = ttp*(mulp + sig2*muep)*rhoi

             bb(1,j) = -c1m
             dd(1,j) = -c1p
             bb(2,j) = -c2m
             dd(2,j) = -c2p

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

             bb(1,j) = bb(1,j) - up
             dd(1,j) = dd(1,j) + um
             bb(2,j) = bb(2,j) - up
             dd(2,j) = dd(2,j) + um

             ! Store the central jacobian and rhs in cc and ff.
             ! Multiply the off-diagonal terms and rhs by the iblank
             ! value so the update determined for iblank = 0 is zero.

             rblank = real(iblank(i,j,k), realType)

             cc(1,1,j) = qq(i,j,k,1,1)
             cc(1,2,j) = qq(i,j,k,1,2)*rblank
             cc(2,1,j) = qq(i,j,k,2,1)*rblank
             cc(2,2,j) = qq(i,j,k,2,2)

             ff(1,j) = dvt(i,j,k,1)*rblank
             ff(2,j) = dvt(i,j,k,2)*rblank

             bb(:,j) = bb(:,j)*rblank
             dd(:,j) = dd(:,j)*rblank

             ! Set off diagonal terms to zero if wall function are used.

             if((i == 2  .and. flagI2(j,k)) .or. &
                (i == il .and. flagIl(j,k)) .or. &
                (j == 2  .and. flagJ2(i,k)) .or. &
                (j == jl .and. flagJl(i,k)) .or. &
                (k == 2  .and. flagK2(i,j)) .or. &
                (k == kl .and. flagKl(i,j))) then
                  bb(1,j) = zero
                  dd(1,j) = zero
                  bb(2,j) = zero
                  dd(2,j) = zero
             endif

           enddo

           ! Solve the tri-diagonal system in j-direction.

           call tdia3(2_intType, jl, bb, cc, dd, ff)

           ! Determine the new rhs for the next direction.

           do j=2,jl
             dvt(i,j,k,1) = qq(i,j,k,1,1)*ff(1,j) + qq(i,j,k,1,2)*ff(2,j)
             dvt(i,j,k,2) = qq(i,j,k,2,1)*ff(1,j) + qq(i,j,k,2,2)*ff(2,j)
           enddo

         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * dd-ADI step in i-direction. As we solve in i-direction, the    *
!      * i-loop is the innermost loop.                                  *
!      *                                                                *
!      ******************************************************************
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

             ! Off-diagonal terms due to the diffusion terms
             ! in i-direction.

             rhoi = one/w(i,j,k,irho)
             mulm = half*(rlv(i-1,j,k) + rlv(i,j,k))
             mulp = half*(rlv(i+1,j,k) + rlv(i,j,k))
             muem = half*(rev(i-1,j,k) + rev(i,j,k))
             muep = half*(rev(i+1,j,k) + rev(i,j,k))

             c1m = ttm*(mulm + sig1*muem)*rhoi
             c1p = ttp*(mulp + sig1*muep)*rhoi

             c2m = ttm*(mulm + sig2*muem)*rhoi
             c2p = ttp*(mulp + sig2*muep)*rhoi

             bb(1,i) = -c1m
             dd(1,i) = -c1p
             bb(2,i) = -c2m
             dd(2,i) = -c2p

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

             bb(1,i) = bb(1,i) - up
             dd(1,i) = dd(1,i) + um
             bb(2,i) = bb(2,i) - up
             dd(2,i) = dd(2,i) + um

             ! Store the central jacobian and rhs in cc and ff.
             ! Multiply the off-diagonal terms and rhs by the iblank
             ! value so the update determined for iblank = 0 is zero.

             rblank = real(iblank(i,j,k), realType)

             cc(1,1,i) = qq(i,j,k,1,1)
             cc(1,2,i) = qq(i,j,k,1,2)*rblank
             cc(2,1,i) = qq(i,j,k,2,1)*rblank
             cc(2,2,i) = qq(i,j,k,2,2)

             ff(1,i) = dvt(i,j,k,1)*rblank
             ff(2,i) = dvt(i,j,k,2)*rblank

             bb(:,i) = bb(:,i)*rblank
             dd(:,i) = dd(:,i)*rblank

             ! Set off diagonal terms to zero if wall function are used.

             if((i == 2  .and. flagI2(j,k)) .or. &
                (i == il .and. flagIl(j,k)) .or. &
                (j == 2  .and. flagJ2(i,k)) .or. &
                (j == jl .and. flagJl(i,k)) .or. &
                (k == 2  .and. flagK2(i,j)) .or. &
                (k == kl .and. flagKl(i,j))) then
                  bb(1,i) = zero
                  dd(1,i) = zero
                  bb(2,i) = zero
                  dd(2,i) = zero
             endif

           enddo

           ! Solve the tri-diagonal system in i-direction.

           call tdia3(2_intType, il, bb, cc, dd, ff)

           ! Determine the new rhs for the next direction.

           do i=2,il
             dvt(i,j,k,1) = qq(i,j,k,1,1)*ff(1,i) + qq(i,j,k,1,2)*ff(2,i)
             dvt(i,j,k,2) = qq(i,j,k,2,1)*ff(1,i) + qq(i,j,k,2,2)*ff(2,i)
           enddo

         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * dd-adi step in k-direction. As we solve in k-direction, the    *
!      * k-loop is the innermost loop.                                  *
!      *                                                                *
!      ******************************************************************
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

             ! Off-diagonal terms due to the diffusion terms
             ! in k-direction.

             rhoi = one/w(i,j,k,irho)
             mulm = half*(rlv(i,j,k-1) + rlv(i,j,k))
             mulp = half*(rlv(i,j,k+1) + rlv(i,j,k))
             muem = half*(rev(i,j,k-1) + rev(i,j,k))
             muep = half*(rev(i,j,k+1) + rev(i,j,k))

             c1m = ttm*(mulm + sig1*muem)*rhoi
             c1p = ttp*(mulp + sig1*muep)*rhoi

             c2m = ttm*(mulm + sig2*muem)*rhoi
             c2p = ttp*(mulp + sig2*muep)*rhoi

             bb(1,k) = -c1m
             dd(1,k) = -c1p
             bb(2,k) = -c2m
             dd(2,k) = -c2p

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

             bb(1,k) = bb(1,k) - up
             dd(1,k) = dd(1,k) + um
             bb(2,k) = bb(2,k) - up
             dd(2,k) = dd(2,k) + um

             ! Store the central jacobian and rhs in cc and ff.
             ! Multiply the off-diagonal terms and rhs by the iblank
             ! value so the update determined for iblank = 0 is zero.

             rblank = real(iblank(i,j,k), realType)

             cc(1,1,k) = qq(i,j,k,1,1)
             cc(1,2,k) = qq(i,j,k,1,2)*rblank
             cc(2,1,k) = qq(i,j,k,2,1)*rblank
             cc(2,2,k) = qq(i,j,k,2,2)

             ff(1,k) = dvt(i,j,k,1)*rblank
             ff(2,k) = dvt(i,j,k,2)*rblank

             bb(:,k) = bb(:,k)*rblank
             dd(:,k) = dd(:,k)*rblank

             ! Set off diagonal terms to zero if wall function are used.

             if((i == 2  .and. flagI2(j,k)) .or. &
                (i == il .and. flagIl(j,k)) .or. &
                (j == 2  .and. flagJ2(i,k)) .or. &
                (j == jl .and. flagJl(i,k)) .or. &
                (k == 2  .and. flagK2(i,j)) .or. &
                (k == kl .and. flagKl(i,j))) then
                  bb(1,k) = zero
                  dd(1,k) = zero
                  bb(2,k) = zero
                  dd(2,k) = zero
             endif

           enddo

           ! Solve the tri-diagonal system in k-direction.

           call tdia3(2_intType, kl, bb, cc, dd, ff)

           ! Store the update in dvt.

           do k=2,kl
             dvt(i,j,k,1) = ff(1,k)
             dvt(i,j,k,2) = ff(2,k)
           enddo

         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Update the turbulent variables. For explicit relaxation the    *
!      * update must be relaxed; for implicit relaxation this has been  *
!      * done via the time step.                                        *
!      *                                                                *
!      ******************************************************************
!
       factor = one
       if(turbRelax == turbRelaxExplicit) factor = alfaTurb

       do k=2,kl
         do j=2,jl
           do i=2,il
             w(i,j,k,itu1) = w(i,j,k,itu1) + factor*dvt(i,j,k,1)
             w(i,j,k,itu1) = max(w(i,j,k,itu1), zero)

             w(i,j,k,itu2) = w(i,j,k,itu2) + factor*dvt(i,j,k,2)
             w(i,j,k,itu2) = max(w(i,j,k,itu2), 1.e-5_realType*wInf(itu2))
           enddo
         enddo
       enddo

       end subroutine kwSolve
