!
!      ******************************************************************
!      *                                                                *
!      * File:          vfSolve.F90                                     *
!      * Author:        Georgi Kalitzin, Edwin van der Weide,           *
!      *                Steve Repsher (blanking)                        *
!      * Starting date: 04-14-2004                                      *
!      * Last modified: 07-05-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine vfSolve(resOnly)
!
!      ******************************************************************
!      *                                                                *
!      * vfSolve solves the v2 transport equation and the               *
!      * f elliptic relaxation equation of the v2-f model               *
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

       real(kind=realType) :: rhoi, ss, spk, ff1, ff2, ff3, ss1, ss2, ss3
       real(kind=realType) :: voli, volmi, volpi
       real(kind=realType) :: xm, ym, zm, xp, yp, zp, xa, ya, za
       real(kind=realType) :: ttm, ttp, mulm, mulp, muem, muep
       real(kind=realType) :: c1m, c1p, c10, c2m, c2p, c20
       real(kind=realType) :: b1, b2, c1, c2, d1, d2
       real(kind=realType) :: qs, uu, um, up, utau
       real(kind=realType) :: tke, tep, tv2, tf2, tkea, tepa, tv2a
       real(kind=realType) :: tkel, tv2l, sle2i, stei
       real(kind=realType) :: rsct, rnu, rn2
       real(kind=realType) :: tu12, tu22, tu32, tu42, tu52
       real(kind=realType) :: rnu23, dtu23, rblank

       real(kind=realType), dimension(itu1:itu5) :: tup

       real(kind=realType), dimension(2:il,2:jl,2:kl,2,2)  :: qq
       real(kind=realType), dimension(2,2:max(il,jl,kl))   :: bb, dd, ff
       real(kind=realType), dimension(2,2,2:max(il,jl,kl)) :: cc

       real(kind=realType), dimension(:,:,:), pointer :: dw2, dvt2, w2, w3
       real(kind=realType), dimension(:,:),   pointer :: rlv2, rlv3
       real(kind=realType), dimension(:,:),   pointer :: rev2, rev3
       real(kind=realType), dimension(:,:),   pointer :: d2Wall2, d2Wall3

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

       sig1 = rvfSigv1
       sig2 = one
!
!      ******************************************************************
!      *                                                                *
!      * Source terms.                                                  *
!      *                                                                *
!      * Determine the source term and its derivative w.r.t. v2 and f   *
!      * for all internal cells of the block.                           *
!      *                                                                *
!      ******************************************************************
!
       do k=2,kl
         do j=2,jl
           do i=2,il

             ! Compute the source terms for both the k and the epsilon
             ! equation. Note that dw(i,j,k,iprod) contains the unscaled
             ! production term.

             rhoi = one/w(i,j,k,irho)
             tke  = w(i,j,k,itu1)
             tep  = w(i,j,k,itu2)
             tv2  = w(i,j,k,itu3)
             tf2  = w(i,j,k,itu4)
             tkea = abs(tke)
             tepa = abs(tep)
             tv2a = abs(tv2)
             tkel = max(tkea,rvfLimitK)
             tv2l = max(tv2a,rvfLimitK*0.666666666_realType)
             stei = tepa/sct(i,j,k)
             sle2i= tepa**2/scl2(i,j,k)

             if(rvfN == 6) then
               rnu  = rlv(i,j,k)*rhoi
               rsct = max(tkea,6.*sqrt(rnu*tepa))
               stei = tepa/rsct

               ! rn2  = rvfCn**2*(rnu*tepa)**1.5
               ! sle2i= rvfCl**2*max(tkea**3,rn2)
             endif

             ss   = prod(i,j,k)
             spk  = rev(i,j,k)*ss*rhoi
             spk  = min(spk, pklim*tepa)

             ff1  =-tepa/tkel
             ff2  = tkel
             ff3  = zero
 
             ss1  =-(rvfC1-1.)/tkel*stei*sle2i
             ss2  =-sle2i
             ss3  = (rvfC1-1.)*2./3.*stei*sle2i + rvfC2*spk/tkel*sle2i

             if(rvfN == 6) then
                ff1 = ff1 - 5.0*tepa/tkel
                ss1 = ss1 + 5.0/tkel*stei*sle2i
             endif

             dvt(i,j,k,1) = ff1*tv2 + ff2*tf2 + ff3
             dvt(i,j,k,2) = ss1*tv2 + ss2*tf2 + ss3

             ! Compute the source term jacobian. Note that only the
             ! destruction terms are linearized to increase the diagonal
             ! dominance of the matrix. Furthermore minus the source
             ! term jacobian is stored.

             qq(i,j,k,1,1) = -ff1
             qq(i,j,k,1,2) = -ff2
             qq(i,j,k,2,1) = -ss1
             qq(i,j,k,2,2) = -ss2

           enddo
         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Advection and unsteady terms.                                  *
!      *                                                                *
!      ******************************************************************
!
       nn = itu3 - 1
       call turbAdvection(2_intType, 1_intType, nn, qq)

       call unsteadyTurbTerm(2_intType, 1_intType, nn, qq)
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
             ! obtain the diffusion term for v2 and f.

             rhoi = one/w(i,j,k,irho)
             mulm = half*(rlv(i,j,k-1) + rlv(i,j,k))
             mulp = half*(rlv(i,j,k+1) + rlv(i,j,k))
             muem = half*(rev(i,j,k-1) + rev(i,j,k))
             muep = half*(rev(i,j,k+1) + rev(i,j,k))

             c1m = ttm*(mulm + sig1*muem)*rhoi
             c1p = ttp*(mulp + sig1*muep)*rhoi
             c10 = c1m + c1p

             c2m = ttm
             c2p = ttp
             c20 = c2m + c2p

             ! Update the residual for this cell and store the possible
             ! coefficients for the matrix in b1, b2, c1, c2, d1 and d2.

             dvt(i,j,k,1) = dvt(i,j,k,1)      + c1m*w(i,j,k-1,itu3) &
                          - c10*w(i,j,k,itu3) + c1p*w(i,j,k+1,itu3)
             dvt(i,j,k,2) = dvt(i,j,k,2)      + c2m*w(i,j,k-1,itu4) &
                          - c20*w(i,j,k,itu4) + c2p*w(i,j,k+1,itu4)

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
                             - b1*max(bmtk1(i,j,itu3,itu3),zero)
               qq(i,j,k,1,2) = qq(i,j,k,1,2) - b1*bmtk1(i,j,itu3,itu4)
               qq(i,j,k,2,1) = qq(i,j,k,2,1) - b2*bmtk1(i,j,itu4,itu3)
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                             - b2*max(bmtk1(i,j,itu4,itu4),zero)
             else if(k == kl) then
               qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                             - d1*max(bmtk2(i,j,itu3,itu3),zero)
               qq(i,j,k,1,2) = qq(i,j,k,1,2) - d1*bmtk2(i,j,itu3,itu4)
               qq(i,j,k,2,1) = qq(i,j,k,2,1) - d2*bmtk2(i,j,itu4,itu3)
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                             - d2*max(bmtk2(i,j,itu4,itu4),zero)
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
             ! obtain the diffusion term for v2 and f.

             rhoi = one/w(i,j,k,irho)
             mulm = half*(rlv(i,j-1,k) + rlv(i,j,k))
             mulp = half*(rlv(i,j+1,k) + rlv(i,j,k))
             muem = half*(rev(i,j-1,k) + rev(i,j,k))
             muep = half*(rev(i,j+1,k) + rev(i,j,k))

             c1m = ttm*(mulm + sig1*muem)*rhoi
             c1p = ttp*(mulp + sig1*muep)*rhoi
             c10 = c1m + c1p

             c2m = ttm
             c2p = ttp
             c20 = c2m + c2p

             ! Update the residual for this cell and store the possible
             ! coefficients for the matrix in b1, b2, c1, c2, d1 and d2.

             dvt(i,j,k,1) = dvt(i,j,k,1)      + c1m*w(i,j-1,k,itu3) &
                          - c10*w(i,j,k,itu3) + c1p*w(i,j+1,k,itu3)
             dvt(i,j,k,2) = dvt(i,j,k,2)      + c2m*w(i,j-1,k,itu4) &
                          - c20*w(i,j,k,itu4) + c2p*w(i,j+1,k,itu4)

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
                             - b1*max(bmtj1(i,k,itu3,itu3),zero)
               qq(i,j,k,1,2) = qq(i,j,k,1,2) - b1*bmtj1(i,k,itu3,itu4)
               qq(i,j,k,2,1) = qq(i,j,k,2,1) - b2*bmtj1(i,k,itu4,itu3)
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                             - b2*max(bmtj1(i,k,itu4,itu4),zero)
             else if(j == jl) then
               qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                             - d1*max(bmtj2(i,k,itu3,itu3),zero)
               qq(i,j,k,1,2) = qq(i,j,k,1,2) - d1*bmtj2(i,k,itu3,itu4)
               qq(i,j,k,2,1) = qq(i,j,k,2,1) - d2*bmtj2(i,k,itu4,itu3)
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                             - d2*max(bmtj2(i,k,itu4,itu4),zero)
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
             ! obtain the diffusion term for v2 and f.

             rhoi = one/w(i,j,k,irho)
             mulm = half*(rlv(i-1,j,k) + rlv(i,j,k))
             mulp = half*(rlv(i+1,j,k) + rlv(i,j,k))
             muem = half*(rev(i-1,j,k) + rev(i,j,k))
             muep = half*(rev(i+1,j,k) + rev(i,j,k))

             c1m = ttm*(mulm + sig1*muem)*rhoi
             c1p = ttp*(mulp + sig1*muep)*rhoi
             c10 = c1m + c1p

             c2m = ttm
             c2p = ttp
             c20 = c2m + c2p

             ! Update the residual for this cell and store the possible
             ! coefficients for the matrix in b1, b2, c1, c2, d1 and d2.

             dvt(i,j,k,1) = dvt(i,j,k,1)      + c1m*w(i-1,j,k,itu3) &
                          - c10*w(i,j,k,itu3) + c1p*w(i+1,j,k,itu3)
             dvt(i,j,k,2) = dvt(i,j,k,2)      + c2m*w(i-1,j,k,itu4) &
                          - c20*w(i,j,k,itu4) + c2p*w(i+1,j,k,itu4)

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
                             - b1*max(bmti1(j,k,itu3,itu3),zero)
               qq(i,j,k,1,2) = qq(i,j,k,1,2) - b1*bmti1(j,k,itu3,itu4)
               qq(i,j,k,2,1) = qq(i,j,k,2,1) - b2*bmti1(j,k,itu4,itu3)
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                             - b2*max(bmti1(j,k,itu4,itu4),zero)
             else if(i == il) then
               qq(i,j,k,1,1) = qq(i,j,k,1,1) + c1 &
                             - d1*max(bmti2(j,k,itu3,itu3),zero)
               qq(i,j,k,1,2) = qq(i,j,k,1,2) - d1*bmti2(j,k,itu3,itu4)
               qq(i,j,k,2,1) = qq(i,j,k,2,1) - d2*bmti2(j,k,itu4,itu3)
               qq(i,j,k,2,2) = qq(i,j,k,2,2) + c2 &
                             - d2*max(bmti2(j,k,itu4,itu4),zero)
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
             dw(i,j,k,itu3) = -vol(i,j,k)*dvt(i,j,k,1)*rblank
             dw(i,j,k,itu4) = -vol(i,j,k)*dvt(i,j,k,2)*rblank
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
               dw2    => dw(2,1:,1:,1:);    dvt2 => dvt(2,1:,1:,1:)
               w2     => w (2,1:,1:,1:);    rlv2 => rlv(2,1:,1:)
               w3     => w (3,1:,1:,1:);    rlv3 => rlv(3,1:,1:)
               d2Wall2=> d2Wall(2,:,:);     rev2 => rev(2,1:,1:)
               d2Wall3=> d2Wall(3,:,:);     rev3 => rev(3,1:,1:)

             case (iMax)
               flag    => flagIl
               dw2    => dw(il  ,1:,1:,1:); dvt2 => dvt(il  ,1:,1:,1:)
               w2     => w (il  ,1:,1:,1:); rlv2 => rlv(il  ,1:,1:)
               w3     => w (il-1,1:,1:,1:); rlv3 => rlv(il-1,1:,1:)
               d2Wall2=> d2Wall(il  ,:,:);  rev2 => rev(il  ,1:,1:)
               d2Wall3=> d2Wall(il-1,:,:);  rev3 => rev(il-1,1:,1:)

             case (jMin)
               flag    => flagJ2
               dw2    => dw(1:,2,1:,1:);    dvt2 => dvt(1:,2,1:,1:)
               w2     => w (1:,2,1:,1:);    rlv2 => rlv(1:,2,1:)
               w3     => w (1:,3,1:,1:);    rlv3 => rlv(1:,3,1:)
               d2Wall2=> d2Wall(:,2,:);     rev2 => rev(1:,2,1:)
               d2Wall3=> d2Wall(:,3,:);     rev3 => rev(1:,3,1:)

             case (jMax)
               flag    => flagJl
               dw2    => dw(1:,jl  ,1:,1:); dvt2 => dvt(1:,jl  ,1:,1:)
               w2     => w (1:,jl  ,1:,1:); rlv2 => rlv(1:,jl  ,1:)
               w3     => w (1:,jl-1,1:,1:); rlv3 => rlv(1:,jl-1,1:)
               d2Wall2=> d2Wall(:,jl  ,:);  rev2 => rev(1:,jl  ,1:)
               d2Wall3=> d2Wall(:,jl-1,:);  rev3 => rev(1:,jl-1,1:)

             case (kMin)
               flag    => flagK2
               dw2    => dw(1:,1:,2,1:);    dvt2 => dvt(1:,1:,2,1:)
               w2     => w (1:,1:,2,1:);    rlv2 => rlv(1:,1:,2)
               w3     => w (1:,1:,3,1:);    rlv3 => rlv(1:,1:,3)
               d2Wall2=> d2Wall(:,:,2);     rev2 => rev(1:,1:,2)
               d2Wall3=> d2Wall(:,:,3);     rev3 => rev(1:,1:,3)

             case (kMax)
               flag    => flagKl
               dw2    => dw(1:,1:,kl  ,1:); dvt2 => dvt(1:,1:,kl,1:)
               w2     => w (1:,1:,kl  ,1:); rlv2 => rlv(1:,1:,kl)
               w3     => w (1:,1:,kl-1,1:); rlv3 => rlv(1:,1:,kl-1)
               d2Wall2=> d2Wall(:,:,kl  );  rev2 => rev(1:,1:,kl)
               d2Wall3=> d2Wall(:,:,kl-1);  rev3 => rev(1:,1:,kl-1)

           end select

           ! Loop over the owned faces of this subface. Therefore the
           ! nodal range of BCData must be used. The offset of +1 is
           ! present, because the starting index of the cell range is
           ! 1 larger than the starting index of the nodal range.

           do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
             do i=(BCData(nn)%inBeg+1),BCData(nn)%inEnd

               ! Enforce v and f in the 1st internal cell from
               ! the wall function table. There is an offset of -1 in
               ! the wall distance. Note that the offset compared to
               ! the current value must be stored. Also note that the
               ! curve fits contain the non-dimensional values.

               utau = viscSubface(nn)%utau(i,j)
               yp = w2(i,j,irho)*d2Wall2(i-1,j-1)*utau/rlv2(i,j)

               ! Set dw2 to zero for proper monitoring of the
               ! convergence.

               dw2(i,j,itu3) = zero
               dw2(i,j,itu4) = zero

               ! Get table values

               call curveTupYp(tup(itu3:itu4), yp, itu3, itu4)
               tu32  = tup(itu3)*utau**2
               tu42  = tup(itu4)*utau**2/rlv2(i,j)*w2(i,j,irho)

               ! Compute f from balance

               if(rvfN .eq. 1) then
                  call curveTupYp(tup(itu1:itu2), yp, itu1, itu2)
                  tu12  = tup(itu1)*utau**2
                  tu22  = tup(itu2)*utau**4/rlv2(i,j)*w2(i,j,irho)
                  call curveTupYp(tup(itu5:itu5), yp, itu5, itu5)
                  tu52  = tup(itu5)*rlv2(i,j)
                  dtu23 = (w3(i,j,itu3)-tu32) &
                        / (d2Wall3(i-1,j-1)-d2Wall2(i-1,j-1))
                  rnu23 = half*( (tu52     +rlv2(i,j))/w2(i,j,irho) + &
                                 (rev3(i,j)+rlv3(i,j))/w3(i,j,irho) )
                  tu42  = (tu22*tu32/tu12-rnu23*dtu23 &
                        / (two*d2Wall2(i-1,j-1)))/tu12
               endif

               ! Set rhs to turbulence variables divide by betaTurb
               ! because the update is scaled by betaTurb.
               ! (see end of routine)

               dvt2(i,j,1) = (tu32 - w2(i,j,itu3))/betaTurb
               dvt2(i,j,2) = (tu42 - w2(i,j,itu4))/betaTurb
               if(rvfN .eq. 1) dvt2(i,j,2) = (tu42 - w2(i,j,itu4))*0.01

               ! Set the wall flag to .true.

               flag(i,j) = .true.

             enddo
           enddo

         enddo bocos
       endif testWallFunctions

       ! Return if only the residual must be computed.

       if( resOnly ) return

       do k=2,kl
         do j=2,jl
           do i=2,il

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

             c2m = ttm
             c2p = ttp

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
             bb(2,j) = bb(2,j)
             dd(2,j) = dd(2,j)

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

             c2m = ttm
             c2p = ttp

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
             bb(2,i) = bb(2,i)
             dd(2,i) = dd(2,i)

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
!      * dd-ADI step in k-direction. As we solve in k-direction, the    *
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

             c2m = ttm
             c2p = ttp

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
             bb(2,k) = bb(2,k)
             dd(2,k) = dd(2,k)

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
!      * Update the turbulent variables.                                *
!      *                                                                *
!      ******************************************************************
!
       do k=2,kl
         do j=2,jl
           do i=2,il
             w(i,j,k,itu3) = w(i,j,k,itu3) + betaTurb*dvt(i,j,k,1)
             w(i,j,k,itu4) = w(i,j,k,itu4) + betaTurb*dvt(i,j,k,2)
           enddo
         enddo
       enddo

       end subroutine vfSolve
