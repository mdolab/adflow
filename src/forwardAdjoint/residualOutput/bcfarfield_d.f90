   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
   !
   !  Differentiation of bcfarfield in forward (tangent) mode:
   !   variations   of useful results: *p *gamma *w
   !   with respect to varying inputs: *p *w gammaconstant
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          bcFarfield.f90                                  *
   !      * Author:        Edwin van der Weide                             *
   !      * Starting date: 03-10-2003                                      *
   !      * Last modified: 06-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE BCFARFIELD_D(secondhalo, correctfork)
   USE FLOWVARREFSTATE
   USE BLOCKPOINTERS_D
   USE BCTYPES
   USE INPUTPHYSICS
   USE CONSTANTS
   USE ITERATION
   IMPLICIT NONE
   !close (UNIT=unitx)
   !
   !      ******************************************************************
   !      *                                                                *
   !      * bcFarfield applies the farfield boundary condition to a block. *
   !      * It is assumed that the pointers in blockPointers are already   *
   !      * set to the correct block on the correct grid level.            *
   !      *                                                                *
   !      ******************************************************************
   !
   !
   !      Subroutine arguments.
   !
   LOGICAL, INTENT(IN) :: secondhalo, correctfork
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: nn, i, j, l
   REAL(kind=realtype) :: nnx, nny, nnz
   REAL(kind=realtype) :: gm1, ovgm1, ac1, ac2
   REAL(kind=realtype) :: ac1d, ac2d
   REAL(kind=realtype) :: r0, u0, v0, w0, qn0, vn0, c0, s0
   REAL(kind=realtype) :: re, ue, ve, we, qne, ce
   REAL(kind=realtype) :: red, ued, ved, wed, qned, ced
   REAL(kind=realtype) :: qnf, cf, uf, vf, wf, sf, cc, qq
   REAL(kind=realtype) :: qnfd, cfd, ufd, vfd, wfd, sfd, ccd
   ! Variables Added for forward AD
   REAL(kind=realtype) :: rho, sf2
   REAL(kind=realtype) :: rhod
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1d, ww2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1d, pp2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: gamma2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: gamma2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   REAL(kind=realtype) :: arg1
   REAL(kind=realtype) :: arg1d
   REAL(kind=realtype) :: pwr1
   REAL(kind=realtype) :: pwr1d
   REAL(kind=realtype) :: pwx1
   REAL(kind=realtype) :: pwx1d
   INTRINSIC SQRT
   INTERFACE 
   SUBROUTINE SETBCPOINTERS2(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
   &        rev1, rev2, offset)
   USE BLOCKPOINTERS_D
   INTEGER(kind=inttype), INTENT(IN) :: nn, offset
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   END SUBROUTINE SETBCPOINTERS2
   END INTERFACE
      INTERFACE 
   SUBROUTINE SETBCPOINTERS_D(nn, ww1, ww1d, ww2, ww2d, pp1, pp1d, &
   &        pp2, pp2d, rlv1, rlv2, rev1, rev2, offset)
   USE BLOCKPOINTERS_D
   INTEGER(kind=inttype), INTENT(IN) :: nn, offset
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1, ww2
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: ww1d, ww2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1, pp2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: pp1d, pp2d
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rlv1, rlv2
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: rev1, rev2
   END SUBROUTINE SETBCPOINTERS_D
   END INTERFACE
      !!$       !File Parameters remove for AD
   !!$       integer :: unitx = 12,ierror
   !!$      integer ::iii,iiii,jjj,jjjj,kkk,kkkk,nnnn,istart2,jstart2,kstart2,iend2,jend2,kend2,n,ii,jj
   !!$      character(len = 16)::outfile
   !!$      
   !!$      outfile = "xoriginal.txt"
   !!$      
   !!$      open (UNIT=unitx,File=outfile,status='old',position='append',action='write',iostat=ierror)
   !!$      if(ierror /= 0)                        &
   !!$           call terminate("verifyResiduals", &
   !!$           "Something wrong when &
   !!$           &calling open")
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Some constants needed to compute the riemann inVariants.
   gm1 = gammainf - one
   ovgm1 = one/gm1
   ! Compute the three velocity components, the speed of sound and
   ! the entropy of the free stream.
   r0 = one/winf(irho)
   u0 = winf(ivx)
   v0 = winf(ivy)
   w0 = winf(ivz)
   arg1 = gammainf*pinfcorr*r0
   c0 = SQRT(arg1)
   pwr1 = winf(irho)**gammainf
   s0 = pwr1/pinfcorr
   gammad = 0.0
   ! Loop over the boundary condition subfaces of this block.
   bocos:DO nn=1,nbocos
   ! Check for farfield boundary conditions.
   IF (bctype(nn) .EQ. farfield) THEN
   ! Nullify the pointers and set them to the correct subface.
   ! They are nullified first, because some compilers require
   ! that.
   !nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)
   CALL SETBCPOINTERS_D(nn, ww1, ww1d, ww2, ww2d, pp1, pp1d, pp2, &
   &                     pp2d, rlv1, rlv2, rev1, rev2, 0_intType)
   ! Set the additional pointer for gamma2.
   SELECT CASE  (bcfaceid(nn)) 
   CASE (imin) 
   gamma2d => gammad(2, 1:, 1:)
   gamma2 => gamma(2, 1:, 1:)
   CASE (imax) 
   gamma2d => gammad(il, 1:, 1:)
   gamma2 => gamma(il, 1:, 1:)
   CASE (jmin) 
   gamma2d => gammad(1:, 2, 1:)
   gamma2 => gamma(1:, 2, 1:)
   CASE (jmax) 
   gamma2d => gammad(1:, jl, 1:)
   gamma2 => gamma(1:, jl, 1:)
   CASE (kmin) 
   gamma2d => gammad(1:, 1:, 2)
   gamma2 => gamma(1:, 1:, 2)
   CASE (kmax) 
   gamma2d => gammad(1:, 1:, kl)
   gamma2 => gamma(1:, 1:, kl)
   END SELECT
   ! Loop over the generic subface to set the state in the
   ! halo cells.
   DO j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
   DO i=bcdata(nn)%icbeg,bcdata(nn)%icend
   ! Store the three components of the unit normal a
   ! bit easier.
   nnx = bcdata(nn)%norm(i, j, 1)
   nny = bcdata(nn)%norm(i, j, 2)
   nnz = bcdata(nn)%norm(i, j, 3)
   !!$       !print out pAdj
   !!$       istart2 = -1!2
   !!$       jstart2 = -1!2
   !!$       kstart2 = -1!2
   !!$       iend2 = 1!2
   !!$       jend2 = 1!2
   !!$       kend2 = 1!2 
   !!$       if(i==BCData(nn)%icBeg) istart2=0
   !!$       if(j==BCData(nn)%jcBeg) jstart2= 0
   !!$       !if(i==BCData(nn)%icBeg+1) istart2=-1
   !!$       !if(j==BCData(nn)%jcBeg+1) jstart2=-1
   !!$!       if(kcell==2) kstart2=-1
   !!$       !if(i==BCData(nn)%icEnd-1) iend2=1
   !!$       !if(j==BCData(nn)%jcEnd-1) jend2=1
   !!$       if(i==BCData(nn)%icEnd) iend2=0
   !!$       if(j==BCData(nn)%jcEnd) jend2=0
   !!$!       if(kcell==kl) kend2=1
   !!$       do jjjj = jstart2,jend2
   !!$          do iiii = istart2,iend2
   !!$             !do kkkk = kstart2,kend2
   !!$               ! do n = 1,3!nw
   !!$                   !do n = 1,1!nw
   !!$                   !do n = 1,nw 
   !!$                   !do sps2 = 1,nTimeIntervalsSpectral
   !!$                   ii = i+iiii
   !!$                   jj = j+jjjj
   !!$                   !k = kcell+kkkk
   !!$                   !print *,'indices',i,j,iiii,jjjj,ii,jj,BCData(nn)%jcEnd, BCData(nn)%icEnd
   !!$
   !!$                   write(unitx,11)i,j,ii,jj,nn,BCData(nn)%norm(ii,jj,1), BCData(nn)%rface(ii,jj)
   !!$11               format(1x,'wadj',5I8,2f20.14) 
   !!$               ! end do
   !!$             end do
   !!$          end do
   ! Compute the normal velocity of the free stream and
   ! substract the normal velocity of the mesh.
   qn0 = u0*nnx + v0*nny + w0*nnz
   vn0 = qn0 - bcdata(nn)%rface(i, j)
   ! Compute the three velocity components, the normal
   ! velocity and the speed of sound of the current state
   ! in the internal cell.
   red = -(one*ww2d(i, j, irho)/ww2(i, j, irho)**2)
   re = one/ww2(i, j, irho)
   ued = ww2d(i, j, ivx)
   ue = ww2(i, j, ivx)
   ved = ww2d(i, j, ivy)
   ve = ww2(i, j, ivy)
   wed = ww2d(i, j, ivz)
   we = ww2(i, j, ivz)
   qned = nnx*ued + nny*ved + nnz*wed
   qne = ue*nnx + ve*nny + we*nnz
   arg1d = (gamma2d(i, j)*re+gamma2(i, j)*red)*pp2(i, j) + gamma2&
   &            (i, j)*re*pp2d(i, j)
   arg1 = gamma2(i, j)*pp2(i, j)*re
   IF (arg1 .EQ. 0.0) THEN
   ced = 0.0
   ELSE
   ced = arg1d/(2.0*SQRT(arg1))
   END IF
   ce = SQRT(arg1)
   ! Compute the new values of the riemann inVariants in
   ! the halo cell. Either the value in the internal cell
   ! is taken (positive sign of the corresponding
   ! eigenvalue) or the free stream value is taken
   ! (otherwise).
   IF (vn0 .GT. -c0) THEN
   ! Outflow or subsonic inflow.
   ac1d = qned + two*ovgm1*ced
   ac1 = qne + two*ovgm1*ce
   ELSE
   ! Supersonic inflow.
   ac1 = qn0 + two*ovgm1*c0
   ac1d = 0.0
   END IF
   IF (vn0 .GT. c0) THEN
   ! Supersonic outflow.
   ac2d = qned - two*ovgm1*ced
   ac2 = qne - two*ovgm1*ce
   ELSE
   ! Inflow or subsonic outflow.
   ac2 = qn0 - two*ovgm1*c0
   ac2d = 0.0
   END IF
   qnfd = half*(ac1d+ac2d)
   qnf = half*(ac1+ac2)
   cfd = fourth*gm1*(ac1d-ac2d)
   cf = fourth*(ac1-ac2)*gm1
   IF (vn0 .GT. zero) THEN
   ! Outflow.
   ufd = ued + nnx*(qnfd-qned)
   uf = ue + (qnf-qne)*nnx
   vfd = ved + nny*(qnfd-qned)
   vf = ve + (qnf-qne)*nny
   wfd = wed + nnz*(qnfd-qned)
   wf = we + (qnf-qne)*nnz
   !Intermediate rho variable added to fix AD bug,ww2 
   ! was not getting picked up here.
   rhod = ww2d(i, j, irho)
   rho = ww2(i, j, irho)
   IF (rho .GT. 0.0 .OR. (rho .LT. 0.0 .AND. gamma2(i, j) .EQ. &
   &                INT(gamma2(i, j)))) THEN
   pwr1d = gamma2(i, j)*rho**(gamma2(i, j)-1)*rhod
   ELSE IF (rho .EQ. 0.0 .AND. gamma2(i, j) .EQ. 1.0) THEN
   pwr1d = rhod
   ELSE
   pwr1d = 0.0
   END IF
   pwr1 = rho**gamma2(i, j)
   sfd = (pwr1d*pp2(i, j)-pwr1*pp2d(i, j))/pp2(i, j)**2
   sf = pwr1/pp2(i, j)
   !print *,'sf1',sfd,ww2d(i, j, irho),gamma2d(i, j),pp2d(i, j)
   !print *,'pwr1d',pwr1d,gamma2(i,j)*ww2(i, j, irho)**(gamma2(i, j)-1)*ww2d(i, j, irho)
   !old version
   !sf2 = ww2(i,j,irho)**gamma2(i,j)/pp2(i,j)
   DO l=nt1mg,nt2mg
   ww1d(i, j, l) = ww2d(i, j, l)
   ww1(i, j, l) = ww2(i, j, l)
   END DO
   ELSE
   ! Inflow
   ufd = nnx*qnfd
   uf = u0 + (qnf-qn0)*nnx
   vfd = nny*qnfd
   vf = v0 + (qnf-qn0)*nny
   wfd = nnz*qnfd
   wf = w0 + (qnf-qn0)*nnz
   sf = s0
   DO l=nt1mg,nt2mg
   ww1d(i, j, l) = 0.0
   ww1(i, j, l) = winf(l)
   END DO
   sfd = 0.0
   END IF
   ! Compute the density, velocity and pressure in the
   ! halo cell.
   ccd = ((cfd*cf+cf*cfd)*gamma2(i, j)-cf**2*gamma2d(i, j))/&
   &            gamma2(i, j)**2
   cc = cf*cf/gamma2(i, j)
   qq = uf*uf + vf*vf + wf*wf
   pwx1d = sfd*cc + sf*ccd
   pwx1 = sf*cc
   IF (pwx1 .GT. 0.0 .OR. (pwx1 .LT. 0.0 .AND. ovgm1 .EQ. INT(&
   &              ovgm1))) THEN
   ww1d(i, j, irho) = ovgm1*pwx1**(ovgm1-1)*pwx1d
   ELSE IF (pwx1 .EQ. 0.0 .AND. ovgm1 .EQ. 1.0) THEN
   ww1d(i, j, irho) = pwx1d
   ELSE
   ww1d(i, j, irho) = 0.0
   END IF
   ww1(i, j, irho) = pwx1**ovgm1
   ww1d(i, j, ivx) = ufd
   ww1(i, j, ivx) = uf
   ww1d(i, j, ivy) = vfd
   ww1(i, j, ivy) = vf
   ww1d(i, j, ivz) = wfd
   ww1(i, j, ivz) = wf
   pp1d(i, j) = ww1d(i, j, irho)*cc + ww1(i, j, irho)*ccd
   pp1(i, j) = ww1(i, j, irho)*cc
   !print *,'state1',ww1d(i,j,irho),i,j
   ! Simply set the laminar and eddy viscosity to
   ! the value in the donor cell. Their values do
   ! not matter too much in the far field.
   IF (viscous) rlv1(i, j) = rlv2(i, j)
   IF (eddymodel) rev1(i, j) = rev2(i, j)
   END DO
   END DO
   ! Compute the energy for these halo's.
   CALL COMPUTEETOT_D(icbeg(nn), icend(nn), jcbeg(nn), jcend(nn), &
   &                   kcbeg(nn), kcend(nn), correctfork)
   ! Extrapolate the state vectors in case a second halo
   ! is needed.
   IF (secondhalo) CALL EXTRAPOLATE2NDHALO_D(nn, correctfork)
   END IF
   END DO bocos
   END SUBROUTINE BCFARFIELD_D
