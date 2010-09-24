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
       subroutine bcFarfield(secondHalo, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * bcFarfield applies the farfield boundary condition to a block. *
!      * It is assumed that the pointers in blockPointers are already   *
!      * set to the correct block on the correct grid level.            *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
       use constants
       use flowVarRefState
       use inputPhysics
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: secondHalo, correctForK
!
!      Local variables.
!
       integer(kind=intType) :: nn, i, j, l

       real(kind=realType) :: nnx, nny, nnz
       real(kind=realType) :: gm1, ovgm1, ac1, ac2
       real(kind=realType) :: r0, u0, v0, w0, qn0, vn0, c0, s0
       real(kind=realType) :: re, ue, ve, we, qne, ce
       real(kind=realType) :: qnf, cf, uf, vf, wf, sf, cc, qq

       real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
       real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
       real(kind=realType), dimension(:,:),   pointer :: gamma2
       real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
       real(kind=realType), dimension(:,:),   pointer :: rev1, rev2
!
!      Interfaces
!
       interface
         subroutine setBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                                  rev1, rev2, offset)
           use blockPointers
           implicit none

           integer(kind=intType), intent(in) :: nn, offset
           real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
           real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
           real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
           real(kind=realType), dimension(:,:),   pointer :: rev1, rev2
         end subroutine setBCPointers
       end interface

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

       gm1   = gammaInf -one
       ovgm1 = one/gm1

       ! Compute the three velocity components, the speed of sound and
       ! the entropy of the free stream.

       r0  = one/wInf(irho)
       u0  = wInf(ivx)
       v0  = wInf(ivy)
       w0  = wInf(ivz)
       c0  = sqrt(gammaInf*pInfCorr*r0)
       s0  = wInf(irho)**gammaInf/pInfCorr

       ! Loop over the boundary condition subfaces of this block.

       bocos: do nn=1,nBocos

         ! Check for farfield boundary conditions.

         testFarfield: if(BCType(nn) == FarField) then

           ! Nullify the pointers and set them to the correct subface.
           ! They are nullified first, because some compilers require
           ! that.

           nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)
           call setBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                              rev1, rev2, 0_intType)

           ! Set the additional pointer for gamma2.

           select case (BCFaceID(nn))
             case (iMin)
               gamma2 => gamma(2,1:,1:)
             case (iMax)
               gamma2 => gamma(il,1:,1:)
             case (jMin)
               gamma2 => gamma(1:,2,1:)
             case (jMax)
               gamma2 => gamma(1:,jl,1:)
             case (kMin)
               gamma2 => gamma(1:,1:,2)
             case (kMax)
               gamma2 => gamma(1:,1:,kl)
           end select

           ! Loop over the generic subface to set the state in the
           ! halo cells.

           do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd

               ! Store the three components of the unit normal a
               ! bit easier.

               nnx = BCData(nn)%norm(i,j,1)
               nny = BCData(nn)%norm(i,j,2)
               nnz = BCData(nn)%norm(i,j,3)

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
               vn0 = qn0 - BCData(nn)%rface(i,j)

               ! Compute the three velocity components, the normal
               ! velocity and the speed of sound of the current state
               ! in the internal cell.

               re  = one/ww2(i,j,irho)
               ue  = ww2(i,j,ivx)
               ve  = ww2(i,j,ivy)
               we  = ww2(i,j,ivz)
               qne = ue*nnx + ve*nny + we*nnz
               ce  = sqrt(gamma2(i,j)*pp2(i,j)*re)

               ! Compute the new values of the riemann inVariants in
               ! the halo cell. Either the value in the internal cell
               ! is taken (positive sign of the corresponding
               ! eigenvalue) or the free stream value is taken
               ! (otherwise).

               if(vn0 > -c0) then       ! Outflow or subsonic inflow.
                 ac1 = qne + two*ovgm1*ce
               else                     ! Supersonic inflow.
                 ac1 = qn0 + two*ovgm1*c0
               endif

               if(vn0 > c0) then        ! Supersonic outflow.
                 ac2 = qne - two*ovgm1*ce
               else                     ! Inflow or subsonic outflow.
                 ac2 = qn0 - two*ovgm1*c0
               endif

               qnf = half*  (ac1 + ac2)
               cf  = fourth*(ac1 - ac2)*gm1

               if(vn0 > zero) then            ! Outflow.

                 uf = ue + (qnf - qne)*nnx
                 vf = ve + (qnf - qne)*nny
                 wf = we + (qnf - qne)*nnz
                 sf = ww2(i,j,irho)**gamma2(i,j)/pp2(i,j)

                 do l=nt1MG,nt2MG
                   ww1(i,j,l) = ww2(i,j,l)
                 enddo

               else                           ! Inflow

                 uf = u0 + (qnf - qn0)*nnx
                 vf = v0 + (qnf - qn0)*nny
                 wf = w0 + (qnf - qn0)*nnz
                 sf = s0

                 do l=nt1MG,nt2MG
                   ww1(i,j,l) = wInf(l)
                 enddo

               endif

               ! Compute the density, velocity and pressure in the
               ! halo cell.

               cc = cf*cf/gamma2(i,j)
               qq = uf*uf + vf*vf + wf*wf
               ww1(i,j,irho) = (sf*cc)**ovgm1
               ww1(i,j,ivx)  = uf
               ww1(i,j,ivy)  = vf
               ww1(i,j,ivz)  = wf
               pp1(i,j)      = ww1(i,j,irho)*cc

               ! Simply set the laminar and eddy viscosity to
               ! the value in the donor cell. Their values do
               ! not matter too much in the far field.

               if( viscous )    rlv1(i,j) = rlv2(i,j)
               if( eddyModel ) rev1(i,j) = rev2(i,j)

             enddo
           enddo

           ! Compute the energy for these halo's.

           call computeEtot(icBeg(nn),icEnd(nn), jcBeg(nn),jcEnd(nn), &
                            kcBeg(nn),kcEnd(nn), correctForK)

           ! Extrapolate the state vectors in case a second halo
           ! is needed.

           if( secondHalo ) call extrapolate2ndHalo(nn, correctForK)

         endif testFarfield
       enddo bocos
!close (UNIT=unitx)
       end subroutine bcFarfield
