!
!      ******************************************************************
!      *                                                                *
!      * File:          bcSubsonicOutflow.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-03-2003                                      *
!      * Last modified: 09-10-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcSubsonicOutflow(secondHalo, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * bcSubsonicOutflow applies the subsonic outflow boundary        *
!      * condition, static pressure prescribed, to a block. It is       *
!      * assumed that the pointers in blockPointers are already set to  *
!      * the correct block on the correct grid level.                   *
!      * Exactly the same boundary condition is also applied for an     *
!      * outflow mass bleed. Therefore the test is for both a subsonic  *
!      * outflow and an bleed outflow.                                  *
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
!      Local parameter.
!
       real(kind=realType), parameter :: twothird = two*third
!
!      Local variables.
!
       integer(kind=intType) :: i, j, l, nn

       real(kind=realType) :: ovg, ovgm1, nnx, nny, nnz
       real(kind=realType) :: pExit, pInt, r, a2, a, ac, ss
       real(kind=realType) :: ue, ve, we, qne, qnh

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
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the boundary condition subfaces of this block.

       bocos: do nn=1,nBocos

         ! Check for the subsonic outflow or outflow bleed
         ! boundary condition.

         outflowSubsonic: if(BCType(nn) == SubsonicOutflow .or. &
                             BCType(nn) == MassBleedOutflow) then

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

               ! Store a couple of variables, such as the static
               ! pressure and grid unit outward normal, a bit easier.

               pExit = BCData(nn)%ps(i,j)

               nnx = BCData(nn)%norm(i,j,1)
               nny = BCData(nn)%norm(i,j,2)
               nnz = BCData(nn)%norm(i,j,3)

               ! Abbreviate 1/gamma and 1/(gamma -1) a bit easier.

               ovg   = one/gamma2(i,j)
               ovgm1 = one/(gamma2(i,j)-one)

               ! Store the internal pressure and correct for the
               ! possible presence of a k-equation.

               pInt = pp2(i,j)
               if( correctForK ) &
                 pInt = pInt - twothird*ww2(i,j,irho)*ww2(i,j,itu1)

               ! Compute the velocity components, the normal velocity
               ! and the speed of sound for the internal cell.

               r   = one/ww2(i,j,irho)
               a2  = gamma2(i,j)*pInt*r
               a   = sqrt(a2)
               ue  = ww2(i,j,ivx)
               ve  = ww2(i,j,ivy)
               we  = ww2(i,j,ivz)
               qne = ue*nnx + ve*nny + we*nnz

               ! Compute the entropy and the acoustic variable.
               ! These riemann inVariants, as well as the tangential
               ! velocity components, are extrapolated.

               ss = pInt*(r**gamma2(i,j))
               ac = qne + two*a*ovgm1

               ! Compute the state in the halo.

               ww1(i,j,irho) = (pExit/ss)**ovg
               pp1(i,j)      = pExit
               a             = sqrt(gamma2(i,j)*pExit/ww1(i,j,irho))
               qnh           = ac - two*a*ovgm1
               ww1(i,j,ivx)  = ue + (qnh - qne)*nnx
               ww1(i,j,ivy)  = ve + (qnh - qne)*nny
               ww1(i,j,ivz)  = we + (qnh - qne)*nnz

               ! Extrapolate the primitive turbulent variables.

               do l=nt1MG,nt2MG
                 ww1(i,j,l) = ww2(i,j,l)
               enddo

               ! Correct the pressure if a k-equation is present.

               if( correctForK )   &
                 pp1(i,j) = pp1(i,j) &
                          + twothird*ww1(i,j,irho)*ww1(i,j,itu1)

               ! Set the viscosities in the halo to the viscosities
               ! in the donor cell.

               if( viscous )   rlv1(i,j) = rlv2(i,j)
               if( eddyModel ) rev1(i,j) = rev2(i,j)

             enddo
           enddo

           ! Compute the energy for these halo's.

           call computeEtot(icBeg(nn),icEnd(nn), jcBeg(nn),jcEnd(nn), &
                             kcBeg(nn),kcEnd(nn), correctForK)

           ! Extrapolate the state vectors in case a second halo
           ! is needed.

           if( secondHalo ) call extrapolate2ndHalo(nn, correctForK)

         endif outflowSubsonic
       enddo bocos

       end subroutine bcSubsonicOutflow
