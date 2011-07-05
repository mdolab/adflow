!
!      ******************************************************************
!      *                                                                *
!      * File:          extrapolate2ndHalo.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-10-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine extrapolate2ndHalo(nn, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * extrapolate2ndHalo determines the states of the second layer   *
!      * halo cells for the given subface of the block. It is assumed   *
!      * that the pointers in blockPointers are already set to the      *
!      * correct block on the correct grid level.                       *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use constants
       use flowVarRefState
       use iteration
       use inputPhysics
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn

       logical, intent(in) :: correctForK
!
!      Local parameter.
!
       real(kind=realType), parameter :: factor = 0.5_realType
!
!      Local variables.
!
       integer(kind=intType) :: i, j, l, idim, ddim

       integer(kind=intType), dimension(3,2) :: crange

       real(kind=realType), dimension(:,:,:), pointer :: ww0, ww1, ww2
       real(kind=realType), dimension(:,:),   pointer :: pp0, pp1, pp2
       real(kind=realType), dimension(:,:),   pointer :: rlv0, rlv1
       real(kind=realType), dimension(:,:),   pointer :: rev0, rev1
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
         end subroutine setBcPointers
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Nullify the pointers and set them to the correct subface.
       ! They are nullified first, because some compilers require that.
       ! Note that rlv0 and rev0 are used here as dummies.

       nullify(ww1, ww2, pp1, pp2, rlv1, rlv0, rev1, rev0)
       call setBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv0, &
                          rev1, rev0, 0_intType)

       ! Set a couple of additional variables needed for the
       ! extrapolation. This depends on the block face on which the
       ! subface is located.

       select case (BCFaceID(nn))

         case (iMin)
           ww0 => w(0,1:,1:,:); pp0 => p(0,1:,1:)
           if( viscous )   rlv0 => rlv(0,1:,1:)
           if( eddyModel ) rev0 => rev(0,1:,1:)
           idim = 1; ddim = 0

         case (iMax)
           ww0 => w(ib,1:,1:,:); pp0 => p(ib,1:,1:)
           if( viscous )   rlv0 => rlv(ib,1:,1:)
           if( eddyModel ) rev0 => rev(ib,1:,1:)
           idim = 1; ddim = ib

         case (jMin)
           ww0 => w(1:,0,1:,:); pp0 => p(1:,0,1:)
           if( viscous )   rlv0 => rlv(1:,0,1:)
           if( eddyModel ) rev0 => rev(1:,0,1:)
           idim = 2; ddim = 0

         case (jMax)
           ww0 => w(1:,jb,1:,:); pp0 => p(1:,jb,1:)
           if( viscous )   rlv0 => rlv(1:,jb,1:)
           if( eddyModel ) rev0 => rev(1:,jb,1:)
           idim = 2; ddim = jb

         case (kMin)
           ww0 => w(1:,1:,0,:); pp0 => p(1:,1:,0)
           if( viscous )   rlv0 => rlv(1:,1:,0)
           if( eddyModel ) rev0 => rev(1:,1:,0)
           idim = 3; ddim = 0

         case (kMax)
           ww0 => w(1:,1:,kb,:); pp0 => p(1:,1:,kb)
           if( viscous )   rlv0 => rlv(1:,1:,kb)
           if( eddyModel ) rev0 => rev(1:,1:,kb)
           idim = 3; ddim = kb

       end select

       ! Loop over the generic subface to set the state in the halo's.

       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
         do i=BCData(nn)%icBeg, BCData(nn)%icEnd

           ! Extrapolate the density, momentum and pressure.
           ! Make sure that a certain threshold is kept.

           ww0(i,j,irho) = two*ww1(i,j,irho) - ww2(i,j,irho)
           ww0(i,j,irho) = max(factor*ww1(i,j,irho),ww0(i,j,irho))

           ww0(i,j,ivx) = two*ww1(i,j,ivx) - ww2(i,j,ivx)
           ww0(i,j,ivy) = two*ww1(i,j,ivy) - ww2(i,j,ivy)
           ww0(i,j,ivz) = two*ww1(i,j,ivz) - ww2(i,j,ivz)

           pp0(i,j) = max(factor*pp1(i,j),two*pp1(i,j) - pp2(i,j))

           ! Extrapolate the turbulent variables. Use constant
           ! extrapolation.

           do l=nt1MG,nt2MG
             ww0(i,j,l) = ww1(i,j,l)
           enddo

           ! The laminar and eddy viscosity, if present. These values
           ! are simply taken constant. Their values do not matter.

           if( viscous )   rlv0(i,j) = rlv1(i,j)
           if( eddyModel ) rev0(i,j) = rev1(i,j)

         enddo
       enddo

       ! Set the range for the halo cells for the energy computation.

       crange(1,1) = icBeg(nn); crange(1,2) = icEnd(nn)
       crange(2,1) = jcBeg(nn); crange(2,2) = jcEnd(nn)
       crange(3,1) = kcBeg(nn); crange(3,2) = kcEnd(nn)

       crange(idim,1) = ddim; crange(idim,2) = ddim

       ! Compute the energy for this halo range.

       call computeEtot(crange(1,1), crange(1,2), crange(2,1), &
                        crange(2,2), crange(3,1), crange(3,2), &
                        correctForK)

       end subroutine extrapolate2ndHalo
