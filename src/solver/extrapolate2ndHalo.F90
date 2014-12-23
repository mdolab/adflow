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

       select case (BCFaceID(nn))
       case (iMin)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd

                ! Extrapolate the density, momentum and pressure.
                ! Make sure that a certain threshold is kept.

                w(0,i,j,irho) = two*w(1,i,j,irho) - w(2,i,j,irho)
                w(0,i,j,irho) = max(factor*w(1,i,j,irho),w(0,i,j,irho))
                
                w(0,i,j,ivx) = two*w(1,i,j,ivx) - w(2,i,j,ivx)
                w(0,i,j,ivy) = two*w(1,i,j,ivy) - w(2,i,j,ivy)
                w(0,i,j,ivz) = two*w(1,i,j,ivz) - w(2,i,j,ivz)
                
                p(0,i,j) = max(factor*p(1,i,j),two*p(1,i,j) - p(2,i,j))

                ! Extrapolate the turbulent variables. Use constant
                ! extrapolation.

                do l=nt1MG,nt2MG
                   w(0,i,j,l) = w(1,i,j,l)
                enddo

                ! The laminar and eddy viscosity, if present. These values
                ! are simply taken constant. Their values do not matter.

                if( viscous )   rlv(0,i,j) = rlv(1,i,j)
                if( eddyModel ) rev(0,i,j) = rev(1,i,j)
             enddo
          enddo
          idim = 1; ddim = 0

       case (iMax)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                w(ib,i,j,irho) = two*w(ie,i,j,irho) - w(il,i,j,irho)
                w(ib,i,j,irho) = max(factor*w(ie,i,j,irho),w(ib,i,j,irho))
                w(ib,i,j,ivx) = two*w(ie,i,j,ivx) - w(il,i,j,ivx)
                w(ib,i,j,ivy) = two*w(ie,i,j,ivy) - w(il,i,j,ivy)
                w(ib,i,j,ivz) = two*w(ie,i,j,ivz) - w(il,i,j,ivz)
                p(ib,i,j) = max(factor*p(ie,i,j),two*p(ie,i,j) - p(il,i,j))
                do l=nt1MG,nt2MG
                   w(ib,i,j,l) = w(ie,i,j,l)
                enddo
                if( viscous )   rlv(ib,i,j) = rlv(ie,i,j)
                if( eddyModel ) rev(ib,i,j) = rev(ie,i,j)
             enddo 
          enddo
          idim = 1; ddim = ib
       case (jMin)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                w(i,0,j,irho) = two*w(i,1,j,irho) - w(i,2,j,irho)
                w(i,0,j,irho) = max(factor*w(i,1,j,irho),w(i,0,j,irho))
                w(i,0,j,ivx) = two*w(i,1,j,ivx) - w(i,2,j,ivx)
                w(i,0,j,ivy) = two*w(i,1,j,ivy) - w(i,2,j,ivy)
                w(i,0,j,ivz) = two*w(i,1,j,ivz) - w(i,2,j,ivz)
                p(i,0,j) = max(factor*p(i,1,j),two*p(i,1,j) - p(i,2,j))
                do l=nt1MG,nt2MG
                   w(i,0,j,l) = w(i,1,j,l)
                end do
                if( viscous )   rlv(i,0,j) = rlv(i,1,j)
                if( eddyModel ) rev(i,0,j) = rev(i,1,j)
             enddo
          enddo
          idim = 2; ddim = 0
       case (jMax)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                w(i,jb,j,irho) = two*w(i,je,j,irho) - w(i,jl,j,irho)
                w(i,jb,j,irho) = max(factor*w(i,je,j,irho),w(i,jb,j,irho))
                w(i,jb,j,ivx) = two*w(i,je,j,ivx) - w(i,jl,j,ivx)
                w(i,jb,j,ivy) = two*w(i,je,j,ivy) - w(i,jl,j,ivy)
                w(i,jb,j,ivz) = two*w(i,je,j,ivz) - w(i,jl,j,ivz)
                p(i,jb,j) = max(factor*p(i,je,j),two*p(i,je,j) - p(i,jl,j))
                do l=nt1MG,nt2MG
                   w(i,jb,j,l) = w(i,je,j,l)
                end do
                if( viscous )   rlv(i,jb,j) = rlv(i,je,j)
                if( eddyModel ) rev(i,jb,j) = rev(i,je,j)
             enddo
          enddo
          idim = 2; ddim = jb
       case (kMin)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                w(i,j,0,irho) = two*w(i,j,1,irho) - w(i,j,2,irho)
                w(i,j,0,irho) = max(factor*w(i,j,1,irho),w(i,j,0,irho))
                w(i,j,0,ivx) = two*w(i,j,1,ivx) - w(i,j,2,ivx)
                w(i,j,0,ivy) = two*w(i,j,1,ivy) - w(i,j,2,ivy)
                w(i,j,0,ivz) = two*w(i,j,1,ivz) - w(i,j,2,ivz)
                p(i,j,0) = max(factor*p(i,j,1),two*p(i,j,1) - p(i,j,2))
                do l=nt1MG,nt2MG
                   w(i,j,0,l) = w(i,j,1,l)
                end do
                if( viscous )   rlv(i,j,0) = rlv(i,j,1)
                if( eddyModel ) rev(i,j,0) = rev(i,j,1)
             enddo
          enddo
          idim = 3; ddim = 0
       case (kMax)
          do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                w(i,j,kb,irho) = two*w(i,j,ke,irho) - w(i,j,kl,irho)
                w(i,j,kb,irho) = max(factor*w(i,j,ke,irho),w(i,j,kb,irho))
                w(i,j,kb,ivx) = two*w(i,j,ke,ivx) - w(i,j,kl,ivx)
                w(i,j,kb,ivy) = two*w(i,j,ke,ivy) - w(i,j,kl,ivy)
                w(i,j,kb,ivz) = two*w(i,j,ke,ivz) - w(i,j,kl,ivz)
                p(i,j,kb) = max(factor*p(i,j,ke),two*p(i,j,ke) - p(i,j,kl))
                do l=nt1MG,nt2MG
                   w(i,j,kb,l) = w(i,j,ke,l)
                end do
                if( viscous )   rlv(i,j,kb) = rlv(i,j,ke)
                if( eddyModel ) rev(i,j,kb) = rev(i,j,ke)
                enddo
             enddo
             idim = 3; ddim = kb
          end select

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
