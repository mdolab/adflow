!
!      ******************************************************************
!      *                                                                *
!      * File:          bcSymm.f90                                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-07-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcSymm(secondHalo)
!
!      ******************************************************************
!      *                                                                *
!      * bcSymm applies the symmetry boundary conditions to a block.    *
!      * It is assumed that the pointers in blockPointers are already   *
!      * set to the correct block on the correct grid level.            *
!      *                                                                *
!      * In case also the second halo must be set the loop over the     *
!      * boundary subfaces is executed twice. This is the only correct  *
!      * way in case the block contains only 1 cell between two         *
!      * symmetry planes, i.e. a 2D problem.                            *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
       use constants
       use flowVarRefState
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: secondHalo
!
!      Local variables.
!
       integer(kind=intType) :: kk, mm, nn, i, j, l

       real(kind=realType) :: vn, nnx, nny, nnz

       real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
       real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
       real(kind=realType), dimension(:,:),   pointer :: gamma1, gamma2
       real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
       real(kind=realType), dimension(:,:),   pointer :: rev1, rev2
!
!      Interfaces
!
       interface
         subroutine setBcPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
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
       ! Set the value of kk; kk == 0 means only single halo, kk == 1
       ! double halo.

       kk = 0
       if( secondHalo ) kk = 1

       ! Loop over the number of times the halo computation must be done.

       nHalo: do mm=0,kk

         ! Loop over the boundary condition subfaces of this block.

         bocos: do nn=1,nBocos

           ! Check for symmetry boundary condition.

           symmetry: if(BCType(nn) == symm) then

             ! Nullify the pointers, because some compilers require that.

             nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)

             ! Set the pointers to the correct subface.

             call setBcPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                                rev1, rev2, mm)

             ! Set the additional pointers for gamma1 and gamma2.

             select case (BCFaceID(nn))
               case (iMin)
                 gamma1 => gamma(1, 1:,1:); gamma2 => gamma(2, 1:,1:)
               case (iMax)
                 gamma1 => gamma(ie,1:,1:); gamma2 => gamma(il,1:,1:)
               case (jMin)
                 gamma1 => gamma(1:,1, 1:); gamma2 => gamma(1:,2, 1:)
               case (jMax)
                 gamma1 => gamma(1:,je,1:); gamma2 => gamma(1:,jl,1:)
               case (kMin)
                 gamma1 => gamma(1:,1:,1 ); gamma2 => gamma(1:,1:,2 )
               case (kMax)
                 gamma1 => gamma(1:,1:,ke); gamma2 => gamma(1:,1:,kl)
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

                 ! Determine twice the normal velocity component,
                 ! which must be substracted from the donor velocity
                 ! to obtain the halo velocity.

                 vn = two*(ww2(i,j,ivx)*nnx + ww2(i,j,ivy)*nny &
                    +      ww2(i,j,ivz)*nnz)

                 ! Determine the flow variables in the halo cell.

                 ww1(i,j,irho) = ww2(i,j,irho)

                 ww1(i,j,ivx) = ww2(i,j,ivx) - vn*nnx
                 ww1(i,j,ivy) = ww2(i,j,ivy) - vn*nny
                 ww1(i,j,ivz) = ww2(i,j,ivz) - vn*nnz

                 ww1(i,j,irhoE) = ww2(i,j,irhoE)

                 ! Simply copy the turbulent variables.

                 do l=nt1MG,nt2MG
                   ww1(i,j,l) = ww2(i,j,l)
                 enddo

                 ! Set the pressure and gamma and possibly the
                 ! laminar and eddy viscosity in the halo.

                 gamma1(i,j) = gamma2(i,j)
                 pp1(i,j)    = pp2(i,j)
                 if( viscous )   rlv1(i,j) = rlv2(i,j)
                 if( eddyModel ) rev1(i,j) = rev2(i,j)

               enddo
             enddo

           endif symmetry
         enddo bocos
       enddo nHalo

       end subroutine bcSymm
