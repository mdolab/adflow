!
!      ******************************************************************
!      *                                                                *
!      * File:          bcSymmPolar.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-02-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcSymmPolar(secondHalo)
!
!      ******************************************************************
!      *                                                                *
!      * bcSymmPolar applies the polar symmetry boundary conditions     *
!      * to a singular line of a block. It is assumed that the pointers *
!      * in blockPointers are already set to the correct block on the   *
!      * correct grid level.                                            *
!      * The polar symmetry condition is a special case of a degenerate *
!      * line, as this line is the axi-symmetric centerline.            *
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
       integer(kind=intType) :: i, j, l, kk, mm, nn

       real(kind=realType) :: nnx, nny, nnz, tmp, vtx, vty, vtz

       real(kind=realType), dimension(:,:,:), pointer :: xline
       real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
       real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
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
       ! Set the value of kk; kk == 0 means only single halo, kk == 1
       ! double halo.

       kk = 0
       if( secondHalo ) kk = 1

       ! Loop over the boundary condition subfaces of this block.

       bocos: do nn=1,nBocos

         ! Check for the polar symmetry boundary condition.

         symmetryPolar: if(BCType(nn) == symmPolar) then

           ! Nullify the pointers, because some compilers require that.

           nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)

           ! Set the pointers for for the coordinates of the centerline.
           ! Needed to determine the direction of the velocity.
           ! This depends on the block face on which this subface is
           ! located.

           select case (BCFaceID(nn))
             case (iMin)
               xline => x(1,:,:,:)
             case (iMax)
               xline => x(il,:,:,:)
             case (jMin)
               xline => x(:,1,:,:)
             case (jMax)
               xline => x(:,jl,:,:)
             case (kMin)
               xline => x(:,:,1,:)
             case (kMax)
               xline => x(:,:,kl,:)
           end select

           ! Loop over the number of times the halo computation must
           ! be done.

           nHalo: do mm=0,kk

             ! Set the pointers to the correct subface.

             call setBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                                rev1, rev2, mm)

             ! Loop over the generic subface to set the state in the
             ! halo cells.

             do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
               do i=BCData(nn)%icBeg, BCData(nn)%icEnd

                 ! Determine the unit vector along the degenerated face.
                 ! However it is not known which is the singular
                 ! direction and therefore determine the direction along
                 ! the diagonal (i,j) -- (i-1,j-1), which is correct for
                 ! both singular i and j-direction. Note that due to the
                 ! usage of the pointer xline there is an offset of +1
                 ! in the indices and therefore (i+1,j+1) - (i,j) must
                 ! be used to determine this vector.

                 nnx = xline(i+1,j+1,1) - xline(i,j,1)
                 nny = xline(i+1,j+1,2) - xline(i,j,2)
                 nnz = xline(i+1,j+1,3) - xline(i,j,3)

                 ! Determine the unit vector in this direction.

                 tmp = one/sqrt(nnx*nnx + nny*nny + nnz*nnz)
                 nnx = nnx*tmp
                 nny = nny*tmp
                 nnz = nnz*tmp

                 ! Determine twice the tangential velocity vector of the
                 ! internal cell.

                 tmp = two*(ww2(i,j,ivx)*nnx + ww2(i,j,ivy)*nny &
                     +      ww2(i,j,ivz)*nnz)
                 vtx = tmp*nnx
                 vty = tmp*nny
                 vtz = tmp*nnz

                 ! Determine the flow variables in the halo cell. The
                 ! velocity is constructed such that the average of the
                 ! internal and the halo cell is along the centerline.
                 ! Note that the magnitude of the velocity does not
                 ! change and thus the energy is identical.

                 ww1(i,j,irho)  = ww2(i,j,irho)
                 ww1(i,j,ivx)   = vtx - ww2(i,j,ivx)
                 ww1(i,j,ivy)   = vty - ww2(i,j,ivy)
                 ww1(i,j,ivz)   = vtz - ww2(i,j,ivz)
                 ww1(i,j,irhoE) = ww2(i,j,irhoE)

                 ! Simply copy the turbulent variables.

                 do l=nt1MG,nt2MG
                   ww1(i,j,l) = ww2(i,j,l)
                 enddo

                 ! Set the pressure and possibly the laminar and
                 ! eddy viscosity in the halo.

                 pp1(i,j) = pp2(i,j)
                 if( viscous )   rlv1(i,j) = rlv2(i,j)
                 if( eddyModel ) rev1(i,j) = rev2(i,j)

               enddo
             enddo

           enddo nHalo

         endif symmetryPolar
       enddo bocos

       end subroutine bcSymmPolar
