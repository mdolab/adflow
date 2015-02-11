!
!      ******************************************************************
!      *                                                                *
!      * File:          bcExtrap.f90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-10-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcExtrap(secondHalo, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * ccExtrap applies the extrapolation boundary condition to a     *
!      * block. It is assumed that the pointers in blockPointers are    *
!      * already set to the correct block on the correct grid level.    *
!      * Extrapolation boundaries are applied to both singular lines or *
!      * points of a block face and to supersonic outlets. They are     *
!      * marked differently because of postprocessing reasons, but      *
!      * their numerical treatment is identical.                        *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
       use constants
       use flowVarRefState
       use inputDiscretization
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
       real(kind=realType), parameter :: factor = 0.5
!
!      Local variables.
!
       integer(kind=intType) :: i, j, l, nn

       real(kind=realType) :: fw2, fw3

       real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2, ww3
       real(kind=realType), dimension(:,:),   pointer :: pp1, pp2, pp3
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

         ! Check for extrapolation or supersonic outlet boundary
         ! conditions.

         extraPolation: if(BCType(nn) == extrap .or. &
                           BCType(nn) == SupersonicOutflow) then

           ! Set the extrapolation weights, depending on the situation.

           if(BCType(nn) == SupersonicOutflow) then

             ! A physical outflow face. Set the weights depending
             ! on the input parameter.

             select case (outflowTreatment)
               case (constantExtrapol)
                 fw2 = one; fw3 = zero
               case (linExtrapol)
                 fw2 = two; fw3 = -one
             end select

           else

             ! Singular block boundary. Use linear extrapolation.

             fw2 = two; fw3 = -one

           endif

           ! Nullify the pointers and set them to the correct subface.
           ! They are nullified first, because some compilers require
           ! that.

           nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)
           call setBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                              rev1, rev2, 0_intType)

           ! Set the additional pointers for ww3 and pp3, depending
           ! on the block face.

           select case (BCFaceID(nn))
             case (iMin)
               ww3 => w(3,1:,1:,:);  pp3 => p(3,1:,1:)
             case (iMax)
               ww3 => w(nx,1:,1:,:); pp3 => p(nx,1:,1:)
             case (jMin)
               ww3 => w(1:,3,1:,:);  pp3 => p(1:,3,1:)
             case (jMax)
               ww3 => w(1:,ny,1:,:); pp3 => p(1:,ny,1:)
             case (kMin)
               ww3 => w(1:,1:,3,:);  pp3 => p(1:,1:,3)
             case (kMax)
               ww3 => w(1:,1:,nz,:); pp3 => p(1:,1:,nz)
           end select

           ! Determine the state in the halo cells for this subface.

           do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd

               ! Extrapolate the density, velocities and pressure.
               ! Make sure that a certain threshold is kept for the
               ! density and pressure.

               ww1(i,j,irho) = fw2*ww2(i,j,irho) + fw3*ww3(i,j,irho)
               ww1(i,j,irho) = max(factor*ww2(i,j,irho), ww1(i,j,irho))

               ww1(i,j,ivx) = fw2*ww2(i,j,ivx) + fw3*ww3(i,j,ivx)
               ww1(i,j,ivy) = fw2*ww2(i,j,ivy) + fw3*ww3(i,j,ivy)
               ww1(i,j,ivz) = fw2*ww2(i,j,ivz) + fw3*ww3(i,j,ivz)

               pp1(i,j) = fw2*pp2(i,j) + fw3*pp3(i,j)
               pp1(i,j) = max(factor*pp2(i,j), pp1(i,j))

               ! Extrapolate the turbulent variables.

               do l=nt1MG,nt2MG
                 ww1(i,j,l) = fw2*ww2(i,j,l) + fw3*ww3(i,j,l)
               enddo

               ! The laminar and eddy viscosity, if present. These
               ! values are simply taken constant. Their values do
               ! not really matter.

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

         endif extraPolation
       enddo bocos

       end subroutine bcExtrap
