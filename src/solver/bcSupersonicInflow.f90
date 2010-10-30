!
!      ******************************************************************
!      *                                                                *
!      * File:          bcSupersonicInflow.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-02-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcSupersonicInflow(secondHalo, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * bcSupersonicInflow applies the supersonic inflow boundary      *
!      * conditions, entire state vector is prescribed, to a block. It  *
!      * is assumed that the pointers in blockPointers are already set  *
!      * to the correct block on the correct grid level.                *
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
       logical, intent(in) :: secondHalo, correctForK
!
!      Local variables.
!
       integer(kind=intType) :: i, j, l, kk, mm, nn
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

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

         ! Check for the supersonic inflow boundary condition.

         inflowSupersonic: if(BCType(nn) == SupersonicInflow) then

           ! Loop over the number of halo cells.

           nHalo: do mm=0,kk

             ! Nullify the pointers, because some compilers require that.

             nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)

             ! Set the pointers to the correct subface.

             call setBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                                rev1, rev2, mm)

             ! Loop over the generic subface to set the state in the
             ! halo cells.

             do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
               do i=BCData(nn)%icBeg, BCData(nn)%icEnd

                 ww1(i,j,irho) = BCData(nn)%rho(i,j)
                 ww1(i,j,ivx)  = BCData(nn)%velx(i,j)
                 ww1(i,j,ivy)  = BCData(nn)%vely(i,j)
                 ww1(i,j,ivz)  = BCData(nn)%velz(i,j)
                 pp1(i,j)      = BCData(nn)%ps(i,j)

                 ! The turbulent variables.

                 do l=nt1MG,nt2MG
                   ww1(i,j,l) = BCData(nn)%turbInlet(i,j,l)
                 enddo

                 ! Set the laminar and eddy viscosity in the halo
                 ! if needed.

                 if( viscous )   rlv1(i,j) = rlv2(i,j)
                 if( eddyModel ) rev1(i,j) = rev2(i,j)

               enddo
             enddo

           enddo nHalo

           ! Set the halo range for the energy to be computed.

           iBeg = icBeg(nn); jBeg = jcBeg(nn); kBeg = kcBeg(nn)
           iEnd = icEnd(nn); jEnd = jcEnd(nn); kEnd = kcEnd(nn)

           if( secondHalo ) then
             select case (BCFaceID(nn))
               case (iMin); iBeg = 0
               case (iMax); iEnd = ib
               case (jMin); jBeg = 0
               case (jMax); jEnd = jb
               case (kMin); kBeg = 0
               case (kMax); kEnd = kb
             end select
           endif

           ! Compute the energy for these halo's.

           call computeEtot(iBeg, iEnd, jBeg, jEnd, kBeg,kEnd, &
                            correctForK)

         endif inflowSupersonic
       enddo bocos

       end subroutine bcSupersonicInflow
