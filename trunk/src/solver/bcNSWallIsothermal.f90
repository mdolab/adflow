!
!      ******************************************************************
!      *                                                                *
!      * File:          bcNSWallIsothermal.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-10-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcNSWallIsothermal(secondHalo, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * bcNSWallAdiabatic applies the viscous isothermal wall          *
!      * boundary condition to a block. It is assumed that the pointers *
!      * in blockPointers are already set to the correct block on the   *
!      * correct grid level.                                            *
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
       integer(kind=intType) :: nn, i, j

       real(kind=realType) :: rhok, t2, t1

       real(kind=realType), dimension(:,:,:), pointer :: uSlip
       real(kind=realType), dimension(:,:),   pointer :: TNS_Wall

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
       ! In case the turbulent transport equations are solved
       ! together with the mean flow equations, aplly the viscous
       ! wall boundary conditions for the turbulent variables.
       ! No need to extrapolate the secondary halo's, because this
       ! is done in extrapolate2ndHalo.

       if( turbCoupled ) call turbBCNSWall(.false.)

       ! Loop over the viscous subfaces of this block. Note that
       ! these are numbered first.

       bocos: do nn=1,nViscBocos

         ! Check for isothermal viscous wall boundary conditions.

         viscWall: if(BCType(nn) == NSWallIsothermal) then

           ! Set the pointers for uSlip and TNSWall to make
           ! the code more readable.

           uSlip    => BCData(nn)%uSlip
           TNS_Wall => BCData(nn)%TNS_Wall

           ! Nullify the pointers and set them to the correct subface.
           ! They are nullified first, because some compilers require
           ! that.

           nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)
           call setBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                              rev1, rev2, 0_intType)

           ! Initialize rhok to zero. This will be overwritten if a
           ! correction for k must be applied.

           rhok = zero

           ! Loop over the generic subface to set the state in the
           ! halo cells.

           do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd

               ! Set the value of rhok if a correcton must be applied.
               ! It probably does not matter too much, because k is very
               ! small near the wall.

               if( correctForK ) rhok = ww2(i,j,irho)*ww2(i,j,itu1)

               ! Compute the temperature in the internal cell and in the
               ! halo cell such that the average is the wall temperature.

               t2 = pp2(i,j)/(RGas*ww2(i,j,irho))
               t1 = two*TNS_Wall(i,j) - t2

               ! Make sure that t1 is within reasonable bounds. These
               ! bounds are such that the clipping is never active in the
               ! converged solution; it is only to avoid instabilities
               ! during the convergence.

               t1 = max(half*TNS_Wall(i,j), t1)
               t1 = min(two *TNS_Wall(i,j), t1)

               ! Determine the variables in the halo. As the spacing
               ! is very small a constant pressure boundary condition
               ! (except for the k correction) is okay. Take the slip
               ! velocity into account.

               pp1(i,j)      =  pp2(i,j) - four*third*rhok
               ww1(i,j,irho) =  pp1(i,j)/(RGas*t1)
               ww1(i,j,ivx)  = -ww2(i,j,ivx) + two*uSlip(i,j,1)
               ww1(i,j,ivy)  = -ww2(i,j,ivy) + two*uSlip(i,j,2)
               ww1(i,j,ivz)  = -ww2(i,j,ivz) + two*uSlip(i,j,3)

               ! Set the viscosities. There is no need to test for a
               ! viscous problem of course. The eddy viscosity is
               ! set to the negative value, as it should be zero on
               ! the wall.

               rlv1(i,j) = rlv2(i,j)
               if( eddyModel ) rev1(i,j) = -rev2(i,j)

             enddo
           enddo

           ! Compute the energy for these halo's.

           call computeEtot(icBeg(nn),icEnd(nn), jcBeg(nn),jcEnd(nn), &
                            kcBeg(nn),kcEnd(nn), correctForK)

           ! Extrapolate the state vectors in case a second halo
           ! is needed.

           if( secondHalo ) call extrapolate2ndHalo(nn, correctForK)

         endif viscWall
       enddo bocos

       end subroutine bcNSWallIsothermal
