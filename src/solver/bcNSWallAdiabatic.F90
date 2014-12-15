!
!      ******************************************************************
!      *                                                                *
!      * File:          bcNsWallAdiabatic.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-10-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcNSWallAdiabatic(secondHalo, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * bcNSWallAdiabatic applies the viscous adiabatic wall           *
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

       real(kind=realType) :: rhok
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

         ! Check for adiabatic viscous wall boundary conditions.

         adiabaticWall: if(BCType(nn) == NSWallAdiabatic) then

           ! Initialize rhok to zero. This will be overwritten if a
           ! correction for k must be applied.

           rhok = zero

           ! Loop over the generic subface to set the state in the
           ! halo cells.
           
           select case (BCFaceID(nn))

           case (iMin)
              do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                 do i=BCData(nn)%icBeg, BCData(nn)%icEnd

                    ! Set the value of rhok if a correcton must be applied.
                    ! It probably does not matter too much, because k is very
                    ! small near the wall.

                    if( correctForK ) rhok = w(2,i,j,irho)*w(2,i,j,itu1)

                    ! Determine the variables in the halo. As the spacing
                    ! is very small a constant pressure boundary condition
                    ! (except for the k correction) is okay. Take the slip
                    ! velocity into account.

                    w(1,i,j,irho) =  w(2,i,j,irho)
                    w(1,i,j,ivx)  = -w(2,i,j,ivx) + two*BCData(nn)%uSlip(i,j,1)
                    w(1,i,j,ivy)  = -w(2,i,j,ivy) + two*BCData(nn)%uSlip(i,j,2)
                    w(1,i,j,ivz)  = -w(2,i,j,ivz) + two*BCData(nn)%uSlip(i,j,3)
                    p(1,i,j)      =  p(2,i,j) - four*third*rhok

                    ! Set the viscosities. There is no need to test for a
                    ! viscous problem of course. The eddy viscosity is
                    ! set to the negative value, as it should be zero on
                    ! the wall.
                    
                    rlv(1,i,j) = rlv(2,i,j)
                    if( eddyModel ) rev(1,i,j) = -rev(2,i,j)
                 enddo
              enddo

           case (iMax)
              do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                 do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                    if( correctForK ) rhok = w(il,i,j,irho)*w(il,i,j,itu1)
                    w(ie,i,j,irho) =  w(il,i,j,irho)
                    w(ie,i,j,ivx)  = -w(il,i,j,ivx) + two*BCData(nn)%uSlip(i,j,1)
                    w(ie,i,j,ivy)  = -w(il,i,j,ivy) + two*BCData(nn)%uSlip(i,j,2)
                    w(ie,i,j,ivz)  = -w(il,i,j,ivz) + two*BCData(nn)%uSlip(i,j,3)
                    p(ie,i,j)      =  p(il,i,j) - four*third*rhok
                    rlv(ie,i,j) = rlv(il,i,j)
                    if( eddyModel ) rev(ie,i,j) = -rev(il,i,j)
                 enddo
              enddo
           case (jMin)
              do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                 do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                    if( correctForK ) rhok = w(i,2,j,irho)*w(i,2,j,itu1)
                    w(i,1,j,irho) =  w(i,2,j,irho)
                    w(i,1,j,ivx)  = -w(i,2,j,ivx) + two*BCData(nn)%uSlip(i,j,1)
                    w(i,1,j,ivy)  = -w(i,2,j,ivy) + two*BCData(nn)%uSlip(i,j,2)
                    w(i,1,j,ivz)  = -w(i,2,j,ivz) + two*BCData(nn)%uSlip(i,j,3)
                    p(i,1,j)      =  p(i,2,j) - four*third*rhok
                    rlv(i,1,j) = rlv(i,2,j)
                    if( eddyModel ) rev(i,1,j) = -rev(i,2,j)
                 enddo
              enddo
           case (jMax)
              do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                 do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                    if( correctForK ) rhok = w(i,jl,j,irho)*w(i,jl,j,itu1)
                    w(i,je,j,irho) =  w(i,jl,j,irho)
                    w(i,je,j,ivx)  = -w(i,jl,j,ivx) + two*BCData(nn)%uSlip(i,j,1)
                    w(i,je,j,ivy)  = -w(i,jl,j,ivy) + two*BCData(nn)%uSlip(i,j,2)
                    w(i,je,j,ivz)  = -w(i,jl,j,ivz) + two*BCData(nn)%uSlip(i,j,3)
                    p(i,je,j)      =  p(i,jl,j) - four*third*rhok
                    rlv(i,je,j) = rlv(i,jl,j)
                    if( eddyModel ) rev(i,je,j) = -rev(i,jl,j)
                 enddo
              enddo
           case (kMin)
              do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                 do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                    if( correctForK ) rhok = w(i,j,1,irho)*w(i,j,1,itu1)
                    w(i,j,1,irho) =  w(i,j,2,irho)
                    w(i,j,1,ivx)  = -w(i,j,2,ivx) + two*BCData(nn)%uSlip(i,j,1)
                    w(i,j,1,ivy)  = -w(i,j,2,ivy) + two*BCData(nn)%uSlip(i,j,2)
                    w(i,j,1,ivz)  = -w(i,j,2,ivz) + two*BCData(nn)%uSlip(i,j,3)
                    p(i,j,1)      =  p(i,j,2) - four*third*rhok
                    rlv(i,j,1) = rlv(i,j,2)
                    if( eddyModel ) rev(i,j,1) = -rev(i,j,2)
                 enddo
              enddo
           case (kMax)
              do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                 do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                    if( correctForK ) rhok = w(i,j,1,irho)*w(i,j,1,itu1)
                    w(i,j,ke,irho) =  w(i,j,kl,irho)
                    w(i,j,ke,ivx)  = -w(i,j,kl,ivx) + two*BCData(nn)%uSlip(i,j,1)
                    w(i,j,ke,ivy)  = -w(i,j,kl,ivy) + two*BCData(nn)%uSlip(i,j,2)
                    w(i,j,ke,ivz)  = -w(i,j,kl,ivz) + two*BCData(nn)%uSlip(i,j,3)
                    p(i,j,ke)      =  p(i,j,kl) - four*third*rhok
                    rlv(i,j,ke) = rlv(i,j,jl)
                    if( eddyModel ) rev(i,j,ke) = -rev(i,j,jl)
                 enddo
              enddo
           end select

           ! Compute the energy for these halo's.

           call computeEtot(icBeg(nn),icEnd(nn), jcBeg(nn),jcEnd(nn), &
                            kcBeg(nn),kcEnd(nn), correctForK)

           ! Extrapolate the state vectors in case a second halo
           ! is needed.

           if( secondHalo ) call extrapolate2ndHalo(nn, correctForK)

         endif adiabaticWall
       enddo bocos

       end subroutine bcNSWallAdiabatic
