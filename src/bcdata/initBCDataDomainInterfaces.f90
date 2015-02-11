!
!      ******************************************************************
!      *                                                                *
!      * File:          initBCDataDomainInterfaces.f90                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-09-2005                                      *
!      * Last modified: 08-17-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initBCDataDomainInterfaces
!
!      ******************************************************************
!      *                                                                *
!      * initBCDataDomainInterfaces initializes the prescribed boundary *
!      * data for domain interfaces. It is just an initialization such  *
!      * the initial halo computations behave normally. During the      *
!      * actual computation this data should be renewed constantly by   *
!      * the coupler.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use flowVarRefState
       use inputIteration
       use inputTimeSpectral
       use BCRoutines
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, nn, mm, l, sps

       real(kind=realType) :: rho, vvx, vvy, vvz, pres, vtotInv

!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of spectral solutions and blocks of
       ! the multigrid start level.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set the pointers for this block.

           call setPointers(nn, mgStartlevel, sps)

           ! Loop over the boundary condition subfaces of this block.

           bocos: do mm=1,nBocos

             ! Check for interface boundary condition types.

             select case (BCType(mm))

               case (DomainInterfaceAll)

                 ! All data must be prescribed. Nullify the pointers
                 ! and set them to the correct subface.

                  call setBCPointers(mm, .False.)

                 ! Loop over the generic subface and simply extrapolate
                 ! the state vector to set the prescribed state.

                 do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                     BCData(mm)%rho(i,j)  = ww2(i,j,irho)
                     BCData(mm)%velx(i,j) = ww2(i,j,ivx)
                     BCData(mm)%vely(i,j) = ww2(i,j,ivy)
                     BCData(mm)%velz(i,j) = ww2(i,j,ivz)
                     BCData(mm)%ps(i,j)   = pp2(i,j)

                     do l=nt1,nt2
                       BCData(mm)%turbInlet(i,j,l) = ww2(i,j,l)
                     enddo
                   enddo
                 enddo

               !=========================================================

               case (DomainInterfaceRhoUVW)

                 ! Density, velocity and turbulent variables are
                 ! prescribed. Nullify the pointers and set them to the
                 ! correct subface.

                 call setBCPointers(mm, .False.)

                 ! Loop over the generic subface and simply extrapolate
                 ! the state vector to set the prescribed state.

                 do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                     BCData(mm)%rho(i,j)  = ww2(i,j,irho)
                     BCData(mm)%velx(i,j) = ww2(i,j,ivx)
                     BCData(mm)%vely(i,j) = ww2(i,j,ivy)
                     BCData(mm)%velz(i,j) = ww2(i,j,ivz)

                     do l=nt1,nt2
                       BCData(mm)%turbInlet(i,j,l) = ww2(i,j,l)
                     enddo
                   enddo
                 enddo

               !=========================================================

               case (DomainInterfaceP)

                 ! Pressure must be prescribed. Nullify the pointers
                 ! and set them to the correct subface.

                 call setBCPointers(mm, .False.)

                 ! Loop over the generic subface and simply extrapolate
                 ! the pressure to set the prescribed value.

                 do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                     BCData(mm)%ps(i,j) = pp2(i,j)
                   enddo
                 enddo

               !=========================================================

               case (DomainInterfaceRho)

                 ! Density must be prescribed. Nullify the pointers
                 ! and set them to the correct subface.

                  call setBCPointers(mm, .False.)

                 ! Loop over the generic subface and simply extrapolate
                 ! the density to set the prescribed value.

                 do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd
                     BCData(mm)%rho(i,j) = ww2(i,j,irho)
                   enddo
                 enddo

               !=========================================================

               case (DomainInterfaceTotal)

                 ! Total conditions must be prescribed. Nullify the
                 ! pointers and set them to the correct subface.

                  call setBCPointers(mm, .False.)

                 ! Loop over the generic subface and simply extrapolate
                 ! the total conditions and the turbulence variables
                 ! to set the prescribed value.

                 do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                     ! Store the variables a bit easier.

                     rho  = ww2(i,j,irho)
                     vvx   = ww2(i,j,ivx)
                     vvy   = ww2(i,j,ivy)
                     vvz   = ww2(i,j,ivz)
                     pres = pp2(i,j)

                     ! Compute the total pressure, total temperature
                     ! and total entahlpy.

                     call computePtot(rho, vvx, vvy, vvz, pres, &
                                      BCData(mm)%ptInlet(i,j), 1_intType)
                     call computeTtot(rho, vvx, vvy, vvz, pres, &
                                      BCData(mm)%ttInlet(i,j), 1_intType)
                     BCData(mm)%htInlet(i,j) = (ww2(i,j,irhoE) + pres) &
                                             / rho

                     ! Determine the velocity direction.

                     vtotInv = one/max(eps,sqrt(vvx*vvx + vvy*vvy + vvz*vvz))
                     BCData(mm)%flowXdirInlet(i,j) = vvx*vtotInv
                     BCData(mm)%flowYdirInlet(i,j) = vvy*vtotInv
                     BCData(mm)%flowZdirInlet(i,j) = vvz*vtotInv

                     ! Simply extrapolate the turbulence variables.

                     do l=nt1,nt2
                       BCData(mm)%turbInlet(i,j,l) = ww2(i,j,l)
                     enddo
                   enddo
                 enddo

             end select

           enddo bocos
         enddo domains
       enddo spectralLoop

       end subroutine initBCDataDomainInterfaces
