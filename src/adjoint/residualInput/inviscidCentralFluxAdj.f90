!
!      ******************************************************************
!      *                                                                *
!      * File:          inviscidCentralFluxAdj.f90                      *
!      * Author:        Edwin van der Weide, C.A.(Sandy) Mader          *
!      *                Seongim Choi                                    *
!      * Starting date: 11-21-2007                                      *
!      * Last modified: 10-22-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine inviscidCentralFluxAdj(wAdj,  pAdj,  dwAdj,         &
                                         siAdj, sjAdj, skAdj, volAdj, &
                                         sFaceIAdj,sFaceJAdj,sFaceKAdj,&
                                         rotRateAdj,                  &
                                         iCell, jCell, kCell,nn,level,sps)
!
!      ******************************************************************
!      *                                                                *
!      * inviscidCentralFluxAdj computes the Euler fluxes using a       *
!      * central discretization for the cell iCell, jCell, kCell of the *
!      * block to which the variables in blockPointers currently point  *
!      * to.                                                            *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers   ! sFaceI,sFaceJ,sFaceK,sI,sJ,sK,blockismoving,addgridvelocities
                           ! vol, nbkGlobal
       use flowVarRefState ! constants (irho, ivx, ivy, imx,..), timeRef
       use inputPhysics    ! equationMode, steady
       use cgnsGrid        !
       use inputTimeSpectral !nTimeIntervalsSpectral
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType) :: iCell, jCell, kCell,nn,level,sps

       real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral), &
                                                      intent(in) :: wAdj
       real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral),    &
                                                      intent(in) :: pAdj
       real(kind=realType), dimension(nw,nTimeIntervalsSpectral), intent(inout) :: dwAdj
       real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral), intent(in) :: siAdj, sjAdj, skAdj
       real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral), intent(in) ::sFaceIAdj,sFaceJAdj,sFaceKAdj
       real(kind=realType), dimension(0:0,0:0,0:0,nTimeIntervalsSpectral), intent(in) :: volAdj
       real(kind=realType), dimension(3),intent(in) ::rotRateAdj
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ii, jj, kk

       real(kind=realType) :: qsp, qsm, rqsp, rqsm, porVel, porFlux
       real(kind=realType) :: pa, fs, sFace, vnp, vnm, fact
       real(kind=realType) :: wx, wy, wz, rvol

!     testing vars
       real(kind=realType) :: wx2, wy2, wz2, rvol2

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!

       ! Initialize sFace to zero. This value will be used if the
       ! block is not moving.

       sFace = 0.0
!
!      ******************************************************************
!      *                                                                *
!      * Advective fluxes in the i-direction.                           *
!      *                                                                *
!      ******************************************************************
!
       i    = iCell-1; j = jCell; k = kCell
       
       fact = -one

       ! Loop over the two faces which contribute to the residual of
       ! the cell considered.

       do ii=-1,0

         ! Set the dot product of the grid velocity and the
         ! normal in i-direction for a moving face.
          
         if( addGridVelocities ) sFace = sFaceIAdj(ii,0,0,sps)

         ! Compute the normal velocities of the left and right state.

         vnp = wAdj(ii+1,0,0,ivx,sps)*sIAdj(ii,0,0,1,sps) &
             + wAdj(ii+1,0,0,ivy,sps)*sIAdj(ii,0,0,2,sps) &
             + wAdj(ii+1,0,0,ivz,sps)*sIAdj(ii,0,0,3,sps)
         vnm = wAdj(ii,  0,0,ivx,sps)*sIAdj(ii,0,0,1,sps) &
             + wAdj(ii,  0,0,ivy,sps)*sIAdj(ii,0,0,2,sps) &
             + wAdj(ii,  0,0,ivz,sps)*sIAdj(ii,0,0,3,sps)
         !print *,'vnp',wAdj(ii+1,0,0,ivx,sps),sIAdj(ii,0,0,1,sps),sps

         ! Set the values of the porosities for this face.
         ! porVel defines the porosity w.r.t. velocity;
         ! porFlux defines the porosity w.r.t. the entire flux.
         ! The latter is only zero for a discontinuous block
         ! boundary that must be treated conservatively.
         ! The default value of porFlux is 0.5, such that the
         ! correct central flux is scattered to both cells.
         ! In case of a boundFlux the normal velocity is set
         ! to sFace.

         porVel  = one
         porFlux = half
         if(porI(i,j,k) == noFlux)    porFlux = 0.0
         if(porI(i,j,k) == boundFlux) then
            porVel = 0.0
            vnp    = sFace
            vnm    = sFace
         endif
          
         ! Incorporate porFlux in porVel.

         porVel = porVel*porFlux

         ! Compute the normal velocities relative to the grid for
         ! the face as well as the mass fluxes.

         qsp = (vnp -sFace)*porVel
         qsm = (vnm -sFace)*porVel

         rqsp = qsp*wAdj(ii+1,0,0,irho,sps)
         rqsm = qsm*wAdj(ii,  0,0,irho,sps)
 
         ! Compute the sum of the pressure multiplied by porFlux.
         ! For the default value of porFlux, 0.5, this leads to
         ! the average pressure.
          
         pa = porFlux*(pAdj(ii+1,0,0,sps) + pAdj(ii,0,0,sps))
         
         ! Compute the fluxes through this face.
          
         fs = rqsp + rqsm
          
         dwAdj(irho,sps) = dwAdj(irho,sps) + fact*fs
          
         fs = rqsp*wAdj(ii+1,0,0,ivx,sps) + rqsm*wAdj(ii,0,0,ivx,sps) &
            + pa*sIAdj(ii,0,0,1,sps)
         dwAdj(imx,sps) = dwAdj(imx,sps) + fact*fs
          
         fs = rqsp*wAdj(ii+1,0,0,ivy,sps) + rqsm*wAdj(ii,0,0,ivy,sps) &
            + pa*sIAdj(ii,0,0,2,sps)
         dwAdj(imy,sps) = dwAdj(imy,sps) + fact*fs
          
         fs = rqsp*wAdj(ii+1,0,0,ivz,sps) + rqsm*wAdj(ii,0,0,ivz,sps) &
            + pa*sIAdj(ii,0,0,3,sps)
         dwAdj(imz,sps) = dwAdj(imz,sps) + fact*fs
          
         fs = qsp*wAdj(ii+1,0,0,irhoE,sps) + qsm*wAdj(ii,0,0,irhoE,sps) &
            + porFlux*(vnp*pAdj(ii+1,0,0,sps) + vnm*pAdj(ii,0,0,sps))
         dwAdj(irhoE,sps) = dwAdj(irhoE,sps) + fact*fs

         ! Update i and set fact to 1 for the second face.

         i    = i + 1
         fact = one

       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Advective fluxes in the j-direction.                           *
!      *                                                                *
!      ******************************************************************
!
       i    = iCell; j = jCell-1; k = kCell
       fact = -one

       ! Loop over the two faces which contribute to the residual of
       ! the cell considered.

       do jj=-1,0

         ! Set the dot product of the grid velocity and the
         ! normal in j-direction for a moving face.
          
         if( addGridVelocities ) sFace = sFaceJAdj(0,jj,0,sps)

         ! Compute the normal velocities of the left and right state.

         vnp = wAdj(0,jj+1,0,ivx,sps)*sJAdj(0,jj,0,1,sps) &
             + wAdj(0,jj+1,0,ivy,sps)*sJAdj(0,jj,0,2,sps) &
             + wAdj(0,jj+1,0,ivz,sps)*sJAdj(0,jj,0,3,sps)
         vnm = wAdj(0,jj,  0,ivx,sps)*sJAdj(0,jj,0,1,sps) &
             + wAdj(0,jj,  0,ivy,sps)*sJAdj(0,jj,0,2,sps) &
             + wAdj(0,jj,  0,ivz,sps)*sJAdj(0,jj,0,3,sps)

         ! Set the values of the porosities for this face.
         ! porVel defines the porosity w.r.t. velocity;
         ! porFlux defines the porosity w.r.t. the entire flux.
         ! The latter is only zero for a discontinuous block
         ! boundary that must be treated conservatively.
         ! The default value of porFlux is 0.5, such that the
         ! correct central flux is scattered to both cells.
         ! In case of a boundFlux the normal velocity is set
         ! to sFace.

         porVel  = one
         porFlux = half
         if(porJ(i,j,k) == noFlux)    porFlux = 0.0
         if(porJ(i,j,k) == boundFlux) then
            porVel = 0.0
            vnp    = sFace
            vnm    = sFace
         endif
          
         ! Incorporate porFlux in porVel.

         porVel = porVel*porFlux

         ! Compute the normal velocities relative to the grid for
         ! the face as well as the mass fluxes.

         qsp = (vnp - sFace)*porVel
         qsm = (vnm - sFace)*porVel

         rqsp = qsp*wAdj(0,jj+1,0,irho,sps)
         rqsm = qsm*wAdj(0,jj,  0,irho,sps)
 
         ! Compute the sum of the pressure multiplied by porFlux.
         ! For the default value of porFlux, 0.5, this leads to
         ! the average pressure.
          
         pa = porFlux*(pAdj(0,jj+1,0,sps) + pAdj(0,jj,0,sps))
          
         ! Compute the fluxes through this face.
          
         fs = rqsp + rqsm
          
         dwAdj(irho,sps) = dwAdj(irho,sps) + fact*fs
          
         fs = rqsp*wAdj(0,jj+1,0,ivx,sps) + rqsm*wAdj(0,jj,0,ivx,sps) &
            + pa*sJAdj(0,jj,0,1,sps)
         dwAdj(imx,sps) = dwAdj(imx,sps) + fact*fs
          
         fs = rqsp*wAdj(0,jj+1,0,ivy,sps) + rqsm*wAdj(0,jj,0,ivy,sps) &
            + pa*sJAdj(0,jj,0,2,sps)
         dwAdj(imy,sps) = dwAdj(imy,sps) + fact*fs
          
         fs = rqsp*wAdj(0,jj+1,0,ivz,sps) + rqsm*wAdj(0,jj,0,ivz,sps) &
            + pa*sJAdj(0,jj,0,3,sps)
         dwAdj(imz,sps) = dwAdj(imz,sps) + fact*fs
          
         fs = qsp*wAdj(0,jj+1,0,irhoE,sps) + qsm*wAdj(0,jj,0,irhoE,sps) &
            + porFlux*(vnp*pAdj(0,jj+1,0,sps) + vnm*pAdj(0,jj,0,sps))
         dwAdj(irhoE,sps) = dwAdj(irhoE,sps) + fact*fs

         ! Update j and set fact to 1 for the second face.

         j    = j + 1
         fact = one
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Advective fluxes in the k-direction.                           *
!      *                                                                *
!      ******************************************************************
!
!       should this be inside, have k=kCell+kk?
       i    = iCell; j = jCell; k = kCell-1
       fact = -one

       ! Loop over the two faces which contribute to the residual of
       ! the cell considered.

       do kk=-1,0

         ! Set the dot product of the grid velocity and the
         ! normal in k-direction for a moving face.
          
         if( addGridVelocities ) sFace = sFaceKAdj(0,0,kk,sps)

         ! Compute the normal velocities of the left and right state.

         vnp = wAdj(0,0,kk+1,ivx,sps)*sKAdj(0,0,kk,1,sps) &
             + wAdj(0,0,kk+1,ivy,sps)*sKAdj(0,0,kk,2,sps) &
             + wAdj(0,0,kk+1,ivz,sps)*sKAdj(0,0,kk,3,sps)
         vnm = wAdj(0,0,kk,  ivx,sps)*sKAdj(0,0,kk,1,sps) &
             + wAdj(0,0,kk,  ivy,sps)*sKAdj(0,0,kk,2,sps) &
             + wAdj(0,0,kk,  ivz,sps)*sKAdj(0,0,kk,3,sps)

         ! Set the values of the porosities for this face.
         ! porVel defines the porosity w.r.t. velocity;
         ! porFlux defines the porosity w.r.t. the entire flux.
         ! The latter is only zero for a discontinuous block
         ! boundary that must be treated conservatively.
         ! The default value of porFlux is 0.5, such that the
         ! correct central flux is scattered to both cells.
         ! In case of a boundFlux the normal velocity is set
         ! to sFace.

         porVel  = one
         porFlux = half
         if(porK(i,j,k) == noFlux)    porFlux = 0.0
         if(porK(i,j,k) == boundFlux) then
            porVel = 0.0
            vnp    = sFace
            vnm    = sFace
         endif
          
         ! Incorporate porFlux in porVel.

         porVel = porVel*porFlux

         ! Compute the normal velocities relative to the grid for
         ! the face as well as the mass fluxes.

         qsp = (vnp - sFace)*porVel
         qsm = (vnm - sFace)*porVel

         rqsp = qsp*wAdj(0,0,kk+1,irho,sps)
         rqsm = qsm*wAdj(0,0,kk,  irho,sps)
 
         ! Compute the sum of the pressure multiplied by porFlux.
         ! For the default value of porFlux, 0.5, this leads to
         ! the average pressure.
          
         pa = porFlux*(pAdj(0,0,kk+1,sps) + pAdj(0,0,kk,sps))
          
         ! Compute the fluxes through this face.
          
         fs = rqsp + rqsm
          
         dwAdj(irho,sps) = dwAdj(irho,sps) + fact*fs
          
         fs = rqsp*wAdj(0,0,kk+1,ivx,sps) + rqsm*wAdj(0,0,kk,ivx,sps) &
            + pa*sKAdj(0,0,kk,1,sps)
         dwAdj(imx,sps) = dwAdj(imx,sps) + fact*fs
          
         fs = rqsp*wAdj(0,0,kk+1,ivy,sps) + rqsm*wAdj(0,0,kk,ivy,sps) &
            + pa*sKAdj(0,0,kk,2,sps)
         dwAdj(imy,sps) = dwAdj(imy,sps) + fact*fs
          
         fs = rqsp*wAdj(0,0,kk+1,ivz,sps) + rqsm*wAdj(0,0,kk,ivz,sps) &
            + pa*sKAdj(0,0,kk,3,sps)
         dwAdj(imz,sps) = dwAdj(imz,sps) + fact*fs
          
         fs = qsp*wAdj(0,0,kk+1,irhoE,sps) + qsm*wAdj(0,0,kk,irhoE,sps) &
            + porFlux*(vnp*pAdj(0,0,kk+1,sps) + vnm*pAdj(0,0,kk,sps))
         dwAdj(irhoE,sps) = dwAdj(irhoE,sps) + fact*fs

         ! Update k and set fact to 1 for the second face.

         k    = k + 1
         fact = one
       enddo

       ! Add the rotational source terms for a moving block in a
       ! steady state computation. These source terms account for the
       ! centrifugal acceleration and the coriolis term. However, as
       ! the equations are solved in the inertial frame and not
       ! in the moving frame, the form is different than what you
       ! normally find in a text book.
       
       rotation: if(blockIsMoving .and. equationMode == steady) then

!          wx = timeRef*rotRateAdj(1)
!          wy = timeRef*rotRateAdj(2)
!          wz = timeRef*rotRateAdj(3)
          !timeref is taken into account in copyAdjointStencil...
          wx = rotRateAdj(1)
          wy = rotRateAdj(2)
          wz = rotRateAdj(3)

          rvol = wAdj(0,0,0,irho,sps)*volAdj(0,0,0,sps)
         
          dwAdj(imx,sps) = dwAdj(imx,sps) &
               + rvol*(wy*wAdj(0,0,0,ivz,sps) - wz*wAdj(0,0,0,ivy,sps))
          dwAdj(imy,sps) = dwAdj(imy,sps) &
               + rvol*(wz*wAdj(0,0,0,ivx,sps) - wx*wAdj(0,0,0,ivz,sps))
          dwAdj(imz,sps) = dwAdj(imz,sps) &
               + rvol*(wx*wAdj(0,0,0,ivy,sps) - wy*wAdj(0,0,0,ivx,sps))


       end if rotation

       end subroutine inviscidCentralFluxAdj
