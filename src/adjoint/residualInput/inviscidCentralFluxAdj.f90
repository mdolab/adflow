!
!      ******************************************************************
!      *                                                                *
!      * File:          inviscidCentralFluxAdj.f90                      *
!      * Author:        Edwin van der Weide                             *
!      *                Seongim Choi                                    *
!      * Starting date: 11-21-2007                                      *
!      * Last modified: 12-17-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine inviscidCentralFluxAdj(wAdj,  pAdj,  dwAdj,         &
                                         siAdj, sjAdj, skAdj, volAdj, &
                                         iCell, jCell, kCell)
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
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType) :: iCell, jCell, kCell

       real(kind=realType), dimension(-2:2,-2:2,-2:2,nw), &
                                                      intent(in) :: wAdj
       real(kind=realType), dimension(-2:2,-2:2,-2:2),    &
                                                      intent(in) :: pAdj
       real(kind=realType), dimension(nw), intent(inout) :: dwAdj
       real(kind=realType), dimension(-2:2,-2:2,-2:2,3), intent(in) :: siAdj, sjAdj, skAdj
       real(kind=realType), dimension(0:0,0:0,0:0), intent(in) :: volAdj
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ii, jj, kk

       real(kind=realType) :: qsp, qsm, rqsp, rqsm, porVel, porFlux
       real(kind=realType) :: pa, fs, sFace, vnp, vnm, fact
       real(kind=realType) :: wx, wy, wz, rvol

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
          
         if( addGridVelocities ) sFace = sFaceI(i,j,k)

         ! Compute the normal velocities of the left and right state.

         vnp = wAdj(ii+1,0,0,ivx)*sIAdj(ii,0,0,1) &
             + wAdj(ii+1,0,0,ivy)*sIAdj(ii,0,0,2) &
             + wAdj(ii+1,0,0,ivz)*sIAdj(ii,0,0,3)
         vnm = wAdj(ii,  0,0,ivx)*sIAdj(ii,0,0,1) &
             + wAdj(ii,  0,0,ivy)*sIAdj(ii,0,0,2) &
             + wAdj(ii,  0,0,ivz)*sIAdj(ii,0,0,3)


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

         qsp = (vnp - sFace)*porVel
         qsm = (vnm - sFace)*porVel

         rqsp = qsp*wAdj(ii+1,0,0,irho)
         rqsm = qsm*wAdj(ii,  0,0,irho)
 
         ! Compute the sum of the pressure multiplied by porFlux.
         ! For the default value of porFlux, 0.5, this leads to
         ! the average pressure.
          
         pa = porFlux*(pAdj(ii+1,0,0) + pAdj(ii,0,0))
         
         ! Compute the fluxes through this face.
          
         fs = rqsp + rqsm
          
         dwAdj(irho) = dwAdj(irho) + fact*fs
          
         fs = rqsp*wAdj(ii+1,0,0,ivx) + rqsm*wAdj(ii,0,0,ivx) &
            + pa*sIAdj(ii,0,0,1)
         dwAdj(imx) = dwAdj(imx) + fact*fs
          
         fs = rqsp*wAdj(ii+1,0,0,ivy) + rqsm*wAdj(ii,0,0,ivy) &
            + pa*sIAdj(ii,0,0,2)
         dwAdj(imy) = dwAdj(imy) + fact*fs
          
         fs = rqsp*wAdj(ii+1,0,0,ivz) + rqsm*wAdj(ii,0,0,ivz) &
            + pa*sIAdj(ii,0,0,3)
         dwAdj(imz) = dwAdj(imz) + fact*fs
          
         fs = qsp*wAdj(ii+1,0,0,irhoE) + qsm*wAdj(ii,0,0,irhoE) &
            + porFlux*(vnp*pAdj(ii+1,0,0) + vnm*pAdj(ii,0,0))
         dwAdj(irhoE) = dwAdj(irhoE) + fact*fs

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
          
         if( addGridVelocities ) sFace = sFaceJ(i,j,k)

         ! Compute the normal velocities of the left and right state.

         vnp = wAdj(0,jj+1,0,ivx)*sJAdj(0,jj,0,1) &
             + wAdj(0,jj+1,0,ivy)*sJAdj(0,jj,0,2) &
             + wAdj(0,jj+1,0,ivz)*sJAdj(0,jj,0,3)
         vnm = wAdj(0,jj,  0,ivx)*sJAdj(0,jj,0,1) &
             + wAdj(0,jj,  0,ivy)*sJAdj(0,jj,0,2) &
             + wAdj(0,jj,  0,ivz)*sJAdj(0,jj,0,3)

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

         rqsp = qsp*wAdj(0,jj+1,0,irho)
         rqsm = qsm*wAdj(0,jj,  0,irho)
 
         ! Compute the sum of the pressure multiplied by porFlux.
         ! For the default value of porFlux, 0.5, this leads to
         ! the average pressure.
          
         pa = porFlux*(pAdj(0,jj+1,0) + pAdj(0,jj,0))
          
         ! Compute the fluxes through this face.
          
         fs = rqsp + rqsm
          
         dwAdj(irho) = dwAdj(irho) + fact*fs
          
         fs = rqsp*wAdj(0,jj+1,0,ivx) + rqsm*wAdj(0,jj,0,ivx) &
            + pa*sJAdj(0,jj,0,1)
         dwAdj(imx) = dwAdj(imx) + fact*fs
          
         fs = rqsp*wAdj(0,jj+1,0,ivy) + rqsm*wAdj(0,jj,0,ivy) &
            + pa*sJAdj(0,jj,0,2)
         dwAdj(imy) = dwAdj(imy) + fact*fs
          
         fs = rqsp*wAdj(0,jj+1,0,ivz) + rqsm*wAdj(0,jj,0,ivz) &
            + pa*sJAdj(0,jj,0,3)
         dwAdj(imz) = dwAdj(imz) + fact*fs
          
         fs = qsp*wAdj(0,jj+1,0,irhoE) + qsm*wAdj(0,jj,0,irhoE) &
            + porFlux*(vnp*pAdj(0,jj+1,0) + vnm*pAdj(0,jj,0))
         dwAdj(irhoE) = dwAdj(irhoE) + fact*fs

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
          
         if( addGridVelocities ) sFace = sFaceK(i,j,k)

         ! Compute the normal velocities of the left and right state.

         vnp = wAdj(0,0,kk+1,ivx)*sKAdj(0,0,kk,1) &
             + wAdj(0,0,kk+1,ivy)*sKAdj(0,0,kk,2) &
             + wAdj(0,0,kk+1,ivz)*sKAdj(0,0,kk,3)
         vnm = wAdj(0,0,kk,  ivx)*sKAdj(0,0,kk,1) &
             + wAdj(0,0,kk,  ivy)*sKAdj(0,0,kk,2) &
             + wAdj(0,0,kk,  ivz)*sKAdj(0,0,kk,3)

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

         rqsp = qsp*wAdj(0,0,kk+1,irho)
         rqsm = qsm*wAdj(0,0,kk,  irho)
 
         ! Compute the sum of the pressure multiplied by porFlux.
         ! For the default value of porFlux, 0.5, this leads to
         ! the average pressure.
          
         pa = porFlux*(pAdj(0,0,kk+1) + pAdj(0,0,kk))
          
         ! Compute the fluxes through this face.
          
         fs = rqsp + rqsm
          
         dwAdj(irho) = dwAdj(irho) + fact*fs
          
         fs = rqsp*wAdj(0,0,kk+1,ivx) + rqsm*wAdj(0,0,kk,ivx) &
            + pa*sKAdj(0,0,kk,1)
         dwAdj(imx) = dwAdj(imx) + fact*fs
          
         fs = rqsp*wAdj(0,0,kk+1,ivy) + rqsm*wAdj(0,0,kk,ivy) &
            + pa*sKAdj(0,0,kk,2)
         dwAdj(imy) = dwAdj(imy) + fact*fs
          
         fs = rqsp*wAdj(0,0,kk+1,ivz) + rqsm*wAdj(0,0,kk,ivz) &
            + pa*sKAdj(0,0,kk,3)
         dwAdj(imz) = dwAdj(imz) + fact*fs
          
         fs = qsp*wAdj(0,0,kk+1,irhoE) + qsm*wAdj(0,0,kk,irhoE) &
            + porFlux*(vnp*pAdj(0,0,kk+1) + vnm*pAdj(0,0,kk))
         dwAdj(irhoE) = dwAdj(irhoE) + fact*fs

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

          wx = timeRef*cgnsDoms(nbkGlobal)%rotRate(1)
          wy = timeRef*cgnsDoms(nbkGlobal)%rotRate(2)
          wz = timeRef*cgnsDoms(nbkGlobal)%rotRate(3)

          rvol = wAdj(0,0,0,irho)*volAdj(0,0,0)
          
          dwAdj(imx) = dwAdj(imx) &
               + rvol*(wy*wAdj(0,0,0,ivz) - wz*wAdj(0,0,0,ivy))
          dwAdj(imy) = dwAdj(imy) &
               + rvol*(wz*wAdj(0,0,0,ivx) - wx*wAdj(0,0,0,ivz))
          dwAdj(imz) = dwAdj(imz) &
               + rvol*(wx*wAdj(0,0,0,ivy) - wy*wAdj(0,0,0,ivx))


       end if rotation

       end subroutine inviscidCentralFluxAdj
