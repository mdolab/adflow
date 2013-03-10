!
!      ******************************************************************
!      *                                                                *
!      * File:          utauWF.f90                                      *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 10-01-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine utauWF(rFilv)
!
!      ******************************************************************
!      *                                                                *
!      * utauWF substitutes the wall shear stress with values from a    *
!      * look-up table, if desired.                                     *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use flowVarRefState
       use inputPhysics
       implicit none
!
!      Subroutine argument.
!
       real(kind=realType), intent(in) :: rFilv
!
!      Local variables.
!
       integer(kind=intType) :: i, j, nn

       real(kind=realType) :: fact
       real(kind=realType) :: tauxx, tauyy, tauzz
       real(kind=realType) :: tauxy, tauxz, tauyz
       real(kind=realType) :: rbar, ubar, vbar, wbar, vx, vy, vz
       real(kind=realType) :: fmx, fmy, fmz, frhoe
       real(kind=realType) :: veln, velnx, velny, velnz, tx, ty, tz
       real(kind=realType) :: veltx, velty, veltz, veltmag
       real(kind=realType) :: txnx, txny, txnz, tynx, tyny, tynz
       real(kind=realType) :: tznx, tzny, tznz
       real(kind=realType) :: tautn, tauWall, utau, re

       real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
       real(kind=realType), dimension(:,:,:), pointer :: ss, rres
       real(kind=realType), dimension(:,:,:), pointer :: norm
       real(kind=realType), dimension(:,:),   pointer :: rrlv2, dd2Wall2
!
!      Function definition.
!
       real(kind=realType) :: curveUpRe
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if no wall functions must be used.

       if(.not. wallFunctions) return

       ! Loop over the viscous subfaces of this block.

       viscSubfaces: do nn=1,nViscBocos

         ! Set a bunch of variables depending on the face id to make
         ! a generic treatment possible.

         select case (BCFaceID(nn))

           case (iMin)
             fact = -one

             ss  => si(1,:,:,:);  rres => fw(2,1:,1:,:)
             ww2 => w(2,1:,1:,:); ww1  => w(1,1:,1:,:)
             dd2Wall2 => d2Wall(2,:,:); rrlv2 => rlv(2,1:,1:)

           !===========================================================

           case (iMax)
             fact = one

             ss  => si(il,:,:,:);  rres => fw(il,1:,1:,:)
             ww2 => w(il,1:,1:,:); ww1  => w(ie,1:,1:,:)
             dd2Wall2 => d2Wall(il,:,:); rrlv2 => rlv(il,1:,1:)

           !===========================================================

           case (jMin)
             fact = -one

             ss  => sj(:,1,:,:);  rres => fw(1:,2,1:,:)
             ww2 => w(1:,2,1:,:); ww1  => w(1:,1,1:,:)
             dd2Wall2 => d2Wall(:,2,:); rrlv2 => rlv(1:,2,1:)

           !===========================================================

           case (jMax)
             fact = one

             ss  => sj(:,jl,:,:);  rres => fw(1:,jl,1:,:)
             ww2 => w(1:,jl,1:,:); ww1  => w(1:,je,1:,:)
             dd2Wall2 => d2Wall(:,jl,:); rrlv2 => rlv(1:,jl,1:)

           !===========================================================

           case (kMin)
             fact = -one

             ss  => sk(:,:,1,:);  rres => fw(1:,1:,2,:)
             ww2 => w(1:,1:,2,:); ww1  => w(1:,1:,1,:)
             dd2Wall2 => d2Wall(:,:,2); rrlv2 => rlv(1:,1:,2)

           !===========================================================

           case (kMax)
             fact = one

             ss  => sk(:,:,kl,:);  rres => fw(1:,1:,kl,:)
             ww2 => w(1:,1:,kl,:); ww1  => w(1:,1:,ke,:)
             dd2Wall2 => d2Wall(:,:,kl); rrlv2 => rlv(1:,1:,kl)

         end select

         ! Set the pointer for the unit outward normals.

         norm => BCData(nn)%norm

         ! Loop over the quadrilateral faces of the subface. Note
         ! that the nodal range of BCData must be used and not the
         ! cell range, because the latter may include the halo's in i
         ! and j-direction. The offset +1 is there, because inBeg and
         ! jnBeg refer to nodal ranges and not to cell ranges.

         do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
           do i=(BCData(nn)%inBeg+1),BCData(nn)%inEnd

             ! Store the viscous stress tensor a bit easier.

             tauxx = viscSubface(nn)%tau(i,j,1)
             tauyy = viscSubface(nn)%tau(i,j,2)
             tauzz = viscSubface(nn)%tau(i,j,3)
             tauxy = viscSubface(nn)%tau(i,j,4)
             tauxz = viscSubface(nn)%tau(i,j,5)
             tauyz = viscSubface(nn)%tau(i,j,6)

             ! Compute the velocities at the wall face; these are only
             ! non-zero for moving a block. Also compute the density,
             ! which is needed to compute the wall shear stress via
             ! wall functions.

             rbar = half*(ww2(i,j,irho) + ww1(i,j,irho))
             ubar = half*(ww2(i,j,ivx)  + ww1(i,j,ivx))
             vbar = half*(ww2(i,j,ivy)  + ww1(i,j,ivy))
             wbar = half*(ww2(i,j,ivz)  + ww1(i,j,ivz))

             ! Compute the velocity difference between the internal cell
             ! and the wall.

             vx = ww2(i,j,ivx) - ubar
             vy = ww2(i,j,ivy) - vbar
             vz = ww2(i,j,ivz) - wbar

             ! Compute the normal velocity of the internal cell.

             veln  = vx*norm(i,j,1) + vy*norm(i,j,2) + vz*norm(i,j,3)
             velnx = veln*norm(i,j,1)
             velny = veln*norm(i,j,2)
             velnz = veln*norm(i,j,3)

             ! Compute the tangential velocity, its magnitude and its
             ! unit vector of the internal cell.

             veltx = vx - velnx
             velty = vy - velny
             veltz = vz - velnz

             veltmag = max(eps,sqrt(veltx**2 + velty**2 + veltz**2))

             tx = veltx/veltmag
             ty = velty/veltmag
             tz = veltz/veltmag

             ! Compute some coefficients needed for the transformation
             ! between the cartesian frame and the frame defined by the
             ! tangential direction (tx,ty,tz) and the normal direction.
             ! The minus sign is present, because for this transformation
             ! the normal direction should be inward pointing and norm
             ! is outward pointing.

             txnx = -tx*norm(i,j,1)
             txny = -tx*norm(i,j,2)
             txnz = -tx*norm(i,j,3)

             tynx = -ty*norm(i,j,1)
             tyny = -ty*norm(i,j,2)
             tynz = -ty*norm(i,j,3)

             tznx = -tz*norm(i,j,1)
             tzny = -tz*norm(i,j,2)
             tznz = -tz*norm(i,j,3)

             ! Compute the tn component of the wall shear stress
             ! tensor. Normally this is the only nonzero shear
             ! stress component in the t-n frame.

             tautn = tauxx*txnx + tauyy*tyny + tauzz*tznz &
                   + tauxy*(txny + tynx)                  &
                   + tauxz*(txnz + tznx)                  &
                   + tauyz*(tynz + tzny)

             ! Compute the Reynolds number using the velocity, density,
             ! laminar viscosity and wall distance. Note that an offset
             ! of -1 must be used in dd2Wall2, because the original array
             ! d2Wall starts at 2.

             re = ww2(i,j,irho)*veltmag*dd2Wall2(i-1,j-1)/rrlv2(i,j)

             ! Determine the friction velocity from the table and
             ! compute the wall shear stress from it.

             utau    = veltmag/max(curveUpRe(re),eps)
             tauWall = rbar*utau*utau

             ! Compute the correction to the wall shear stress tautn and
             ! transform this correction back to the cartesian frame.
             ! Take rFilv into account, such that the correction to the
             ! stress tensor is computed correctly.

             tautn = rFilv*tauWall - tautn

             tauxx = two*tautn*txnx
             tauyy = two*tautn*tyny
             tauzz = two*tautn*tznz

             tauxy = tautn*(txny + tynx)
             tauxz = tautn*(txnz + tznx)
             tauyz = tautn*(tynz + tzny)

             ! Compute the correction to the viscous flux at the wall.

             fmx   = tauxx*ss(i,j,1) + tauxy*ss(i,j,2) &
                   + tauxz*ss(i,j,3)
             fmy   = tauxy*ss(i,j,1) + tauyy*ss(i,j,2) &
                   + tauyz*ss(i,j,3)
             fmz   = tauxz*ss(i,j,1) + tauyz*ss(i,j,2) &
                   + tauzz*ss(i,j,3)
             frhoE = (ubar*tauxx + vbar*tauxy + wbar*tauxz)*ss(i,j,1) &
                   + (ubar*tauxy + vbar*tauyy + wbar*tauyz)*ss(i,j,2) &
                   + (ubar*tauxz + vbar*tauyz + wbar*tauzz)*ss(i,j,3)

             ! Add them to the residual. Note that now the factor rFilv
             ! is already taken into account via tau. Fact is present to
             ! take inward/outward pointing normals into account

             rres(i,j,imx)   = rres(i,j,imx)   - fact*fmx
             rres(i,j,imy)   = rres(i,j,imy)   - fact*fmy
             rres(i,j,imz)   = rres(i,j,imz)   - fact*fmz
             rres(i,j,irhoE) = rres(i,j,irhoE) - fact*frhoE

             ! Store the friction velocity for later use.

             viscSubface(nn)%utau(i,j) = utau

             ! Also add the correction to the wall stress tensor.

             viscSubface(nn)%tau(i,j,1) = &
                                  viscSubface(nn)%tau(i,j,1) + tauxx
             viscSubface(nn)%tau(i,j,2) = &
                                  viscSubface(nn)%tau(i,j,2) + tauyy
             viscSubface(nn)%tau(i,j,3) = &
                                  viscSubface(nn)%tau(i,j,3) + tauzz
             viscSubface(nn)%tau(i,j,4) = &
                                  viscSubface(nn)%tau(i,j,4) + tauxy
             viscSubface(nn)%tau(i,j,5) = &
                                  viscSubface(nn)%tau(i,j,5) + tauxz
             viscSubface(nn)%tau(i,j,6) = &
                                  viscSubface(nn)%tau(i,j,6) + tauyz
           enddo
         enddo

       enddo viscSubfaces

       end subroutine utauWF
