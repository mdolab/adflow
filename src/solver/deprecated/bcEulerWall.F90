!
!      ******************************************************************
!      *                                                                *
!      * File:          bcEulerWall.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-07-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcEulerWall(secondHalo, correctForK)
!
!      ******************************************************************
!      *                                                                *
!      * bcEulerWall applies the inviscid wall boundary condition to    *
!      * a block. It is assumed that the pointers in blockPointers are  *
!      * already set to the correct block on the correct grid level.    *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
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
!      Local variables.
!
       integer(kind=intType) :: nn, j, k, l
       integer(kind=intType) :: jm1, jp1, km1, kp1
       integer(kind=intType) :: wallTreatment

       real(kind=realType) :: sixa, siya, siza, sjxa, sjya, sjza
       real(kind=realType) :: skxa, skya, skza, a1, b1
       real(kind=realType) :: rxj, ryj, rzj, rxk, ryk, rzk
       real(kind=realType) :: dpj, dpk, ri, rj, rk, qj, qk, vn
       real(kind=realType) :: uux, uuy, uuz

#ifndef TAPENADE_REVERSE
       real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
       real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
       real(kind=realType), dimension(:,:),   pointer :: pp3, pp4
       real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
       real(kind=realType), dimension(:,:),   pointer :: rev1, rev2
       real(kind=realType), dimension(:,:,:), pointer :: ssi, ssj, ssk
       real(kind=realType), dimension(:,:,:), pointer :: ss

!
!      Interfaces
!
       interface
         subroutine setBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                                  rev1, rev2, offset)
           use blockPointers
           use flowVarRefState
           implicit none

           integer(kind=intType), intent(in) :: nn, offset
           real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
           real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
           real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
           real(kind=realType), dimension(:,:),   pointer :: rev1, rev2
         end subroutine setBCPointers

         subroutine resetBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                                  rev1, rev2, offset)
           use blockPointers
           use flowVarRefState
           implicit none

           integer(kind=intType), intent(in) :: nn, offset
           real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
           real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
           real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
           real(kind=realType), dimension(:,:),   pointer :: rev1, rev2
         end subroutine resetBCPointers


        subroutine setpp3pp4(nn, pp3, pp4)

           use blockPointers
           implicit none

           integer(kind=intType), intent(in) :: nn
           real(kind=realType), dimension(:,:),   pointer :: pp3, pp4
         end subroutine setpp3pp4

        subroutine resetpp3pp4(nn, pp3, pp4)

           use blockPointers
           implicit none

           integer(kind=intType), intent(in) :: nn
           real(kind=realType), dimension(:,:),   pointer :: pp3, pp4
         end subroutine resetpp3pp4

        subroutine setss(nn, ssi, ssj, ssk, ss)

           use blockPointers
           implicit none

           integer(kind=intType), intent(in) :: nn
           real(kind=realType), dimension(:,:,:), pointer :: ssi, ssj, ssk
           real(kind=realType), dimension(:,:,:), pointer :: ss
         end subroutine setss

        subroutine resetss(nn, ssi, ssj, ssk, ss)

           use blockPointers
           implicit none

           integer(kind=intType), intent(in) :: nn
           real(kind=realType), dimension(:,:,:), pointer :: ssi, ssj, ssk
           real(kind=realType), dimension(:,:,:), pointer :: ss
         end subroutine resetss

       end interface

#else
       real(kind=realType), dimension(imaxDim,jmaxDim,nw) :: ww1, ww2
       real(kind=realType), dimension(imaxDim,jmaxDim) :: pp1, pp2
       real(kind=realType), dimension(imaxDim,jmaxDim) :: pp3, pp4
       real(kind=realType), dimension(imaxDim,jmaxDim) :: rlv1, rlv2
       real(kind=realType), dimension(imaxDim,jmaxDim) :: rev1, rev2
       real(kind=realType), dimension(imaxDim,jmaxDim,3) :: ssi, ssj, ssk
       real(kind=realType), dimension(imaxDim,jmaxDim,3) :: ss
#endif
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Make sure that on the coarser grids the constant pressure
       ! boundary condition is used.

       wallTreatment = eulerWallBcTreatment
       if(currentLevel > groundLevel) wallTreatment = constantPressure

       ! Loop over the boundary condition subfaces of this block.

       bocos: do nn=1,nBocos

         ! Check for Euler wall boundary condition.

         invWall: if(BCType(nn) == EulerWall) then

           ! Set the pointers for the unit normal and the normal
           ! velocity to make the code more readable.
           
           ! Modify to use actual pointer - Peter Lyu
           !norm  => BCData(nn)%norm
           !rface => BCData(nn)%rface

           ! Nullify the pointers and set them to the correct subface.
           ! They are nullified first, because some compilers require
           ! that.

           !nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)
#ifndef TAPENADE_REVERSE
           call setBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                              rev1, rev2, 0)
#else
           call setBCPointersBwd(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                rev1, rev2, 0)
#endif

!
!          **************************************************************
!          *                                                            *
!          * Determine the boundary condition treatment and compute the *
!          * undivided pressure gradient accordingly. This gradient is  *
!          * temporarily stored in the halo pressure.                   *
!          *                                                            *
!          **************************************************************
!
           BCTreatment: select case (wallTreatment)

             case (constantPressure)
               ! Constant pressure. Set the gradient to zero.

               do k=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                 do j=BCData(nn)%icBeg, BCData(nn)%icEnd
                   pp1(j,k) = zero
                 enddo
               enddo

             !===========================================================

             case (linExtrapolPressure)

               ! Linear extrapolation. First set the additional pointer
               ! for pp3, depending on the block face.

#ifndef TAPENADE_REVERSE
               call setpp3pp4(nn, pp3, pp4)
#else
               call setpp3pp4Bwd(nn, pp3, pp4)
#endif
               ! Compute the gradient.

               do k=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                 do j=BCData(nn)%icBeg, BCData(nn)%icEnd
                   pp1(j,k) = pp3(j,k) - pp2(j,k)
                 enddo
               enddo

#ifndef TAPENADE_REVERSE
               call resetpp3pp4(nn, pp3, pp4)
#else
               call resetpp3pp4Bwd(nn, pp3, pp4)
#endif


             !===========================================================

             case (quadExtrapolPressure)

               ! Quadratic extrapolation. First set the additional
               ! pointers for pp3 and pp4, depending on the block face.

#ifndef TAPENADE_REVERSE
               call setpp3pp4(nn, pp3, pp4)
#else
               call setpp3pp4Bwd(nn, pp3, pp4)
#endif

               ! Compute the gradient.

               do k=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                 do j=BCData(nn)%icBeg, BCData(nn)%icEnd
                   pp1(j,k) = two*pp3(j,k) - 1.5_realType*pp2(j,k) &
                            - half*pp4(j,k)
                 enddo
               enddo

#ifndef TAPENADE_REVERSE
               call resetpp3pp4(nn, pp3, pp4)
#else
               call resetpp3pp4Bwd(nn, pp3, pp4)
#endif

             !===========================================================

             case (normalMomentum)

               ! Pressure gradient is computed using the normal momentum
               ! equation. First set a couple of additional variables for
               ! the normals, depending on the block face. Note that the
               ! construction 1: should not be used in these pointers,
               ! because element 0 is needed. Consequently there will be
               ! an offset of 1 for these normals. This is commented in
               ! the code. For moving faces also the grid velocity of
               ! the 1st cell center from the wall is needed.

#ifndef TAPENADE_REVERSE
               call setss(nn, ssi, ssj, ssk, ss)
#else
               call setssBwd(nn, ssi, ssj, ssk, ss)
#endif


               ! Loop over the faces of the generic subface.

               do k=BCData(nn)%jcBeg, BCData(nn)%jcEnd

                 ! Store the indices k+1, k-1 a bit easier and make
                 ! sure that they do not exceed the range of the arrays.

                 km1 = k-1; km1 = max(BCData(nn)%jcBeg,km1)
                 kp1 = k+1; kp1 = min(BCData(nn)%jcEnd,kp1)

                 ! Compute the scaling factor for the central difference
                 ! in the k-direction.

                 b1 = one/max(1_intType,(kp1-km1))

                 ! The j-loop.

                 do j=BCData(nn)%icBeg, BCData(nn)%icEnd

                   ! The indices j+1 and j-1. Make sure that they
                   ! do not exceed the range of the arrays.

                   jm1 = j-1; jm1 = max(BCData(nn)%icBeg,jm1)
                   jp1 = j+1; jp1 = min(BCData(nn)%icEnd,jp1)

                   ! Compute the scaling factor for the central
                   ! difference in the j-direction.

                   a1 = one/max(1_intType,(jp1-jm1))

                   ! Compute (twice) the average normal in the generic i,
                   ! j and k-direction. Note that in j and k-direction
                   ! the average in the original indices should be taken
                   ! using j-1 and j (and k-1 and k). However due to the
                   ! usage of pointers ssj and ssk there is an offset in
                   ! the indices of 1 and therefore now the correct
                   ! average is obtained with the indices j and j+1
                   ! (k and k+1).

                   sixa = two*ssi(j,k,1)
                   siya = two*ssi(j,k,2)
                   siza = two*ssi(j,k,3)

                   sjxa = ssj(j,k,1) + ssj(j+1,k,1)
                   sjya = ssj(j,k,2) + ssj(j+1,k,2)
                   sjza = ssj(j,k,3) + ssj(j+1,k,3)

                   skxa = ssk(j,k,1) + ssk(j,k+1,1)
                   skya = ssk(j,k,2) + ssk(j,k+1,2)
                   skza = ssk(j,k,3) + ssk(j,k+1,3)

                   ! Compute the difference of the normal vector and
                   ! pressure in j and k-direction. As the indices are
                   ! restricted to the 1st halo-layer, the computation
                   ! of the internal halo values is not consistent;
                   ! however this is not really a problem, because these
                   ! values are overwritten in the communication pattern.

                   rxj = a1*(BCData(nn)%norm(jp1,k,1) - BCData(nn)%norm(jm1,k,1))
                   ryj = a1*(BCData(nn)%norm(jp1,k,2) - BCData(nn)%norm(jm1,k,2))
                   rzj = a1*(BCData(nn)%norm(jp1,k,3) - BCData(nn)%norm(jm1,k,3))
                   dpj = a1*(pp2(jp1,k)    - pp2(jm1,k))

                   rxk = b1*(BCData(nn)%norm(j,kp1,1) - BCData(nn)%norm(j,km1,1))
                   ryk = b1*(BCData(nn)%norm(j,kp1,2) - BCData(nn)%norm(j,km1,2))
                   rzk = b1*(BCData(nn)%norm(j,kp1,3) - BCData(nn)%norm(j,km1,3))
                   dpk = b1*(pp2(j,kp1)    - pp2(j,km1))

                   ! Compute the dot product between the unit vector
                   ! and the normal vectors in i, j and k-direction.

                   ri = BCData(nn)%norm(j,k,1)*sixa + BCData(nn)%norm(j,k,2)*siya &
                      + BCData(nn)%norm(j,k,3)*siza
                   rj = BCData(nn)%norm(j,k,1)*sjxa + BCData(nn)%norm(j,k,2)*sjya &
                      + BCData(nn)%norm(j,k,3)*sjza
                   rk = BCData(nn)%norm(j,k,1)*skxa + BCData(nn)%norm(j,k,2)*skya &
                      + BCData(nn)%norm(j,k,3)*skza

                   ! Store the velocity components in ux, uy and uz and
                   ! subtract the mesh velocity if the face is moving.

                   uux = ww2(j,k,ivx)
                   uuy = ww2(j,k,ivy)
                   uuz = ww2(j,k,ivz)

                   if( addGridVelocities ) then
                     uux = uux - ss(j,k,1)
                     uuy = uuy - ss(j,k,2)
                     uuz = uuz - ss(j,k,3)
                   endif

                   ! Compute the velocity components in j and
                   ! k-direction.

                   qj = uux*sjxa + uuy*sjya + uuz*sjza
                   qk = uux*skxa + uuy*skya + uuz*skza

                   ! Compute the pressure gradient, which is stored
                   ! in pp1. I'm not entirely sure whether this
                   ! formulation is correct for moving meshes. It could
                   ! be that an additional term is needed there.

                   pp1(j,k) = ((qj*(uux*rxj + uuy*ryj + uuz*rzj)      &
                            +   qk*(uux*rxk + uuy*ryk + uuz*rzk))     &
                            *  ww2(j,k,irho) - rj*dpj - rk*dpk)/ri
                 enddo
               enddo

#ifndef TAPENADE_REVERSE
               call resetss(nn, ssi, ssj, ssk, ss)
#else
               call resetssBwd(nn, ssi, ssj, ssk, ss)
#endif

           end select BCTreatment

           ! Determine the state in the halo cell. Again loop over
           ! the cell range for this subface.

           do k=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do j=BCData(nn)%icBeg, BCData(nn)%icEnd

               ! Compute the pressure density and velocity in the
               ! halo cell. Note that rface is the grid velocity
               ! component in the direction of norm, i.e. outward
               ! pointing.

               pp1(j,k) = dim(pp2(j,k),pp1(j,k))

               vn = two*(BCData(nn)%rface(j,k) - ww2(j,k,ivx)*BCData(nn)%norm(j,k,1) &
                                    - ww2(j,k,ivy)*BCData(nn)%norm(j,k,2) &
                                    - ww2(j,k,ivz)*BCData(nn)%norm(j,k,3))
               
               ww1(j,k,irho) = ww2(j,k,irho)
               ww1(j,k,ivx)  = ww2(j,k,ivx) + vn*BCData(nn)%norm(j,k,1)
               ww1(j,k,ivy)  = ww2(j,k,ivy) + vn*BCData(nn)%norm(j,k,2)
               ww1(j,k,ivz)  = ww2(j,k,ivz) + vn*BCData(nn)%norm(j,k,3)

               ! Just copy the turbulent variables.

               do l=nt1MG,nt2MG
                 ww1(j,k,l) = ww2(j,k,l)
               enddo

               ! The laminar and eddy viscosity, if present.

               if( viscous )    rlv1(j,k) = rlv2(j,k)
               if( eddyModel ) rev1(j,k) = rev2(j,k)

             enddo
           enddo

           ! deallocation all pointer
#ifndef TAPENADE_REVERSE
           call resetBCPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                              rev1, rev2, 0)
#else
           call resetBCPointersBwd(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
                rev1, rev2, 0)
#endif

           ! Compute the energy for these halo's.

           call computeEtot(icBeg(nn),icEnd(nn), jcBeg(nn),jcEnd(nn), &
                            kcBeg(nn),kcEnd(nn), correctForK)

           ! Extrapolate the state vectors in case a second halo
           ! is needed.

           if( secondHalo ) call extrapolate2ndHalo(nn, correctForK)
  
        endif invWall
       enddo bocos

       end subroutine bcEulerWall
