!
!      ******************************************************************
!      *                                                                *
!      * File:          bcEulerWallAdj.f90                              *
!      * Author:        Edwin van der Weide                             *
!      *                Seongim Choi,C.A.(Sandy) Mader                  *
!      * Starting date: 03-21-2006                                      *
!      * Last modified: 10-27-2008                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine bcEulerWallForcesAdj(secondHalo, wAdj,pAdj,sAdj,      &
     siAdj, sjAdj, skAdj, normAdj,rFaceAdj,mm,&
     iiBeg,iiEnd,jjBeg,jjEnd,i2Beg,i2End,j2Beg,j2End)

  !
  !      ******************************************************************
  !      *                                                                *
  !      * bcEulerWallAdj pplies the inviscid wall boundary condition to  *
  !      * subface nn of the block to which the pointers in blockPointers *
  !      * currently point.                                               *
  !      *                                                                *
  !      ******************************************************************
  !
  use BCTypes
  !       use blockPointers, only : il, jl, kl, BCData, BCFaceID, &
  !                                 si, sj, sk, s, addGridVelocities
  use blockPointers, only : il, jl, kl, BCData, BCFaceID,s, addGridVelocities, nBocos, BCType,ib,jb,kb,ie,je,ke
  use constants
  use flowVarRefState
  use inputDiscretization
  use iteration
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType):: mm
  
  integer(kind=intType) :: wallTreatment

  logical, intent(in) :: secondHalo

  integer(kind=intType), intent(in) :: iiBeg,iiEnd,jjBeg,jjEnd
  integer(kind=intType), intent(in) :: i2Beg,i2End,j2Beg,j2End
  
  !       integer(kind=intType) :: iCell, jCell, kCell
  integer(kind=intType):: icBeg,icEnd,jcBeg,jcEnd
  
  real(kind=realType), dimension(0:ib,0:jb,0:kb,nw), intent(in) :: wAdj
  real(kind=realType), dimension(0:ib,0:jb,0:kb), intent(in) :: pAdj
  
  real(kind=realType), dimension(0:2,iiBeg:iiEnd,jjBeg:jjEnd,3) :: siAdj 
  ! notice the range of y dim is set 1:2 which corresponds to 1/jl
  real(kind=realType), dimension(iiBeg:iiEnd,0:2,jjBeg:jjEnd,3) :: sjAdj
  ! notice the range of z dim is set 1:2 which corresponds to 1/kl
  real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,0:2,3) :: skAdj
  real(kind=realType), dimension(0:ie,0:je,0:ke,3) :: sAdj
  real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3) :: ssAdj
  real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3) :: normAdj
  real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd) :: rFaceAdj

  real(kind=realType), dimension(iibeg:iiend,jjbeg:jjend,nw):: wAdj0, wAdj1
  real(kind=realType), dimension(iibeg:iiend,jjbeg:jjend,nw):: wAdj2, wAdj3
  real(kind=realType), dimension(iibeg:iiend,jjbeg:jjend) :: pAdj0, pAdj1
  real(kind=realType), dimension(iibeg:iiend,jjbeg:jjend) :: pAdj2, pAdj3
  
  real(kind=realType), dimension(0:ib,0:jb,0:kb)::rlvAdj, revAdj
  real(kind=realType), dimension(iibeg:iiend,jjbeg:jjend)::rlvAdj1, rlvAdj2
  real(kind=realType), dimension(iibeg:iiend,jjbeg:jjend)::revAdj1, revAdj2

!!$  real(kind=realType), dimension(-2:2,-2:2,nw) :: wAdj0, wAdj1
!!$  real(kind=realType), dimension(-2:2,-2:2,nw) :: wAdj2, wAdj3
!!$  real(kind=realType), dimension(-2:2,-2:2)    :: pAdj0, pAdj1
!!$  real(kind=realType), dimension(-2:2,-2:2)    :: pAdj2, pAdj3
!!$
!!$  real(kind=realType), dimension(-2:2,-2:2,-2:2)::rlvAdj, revAdj
!!$  real(kind=realType), dimension(-2:2,-2:2)::rlvAdj1, rlvAdj2
!!$  real(kind=realType), dimension(-2:2,-2:2)::revAdj1, revAdj2
!!$
!!$
!!$  real(kind=realType), dimension(-2:2,-2:2,-2:2,3), intent(in) ::siAdj, sjAdj, skAdj
!!$  real(kind=realType), dimension(nBocos,-2:2,-2:2,3), intent(in) :: normAdj
!!$  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw),intent(in) :: wAdj
!!$  real(kind=realType), dimension(-2:2,-2:2,-2:2),intent(in) :: pAdj

  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, k, l, ii, jj, kk
  integer(kind=intType) :: jm1,  jp1,  km1,  kp1
  integer(kind=intType) :: jjm1, jjp1, kkm1, kkp1

  real(kind=realType) :: sixa, siya, siza, sjxa, sjya, sjza
  real(kind=realType) :: skxa, skya, skza, a1, b1
  real(kind=realType) :: rxj, ryj, rzj, rxk, ryk, rzk
  real(kind=realType) :: dpj, dpk, ri, rj, rk, qj, qk, vn
  real(kind=realType) :: ux, uy, uz, ovgm1, gm53, factK
  real(kind=realType) :: rface

  real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3) :: ssi, ssj, ssk
  !  real(kind=realType), dimension(:,:,:), pointer :: ss
  !

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !

  ! Make sure that on the coarser grids the constant pressure
  ! boundary condition is used.

  wallTreatment = wallBcTreatment
  if(currentLevel > groundLevel) wallTreatment = constantPressure

  icbeg = iibeg
  icend = iiend
  jcbeg = jjbeg
  jcend = jjend

  ! Check for Euler wall boundary condition.
  
  invWall: if(BCType(mm) == EulerWall) then

     ! Set the pointers for the unit normal and the normal
     ! velocity to make the code more readable.
     
     !!?norm  => BCData(nn)%norm
     !!?rface => BCData(nn)%rface
     
     !print *,'extracting states'
     !Copy the states and other parameters to subfaces
     call extractBCStatesForcesAdj(mm,wAdj,pAdj,wAdj0, wAdj1, wAdj2,wAdj3,&
          pAdj0,pAdj1, pAdj2,pAdj3,&
          rlvAdj, revAdj,rlvAdj1, rlvAdj2,revAdj1, revAdj2,&
          iibeg,jjbeg,iiend,jjend,secondHalo)
     !print *,'states extracted'
     ! Some initialization
           ssi = zero
           ssj = zero
           ssk = zero

           ! Easier storage of variables involving gamma.

           ovgm1 = one/(gammaInf - one)
           gm53  =  gammaInf - five*third
           factK = -ovgm1*gm53

           ! Determine the boundary condition treatment and compute the
           ! undivided pressure gradient accordingly. This gradient is
           ! temporarily stored in the halo pressure.
           !print *,'selecting case',wallBCTreatment
           BCTreatment: select case (wallBCTreatment)

           case (constantPressure)
              !print *,'constant pressure'
              ! Constant pressure. Set the gradient to zero.
                            
              do j=jcBeg, jcEnd
                 do i=icBeg, icEnd
                    ii = i 
                    jj = j 
                    pAdj1(ii,jj) = zero
                 enddo
              enddo

              !===============================================================

           case (linExtrapolPressure)
              !print *,'linear extrapolation'
              ! Linear extrapolation. Compute the gradient.
              
              do j=jcBeg, jcEnd
                 do i=icBeg, icEnd
                    ii = i 
                    jj = j 
                    
                    pAdj1(ii,jj) = pAdj3(ii,jj) - pAdj2(ii,jj)
                 enddo
              enddo
!!$
!!$
!!$              !===============================================================
!!$
!!$           case (quadExtrapolPressure)
!!$
!!$              ! Quadratic extrapolation. Does not fit within the
!!$              ! current data structures.
!!$              
!!$              !call terminate("bcEulerWallAdj", "Quadratic extrapolation does not fit within the current data structure for the boundary stuff")
!!$              call terminate("bcEulerWallAdj", "Quadratic")
!!$              !print *, "bcEulerWallAdj: quadExtrapolPressure: STOP"
!!$
              !===============================================================

           case (normalMomentum)
              !print *,' normal momentum'
              !call terminate("bcEulerWallAdj", &
              !               "No normal momentum in this version.")
              !print *, "bcEulerWallAdj: STOP"

              ! Pressure gradient is computed using the normal momentum
              ! equation. First set a couple of additional variables for
              ! the normals, depending on the block face. Note that the
              ! construction 1: should not be used in these pointers,
              ! because element 0 is needed. Consequently there will be
              ! an offset of 1 for these normals. This is commented in
              ! the code. For moving faces also the grid velocity of
              ! the 1st cell center from the wall is needed.

              select case (BCFaceID(mm))
              case (iMin)
                 if(secondHalo) then
                    ssi(:,:,:) = siAdj(1,:,:,:)
                    ssj(:,:,:) = sjAdj(2,:,:,:)
                    ssk(:,:,:) = skAdj(2,:,:,:)
                    if( addGridVelocities ) ssAdj(:,:,:) = sAdj(2,:,:,:)
                 else
                    ssi(:,:,:) = siAdj(0,:,:,:)
                    ssj(:,:,:) = sjAdj(1,:,:,:)
                    ssk(:,:,:) = skAdj(1,:,:,:)     
                    if( addGridVelocities ) ssAdj(:,:,:) = sAdj(1,:,:,:)
                 end if
                 

                 !===========================================================

              case (iMax)
                 if(secondHalo)then
                    ssi(:,:,:) = siAdj(ib-2,:,:,:)
                    ssj(:,:,:) = sjAdj(ib-2,:,:,:)
                    ssk(:,:,:) = skAdj(ib-2,:,:,:)
                    if( addGridVelocities ) ssAdj(:,:,:) = sAdj(ib-2,:,:,:)
                 else
                    ssi(:,:,:) = siAdj(ib-1,:,:,:)
                    ssj(:,:,:) = sjAdj(ib-1,:,:,:)
                    ssk(:,:,:) = skAdj(ib-1,:,:,:)
                    if( addGridVelocities ) ssAdj(:,:,:) = sAdj(ib-1,:,:,:)
                 end if

                 

                 !===========================================================

              case (jMin)
                 if(secondHalo) then
                    ssi(:,:,:) = sjAdj(:,1,:,:)
                    ssj(:,:,:) = siAdj(:,2,:,:)
                    ssk(:,:,:) = skAdj(:,2,:,:)
                    if( addGridVelocities ) ssAdj(:,:,:) = sAdj(:,2,:,:)
                 else
                    ssi(:,:,:) = sjAdj(:,0,:,:)
                    ssj(:,:,:) = siAdj(:,1,:,:)
                    ssk(:,:,:) = skAdj(:,1,:,:)
                    if( addGridVelocities ) ssAdj(:,:,:) = sAdj(:,1,:,:)
                 end if

                 !===========================================================

              case (jMax)
                 if(secondHalo) then
                    ssi(:,:,:) = sjAdj(:,jb-2,:,:)
                    ssj(:,:,:) = siAdj(:,jb-2,:,:)
                    ssk(:,:,:) = skAdj(:,jb-2,:,:)
                    if( addGridVelocities ) ssAdj(:,:,:) = sAdj(:,jb-2,:,:)
                 else
                    ssi(:,:,:) = sjAdj(:,jb-1,:,:)
                    ssj(:,:,:) = siAdj(:,jb-1,:,:)
                    ssk(:,:,:) = skAdj(:,jb-1,:,:)
                    if( addGridVelocities ) ssAdj(:,:,:) = sAdj(:,jb-1,:,:)
                 end if

                 

                 !===========================================================

              case (kMin)
                 if(secondHalo) then
                    ssi(:,:,:) = skAdj(:,:,1,:)
                    ssj(:,:,:) = siAdj(:,:,2,:)
                    ssk(:,:,:) = sjAdj(:,:,2,:)
                    if( addGridVelocities ) ssAdj(:,:,:) = sAdj(:,:,2,:)
                 else
                    ssi(:,:,:) = skAdj(:,:,0,:)
                    ssj(:,:,:) = siAdj(:,:,1,:)
                    ssk(:,:,:) = sjAdj(:,:,1,:)
                    if( addGridVelocities ) ssAdj(:,:,:) = sAdj(:,:,1,:)
                 end if

                

                 !===========================================================

              case (kMax)
                 if(secondHalo) then
                    ssi(:,:,:) = skAdj(:,:,kb-2,:)
                    ssj(:,:,:) = siAdj(:,:,kb-2,:)
                    ssk(:,:,:) = sjAdj(:,:,kb-2,:)
                    if( addGridVelocities ) ssAdj(:,:,:) = sAdj(:,:,kb-2,:)

                 else
                    ssi(:,:,:) = skAdj(:,:,kb-1,:)
                    ssj(:,:,:) = siAdj(:,:,kb-1,:)
                    ssk(:,:,:) = sjAdj(:,:,kb-1,:)
                    if( addGridVelocities ) ssAdj(:,:,:) = sAdj(:,:,kb-1,:)

                 end if

                 
              end select

              ! Loop over the faces of the generic subface.
              ! Note that now the running indices are j and k. This is
              ! done, because the generic i-direction is assumed to
              ! be the normal direction.

              do k=jcBeg, jcEnd

                 ! Store the indices k+1, k-1 a bit easier and make
                 ! sure that they do not exceed the range of the arrays.

                 km1 = k-1; km1 = max(BCData(mm)%jcBeg,km1)
                 kp1 = k+1; kp1 = min(BCData(mm)%jcEnd,kp1)

                 ! Compute the scaling factor for the central difference
                 ! in the k-direction.

                 b1 = one/max(1,(kp1-km1))

                 ! Compute the offset indices.

                 kk   = k   
                 kkm1 = km1 
                 kkp1 = kp1 

                 ! The generic j-direction.

                 do j=icBeg, icEnd

                    ! The indices j+1 and j-1. Make sure that they
                    ! do not exceed the range of the arrays.

                    jm1 = j-1; jm1 = max(BCData(mm)%icBeg,jm1)
                    jp1 = j+1; jp1 = min(BCData(mm)%icEnd,jp1)

                    ! Compute the scaling factor for the central
                    ! difference in the j-direction.

                    a1 = one/max(1,(jp1-jm1))

                    jj   = j   
                    jjm1 = jm1 
                    jjp1 = jp1 

                    ! Compute (twice) the average normal in the generic i,
                    ! j and k-direction. Note that in j and k-direction
                    ! the average in the original indices should be taken
                    ! using j-1 and j (and k-1 and k). However due to the
                    ! usage of pointers ssj and ssk there is an offset in
                    ! the indices of 1 and therefore now the correct
                    ! average is obtained with the indices j and j+1
                    ! (k and k+1).

                    sixa = two*ssi(jj,kk,1)
                    siya = two*ssi(jj,kk,2)
                    siza = two*ssi(jj,kk,3)

                    sjxa = ssj(jj-1,kk,1) + ssj(jj,kk,1)! it was ssj(j,k,1) + ssj(j+1,k,1)
                    sjya = ssj(jj-1,kk,2) + ssj(jj,kk,2)! it was ssj(j,k,2) + ssj(j+1,k,2)
                    sjza = ssj(jj-1,kk,3) + ssj(jj,kk,3)! it was ssj(j,k,3) + ssj(j+1,k,3)

                    skxa = ssk(jj,kk-1,1) + ssk(jj,kk,1) ! it was ssk(j,k,1) + ssk(j,k+1,1)
                    skya = ssk(jj,kk-1,2) + ssk(jj,kk,2) ! it was ssk(j,k,2) + ssk(j,k+1,2)
                    skza = ssk(jj,kk-1,3) + ssk(jj,kk,3) ! it was ssk(j,k,3) + ssk(j,k+1,3)

                    ! Compute the difference of the normal vector and
                    ! pressure in j and k-direction. As the indices are
                    ! restricted to the 1st halo-layer, the computation
                    ! of the internal halo values is not consistent;
                    ! however this is not really a problem, because these
                    ! values are overwritten in the communication pattern.

                    rxj = a1*(normAdj(jjp1,kk,1) - normAdj(jjm1,kk,1))
                    ryj = a1*(normAdj(jjp1,kk,2) - normAdj(jjm1,kk,2))
                    rzj = a1*(normAdj(jjp1,kk,3) - normAdj(jjm1,kk,3))
                    !print *, "jjp1,jjm1, kk =", jjp1,jjm1, kk
                    dpj = a1*(pAdj2(jjp1,kk) - pAdj2(jjm1,kk))

                    rxk = b1*(normAdj(jj,kkp1,1) - normAdj(jj,kkm1,1))
                    ryk = b1*(normAdj(jj,kkp1,2) - normAdj(jj,kkm1,2))
                    rzk = b1*(normAdj(jj,kkp1,3) - normAdj(jj,kkm1,3))
                    !print *, "jj, kkp1, kkm1 =", jj, kkp1, kkm1
                    dpk = b1*(pAdj2(jj,kkp1) - pAdj2(jj,kkm1))

                    ! Compute the dot product between the unit vector
                    ! and the normal vectors in i, j and k-direction.

                    ri = normAdj(jj,kk,1)*sixa + normAdj(jj,kk,2)*siya &
                         + normAdj(jj,kk,3)*siza
                    rj = normAdj(jj,kk,1)*sjxa + normAdj(jj,kk,2)*sjya &
                         + normAdj(jj,kk,3)*sjza
                    rk = normAdj(jj,kk,1)*skxa + normAdj(jj,kk,2)*skya &
                         + normAdj(jj,kk,3)*skza

                    ! Store the velocity components in ux, uy and uz and
                    ! subtract the mesh velocity if the face is moving.

                    ux = wAdj2(jj,kk,ivx)
                    uy = wAdj2(jj,kk,ivy)
                    uz = wAdj2(jj,kk,ivz)

                    if( addGridVelocities ) then
                       ux = ux - ssAdj(jj,kk,1)
                       uy = uy - ssAdj(jj,kk,2)
                       uz = uz - ssAdj(jj,kk,3)
                    endif

                    ! Compute the velocity components in j and
                    ! k-direction.

                    qj = ux*sjxa + uy*sjya + uz*sjza
                    qk = ux*skxa + uy*skya + uz*skza

                    ! Compute the pressure gradient, which is stored
                    ! in pAdj1. I'm not entirely sure whether this
                    ! formulation is correct for moving meshes. It could
                    ! be that an additional term is needed there.

                    pAdj1(jj,kk) = ((qj*(ux*rxj + uy*ryj + uz*rzj)      &
                         +   qk*(ux*rxk + uy*ryk + uz*rzk))     &
                         *  wAdj2(jj,kk,irho) - rj*dpj - rk*dpk)/ri


                 end do
              end do
              
           end select BCTreatment

           ! Determine the state in the halo cell. Again loop over
           ! the cell range for this subface.
           !print *,' determining halo cells'

           do j=jcBeg, jcEnd
              do i=icBeg, icEnd
                 ii = i 
                 jj = j 
                 
                 rface = rFaceAdj(ii,jj) !BCData(mm)%rface(i,j)

                 ! Compute the pressure density and velocity in the
                 ! halo cell. Note that rface is the grid velocity
                 ! component in the direction of norm, i.e. outward
                 ! pointing.

                 pAdj1(ii,jj) = max(zero, pAdj2(ii,jj)-pAdj1(ii,jj) )
!!$
!!$!                 vn = two*(BCData(nn)%rface(i,j)              &
!!$!                      - wAdj2(ii,jj,ivx)*normAdj(nn,ii,jj,1) &
!!$!                      - wAdj2(ii,jj,ivy)*normAdj(nn,ii,jj,2) &
!!$ !                     - wAdj2(ii,jj,ivz)*normAdj(nn,ii,jj,3))
!!$
                 vn = two*(rface - wAdj2(ii,jj,ivx)*normAdj(ii,jj,1) &
                      - wAdj2(ii,jj,ivy)*normAdj(ii,jj,2) &
                      - wAdj2(ii,jj,ivz)*normAdj(ii,jj,3))

                 wAdj1(ii,jj,irho) = wAdj2(ii,jj,irho)
                 wAdj1(ii,jj,ivx)  = wAdj2(ii,jj,ivx) + vn*normAdj(ii,jj,1)
                 wAdj1(ii,jj,ivy)  = wAdj2(ii,jj,ivy) + vn*normAdj(ii,jj,2)
                 wAdj1(ii,jj,ivz)  = wAdj2(ii,jj,ivz) + vn*normAdj(ii,jj,3)

                 ! Just copy the turbulent variables.

                 do l=nt1MG,nt2MG
                    wAdj1(ii,jj,l) = wAdj2(ii,jj,l)
                 enddo


                 !
                 !        Input the viscous effects - rlv1(), and rev1()
                 !

                 ! Compute the total energy.

                 wAdj1(ii,jj,irhoE) = ovgm1*pAdj1(ii,jj)     &
                      + half*wAdj1(ii,jj,irho) &
                      *     (wAdj1(ii,jj,ivx)**2 &
                      +      wAdj1(ii,jj,ivy)**2 &
                      +      wAdj1(ii,jj,ivz)**2)

                 if( kPresent )                            &
                      wAdj1(ii,jj,irhoE) = wAdj1(ii,jj,irhoE) &
                      - factK*wAdj1(ii,jj,irho) &
                      *       wAdj1(ii,jj,itu1)
              enddo
           enddo

 
           ! Extrapolate the state vectors in case a second halo
           ! is needed.

           
           if( secondHalo )   &
                call extrapolate2ndHaloForcesAdj(mm,iiBeg, iiEnd, jjBeg, jjEnd,  &
                wAdj0, wAdj1, &
                wAdj2, pAdj0, pAdj1, pAdj2)
           
           call replaceBCStatesForcesAdj(mm,  wAdj0,wAdj1, wAdj2, wAdj3,&
                pAdj0,pAdj1, pAdj2, pAdj3,rlvAdj1, rlvAdj2,revAdj1, revAdj2,&
                wAdj,pAdj,rlvAdj,revAdj,&
                iibeg,jjbeg,iiend,jjend,secondHalo)
            
        endif invWall

      end subroutine bcEulerWallForcesAdj
