subroutine bcTurbWall(nn)
  !
  !       bcTurbWall applies the implicit treatment of the viscous       
  !       wall boundary condition for the turbulence model used to the   
  !       given subface nn.                                              
  !       It is assumed that the pointers in blockPointers are           
  !       already set to the correct block.                              
  !
  use blockPointers
  use flowVarRefState
  use inputPhysics
  use constants
  use paramTurb
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, ii, jj, iiMax, jjMax

  real(kind=realType) :: tmpd, tmpe, tmpf, nu

  real(kind=realType), dimension(:,:,:,:), pointer :: bmt
  real(kind=realType), dimension(:,:,:),   pointer :: bvt, ww2
  real(kind=realType), dimension(:,:),     pointer :: rlv2, dd2Wall



  ! Determine the turbulence model used and loop over the faces
  ! of the subface and set the values of bmt and bvt for an
  ! implicit treatment.

  select case (turbModel)

  case (spalartAllmaras, spalartAllmarasEdwards)

     ! Spalart-allmaras type of model. Value at the wall is zero,
     ! so simply negate the internal value.
     select case (BCFaceID(nn))
     case (iMin)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              bmti1(i,j,itu1,itu1) = one
           enddo
        enddo
     case (iMax)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              bmti2(i,j,itu1,itu1) = one
           enddo
        enddo
     case (jMin)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              bmtj1(i,j,itu1,itu1) = one
           enddo
        enddo
     case (jMax)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              bmtj2(i,j,itu1,itu1) = one
           enddo
        enddo

     case (kMin)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              bmtk1(i,j,itu1,itu1) = one
           enddo
        enddo

     case (kMax)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              bmtk2(i,j,itu1,itu1) = one
           enddo
        enddo
     end select

     !        ================================================================
  case (komegaWilcox, komegaModified, menterSST)

     ! K-omega type of models. K is zero on the wall and thus the
     ! halo value is the negative of the first internal cell.
     ! For omega the situation is a bit more complicated.
     ! Theoretically omega is infinity, but it is set to a large
     ! value, see menter's paper. The halo value is constructed
     ! such that the wall value is correct. Make sure that i and j
     ! are limited to physical dimensions of the face for the wall
     ! distance. Due to the usage of the dd2Wall pointer and the
     ! fact that the original d2Wall array starts at 2, there is
     ! an offset of -1 present in dd2Wall.

     select case (BCFaceID(nn))
     case (iMin)
        iiMax = jl; jjMax = kl

        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           jj = max(2,min(j,jjMax))

           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              ii = max(2,min(i,iiMax))

              nu   = rlv(2,i,j)/w(2,i,j,irho)
              tmpd = one/(rkwBeta1*(d2Wall(2,ii,jj)**2))
              
              bmti1(i,j,itu1,itu1) = one
              bmti1(i,j,itu2,itu2) = one

              bvti1(i,j,itu2) = two*60.0_realType*nu*tmpd
           enddo
        enddo 

     case (iMax)
        iiMax = jl; jjMax = kl

        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           jj = max(2,min(j,jjMax))

           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              ii = max(2,min(i,iiMax))

              nu   = rlv(jl,i,j)/w(il,i,j,irho)
              tmpd = one/(rkwBeta1*(d2Wall(il,ii,jj)**2))
              
              bmti2(i,j,itu1,itu1) = one
              bmti2(i,j,itu2,itu2) = one

              bvti2(i,j,itu2) = two*60.0_realType*nu*tmpd
           enddo
        enddo 

     case (jMin)
        iiMax = il; jjMax = kl

        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           jj = max(2,min(j,jjMax))

           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              ii = max(2,min(i,iiMax))

              nu   = rlv(i,2,j)/w(i,2,j,irho)
              tmpd = one/(rkwBeta1*(d2Wall(ii,2,jj)**2))
              
              bmtj1(i,j,itu1,itu1) = one
              bmtj1(i,j,itu2,itu2) = one

              bvtj1(i,j,itu2) = two*60.0_realType*nu*tmpd
           enddo
        enddo 

     case (jMax)
        iiMax = il; jjMax = kl

        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           jj = max(2,min(j,jjMax))

           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              ii = max(2,min(i,iiMax))

              nu   = rlv(i,jl,j)/w(i,jl,j,irho)
              tmpd = one/(rkwBeta1*(d2Wall(ii,jl,jj)**2))
              
              bmtj2(i,j,itu1,itu1) = one
              bmtj2(i,j,itu2,itu2) = one

              bvtj2(i,j,itu2) = two*60.0_realType*nu*tmpd
           enddo
        enddo 

     case (kMin)
        iiMax = il; jjMax = jl

        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           jj = max(2,min(j,jjMax))

           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              ii = max(2,min(i,iiMax))

              nu   = rlv(i,j,2)/w(i,j,2,irho)
              tmpd = one/(rkwBeta1*(d2Wall(ii,jj,2)**2))
              
              bmtk1(i,j,itu1,itu1) = one
              bmtk1(i,j,itu2,itu2) = one

              bvtk1(i,j,itu2) = two*60.0_realType*nu*tmpd
           enddo
        enddo 

     case (kMax)
        iiMax = il; jjMax = jl

        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           jj = max(2,min(j,jjMax))

           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              ii = max(2,min(i,iiMax))

              nu   = rlv(i,j,kl)/w(i,j,kl,irho)
              tmpd = one/(rkwBeta1*(d2Wall(ii,jj,kl)**2))
              
              bmtk2(i,j,itu1,itu1) = one
              bmtk2(i,j,itu2,itu2) = one

              bvtk2(i,j,itu2) = two*60.0_realType*nu*tmpd
           enddo
        enddo 
     end select

     !        ================================================================

  case (ktau)

     ! K-tau model. Both k and tau are zero at the wall, so the
     ! negative value of the internal cell is taken for the halo.
     select case (BCFaceID(nn))
     case (iMin)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              bmti1(i,j,itu1,itu1) = one
              bmti1(i,j,itu2,itu2) = one
           enddo
        enddo
     case (iMax)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              bmti2(i,j,itu1,itu1) = one
              bmti2(i,j,itu2,itu2) = one
           enddo
        enddo
     case (jMin)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              bmtj1(i,j,itu1,itu1) = one
              bmtj1(i,j,itu2,itu2) = one
           enddo
        enddo
     case (jMax)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              bmtj2(i,j,itu1,itu1) = one
              bmtj2(i,j,itu2,itu2) = one
           enddo
        enddo

     case (kMin)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              bmtk1(i,j,itu1,itu1) = one
              bmtk1(i,j,itu2,itu2) = one
           enddo
        enddo

     case (kMax)
        do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
           do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              bmtk2(i,j,itu1,itu1) = one
              bmtk2(i,j,itu2,itu2) = one
           enddo
        enddo
     end select

     !        ================================================================
#ifndef USE_TAPENADE
  case (v2f)
     
     ! Set some variables depending on the block face on which the
     ! subface is located. Needed for a general treatment.

     select case (BCFaceID(nn))
     case (iMin)
        iiMax = jl; jjMax = kl
        bmt => bmti1; bvt => bvti1; ww2 => w(2 ,1:,1:,:)
        rlv2 => rlv(2, 1:,1:); dd2Wall => d2Wall(2, :,:)

     case (iMax)
        iiMax = jl; jjMax = kl
        bmt => bmti2; bvt => bvti2; ww2 => w(il,1:,1:,:)
        rlv2 => rlv(il,1:,1:); dd2Wall => d2Wall(il,:,:)

     case (jMin)
        iiMax = il; jjMax = kl
        bmt => bmtj1; bvt => bvtj1; ww2 => w(1:,2 ,1:,:)
        rlv2 => rlv(1:,2 ,1:); dd2Wall => d2Wall(:,2 ,:)

     case (jMax)
        iiMax = il; jjMax = kl
        bmt => bmtj2; bvt => bvtj2; ww2 => w(1:,jl,1:,:)
        rlv2 => rlv(1:,jl,1:); dd2Wall => d2Wall(:,jl,:)

     case (kMin)
        iiMax = il; jjMax = jl
        bmt => bmtk1; bvt => bvtk1; ww2 => w(1:,1:,2 ,:)
        rlv2 => rlv(1:,1:,2 ); dd2Wall => d2Wall(:,:,2 )

     case (kMax)
        iiMax = il; jjMax = jl
        bmt => bmtk2; bvt => bvtk2; ww2 => w(1:,1:,kl,:)
        rlv2 => rlv(1:,1:,kl); dd2Wall => d2Wall(:,:,kl)
     end select

     ! V2f turbulence model. Same story for the wall distance as
     ! for k-omega. For this model there is a coupling between the
     ! equations via the boundary conditions.

     do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
        jj = max(2,min(j,jjMax))

        do i=BCData(nn)%icBeg, BCData(nn)%icEnd
           ii = max(2,min(i,iiMax))

           nu   = rlv2(i,j)/ww2(i,j,irho)
           tmpd = one/(dd2Wall(ii-1,jj-1)**2)
           tmpe = two*nu*tmpd
           tmpf =-20.0_realType*(nu*tmpd)**2 &
                / abs(tmpe*ww2(i,j,itu1))
           if(rvfN == 6) tmpf = zero

           bmt(i,j,itu1,itu1) = one
           bmt(i,j,itu2,itu2) = one
           bmt(i,j,itu3,itu3) = one
           bmt(i,j,itu4,itu4) = one

           bmt(i,j,itu2,itu1) = -two*tmpe
           bmt(i,j,itu4,itu3) = -two*tmpf
        enddo
     enddo
#endif
  end select
end subroutine bcTurbWall
