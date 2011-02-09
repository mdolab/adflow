!
!      ******************************************************************
!      *                                                                *
!      * File:          bcTurbWall.F90                                  *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 06-26-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcTurbWall(nn)
!
!      ******************************************************************
!      *                                                                *
!      * bcTurbWall applies the implicit treatment of the viscous       *
!      * wall boundary condition for the turbulence model used to the   *
!      * given subface nn.                                              *
!      * It is assumed that the pointers in blockPointers are           *
!      * already set to the correct block.                              *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
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
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
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

       ! Determine the turbulence model used and loop over the faces
       ! of the subface and set the values of bmt and bvt for an
       ! implicit treatment.

       select case (turbModel)

         case (spalartAllmaras, spalartAllmarasEdwards)

           ! Spalart-allmaras type of model. Value at the wall is zero,
           ! so simply negate the internal value.

           do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
               bmt(i,j,itu1,itu1) = one
             enddo
           enddo

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

           do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             jj = max(2_intType,min(j,jjMax))

             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
               ii = max(2_intType,min(i,iiMax))

               nu   = rlv2(i,j)/ww2(i,j,irho)
               tmpd = one/(rkwBeta1*(dd2Wall(ii-1,jj-1)**2))

               bmt(i,j,itu1,itu1) = one
               bmt(i,j,itu2,itu2) = one

               bvt(i,j,itu2) = two*60.0_realType*nu*tmpd
             enddo
           enddo

!        ================================================================

         case (ktau)

           ! K-tau model. Both k and tau are zero at the wall, so the
           ! negative value of the internal cell is taken for the halo.

           do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
               bmt(i,j,itu1,itu1) = one
               bmt(i,j,itu2,itu2) = one
             enddo
           enddo

!        ================================================================

         case (v2f)

           ! V2f turbulence model. Same story for the wall distance as
           ! for k-omega. For this model there is a coupling between the
           ! equations via the boundary conditions.

           do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
             jj = max(2_intType,min(j,jjMax))

             do i=BCData(nn)%icBeg, BCData(nn)%icEnd
               ii = max(2_intType,min(i,iiMax))

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

         case default
           call terminate("bcTurbWall", &
                          "Turbulence model not implemented yet")

       end select

       end subroutine bcTurbWall
