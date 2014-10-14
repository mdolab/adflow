
!
!      ******************************************************************
!      *                                                                *
!      * File:          setss.f90                                       *
!      * Author:        Eirikur Jonsson                                 *
!      * Starting date: 10-14-2014                                      *
!      * Last modified: 10-14-2014                                      *
!      *                                                                *
!      ******************************************************************

       subroutine setss(nn, ssi, ssj, ssk, ss)
       
       use BCTypes
       use blockPointers
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
       real(kind=realType), dimension(:,:,:), pointer :: ssi, ssj, ssk
       real(kind=realType), dimension(:,:,:), pointer :: ss

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the face id on which the subface is located and set
       ! the pointers accordinly.

		   select case (BCFaceID(nn))
			 case (iMin)
			   ssi => si(1,:,:,:)
			   ssj => sj(2,:,:,:)
			   ssk => sk(2,:,:,:)

			   if( addGridVelocities ) ss => s(2,:,:,:)

			 !=======================================================

			 case (iMax)
			   ssi => si(il,:,:,:)
			   ssj => sj(il,:,:,:)
			   ssk => sk(il,:,:,:)

			   if( addGridVelocities ) ss => s(il,:,:,:)

			 !=======================================================

			 case (jMin)
			   ssi => sj(:,1,:,:)
			   ssj => si(:,2,:,:)
			   ssk => sk(:,2,:,:)

			   if( addGridVelocities ) ss => s(:,2,:,:)

			 !=======================================================

			 case (jMax)
			   ssi => sj(:,jl,:,:)
			   ssj => si(:,jl,:,:)
			   ssk => sk(:,jl,:,:)

			   if( addGridVelocities ) ss => s(:,jl,:,:)

			 !=======================================================

			 case (kMin)
			   ssi => sk(:,:,1,:)
			   ssj => si(:,:,2,:)
			   ssk => sj(:,:,2,:)

			   if( addGridVelocities ) ss => s(:,:,2,:)

			 !=======================================================

			 case (kMax)
			   ssi => sk(:,:,kl,:)
			   ssj => si(:,:,kl,:)
			   ssk => sj(:,:,kl,:)

			   if( addGridVelocities ) ss => s(:,:,kl,:)

		   end select
       end subroutine setss


