
!
!      ******************************************************************
!      *                                                                *
!      * File:          resetSSBwd.f90                                  *
!      * Author:        Eirikur Jonsson                                 *
!      * Starting date: 10-14-2014                                      *
!      * Last modified: 10-14-2014                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine resetssBwd(nn, ssi, ssj, ssk, ss)
       
       use BCTypes
       use blockPointers
!       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
       real(kind=realType), dimension(imaxDim,jmaxDim,3) :: ssi, ssj, ssk
       real(kind=realType), dimension(imaxDim,jmaxDim,3) :: ss

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
			   si(1,:,:,:) = ssi 
			   sj(2,:,:,:) = ssj
			   sk(2,:,:,:) = ssk

			   if( addGridVelocities ) s(2,:,:,:) = ss

			 !=======================================================

			 case (iMax)
			   si(il,:,:,:) = ssi
			   sj(il,:,:,:) = ssj
			   sk(il,:,:,:) = ssk

			   if( addGridVelocities ) s(il,:,:,:) = ss

			 !=======================================================

			 case (jMin)
			   sj(:,1,:,:) = ssi
			   si(:,2,:,:) = ssj
			   sk(:,2,:,:) = ssk

			   if( addGridVelocities ) s(:,2,:,:) = ss

			 !=======================================================

			 case (jMax)
			   sj(:,jl,:,:) = ssi
			   si(:,jl,:,:) = ssj
			   sk(:,jl,:,:) = ssk

			   if( addGridVelocities ) s(:,jl,:,:) = ss

			 !=======================================================

			 case (kMin)
			   sk(:,:,1,:) = ssi
			   si(:,:,2,:) = ssj
			   sj(:,:,2,:) = ssk

			   if( addGridVelocities ) s(:,:,2,:) = ss

			 !=======================================================

			 case (kMax)
			   sk(:,:,kl,:) = ssi
			   si(:,:,kl,:) = ssj
			   sj(:,:,kl,:) = ssk

			   if( addGridVelocities ) s(:,:,kl,:) = ss

		   end select

       end subroutine resetssBwd


