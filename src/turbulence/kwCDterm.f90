!
!      ******************************************************************
!      *                                                                *
!      * File:          kwCDterm.f90                                    *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 07-09-2003                                      *
!      * Last modified: 04-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine kwCDterm
!
!      ******************************************************************
!      *                                                                *
!      * kwCDterm computes the cross-diffusion term in the omega-eqn    *
!      * for the SST version as well as the modified k-omega turbulence *
!      * model. It is assumed that the pointers in blockPointers and    *
!      * turbMod are already set.                                       *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use turbMod
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k
       real(kind=realType)   :: kx, ky, kz, wx, wy, wz
       real(kind=realType)   :: lnwip1, lnwim1, lnwjp1, lnwjm1
       real(kind=realType)   :: lnwkp1, lnwkm1
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the cell centers of the given block. It may be more
       ! efficient to loop over the faces and to scatter the gradient,
       ! but in that case the gradients for k and omega must be stored.
       ! In the current approach no extra memory is needed.

       do k=2,kl
         do j=2,jl
           do i=2,il

             ! Compute the gradient of k in the cell center. Use is made
             ! of the fact that the surrounding normals sum up to zero,
             ! such that the cell i,j,k does not give a contribution.
             ! The gradient is scaled by a factor 1/2vol.

             kx = w(i+1,j,k,itu1)*si(i,j,k,1) - w(i-1,j,k,itu1)*si(i-1,j,k,1) &
                + w(i,j+1,k,itu1)*sj(i,j,k,1) - w(i,j-1,k,itu1)*sj(i,j-1,k,1) &
                + w(i,j,k+1,itu1)*sk(i,j,k,1) - w(i,j,k-1,itu1)*sk(i,j,k-1,1)
             ky = w(i+1,j,k,itu1)*si(i,j,k,2) - w(i-1,j,k,itu1)*si(i-1,j,k,2) &
                + w(i,j+1,k,itu1)*sj(i,j,k,2) - w(i,j-1,k,itu1)*sj(i,j-1,k,2) &
                + w(i,j,k+1,itu1)*sk(i,j,k,2) - w(i,j,k-1,itu1)*sk(i,j,k-1,2)
             kz = w(i+1,j,k,itu1)*si(i,j,k,3) - w(i-1,j,k,itu1)*si(i-1,j,k,3) &
                + w(i,j+1,k,itu1)*sj(i,j,k,3) - w(i,j-1,k,itu1)*sj(i,j-1,k,3) &
                + w(i,j,k+1,itu1)*sk(i,j,k,3) - w(i,j,k-1,itu1)*sk(i,j,k-1,3)

             ! Compute the logarithm of omega in the points that
             ! contribute to the gradient in this cell.

             lnwip1 = log(abs(w(i+1,j,k,itu2)))
             lnwim1 = log(abs(w(i-1,j,k,itu2)))
             lnwjp1 = log(abs(w(i,j+1,k,itu2)))
             lnwjm1 = log(abs(w(i,j-1,k,itu2)))
             lnwkp1 = log(abs(w(i,j,k+1,itu2)))
             lnwkm1 = log(abs(w(i,j,k-1,itu2)))

             ! Compute the scaled gradient of ln omega.

             wx = lnwip1*si(i,j,k,1) - lnwim1*si(i-1,j,k,1) &
                + lnwjp1*sj(i,j,k,1) - lnwjm1*sj(i,j-1,k,1) &
                + lnwkp1*sk(i,j,k,1) - lnwkm1*sk(i,j,k-1,1)
             wy = lnwip1*si(i,j,k,2) - lnwim1*si(i-1,j,k,2) &
                + lnwjp1*sj(i,j,k,2) - lnwjm1*sj(i,j-1,k,2) &
                + lnwkp1*sk(i,j,k,2) - lnwkm1*sk(i,j,k-1,2)
             wz = lnwip1*si(i,j,k,3) - lnwim1*si(i-1,j,k,3) &
                + lnwjp1*sj(i,j,k,3) - lnwjm1*sj(i,j-1,k,3) &
                + lnwkp1*sk(i,j,k,3) - lnwkm1*sk(i,j,k-1,3)

             ! Compute the dot product grad k grad ln omega.
             ! Multiply it by the correct scaling factor and store it.

             kwCD(i,j,k) = fourth*(kx*wx + ky*wy + kz*wz)/(vol(i,j,k)**2)

           enddo
         enddo
       enddo

       end subroutine kwCDterm
