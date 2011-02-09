!
!      ******************************************************************
!      *                                                                *
!      * File:          ktCDterm.f90                                    *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 07-09-2003                                      *
!      * Last modified: 04-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine ktCDterm
!
!      ******************************************************************
!      *                                                                *
!      * ktCdterm computes the cross-diffusion term in the tau-eqn      *
!      * for the k-tau turbulence model for the given block.            *
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
       real(kind=realType)   :: kx, ky, kz, tx, ty, tz, cd
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the cell centers of the given block. It may be more
       ! efficient to loop over the faces and to scatter the gradient,
       ! but in that case the gradients for k and tau must be stored.
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

             ! Idem for tau.

             tx = w(i+1,j,k,itu2)*si(i,j,k,1) - w(i-1,j,k,itu2)*si(i-1,j,k,1) &
                + w(i,j+1,k,itu2)*sj(i,j,k,1) - w(i,j-1,k,itu2)*sj(i,j-1,k,1) &
                + w(i,j,k+1,itu2)*sk(i,j,k,1) - w(i,j,k-1,itu2)*sk(i,j,k-1,1)
             ty = w(i+1,j,k,itu2)*si(i,j,k,2) - w(i-1,j,k,itu2)*si(i-1,j,k,2) &
                + w(i,j+1,k,itu2)*sj(i,j,k,2) - w(i,j-1,k,itu2)*sj(i,j-1,k,2) &
                + w(i,j,k+1,itu2)*sk(i,j,k,2) - w(i,j,k-1,itu2)*sk(i,j,k-1,2)
             tz = w(i+1,j,k,itu2)*si(i,j,k,3) - w(i-1,j,k,itu2)*si(i-1,j,k,3) &
                + w(i,j+1,k,itu2)*sj(i,j,k,3) - w(i,j-1,k,itu2)*sj(i,j-1,k,3) &
                + w(i,j,k+1,itu2)*sk(i,j,k,3) - w(i,j,k-1,itu2)*sk(i,j,k-1,3)

             ! Compute the dot product grad k grad tau. Multiply it by
             ! the correct scaling factor.

             cd = fourth*(kx*tx + ky*ty + kz*tz)/(vol(i,j,k)**2)

             ! Cd must only be taken into account when it is negative.
             ! Take care of that here.

             ktCD(i,j,k) = min(zero, cd)

           enddo
         enddo
       enddo

       end subroutine ktCDterm
