!
!      ******************************************************************
!      *                                                                *
!      * File:          prodWmag2.f90                                   *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 06-23-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine prodWmag2
!
!      ******************************************************************
!      *                                                                *
!      * prodWmag2 computes the term:                                   *
!      *    2*oij*oij  with oij=0.5*(duidxj - dujdxi).                  *
!      * This is equal to the magnitude squared of the vorticity.       *
!      * It is assumed that the pointer vort, stored in turbMod, is     *
!      * already set to the correct entry.                              *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use flowVarRefState
       use section
       use turbMod
       implicit none
!
!      Local variables.
!
       integer :: i, j, k

       real(kind=realType) :: uy, uz, vx, vz, wx, wy
       real(kind=realType) :: fact, vortx, vorty, vortz
       real(kind=realType) :: omegax, omegay, omegaz
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the non-dimensional wheel speed of this block.

       omegax = timeRef*sections(sectionID)%rotRate(1)
       omegay = timeRef*sections(sectionID)%rotRate(2)
       omegaz = timeRef*sections(sectionID)%rotRate(3)

       ! Loop over the cell centers of the given block. It may be more
       ! efficient to loop over the faces and to scatter the gradient,
       ! but in that case the gradients for u, v and w must be stored.
       ! In the current approach no extra memory is needed.

       do k=2,kl
         do j=2,jl
           do i=2,il

             ! Compute the necessary derivatives of u in the cell center.
             ! Use is made of the fact that the surrounding normals sum up
             ! to zero, such that the cell i,j,k does not give a
             ! contribution. The gradient is scaled by a factor 2*vol.

             uy = w(i+1,j,k,ivx)*si(i,j,k,2) - w(i-1,j,k,ivx)*si(i-1,j,k,2) &
                + w(i,j+1,k,ivx)*sj(i,j,k,2) - w(i,j-1,k,ivx)*sj(i,j-1,k,2) &
                + w(i,j,k+1,ivx)*sk(i,j,k,2) - w(i,j,k-1,ivx)*sk(i,j,k-1,2)
             uz = w(i+1,j,k,ivx)*si(i,j,k,3) - w(i-1,j,k,ivx)*si(i-1,j,k,3) &
                + w(i,j+1,k,ivx)*sj(i,j,k,3) - w(i,j-1,k,ivx)*sj(i,j-1,k,3) &
                + w(i,j,k+1,ivx)*sk(i,j,k,3) - w(i,j,k-1,ivx)*sk(i,j,k-1,3)

             ! Idem for the gradient of v.

             vx = w(i+1,j,k,ivy)*si(i,j,k,1) - w(i-1,j,k,ivy)*si(i-1,j,k,1) &
                + w(i,j+1,k,ivy)*sj(i,j,k,1) - w(i,j-1,k,ivy)*sj(i,j-1,k,1) &
                + w(i,j,k+1,ivy)*sk(i,j,k,1) - w(i,j,k-1,ivy)*sk(i,j,k-1,1)
             vz = w(i+1,j,k,ivy)*si(i,j,k,3) - w(i-1,j,k,ivy)*si(i-1,j,k,3) &
                + w(i,j+1,k,ivy)*sj(i,j,k,3) - w(i,j-1,k,ivy)*sj(i,j-1,k,3) &
                + w(i,j,k+1,ivy)*sk(i,j,k,3) - w(i,j,k-1,ivy)*sk(i,j,k-1,3)

             ! And for the gradient of w.

             wx = w(i+1,j,k,ivz)*si(i,j,k,1) - w(i-1,j,k,ivz)*si(i-1,j,k,1) &
                + w(i,j+1,k,ivz)*sj(i,j,k,1) - w(i,j-1,k,ivz)*sj(i,j-1,k,1) &
                + w(i,j,k+1,ivz)*sk(i,j,k,1) - w(i,j,k-1,ivz)*sk(i,j,k-1,1)
             wy = w(i+1,j,k,ivz)*si(i,j,k,2) - w(i-1,j,k,ivz)*si(i-1,j,k,2) &
                + w(i,j+1,k,ivz)*sj(i,j,k,2) - w(i,j-1,k,ivz)*sj(i,j-1,k,2) &
                + w(i,j,k+1,ivz)*sk(i,j,k,2) - w(i,j,k-1,ivz)*sk(i,j,k-1,2)

             ! Compute the three components of the vorticity vector.
             ! Substract the part coming from the rotating frame.

             fact = half/vol(i,j,k)

             vortx = fact*(wy - vz) - two*omegax
             vorty = fact*(uz - wx) - two*omegay
             vortz = fact*(vx - uy) - two*omegaz

             ! Compute the magnitude squared of the vorticity.

             vort(i,j,k) = vortx**2 + vorty**2 + vortz**2

           enddo
         enddo
       enddo

       end subroutine prodWmag2
