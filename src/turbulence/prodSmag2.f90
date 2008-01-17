!
!      ******************************************************************
!      *                                                                *
!      * File:          prodSmag2.f90                                   *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 08-01-2003                                      *
!      * Last modified: 04-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine prodSmag2
!
!      ******************************************************************
!      *                                                                *
!      * prodSmag2 computes the term:                                   *
!      *        2*sij*sij - 2/3 div(u)**2 with  sij=0.5*(duidxj+dujdxi) *
!      * which is used for the turbulence equations.                    *
!      * It is assumed that the pointer prod, stored in turbMod, is     *
!      * already set to the correct entry.                              *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use turbMod
       implicit none
!
!      Local parameter
!
       real(kind=realType), parameter :: f23 = two*third
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k
       real(kind=realType)   :: ux, uy, uz, vx, vy, vz, wx, wy, wz
       real(kind=realType)   :: div2, fact, sxx, syy, szz, sxy, sxz, syz
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the cell centers of the given block. It may be more
       ! efficient to loop over the faces and to scatter the gradient,
       ! but in that case the gradients for u, v and w must be stored.
       ! In the current approach no extra memory is needed.

       do k=2,kl
         do j=2,jl
           do i=2,il

             ! Compute the gradient of u in the cell center. Use is made
             ! of the fact that the surrounding normals sum up to zero,
             ! such that the cell i,j,k does not give a contribution.
             ! The gradient is scaled by the factor 2*vol.

             ux = w(i+1,j,k,ivx)*si(i,j,k,1) - w(i-1,j,k,ivx)*si(i-1,j,k,1) &
                + w(i,j+1,k,ivx)*sj(i,j,k,1) - w(i,j-1,k,ivx)*sj(i,j-1,k,1) &
                + w(i,j,k+1,ivx)*sk(i,j,k,1) - w(i,j,k-1,ivx)*sk(i,j,k-1,1)
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
             vy = w(i+1,j,k,ivy)*si(i,j,k,2) - w(i-1,j,k,ivy)*si(i-1,j,k,2) &
                + w(i,j+1,k,ivy)*sj(i,j,k,2) - w(i,j-1,k,ivy)*sj(i,j-1,k,2) &
                + w(i,j,k+1,ivy)*sk(i,j,k,2) - w(i,j,k-1,ivy)*sk(i,j,k-1,2)
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
             wz = w(i+1,j,k,ivz)*si(i,j,k,3) - w(i-1,j,k,ivz)*si(i-1,j,k,3) &
                + w(i,j+1,k,ivz)*sj(i,j,k,3) - w(i,j-1,k,ivz)*sj(i,j-1,k,3) &
                + w(i,j,k+1,ivz)*sk(i,j,k,3) - w(i,j,k-1,ivz)*sk(i,j,k-1,3)

             ! Compute the components of the stress tensor.
             ! The combination of the current scaling of the velocity
             ! gradients (2*vol) and the definition of the stress tensor,
             ! leads to the factor 1/(4*vol).

             fact = fourth/vol(i,j,k)

             sxx = two*fact*ux
             syy = two*fact*vy
             szz = two*fact*wz

             sxy = fact*(uy + vx)
             sxz = fact*(uz + wx)
             syz = fact*(vz + wy)

             ! Compute 2/3 * divergence of velocity squared

             div2 = f23*(sxx+syy+szz)**2

             ! Store the square of strain as the production term.

             prod(i,j,k) = two*(two*(sxy**2 + sxz**2 + syz**2) &
                         +           sxx**2 + syy**2 + szz**2) - div2

           enddo
         enddo
       enddo

       end subroutine prodSmag2
