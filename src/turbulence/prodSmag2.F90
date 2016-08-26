!
!       File:          prodSmag2.f90                                   
!       Author:        Georgi Kalitzin, Edwin van der Weide            
!       Starting date: 08-01-2003                                      
!       Last modified: 04-12-2005                                      
!
subroutine prodSmag2
  !
  !       prodSmag2 computes the term:                                   
  !              2*sij*sij - 2/3 div(u)**2 with  sij=0.5*(duidxj+dujdxi) 
  !       which is used for the turbulence equations.                    
  !       It is assumed that the pointer prod, stored in turbMod, is     
  !       already set to the correct entry.                              
  !
  use constants
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
  integer(kind=intType) :: i, j, k, ii
  real(kind=realType)   :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
  real(kind=realType)   :: div2, fact, sxx, syy, szz, sxy, sxz, syz

  ! Loop over the cell centers of the given block. It may be more
  ! efficient to loop over the faces and to scatter the gradient,
  ! but in that case the gradients for u, v and w must be stored.
  ! In the current approach no extra memory is needed.

#ifdef TAPENADE_FAST
  !$AD II-LOOP
  do ii=0,nx*ny*nz-1
     i = mod(ii, nx) + 2
     j = mod(ii/nx, ny) + 2
     k = ii/(nx*ny) + 2
#else
     do k=2, kl
        do j=2, jl
           do i=2, il
#endif   

              ! Compute the gradient of u in the cell center. Use is made
              ! of the fact that the surrounding normals sum up to zero,
              ! such that the cell i,j,k does not give a contribution.
              ! The gradient is scaled by the factor 2*vol.

              uux = w(i+1,j,k,ivx)*si(i,j,k,1) - w(i-1,j,k,ivx)*si(i-1,j,k,1) &
                   + w(i,j+1,k,ivx)*sj(i,j,k,1) - w(i,j-1,k,ivx)*sj(i,j-1,k,1) &
                   + w(i,j,k+1,ivx)*sk(i,j,k,1) - w(i,j,k-1,ivx)*sk(i,j,k-1,1)
              uuy = w(i+1,j,k,ivx)*si(i,j,k,2) - w(i-1,j,k,ivx)*si(i-1,j,k,2) &
                   + w(i,j+1,k,ivx)*sj(i,j,k,2) - w(i,j-1,k,ivx)*sj(i,j-1,k,2) &
                   + w(i,j,k+1,ivx)*sk(i,j,k,2) - w(i,j,k-1,ivx)*sk(i,j,k-1,2)
              uuz = w(i+1,j,k,ivx)*si(i,j,k,3) - w(i-1,j,k,ivx)*si(i-1,j,k,3) &
                   + w(i,j+1,k,ivx)*sj(i,j,k,3) - w(i,j-1,k,ivx)*sj(i,j-1,k,3) &
                   + w(i,j,k+1,ivx)*sk(i,j,k,3) - w(i,j,k-1,ivx)*sk(i,j,k-1,3)

              ! Idem for the gradient of v.

              vvx = w(i+1,j,k,ivy)*si(i,j,k,1) - w(i-1,j,k,ivy)*si(i-1,j,k,1) &
                   + w(i,j+1,k,ivy)*sj(i,j,k,1) - w(i,j-1,k,ivy)*sj(i,j-1,k,1) &
                   + w(i,j,k+1,ivy)*sk(i,j,k,1) - w(i,j,k-1,ivy)*sk(i,j,k-1,1)
              vvy = w(i+1,j,k,ivy)*si(i,j,k,2) - w(i-1,j,k,ivy)*si(i-1,j,k,2) &
                   + w(i,j+1,k,ivy)*sj(i,j,k,2) - w(i,j-1,k,ivy)*sj(i,j-1,k,2) &
                   + w(i,j,k+1,ivy)*sk(i,j,k,2) - w(i,j,k-1,ivy)*sk(i,j,k-1,2)
              vvz = w(i+1,j,k,ivy)*si(i,j,k,3) - w(i-1,j,k,ivy)*si(i-1,j,k,3) &
                   + w(i,j+1,k,ivy)*sj(i,j,k,3) - w(i,j-1,k,ivy)*sj(i,j-1,k,3) &
                   + w(i,j,k+1,ivy)*sk(i,j,k,3) - w(i,j,k-1,ivy)*sk(i,j,k-1,3)

              ! And for the gradient of w.

              wwx = w(i+1,j,k,ivz)*si(i,j,k,1) - w(i-1,j,k,ivz)*si(i-1,j,k,1) &
                   + w(i,j+1,k,ivz)*sj(i,j,k,1) - w(i,j-1,k,ivz)*sj(i,j-1,k,1) &
                   + w(i,j,k+1,ivz)*sk(i,j,k,1) - w(i,j,k-1,ivz)*sk(i,j,k-1,1)
              wwy = w(i+1,j,k,ivz)*si(i,j,k,2) - w(i-1,j,k,ivz)*si(i-1,j,k,2) &
                   + w(i,j+1,k,ivz)*sj(i,j,k,2) - w(i,j-1,k,ivz)*sj(i,j-1,k,2) &
                   + w(i,j,k+1,ivz)*sk(i,j,k,2) - w(i,j,k-1,ivz)*sk(i,j,k-1,2)
              wwz = w(i+1,j,k,ivz)*si(i,j,k,3) - w(i-1,j,k,ivz)*si(i-1,j,k,3) &
                   + w(i,j+1,k,ivz)*sj(i,j,k,3) - w(i,j-1,k,ivz)*sj(i,j-1,k,3) &
                   + w(i,j,k+1,ivz)*sk(i,j,k,3) - w(i,j,k-1,ivz)*sk(i,j,k-1,3)

              ! Compute the components of the stress tensor.
              ! The combination of the current scaling of the velocity
              ! gradients (2*vol) and the definition of the stress tensor,
              ! leads to the factor 1/(4*vol).

              fact = fourth/vol(i,j,k)

              sxx = two*fact*uux
              syy = two*fact*vvy
              szz = two*fact*wwz

              sxy = fact*(uuy + vvx)
              sxz = fact*(uuz + wwx)
              syz = fact*(vvz + wwy)

              ! Compute 2/3 * divergence of velocity squared

              div2 = f23*(sxx+syy+szz)**2

              ! Store the square of strain as the production term.

              scratch(i,j,k, iprod) = two*(two*(sxy**2 + sxz**2 + syz**2) &
                   +           sxx**2 + syy**2 + szz**2) - div2
#ifdef TAPENADE_FAST
           end do
#else
        enddo
     enddo
  enddo
#endif  
end subroutine prodSmag2
