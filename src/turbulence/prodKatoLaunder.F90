!
!      ******************************************************************
!      *                                                                *
!      * File:          prodKatoLaunder.f90                             *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 08-01-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine prodKatoLaunder
  !
  !      ******************************************************************
  !      *                                                                *
  !      * prodKatoLaunder computes the turbulent production term using   *
  !      * the Kato-Launder formulation.                                  *
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
  integer(kind=intType) :: i, j, k

  real(kind=realType) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
  real(kind=realType) :: qxx, qyy, qzz, qxy, qxz, qyz, sijsij
  real(kind=realType) :: oxy, oxz, oyz, oijoij
  real(kind=realType) :: fact, omegax, omegay, omegaz
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Determine the non-dimensional wheel speed of this block.
  ! The vorticity term, which appears in Kato-Launder is of course
  ! not frame invariant. To approximate frame invariance the wheel
  ! speed should be substracted from oxy, oxz and oyz, which results
  ! in the vorticity in the rotating frame. However some people
  ! claim that the absolute vorticity should be used to obtain the
  ! best results. In that omega should be set to zero.

  omegax = timeRef*sections(sectionID)%rotRate(1)
  omegay = timeRef*sections(sectionID)%rotRate(2)
  omegaz = timeRef*sections(sectionID)%rotRate(3)

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
              ! The gradient is scaled by a factor 2*vol.

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

              ! Compute the strain and vorticity terms. The multiplication
              ! is present to obtain the correct gradients. Note that
              ! the wheel speed is substracted from the vorticity terms.

              fact = half/vol(i,j,k)

              qxx = fact*uux
              qyy = fact*vvy
              qzz = fact*wwz

              qxy = fact*half*(uuy + vvx)
              qxz = fact*half*(uuz + wwx)
              qyz = fact*half*(vvz + wwy)

              oxy = fact*half*(vvx - uuy) - omegaz
              oxz = fact*half*(uuz - wwx) - omegay
              oyz = fact*half*(wwy - vvz) - omegax

              ! Compute the summation of the strain and vorticity tensors.

              sijsij = two*(qxy**2 + qxz**2 + qyz**2) &
                   +      qxx**2 + qyy**2 + qzz**2
              oijoij = two*(oxy**2 + oxz**2 + oyz**2)

              ! Compute the production term.

              dw(i,j,k, iprod) = two*sqrt(sijsij*oijoij)
#ifdef TAPENADE_FAST
           end do
#else
        enddo
     enddo
  enddo
#endif  

end subroutine prodKatoLaunder
