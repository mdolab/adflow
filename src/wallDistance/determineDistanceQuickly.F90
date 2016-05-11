subroutine updateWallDistancesQuickly(nn, level, sps)

  ! This is the actual update routine that uses xSurf. It is done on
  ! block-level-sps basis.  This is the used to update the wall
  ! distance. Most importantly, this routine is included in the
  ! reverse mode AD routines, but NOT the forward mode. Since it is
  ! done on a per-block basis, it is assumed that the required block
  ! pointers are already set. 

  use blockPointers
  use wallDistanceData
  implicit none

  ! Subroutine arguments
  integer(kind=intType) :: nn, level, sps

  ! Local Variables
  integer(kind=intType) :: i, j, k, ii, ind(4)
  real(kind=realType) :: xp(3), xc(3), u, v

#ifdef TAPENADE_FAST
  !$AD II-LOOP
  do ii=0,nx*ny*nz-1
     i = mod(ii, nx) + 2
     j = mod(ii/nx, ny) + 2
     k = ii/(nx*ny) + 2
#else
     do k=2,kl
        do j=2,jl
           do i=2,il
#endif

              if (flowDoms(nn, level, sps)%surfNodeIndices(1, i, j, k) == 0) then 
                 ! This node is too far away and has no
                 ! association. Set the distance to a large constant. 
                 d2wall(i, j, k) = large
                 cycle
              end if

              ! Extract elemID and u-v position for the association of
              ! this cell:

              ind = flowDoms(nn,level,sps)%surfNodeIndices(:, i, j, k)
              u   = flowDoms(nn,level,sps)%uv(1,i,j,k)
              v    = flowDoms(nn,level,sps)%uv(2,i,j,k)

              ! Now we have the 4 corners, use bi-linear shape
              ! functions o to get target: (CCW ordering remember!)
            
              xp(:) = &
                   (one-u)*(one-v)*xSurf(3*(ind(1)-1)+1:3*ind(1)) + &
                   (    u)*(one-v)*xSurf(3*(ind(2)-1)+1:3*ind(2)) + &
                   (    u)*(    v)*xSurf(3*(ind(3)-1)+1:3*ind(3)) + &
                   (one-u)*(    v)*xSurf(3*(ind(4)-1)+1:3*ind(4)) 
  
              ! Get the cell center
              xc(1) = eighth*(x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
                   +         x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
                   +         x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                   +         x(i-1,j,  k,  1) + x(i,j,  k,  1))
              
              xc(2) = eighth*(x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
                   +         x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
                   +         x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                   +         x(i-1,j,  k,  2) + x(i,j,  k,  2))
              
              xc(3) = eighth*(x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
                   +         x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
                   +         x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                   +         x(i-1,j,  k,  3) + x(i,j,  k,  3))
              

              ! Now we have the two points...just take the norm of the
              ! distance between them
              
              d2wall(i,j,k) = sqrt(&
                   (xc(1)-xp(1))**2 + (xc(2)-xp(2))**2 + (xc(3)-xp(3))**2)
#ifdef TAPENADE_FAST
           end do
#else
        enddo
     enddo
  enddo
#endif
  
end subroutine updateWallDistancesQuickly
