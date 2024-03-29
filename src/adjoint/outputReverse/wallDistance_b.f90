!        generated by tapenade     (inria, ecuador team)
!  tapenade 3.16 (develop) - 22 aug 2023 15:51
!
module walldistance_b
  use constants, only : inttype, realtype
  use walldistancedata
  implicit none
  save 

contains
!  differentiation of updatewalldistancesquickly in reverse (adjoint) mode (with options noisize i4 dr8 r8):
!   gradient     of useful results: *x *d2wall *xsurf
!   with respect to varying inputs: *x *d2wall *xsurf
!   rw status of diff variables: *x:incr *d2wall:in-out *xsurf:incr
!   plus diff mem management of: x:in d2wall:in xsurf:in
  subroutine updatewalldistancesquickly_b(nn, level, sps)
! this is the actual update routine that uses xsurf. it is done on
! block-level-sps basis.  this is the used to update the wall
! distance. most importantly, this routine is included in the
! reverse mode ad routines, but not the forward mode. since it is
! done on a per-block basis, it is assumed that the required block
! pointers are already set.
    use constants
    use blockpointers, only : nx, ny, nz, il, jl, kl, x, xd, flowdoms,&
&   flowdomsd, d2wall, d2walld
    implicit none
! subroutine arguments
    integer(kind=inttype) :: nn, level, sps
! local variables
    integer(kind=inttype) :: i, j, k, ii, ind(4)
    real(kind=realtype) :: xp(3), xc(3), u, v
    real(kind=realtype) :: xpd(3), xcd(3)
    intrinsic mod
    intrinsic sqrt
    real(kind=realtype) :: tempd
    real(kind=realtype) :: tempd0
    real(kind=realtype) :: tempd1
    real(kind=realtype) :: tempd2
    xcd = 0.0_8
!$bwd-of ii-loop 
    do ii=0,nx*ny*nz-1
      i = mod(ii, nx) + 2
      j = mod(ii/nx, ny) + 2
      k = ii/(nx*ny) + 2
      if (flowdoms(nn, level, sps)%surfnodeindices(1, i, j, k) .eq. 0) &
&     then
        d2walld(i, j, k) = 0.0_8
      else
! extract elemid and u-v position for the association of
! this cell:
        ind = flowdoms(nn, level, sps)%surfnodeindices(:, i, j, k)
        u = flowdoms(nn, level, sps)%uv(1, i, j, k)
        v = flowdoms(nn, level, sps)%uv(2, i, j, k)
! now we have the 4 corners, use bi-linear shape
! functions o to get target: (ccw ordering remember!)
        xp(:) = (one-u)*(one-v)*xsurf(3*(ind(1)-1)+1:3*ind(1)) + u*(one-&
&         v)*xsurf(3*(ind(2)-1)+1:3*ind(2)) + u*v*xsurf(3*(ind(3)-1)+1:3&
&         *ind(3)) + (one-u)*v*xsurf(3*(ind(4)-1)+1:3*ind(4))
! get the cell center
        xc(1) = eighth*(x(i-1, j-1, k-1, 1)+x(i, j-1, k-1, 1)+x(i-1, j, &
&         k-1, 1)+x(i, j, k-1, 1)+x(i-1, j-1, k, 1)+x(i, j-1, k, 1)+x(i-&
&         1, j, k, 1)+x(i, j, k, 1))
        xc(2) = eighth*(x(i-1, j-1, k-1, 2)+x(i, j-1, k-1, 2)+x(i-1, j, &
&         k-1, 2)+x(i, j, k-1, 2)+x(i-1, j-1, k, 2)+x(i, j-1, k, 2)+x(i-&
&         1, j, k, 2)+x(i, j, k, 2))
        xc(3) = eighth*(x(i-1, j-1, k-1, 3)+x(i, j-1, k-1, 3)+x(i-1, j, &
&         k-1, 3)+x(i, j, k-1, 3)+x(i-1, j-1, k, 3)+x(i, j-1, k, 3)+x(i-&
&         1, j, k, 3)+x(i, j, k, 3))
! now we have the two points...just take the norm of the
! distance between them
        xpd = 0.0_8
        if ((xc(1)-xp(1))**2 + (xc(2)-xp(2))**2 + (xc(3)-xp(3))**2 .eq. &
&           0.0_8) then
          tempd = 0.0_8
        else
          tempd = d2walld(i, j, k)/(2.0*sqrt((xc(1)-xp(1))**2+(xc(2)-xp(&
&           2))**2+(xc(3)-xp(3))**2))
        end if
        d2walld(i, j, k) = 0.0_8
        tempd0 = 2*(xc(1)-xp(1))*tempd
        tempd1 = 2*(xc(2)-xp(2))*tempd
        tempd2 = 2*(xc(3)-xp(3))*tempd
        xcd(3) = xcd(3) + tempd2
        xpd(3) = xpd(3) - tempd2
        xcd(2) = xcd(2) + tempd1
        xpd(2) = xpd(2) - tempd1
        xcd(1) = xcd(1) + tempd0
        xpd(1) = xpd(1) - tempd0
        tempd = eighth*xcd(3)
        xcd(3) = 0.0_8
        xd(i-1, j-1, k-1, 3) = xd(i-1, j-1, k-1, 3) + tempd
        xd(i, j-1, k-1, 3) = xd(i, j-1, k-1, 3) + tempd
        xd(i-1, j, k-1, 3) = xd(i-1, j, k-1, 3) + tempd
        xd(i, j, k-1, 3) = xd(i, j, k-1, 3) + tempd
        xd(i-1, j-1, k, 3) = xd(i-1, j-1, k, 3) + tempd
        xd(i, j-1, k, 3) = xd(i, j-1, k, 3) + tempd
        xd(i-1, j, k, 3) = xd(i-1, j, k, 3) + tempd
        xd(i, j, k, 3) = xd(i, j, k, 3) + tempd
        tempd = eighth*xcd(2)
        xcd(2) = 0.0_8
        xd(i-1, j-1, k-1, 2) = xd(i-1, j-1, k-1, 2) + tempd
        xd(i, j-1, k-1, 2) = xd(i, j-1, k-1, 2) + tempd
        xd(i-1, j, k-1, 2) = xd(i-1, j, k-1, 2) + tempd
        xd(i, j, k-1, 2) = xd(i, j, k-1, 2) + tempd
        xd(i-1, j-1, k, 2) = xd(i-1, j-1, k, 2) + tempd
        xd(i, j-1, k, 2) = xd(i, j-1, k, 2) + tempd
        xd(i-1, j, k, 2) = xd(i-1, j, k, 2) + tempd
        xd(i, j, k, 2) = xd(i, j, k, 2) + tempd
        tempd = eighth*xcd(1)
        xcd(1) = 0.0_8
        xd(i-1, j-1, k-1, 1) = xd(i-1, j-1, k-1, 1) + tempd
        xd(i, j-1, k-1, 1) = xd(i, j-1, k-1, 1) + tempd
        xd(i-1, j, k-1, 1) = xd(i-1, j, k-1, 1) + tempd
        xd(i, j, k-1, 1) = xd(i, j, k-1, 1) + tempd
        xd(i-1, j-1, k, 1) = xd(i-1, j-1, k, 1) + tempd
        xd(i, j-1, k, 1) = xd(i, j-1, k, 1) + tempd
        xd(i-1, j, k, 1) = xd(i-1, j, k, 1) + tempd
        xd(i, j, k, 1) = xd(i, j, k, 1) + tempd
        xsurfd(3*(ind(1)-1)+1:3*ind(1)) = xsurfd(3*(ind(1)-1)+1:3*ind(1)&
&         ) + (one-u)*(one-v)*xpd
        xsurfd(3*(ind(2)-1)+1:3*ind(2)) = xsurfd(3*(ind(2)-1)+1:3*ind(2)&
&         ) + u*(one-v)*xpd
        xsurfd(3*(ind(3)-1)+1:3*ind(3)) = xsurfd(3*(ind(3)-1)+1:3*ind(3)&
&         ) + u*v*xpd
        xsurfd(3*(ind(4)-1)+1:3*ind(4)) = xsurfd(3*(ind(4)-1)+1:3*ind(4)&
&         ) + (one-u)*v*xpd
      end if
    end do
  end subroutine updatewalldistancesquickly_b

  subroutine updatewalldistancesquickly(nn, level, sps)
! this is the actual update routine that uses xsurf. it is done on
! block-level-sps basis.  this is the used to update the wall
! distance. most importantly, this routine is included in the
! reverse mode ad routines, but not the forward mode. since it is
! done on a per-block basis, it is assumed that the required block
! pointers are already set.
    use constants
    use blockpointers, only : nx, ny, nz, il, jl, kl, x, flowdoms, &
&   d2wall
    implicit none
! subroutine arguments
    integer(kind=inttype) :: nn, level, sps
! local variables
    integer(kind=inttype) :: i, j, k, ii, ind(4)
    real(kind=realtype) :: xp(3), xc(3), u, v
    intrinsic mod
    intrinsic sqrt
!$ad ii-loop
    do ii=0,nx*ny*nz-1
      i = mod(ii, nx) + 2
      j = mod(ii/nx, ny) + 2
      k = ii/(nx*ny) + 2
      if (flowdoms(nn, level, sps)%surfnodeindices(1, i, j, k) .eq. 0) &
&     then
! this node is too far away and has no
! association. set the distance to a large constant.
        d2wall(i, j, k) = large
      else
! extract elemid and u-v position for the association of
! this cell:
        ind = flowdoms(nn, level, sps)%surfnodeindices(:, i, j, k)
        u = flowdoms(nn, level, sps)%uv(1, i, j, k)
        v = flowdoms(nn, level, sps)%uv(2, i, j, k)
! now we have the 4 corners, use bi-linear shape
! functions o to get target: (ccw ordering remember!)
        xp(:) = (one-u)*(one-v)*xsurf(3*(ind(1)-1)+1:3*ind(1)) + u*(one-&
&         v)*xsurf(3*(ind(2)-1)+1:3*ind(2)) + u*v*xsurf(3*(ind(3)-1)+1:3&
&         *ind(3)) + (one-u)*v*xsurf(3*(ind(4)-1)+1:3*ind(4))
! get the cell center
        xc(1) = eighth*(x(i-1, j-1, k-1, 1)+x(i, j-1, k-1, 1)+x(i-1, j, &
&         k-1, 1)+x(i, j, k-1, 1)+x(i-1, j-1, k, 1)+x(i, j-1, k, 1)+x(i-&
&         1, j, k, 1)+x(i, j, k, 1))
        xc(2) = eighth*(x(i-1, j-1, k-1, 2)+x(i, j-1, k-1, 2)+x(i-1, j, &
&         k-1, 2)+x(i, j, k-1, 2)+x(i-1, j-1, k, 2)+x(i, j-1, k, 2)+x(i-&
&         1, j, k, 2)+x(i, j, k, 2))
        xc(3) = eighth*(x(i-1, j-1, k-1, 3)+x(i, j-1, k-1, 3)+x(i-1, j, &
&         k-1, 3)+x(i, j, k-1, 3)+x(i-1, j-1, k, 3)+x(i, j-1, k, 3)+x(i-&
&         1, j, k, 3)+x(i, j, k, 3))
! now we have the two points...just take the norm of the
! distance between them
        d2wall(i, j, k) = sqrt((xc(1)-xp(1))**2 + (xc(2)-xp(2))**2 + (xc&
&         (3)-xp(3))**2)
      end if
    end do
  end subroutine updatewalldistancesquickly
! ----------------------------------------------------------------------
!                                                                      |
!                    no tapenade routine below this line               |
!                                                                      |
! ----------------------------------------------------------------------

end module walldistance_b

