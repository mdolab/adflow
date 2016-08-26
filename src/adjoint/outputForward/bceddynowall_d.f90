!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of bceddynowall in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: *rev
!   with respect to varying inputs: *rev
!   plus diff mem management of: rev:in bcdata:in
subroutine bceddynowall_d(nn)
!
!       bceddynowall sets the eddy viscosity in the halo cells of      
!       subface nn of the block given in blockpointers. the boundary   
!       condition on the subface can be anything but a viscous wall.   
!       a homogeneous neumann condition is applied, which means that   
!       the eddy viscosity is simply copied from the interior cell.    
!
  use constants
  use blockpointers
  implicit none
!
!      subroutine arguments.
!
  integer(kind=inttype), intent(in) :: nn
!
!      local variables.
!
  integer(kind=inttype) :: i, j
! determine the face id on which the subface and copy
  select case  (bcfaceid(nn)) 
  case (imin) 
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        revd(1, i, j) = revd(2, i, j)
        rev(1, i, j) = rev(2, i, j)
      end do
    end do
  case (imax) 
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        revd(ie, i, j) = revd(il, i, j)
        rev(ie, i, j) = rev(il, i, j)
      end do
    end do
  case (jmin) 
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        revd(i, 1, j) = revd(i, 2, j)
        rev(i, 1, j) = rev(i, 2, j)
      end do
    end do
  case (jmax) 
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        revd(i, je, j) = revd(i, jl, j)
        rev(i, je, j) = rev(i, jl, j)
      end do
    end do
  case (kmin) 
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        revd(i, j, 1) = revd(i, j, 2)
        rev(i, j, 1) = rev(i, j, 2)
      end do
    end do
  case (kmax) 
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        revd(i, j, ke) = revd(i, j, kl)
        rev(i, j, ke) = rev(i, j, kl)
      end do
    end do
  end select
end subroutine bceddynowall_d
