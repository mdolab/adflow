!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of turb2ndhalo in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *rev *w
!   with respect to varying inputs: *rev *w
!   plus diff mem management of: rev:in w:in bcdata:in
!
!      ******************************************************************
!      *                                                                *
!      * file:          turb2ndhalo.f90                                 *
!      * author:        edwin van der weide                             *
!      * starting date: 06-16-2003                                      *
!      * last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine turb2ndhalo_b(nn)
!
!      ******************************************************************
!      *                                                                *
!      * turb2ndhalo sets the turbulent variables in the second halo    *
!      * cell for the given subface. simple constant extrapolation is   *
!      * used to avoid problems.                                        *
!      *                                                                *
!      ******************************************************************
!
  use blockpointers
  use bctypes
  use flowvarrefstate
  implicit none
!
!      subroutine arguments.
!
  integer(kind=inttype), intent(in) :: nn
!
!      local variables.
!
  integer(kind=inttype) :: i, j, l
  real(kind=realtype) :: tmp
  real(kind=realtype) :: tmp0
  real(kind=realtype) :: tmp1
  real(kind=realtype) :: tmp2
  real(kind=realtype) :: tmp3
  real(kind=realtype) :: tmp4
  integer :: branch
  real(kind=realtype) :: tmpd
  real(kind=realtype) :: tmpd4
  real(kind=realtype) :: tmpd3
  real(kind=realtype) :: tmpd2
  real(kind=realtype) :: tmpd1
  real(kind=realtype) :: tmpd0
!
!      ******************************************************************
!      *                                                                *
!      * begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
! determine the face on which this subface is located and set
! some pointers accordingly.
! loop over the turbulent variables and set the second halo
! value. if this is an eddy model, also set the eddy viscosity.
  select case  (bcfaceid(nn)) 
  case (imin) 
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        if (eddymodel) then
          call pushcontrol1b(1)
        else
          call pushcontrol1b(0)
        end if
      end do
    end do
    do j=bcdata(nn)%jcend,bcdata(nn)%jcbeg,-1
      do i=bcdata(nn)%icend,bcdata(nn)%icbeg,-1
        call popcontrol1b(branch)
        if (branch .ne. 0) then
          revd(1, i, j) = revd(1, i, j) + revd(0, i, j)
          revd(0, i, j) = 0.0_8
        end if
        do l=nt2,nt1,-1
          wd(1, i, j, l) = wd(1, i, j, l) + wd(0, i, j, l)
          wd(0, i, j, l) = 0.0_8
        end do
      end do
    end do
  case (imax) 
!===============================================================
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        if (eddymodel) then
          call pushcontrol1b(1)
        else
          call pushcontrol1b(0)
        end if
      end do
    end do
    do j=bcdata(nn)%jcend,bcdata(nn)%jcbeg,-1
      do i=bcdata(nn)%icend,bcdata(nn)%icbeg,-1
        call popcontrol1b(branch)
        if (branch .ne. 0) then
          tmpd0 = revd(ib, i, j)
          revd(ib, i, j) = 0.0_8
          revd(ie, i, j) = revd(ie, i, j) + tmpd0
        end if
        do l=nt2,nt1,-1
          tmpd = wd(ib, i, j, l)
          wd(ib, i, j, l) = 0.0_8
          wd(ie, i, j, l) = wd(ie, i, j, l) + tmpd
        end do
      end do
    end do
  case (jmin) 
!===============================================================
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        if (eddymodel) then
          call pushcontrol1b(1)
        else
          call pushcontrol1b(0)
        end if
      end do
    end do
    do j=bcdata(nn)%jcend,bcdata(nn)%jcbeg,-1
      do i=bcdata(nn)%icend,bcdata(nn)%icbeg,-1
        call popcontrol1b(branch)
        if (branch .ne. 0) then
          revd(i, 1, j) = revd(i, 1, j) + revd(i, 0, j)
          revd(i, 0, j) = 0.0_8
        end if
        do l=nt2,nt1,-1
          wd(i, 1, j, l) = wd(i, 1, j, l) + wd(i, 0, j, l)
          wd(i, 0, j, l) = 0.0_8
        end do
      end do
    end do
  case (jmax) 
!===============================================================
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        if (eddymodel) then
          call pushcontrol1b(1)
        else
          call pushcontrol1b(0)
        end if
      end do
    end do
    do j=bcdata(nn)%jcend,bcdata(nn)%jcbeg,-1
      do i=bcdata(nn)%icend,bcdata(nn)%icbeg,-1
        call popcontrol1b(branch)
        if (branch .ne. 0) then
          tmpd2 = revd(i, jb, j)
          revd(i, jb, j) = 0.0_8
          revd(i, je, j) = revd(i, je, j) + tmpd2
        end if
        do l=nt2,nt1,-1
          tmpd1 = wd(i, jb, j, l)
          wd(i, jb, j, l) = 0.0_8
          wd(i, je, j, l) = wd(i, je, j, l) + tmpd1
        end do
      end do
    end do
  case (kmin) 
!===============================================================
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        if (eddymodel) then
          call pushcontrol1b(1)
        else
          call pushcontrol1b(0)
        end if
      end do
    end do
    do j=bcdata(nn)%jcend,bcdata(nn)%jcbeg,-1
      do i=bcdata(nn)%icend,bcdata(nn)%icbeg,-1
        call popcontrol1b(branch)
        if (branch .ne. 0) then
          revd(i, j, 1) = revd(i, j, 1) + revd(i, j, 0)
          revd(i, j, 0) = 0.0_8
        end if
        do l=nt2,nt1,-1
          wd(i, j, 1, l) = wd(i, j, 1, l) + wd(i, j, 0, l)
          wd(i, j, 0, l) = 0.0_8
        end do
      end do
    end do
  case (kmax) 
!===============================================================
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        if (eddymodel) then
          call pushcontrol1b(1)
        else
          call pushcontrol1b(0)
        end if
      end do
    end do
    do j=bcdata(nn)%jcend,bcdata(nn)%jcbeg,-1
      do i=bcdata(nn)%icend,bcdata(nn)%icbeg,-1
        call popcontrol1b(branch)
        if (branch .ne. 0) then
          tmpd4 = revd(i, j, kb)
          revd(i, j, kb) = 0.0_8
          revd(i, j, ke) = revd(i, j, ke) + tmpd4
        end if
        do l=nt2,nt1,-1
          tmpd3 = wd(i, j, kb, l)
          wd(i, j, kb, l) = 0.0_8
          wd(i, j, ke, l) = wd(i, j, ke, l) + tmpd3
        end do
      end do
    end do
  end select
end subroutine turb2ndhalo_b
