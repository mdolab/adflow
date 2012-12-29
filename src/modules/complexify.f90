!******************************************************************************
!	Written for 'complexify.py 1.3'
!	J.R.R.A.Martins 1999
!       21-Apr-00  Fixed tan, sinh, cosh
!                  sign now returns complex
!                  added log10 and nint
!                  changed ==, /= and >= -- see comments below
!       20-May-00  added cosd, sind, and epsilon
!       11-Jul-00  took away cosd, sind (they are reserved, but not
!                  intrinsic functions in F90)
!       21-Jul-00  converted all trig functions to the value/derivative
!                  formulas -- not general complex number formulas
!       15-Aug-00  Fixed bug in atan2 formula and added the rest of the
!                  _ci and _ic cominations to the relational operators.
!                  P. Sturdza
!                  
!******************************************************************************
!
! Assume all code is compiled with double precision (-r8 compiler flag)
!

!TODO:
!     more typ combinations: cc, cr, rc, ic ?
!     check all fcns
!

module complexify


  implicit none
  
! ABS
  interface abs
     module procedure abs_c
  end interface

! COSD
!  interface cosd
!     module procedure cosd_c
!  end interface

! ACOS
  interface acos
     module procedure acos_c
  end interface

! SIND
!  interface sind
!     module procedure sind_c
!  end interface

! ASIN
  interface asin
     module procedure asin_c
  end interface

! ATAN
  interface atan
     module procedure atan_c
  end interface

! ATAN2
  interface atan2
     module procedure atan2_cc
  end interface

! COSH
  interface cosh
     module procedure cosh_c
  end interface

! MAX (limited to 2-4 complex args, 2 mixed args)
  interface max
     module procedure max_cc
     module procedure max_cr
     module procedure max_rc
     module procedure max_ccc     ! added because of DFLUX.f
     module procedure max_cccc     ! added because of DFLUX.f
  end interface

! MIN (limited to 2-4 complex args, 2 mixed args)
  interface min
     module procedure min_cc
     module procedure min_cr
     module procedure min_rc
     module procedure min_ccc
     module procedure min_cccc
  end interface

! MINVAL
  interface minval
      module procedure minval_c
  end interface minval

! MAXVAL
  interface maxval
      module procedure maxval_c
  end interface maxval

! SIGN
  interface sign
     module procedure sign_cc
     module procedure sign_cca
     module procedure sign_cr
     module procedure sign_rc
  end interface

! DIM
  interface dim
     module procedure dim_cc
     module procedure dim_cr
     module procedure dim_rc
  end interface

! SINH
  interface sinh
     module procedure sinh_c
  end interface
  
! TAN
  interface tan
     module procedure tan_c
  end interface
  
! TANH
  interface tanh
     module procedure tanh_c
  end interface

! LOG10
  interface log10
     module procedure log10_c
  end interface

! NINT
  interface nint
     module procedure nint_c
  end interface

! EPSILON
  interface epsilon
     module procedure epsilon_c
  end interface

! <
  interface operator (<)
     module procedure lt_cc
     module procedure lt_cr
     module procedure lt_rc
     module procedure lt_ci
     module procedure lt_ic
  end interface

! <=
  interface operator (<=)
     module procedure le_cc
     module procedure le_cr
     module procedure le_rc
     module procedure le_ci
     module procedure le_ic
  end interface

! >
  interface operator (>)
     module procedure gt_cc
     module procedure gt_cr
     module procedure gt_rc
     module procedure gt_ci
     module procedure gt_ic
     module procedure gt_cac
     module procedure gt_car
  end interface

!! MIPSpro Compilers: Version 7.30 won't take .ge. and .eq..
!! But pgf90 on Linux doesn't complain, go figure.
!! It looks like a strict interpretation of FORTRAN should
!! not allow overloading of .eq. and .ne. since they already
!! have a definition for type complex, so define new operators
!! called .ceq., .cne. and, for MIPS, .cge.
!!
!! comment out (and uncomment) the appropriate versions for
!! your compiler
!!
! >= 
  interface operator (>=)
     module procedure ge_cc
     module procedure ge_cr
     module procedure ge_rc
     module procedure ge_ci
     module procedure ge_ic
  end interface
! interface operator (.cge.)
!    module procedure ge_cc
!    module procedure ge_rr
!    module procedure ge_ii
!    module procedure ge_aa
!    module procedure ge_cr
!    module procedure ge_rc
!    module procedure ge_ci
!    module procedure ge_ic
!    module procedure ge_ir
!    module procedure ge_ri
! end interface

! ==
!  interface operator (==)
!     module procedure eq_cc
!     module procedure eq_cr
!     module procedure eq_rc
!     module procedure eq_ci
!     module procedure eq_ic
!  end interface
  interface operator (.ceq.)
     module procedure eq_cc
     module procedure eq_rr
     module procedure eq_ii
     module procedure eq_iai
     module procedure eq_iaia
     module procedure eq_i8i8
     module procedure eq_i1i1
     module procedure eq_aa
     module procedure eq_cr
     module procedure eq_rc
     module procedure eq_ci
     module procedure eq_ic
     module procedure eq_ir
     module procedure eq_ri
  end interface

! /=
!  interface operator (/=)
!     module procedure ne_cc
!     module procedure ne_cr
!     module procedure ne_rc
!     module procedure ne_ci
!     module procedure ne_ic
!  end interface
  interface operator (.cne.)
     module procedure ne_cc
     module procedure ne_rr
     module procedure ne_ii
     module procedure ne_aa
     module procedure ne_cr
     module procedure ne_rc
     module procedure ne_ci
     module procedure ne_ic
     module procedure ne_ir
     module procedure ne_ri
  end interface

! floor
  interface floor
     module procedure floor_c
  end interface

contains

!******************************************************************************
!
!   Function definitions
!
!******************************************************************************

! ABS, intrinsic
  complex*16 function abs_c(val)
    complex*16, intent(in) :: val
    abs_c = val
    if (real(val) < 0) abs_c = cmplx(-real(val),-aimag(val))
    return
  end function abs_c

! COSD
!  complex*16 function cosd_c(z)
!    complex*16, intent(in) :: z
!    cosd_c = cos(z*3.14159265358979323846/180.)
!  end function cosd_c

! SIND
!  complex*16 function sind_c(z)
!    complex*16, intent(in) :: z
!    sind_c = sin(z*3.14159265358979323846/180.)
!  end function sind_c

! ACOS
  complex*16 function acos_c(z)
    complex*16, intent(in) :: z
!   acos_c = - cmplx(0., 1.)*log(z+sqrt(z**2-1.))
!   not general complex valued formula:
    acos_c = cmplx(acos(real(z)),-aimag(z)/sqrt(1.-real(z)**2))
    return
  end function acos_c

! ASIN
  complex*16 function asin_c(z)
    complex*16, intent(in) :: z
!   asin_c = - cmplx(0., 1.)*log(cmplx(0.,1.)*z+sqrt(1.-z**2))
!   not general complex valued formula:
    asin_c = cmplx(asin(real(z)),aimag(z)/sqrt(1.-real(z)**2))
    return
  end function asin_c

! ATAN
  complex*16 function atan_c(z)
    complex*16, intent(in) :: z
!   complex*16 z2
!   real*8 pi2, xans, yans, r, r2, x, y
!   pi2 = 2.0*atan(1.0)
!   r      = sqrt(real(z)**2+aimag(z)**2)
!   x      = real(z)
!   y      = aimag(z)
!   r2     = r*r
!   xans   = 0.5*atan2 (2.0*x, 1.0-r2)
!   yans   = 0.25*log((r2+2.0*y+1.0)/(r2-2.0*y+1.0))
!   atan_c = cmplx (xans, yans)
!   not general complex valued formula:
    atan_c = cmplx(atan(real(z)),aimag(z)/(1.+real(z)**2))
    return
  end function atan_c
  
! ATAN2
  complex*16 function atan2_cc(csn, ccs)
    complex*16, intent(in) :: csn, ccs
!   real*8 pi
!   pi = 4.0*atan(1.0)
!   if (sqrt(real(ccs)**2 + aimag(ccs)**2).eq.0.) then  ! abs orig
!     if (sqrt(real(csn)**2+aimag(csn)**2).eq.0.) then
!       atan2_cc = cmplx(0.0)
!     else
!       atan2_cc = cmplx(sign(0.5*pi,real(csn)), 0.0)
!     end if
!   else
!     atan2_cc = atan(csn/ccs)
!     if (real(ccs).lt.0.) atan2_cc = atan2_cc + pi
!     if (real(atan2_cc).gt.pi) atan2_cc = atan2_cc - 2.0*pi
!   end if
!   not general complex valued formula:
    real*8 a,b,c,d
    a=real(csn)
    b=aimag(csn)
    c=real(ccs)
    d=aimag(ccs)
    atan2_cc=cmplx(atan2(a,c),(c*b-a*d)/(a**2+c**2))
    return
  end function atan2_cc

! COSH
  complex*16 function cosh_c(z)
    complex*16, intent(in) :: z
!   complex*16 eplus, eminus
!   eplus = exp(z)
!   eminus = exp(z)
!   cosh_c = (eplus + eminus)/2.
!   not general complex valued formula:
    cosh_c=cmplx(cosh(real(z)),aimag(z)*sinh(real(z)))
    return
  end function cosh_c

! SINH
  complex*16 function sinh_c(z)
    complex*16, intent(in) :: z
!   complex*16 eplus, eminus
!   eplus = exp(z)
!   eminus = exp(z)
!   sinh_c = (eplus - eminus)/2.
!   not general complex valued formula:
    sinh_c=cmplx(sinh(real(z)),aimag(z)*cosh(real(z)))
    return
  end function sinh_c

! TAN
  complex*16 function tan_c(z)
    complex*16, intent(in) :: z
!   complex*16 eiplus, eiminus
!   eiplus = exp(cmplx(0.,1.)*z)
!   eiminus = exp(-cmplx(0.,1.)*z)
!   tan_c = cmplx(0.,1.)*(eiminus - eiplus)/(eiplus + eiminus)
!   not general complex valued formula:
    tan_c=cmplx(tan(real(z)),aimag(z)/cos(real(z))**2)
    return
  end function tan_c
  
! TANH
  complex*16 function tanh_c(a)
    complex*16, intent(in) :: a
!   complex*16 eplus, eminus
!   if(real(a) > 50)then
!     tanh_c = 1.
!   else
!     eplus = exp(a)
!     eminus = exp(-a)
!     tanh_c = (eplus - eminus)/(eplus + eminus)
!   end if
!   not general complex valued formula:
    tanh_c=cmplx(tanh(real(a)),aimag(a)/cosh(real(a))**2)
    return
  end function tanh_c

! MAX, intrinsic
  complex*16 function max_cc(val1, val2)
    complex*16, intent(in) :: val1, val2
    if (real(val1) > real(val2)) then
      max_cc = val1
    else
      max_cc = val2
    endif
    return
  end function max_cc
  complex*16 function max_cr(val1, val2)
    complex*16, intent(in) :: val1    
    real*8, intent(in) :: val2    
    if (real(val1) > val2) then
      max_cr = val1
    else
      max_cr = cmplx(val2, 0.)
    endif
    return
  end function max_cr
  complex*16 function max_rc(val1, val2)
    real*8, intent(in) :: val1
    complex*16, intent(in) :: val2
    if (val1 > real(val2)) then
      max_rc = cmplx(val1, 0.)
    else
      max_rc = val2
    endif
    return
  end function max_rc
  complex*16 function max_ccc(val1, val2, val3)
    complex*16, intent(in) :: val1, val2, val3
    if (real(val1) > real(val2)) then
      max_ccc = val1
    else
      max_ccc = val2
    endif
    if (real(val3) > real(max_ccc)) then
      max_ccc = val3
    endif
    return
  end function max_ccc
  function max_cccc(val1, val2, val3, val4)
    complex*16, intent(in) :: val1, val2, val3, val4
    complex*16 max_cccc
    complex*16 max_cccc2
    if (real(val1) > real(val2)) then
      max_cccc = val1
    else
      max_cccc = val2
    endif
    if (real(val3) > real(val4)) then
      max_cccc2 = val3
    else
      max_cccc2 = val4
    endif
    if (real(max_cccc2) > real(max_cccc)) then
      max_cccc = max_cccc2
    endif
    return
  end function max_cccc

! MIN, intrinsic
  complex*16 function min_cc(val1, val2)
    complex*16, intent(in) :: val1, val2
    if (real(val1) < real(val2)) then
      min_cc = val1
    else
      min_cc = val2
    endif
    return
  end function min_cc
  complex*16 function min_cr(val1, val2)
    complex*16, intent(in) :: val1    
    real*8, intent(in) :: val2    
    if (real(val1) < val2) then
      min_cr = val1
    else
      min_cr = cmplx(val2, 0.)
    endif
    return
  end function min_cr
  complex*16 function min_rc(val1, val2)
    real*8, intent(in) :: val1
    complex*16, intent(in) :: val2
    if (val1 < real(val2)) then
      min_rc = cmplx(val1, 0.)
    else
      min_rc = val2
    endif
    return
  end function min_rc
  complex*16 function min_ccc(val1, val2, val3)
    complex*16, intent(in) :: val1, val2, val3
    if (real(val1) < real(val2)) then
      min_ccc = val1
    else
      min_ccc = val2
    endif
    if (real(val3) < real(min_ccc)) then
      min_ccc = val3
    endif
    return
  end function min_ccc
  function min_cccc(val1, val2, val3, val4)
    complex*16, intent(in) :: val1, val2, val3, val4
    complex*16 min_cccc
    complex*16 min_cccc2
    if (real(val1) < real(val2)) then
      min_cccc = val1
    else
      min_cccc = val2
    endif
    if (real(val3) < real(val4)) then
      min_cccc2 = val3
    else
      min_cccc2 = val4
    endif
    if (real(min_cccc2) < real(min_cccc)) then
      min_cccc = min_cccc2
    endif
    return
  end function min_cccc

! MINVAL: minimum of an array
! Assumes a 1D array!
  complex*16 function minval_c(z)
    complex*16, intent(in) :: z(:)
    minval_c = cmplx(minval(real(z)), aimag(z(minloc(real(z),dim=1))))
  end function minval_c

! MAXVAL: maximum of an array
! Assumes a 1D array!
  complex*16 function maxval_c(z)
    complex*16, intent(in) :: z(:)
    maxval_c = cmplx(maxval(real(z)), aimag(z(maxloc(real(z),dim=1))))
  end function maxval_c

!! MINLOC: location of minimum in an array
!  complex*16 function minloc_c(z)
!    complex*16, intent(in) :: z(:)
!    !integer n
!    !n = size(z)
!    minloc_c = minloc(real(z))
!  end function minval_c

  
! SIGN, intrinsic, assume that val1 is always a complex*16
!                  in reality could be int
  complex*16 function sign_cc(val1, val2)
    complex*16, intent(in) :: val1, val2
    real*8  sign
    if (real(val2) < 0.) then
      sign = -1.
    else
      sign = 1.
    endif
    sign_cc = sign * val1
    return
  end function sign_cc
  function sign_cca(val1, val2) ! NEW, not verified
    complex*16, intent(in) :: val1
    complex*16, intent(in) :: val2(:)
    complex*16 sign_cca(size(val2))
    real*8 sign
    integer i, n
    n = size(val2)
    do i = 1, n
       if (real(val2(i)) < 0.) then
          sign = -1.
       else
          sign = 1.
       endif
       sign_cca(i) = sign * val1
    enddo
    return
  end function sign_cca
  complex*16 function sign_cr(val1, val2)
    complex*16, intent(in) :: val1
    real*8, intent(in) :: val2
    real*8 sign
    if (real(val2) < 0.) then
      sign = -1.
    else
      sign = 1.
    endif
    sign_cr = sign * val1
    return
  end function sign_cr
  complex*16 function sign_rc(val1, val2)
    real*8, intent(in) :: val1
    complex*16, intent(in) :: val2
    real*8 sign
    if (real(val2) < 0.) then
      sign = -1.
    else
      sign = 1.
    endif
    sign_rc = sign * val1
    return
  end function sign_rc

! DIM, intrinsic
  complex*16 function dim_cc(val1, val2)
    complex*16, intent(in) :: val1, val2
    if (val1 > val2) then
      dim_cc = val1 - val2
    else
      dim_cc = cmplx(0., 0.)
    endif
    return
  end function dim_cc
  complex*16 function dim_cr(val1, val2)
    complex*16, intent(in) :: val1
    real*8, intent(in) :: val2
    if (val1 > val2) then
      dim_cr = val1 - cmplx(val2, 0.)
    else
      dim_cr = cmplx(0., 0.)
    endif
    return
  end function dim_cr
  complex*16 function dim_rc(val1, val2)
    real*8, intent(in) :: val1
    complex*16, intent(in) :: val2
    if (val1 > val2) then
      dim_rc = cmplx(val1, 0.) - val2
    else
      dim_rc = cmplx(0., 0.)
    endif
    return
  end function dim_rc
  
! LOG10
  complex*16 function log10_c(z)
    complex*16, intent(in) :: z
    log10_c=log(z)/log((10.0,0.0))
  end function log10_c

! NINT
  integer function nint_c(z)
    complex*16, intent(in) :: z
    nint_c = nint(real(z))
  end function nint_c

! EPSILON !! bad news ulness compiled with -r8
  complex*16 function epsilon_c(z)
    complex*16, intent(in) :: z
    epsilon_c=epsilon(real(z))
  end function epsilon_c

! <, .lt.
  logical function lt_cc(lhs, rhs)
    complex*16, intent(in) :: lhs, rhs
    lt_cc = real(lhs) < real(rhs)
  end function lt_cc
  logical function lt_cr(lhs, rhs)
    complex*16, intent(in) :: lhs
    real*8, intent(in) :: rhs
    lt_cr = real(lhs) < rhs
  end function lt_cr
  logical function lt_rc(lhs, rhs)
    real*8, intent(in) :: lhs
    complex*16, intent(in) :: rhs
    lt_rc = lhs < real(rhs)
  end function lt_rc
  logical function lt_ci(lhs, rhs)
    complex*16, intent(in) :: lhs
    integer, intent(in) :: rhs
    lt_ci = real(lhs) < rhs
  end function lt_ci
  logical function lt_ic(lhs, rhs)
    integer, intent(in) :: lhs
    complex*16, intent(in) :: rhs
    lt_ic = lhs < real(rhs)
  end function lt_ic

! <=, .le.
  logical function le_cc(lhs, rhs)
    complex*16, intent(in) :: lhs, rhs
    le_cc = real(lhs) <= real(rhs)
  end function le_cc
  logical function le_cr(lhs, rhs)
    complex*16, intent(in) :: lhs
    real*8, intent(in) :: rhs
    le_cr = real(lhs) <= rhs
  end function le_cr
  logical function le_rc(lhs, rhs)
    real*8, intent(in) :: lhs
    complex*16, intent(in) :: rhs
    le_rc = lhs <= real(rhs)
  end function le_rc
  logical function le_ci(lhs, rhs)
    complex*16, intent(in) :: lhs
    integer, intent(in) :: rhs
    le_ci = real(lhs) <= rhs
  end function le_ci
  logical function le_ic(lhs, rhs)
    integer, intent(in) :: lhs
    complex*16, intent(in) :: rhs
    le_ic = lhs <= real(rhs)
  end function le_ic

! >, .gt.
  logical function gt_cc(lhs, rhs)
    complex*16, intent(in) :: lhs, rhs
    gt_cc = real(lhs) > real(rhs)
  end function gt_cc
  logical function gt_cr(lhs, rhs)
    complex*16, intent(in) :: lhs
    real*8, intent(in) :: rhs
    gt_cr = real(lhs) > rhs
  end function gt_cr
  logical function gt_rc(lhs, rhs)
    real*8, intent(in) :: lhs
    complex*16, intent(in) :: rhs
    gt_rc = lhs > real(rhs)
  end function gt_rc
  logical function gt_ci(lhs, rhs)
    complex*16, intent(in) :: lhs
    integer, intent(in) :: rhs
    gt_ci = real(lhs) > rhs
  end function gt_ci
  logical function gt_ic(lhs, rhs)
    integer, intent(in) :: lhs
    complex*16, intent(in) :: rhs
    gt_ic = lhs > real(rhs)
  end function gt_ic
!  function gt_caca(lhs, rhs) ! Arrays
!    complex*16, intent(in) :: lhs(:), rhs(:)
!    logical gt_caca(size(lhs))
!    integer n
!    n = size(lhs)
!    gt_caca = real(lhs) > real(rhs)
!  end function gt_caca
  function gt_cac(lhs, rhs) ! Arrays
    complex*16, intent(in) :: lhs(:)
    complex*16, intent(in) :: rhs
    logical gt_cac(size(lhs))
    integer n
    n = size(lhs)
    gt_cac = real(lhs) > real(rhs)
  end function gt_cac
  function gt_car(lhs, rhs) ! Arrays
    complex*16, intent(in) :: lhs(:)
    real*8, intent(in) :: rhs
    logical gt_car(size(lhs))
    integer n
    n = size(lhs)
    gt_car = real(lhs) > rhs
  end function gt_car

!! here are the redefined ones:
! >=, .ge.
  logical function ge_cc(lhs, rhs)
    complex*16, intent(in) :: lhs, rhs
    ge_cc = real(lhs) >= real(rhs)
  end function ge_cc
  logical function ge_rr(lhs, rhs)
    real*8, intent(in) :: lhs, rhs
    ge_rr = lhs >= rhs
  end function ge_rr
  logical function ge_ii(lhs, rhs)
    integer, intent(in) :: lhs, rhs
    ge_ii = lhs >= rhs
  end function ge_ii
  logical function ge_aa(lhs, rhs)
    character(len=*), intent(in) :: lhs, rhs
    ge_aa = lhs >= rhs
  end function ge_aa
  logical function ge_cr(lhs, rhs)
    complex*16, intent(in) :: lhs
    real*8, intent(in) :: rhs
    ge_cr = real(lhs) >= rhs
  end function ge_cr
  logical function ge_rc(lhs, rhs)
    real*8, intent(in) :: lhs
    complex*16, intent(in) :: rhs
    ge_rc = lhs >= real(rhs)
  end function ge_rc
  logical function ge_ci(lhs, rhs)
    complex*16, intent(in) :: lhs
    integer, intent(in) :: rhs
    ge_ci = real(lhs) >= rhs
  end function ge_ci
  logical function ge_ic(lhs, rhs)
    integer, intent(in) :: lhs
    complex*16, intent(in) :: rhs
    ge_ic = lhs >= real(rhs)
  end function ge_ic
  logical function ge_ir(lhs, rhs)
    integer, intent(in) :: lhs
    real*8, intent(in) :: rhs
    ge_ir = lhs >= rhs
  end function ge_ir
  logical function ge_ri(lhs, rhs)
    real*8, intent(in) :: lhs
    integer, intent(in) :: rhs
    ge_ri = lhs >= rhs
  end function ge_ri

! ==, .eq.
  logical function eq_cc(lhs, rhs)
    complex*16, intent(in) :: lhs, rhs
    eq_cc = real(lhs) == real(rhs)
  end function eq_cc
  logical function eq_rr(lhs, rhs)
    real*8, intent(in) :: lhs, rhs
    eq_rr = lhs == rhs
  end function eq_rr
  logical function eq_ii(lhs, rhs)
    integer, intent(in) :: lhs, rhs
    eq_ii = lhs == rhs
  end function eq_ii
  ! lhs and rhs are rank 1 integer arrays
  function eq_iaia(lhs, rhs) 
    integer, intent(in) :: lhs(:), rhs(:)
    logical eq_iaia(size(lhs))
    eq_iaia = lhs == rhs
  end function eq_iaia
  ! lhs is a rank 3 integer array
  function eq_iai(lhs, rhs) 
    integer, intent(in) :: lhs(:,:,:)
    integer, intent(in) :: rhs
    logical eq_iai(size(lhs,1), size(lhs,2), size(lhs,3))
    eq_iai = lhs == rhs
  end function eq_iai
  logical function eq_i8i8(lhs, rhs)
    integer(kind = 8), intent(in) :: lhs, rhs
    eq_i8i8 = lhs == rhs
  end function eq_i8i8
  logical function eq_i1i1(lhs, rhs)
    integer(kind = 1), intent(in) :: lhs, rhs
    eq_i1i1 = lhs == rhs
  end function eq_i1i1
  logical function eq_aa(lhs, rhs)
    character(len=*), intent(in) :: lhs, rhs
    eq_aa = lhs == rhs
  end function eq_aa
  logical function eq_cr(lhs, rhs)
    complex*16, intent(in) :: lhs
    real*8, intent(in) :: rhs
    eq_cr = real(lhs) == rhs
  end function eq_cr
  logical function eq_rc(lhs, rhs)
    real*8, intent(in) :: lhs
    complex*16, intent(in) :: rhs
    eq_rc = lhs == real(rhs)
  end function eq_rc
  logical function eq_ci(lhs, rhs)
    complex*16, intent(in) :: lhs
    integer, intent(in) :: rhs
    eq_ci = real(lhs) == rhs
  end function eq_ci
  logical function eq_ic(lhs, rhs)
    integer, intent(in) :: lhs
    complex*16, intent(in) :: rhs
    eq_ic = lhs == real(rhs)
  end function eq_ic
  logical function eq_ir(lhs, rhs)
    integer, intent(in) :: lhs
    real*8, intent(in) :: rhs
    eq_ir = lhs == rhs
  end function eq_ir
  logical function eq_ri(lhs, rhs)
    real*8, intent(in) :: lhs
    integer, intent(in) :: rhs
    eq_ri = lhs == rhs
  end function eq_ri

! /=, .ne.
  logical function ne_cc(lhs, rhs)
    complex*16, intent(in) :: lhs, rhs
    ne_cc = real(lhs) /= real(rhs)
  end function ne_cc
  logical function ne_rr(lhs, rhs)
    real*8, intent(in) :: lhs, rhs
    ne_rr = lhs /= rhs
  end function ne_rr
  logical function ne_ii(lhs, rhs)
    integer, intent(in) :: lhs, rhs
    ne_ii = lhs /= rhs
  end function ne_ii
  logical function ne_aa(lhs, rhs)
    character(len=*), intent(in) :: lhs, rhs
    ne_aa = lhs /= rhs
  end function ne_aa
  logical function ne_cr(lhs, rhs)
    complex*16, intent(in) :: lhs
    real*8, intent(in) :: rhs
    ne_cr = real(lhs) /= rhs
  end function ne_cr
  logical function ne_rc(lhs, rhs)
    real*8, intent(in) :: lhs
    complex*16, intent(in) :: rhs
    ne_rc = lhs /= real(rhs)
  end function ne_rc
  logical function ne_ci(lhs, rhs)
    complex*16, intent(in) :: lhs
    integer, intent(in) :: rhs
    ne_ci = real(lhs) /= rhs
  end function ne_ci
  logical function ne_ic(lhs, rhs)
    integer, intent(in) :: lhs
    complex*16, intent(in) :: rhs
    ne_ic = lhs /= real(rhs)
  end function ne_ic
  logical function ne_ir(lhs, rhs)
    integer, intent(in) :: lhs
    real*8, intent(in) :: rhs
    ne_ir = lhs /= rhs
  end function ne_ir
  logical function ne_ri(lhs, rhs)
    real*8, intent(in) :: lhs
    integer, intent(in) :: rhs
    ne_ri = lhs /= rhs
  end function ne_ri

! floor: the largest integer less than or equal to the argument
  function floor_c(z)
    complex*16, intent(in) :: z(:)
    complex*16 floor_c(size(z))
    integer n
    n = size(z)
    floor_c = floor(real(z(1:n)))
  end function floor_c

end module complexify

