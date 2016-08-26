!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of getdirvector in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: alpha beta winddirection
!   with respect to varying inputs: alpha beta
!
!      file:          getdirvector.f90                                
!      author:        andre c. marta                                  
!      starting date: 10-25-2005                                      
!      last modified: 10-26-2006                                      
!
subroutine getdirvector_b(refdirection, alpha, alphad, beta, betad, &
& winddirection, winddirectiond, liftindex)
!(xb,yb,zb,alpha,beta,xw,yw,zw)
!
!      convert the angle of attack and side slip angle to wind axes.  
!      the components of the wind direction vector (xw,yw,zw) are     
!      computed given the direction angles in radians and the body    
!      direction by performing two rotations on the original          
!      direction vector:                                              
!        1) rotation about the zb or yb-axis: alpha clockwise (cw)    
!           (xb,yb,zb) -> (x1,y1,z1)                                  
!        2) rotation about the yl or z1-axis: beta counter-clockwise  
!           (ccw)  (x1,y1,z1) -> (xw,yw,zw)                           
!         input arguments:                                            
!            alpha    = angle of attack in radians                    
!            beta     = side slip angle in radians                    
!            refdirection = reference direction vector                
!         output arguments:                                           
!            winddirection = unit wind vector in body axes            
!
  use constants
  use utils_b, only : terminate
  implicit none
!
!     subroutine arguments.
!
  real(kind=realtype), dimension(3), intent(in) :: refdirection
  real(kind=realtype) :: alpha, beta
  real(kind=realtype) :: alphad, betad
  real(kind=realtype), dimension(3) :: winddirection
  real(kind=realtype), dimension(3) :: winddirectiond
  integer(kind=inttype) :: liftindex
!
!     local variables.
!
  real(kind=realtype) :: rnorm, x1, y1, z1, xbn, ybn, zbn, xw, yw, zw
  real(kind=realtype) :: x1d, y1d, z1d, xbnd, ybnd, zbnd, xwd, ywd, zwd
  real(kind=realtype) :: tmp
  real(kind=realtype) :: tmpd
  intrinsic sqrt
  integer :: branch
!      begin execution.                                               
!
! normalize the input vector.
  rnorm = sqrt(refdirection(1)**2 + refdirection(2)**2 + refdirection(3)&
&   **2)
  xbn = refdirection(1)/rnorm
  ybn = refdirection(2)/rnorm
  zbn = refdirection(3)/rnorm
!!$      ! compute the wind direction vector.
!!$
!!$      ! 1) rotate alpha radians cw about y-axis
!!$      !    ( <=> rotate y-axis alpha radians ccw)
!!$
!!$      call vectorrotation(x1, y1, z1, 2, alpha, xbn, ybn, zbn)
!!$
!!$      ! 2) rotate beta radians ccw about z-axis
!!$      !    ( <=> rotate z-axis -beta radians ccw)
!!$
!!$      call vectorrotation(xw, yw, zw, 3, -beta, x1, y1, z1)
  if (liftindex .eq. 2) then
! compute the wind direction vector.aerosurf axes different!!
! 1) rotate alpha radians cw about z-axis
!    ( <=> rotate z-axis alpha radians ccw)
    tmp = -alpha
    call vectorrotation(x1, y1, z1, 3, tmp, xbn, ybn, zbn)
! 2) rotate beta radians ccw about y-axis
!    ( <=> rotate z-axis -beta radians ccw)
    tmp = -beta
    call pushcontrol2b(0)
  else if (liftindex .eq. 3) then
! compute the wind direction vector.aerosurf axes different!!
! 1) rotate alpha radians cw about z-axis
!    ( <=> rotate z-axis alpha radians ccw)
    call vectorrotation(x1, y1, z1, 2, alpha, xbn, ybn, zbn)
! 2) rotate beta radians ccw about y-axis
!    ( <=> rotate z-axis -beta radians ccw)
    call pushcontrol2b(1)
  else
    call pushcontrol2b(2)
  end if
  zwd = winddirectiond(3)
  winddirectiond(3) = 0.0_8
  ywd = winddirectiond(2)
  winddirectiond(2) = 0.0_8
  xwd = winddirectiond(1)
  call popcontrol2b(branch)
  if (branch .eq. 0) then
    tmpd = 0.0_8
    call vectorrotation_b(xw, xwd, yw, ywd, zw, zwd, 2, tmp, tmpd, x1, &
&                   x1d, y1, y1d, z1, z1d)
    betad = betad - tmpd
    tmp = -alpha
    tmpd = 0.0_8
    call vectorrotation_b(x1, x1d, y1, y1d, z1, z1d, 3, tmp, tmpd, xbn, &
&                   xbnd, ybn, ybnd, zbn, zbnd)
    alphad = alphad - tmpd
  else if (branch .eq. 1) then
    call vectorrotation_b(xw, xwd, yw, ywd, zw, zwd, 3, beta, betad, x1&
&                   , x1d, y1, y1d, z1, z1d)
    call vectorrotation_b(x1, x1d, y1, y1d, z1, z1d, 2, alpha, alphad, &
&                   xbn, xbnd, ybn, ybnd, zbn, zbnd)
  end if
end subroutine getdirvector_b
