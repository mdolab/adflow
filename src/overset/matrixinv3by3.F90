!
!     ******************************************************************  
!     * Plain invert of A to solve Ax=B                                *
!     *                                                                *
!     * a(3,3) --- matrix to be inverted                               *
!     * a(:,4) --- contains the right hand side vector, also stores    *
!     *            the final solution                                  *
!     *                                                                *
!     ******************************************************************

subroutine matrixinv3by3(a, ok_flag)

  use precision

  implicit none

  ! Input variables
  real(kind=realType), intent(inout) :: a(3, 4)
  logical,             intent(out)   :: ok_flag   

  ! Working variables
  integer(kind=intType) :: n
  real(kind=realType) :: det, cofactor(3,3), ainv(3, 3), rin(3,1), rout(3,1), eps

  eps = 1.0e-10
  ainv = 0.d0

  det =   a(1,1)*a(2,2)*a(3,3)  &
        - a(1,1)*a(2,3)*a(3,2)  &
        - a(1,2)*a(2,1)*a(3,3)  &
        + a(1,2)*a(2,3)*a(3,1)  &
        + a(1,3)*a(2,1)*a(3,2)  &
        - a(1,3)*a(2,2)*a(3,1)
 
   
   if (abs(det) <= eps) then
      ainv = 0.0d0
      ok_flag = .false.
      return
   end if
   
   cofactor(1,1) = +(a(2,2)*a(3,3)-a(2,3)*a(3,2))
   cofactor(1,2) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
   cofactor(1,3) = +(a(2,1)*a(3,2)-a(2,2)*a(3,1))
   cofactor(2,1) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
   cofactor(2,2) = +(a(1,1)*a(3,3)-a(1,3)*a(3,1))
   cofactor(2,3) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
   cofactor(3,1) = +(a(1,2)*a(2,3)-a(1,3)*a(2,2))
   cofactor(3,2) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
   cofactor(3,3) = +(a(1,1)*a(2,2)-a(1,2)*a(2,1))
   
   ainv = transpose(cofactor) / det
   
   ok_flag = .true.
   
   do n=1, 3
      rin(n,1) = a(n, 4)
   end do

   !save solution on a(:,4)
   rout = matmul(ainv,rin)

   do n=1, 3
      a(n, 4) = rout(n,1)
   end do

end subroutine

