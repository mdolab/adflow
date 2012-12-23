!
!     ******************************************************************
!     *                                                                *
!     * File:          newtonStep.f90                                  *
!     * Author:        Juan J. Alonso                                  *
!     * Starting date: 04-20-2006                                      *
!     * Last modified: 04-25-2006                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine newtonStep(hess,grad,step,iErr)
!
!     ******************************************************************
!     *                                                                *
!     * Compute the Newton step given by the Hessian matrix and the    *
!     * gradient vector of the distance squared function.              *
!     *                                                                *
!     ******************************************************************
!
      use precision

      implicit none
!
!     Subroutine arguments.
!
      real(kind=realType), dimension(3),   intent(in)  :: grad
      real(kind=realType), dimension(3),  intent(out)  :: step
      real(kind=realType), dimension(3,3), intent(in)  :: hess

      integer(kind=intType), intent(out) :: iErr

!
!     Local variables.
!
      real(kind=realType) :: determinant

!
!     ******************************************************************
!     *                                                                *
!     * Compute the Newton step as the solution of the problem         *
!     *                                                                *
!     * Hessian . step = -grad                                         *
!     *                                                                *
!     * using simple Cramer's rule                                     *
!     *                                                                *
!     ******************************************************************
!
      ! Compute the determinant of the Hessian Matrix

      determinant = hess(1,1)*(hess(2,2)*hess(3,3)-hess(2,3)*hess(3,2)) &
                   +hess(1,2)*(hess(2,3)*hess(3,1)-hess(2,1)*hess(3,3)) &
                   +hess(1,3)*(hess(2,1)*hess(3,2)-hess(2,2)*hess(3,1))

      ! First component of the step

      step(1) = grad(1)  *(hess(2,2)*hess(3,3)-hess(2,3)*hess(3,2)) &
               +hess(1,2)*(hess(2,3)*grad(3)  -grad(2)  *hess(3,3)) &
               +hess(1,3)*(grad(2)  *hess(3,2)-hess(2,2)*grad(3)  )
      step(1) = -step(1)/determinant

      ! Second component of the step

      step(2) = hess(1,1)*(grad(2)  *hess(3,3)-hess(2,3)*grad(3)  ) &
               +grad(1)  *(hess(2,3)*hess(3,1)-hess(2,1)*hess(3,3)) &
               +hess(1,3)*(hess(2,1)*grad(3)  -grad(2)  *hess(3,1))
      step(2) = -step(2)/determinant

      ! First component of the step

      step(3) = hess(1,1)*(hess(2,2)*grad(3)  -grad(2)  *hess(3,2)) &
               +hess(1,2)*(grad(2)  *hess(3,1)-hess(2,1)*grad(3)  ) &
               +grad(1)  *(hess(2,1)*hess(3,2)-hess(2,2)*hess(3,1))
      step(3) = -step(3)/determinant

      ! Return iErr = 0 for the time being.

      iErr = 0

      return

      end subroutine newtonStep
