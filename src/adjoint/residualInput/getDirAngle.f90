!
!     ******************************************************************
!     *                                                                *
!     * File:          getDirAngle.f90                                 *
!     * Author:        Andre C. Marta,C.A.(Sandy) Mader                *
!     * Starting date: 10-25-2005                                      *
!     * Last modified: 06-13-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine getDirAngle(freeStreamAxis,liftAxis,liftIndex,alpha,beta)
!
!     ******************************************************************
!     *                                                                *
!     * Convert the wind axes to angle of attack and side slip angle.  *
!     * The direction angles alpha and beta are computed given the     *
!     * components of the wind direction vector (freeStreamAxis), the  *
!     * lift direction vector (liftAxis) and assuming that the         *
!     * body direction (xb,yb,zb) is in the default ijk coordinate     *
!     * system. The rotations are determined by first determining      *
!     * whether the lift is primarily in the j or k direction and then *
!     * determining the angles accordingly.                            *
!     * direction vector:                                              *
!     *   1) Rotation about the zb or yb -axis: alpha clockwise (CW)   *
!     *      (xb,yb,zb) -> (x1,y1,z1)                                  *
!     *                                                                *
!     *   2) Rotation about the yl or z1 -axis: beta counter-clockwise *
!     *      (CCW) (x1,y1,z1) -> (xw,yw,zw)                            *
!     *                                                                *
!     *    input arguments:                                            *
!     *       freeStreamAxis = wind vector in body axes                *
!     *       liftAxis       = lift direction vector in body axis      *       
!     *    output arguments:                                           *
!     *       alpha    = angle of attack in radians                    *
!     *       beta     = side slip angle in radians                    *
!     *                                                                *
!     ******************************************************************
!
      use constants
      implicit none
!
!     Subroutine arguments.
!
!      real(kind=realType), intent(in)  :: xw, yw, zw
      real(kind=realType), dimension(3),intent(in) :: freeStreamAxis
      real(kind=realType), dimension(3),intent(in) :: liftAxis
      real(kind=realType), intent(out) :: alpha, beta
      integer(kind=intType), intent(out)::liftIndex
!
!     Local variables.
!
      real(kind=realType) :: rnorm
      integer(kind=intType)::flowIndex,i
      real(kind=realType), dimension(3) :: freeStreamAxisNorm
      integer(kind=intType), dimension(1) :: temp
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
      
      !print *,'input',freeStreamAxis,liftAxis,liftIndex,alpha,beta
      ! Check that the dominant free stream direction is x 
      !print *,'in get angle'
      temp = maxloc(freeStreamAxis)
      flowIndex = temp(1)
      if (flowIndex .ne. 1) then
         call terminate('getDirAngle', 'Dominant Flow not in +ve x direction')
      endif
      !print *,'flow Index'
      ! Determine the dominant lift direction
      temp = maxloc(liftAxis)
      !print *,'temp',temp,temp(1)
      liftIndex = temp(1)
      
      ! Normalize the freeStreamDirection vector.
      !print *,'liftIndex1'
      rnorm = sqrt( freeStreamAxis(1)**2 + freeStreamAxis(2)**2 + freeStreamAxis(3)**2 )
      do i =1,3
         freeStreamAxisNorm(i) = freeStreamAxis(i)/rnorm
      enddo
      !print *,'liftIndex',liftIndex
      if (liftIndex == 2) then
         ! different coordinate system for aerosurf
         ! Wing is in z- direction
         ! Compute angle of attack alpha.
      
         alpha = asin(freeStreamAxisNorm(2))
         
         ! Compute side-slip angle beta.
         
         beta  = -atan2(freeStreamAxisNorm(3),freeStreamAxisNorm(1))


      elseif (liftIndex ==3) then
         ! Wing is in y- direction
         
         ! Compute angle of attack alpha.
         
         alpha = asin(freeStreamAxisNorm(3))
         
         ! Compute side-slip angle beta.
         
         beta  = atan2(freeStreamAxisNorm(2),freeStreamAxisNorm(1))
      else
         call terminate('getDirAngle', 'Invalid Lift Direction')
      endif

      !print *,'end get angle'
      

      end subroutine getDirAngle
