!
!      ******************************************************************
!      *                                                                *
!      * File:          computeLeastSquaresRegression.f90               *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 10-14-2009                                      *
!      * Last modified: 10-14-2009                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine computeLeastSquaresRegression(y,x,npts,m,b)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Computes the slope of best fit for a set of x,y data of length *
  !      * npts                                                           *
  !      *                                                                *
  !      ******************************************************************
  !
  use precision
  implicit none
  !Subroutine arguments 
  integer(kind=intType)::npts
  real(kind=realType),dimension(npts)  :: x,y
  real(kind=realType)::m,b

  !local variables
  real(kind=realType)::sumx,sumy,sumx2,sumxy
  integer(kind=intType)::i

  !begin execution
  sumx=0.0
  sumy=0.0
  sumx2=0.0
  sumxy=0.0
  do i = 1,npts

     sumx=sumx+x(i)
     sumy=sumy+y(i)
     sumx2=sumx2+x(i)*x(i)
     sumxy=sumxy+x(i)*y(i)
  enddo

  m = ((npts*sumxy)-(sumy*sumx))/((npts*sumx2)-(sumx)**2)
  b = (sumy*sumx2-(sumx*sumxy))/((npts*sumx2)-(sumx)**2)

end subroutine computeLeastSquaresRegression
