!
!      ******************************************************************
!      *                                                                *
!      * File:          forcesCouplingAdj.f90                           *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader           *
!      * Starting date: 08-17-2008                                      *
!      * Last modified: 08-17-2008                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine forcesCouplingAdj(pAdj,pts,normAdj,force,fact,&
     iBeg,iEnd,jBeg,jEnd,iNode,jNode)
  !
  
  use flowVarRefState


  implicit none

  ! Subroutine Arguments
  real(kind=realType), intent(in) :: pAdj(3,2,2)
  real(kind=realType), intent(in) :: pts(3,3,3)
  real(kind=realType), intent(in) :: normAdj(3,2,2)
  real(kind=realType), intent(out):: force(3)
  real(kind=realType), intent(in) :: fact
  integer(kind=intType),intent(in):: iBeg,iEnd,jBeg,jEnd,iNode,jNode

  ! Local Variables
  integer(kind=intTYpe) :: i,j
  real(kind=realType) :: pp,scaleDim
  scaleDim = pRef
  force(:) = 0.0

  ! Force is the contribution of each of 4 cells

  do j=1,2
     do i=1,2

        if (.not.(iNode+i-2 < iBeg .or. iNode+i-1 > iEnd .or. &
                  jNode+j-2 < jBeg .or. jNode+j-1 > jEnd)) then

           ! Calculate the Pressure
           pp = half*(pAdj(1,i,j) + pAdj(2,i,j)) - Pinf
           pp = fourth*fact*scaleDim*pp
           
           force(:) = force(:) + pp*normAdj(:,i,j)

         end if
      end do
   end do
 end subroutine forcesCouplingAdj
