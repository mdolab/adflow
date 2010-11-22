!
!     ******************************************************************
!     *                                                                *
!     * File:          computeForceCouplingAdj.f90                     *
!     * Author:        Gaetan Kenway                                   *
!     * Starting date: 09-27-2010                                      *
!     * Last modified: 09-27-2010                                      *
!     *                                                                *
!     ******************************************************************

subroutine computeForceCouplingAdj(force,moment,pts,wAdj,refPoint,fact,&
     iBeg,iEnd,jBeg,jEnd,iNode,jNode,righthanded)
  
  use constants
  use flowvarrefstate      !nw


  implicit none

  !     Subroutine arguments.

  real(kind=realType), intent(in) :: pts(3,3,3)
  real(kind=realType), intent(in) :: wAdj(2,2,2,nw)
  real(kind=realType), intent(in) :: fact
  real(kind=realType), intent(in) :: refPoint(3)
  integer(kind=intType),intent(in):: iBeg,iEnd,jBeg,jEnd,iNode,jNode
  logical , intent(in) :: righthanded
  real(kind=realType), intent(out) :: force(3)
  real(kind=realType), intent(out) :: moment(3)

  ! Local Variables
  integer(kind=intType):: i,j,k,l,kk
  real(kind=realType)  :: pAdj(3,2,2)
  real(kind=realType)  :: normAdj(3,2,2)

 ! Compute the pressure in the 2x2x2 block
  padj = 0.0

  call computeForceCouplingPressureAdj(wAdj,pAdj(2:3,:,:))
  call getSurfaceNormalsCouplingAdj(pts,normAdj,rightHanded)

  ! Only BC Euler wall boundary conditions are implemented. 
  call bcEulerWallForceCouplingAdj(wAdj,pAdj)

  ! For NS probably will need to call 
  ! call bcNSWallAdiabaticForceCouplingAdj(wAdj,pAdj)

  ! Also for NS calculation, will need to compute the tau on the
  ! subfaces, with 

  ! call viscousFluxCouplingAdj(     )

  call forcesCouplingAdj(pAdj,pts,normAdj,refPoint,force,moment,fact,&
       iBeg,iEnd,jBeg,jEnd,iNode,jNode)

end subroutine computeForceCouplingAdj
