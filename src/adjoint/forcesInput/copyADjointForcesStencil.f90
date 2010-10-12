!
!     ******************************************************************
!     *                                                                *
!     * File:          copyADjointForcesStencil.f90                    *
!     * Author:        C.A.(Sandy) Mader                               *
!     *                Seongim Choi                                    *
!     * Starting date: 01-15-2007                                      *
!     * Last modified: 05-06-2008                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine copyADjointForcesStencil(wAdj,xAdj,alphaAdj,betaAdj,&
     MachAdj,machCoefAdj,machGridAdj,prefAdj,rhorefAdj, pinfdimAdj,&
     rhoinfdimAdj,rhoinfAdj, pinfAdj,rotRateAdj,rotCenterAdj,murefAdj,&
     timerefAdj,pInfCorrAdj,pointRefAdj,rotPointAdj,nn,level,sps,&
     liftIndex)

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Transfer the variable  x to the auxiliary stencil              *
  !     * xAdj is used by the Tapenade differentiated routines.          *
  !     *                                                                *
  !     * And compute boundary face normals (siAdj, sjAdj, skAdj)        *
  !     *                                                                *
  !     * It is assumed that the pointers in blockPointers have already  *
  !     * been set.                                                      *
  !     *                                                                *
  !     ******************************************************************
  !
  use blockPointers   ! w,il,jl,kl,ie,je,ke
  use communication   ! myID for debug
  use flowvarrefstate ! nw
  use inputPhysics    ! Mach,veldirfreestream
  use cgnsgrid        ! cgnsdoms
  use inputMotion     ! rotPoint
  implicit none
  !
  !     Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn,level,sps    
  real(kind=realType), dimension(0:ie,0:je,0:ke,3), intent(out) :: xAdj
  real(kind=realType), dimension(0:ib,0:jb,0:kb,1:nw), intent(out) :: wAdj

  real(kind=realType) :: alphaAdj, betaAdj,MachAdj,MachCoefAdj,MachGridAdj
  REAL(KIND=REALTYPE) :: prefAdj, rhorefAdj,pInfCorrAdj
  REAL(KIND=REALTYPE) :: pinfdimAdj, rhoinfdimAdj
  REAL(KIND=REALTYPE) :: rhoinfAdj, pinfAdj
  REAL(KIND=REALTYPE) :: murefAdj, timerefAdj
  integer(kind=intType)::liftIndex

  real(kind=realType), dimension(3),intent(out) ::rotRateAdj,rotCenterAdj
  real(kind=realType), dimension(3),intent(out) ::pointRefAdj,rotPointAdj
  !
  !     Local variables.
  !
  integer(kind=intType) :: ii, jj, kk, i, j, k, l, m, n

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !

  ! Initialize the auxiliar array wAdj.
  wAdj = zero
  wAdj(0:ib,0:jb,0:kb,1:nw) = w(0:ib,0:jb,0:kb,1:nw) 
  xAdj = zero
  xAdj(0:ie,0:je,0:ke,:) = x(0:ie,0:je,0:k,:)

  MachAdj = Mach
  MachCoefAdj = MachCoef
  MachGridAdj = MachGrid
  call getDirAngle(velDirFreestream,LiftDirection,liftIndex,alphaAdj,betaAdj)

  prefAdj = pRef
  rhorefAdj = rhoref
  pinfdimAdj = pinfdim
  rhoinfdimAdj = rhoinfdim
  rhoinfAdj = rhoinf
  pinfAdj = pInf
  murefAdj = muref
  timerefAdj = timeref
  pInfCorrAdj = pInfCorr

  ! Store the rotation center and determine the
  ! nonDimensional rotation rate of this block. As the
  ! reference length is 1 timeRef == 1/uRef and at the end
  ! the nonDimensional velocity is computed.

  j = nbkGlobal

  rotCenterAdj = cgnsDoms(j)%rotCenter
  rotRateAdj   = timeRef*cgnsDoms(j)%rotRate
  pointRefAdj = pointRef
  rotPointAdj = rotPoint
end subroutine copyADjointForcesStencil
