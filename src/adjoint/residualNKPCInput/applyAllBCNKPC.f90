!
!      ******************************************************************
!      *                                                                *
!      * File:          applyAllBCAdj.f90                               *
!      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
!      *                C.A.(Sandy) Mader                               *
!      * Starting date: 04-16-2008                                      *
!      * Last modified: 04-17-2008                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine applyAllBCNKPC(wInfAdj,pInfCorrAdj,wAdj, pAdj,sAdj, &
     siAdj, sjAdj, skAdj, volAdj, normAdj, rFaceAdj,&
     iCell, jCell, kCell,secondHalo,nn,level,sps,sps2)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * applyAllBCAdj applies the possible boundary conditions for the *
  !      * halo cells adjacent to the cell for which the residual needs   *
  !      * to be computed.                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  use BCTypes
  use blockPointers!, only : ie, ib, je, jb, ke, kb, nBocos, &
  !         BCFaceID, BCType, BCData,p,w
  use flowVarRefState
  use inputDiscretization !precond,choimerkle, etc...
  use inputTimeSpectral !nIntervalTimespectral
  implicit none
  !
  !      Subroutine arguments.
  !
  logical:: secondHalo

  integer(kind=intType) :: iCell, jCell, kCell
  integer(kind=intType) :: nn,level, sps,sps2

  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral), &
       intent(inout) :: wAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral),    &
       intent(inout) :: pAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,3,nTimeIntervalsSpectral),intent(in) :: sAdj
  real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral), &
       intent(in) :: siAdj, sjAdj, skAdj

  real(kind=realType),dimension(nTimeIntervalsSpectral), intent(in) :: volAdj
  real(kind=realType), intent(in) ::pInfCorrAdj
  real(kind=realType), dimension(nBocos,-2:2,-2:2,3,nTimeIntervalsSpectral), intent(in) :: normAdj
  real(kind=realType), dimension(nBocos,-2:2,-2:2,nTimeIntervalsSpectral), intent(in) ::rFaceAdj
  real(kind=realType), dimension(nw),intent(in)::wInfAdj

  !
  !      Local variables.
  !

  integer(kind=intType)::i,j,k,ii,jj,kk,l
  integer(kind=intType) :: iStart,iEnd,jStart,jEnd,kStart,kEnd

  logical :: correctForK
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !

  call bcSymmNKPC(wAdj,pAdj,normAdj,iCell,jCell,kCell,secondHalo,nn,level,sps,sps2)

  select case (precond)

  case (noPrecond)
     call bcFarfieldNKPC(secondHalo,wInfAdj,pInfCorrAdj, wAdj,pAdj,      &
          siAdj, sjAdj, skAdj, normAdj,rFaceAdj,iCell,jCell,kCell,nn,level,sps,sps2)
  case (Turkel)
     call terminate("applyAllBC", "Farfield boundary conditions for Turkel preconditioner not implemented")
  case (ChoiMerkle)
     call terminate("applyAllBC", "Farfield boundary conditions for Choi and Merkle preconditioner not implemented")
  end select

  call bcEulerWallNKPC(secondHalo, wAdj,pAdj,sAdj,      &
       siAdj, sjAdj, skAdj, normAdj,rFaceAdj,iCell,jCell,kCell,nn,level,sps,sps2)

end subroutine applyAllBCNKPC
