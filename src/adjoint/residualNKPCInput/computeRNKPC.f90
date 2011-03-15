!
!      ******************************************************************
!      *                                                                *
!      * File:          computeRAdj.f90                                 *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 02-01-2008                                      *
!      * Last modified: 04-23-2008                                      *
!      *                                                                *
!      ******************************************************************

subroutine computeRNKPC(wAdj, dwAdj, siAdj, sjAdj, skAdj, sAdj, volAdj, &
     sfaceIAdj, sfaceJAdj, sfaceKadj, rotRateAdj, iCell, jCell, kCell, &
     nn, level, sps)

  use blockPointers
  use flowVarRefState
  use inputTimeSpectral
  use inputPhysics
  use section 
  implicit none

  ! Passed in Variables

  ! Input Variables --- Note no intent(in) --- this can cause problems
  ! with reverse mode AD
  integer(kind=intType), intent(in) :: iCell, jCell, kCell,nn,level,sps
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral) &
       :: wAdj
  real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral) :: &
       siAdj, sjAdj, skAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,3, &
       nTimeIntervalsSpectral) :: sAdj
  real(kind=realType),dimension(nTimeIntervalsSpectral):: volAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral) :: &
       sFaceIAdj,sFaceJAdj,sFaceKAdj
  real(kind=realType), dimension(3) :: rotRateAdj 
  ! Ouptut Variables
  real(kind=realType), dimension(nw,nTimeIntervalsSpectral)  :: dwAdj

  !      Set Local Variables
 
  real(kind=realType), dimension(-2:2,-2:2,-2:2,&
       nTimeIntervalsSpectral) :: pAdj
   real(kind=realType), dimension(nBocos,-2:2,-2:2,3,&
       nTimeIntervalsSpectral) :: normAdj
  real(kind=realType), dimension(nBocos,-2:2,-2:2,&
       nTimeIntervalsSpectral) ::rFaceAdj

  real(kind=realType), dimension(-1:1,-1:1,-1:1,&
       nTimeIntervalsSpectral) :: radIAdj,radJAdj,radKAdj

  real(kind=realType), dimension(nSections) :: t

  ! Integer variables
  integer(kind=intType) :: sps2
 
  ! Logical 
  logical :: secondHalo,correctForK
  secondHalo = .True. 
  correctForK = .false.
 
  do sps2 = 1,nTimeIntervalsSpectral
     call computeNormNKPC(siAdj,sjAdj,skAdj,normAdj,iCell,jCell,kCell,&
          nn,level,sps,sps2)

     !Compute the Pressure in the stencil based on the current 
     !States
  
     call normalVelocitiesAllLevelsNKPC(sps,iCell, jCell, kCell,sFaceIAdj,&
          sFaceJAdj,sFaceKAdj,siAdj, sjAdj, skAdj,rFaceAdj,nn,level,sps2)

     call computePressureNKPC(wAdj, pAdj,nn,level,sps,sps2)  

     ! Apply all boundary conditions to stencil.
     ! In case of a full mg mode, and a segegated turbulent solver,
     ! first call the turbulent boundary conditions, such that the
     ! turbulent kinetic energy is properly initialized in the halo's.

     call applyAllBCNKPC(wInf,pInfCorr,wAdj, pAdj,sAdj, &
          siAdj, sjAdj, skAdj, volAdj, normAdj, &
          rFaceAdj,iCell, jCell, kCell,secondHalo,nn,level,sps,sps2)
     call timeStepNKPC(.true.,wAdj,pAdj,siAdj, sjAdj, skAdj,&
          sFaceIAdj,sFaceJAdj,sFaceKAdj,volAdj,radIAdj,radJAdj,radKAdj,&
          iCell, jCell, kCell,pInfCorr,rhoInf,nn,level,sps,sps2)
  end do

  call initresNKPC(1, nwf,wAdj,volAdj,dwAdj,nn,level,sps)
  call residualNKPC(wAdj,pAdj,siAdj,sjAdj,skAdj,volAdj,normAdj,&
       sFaceIAdj,sFaceJAdj,sFaceKAdj,radIAdj,radJAdj,radKAdj,&
       dwAdj, iCell, jCell, kCell,rotRateAdj,correctForK,nn,level,sps)

  ! Finally scale the residual by the inverse of the volume.
  do sps2=1,ntimeintervalsspectral
     dwAdj(:,sps2) =dwAdj(:,sps2)/voladj(sps2)
  end do
  
end subroutine computeRNKPC
