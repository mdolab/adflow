!
!      ******************************************************************
!      *                                                                *
!      * File:          blockPointers.f90                               *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 03-07-2003                                      *
!      * Last modified: 11-21-2007                                      *
!      *                                                                *
!      ******************************************************************
!
module blockPointers
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This module contains the pointers for all variables inside a   *
  !      * block. The pointers are set via the subroutine setPointers,    *
  !      * which can be found in the utils directory. In this way the     *
  !      * code becomes much more readable. The relation to the original  *
  !      * multiblock grid is not copied, because it does not affect the  *
  !      * computation.                                                   *
  !      *                                                                *
  !      * See the module block for the meaning of the variables.         *
  !      *                                                                *
  !      * Note that the dimensions are not pointers, but integers.       *
  !      * Consequently changing dimensions of a block must be done only  *
  !      * with the variables of floDoms.                                 *
  !      *                                                                *
  !      ******************************************************************
  !
  use constants, only : intType, realType, porType
  use block, only : fringeType, BCDataType, viscSubFaceType, flowDoms, nDom
#ifndef USE_TAPENADE
  use block, only : flowDomsd
#endif
  implicit none
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Additional info, such that it is known to which block the data *
  !      * inside this module belongs.                                    *
  !      *                                                                *
  !      ******************************************************************
  !
  ! sectionID:   the section to which this block belongs.
  ! nbkLocal :   local block number.
  ! nbkGlobal:   global block number in the original cgns grid.
  ! mgLevel:     the multigrid level.
  ! spectralSol: the spectral solution index of this block.

  integer(kind=intType) :: sectionID
  integer(kind=intType) :: nbkLocal, nbkGlobal, mgLevel
  integer(kind=intType) :: spectralSol
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Variables, which are either copied or the pointer is set to    *
  !      * the correct variable in the block. See the module block for    *
  !      * meaning of the variables.                                      *
  !      *                                                                *
  !      ******************************************************************
  !
  integer(kind=intType) :: nx, ny, nz, il, jl, kl
  integer(kind=intType) :: ie, je, ke, ib, jb, kb
  integer(kind=intType) :: maxDim, imaxDim, jmaxDim

  logical :: rightHanded

  integer(kind=intType) :: iBegOr, iEndOr, jBegOr, jEndOr
  integer(kind=intType) :: kBegOr, kEndOr

  integer(kind=intType) :: nSubface, n1to1, nBocos, nViscBocos

  integer(kind=intType), dimension(:), pointer :: BCType
  integer(kind=intType), dimension(:), pointer :: BCFaceID

  integer(kind=intType), dimension(:), pointer :: cgnsSubface

  integer(kind=intType), dimension(:), pointer :: inBeg, inEnd
  integer(kind=intType), dimension(:), pointer :: jnBeg, jnEnd
  integer(kind=intType), dimension(:), pointer :: knBeg, knEnd

  integer(kind=intType), dimension(:), pointer :: dinBeg, dinEnd
  integer(kind=intType), dimension(:), pointer :: djnBeg, djnEnd
  integer(kind=intType), dimension(:), pointer :: dknBeg, dknEnd

  integer(kind=intType), dimension(:), pointer :: icBeg, icEnd
  integer(kind=intType), dimension(:), pointer :: jcBeg, jcEnd
  integer(kind=intType), dimension(:), pointer :: kcBeg, kcEnd

  integer(kind=intType), dimension(:), pointer :: neighBlock
  integer(kind=intType), dimension(:), pointer :: neighProc
  integer(kind=intType), dimension(:), pointer :: l1, l2, l3
  integer(kind=intType), dimension(:), pointer :: groupNum

  integer(kind=intType), dimension(:,:,:), pointer :: iblank
  type(fringeType), dimension(:, :, :), pointer :: fringes
  integer(kind=intType), dimension(:, :), pointer :: orphans
  integer(kind=intType) :: nOrphans

  integer(kind=intType), dimension(:), pointer :: neighBlockOver
  integer(kind=intType), dimension(:), pointer :: neighProcOver

  type(BCDataType),      dimension(:), pointer :: BCData
  type(viscSubfaceType), dimension(:), pointer :: viscSubface

  integer(kind=intType), dimension(:,:), pointer :: viscIMinPointer
  integer(kind=intType), dimension(:,:), pointer :: viscIMaxPointer
  integer(kind=intType), dimension(:,:), pointer :: viscJMinPointer
  integer(kind=intType), dimension(:,:), pointer :: viscJMaxPointer
  integer(kind=intType), dimension(:,:), pointer :: viscKMinPointer
  integer(kind=intType), dimension(:,:), pointer :: viscKMaxPointer

  real(kind=realType), dimension(:,:,:,:),   pointer :: x
  real(kind=realType), dimension(:,:,:,:,:), pointer :: xOld
  real(kind=realType), dimension(:,:,:,:), pointer :: sI, sJ, sK
  real(kind=realType), dimension(:,:,:),   pointer :: vol
  real(kind=realType), dimension(:,:,:),   pointer :: volref
  real(kind=realType), dimension(:,:,:,:), pointer :: volOld
  real(kind=realType), dimension(:,:,:,:),   pointer :: dadidata

  integer(kind=porType), dimension(:,:,:), pointer :: porI, porJ, porK

  integer(kind=intType), dimension(:,:,:), pointer :: indFamilyI
  integer(kind=intType), dimension(:,:,:), pointer :: indFamilyJ
  integer(kind=intType), dimension(:,:,:), pointer :: indFamilyK

  integer(kind=intType), dimension(:,:,:), pointer :: factFamilyI
  integer(kind=intType), dimension(:,:,:), pointer :: factFamilyJ
  integer(kind=intType), dimension(:,:,:), pointer :: factFamilyK

  real(kind=realType), dimension(:,:,:,:,:), pointer :: rotMatrixI
  real(kind=realType), dimension(:,:,:,:,:), pointer :: rotMatrixJ
  real(kind=realType), dimension(:,:,:,:,:), pointer :: rotMatrixK

  logical :: blockIsMoving, addGridVelocities

  real(kind=realType), dimension(:,:,:), pointer :: sFaceI, sFaceJ, sfaceK
  real(kind=realType), dimension(:,:,:,:),   pointer :: w
  real(kind=realType), dimension(:,:,:,:,:), pointer :: wOld

  real(kind=realType), dimension(:,:,:),     pointer :: p, gamma, aa
  real(kind=realType), dimension(:,:,:),     pointer :: shockSensor
  real(kind=realType), dimension(:,:,:),     pointer :: rlv, rev
  real(kind=realType), dimension(:,:,:,:),   pointer :: s
  real(kind=realType), dimension(:,:,:),     pointer :: p1
  real(kind=realType), dimension(:,:,:,:),   pointer :: dw, fw
  real(kind=realType), dimension(:,:,:,:),   pointer :: scratch
  real(kind=realType), dimension(:,:,:,:,:), pointer :: dwOldRK
  real(kind=realType), dimension(:,:,:,:),   pointer :: w1, wr
  real(kind=realType), dimension(:, :, :), pointer:: ux, uy, uz
  real(kind=realType), dimension(:, :, :), pointer:: vx, vy, vz 
  real(kind=realType), dimension(:, :, :), pointer:: wx, wy, wz 
  real(kind=realType), dimension(:, :, :), pointer:: qx, qy, qz 
  
  integer(kind=intType), dimension(:,:), pointer :: mgIFine
  integer(kind=intType), dimension(:,:), pointer :: mgJFine
  integer(kind=intType), dimension(:,:), pointer :: mgKFine

  real(kind=realType),   dimension(:),   pointer :: mgIWeight
  real(kind=realType),   dimension(:),   pointer :: mgJWeight
  real(kind=realType),   dimension(:),   pointer :: mgKWeight

  integer(kind=intType), dimension(:,:), pointer :: mgICoarse
  integer(kind=intType), dimension(:,:), pointer :: mgJCoarse
  integer(kind=intType), dimension(:,:), pointer :: mgKCoarse

  real(kind=realType), dimension(:,:,:,:), pointer :: wn
  real(kind=realType), dimension(:,:,:),   pointer :: pn
  real(kind=realType), dimension(:,:,:),   pointer :: dtl
  real(kind=realType), dimension(:,:,:),   pointer :: radI, radJ, radK

  real(kind=realType), dimension(:,:,:), pointer :: d2Wall
  real(kind=realType), dimension(:,:,:),   pointer :: filterDES  ! eran-des
  real(kind=realType), dimension(:,:,:,:), pointer :: bmti1
  real(kind=realType), dimension(:,:,:,:), pointer :: bmti2
  real(kind=realType), dimension(:,:,:,:), pointer :: bmtj1
  real(kind=realType), dimension(:,:,:,:), pointer :: bmtj2
  real(kind=realType), dimension(:,:,:,:), pointer :: bmtk1
  real(kind=realType), dimension(:,:,:,:), pointer :: bmtk2
  real(kind=realType), dimension(:,:,:), pointer :: bvti1, bvti2
  real(kind=realType), dimension(:,:,:), pointer :: bvtj1, bvtj2
  real(kind=realType), dimension(:,:,:), pointer :: bvtk1, bvtk2

  integer(kind=intType), dimension(:,:,:), pointer :: globalNode
  integer(kind=intType), dimension(:,:,:), pointer :: globalCell
  real(kind=realType), dimension(:, :, :, :), pointer :: xSeed
  integer(kind=intType), dimension(:, :, :), pointer :: wallInd

  integer(kind=intType), dimension(:), pointer ::ifaceptb,iedgeptb

  real(kind=realType), dimension(:,:,:,:), pointer :: w_offTimeInstance
  real(kind=realType), dimension(:,:,:), pointer :: vol_offTimeInstance

  ! *******************************
  ! Added by HDN
  ! *******************************
  real(kind=realType), dimension(:,:,:,:),   pointer :: xALE
  real(kind=realType), dimension(:,:,:,:),   pointer :: sVeloIALE, sVeloJALE, sVeloKALE
  real(kind=realType), dimension(:,:,:,:,:), pointer :: sIALE, sJALE, sKALE
  real(kind=realType), dimension(:,:,:,:),   pointer :: sFaceIALE, sFaceJALE, sFaceKALE
  real(kind=realType), dimension(:,:,:,:,:), pointer :: dwALE, fwALE
  

#ifndef USE_TAPENADE
  TYPE(VISCSUBFACETYPE), DIMENSION(:), POINTER :: viscsubfaced

  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: xd
  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: sid, sjd, skd

  real(kind=realType), dimension(:,:,:),   pointer ::vold
  
  REAL(kind=realtype), DIMENSION(:, :, :, :, :), POINTER :: rotmatrixid
  REAL(kind=realtype), DIMENSION(:, :, :, :, :), POINTER :: rotmatrixjd
  REAL(kind=realtype), DIMENSION(:, :, :, :, :), POINTER :: rotmatrixkd

  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: sfaceid
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: sfacejd
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: sfacekd

  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: wd
  REAL(kind=realtype), DIMENSION(:, :, :, :, :), POINTER :: woldd

  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: uxd
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: uyd
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: uzd

  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: vxd
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: vyd
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: vzd

  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: wxd
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: wyd
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: wzd

  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: qxd
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: qyd
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: qzd

  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: pd, gammad, aad
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: rlvd, revd

  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: sd

  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: dwd, fwd
  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: w1d, wrd
  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: scratchd

  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: dtld
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: radid
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: radjd
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: radkd

  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmti1d
  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmti2d
  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmtj1d
  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmtj2d
  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmtk1d
  REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmtk2d

  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: bvti1d, bvti2d
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: bvtj1d, bvtj2d
  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: bvtk1d, bvtk2d

  REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: d2walld

  real(kind=realType), dimension(:,:,:,:), pointer :: w_offTimeInstanced
  real(kind=realType), dimension(:,:,:), pointer :: vol_offTimeInstanced
  
  type(BCDataType),      dimension(:), pointer :: BCDatad

  real(kind=realType), dimension(:, :, :, :, :), pointer :: PCMat
    real(kind=realType), dimension(:, :, :, :), pointer :: PCvec1, PCvec2  

  real(kind=realType), dimension(:, :, :, :), pointer :: i_D_fact,  j_D_fact,  k_D_fact 
  real(kind=realType), dimension(:, :, :, :), pointer :: i_L_Fact,  j_L_Fact,  k_L_Fact 
  real(kind=realType), dimension(:, :, :, :), pointer :: i_U_Fact,  j_U_Fact,  k_U_Fact 
  real(kind=realType), dimension(:, :, :, :), pointer :: i_U2_Fact, j_U2_Fact, k_U2_Fact
  integer(kind=intType), dimension(:, :, :, :), pointer :: i_ipiv, j_ipiv, k_ipiv
#endif

end module blockPointers
