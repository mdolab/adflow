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
       use block
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

       integer(kind=intType) :: nCellsOverset, nCellsOversetAll
       integer(kind=intType) :: nHoles, nOrphans

       integer(kind=intType), dimension(:,:,:), pointer :: iblank

       integer(kind=intType), dimension(:,:), pointer :: ibndry
       integer(kind=intType), dimension(:,:), pointer :: idonor
       real(kind=realType),   dimension(:,:), pointer :: overint

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
       real(kind=realType), dimension(:,:,:,:), pointer :: volOld

       integer(kind=porType), dimension(:,:,:), pointer :: porI
       integer(kind=porType), dimension(:,:,:), pointer :: porJ
       integer(kind=porType), dimension(:,:,:), pointer :: porK

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

       real(kind=realType), dimension(:,:,:), pointer :: sFaceI
       real(kind=realType), dimension(:,:,:), pointer :: sFaceJ
       real(kind=realType), dimension(:,:,:), pointer :: sFaceK

       real(kind=realType), dimension(:,:,:,:),   pointer :: w
       real(kind=realType), dimension(:,:,:,:,:), pointer :: wOld
       real(kind=realType), dimension(:,:,:),     pointer :: p, gamma
       real(kind=realType), dimension(:,:,:),     pointer :: rlv, rev
       real(kind=realType), dimension(:,:,:,:),   pointer :: s

       real(kind=realType), dimension(:,:,:),     pointer :: p1
       real(kind=realType), dimension(:,:,:,:),   pointer :: dw, fw
       real(kind=realType), dimension(:,:,:,:,:), pointer :: dwOldRK
       real(kind=realType), dimension(:,:,:,:),   pointer :: w1, wr

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
       real(kind=realType), dimension(:,:,:),   pointer :: radI
       real(kind=realType), dimension(:,:,:),   pointer :: radJ
       real(kind=realType), dimension(:,:,:),   pointer :: radK

       real(kind=realType), dimension(:,:,:), pointer :: d2Wall

       real(kind=realType), dimension(:,:,:,:), pointer :: bmti1
       real(kind=realType), dimension(:,:,:,:), pointer :: bmti2
       real(kind=realType), dimension(:,:,:,:), pointer :: bmtj1
       real(kind=realType), dimension(:,:,:,:), pointer :: bmtj2
       real(kind=realType), dimension(:,:,:,:), pointer :: bmtk1
       real(kind=realType), dimension(:,:,:,:), pointer :: bmtk2

       real(kind=realType), dimension(:,:,:), pointer :: bvti1, bvti2
       real(kind=realType), dimension(:,:,:), pointer :: bvtj1, bvtj2
       real(kind=realType), dimension(:,:,:), pointer :: bvtk1, bvtk2

       end module blockPointers
