!
!      ******************************************************************
!      *                                                                *
!      * File:          setPointers.f90                                 *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 03-07-2003                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setPointers(nn,mm,ll)
!
!      ******************************************************************
!      *                                                                *
!      * setPointers makes the variables in blockPointers point to      *
!      * block nn for grid level mm and spectral solution ll.           *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: nn, mm, ll
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the info of the current block, such that inside the
       ! module blockPointers it is known to which block the data
       ! belongs.

       sectionID   = flowDoms(nn,mm,ll)%sectionID
       nbkLocal    = nn
       nbkGlobal   = flowDoms(nn,mm,ll)%cgnsBlockID
       mgLevel     = mm
       spectralSol = ll

       ! Block dimensions.

       nx = flowDoms(nn,mm,ll)%nx
       ny = flowDoms(nn,mm,ll)%ny
       nz = flowDoms(nn,mm,ll)%nz

       il = flowDoms(nn,mm,ll)%il
       jl = flowDoms(nn,mm,ll)%jl
       kl = flowDoms(nn,mm,ll)%kl

       ie = flowDoms(nn,mm,ll)%ie
       je = flowDoms(nn,mm,ll)%je
       ke = flowDoms(nn,mm,ll)%ke

       ib = flowDoms(nn,mm,ll)%ib
       jb = flowDoms(nn,mm,ll)%jb
       kb = flowDoms(nn,mm,ll)%kb

       ! Point range in the corresponding cgns block

       iBegor = flowDoms(nn,mm,ll)%iBegor
       iEndor = flowDoms(nn,mm,ll)%iEndor
       jBegor = flowDoms(nn,mm,ll)%jBegor
       jEndor = flowDoms(nn,mm,ll)%jEndor
       kBegor = flowDoms(nn,mm,ll)%kBegor
       kEndor = flowDoms(nn,mm,ll)%kEndor

       ! Subface info. Note that the pointers point to the 1st spectral
       ! mode, because this is the only one allocated. The info is the
       ! same for all modes.

       nSubface   = flowDoms(nn,mm,ll)%nSubface
       n1to1      = flowDoms(nn,mm,ll)%n1to1
       nBocos     = flowDoms(nn,mm,ll)%nBocos
       nViscBocos = flowDoms(nn,mm,ll)%nViscBocos

       BCType      => flowDoms(nn,mm,1)%BCType
       BCFaceID    => flowDoms(nn,mm,1)%BCFaceID
       cgnsSubface => flowDoms(nn,mm,1)%cgnsSubface

       inBeg => flowDoms(nn,mm,1)%inBeg
       jnBeg => flowDoms(nn,mm,1)%jnBeg
       knBeg => flowDoms(nn,mm,1)%knBeg
       inEnd => flowDoms(nn,mm,1)%inEnd
       jnEnd => flowDoms(nn,mm,1)%jnEnd
       knEnd => flowDoms(nn,mm,1)%knEnd

       dinBeg => flowDoms(nn,mm,1)%dinBeg
       djnBeg => flowDoms(nn,mm,1)%djnBeg
       dknBeg => flowDoms(nn,mm,1)%dknBeg
       dinEnd => flowDoms(nn,mm,1)%dinEnd
       djnEnd => flowDoms(nn,mm,1)%djnEnd
       dknEnd => flowDoms(nn,mm,1)%dknEnd

       icBeg => flowDoms(nn,mm,1)%icBeg
       jcBeg => flowDoms(nn,mm,1)%jcBeg
       kcBeg => flowDoms(nn,mm,1)%kcBeg
       icEnd => flowDoms(nn,mm,1)%icEnd
       jcEnd => flowDoms(nn,mm,1)%jcEnd
       kcEnd => flowDoms(nn,mm,1)%kcEnd

       neighBlock => flowDoms(nn,mm,1)%neighBlock
       neighProc  => flowDoms(nn,mm,1)%neighProc
       l1         => flowDoms(nn,mm,1)%l1
       l2         => flowDoms(nn,mm,1)%l2
       l3         => flowDoms(nn,mm,1)%l3
       groupNum   => flowDoms(nn,mm,1)%groupNum

       ! Overset boundary and hole info.

       nCellsOverset     = flowDoms(nn,mm,ll)%nCellsOverset
       nCellsOversetAll  = flowDoms(nn,mm,ll)%nCellsOversetAll
       nOrphans          = flowDoms(nn,mm,ll)%nOrphans
       nHoles            = flowDoms(nn,mm,ll)%nHoles

       ibndry           => flowDoms(nn,mm,ll)%ibndry
       idonor           => flowDoms(nn,mm,ll)%idonor
       overint          => flowDoms(nn,mm,ll)%overint
       neighProcOver    => flowDoms(nn,mm,ll)%neighProcOver
       neighBlockOver   => flowDoms(nn,mm,ll)%neighBlockOver

       iblank => flowDoms(nn,mm,ll)%iblank

       ! The data for boundary subfaces.

       BCData => flowDoms(nn,mm,ll)%BCData

       ! The stress tensor and heat flux vector at viscous wall faces
       ! as well as the face pointers to these viscous wall faces.
       ! The latter point to the 1st spectral mode, because they are
       ! the only ones allocated. The info is the same for all modes.

       viscSubface => flowDoms(nn,mm,ll)%viscSubface

       viscIminPointer => flowDoms(nn,mm,1)%viscIminPointer
       viscImaxPointer => flowDoms(nn,mm,1)%viscImaxPointer
       viscJminPointer => flowDoms(nn,mm,1)%viscJminPointer
       viscJmaxPointer => flowDoms(nn,mm,1)%viscJmaxPointer
       viscKminPointer => flowDoms(nn,mm,1)%viscKminPointer
       viscKmaxPointer => flowDoms(nn,mm,1)%viscKmaxPointer

       ! Mesh related variables. The porosities point to the 1st
       ! spectral mode, because they are the only ones allocated.
       ! The info is the same for all modes.
       ! Note that xOld and volOld always point to the finest
       ! grid level.

       x    => flowDoms(nn,mm,ll)%x
       xOld => flowDoms(nn,1,ll)%xOld

       si     => flowDoms(nn,mm,ll)%si
       sj     => flowDoms(nn,mm,ll)%sj
       sk     => flowDoms(nn,mm,ll)%sk
       vol    => flowDoms(nn,mm,ll)%vol
       volOld => flowDoms(nn,1,ll)%volOld

       porI => flowDoms(nn,mm,1)%porI
       porJ => flowDoms(nn,mm,1)%porJ
       porK => flowDoms(nn,mm,1)%porK

       indFamilyI => flowDoms(nn,mm,1)%indFamilyI
       indFamilyJ => flowDoms(nn,mm,1)%indFamilyJ
       indFamilyK => flowDoms(nn,mm,1)%indFamilyK

       factFamilyI => flowDoms(nn,mm,1)%factFamilyI
       factFamilyJ => flowDoms(nn,mm,1)%factFamilyJ
       factFamilyK => flowDoms(nn,mm,1)%factFamilyK

       rotMatrixI => flowDoms(nn,mm,ll)%rotMatrixI
       rotMatrixJ => flowDoms(nn,mm,ll)%rotMatrixJ
       rotMatrixK => flowDoms(nn,mm,ll)%rotMatrixK

       blockIsMoving     = flowDoms(nn,mm,ll)%blockIsMoving
       addGridVelocities = flowDoms(nn,mm,ll)%addGridVelocities

       sFaceI => flowDoms(nn,mm,ll)%sFaceI
       sFaceJ => flowDoms(nn,mm,ll)%sFaceJ
       sFaceK => flowDoms(nn,mm,ll)%sFaceK

       ! Flow variables. Note that wOld, gamma and the laminar viscosity
       ! point to the entries on the finest mesh. The reason is that
       ! they are computed from the other variables. For the eddy
       ! viscosity this is not the case because in a segregated solver
       ! its values are obtained from the fine grid level.

       w     => flowDoms(nn,mm,ll)%w
       wOld  => flowDoms(nn,1, ll)%wOld
       p     => flowDoms(nn,mm,ll)%p
       gamma => flowDoms(nn,1, ll)%gamma
       rlv   => flowDoms(nn,1, ll)%rlv
       rev   => flowDoms(nn,mm,ll)%rev
       s     => flowDoms(nn,mm,ll)%s

       ! Residual and multigrid variables. The residual point to the
       ! finest grid entry, the multigrid variables to their own level.

       dw => flowDoms(nn,1,ll)%dw
       fw => flowDoms(nn,1,ll)%fw

       dwOldRK => flowDoms(nn,1,ll)%dwOldRK

       p1 => flowDoms(nn,mm,ll)%p1
       w1 => flowDoms(nn,mm,ll)%w1
       wr => flowDoms(nn,mm,ll)%wr

       ! Variables, which allow a more flexible multigrid treatment.
       ! They are the same for all spectral modes and therefore they
       ! point to the 1st mode.

       mgIFine => flowDoms(nn,mm,1)%mgIFine
       mgJFine => flowDoms(nn,mm,1)%mgJFine
       mgKFine => flowDoms(nn,mm,1)%mgKFine

       mgIWeight => flowDoms(nn,mm,1)%mgIWeight
       mgJWeight => flowDoms(nn,mm,1)%mgJWeight
       mgKWeight => flowDoms(nn,mm,1)%mgKWeight

       mgICoarse => flowDoms(nn,mm,1)%mgICoarse
       mgJCoarse => flowDoms(nn,mm,1)%mgJCoarse
       mgKCoarse => flowDoms(nn,mm,1)%mgKCoarse

       ! Time-stepping variables and spectral radIi.
       ! They all point to the fine mesh entry.

       wn  => flowDoms(nn,1,ll)%wn
       pn  => flowDoms(nn,1,ll)%pn
       dtl => flowDoms(nn,1,ll)%dtl

       radI => flowDoms(nn,1,ll)%radI
       radJ => flowDoms(nn,1,ll)%radJ
       radK => flowDoms(nn,1,ll)%radK

       ! Wall distance for the turbulence models.

       d2Wall => flowDoms(nn,mm,ll)%d2Wall

       ! Arrays used for the implicit treatment of the turbulent wall
       ! boundary conditions. As these variables are only allocated for
       ! the 1st spectral solution of the fine mesh, the pointers point
       ! to those arrays.

       bmti1 => flowDoms(nn,1,1)%bmti1
       bmti2 => flowDoms(nn,1,1)%bmti2
       bmtj1 => flowDoms(nn,1,1)%bmtj1
       bmtj2 => flowDoms(nn,1,1)%bmtj2
       bmtk1 => flowDoms(nn,1,1)%bmtk1
       bmtk2 => flowDoms(nn,1,1)%bmtk2

       bvti1 => flowDoms(nn,1,1)%bvti1
       bvti2 => flowDoms(nn,1,1)%bvti2
       bvtj1 => flowDoms(nn,1,1)%bvtj1
       bvtj2 => flowDoms(nn,1,1)%bvtj2
       bvtk1 => flowDoms(nn,1,1)%bvtk1
       bvtk2 => flowDoms(nn,1,1)%bvtk2

       end subroutine setPointers
