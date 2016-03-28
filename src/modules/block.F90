!
!      ******************************************************************
!      *                                                                *
!      * File:          block.f90                                       *
!      * Author:        Edwin van der Weide, Steve Repsher,             *
!      *                Seonghyeon Hahn                                 *
!      * Starting date: 12-19-2002                                      *
!      * Last modified: 11-21-2007                                      *
!      *                                                                *
!      ******************************************************************
!
module block
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This module contains the definition of the derived data type   *
  !      * for block, which is the basic building block for this code.    *
  !      *                                                                *
  !      * Apart from the derived data type for block, this module also   *
  !      * contains the actual array for storing the blocks and the       *
  !      * number of blocks stored on this processor.                     *
  !      *                                                                *
  !      ******************************************************************
  !
  use constants
  implicit none
  save

  ! Parameters used for coarsening definition.

  integer(kind=porType), parameter :: leftStarted  = -1_porType
  integer(kind=porType), parameter :: regular      =  0_porType
  integer(kind=porType), parameter :: rightStarted =  1_porType

  ! Parameters used for subsonic inlet bc treatment.

  integer(kind=intType), parameter :: noSubInlet      = 0
  integer(kind=intType), parameter :: totalConditions = 1
  integer(kind=intType), parameter :: massFlow        = 2
  !
  !      ******************************************************************
  !      *                                                                *
  !      * The definition of the derived data type visc_subface_type,     *
  !      * which stores the viscous stress tensor and heat flux vector.   *
  !      * In this way it is avoided that these quantities must be        *
  !      * recomputed for the viscous forces and postprocessing. This     *
  !      * saves both time and a considerable amount of code.             *
  !      *                                                                *
  !      ******************************************************************
  !
  type viscSubfaceType

     ! tau(:,:,6): The 6 components of the viscous stress tensor.
     !             The first 2 dimensions of these arrays are equal
     !             to the dimenions of the cell subface without any
     !             halo cell. Consequently the starting index is
     !             arbitrary, such that no offset computation is
     !             needed when the arrays are accessed.
     ! q(:,:,3):   Same story for the heat flux vector.
     ! uTau(:,:):  And for the friction velocity.

     real(kind=realType), dimension(:,:,:), pointer :: tau, q
     real(kind=realType), dimension(:,:),   pointer :: uTau

  end type viscSubfaceType
  !
  !      ******************************************************************
  !      *                                                                *
  !      * The definition of the derived data type BCDataType, which      *
  !      * stores the prescribed data of boundary faces as well as unit   *
  !      * normals. For all the arrays the first two dimensions equal the *
  !      * dimensions of the subface, possibly extended with halo cells.  *
  !      * Consequently the starting index is arbitrary, such that no     *
  !      * offset computation is needed when the array is accessed.       *
  !      *                                                                *
  !      ******************************************************************
  !
  type BCDataType

     ! inBeg, inEnd: Node range in the first direction of the subface
     ! jnBeg, jnEnd: Idem in the second direction.
     ! icBeg, icEnd: Cell range in the first direction of the subface
     ! jcBeg, jcEnd: Idem in the second direction.

     integer(kind=intType) :: inBeg, inEnd, jnBeg, jnEnd
     integer(kind=intType) :: icBeg, icEnd, jcBeg, jcEnd

     ! norm(:,:,3):  The unit normal; it points out of the domain.
     ! rface(:,:):   Velocity of the face in the direction of the
     !               outward pointing normal. only allocated for
     !               the boundary conditions that need this info.

     real(kind=realType), dimension(:,:,:), pointer :: norm
     real(kind=realType), dimension(:,:),   pointer :: rface
     real(kind=realType), dimension(:,:,:), pointer :: F, Fv, Fp
     real(kind=realType), dimension(:,:), pointer :: dualArea

     ! symNorm is the normal for (symmertry) boundary conditions.
     ! symNormSet is set to false until symNorm is computed at the
     ! beginning of a simulation. symNorm then remains constant for
     ! the remainder of the simulation. This is ok, since if the
     ! normal of the symmetry plane is changing, your results are
     ! invalid anyway.  These values are only used on symmetry
     ! plane. They are undefined for other BC's.
     real(kind=realType), dimension(3) :: symNorm

     logical :: symNormSet 

     ! subsonicInletTreatment: which boundary condition treatment
     !                         to use for subsonic inlets; either
     !                         totalConditions or massFlow.

     integer(kind=intType) :: subsonicInletTreatment

     ! uSlip(:,:,3):  the 3 components of the velocity vector on
     !                a viscous wall. 
     ! TNS_Wall(:,:): Wall temperature for isothermal walls.

     real(kind=realType), dimension(:,:,:), pointer :: uSlip
     real(kind=realType), dimension(:,:),   pointer :: TNS_Wall
     integer(kind=intType) :: mask


     ! *******************************
     ! Added by HDN
     ! *******************************
     ! normALE(0:nALEsteps,ie:ib,je:jb,3)
     !                 - Storage of norm for intermediate meshes.
     ! rFaceALE(0:nALEsteps,iBeg:iEnd,jBeg:jEnd)
     !                 - Storage of rface for intermediate meshes.
     ! uSlipALE(0:nALEsteps,iBeg:iEnd,jBeg:jEnd,3)
     !                 - Storage of uSlip for intermediate meshes.
     ! sHeatFlux(0:nALEsteps,iBeg:iEnd,jBeg:jEnd,3)
     !                 - Surface heat flux.
     real(kind=realType), dimension(:,:,:,:), pointer :: normALE
     real(kind=realType), dimension(:,:,:),   pointer :: rFaceALE
     real(kind=realType), dimension(:,:,:,:), pointer :: uSlipALE
     real(kind=realType), dimension(:,:),     pointer :: sHeatFlux


     ! ptInlet(:,:):       Total pressure at subsonic inlets.
     ! ttInlet(:,:):       Total temperature at subsonic inlets.
     ! htInlet(:,:):       Total enthalpy at subsonic inlets.
     ! flowXDirInlet(:,:): X-direction of the flow for subsonic
     !                     inlets.
     ! flowYDirInlet(:,:): Idem in y-direction.
     ! flowZDirInlet(:,:): Idem in z-direction.

     real(kind=realType), dimension(:,:), pointer :: ptInlet, ttInlet, htInlet
     real(kind=realType), dimension(:,:), pointer :: flowXDirInlet
     real(kind=realType), dimension(:,:), pointer :: flowYDirInlet
     real(kind=realType), dimension(:,:), pointer :: flowZDirInlet

     ! turbInlet(:,:,nt1:nt2): Turbulence variables at inlets,
     !                         either subsonic or supersonic.

     real(kind=realType), dimension(:,:,:), pointer :: turbInlet

     ! rho(:,:):  density; used for multiple bc's.
     ! velX(:,:): x-velocity; used for multiple bc's.
     ! velY(:,:): y-velocity; used for multiple bc's.
     ! velZ(:,:): z-velocity; used for multiple bc's.
     ! ps(:,:):   static pressure; used for multiple bc's.

     real(kind=realType), dimension(:,:), pointer :: rho
     real(kind=realType), dimension(:,:), pointer :: velX
     real(kind=realType), dimension(:,:), pointer :: velY
     real(kind=realType), dimension(:,:), pointer :: velZ
     real(kind=realType), dimension(:,:), pointer :: ps

     ! Surface blanking for force integration
     integer(kind=intType), dimension(:,:), pointer :: iblank
      
  end type BCDataType

  type surfaceNodeWeightArray
     real(kind=realType), dimension(:, :, :), pointer :: weight
  end type surfaceNodeWeightArray

 type fringeType

     ! Make everything in here static such that we can easily copy the
     ! datatype and use MPI to communicate them directly. This data
     ! type bears some resemblence to the haloList type used for the
     ! B2B preprocessing.

    ! qualaity is the best quality that has been found from a
    ! DONOR cell. It is initialized to a 'large' value.
    real(kind=realType) :: quality

    ! This is the information regarding where the cell came from. 
    integer(kind=intType) :: myBlock, myI, myJ, myK

    ! This is the information about the donor that was found. Note we
    ! use dI, dJ, dK, short for donorI, etc.
    integer(kind=intType) :: donorProc, donorBlock, dI, dJ, dK
    real(kind=realType) :: donorFrac(3)

    ! gInd are the global indices of the donor cells. We will need
    ! these for forming the PC for the Newton Krylov solver
    integer(kind=intType), dimension(8) :: gInd
    
    ! Status integer: This stores the following status information:
    ! isDonor, isHole, isCompute, isFloodSeed, isFlooded, isWall and
    ! is wallDonor. 

    integer(kind=intType) :: status

  end type fringeType
  
  interface operator(<=)
     module procedure lessEqualFringeType
  end interface operator(<=)

  interface operator(<)
     module procedure lessFringeType
  end interface operator(<)


  !      ******************************************************************
  !      *                                                                *
  !      * The definition of the derived data type block_type, which      *
  !      * stores dimensions, coordinates, solution, etc.                 *
  !      *                                                                *
  !      ******************************************************************
  !
  type blockType
     !
     !        ****************************************************************
     !        *                                                              *
     !        * Block dimensions and orientation.                            *
     !        *                                                              *
     !        ****************************************************************
     !
     !  nx, ny, nz - Block integer dimensions for no halo cell based
     !               quantities.
     !  il, jl, kl - Block integer dimensions for no halo node based
     !               quantities.
     !  ie, je, ke - Block integer dimensions for single halo
     !               cell-centered quantities.
     !  ib, jb, kb - Block integer dimensions for double halo
     !               cell-centered quantities.
     ! rightHanded - Whether or not the block is a right handed.
     !               If not right handed it is left handed of course.

     integer(kind=intType) :: nx, ny, nz
     integer(kind=intType) :: il, jl, kl
     integer(kind=intType) :: ie, je, ke
     integer(kind=intType) :: ib, jb, kb

     logical :: rightHanded
     !
     !        ****************************************************************
     !        *                                                              *
     !        * Block boundary conditions.                                   *
     !        *                                                              *
     !        ****************************************************************
     !
     !  nSubface             - Number of subfaces on this block.
     !  n1to1                - Number of 1 to 1 block boundaries.
     !  nBocos               - Number of physical boundary subfaces.
     !  nViscBocos           - Number of viscous boundary subfaces.
     !  BCType(:)            - Boundary condition type for each
     !                         subface. See the module BCTypes for
     !                         the possibilities.
     !  BCFaceID(:)          - Block face location of each subface.
     !                         possible values are: iMin, iMax, jMin,
     !                         jMax, kMin, kMax. see also module
     !                         BCTypes.
     !  cgnsSubface(:)       - The subface in the corresponding cgns
     !                         block. As cgns distinguishes between
     !                         boundary and internal boundaries, the
     !                         BCType of the subface is needed to
     !                         know which one to take.
     !  inBeg(:), inEnd(:)   - Lower and upper limits for the nodes
     !  jnBeg(:), jnEnd(:)     in each of the index directions on a
     !  knBeg(:), knEnd(:)     given subface. Note that one of these
     !                         indices will not change since we will
     !                         be moving on a face.
     !  dinBeg(:), dinEnd(:) - Lower and upper limits for the nodes
     !  djnBeg(:), djnEnd(:)   in the each of the index directions
     !  dknBeg(:), dknEnd(:)   of the donor subface for this
     !                         particular subface. Note that one of
     !                         these indices will not change since we
     !                         will be moving on a face.
     !  icBeg(:), icEnd(:)   - Lower and upper limits for the cells
     !  jcBeg(:), jcEnd(:)     in each of the index directions for
     !  kcBeg(:), kcEnd(:)     the subface. The cells indicated by
     !                         this range are halo cells (the 
     !                         constant index) adjacent to the face.
     !                         a possible overlap outside the block
     !                         is stored.
     !  neighBlock(:)        - Local block number to which this
     !                         subface connects. This value is set to
     !                         zero if this subface is not connected
     !                         to another block.
     !  neighProc(:)         - Processor number where the neighbor
     !                         block is stored. This value is set to
     !                         -1 if this subface is not connected
     !                         to another block.
     !  l1(:), l2(:),        - Short hand for the transformation
     !  l3(:)                  matrix between this subface and the
     !                         neighbor block. These values are set
     !                         to zero if this subface is not
     !                         connected to another block.
     !  groupNum(:)          - Group number to which this subface
     !                         belongs. If this subface does not
     !                         belong to any group, the corresponding
     !                         entry in this array is zeroed out. If
     !                         the subface belongs to a sliding mesh
     !                         interface the absolute value of 
     !                         groupNum contains the number of the
     !                         sliding mesh interface. One side of
     !                         the interface gets a positive number,
     !                         the other side a negative one.
     !
     !
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

     !
     !        ****************************************************************
     !        *                                                              *
     !        * Overset interpolation information                            *
     !        *                                                              *
     !        ****************************************************************

     integer(kind=intType), dimension(:,:,:), pointer :: iblank
     type(fringeType) , dimension(:, :, :), pointer :: fringes
     integer(kind=intType), dimension(:, :), pointer :: orphans
     integer(kind=intType) :: nOrphans

     !
     !        ****************************************************************
     !        *                                                              *
     !        * Boundary data for the boundary subfaces.                     *
     !        *                                                              *
     !        ****************************************************************
     !
     ! BCData(nBocos): The boundary data for each of the boundary
     !                 subfaces.

     type(BCDataType), dimension(:), pointer :: BCData
     !
     !        ****************************************************************
     !        *                                                              *
     !        * The stress tensor and heat flux vector at viscous wall faces *
     !        * as well as the face pointers to these viscous wall faces.    *
     !        *                                                              *
     !        ****************************************************************
     !
     ! viscSubface(nViscBocos):    Storage for the viscous stress
     !                             tensor and heat flux vector for
     !                             the viscous subfaces.
     ! viscIMinPointer(2:jl,2:kl): Pointer to viscous subface for 
     !                             the iMin block face. If the face
     !                             is not part of a viscous subface
     !                             this value is set to 0.
     ! viscIMaxPointer(2:jl,2:kl): Idem for iMax block face.
     ! viscJMinPointer(2:il,2:kl): Idem for jMin block face.
     ! viscJMaxPointer(2:il,2:kl): Idem for jmax block face.
     ! viscKMinPointer(2:il,2:jl): Idem for kMin block face.
     ! viscKMaxPointer(2:il,2:jl): Idem for kMax block face.

     type(viscSubfaceType), dimension(:), pointer :: viscSubface

     integer(kind=intType), dimension(:,:), pointer :: viscIMinPointer
     integer(kind=intType), dimension(:,:), pointer :: viscIMaxPointer
     integer(kind=intType), dimension(:,:), pointer :: viscJMinPointer
     integer(kind=intType), dimension(:,:), pointer :: viscJMaxPointer
     integer(kind=intType), dimension(:,:), pointer :: viscKMinPointer
     integer(kind=intType), dimension(:,:), pointer :: viscKMaxPointer
     !
     !        ****************************************************************
     !        *                                                              *
     !        * Mesh related variables.                                      *
     !        *                                                              *
     !        ****************************************************************
     !
     !  x(0:ie,0:je,0:ke,3)  - xyz locations of grid points in block.
     !  xOld(nOld,:,:,:,:)   - Coordinates on older time levels;
     !                         only needed for unsteady problems on
     !                         deforming grids. Only allocated on
     !                         the finest grid level. The blank
     !                         dimensions are equal to the dimensions
     !                         of x.
     !  sI(0:ie,1:je,1:ke,3) - Projected areas in the i-coordinate
     !                         direction. Normals point in the
     !                         direction of increasing i.
     !  sJ(1:ie,0:je,1:ke,3) - Projected areas in the j-coordinate
     !                         direction. Normals point in the
     !                         direction of increasing j.
     !  sK(1:ie,1:je,0:ke,3) - Projected areas in the k-coordinate
     !                         direction. Normals point in the
     !                         direction of increasing k.
     !  vol(0:ib,0:jb,0:kb)  - Cell volumes. The second level halo
     !                         is present for a multigrid option.
     !  volOld(nold,2:il,..) - Volumes on older time levels; only
     !                         needed for unsteady problems on
     !                         deforming grids. Only allocated on
     !                         the finest grid level.
     !  uv(2,2:il,2:jl,2:kl) - Parametric location on elemID for each cell. 
     !                         Only used for fast wall distance calcs. 
     ! elemID(2:il,2:jl,2:kl)- Element ID each face is attached it
     !  porI(1:il,2:jl,2:kl) - Porosity in the i direction.
     !  porJ(2:il,1:jl,2:kl) - Porosity in the j direction.
     !  porK(2:il,2:jl,1:kl) - Porosity in the k direction.
     !
     !  indFamilyI(:,:,:)  - Index of the i-face in the arrays
     !                       to compute the local mass flow
     !                       for a family or sliding mesh interface.
     !                       Dimension is (1:il,2:jl,2:kl).
     !  indFamilyJ(:,:,:)  - Idem for the j-faces.
     !                       Dimension is (2:il,1:jl,2:kl).
     !  indFamilyK(:,:,:)  - Idem for the k-faces.
     !                       Dimension is (2:il,2:jl,1:kl)
     !  factFamilyI(:,:,:) - Corresponding factor to make sure
     !                       that the massflow is defined positive
     !                       when it enters the block and to define
     !                       the mass flow of the entire wheel
     !                       instead of a sector. Hence the possible
     !                       values or -nSlices and nSlices, where
     !                       nSlices or the number of sections to
     !                       obtain the full wheel.
     !  factFamilyJ(:,:,:) - Idem for the j-faces.
     !  factFamilyK(:,:,:) - Idem for the k-faces.
     !
     !  rotMatrixI(:,:,:,:,:) - Rotation matrix of the i-faces to
     !                          transform the velocity components
     !                          from Cartesian to local cylindrical.
     !                          This is needed only for problems with
     !                          rotational periodicity in combination
     !                          with an upwind scheme.
     !                          Dimension is (1:il,2:jl,2:kl,3,3).
     !  rotMatrixJ(:,:,:,:,:) - Idem for the j-faces.
     !                          Dimension is (2:il,1:jl,2:kl,3,3).
     !  rotMatrixK(:,:,:,:,:) - Idem for the k-faces.
     !                          Dimension is (2:il,2:jl,1:kl,3,3).
     !
     !  blockIsMoving      - Whether or not the block is moving.
     !  addGridVelocities  - Whether or not the face velocities
     !                       are allocated and set.
     !  sFaceI(0:ie,je,ke) - Dot product of the face velocity and
     !                       the normal in i-direction.
     !  sFaceJ(ie,0:je,ke) - Idem in j-direction.
     !  sFaceK(ie,je,0:ke) - Idem in k-direction.

     real(kind=realType), dimension(:,:,:,:),   pointer :: x, xtmp
     real(kind=realType), dimension(:,:,:,:,:),   pointer :: Xold
     real(kind=realType), dimension(:,:,:,:), pointer :: sI, sJ, sK
     real(kind=realType), dimension(:,:,:),   pointer :: vol
     real(kind=realType), dimension(:,:,:,:), pointer :: volOld
     real(kind=realType), dimension(:,:,:,:), pointer ::dadidata
     real(kind=realType), dimension(:,:,:,:), pointer :: uv
     integer(kind=intType), dimension(:,:,:), pointer :: elemID

     integer(kind=porType), dimension(:,:,:), pointer :: porI, porJ, porK
     integer(kind=intType), dimension(:,:,:), pointer :: indFamilyI, indFamilyJ, indFamilyK
     integer(kind=intType), dimension(:,:,:), pointer :: factFamilyI, factFamilyJ, factFamilyK
     real(kind=realType), dimension(:,:,:,:,:), pointer :: rotMatrixI, rotMatrixJ, rotMatrixK

     logical :: blockIsMoving, addGridVelocities

     real(kind=realType), dimension(:,:,:), pointer :: sFaceI, sFaceJ, sFaceK


     ! *******************************
     ! Added by HDN
     ! *******************************
     ! xALE(0:ie,0:je,0:ke,3)              - Temporary storage of x so that
     !                                       intermediate meshes can be stored in
     !                                       x directly
     ! sVeloIALE(0:ie,1:je,1:ke,3)
     ! sVeloJALE(1:ie,0:je,1:ke,3)
     ! sVeloKALE(1:ie,1:je,0:ke,3)         - Storage of surface velocities at one
     !                                       time step
     ! sIALE(0:nALEsteps,0:ie,1:je,1:ke,3)
     ! sJALE(0:nALEsteps,1:ie,0:je,1:ke,3)
     ! sKALE(0:nALEsteps,1:ie,1:je,0:ke,3) - Storage of sI, sJ, sK for intermediate
     !                                       meshes
     ! sFaceIALE(0:nALEsteps,0:ie,je,ke)
     ! sFaceJALE(0:nALEsteps,ie,0:je,ke)
     ! sFaceKALE(0:nALEsteps,ie,je,0:ke)   - Storage of sFaceI, sFaceJ, sFaceK for
     !                                       intermediate meshes
     real(kind=realType), dimension(:,:,:,:),   pointer :: xALE
     real(kind=realType), dimension(:,:,:,:),   pointer :: sVeloIALE, sVeloJALE, sVeloKALE
     real(kind=realType), dimension(:,:,:,:,:), pointer :: sIALE, sJALE, sKALE
     real(kind=realType), dimension(:,:,:,:),   pointer :: sFaceIALE, sFaceJALE, sFaceKALE




     !
     !        ****************************************************************
     !        *                                                              *
     !        * Flow variables.                                              *
     !        *                                                              *
     !        ****************************************************************
     !
     ! w(0:ib,0:jb,0:kb,1:nw)       - The set of independent variables
     !                                w(i,j,k,1:nwf) flow field
     !                                variables, which are rho, u, 
     !                                v, w and rhoE. In other words
     !                                the velocities  are stored and
     !                                not the momentum!!!!
     !                                w(i,j,k,nt1:nt2) turbulent 
     !                                variables; also the primitive
     !                                variables are stored.
     ! wOld(nOld,2:il,2:jl,2:kl,nw) - Solution on older time levels,
     !                                needed for the time integration
     !                                for unsteady problems. In
     !                                constrast to w, the conservative
     !                                variables are stored in wOld for
     !                                the flow variables; the turbulent
     !                                variables are always the
     !                                primitive ones.
     !                                Only allocated on the finest
     !                                mesh.
     ! p(0:ib,0:jb,0:kb)            - Static pressure.
     ! gamma(0:ib,0:jb,0:kb)        - Specific heat ratio; only
     !                                allocated on the finest grid.
     ! rlv(0:ib,0:jb,0:kb)          - Laminar viscosity; only
     !                                allocated on the finest mesh
     !                                and only for viscous problems.
     ! rev(0:ib,0:jb,0:kb)          - Eddy viscosity; only
     !                                allocated rans problems with
     !                                eddy viscosity models.
     ! s(1:ie,1:je,1:ke,3)          - Mesh velocities of the cell
     !                                centers; only for moving mesh
     !                                problems.

     real(kind=realType), dimension(:,:,:,:),   pointer :: w, wtmp
     real(kind=realType), dimension(:,:,:,:,:), pointer :: dw_deriv
     real(kind=realType), dimension(:,:,:,:,:), pointer :: wOld
     real(kind=realType), dimension(:,:,:),     pointer :: p, gamma, aa
     real(kind=realType), dimension(:,:,:),     pointer :: rlv, rev
     real(kind=realType), dimension(:,:,:,:),   pointer :: s
     real(kind=realType), dimension(:,:,:),   pointer :: shockSensor

     ! Nodal Fluxes: ux,uy,uz,vx,vy,vz,wx,wy,wz,qx,qy,qz(il, jl, kl)
     real(kind=realType), dimension(:, :, :), pointer:: ux, uy, uz
     real(kind=realType), dimension(:, :, :), pointer:: vx, vy, vz 
     real(kind=realType), dimension(:, :, :), pointer:: wx, wy, wz 
     real(kind=realType), dimension(:, :, :), pointer:: qx, qy, qz 
     



     !
     !        ****************************************************************
     !        *                                                              *
     !        * Residual and multigrid variables.                            *
     !        *                                                              *
     !        ****************************************************************
     !
     ! dw(0:ib,0:jb,0:kb,1:nw)   - Values of convective and combined
     !                             flow residuals. Only allocated on
     !                             the finest mesh.
     ! fw(0:ib,0:jb,0:kb,1:nwf)  - values of artificial dissipation
     !                             and viscous residuals.
     !                             Only allocated on the finest mesh.

     ! dwOldRK(:,2:il,2:jl,2:kl,nw) - Old residuals for the time
     !                                accurate Runge-Kutta schemes.
     !                                The first dimension is
     !                                nRKStagesUnsteady - 1.Only
     !                                allocated on the finest level
     !                                and only in unsteady mode for
     !                                Runge-Kutta schemes.

     ! w1(1:ie,1:je,1:ke,1:nMGVar) - Values of the mg variables
     !                               upon first entry to a coarser
     !                               mesh; only allocated on the
     !                               coarser grids. The variables
     !                               used to compute the multigrid
     !                               corrections are rho, u, v, w
     !                               and p; the rhoE value is used
     !                               for unsteady problems only.
     ! p1(1:ie,1:je,1:ke)          - Value of the pressure upon
     !                               first entry to a coarser grid;
     !                               only allocated on the coarser
     !                               grids.
     ! wr(2:il,2:jl,2:kl,1:nMGVar) - Multigrid forcing terms; only 
     !                               allocated on the coarser grids.
     !                               The forcing term of course
     !                               contains conservative residuals,
     !                               at least for the flow variables.
     ! shockSensor(0:ib,0:jb,0:kb)   Precomputed sensor value for shock 
     !                               that is *NOT* differentated. 
     ! scratch(0:ib,0:jb,0:kb,5)     Scratch space for the turbulence
     !                               models. NOMINALLY this could use
     !                               dw and the code was nominally setup
     !                               for this originally. However, this 
     !                               complicates reverse mode sensitivities
     !                               So we use this instead. 

     real(kind=realType), dimension(:,:,:),     pointer :: p1
     real(kind=realType), dimension(:,:,:,:),   pointer :: dw, fw
     real(kind=realType), dimension(:,:,:,:),   pointer :: dwtmp, dwtmp2
     real(kind=realType), dimension(:,:,:,:,:), pointer :: dwOldRK
     real(kind=realType), dimension(:,:,:,:),   pointer :: w1, wr
     real(kind=realType), dimension(:,:,:,:),   pointer :: scratch


     ! *******************************
     ! Added by HDN
     ! Used for ALE. Only allocated on the finest mesh.
     ! Extra dim is used to store initial residuals
     ! *******************************
     ! dwALE(0:nALEsteps,0:ib,0:jb,0:kb,1:nw)   - Values of ONLY the convective flux
     !                                            of intermediate meshes.
     ! fwALE(0:nALEsteps,0:ib,0:jb,0:kb,1:nwf)  - values of ONLY the artificial
     !                                            dissipation of intermediate meshes.
     real(kind=realType), dimension(:,:,:,:,:),   pointer :: dwALE, fwALE



     ! mgIFine(2:il,2) - The two fine grid i-cells used for the
     !                   restriction of the solution and residual to
     !                   the coarse grid. Only on the coarser grids.
     ! mgJFine(2:jl,2) - Idem for j-cells.
     ! mgKFine(2:kl,2) - Idem for k-cells.

     ! mgIWeight(2:il) - Weight for the residual restriction in
     !                   in i-direction. Value is either 0.5 or 1.0,
     !                   depending whether mgIFine(,1) is equal to
     !                   or differs from mgIFine(,2).
     ! mgJWeight(2:jl) - Idem for weights in j-direction.
     ! mgKWeight(2:kl) - Idem for weights in k-direction.

     ! mgICoarse(2:il,2) - The two coarse grid i-cells used for the
     !                     interpolation of the correction to the
     !                     fine grid. Not on the coarsest grid.
     ! mgJCoarse(2:jl,2) - Idem for j-cells.
     ! mgKCoarse(2:kl,2) - Idem for k-cells.

     integer(kind=intType), dimension(:,:), pointer :: mgIFine
     integer(kind=intType), dimension(:,:), pointer :: mgJFine
     integer(kind=intType), dimension(:,:), pointer :: mgKFine

     real(kind=realType),   dimension(:),   pointer :: mgIWeight
     real(kind=realType),   dimension(:),   pointer :: mgJWeight
     real(kind=realType),   dimension(:),   pointer :: mgKWeight

     integer(kind=intType), dimension(:,:), pointer :: mgICoarse
     integer(kind=intType), dimension(:,:), pointer :: mgJCoarse
     integer(kind=intType), dimension(:,:), pointer :: mgKCoarse

     ! iCoarsened - How this block was coarsened in i-direction.
     ! jCoarsened - How this block was coarsened in j-direction.
     ! kCoarsened - How this block was coarsened in k-direction.

     integer(kind=porType) :: iCoarsened, jCoarsened, kCoarsened

     ! iCo: Indicates whether or not i grid lines are present on the
     !      coarse grid; not allocated for the coarsest grid.
     ! jCo: Idem in j-direction.
     ! kCo: Idem in k-direction.

     logical, dimension(:), pointer :: iCo, jCo, kCo
     !
     !        ****************************************************************
     !        *                                                              *
     !        * Time-stepping and spectral radii variables.                  *
     !        * only allocated on the finest grid.                           *
     !        *                                                              *
     !        ****************************************************************
     !
     ! wn(2:il,2:jl,2:kl,1:nMGVar) - Values of the update variables
     !                               at the beginning of the RungeKutta
     !                               iteration. Only allocated for
     !                               RungeKutta smoother.
     ! pn(2:il,2:jl,2:kl)          - The pressure for the RungeKutta
     !                               smoother.
     ! dtl(1:ie,1:je,1:ke)         - Time step
     ! radI(1:ie,1:je,1:ke)        - Spectral radius in i-direction.
     ! radJ(1:ie,1:je,1:ke)        - Spectral radius in j-direction.
     ! radK(1:ie,1:je,1:ke)        - Spectral radius in k-direction.

     real(kind=realType), dimension(:,:,:,:), pointer :: wn
     real(kind=realType), dimension(:,:,:),   pointer :: pn
     real(kind=realType), dimension(:,:,:),   pointer :: dtl
     real(kind=realType), dimension(:,:,:),   pointer :: radI, radJ, radK

     !
     !        ****************************************************************
     !        *                                                              *
     !        * Variables for Iso/Surface Slice generation                   *
     !        *                                                              *
     !        ****************************************************************
     ! fc(1:ie,1:je,1:ke) - cell center values of the function to be iso-valued
     ! fn(1:il,1:jl,1:kl) - node values of the function to be iso-valued
     ! Note these are are only allocated temporaily during solution writing.

     real(kind=realType), dimension(:, :, :), pointer :: fc
     real(kind=realType), dimension(:, :, :), pointer :: fn
     


     !
     !        ****************************************************************
     !        *                                                              *
     !        * Turbulence model variables.                                  *
     !        *                                                              *
     !        ****************************************************************
     !
     ! d2Wall(2:il,2:jl,2:kl) - Distance from the center of the cell
     !                          to the nearest viscous wall.

     real(kind=realType), dimension(:,:,:), pointer :: d2Wall, filterDES


     ! bmti1(je,ke,nt1:nt2,nt1:nt2): Matrix used for the implicit
     !                               boundary condition treatment of
     !                               the turbulence equations at the
     !                               iMin boundary. Only allocated on
     !                               the finest level and for the 1st
     !                               spectral solution.
     ! bmti2(je,ke,nt1:nt2,nt1:nt2): Idem for the iMax boundary.
     ! bmtj1(ie,ke,nt1:nt2,nt1:nt2): Idem for the jMin boundary.
     ! bmtj2(ie,ke,nt1:nt2,nt1:nt2): Idem for the jMax boundary.
     ! bmtk1(ie,je,nt1:nt2,nt1:nt2): Idem for the kMin boundary.
     ! bmtk2(ie,je,nt1:nt2,nt1:nt2): Idem for the kMax boundary.

     real(kind=realType), dimension(:,:,:,:), pointer :: bmti1, bmti2
     real(kind=realType), dimension(:,:,:,:), pointer :: bmtj1, bmtj2
     real(kind=realType), dimension(:,:,:,:), pointer :: bmtk1, bmtk2


     ! bvti1(je,ke,nt1:nt2): RHS vector used for the implicit
     !                       boundary condition treatment of the
     !                       turbulence equations at the iMin
     !                       boundary. Only allocated on the finest
     !                       level and for the 1st spectral solution.
     ! bvti2(je,ke,nt1:nt2): Idem for the iMax boundary.
     ! bvtj1(ie,ke,nt1:nt2): Idem for the jMin boundary.
     ! bvtj2(ie,ke,nt1:nt2): Idem for the jMax boundary.
     ! bvti2(je,ke,nt1:nt2): Idem for the iMax boundary.
     ! bvtk1(ie,ke,nt1:nt2): Idem for the kMin boundary.
     ! bvtk2(ie,ke,nt1:nt2): idem for the kMax boundary.

     real(kind=realType), dimension(:,:,:), pointer :: bvti1, bvti2
     real(kind=realType), dimension(:,:,:), pointer :: bvtj1, bvtj2
     real(kind=realType), dimension(:,:,:), pointer :: bvtk1, bvtk2
     !
     !        ****************************************************************
     !        *                                                              *
     !        * Relation to the original cgns grid.                          *
     !        *                                                              *
     !        ****************************************************************
     !
     ! sectionID      - The section of the grid this block belongs to.
     ! cgnsBlockID    - Block/zone number of the cgns grid to which
     !                  this block is related.
     ! iBegOr, iEndOr - Range of points of this block in the
     ! jBegOr, jEndOr   corresponding cgns block, i.e. for this block
     ! kBegOr, kEndOr   iBegOr <= i <= iEndOr, jBegOr <= j <= jEndOr, 
     !                  kBegOr <= k <= kEndOr.
     !                  It is of course possible that the entire
     !                  block is stored.

     integer(kind=intType) :: cgnsBlockID, sectionID
     integer(kind=intType) :: iBegOr, iEndOr, jBegOr, jEndOr
     integer(kind=intType) :: kBegOr, kEndOr
     type(surfaceNodeWeightArray) , dimension(6) :: nodalWeights 
     !
     !        ****************************************************************
     !        *                                                              *
     !        * Adjoint solver variables.                                    *
     !        *                                                              *
     !        ****************************************************************
     !
     ! globalNode(ib:ie,jb:je,kb:ke):  Global node numbering.
     ! globalCell(0:ib,0:jb,0:kb):     Global cell numbering.
     ! color(0:ib,0:jb,0:kb)     :     Temporary coloring array used for 
     !                                 forward mode AD/FD calculations
     integer(kind=intType), dimension(:,:,:), pointer :: globalNode
     integer(kind=intType), dimension(:,:,:), pointer :: globalCell
     integer(kind=intType), dimension(:,:,:), pointer :: color

     integer(kind=intType),dimension(:),pointer::ifaceptb
     integer(kind=intType),dimension(:),pointer::iedgeptb

  end type blockType

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Array of all blocks at all multigrid levels and spectral sols. *
  !      *                                                                *
  !      ******************************************************************
  !
  ! nDom:            total number of computational blocks.
  ! flowDoms(:,:,:): array of blocks. Dimensions are
  !                  (nDom,nLevels,nTimeIntervalsSpectral)

  integer(kind=intType) :: nDom


#ifdef USE_TAPENADE
  ! This is never actually compiled...just make tapenade think it
  ! isn't allocatable
  type(blockType), dimension(nn:nn,1,ntimeIntervalsSpectral) :: flowDoms
#else
  type(blockType), allocatable, dimension(:,:,:) :: flowDoms
  type(blockType), allocatable, dimension(:,:,:) :: flowDomsd
  type(blockType), allocatable, dimension(:,:,:) :: flowDomsb
#endif

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Additional info needed in the flow solver.                     *
  !      *                                                                *
  !      ******************************************************************
  !
  ! nCellGlobal(nLev) - Global number of cells on every mg level.

  integer(kind=intType), allocatable, dimension(:) :: nCellGlobal

  contains


  logical function lessEqualFringeType(g1, g2)
    
    !        ****************************************************************
    !        *                                                              *
    !        * lessEqual returns .true. if g1 <= g2 and .false. otherwise.  *
    !        * The comparison is firstly based on the processor ID of the   *
    !        * donor, then the block, then then the I, J, K                 *
    !        *                                                              *
    !        ****************************************************************
    !
    implicit none
    !
    !        Function arguments.
    !
    type(fringeType), intent(in) :: g1, g2
    !
    ! Compare the donor processors first. If not equal,
    ! set lessEqual appropriately and return.
    
    if(g1%donorProc < g2%donorProc) then
       lessEqualfringeType = .true.
       return
    else if(g1%donorProc > g2%donorProc) then
       lessEqualfringeType = .false.
       return
    endif

    ! Donor processors are identical. Now we check the block
    
    if(g1%donorBlock < g2%donorBlock) then
       lessEqualfringeType = .true.
       return
    else if(g1%donorBlock > g2%donorBlock) then
       lessEqualfringeType = .false.
       return
    endif
    
    ! Compare the indices of the halo. First k, then j and
    ! finally i.
    
    if(g1%dK < g2%dK) then
       lessEqualfringeType = .true.
       return
    else if(g1%dK > g2%dK) then
       lessEqualfringeType = .false.
       return
    endif
    
    if(g1%dJ < g2%dJ) then
       lessEqualfringeType = .true.
       return
    else if(g1%dJ > g2%dJ) then
       lessEqualfringeType = .false.
       return
    endif
         
    if(g1%dI < g2%dI) then
       lessEqualfringeType = .true.
       return
    else if(g1%dI > g2%dI) then
       lessEqualfringeType = .false.
       return
    endif

    ! Both entities are identical. So set lessEqual to .true.
    
    lessEqualfringeType = .true.
    
  end function lessEqualFringeType

  logical function lessFringeType(g1, g2)
    
    !        ****************************************************************
    !        *                                                              *
    !        * less returns .true. if g1 <= g2 and .false. otherwise.  *
    !        * The comparison is firstly based on the processor ID of the   *
    !        * donor, then the block, then then the I, J, K                 *
    !        *                                                              *
    !        ****************************************************************
    !
    implicit none
    !
    !        Function arguments.
    !
    type(fringeType), intent(in) :: g1, g2
    !
    ! Compare the donor processors first. If not equal,
    ! set less appropriately and return.
    
    if(g1%donorProc < g2%donorProc) then
       lessfringeType = .true.
       return
    else if(g1%donorProc > g2%donorProc) then
       lessfringeType = .false.
       return
    endif

    ! Donor processors are identical. Now we check the block
    
    if(g1%donorBlock < g2%donorBlock) then
       lessfringeType = .true.
       return
    else if(g1%donorBlock > g2%donorBlock) then
       lessfringeType = .false.
       return
    endif
    
    ! Compare the indices of the halo. First k, then j and
    ! finally i.
    
    if(g1%dK < g2%dK) then
       lessfringeType = .true.
       return
    else if(g1%dK > g2%dK) then
       lessfringeType = .false.
       return
    endif
    
    if(g1%dJ < g2%dJ) then
       lessfringeType = .true.
       return
    else if(g1%dJ > g2%dJ) then
       lessfringeType = .false.
       return
    endif
         
    if(g1%dI < g2%dI) then
       lessfringeType = .true.
       return
    else if(g1%dI > g2%dI) then
       lessfringeType = .false.
       return
    endif

    ! Both entities are identical. So set less to .False.
    
    lessFringeType = .False.
    
  end function lessFringeType

end module block
