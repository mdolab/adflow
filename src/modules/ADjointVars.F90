!
!     ******************************************************************
!     *                                                                *
!     * File:          ADjointVars.F90                                 *
!     * Author:        Andre C. Marta,C.A.(Sandy) Mader                *
!     * Starting date: 07-21-2006                                      *
!     * Last modified: 01-17-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      module ADjointVars
!
!     ******************************************************************
!     *                                                                *
!     * This module contains the variables used by discrete adjoint    *
!     * solver, such as logicals, mappings, etc.                       *
!     *                                                                *
!     ******************************************************************
!
      use constants
      implicit none
!
!     ******************************************************************
!     *                                                                *
!     * Cost functions.                                                *
!     *                                                                *
!     ******************************************************************
!
      ! nCostFunction Number of cost/constraint functions.

      integer(kind=intType), parameter :: nCostFunction = 14_intType

      integer(kind=intType), parameter :: costFuncLiftCoef = 1_intType,&
                                          costFuncDragCoef = 2_intType,&
					  costFuncForceXCoef = 3_intType,&
                                          costFuncForceYCoef = 4_intType,&
                                          costFuncForceZCoef = 5_intType,&
                                          costFuncMomXCoef = 6_intType,&
                                          costFuncMomYCoef = 7_intType,&
                                          costFuncMomZCoef = 8_intType,&
                                          costFuncCmzAlpha  = 9_intType,&
                                          costFuncCm0      = 10_intType,&
                                          costFuncClAlpha  = 11_intType,&
                                          costFuncCl0      = 12_intType,&
                                          costFuncCdAlpha  = 13_intType,&
                                          costFuncCd0      = 14_intType

      ! Cost function names.

      character(len=maxStringLen), parameter :: &
                                         costFuncNameLift = "LiftCoef"
      character(len=maxStringLen), parameter :: &
                                         costFuncNameDrag = "DragCoef"
      character(len=maxStringLen), parameter :: &
                                         costFuncNameForceX = "ForceXCoef"
      character(len=maxStringLen), parameter :: &
                                         costFuncNameForceY = "ForceYCoef"
      character(len=maxStringLen), parameter :: &
                                         costFuncNameForceZ = "ForceZCoef"
      character(len=maxStringLen), parameter :: &
                                         costFuncNameMomX = "MomXCoef"
      character(len=maxStringLen), parameter :: &
                                         costFuncNameMomY = "MomYCoef"
      character(len=maxStringLen), parameter :: &
                                         costFuncNameMomZ = "MomZCoef"

      ! functionName  Name of the cost/function. Size  [nCostFunction].
      ! functionValue Value of the cost/function. Size  [nCostFunction].
      ! functionGrad  Gradient values of the cost function w.r.t. the
      !               design variables.
      !               Size [nCostFunction,nDesignExtra].

      character(len=maxStringLen), allocatable, dimension(:) :: &
                                                          functionName
      real(kind=realType), allocatable, dimension(:)   :: functionValue
      real(kind=realType), allocatable, dimension(:,:) :: functionGrad
      real(kind=realType), allocatable, dimension(:,:) :: ADjoint
      real(kind=realType), allocatable, dimension(:,:) :: functionGradSpatial
      real(kind=realType), allocatable, dimension(:,:) :: functionGradSurfaceDV
      real(kind=realType), allocatable, dimension(:,:) :: functionGradSurfaceDisp
      real(kind=realType), allocatable, dimension(:,:) :: functionGradStruct
      real(kind=realType), allocatable, dimension(:,:) :: functionGradCoupling
      real(kind=realType), allocatable, dimension(:,:) :: functionGradCouplingExp
!
!     ******************************************************************
!     *                                                                *
!     * Design variables.                                              *
!     *                                                                *
!     ******************************************************************
!
      ! nDesignAOA     Extra design variables: angle of atack
      ! nDesignSSA     Extra design variables: side-slip angle
      ! nDesignExtra   Extra number of design variables:
      ! nDesignSpatial Spatial design variables in grid (3*nNodesGlobal)
      ! nDesign        Total number of design variables
      !               = nNodesGlobal (for sigma)
      !               + 3*nNodesGlobal (for coord) + nDesignExtra

      integer(kind=intType), parameter :: nDesignAOA     = 1_intType
      integer(kind=intType), parameter :: nDesignSSA     = 2_intType
      integer(kind=intType), parameter :: nDesignMach     = 3_intType
      integer(kind=intType), parameter :: nDesignMachGrid = 4_intType
      integer(kind=intType), parameter :: nDesignRotX     = 5_intType
      integer(kind=intType), parameter :: nDesignRotY     = 6_intType	
      integer(kind=intType), parameter :: nDesignRotZ     = 7_intType
      integer(kind=intType)            :: nDesignExtra
      integer(kind=intType)	       :: nDesignSpatial 
      integer(kind=intType)            :: nDesign

      ! Design variable names.

      character(len=maxStringLen), parameter :: &
                                         varNameAOA   = "AngleOfAttack"
      character(len=maxStringLen), parameter :: &
                                         varNameSSA   = "SideSlipAngle"


      ! xDesignVarName   Name of the extra of design variables
      ! xDesignVar       Value of the extra of design variables
      ! xDesignVarLower  Lower bound value of the extra of design vars.
      ! xDesignVarUpper  Upper bound value of the extra of design vars.

      character(len=maxStringLen), allocatable, dimension(:) :: &
                                         xDesignVarName

      real(kind=realType), allocatable, dimension(:) :: &
                                                     xDesignVar, &
                                                     xDesignVarLower, &
                                                     xDesignVarUpper
!
!     ******************************************************************
!     *                                                                *
!     * Global numbering.                                              *
!     *                                                                *
!     ******************************************************************
!
      ! nNodesGlobal  Total number of nodes.
      ! nNodesLocal   Number of nodes owned by the processor.
      ! nOffsetLocal  Global node number offset per processor.

      integer(kind=intType) :: nNodesGlobal, nNodesLocal, nNodeOffsetLocal
      integer(kind=intType) :: nCellsGlobal, nCellsLocal, nCellOffsetLocal
      integer(kind=intType) :: nSurfNodesGlobal, nSurfNodesLocal

      ! ****************************************************************
      ! Finite-difference approximation step size.

      ! Relative/absolute finite-difference step for design variables
      ! to estimate d J / d alpha(n) and d R / d alpha(n)

      real(kind=realType), parameter :: adjEpsFd = 1.0e-4_realType
      real(kind=realType), parameter :: adjRelFd = 1.0e-5_realType
      real(kind=realType), parameter :: adjAbsFd = 1.0e-5_realType

      end module ADjointVars
