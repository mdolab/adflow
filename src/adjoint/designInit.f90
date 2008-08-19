!
!     ******************************************************************
!     *                                                                *
!     * File:          designInit.f90                                  *
!     * Authors:       C.A.(Sandy) Mader                               *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 01-14-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine designInit
!
!     ******************************************************************
!     *                                                                *
!     * designInit determines the number of design variables,          *
!     *   allocates the arrays that store the cost function names,     *
!     *   values and gradients and fills in the function names. It     *
!     *   also allocates the arrays that store the design variables,   *
!     *   their names, lower and upper bounds, filling their names.    *
!     *                                                                *
!     * The variables are split in the spatial variables (x) and       *
!     *   and extra variables, such as angle of attack.                *
!     *                                                                *
!     * The extra design variables are ordered as:                     *
!     *  - angle of attack                                             *
!     *  - side-slip angle                                             *
!     *  - ...                                                         *
!     *                                                                *
!     * The routines that compute dR/da (setupGradientMatrixExtra)     *
!     * and dJ/da (setupGradientRHSExtra) have to be consistent with   *
!     * this ordering...                                               *
!     *                                                                *
!     ******************************************************************
!
      use ADjointVars     ! nDesign, nDesignExtra, nDesignSpatial
                          ! nNodesGlobal, nDesignDipoles
      use communication   ! myID
      use flowVarRefState ! magnetic
      implicit none
!
!     Local variables.
!
      character(len=maxStringLen) :: varName, tmp
      integer(kind=intType)       :: n, nn,ierr
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
      ! Determine the number of extra variables. This includes:
      ! - angle of attack
      ! - side-slip angle
      ! - Mach Number
      ! - will eventually include Flight Dynamic Derivatives.

      nDesignExtra = 3
  
      ! Determine the total number of design variables.
      !   + nDesignSpatial (= 3 * nNodesGlobal)
      !   + nDesignExtra

      nDesignSpatial = 3 * nNodesGlobal
      nDesign        = nDesignExtra + nDesignSpatial

!
!     ******************************************************************
!     *                                                                *
!     * Cost functions.                                                *
!     *                                                                *
!     ******************************************************************
!
      ! Allocate memory to store the function names, their values
      ! and gradients.
      allocate(functionName(nCostFunction),stat=ierr)
      if(ierr /= 0)                       &
           call terminate("designInit", &
           "Memory allocation failure for functionName")
      allocate(functionValue(nCostFunction))
      allocate(functionGrad (nCostFunction,nDesignExtra))
      allocate(adjoint(nCostFunction,nw * nNodesGlobal))
      allocate(functionGradSpatial(nCostFunction,nDesignSpatial))
      allocate(functionGradCoupling(nCostFunction,nDesignSpatial))
      allocate(functionGradCouplingExp(nCostFunction,nDesignSpatial))
      
      !print *,'function storage allocated'

      ! Set the cost function names.

      functionName(costFuncLiftCoef) = costFuncNameLift
      functionName(costFuncDragCoef) = costFuncNameDrag
      functionName(costFuncForceXCoef) = costFuncNameForceX
      functionName(costFuncForceYCoef) = costFuncNameForceY
      functionName(costFuncForceZCoef) = costFuncNameForceZ
      functionName(costFuncMomXCoef) = costFuncNameMomX
      functionName(costFuncMomYCoef) = costFuncNameMomY
      functionName(costFuncMomZCoef) = costFuncNameMomZ
!
!     ******************************************************************
!     *                                                                *
!     * Design variables.                                              *
!     *                                                                *
!     ******************************************************************
!
      ! Allocate memory to store the design variable names, values,
      ! lower and upper bounds.

      allocate(xDesignVarName (nDesignExtra))
      allocate(xDesignVar     (nDesignExtra))
      allocate(xDesignVarLower(nDesignExtra))
      allocate(xDesignVarUpper(nDesignExtra))

      ! Set the design variable names.

      ! Angle of attack and side-slip angle.

      xDesignVarName(nDesignAOA) = trim(varNameAOA)
      xDesignVarName(nDesignSSA) = trim(varNameSSA)

 

      ! Formats.

  10  format(i5)
  20  format(a,a)

      end subroutine designInit
