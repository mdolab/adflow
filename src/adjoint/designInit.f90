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
      use inputTimeSpectral !nTimeIntervalsSpectral
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
      ! - Mach Number,mach NumberGrid
      ! - Rotational Rate (3)
      ! - Rotational Center (3)
      ! - Moment Reference Point (3)
      ! - will eventually include Flight Dynamic Derivatives.

      nDesignExtra = 13!7
  
      ! Determine the total number of design variables.
      !   + nDesignSpatial (= 3 * nNodesGlobal)
      !   + nDesignExtra

      nDesignSpatial = (3 * nNodesGlobal)*nTimeIntervalsSpectral
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
      if( .not. allocated(functionName))then
         allocate(functionName(nCostFunction),stat=ierr)
      endif
      if(ierr /= 0)                       &
           call terminate("designInit", &
           "Memory allocation failure for functionName")

      if( .not. allocated(functionValue))then
         allocate(functionValue(nCostFunction),stat=ierr)
         if(ierr /= 0)                       &
              call terminate("designInit", &
              "Memory allocation failure for functionValue")
      end if
      if( .not. allocated(functionGrad))then
         allocate(functionGrad (nCostFunction,nDesignExtra),stat=ierr)
         if(ierr /= 0)                       &
              call terminate("designInit", &
              "Memory allocation failure for functionGrad")
      endif
!!$      if( .not. allocated(adjoint))then
!!$         allocate(adjoint(nCostFunction,nw * nCellsGlobal),stat=ierr)
!!$         if(ierr /= 0)                       &
!!$              call terminate("designInit", &
!!$              "Memory allocation failure for adjoint")
!!$      endif 
!      if( .not. allocated(functionGradSpatial))then
!         allocate(functionGradSpatial(nCostFunction,nDesignSpatial),stat=ierr)
!         if(ierr /= 0)                       &
!              call terminate("designInit", &
!              "Memory allocation failure for functionGradSpatial")
!      endif 

!!$      if( .not. allocated(functionGradStruct))then
!!$         allocate(functionGradStruct(nCostFunction,nDesignSpatial),stat=ierr)
!!$ if(ierr /= 0)                       &
!!$              call terminate("designInit", &
!!$              "Memory allocation failure for functionGradStruct")
!!$      endif 
!!$      if( .not. allocated(functionGradCoupling))then
!!$         allocate(functionGradCoupling(nCostFunction,nDesignSpatial),stat=ierr)
!!$         if(ierr /= 0)                       &
!!$              call terminate("designInit", &
!!$              "Memory allocation failure for functionGradCoupling")
!!$      endif
!!$      if( .not. allocated(functionGradCouplingExp))then
!!$         allocate(functionGradCouplingExp(nCostFunction,nDesignSpatial),stat=ierr) 
!!$         if(ierr /= 0)                       &
!!$              call terminate("designInit", &
!!$              "Memory allocation failure for functionGradCouplingExp")
!!$      endif
      
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

!!$      if( .not. allocated(xDesignVarName))then
!!$         allocate(xDesignVarName (nDesignExtra),stat=ierr)
!!$         if(ierr /= 0)                       &
!!$              call terminate("designInit", &
!!$              "Memory allocation failure for xDesignVarName")
!!$      endif
!!$      if( .not. allocated(xDesignVar))then
!!$         allocate(xDesignVar     (nDesignExtra),stat=ierr)
!!$         if(ierr /= 0)                       &
!!$              call terminate("designInit", &
!!$              "Memory allocation failure for xDesignVar")
!!$      endif
!!$      if( .not. allocated(xDesignVarLower))then
!!$         allocate(xDesignVarLower(nDesignExtra),stat=ierr)
!!$         if(ierr /= 0)                       &
!!$              call terminate("designInit", &
!!$              "Memory allocation failure for xDesignVarLower")
!!$      endif
!!$      if( .not. allocated(xDesignVarUpper))then
!!$         allocate(xDesignVarUpper(nDesignExtra),stat=ierr)
!!$         if(ierr /= 0)                       &
!!$              call terminate("designInit", &
!!$              "Memory allocation failure for xDesignVarUpper")
!!$      endif
!!$      
!!$      ! Set the design variable names.
!!$      
!!$      ! Angle of attack and side-slip angle.
!!$
!!$      xDesignVarName(nDesignAOA) = trim(varNameAOA)
!!$      xDesignVarName(nDesignSSA) = trim(varNameSSA)

 

      ! Formats.

  10  format(i5)
  20  format(a,a)

      end subroutine designInit
