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
      use costFunctions
      implicit none
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
      integer(kind=intType), parameter :: nDesignRotCenX     = 8_intType
      integer(kind=intType), parameter :: nDesignRotCenY     = 9_intType
      integer(kind=intType), parameter :: nDesignRotCenZ     = 10_intType
      integer(kind=intType), parameter :: nDesignPointRefX     = 11_intType
      integer(kind=intType), parameter :: nDesignPointRefY     = 12_intType
      integer(kind=intType), parameter :: nDesignPointRefZ     = 13_intType
      
      integer(kind=intType)            :: nDesignExtra
      integer(kind=intType)	       :: nDesignSpatial 
      integer(kind=intType)            :: nDesign

      
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
