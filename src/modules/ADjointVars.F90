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
      ! nDesignExtra   Extra number of "extra" design variables 
      !                (listed below)

      integer(kind=intType) :: nDesignAOA = -1  
      integer(kind=intType) :: nDesignSSA = -1   
      integer(kind=intType) :: nDesignMach = -1   
      integer(kind=intType) :: nDesignMachGrid = -1
      integer(kind=intType) :: nDesignRotX = -1   
      integer(kind=intType) :: nDesignRotY = -1    
      integer(kind=intType) :: nDesignRotZ = -1    
      integer(kind=intType) :: nDesignRotCenX = -1
      integer(kind=intType) :: nDesignRotCenY = -1 
      integer(kind=intType) :: nDesignRotCenZ = -1  
      integer(kind=intType) :: nDesignPointRefX = -1  
      integer(kind=intType) :: nDesignPointRefY = -1   
      integer(kind=intType) :: nDesignPointRefZ = -1   
      
      integer(kind=intType) :: nDesignExtra = 0

      real(kind=realType),dimension(:),allocatable :: dIda
      
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
      real(kind=realType) :: timings(20)
      ! ****************************************************************
      ! Finite-difference approximation step size.

      ! Relative/absolute finite-difference step for design variables
      ! to estimate d J / d alpha(n) and d R / d alpha(n)

      real(kind=realType), parameter :: adjEpsFd = 1.0e-4_realType
      real(kind=realType), parameter :: adjRelFd = 1.0e-5_realType
      real(kind=realType), parameter :: adjAbsFd = 1.0e-5_realType


    end module ADjointVars
