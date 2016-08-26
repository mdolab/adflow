
      module ADjointVars

      use constants, only : intType
      implicit none
!
      ! nDesignExtra Extra number of "extra" design variables (listed
      !                below)

      integer(kind=intType) :: nDesignAOA   
      integer(kind=intType) :: nDesignSSA    
      integer(kind=intType) :: nDesignMach    
      integer(kind=intType) :: nDesignMachGrid
      integer(kind=intType) :: nDesignRotX    
      integer(kind=intType) :: nDesignRotY    
      integer(kind=intType) :: nDesignRotZ    
      integer(kind=intType) :: nDesignRotCenX 
      integer(kind=intType) :: nDesignRotCenY  
      integer(kind=intType) :: nDesignRotCenZ   
      integer(kind=intType) :: nDesignPointRefX  
      integer(kind=intType) :: nDesignPointRefY   
      integer(kind=intType) :: nDesignPointRefZ   
      integer(kind=intType) :: nDesignLengthRef
      integer(kind=intType) :: nDesignSurfaceRef
      integer(kind=intType) :: nDesignDissError
      integer(kind=intType) :: nDesignPressure
      integer(kind=intType) :: nDesignTemperature
      integer(kind=intType) :: nDesignDensity
      integer(kind=intType) :: nDesignExtra = 0
      
      ! nNodesGlobal  Total number of nodes on each level
      ! nNodesLocal   Number of nodes owned by the processor on each level
      ! nOffsetLocal  Global node number offset per processor on each level

      ! Note: We're going to assume no more than 20 multigrid
      ! levels...this really should NEVER be exceed...
      integer(kind=intType), parameter :: maxLevels = 20
      integer(kind=intType), dimension(maxLevels) :: nNodesGlobal, nNodesLocal, nNodeOffsetLocal
      integer(kind=intType), dimension(maxLevels) :: nCellsGlobal, nCellsLocal, nCellOffsetLocal
      logical :: derivVarsAllocated = .False.
    end module ADjointVars
