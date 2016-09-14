
      module ADjointVars

      use constants, only : intType
      implicit none

      ! Indices of the extra deisgn design variables.
      integer(kind=intType), parameter :: iAlpha=1
      integer(kind=intType), parameter :: iBeta=2
      integer(kind=intType), parameter :: iMach=3
      integer(kind=intType), parameter :: iMachGrid=4
      integer(kind=intType), parameter :: iRotX=5
      integer(kind=intType), parameter :: iRotY=6
      integer(kind=intType), parameter :: iRotZ=7
      integer(kind=intType), parameter :: iRotCenX=8
      integer(kind=intType), parameter :: iRotCenY=9
      integer(kind=intType), parameter :: iRotCenZ=10
      integer(kind=intType), parameter :: iPointRefX=11
      integer(kind=intType), parameter :: iPointRefY=12
      integer(kind=intType), parameter :: iPointRefZ=13
      integer(kind=intType), parameter :: iPressure=14
      integer(kind=intType), parameter :: iTemperature=15
      integer(kind=intType), parameter :: iDensity=16

      integer(kind=intType), parameter :: nDesignExtra=16

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
