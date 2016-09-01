module indirectHalo
  !
  !       This local module contains the derived data type used to       
  !       determine the indirect halo's as well as an array of this type.
  !
  use precision
  implicit none
  save

  public
  private :: lessEqualIndirectHaloType
  private :: lessIndirectHaloType
  !
  !       The definition of the derived data type indirectHaloType.      
  !
  type indirectHaloType

     ! myBlock      : Local block ID of the halo.
     ! myI, myJ, myK: i,j,k indices of the halo.
     ! myDirectHalo : Index in the haloListType where the
     !                corresponding direct halo is stored.
     ! levOfInd     : Level of indirectness.
     ! donorProc    : Processor where donor of the direct halo is
     !                stored. In case this halo is a boundary
     !                halo, donorProc is set to -1.

     integer(kind=intType) :: myBlock
     integer(kind=intType) :: myI, myJ, myK
     integer(kind=intType) :: myDirectHalo
     integer(kind=intType) :: levOfInd
     integer(kind=intType) :: donorProc

  end type indirectHaloType

  ! Interface for the extension of the operators <= and <.
  ! These are needed for the sorting of indirectHaloType.

  interface operator(<=)
     module procedure lessEqualIndirectHaloType
  end interface operator(<=)

  interface operator(<)
     module procedure lessIndirectHaloType
  end interface operator(<)

  ! nIndHalo         : Number of indirect halo's to be treated.
  ! indHalo(nIndHalo): The indirect halo's.

  integer(kind=intType) :: nIndHalo
  type(indirectHaloType), dimension(:), allocatable :: indHalo

  ! nLevOfInd               : Number of levels of indirectness
  ! nHaloPerLev(0:nLevOfInd): Number of indirect halo's per level
  !                           of indirectness; stored in
  !                           cumulative storage format.
  ! nHaloPerProc(0:nProc)   : Number of indirect halo's per
  !                           processor for a given level of
  !                           indirectness; cumulative storage.
  !                           nHaloPerProc(0) is not 0,
  !                           because of the presence of boundary
  !                           halo's, which get proc ID -1.

  integer(kind=intType) :: nLevOfInd
  integer(kind=intType), dimension(:), allocatable :: nHaloPerLev
  integer(kind=intType), dimension(:), allocatable :: nHaloPerProc

contains
  !
  !         Functions to simulate the operators <= and <.                
  !
  logical function lessEqualIndirectHaloType(g1, g2)
    !
    !         This function returns .true. if g1 <= g2 and .false.         
    !         otherwise. The comparison is firstly based on the level of   
    !         indirectness followed by the donor processor, the            
    !         corresponding direct halo, my block ID and finally the i, j  
    !         and k indices.                                               
    !
    implicit none
    !
    !        Function arguments.
    !
    type(indirectHaloType), intent(in) :: g1, g2

    ! Compare the level of indirectness. If not equal, set
    ! lessEqual appropriately and return.

    if(g1%levOfInd < g2%levOfInd) then
       lessEqualIndirectHaloType = .true.
       return
    else if(g1%levOfInd > g2%levOfInd) then
       lessEqualIndirectHaloType = .false.
       return
    endif

    ! Compare the donor processors.

    if(g1%donorProc < g2%donorProc) then
       lessEqualIndirectHaloType = .true.
       return
    else if(g1%donorProc > g2%donorProc) then
       lessEqualIndirectHaloType = .false.
       return
    endif

    ! Compare the direct halo.

    if(g1%myDirectHalo < g2%myDirectHalo) then
       lessEqualIndirectHaloType = .true.
       return
    else if(g1%myDirectHalo > g2%myDirectHalo) then
       lessEqualIndirectHaloType = .false.
       return
    endif

    ! Compare my block ID.

    if(g1%myBlock < g2%myBlock) then
       lessEqualIndirectHaloType = .true.
       return
    else if(g1%myBlock > g2%myBlock) then
       lessEqualIndirectHaloType = .false.
       return
    endif

    ! Finally compare the halo indices. Start with k.

    if(g1%myK < g2%myK) then
       lessEqualIndirectHaloType = .true.
       return
    else if(g1%myK > g2%myK) then
       lessEqualIndirectHaloType = .false.
       return
    endif

    ! The j index.

    if(g1%myJ < g2%myJ) then
       lessEqualIndirectHaloType = .true.
       return
    else if(g1%myJ > g2%myJ) then
       lessEqualIndirectHaloType = .false.
       return
    endif

    ! The i index.

    if(g1%myI < g2%myI) then
       lessEqualIndirectHaloType = .true.
       return
    else if(g1%myI > g2%myI) then
       lessEqualIndirectHaloType = .false.
       return
    endif

    ! Both entities are identical. So set lessEqual to .true.

    lessEqualIndirectHaloType = .true.

  end function lessEqualIndirectHaloType

  !        ================================================================

  logical function lessIndirectHaloType(g1, g2)
    !
    !         This function returns .true. If g1 < g2 and .false.          
    !         otherwise. It is basically the same as the lessEqual         
    !         function, except that the equality is now considered as      
    !         .false.                                                      
    !
    implicit none
    !
    !        Function arguments.
    !
    type(indirectHaloType), intent(in) :: g1, g2

    ! Compare the level of indirectness. If not equal, set
    ! lessIndirectHaloType appropriately and return.

    if(g1%levOfInd < g2%levOfInd) then
       lessIndirectHaloType = .true.
       return
    else if(g1%levOfInd > g2%levOfInd) then
       lessIndirectHaloType = .false.
       return
    endif

    ! Compare the donor processors.

    if(g1%donorProc < g2%donorProc) then
       lessIndirectHaloType = .true.
       return
    else if(g1%donorProc > g2%donorProc) then
       lessIndirectHaloType = .false.
       return
    endif

    ! Compare the direct halo.

    if(g1%myDirectHalo < g2%myDirectHalo) then
       lessIndirectHaloType = .true.
       return
    else if(g1%myDirectHalo > g2%myDirectHalo) then
       lessIndirectHaloType = .false.
       return
    endif

    ! Compare my block id.

    if(g1%myBlock < g2%myBlock) then
       lessIndirectHaloType = .true.
       return
    else if(g1%myBlock > g2%myBlock) then
       lessIndirectHaloType = .false.
       return
    endif

    ! Finally compare the halo indices. Start with k.

    if(g1%myK < g2%myK) then
       lessIndirectHaloType = .true.
       return
    else if(g1%myK > g2%myK) then
       lessIndirectHaloType = .false.
       return
    endif

    ! The j index.

    if(g1%myJ < g2%myJ) then
       lessIndirectHaloType = .true.
       return
    else if(g1%myJ > g2%myJ) then
       lessIndirectHaloType = .false.
       return
    endif

    ! The i index.

    if(g1%myI < g2%myI) then
       lessIndirectHaloType = .true.
       return
    else if(g1%myI > g2%myI) then
       lessIndirectHaloType = .false.
       return
    endif

    ! Both entities are identical.
    ! So set lessIndirectHaloType to .false.

    lessIndirectHaloType = .false.

  end function lessIndirectHaloType

end module indirectHalo
module haloList
  !
  !       This local module contains temporary variables to create the   
  !       list of halo cells and nodes.                                  
  !
  use precision
  implicit none
  save

  public
  private :: lessEqualHaloListType
  private :: lessHaloListType
  !
  !       The definition of the variables for the 3 lists.               
  !
  type haloListType

     ! myBlock:           local block ID of the halo.
     ! myI, myJ, myK:     i,j,k indices of the halo.
     ! donorProc :        processor where donor is stored. In case
     !                    the halo is a boundary halo, donorProc
     !                    is set to -1.
     ! donorBlock:        block ID of the donor. In case the halo
     !                    is a boundary halo donorBlock is set
     !                    to the corresponding boundary condition.
     ! dI, dJ, dK:        i,j,k indices of the donor.
     ! levOfInd:          level of indirectness.
     ! interp(..):        interpolants for the donor stencil; only
     !                    allocated for lists requiring this info.
     ! nPeriodicSubfaces: Number of periodic subfaces that are
     !                    crossed when going from the halo to the
     !                    donor. This is at most the level of
     !                    indirectness of the halo.
     ! periodicSubfaces:  The corresponding subfaces ID's according
     !                    to the sequence defined in periodicGlobal.

     integer(kind=intType) :: myBlock
     integer(kind=intType) :: myI, myJ, myK
     integer(kind=intType) :: donorProc, donorBlock
     integer(kind=intType) :: dI, dJ, dK
     integer(kind=intType) :: levOfInd
     real(kind=realType), dimension(:), pointer :: interp

     integer(kind=intType) :: nPeriodicSubfaces
     integer(kind=intType), dimension(:), pointer :: periodicSubfaces

  end type haloListType

  ! Interface for the extension of the operators <= and <.
  ! These are needed for the sorting of haloListType.

  interface operator(<=)
     module procedure lessEqualHaloListType
  end interface operator(<=)

  interface operator(<)
     module procedure lessHaloListType
  end interface operator(<)

  ! nCellHalo1st: # of 1st level cell halo's
  ! nCellHalo2nd: # of 2nd level cell halo's
  ! nNodeHalo1st: # of 1st level node halo's

  integer(kind=intType) :: nCellHalo1st, nCellHalo2nd
  integer(kind=intType) :: nNodeHalo1st

  ! iiCell1st: Counter variable for the 1st level cell halo's
  ! iiCell2nd: Counter variable for the 2nd level cell halo's
  ! iiNode1st: Counter variable for the 1st level node halo's

  integer(kind=intType) :: iiCell1st, iiCell2nd, iiNode1st

  ! cellHalo1st(nCellHalo1st) :: List of halo info for 1st
  !                              level cell halo's.
  ! cellHalo2nd(nCellHalo2nd) :: Idem for 2nd level cell halo's.
  ! nodeHalo1st(nNodeHalo1st) :: Idem for 1st level node halo's.

  type(haloListType), dimension(:), allocatable :: cellHalo1st
  type(haloListType), dimension(:), allocatable :: cellHalo2nd
  type(haloListType), dimension(:), allocatable :: nodeHalo1st

  ! transformCell(nCellHalo1st,3) :: Short hand for the transformation
  !                                  matrix between the halo and
  !                                  the donor for cell based halo's.
  !                                  In principle the size equals the
  !                                  number of faces (i.e. direct)
  !                                  halo's, but the difference is
  !                                  not so large.
  ! transformNode(nNodeHalo1st,3) :: Idem for the nodes.

  integer(kind=intType), dimension(:,:), allocatable :: transformCell
  integer(kind=intType), dimension(:,:), allocatable :: transformNode
  !
  !       The definition of the index variables, which store for each    
  !       i,j,k in the block the index in the corresponding list.        
  !       I know I'm wasting memory here (because only the halo's are    
  !       relevant), but that's not too much of a problem. The reason is 
  !       that neither the metrics nor the variables have been allocated 
  !       yet. So later on, much more memory is needed than the single   
  !       integer for each cell/node used here.                          
  !
  type indexListType

     ! entryList(:,:,:)  :: Corresponding entry in the list.
     !                      Dimensions are either 0:ie,0:je,0:ke for
     !                      the node or 0:ib,0:jb,0:kb for the cell
     !                      based halo's. The latter is then suited
     !                      for the 2nd level halo's.

     integer(kind=intType), dimension(:,:,:), pointer :: entryList

  end type indexListType

  ! nodeIndex(nDom) :: The node indices for every block.
  ! cellIndex(nDom) :: Idem for the cells.

  type(indexListType), allocatable, dimension(:) :: nodeIndex
  type(indexListType), allocatable, dimension(:) :: cellIndex

contains
  !
  !         Functions to simulate the operators <= and <.                
  !
  logical function lessEqualHaloListType(g1, g2)
    !
    !         lessEqual returns .true. if g1 <= g2 and .false. otherwise.  
    !         The comparison is firstly based on the processor ID of the   
    !         donor. After that it depends whether the halo is a boundary  
    !         halo or not. Note that boundary halo's have a donor processor
    !         if of -1, such that they are always first in the list.       
    !
    implicit none
    !
    !        Function arguments.
    !
    type(haloListType), intent(in) :: g1, g2

    ! Compare the donor processors first. If not equal,
    ! set lessEqual appropriately and return.

    if(g1%donorProc < g2%donorProc) then
       lessEqualHaloListType = .true.
       return
    else if(g1%donorProc > g2%donorProc) then
       lessEqualHaloListType = .false.
       return
    endif

    ! Donor processors are identical. Now it depends whether we are
    ! dealing with boundary halo's or not.

    boundary: if(g1%donorProc == -1) then ! And thus
       ! g2%donorProc == -1

       ! Both halo's are boundary halo's. Compare the block ID of
       ! the halo's.

       if(g1%myBlock < g2%myBlock) then
          lessEqualHaloListType = .true.
          return
       else if(g1%myBlock > g2%myBlock) then
          lessEqualHaloListType = .false.
          return
       endif

       ! Compare the boundary conditions, which are stored in
       ! donorBlock. Note that the sequence in BCTypes is such that
       ! the most important BC has the highest number.

       if(g1%donorBlock < g2%donorBlock) then
          lessEqualHaloListType = .true.
          return
       else if(g1%donorBlock > g2%donorBlock) then
          lessEqualHaloListType = .false.
          return
       endif

       ! As it is possible that indirect halo's need donor info from
       ! direct halo's or even indirect halo's with a smaller level
       ! of indirectness, compare the level of indirectness.

       if(g1%levOfInd < g2%levOfInd) then
          lessEqualHaloListType = .true.
          return
       else if(g1%levOfInd > g2%levOfInd) then
          lessEqualHaloListType = .false.
          return
       endif

       ! Compare the indices of the halo. First k, then j and
       ! finally i.

       if(g1%myK < g2%myK) then
          lessEqualHaloListType = .true.
          return
       else if(g1%myK > g2%myK) then
          lessEqualHaloListType = .false.
          return
       endif

       if(g1%myJ < g2%myJ) then
          lessEqualHaloListType = .true.
          return
       else if(g1%myJ > g2%myJ) then
          lessEqualHaloListType = .false.
          return
       endif

       if(g1%myI < g2%myI) then
          lessEqualHaloListType = .true.
          return
       else if(g1%myI > g2%myI) then
          lessEqualHaloListType = .false.
          return
       endif

       ! No need to compare anything else; g1 == g2.

    else boundary

       ! Both halo's are internal halo's, whose donor is stored on
       ! the same processor. Compare the donor blocks.

       if(g1%donorBlock < g2%donorBlock) then
          lessEqualHaloListType = .true.
          return
       else if(g1%donorBlock > g2%donorBlock) then
          lessEqualHaloListType = .false.
          return
       endif

       ! Also the blocks are identical. Compare the donor indices.
       ! First the k index.

       if(g1%dK < g2%dK) then
          lessEqualHaloListType = .true.
          return
       else if(g1%dK > g2%dK) then
          lessEqualHaloListType = .false.
          return
       endif

       ! The j index.

       if(g1%dJ < g2%dJ) then
          lessEqualHaloListType = .true.
          return
       else if(g1%dJ > g2%dJ) then
          lessEqualHaloListType = .false.
          return
       endif

       ! And the i index.

       if(g1%dI < g2%dI) then
          lessEqualHaloListType = .true.
          return
       else if(g1%dI > g2%dI) then
          lessEqualHaloListType = .false.
          return
       endif

       ! The donors are identical. Compare the halo's.
       ! First the block id.

       if(g1%myBlock < g2%myBlock) then
          lessEqualHaloListType = .true.
          return
       else if(g1%myBlock > g2%myBlock) then
          lessEqualHaloListType = .false.
          return
       endif

       ! Halo blocks are also identical. Finally compare the
       ! halo indices. Start with k.

       if(g1%myK < g2%myK) then
          lessEqualHaloListType = .true.
          return
       else if(g1%myK > g2%myK) then
          lessEqualHaloListType = .false.
          return
       endif

       ! The j index.

       if(g1%myJ < g2%myJ) then
          lessEqualHaloListType = .true.
          return
       else if(g1%myJ > g2%myJ) then
          lessEqualHaloListType = .false.
          return
       endif

       ! The i index.

       if(g1%myI < g2%myI) then
          lessEqualHaloListType = .true.
          return
       else if(g1%myI > g2%myI) then
          lessEqualHaloListType = .false.
          return
       endif

    endif boundary

    ! Both entities are identical. So set lessEqual to .true.

    lessEqualHaloListType = .true.

  end function lessEqualHaloListType

  !        ================================================================

  logical function lessHaloListType(g1, g2)
    !
    !         This function returns .true. if g1 < g2 and .false.          
    !         otherwise. It is basically the same as the lessEqual         
    !         function, except that the equality is now considered as      
    !         .false.                                                      
    !
    implicit none
    !
    !        Function arguments.
    !
    type(haloListType), intent(in) :: g1, g2

    ! Compare the donor processors first. If not equal,
    ! set the function appropriately and return.

    if(g1%donorProc < g2%donorProc) then
       lessHaloListType = .true.
       return
    else if(g1%donorProc > g2%donorProc) then
       lessHaloListType = .false.
       return
    endif

    ! Donor processors are identical. Now it depends whether we are
    ! dealing with boundary halo's or not.

    boundary: if(g1%donorProc == -1) then ! And thus
       ! g2%donorProc == -1

       ! Both halo's are boundary halo's. Compare the block ID of
       ! the halo's.

       if(g1%myBlock < g2%myBlock) then
          lessHaloListType = .true.
          return
       else if(g1%myBlock > g2%myBlock) then
          lessHaloListType = .false.
          return
       endif

       ! Compare the boundary conditions, which are stored in
       ! donorBlock. Note that the sequence in BCTypes is such that
       ! the most important bc has the highest number.

       if(g1%donorBlock < g2%donorBlock) then
          lessHaloListType = .true.
          return
       else if(g1%donorBlock > g2%donorBlock) then
          lessHaloListType = .false.
          return
       endif

       ! As it is possible that indirect halo's need donor info from
       ! direct halo's or even indirect halo's with a smaller level
       ! of indirectness, compare the level of indirectness.

       if(g1%levOfInd < g2%levOfInd) then
          lessHaloListType = .true.
          return
       else if(g1%levOfInd > g2%levOfInd) then
          lessHaloListType = .false.
          return
       endif

       ! Compare the indices of the halo. First k, then j and
       ! finally i.

       if(g1%myK < g2%myK) then
          lessHaloListType = .true.
          return
       else if(g1%myK > g2%myK) then
          lessHaloListType = .false.
          return
       endif

       if(g1%myJ < g2%myJ) then
          lessHaloListType = .true.
          return
       else if(g1%myJ > g2%myJ) then
          lessHaloListType = .false.
          return
       endif

       if(g1%myI < g2%myI) then
          lessHaloListType = .true.
          return
       else if(g1%myI > g2%myI) then
          lessHaloListType = .false.
          return
       endif

       ! No need to compare anything else. G1 == g2.

    else boundary

       ! Both halo's are internal halo's, whose donor is stored on
       ! the same processor. Compare the donor blocks.

       if(g1%donorBlock < g2%donorBlock) then
          lessHaloListType = .true.
          return
       else if(g1%donorBlock > g2%donorBlock) then
          lessHaloListType = .false.
          return
       endif

       ! Also the blocks are identical. Compare the donor indices.
       ! First the k index.

       if(g1%dK < g2%dK) then
          lessHaloListType = .true.
          return
       else if(g1%dK > g2%dK) then
          lessHaloListType = .false.
          return
       endif

       ! The j index.

       if(g1%dJ < g2%dJ) then
          lessHaloListType = .true.
          return
       else if(g1%dJ > g2%dJ) then
          lessHaloListType = .false.
          return
       endif

       ! And the i index.

       if(g1%dI < g2%dI) then
          lessHaloListType = .true.
          return
       else if(g1%dI > g2%dI) then
          lessHaloListType = .false.
          return
       endif

       ! The donors are identical. Compare the halo's.
       ! First the block id.

       if(g1%myBlock < g2%myBlock) then
          lessHaloListType = .true.
          return
       else if(g1%myBlock > g2%myBlock) then
          lessHaloListType = .false.
          return
       endif

       ! Halo blocks are also identical. Finally compare the
       ! halo indices. Start with k.

       if(g1%myK < g2%myK) then
          lessHaloListType = .true.
          return
       else if(g1%myK > g2%myK) then
          lessHaloListType = .false.
          return
       endif

       ! The j index.

       if(g1%myJ < g2%myJ) then
          lessHaloListType = .true.
          return
       else if(g1%myJ > g2%myJ) then
          lessHaloListType = .false.
          return
       endif

       ! The i index.

       if(g1%myI < g2%myI) then
          lessHaloListType = .true.
          return
       else if(g1%myI > g2%myI) then
          lessHaloListType = .false.
          return
       endif

    endif boundary

    ! Both entities are identical.
    ! So set lessHaloListType to .false.

    lessHaloListType = .false.

  end function lessHaloListType

end module haloList

module checkVolBlock
  !
  !       Local module, which contains the definition of the derived     
  !       datatype used to test for negative volumes in the grid.        
  !
  implicit none
  save

  type checkVolBlockType

     ! blockHasNegVol:              Whether or not the block
     !                              contains negative volumes.
     ! volumeIsNeg(2:il,2:jl,2:kl): Whether or not the owned volumes
     !                              are negative.

     logical :: blockHasNegVol
     logical, dimension(:,:,:), pointer :: volumeIsNeg

  end type checkVolBlockType

end module checkVolBlock
module periodicInfo
  !
  !       Local module that contains derived datatypes as well as arrays 
  !       of these derived datatypes to store information related to     
  !       periodicity.                                                   
  !
  use precision
  implicit none
  save

  public
  private :: lessCGNSPeriodicType
  private :: equalCGNSPeriodicType
  private :: lessPeriodicSubfacesHaloT
  private :: lessEqualPeriodicSubfacesHaloT
  private :: equalPeriodicSubfacesHaloT

  ! Definition of the derived data type for storing the periodic
  ! faces of the cgns grid a bit easier.

  type cgnsPeriodicType

     ! cgnsBlock   :: the block ID in the cgns grid.
     ! cgnsSubface :: the suface ID in this block.

     integer(kind=intType) :: cgnsBlock, cgnsSubface

  end type cgnsPeriodicType

  ! Interface for the extension of the operators < and ==.
  ! These are needed for the sorting and searching of
  ! cgnsPeriodicType.

  interface operator(<)
     module procedure lessCGNSPeriodicType
  end interface operator(<)

  interface operator(==)
     module procedure equalCGNSPeriodicType
  end interface operator(==)

  ! nPeriodicGlobal :: Total number of periodic faces in cgns grid.
  ! periodicGlobal  :: The corresponding faces.

  integer(kind=intType) :: nPeriodicGlobal
  type(cgnsPeriodicType), dimension(:), allocatable :: periodicGlobal

  ! Definition of the derived data type to store the periodic
  ! subfaces that are crossed when the halo and the donor are
  ! connected. The direction is from the halo to the donor.

  type periodicSubfacesHaloType

     ! internalHalo:      Whether or not the halo is an internal
     !                    halo, i.e. it is stored on the
     !                    same processor as the donor.
     ! indexInHaloList:   The corresponding index in either the
     !                    node or cell halo list.
     ! nPeriodicSubfaces: Number of periodic subfaces that are
     !                    crossed. This is at most the level of
     !                    indirectness of the halo.
     ! periodicSubfaces:  The corresponding subfaces ID's according
     !                    to the sequence defined in periodicGlobal.

     logical               :: internalHalo
     integer(kind=intType) :: indexInHaloList
     integer(kind=intType) :: nPeriodicSubfaces
     integer(kind=intType), dimension(:), pointer :: periodicSubfaces

  end type periodicSubfacesHaloType

  ! Interface for the extension of the operators <, <= and ==.
  ! This is needed for the sorting and comparing of variables of
  ! the type periodicSubfacesHaloType

  interface operator(<)
     module procedure lessPeriodicSubfacesHaloT
  end interface operator(<)

  interface operator(<=)
     module procedure lessEqualPeriodicSubfacesHaloT
  end interface operator(<=)

  interface operator(==)
     module procedure equalPeriodicSubfacesHaloT
  end interface operator(==)

  !=================================================================

contains

  !=================================================================
  !
  !         Functions to simulate the operators < and ==.                
  !
  logical function lessCGNSPeriodicType(g1, g2)
    !
    !         lessCGNSPeriodicType returns .true. if g1 is considered      
    !         smaller than g2. This comparison is first based on the block 
    !         ID followed by the subface id.                               
    !
    implicit none
    !
    !        Function arguments.
    !
    type(cgnsPeriodicType), intent(in) :: g1, g2

    ! Compare the block ID. If not equal set lessCGNSPeriodicType
    ! accordingly.

    if(g1%cgnsBlock < g2%cgnsBlock) then
       lessCGNSPeriodicType = .true.
       return
    else if(g1%cgnsBlock > g2%cgnsBlock) then
       lessCGNSPeriodicType = .false.
       return
    endif

    ! Block ID's are identical. Compare the subfaces.

    if(g1%cgnsSubface < g2%cgnsSubface) then
       lessCGNSPeriodicType = .true.
       return
    else if(g1%cgnsSubface > g2%cgnsSubface) then
       lessCGNSPeriodicType = .false.
       return
    endif

    ! Both objects are identical.
    ! Set lessCGNSPeriodicType to .false.

    lessCGNSPeriodicType = .false.

  end function lessCGNSPeriodicType

  !        ================================================================

  logical function equalCGNSPeriodicType(g1, g2)
    !
    !         equalCGNSPeriodicType returns .true. if g1 is considered     
    !         equal to g2, i.e. both the block and subface ID must match,  
    !
    implicit none
    !
    !        Function arguments.
    !
    type(cgnsPeriodicType), intent(in) :: g1, g2

    equalCGNSPeriodicType = .false.
    if(g1%cgnsBlock   == g2%cgnsBlock .and. &
         g1%cgnsSubface == g2%cgnsSubface)    &
         equalCGNSPeriodicType = .true.

  end function equalCGNSPeriodicType

  !        ================================================================

  logical function lessPeriodicSubfacesHaloT(g1, g2)
    !
    !         lessPeriodicSubfacesHaloT returns .true. if g1 is            
    !         considered smaller than g2.                                  
    !
    implicit none
    !
    !        Function arguments.
    !
    type(periodicSubfacesHaloType), intent(in) :: g1, g2
    !
    !        Local variables.
    !
    integer(kind=intType) :: nn, i1, i2

    ! First compare whether or not both g1 and g2 are internal
    ! halo's. Fortran does not allow a direct comparison of
    ! logicals and therefore the integers i1 and i2 are used.

    i1 = 1; if( g1%internalHalo ) i1 = 0
    i2 = 1; if( g2%internalHalo ) i2 = 0

    if(i1 < i2) then
       lessPeriodicSubfacesHaloT = .true.
       return
    else if(i1 > i2) then
       lessPeriodicSubfacesHaloT = .false.
       return
    endif

    ! Compare the number of periodic subfaces.

    if(g1%nPeriodicSubfaces < g2%nPeriodicSubfaces) then
       lessPeriodicSubfacesHaloT = .true.
       return
    else if(g1%nPeriodicSubfaces > g2%nPeriodicSubfaces) then
       lessPeriodicSubfacesHaloT = .false.
       return
    endif

    ! The number of periodic subfaces is the same. Compare the
    ! subfaces themselves. It is assumed that the subfaces are
    ! sorted in increading order. This can be done, because the
    ! periodic transformations are commuting matrices.

    do nn=1,g1%nPeriodicSubfaces
       if(g1%periodicSubfaces(nn) < g2%periodicSubfaces(nn)) then
          lessPeriodicSubfacesHaloT = .true.
          return
       else if(g1%periodicSubfaces(nn) > g2%periodicSubfaces(nn)) then
          lessPeriodicSubfacesHaloT = .false.
          return
       endif
    enddo

    ! The periodic subfaces are identical as well. Compare the
    ! indices in the list.

    if(g1%indexInHaloList < g2%indexInHaloList) then
       lessPeriodicSubfacesHaloT = .true.
       return
    else if(g1%indexInHaloList > g2%indexInHaloList) then
       lessPeriodicSubfacesHaloT = .false.
       return
    endif

    ! Both objects are the same. Return .false.

    lessPeriodicSubfacesHaloT = .false.

  end function lessPeriodicSubfacesHaloT

  !        ================================================================

  logical function lessEqualPeriodicSubfacesHaloT(g1, g2)
    !
    !         lessEqualPeriodicSubfacesHaloT returns .true. if g1 is       
    !         considered smaller than or equal to g2.                      
    !
    implicit none
    !
    !        Function arguments.
    !
    type(periodicSubfacesHaloType), intent(in) :: g1, g2
    !
    !        Local variables.
    !
    integer(kind=intType) :: nn, i1, i2

    ! First compare whether or not both g1 and g2 are internal
    ! halo's. Fortran does not allow a direct comparison of
    ! logicals and therefore the integers i1 and i2 are used.

    i1 = 1; if( g1%internalHalo ) i1 = 0
    i2 = 1; if( g2%internalHalo ) i2 = 0

    if(i1 < i2) then
       lessEqualPeriodicSubfacesHaloT = .true.
       return
    else if(i1 > i2) then
       lessEqualPeriodicSubfacesHaloT = .false.
       return
    endif

    ! Compare the number of periodic subfaces.

    if(g1%nPeriodicSubfaces < g2%nPeriodicSubfaces) then
       lessEqualPeriodicSubfacesHaloT = .true.
       return
    else if(g1%nPeriodicSubfaces > g2%nPeriodicSubfaces) then
       lessEqualPeriodicSubfacesHaloT = .false.
       return
    endif

    ! The number of periodic subfaces is the same. Compare the
    ! subfaces themselves. It is assumed that the subfaces are
    ! sorted in increading order. This can be done, because the
    ! periodic transformations are commuting matrices.

    do nn=1,g1%nPeriodicSubfaces
       if(g1%periodicSubfaces(nn) < g2%periodicSubfaces(nn)) then
          lessEqualPeriodicSubfacesHaloT = .true.
          return
       else if(g1%periodicSubfaces(nn) > g2%periodicSubfaces(nn)) then
          lessEqualPeriodicSubfacesHaloT = .false.
          return
       endif
    enddo

    ! The periodic subfaces are identical as well. Compare the
    ! indices in the list.

    if(g1%indexInHaloList < g2%indexInHaloList) then
       lessEqualPeriodicSubfacesHaloT = .true.
       return
    else if(g1%indexInHaloList > g2%indexInHaloList) then
       lessEqualPeriodicSubfacesHaloT = .false.
       return
    endif

    ! Both objects are the same. Return .true.

    lessEqualPeriodicSubfacesHaloT = .true.

  end function lessEqualPeriodicSubfacesHaloT

  !        ================================================================

  logical function equalPeriodicSubfacesHaloT(g1, g2)
    !
    !         equalPeriodicSubfacesHaloT returns .true. if g1 is           
    !         considered equal to g2. The equal operator is only used to   
    !         find the different number of periodic transformations in     
    !         determinePeriodicData. Hence only the periodic subfaces of   
    !         the halo's are compared and g1 and g2 are considered equal   
    !         if the subfaces are equal, even if other member variables    
    !         differ.                                                      
    !
    implicit none
    !
    !        Function arguments.
    !
    type(periodicSubfacesHaloType), intent(in) :: g1, g2
    !
    !        Local variables.
    !
    integer(kind=intType) :: nn

    if(g1%nPeriodicSubfaces /= g2%nPeriodicSubfaces) then
       equalPeriodicSubfacesHaloT = .false.
       return
    endif

    do nn=1,g1%nPeriodicSubfaces
       if(g1%periodicSubfaces(nn) /= g2%periodicSubfaces(nn)) then
          equalPeriodicSubfacesHaloT = .false.
          return
       endif
    enddo

    equalPeriodicSubfacesHaloT = .true.

  end function equalPeriodicSubfacesHaloT

end module periodicInfo
module bcHalo
  !
  !       This local module contains the derived datatype bcHaloType,    
  !       which is used to determine the boundary condition for an       
  !       indirect halo when the nearest direct halo's are all boundary  
  !       halo's.                                                        
  !
  use precision
  implicit none
  save

  public
  private :: lessEqualBCHaloType
  !
  !       The definition of the derived datatype.                        
  !
  type bcHaloType

     ! directHalo: Index in the haloListType where the
     !             corresponding direct halo is stored.
     ! BC:         Corresponding boundary condition.

     integer(kind=intType) :: directHalo, BC

  end type bcHaloType

  ! Interface for the extension of the operator <= needed for the
  ! sorting of bcHaloType. Note that the = operator does not
  ! need to be defined, because bcHaloType only contains
  ! primitive types.

  interface operator(<=)
     module procedure lessEqualBCHaloType
  end interface operator(<=)

contains
  !
  logical function lessEqualBCHaloType(g1, g2)
    !
    !         Function to simulate the operator <= for bcHaloType.         
    !         It first compares the boundary condition. If equal the index 
    !         of the direct halo is compared, although this is not really  
    !         important.                                                   
    !         LessEqual returns .true. if g1 <= g2 and .false. otherwise.  
    !
    implicit none
    !
    !        Function arguments.
    !
    type(bcHaloType), intent(in) :: g1, g2

    ! First compare the boundary conditions. Note that the sequence
    ! in BCTypes is such that the most important BC has the
    ! highest number.

    if(g1%BC < g2%BC) then
       lessEqualBCHaloType = .true.
       return
    else if(g1%BC > g2%BC) then
       lessEqualBCHaloType = .false.
       return
    endif

    ! Boundary conditions are equal. Just compare the index.

    if(g1%directHalo < g2%directHalo) then
       lessEqualBCHaloType = .true.
       return
    else if(g1%directHalo > g2%directHalo) then
       lessEqualBCHaloType = .false.
       return
    endif

    ! g1 and g2 are equal. Return .true.

    lessEqualBCHaloType = .true.

  end function lessEqualBCHaloType

  !        ================================================================

  subroutine sortBCHaloType(bcHaloArray, nn)
    !
    !         SortBCHaloType sorts the given number of BCHalo's in         
    !         increasing order. Note that this routine is called sort and  
    !         not qsort, because only an insertion sort is done here. The  
    !         reason is that nn <= 3 and thus an insertion sort is okay.   
    !
    implicit none
    !
    !        Subroutine arguments
    !
    integer(kind=intType), intent(in) :: nn
    type(bcHaloType), dimension(*), intent(inout) :: bcHaloArray
    !
    !        Local variables.
    !
    integer(kind=intType) :: i, j

    type(bcHaloType) :: a

    do j=1,nn
       a = bcHaloArray(j)
       do i=(j-1),1,-1
          if(bcHaloArray(i) <= a) exit
          bcHaloArray(i+1) = bcHaloArray(i)
       enddo
       bcHaloArray(i+1) = a
    enddo

  end subroutine sortBCHaloType

end module bcHalo

module coarse1to1Subface
  !
  !       This local module contains the derived datatype                
  !       coarse1to1SubfaceType, which is used to determine the 1 to 1   
  !       block boundaries for the coarser grids.                        
  !
  use precision
  implicit none
  save
  !
  !       The definition of the derived datatype.                        
  !
  type coarse1to1SubfaceType

     ! Nodal range in the three coordinates directions for the
     ! coarse grid subface.

     integer(kind=intType) :: iBeg, jBeg, kBeg, iEnd, jEnd, kEnd

     ! Processor and block id of the neighboring block.

     integer(kind=intType) :: neighProc, neighBlock

     ! Number of points in the three coordinate directions for the
     ! coarse grid donor subface.

     integer(kind=intType) :: ndi, ndj, ndk

     ! Corresponding i, j and k indices of the fine grid donor block
     ! for each of the coarse grid subface lines.

     integer(kind=intType), dimension(:), pointer :: idfine
     integer(kind=intType), dimension(:), pointer :: jdfine
     integer(kind=intType), dimension(:), pointer :: kdfine

  end type coarse1to1SubfaceType

  ! Number of 1 to 1 fine grid subfaces on this processor.

  integer(kind=intType) :: nSubface1to1

  ! Array of 1 to 1 subfaces.

  type(coarse1to1SubfaceType), dimension(:), allocatable :: subface1to1

end module coarse1to1Subface

!      ==================================================================

module coarseningInfo
  !
  !       This local module contains the derived datatype                
  !       coarseningInfoType, which stores for a given block the grid    
  !       lines to keep for the coarse grid.                             
  !
  type coarseningInfoType

     ! Logical, which indicate whether or not a fine grid 1 to 1
     ! block boundary is still a 1 to 1 block boundary on the
     ! coarse grid.

     logical, dimension(:), pointer :: coarseIs1to1

  end type coarseningInfoType

  ! Array to store the info for all the blocks.

  type(coarseningInfoType), dimension(:), allocatable :: coarseInfo

end module coarseningInfo

