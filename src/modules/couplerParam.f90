!
!      ******************************************************************
!      *                                                                *
!      * File:          couplerParam.f90                                *
!      * Author:        Seonghyeon Hahn, Edwin van der Weide            *
!      * Starting date: 01-31-2005                                      *
!      * Last modified: 02-08-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       module couplerParam
!
!      ******************************************************************
!      *                                                                *
!      * Definition of some parameters which are required for the       *
!      * coupling process.                                              *
!      *                                                                *
!      ******************************************************************
!
       use precision
       use constants
       implicit none
       save

       ! MachIni  : Mach number for flow initialization.
       ! pIni     : Pressure for flow initialization.
       ! rhoIni   : Density for flow initialization.
       ! velDirIni: Velocity direction for initialization.

       real(kind=realType) :: MachIni, pIni, rhoIni
       real(kind=realType), dimension(3) :: velDirIni

       ! maxCplNameLen  : Maximum length of strings for the names
       !                  of codes and interfaces.
       ! codeName       : Name given to a specific instance of SUmb.
       ! cplGetCoarseSol: Logical flag whether coarse-level
       !                  solutions are obtained from the coupler
       !                  or not.

       integer, parameter :: maxCplNameLen = 80
       character(len=maxCplNameLen) :: codeName
       logical :: cplGetCoarseSol

       ! nTetraTrue : Number of tetrahedra in the local domain.
       ! nPyraTrue  : Idem for pyramids.
       ! nPrismTrue : Idem for prisms.
       ! nHexaTrue  : Idem for hexahedra.
       ! nNodesTrue : Idem for nodes.
       ! nTetraAlloc: max(1,nTetraTrue), which is to avoid a zero-size
       !              array that may have some conflict with python.
       ! nPyraAlloc : Idem for pyramids.
       ! nPrismAlloc: Idem for prisms.
       ! nHexaAlloc : Idem for hexahedra.
       ! nNodesAlloc: Idem for nodes. Note that local solutions must be
       !              interpolated before being provided to the coupler,
       !              because SUmb uses the cell-centered configuration.

       integer(kind=intType) :: nTetraTrue, nPyraTrue, nPrismTrue, &
                                nHexaTrue, nNodesTrue
       integer(kind=intType) :: nTetraAlloc, nPyraAlloc, nPrismAlloc, &
                                nHexaAlloc, nNodesAlloc

       ! nDataSUmb       : Number of flow variables which this particular
       !                   instance of SUmb can provide to the coupler.
       !                   See the subroutine sumb_getMeshData.
       ! dataNamesSUmb(:): Names of flow variables which this particular
       !                   instance of SUmb can provide to the coupler.
       !                   See the subroutine sumb_getMeshData.
       ! iwSUmb(:)       : Flow-variable indices of SUmb for the
       !                   corresponding element of dataNamesSUmb(:).
       !                   See the subroutine sumb_getMeshData.

       integer(kind=intType) :: nDataSUmb
       character(len=maxCplNameLen), allocatable, &
                                     dimension(:) :: dataNamesSUmb
       integer(kind=intType), allocatable, dimension(:) :: iwSUmb

       end module couplerParam
