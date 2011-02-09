!
!      ******************************************************************
!      *                                                                *
!      * File:          outputMod.f90                                   *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 08-13-2004                                      *
!      * Last modified: 03-29-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       module outputMod
!
!      ******************************************************************
!      *                                                                *
!      * This local module contains variables used when writing the     *
!      * grid and solution files.                                       *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none

       ! nblocksCGNSblock(0:cgnsNDom): The number of local blocks per
       !                               cgns block in cumulative storage
       !                               format.
       ! blocksCGNSblock(nDom):        The corresponding local block ID's.

       integer(kind=intType), dimension(:), allocatable :: nblocksCGNSblock
       integer(kind=intType), dimension(:), allocatable :: blocksCGNSblock

       ! nDomPerProc(0:nProc):  The number of flow domains for each 
       !                        processor in cumulative storage format.
       ! IDsBegOrAllDoms(4,..): The CGNS block numbers and the beginning
       !                        indices for all of the flow domains on
       !                        every processor.

       integer(kind=intType), dimension(:),   allocatable :: nDomPerProc
       integer(kind=intType), dimension(:,:), allocatable :: IDsBegOrAllDoms

       ! nGridsToWrite:   Number of grid files to write.
       ! nVolSolToWrite:  Number of volume solution files to write.
       !                  For CGNS nVolSolToWrite == nGridsToWrite.
       ! nSurfSolToWrite: Number of surface solution files to write.

       integer(kind=intType) :: nGridsToWrite
       integer(kind=intType) :: nVolSolToWrite
       integer(kind=intType) :: nSurfSolToWrite

       ! gridFileNames(nGridsToWrite):      Names of the grid files to
       !                                    write.
       ! volSolFileNames(nVolSolToWrite):   Names of the volume solution
       !                                    files to write.
       ! surfSolFileNames(nSurfSolToWrite): Names of the surface solution
       !                                    files to write.
       ! fileIDs(nFilesToWrite):            Identifiers for the files to
       !                                    write. As the grids, volume
       !                                    solution and surface solution
       !                                    files are written one after
       !                                    the other, only one set is
       !                                    needed.
       ! cgnsBases(nFilesToWrite):          The CGNS base IDs of the
       !                                    files.

       character(len=maxStringLen), dimension(:), allocatable :: &
                        gridFileNames, volSolFileNames, surfSolFileNames

       integer, dimension(:), allocatable :: fileIDs
       integer, dimension(:), allocatable :: cgnsBases

       ! useLinksInCGNS:     Whether or not to use links in CGNS between
       !                     the grid and volume solution files. If not,
       !                     the grid and solution are written in the
       !                     same file.
       ! writePlot3DConn:    Whether or not to write the Plot3D
       !                     connectivity file.
       ! writeFormatInParam: Whether or not the format must be written
       !                     in the parameter file if this file must
       !                     be updated automatically.

       logical :: useLinksInCGNS
       logical :: writePlot3DConn
       logical :: writeFormatInParam

       end module outputMod
