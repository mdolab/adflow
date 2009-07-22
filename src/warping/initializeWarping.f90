!
! ***********************************
! *  File: initializeWarping..f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 12-12-2008
! *  Modified: 12-12-2008
! ***********************************

subroutine initializeWarping

use blockPointers
use communication
use monitor
use inputIO
use mdDataLocal
use mdData
implicit none


! Local Arguments

integer(kind=intType)::sps = 1,famid = 0,startind,endInd,level=1




!begin execution

!Call old indexing and coordinate functions to populate
!Baseline datastructures in mdData and mdDataLocal
call mdCreateSurfCoorListLocal(sps,famid,startInd,endInd)

call mdCreateSurfIndListLocal(famID,startInd,endInd)

call mdCreateNsurfNodes

!Run the index synchronization and reduction
call synchronizeIndices

!Create the reduced surface coordinate list for python interface
call mdCreateGlobalReducedSurfaceList

!Run the Adjoint preprocessing to setup global volume numbering and
!initialize the various PETSc variables

! Assertion testing, memory allocation and global number indexing.
call preprocessingADjoint(level)

!initialize Petsc
call initializePETSc

!setup additional PETSc Warping Variables
call warpingInitializePETSc

call createWarpingPETScMat


end subroutine initializeWarping


