!
! ***********************************
! *  File: initializeWarping..f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 12-12-2008
! *  Modified: 12-12-2008
! ***********************************

subroutine initializeWarping(famID)

use blockPointers
use communication
use monitor
use inputIO
use mdDataLocal
use mdData
implicit none

!subroutine Arguments
integer(kind=intType)::famid
! Local Arguments

integer(kind=intType)::sps = 1,startind,endInd,level=1

if(myid==0)print *,'in initializewarping...',famid

!Call old indexing and coordinate functions to populate
!Baseline datastructures in mdData and mdDataLocal
call mdCreateSurfCoorListLocal(sps,famid,startInd,endInd)
!print *,'famID, ind.',famID
call mdCreateSurfIndListLocal(famID,startInd,endInd)
!print *,'local indices created'
call mdCreateNsurfNodes
!print *,'nsurface nodes complete'
!Run the index synchronization and reduction
call synchronizeIndices
!print *,'synchronize indices complete'
!Create the reduced surface coordinate list for python interface
!call mdCreateGlobalReducedSurfaceList

!Run the Adjoint preprocessing to setup global volume numbering and
!initialize the various PETSc variables

! Assertion testing, memory allocation and global number indexing.
call preprocessingADjoint(level)

!initialize Petsc
call initializePETSc

!setup additional PETSc Warping Variables
call warpingInitializePETSc

call createWarpingPETScMat

!print *,'indices',mdSurfGlobalIndLocal,'xx',mdGlobalSurfxx
end subroutine initializeWarping


