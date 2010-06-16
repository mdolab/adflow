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
use mddata
use mddatalocal
implicit none

!subroutine Arguments
integer(kind=intType)::i,ierr
! Local Arguments

integer(kind=intType)::sps = 1,startind,endInd,level=1
call mdCreateNsurfNodesLocal
call mdCreateSurfListLocal(sps)
call synchronizeIndices
call flagImplicitEdgesAndFacesSurface

!Run the Adjoint preprocessing to setup global volume numbering and
!initialize the various PETSc variables

! Assertion testing, memory allocation and global number indexing.
call preprocessingADjoint(level)
call warpingInitializePETSc
call createWarpingPETScMat
end subroutine initializeWarping


