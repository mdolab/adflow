subroutine preprocessingCustomOverset
  !
  !      ******************************************************************
  !      *                                                                *
  !      * preprocessing determines the communication patterns between    *
  !      * the processors for all the mg levels, computes the wall        *
  !      * distances and the metrics.                                     *
  !      *                                                                *
  !      ******************************************************************
  !
  use block
  use blockPointers
  use cgnsGrid
  use commMixing
  use commSliding
  use communication
  use inputPhysics
  use inputTimeSpectral
  use interfaceGroups
  use section
  use wallDistanceData
  use overset
  use utils, only : EChk, setPointers, setBufferSizes
  implicit none
  !
  !      Local variables.
  !
  integer :: ierr

  integer(kind=intType) :: nLevels, level, nn, mm, nsMin, nsMax
  integer(kind=intType) :: i, j, k, n, m, l, sps, nVolPos, nVolNeg
  real(kind=realType) :: xp, yp, zp, vp1, vp2, vp3, vp4, vp5, vp6
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Check that for the unsteady modes the number of periodic slices
  ! is identical for all sections.

  nsMin = sections(1)%nSlices
  nsMax = sections(1)%nSlices

  do nn=2,nSections
     nsMin = min(nsMin,sections(nn)%nSlices)
     nsMax = max(nsMax,sections(nn)%nSlices)
  enddo

  if((equationMode == unsteady .or. &
       equationMode == timeSpectral) .and. nsMin < nsMax) then

     if(myID == 0)                     &
          call returnFail("preprocessing", &
          "Different rotational periodicity encountered &
          &for time accurate computation")
     call mpi_barrier(SUmb_comm_world, ierr)

  endif

  ! Determine the number of multigrid levels needed in the
  ! computation, and allocate the memory for the node and cell
  ! communication patterns, including the memory copies for internal
  ! communication, and the global number of cells on each level.
  ! Note that this communication pattern does not change in time.

  nLevels = ubound(flowDoms,2)
  nn      = nLevels
  allocate(commPatternCell_1st(nn), commPatternCell_2nd(nn), &
       commPatternNode_1st(nn), internalCell_1st(nn),    &
       internalCell_2nd(nn),    internalNode_1st(nn),    &
       nCellGlobal(nn),         stat=ierr)
  if(ierr /= 0)                     &
       call returnFail("preprocessing", &
       "Memory allocation failure for commPatterns")

  ! Set the sizes here so that we know how to dealloc the stuff
  ! later on
  do i=1,nLevels
     commPatternCell_1st(i)%nPeriodic = 0
     commPatternCell_1st(i)%nProcSend = 0
     commPatternCell_1st(i)%nProcRecv = 0

     commPatternCell_2nd(i)%nPeriodic = 0
     commPatternCell_2nd(i)%nProcSend = 0
     commPatternCell_2nd(i)%nProcRecv = 0

     commPatternNode_1st(i)%nPeriodic = 0
     commPatternNode_1st(i)%nProcSend = 0
     commPatternNode_1st(i)%nProcRecv = 0

     internalCell_1st(i)%nPeriodic = 0
     internalCell_2nd(i)%nPeriodic = 0
     internalNode_1st(i)%nPeriodic = 0
  end do

  ! Allocate the memory for the sliding mesh communication pattern.
  ! This pattern changes in time and therefore each spectral time
  ! value has its own sliding mesh communication pattern.

  mm = nTimeIntervalsSpectral

  allocate(commSlidingCell_1st(nn,mm), &
       commSlidingCell_2nd(nn,mm), &
       intSlidingCell_1st(nn,mm),  &
       intSlidingCell_2nd(nn,mm),  stat=ierr)
  if(ierr /= 0)                     &
       call returnFail("preprocessing", &
       "Memory allocation failure for &
       &slidingCommPatterns")

  ! Allocate the memory for the communication of the mixing plane
  ! halo cells. As this type of "boundary condition" can only
  ! be applied for steady state computations, there is no
  ! dependency on the number of time spectral solutions.

  mm = nInterfaceGroups
  allocate(commPatternMixing(nn,mm,2), stat=ierr)
  if(ierr /= 0)                     &
       call returnFail("preprocessing", &
       "Memory allocation failure for &
       &commPatternMixing")

  ! Allocate the memory for the overset mesh communication pattern.
  ! This pattern changes in time and therefore each spectral time
  ! value has its own sliding mesh communication pattern.

  mm = nTimeIntervalsSpectral
  allocate(commPatternOverset(nn,mm), internalOverset(nn,mm), &
       overlapMatrix(nn, mm), stat=ierr)
  if(ierr /= 0)                     &
       call returnFail("preprocessing", &
       "Memory allocation failure for commOverset")

  ! Determine the fine grid 1 to 1 matching communication pattern.

  call determineCommPattern(1_intType)

  ! Initialize the send and receive buffer sizes to 0 and determine
  ! its size for the finest grid for the 1 to 1 communication.

  sendBufferSize_1to1  = 0
  recvBufferSize_1to1  = 0
  sendBufferSizeSlide  = 0
  recvBufferSizeSlide  = 0
  sendBufferSizeOver   = 0
  recvBufferSizeOver   = 0

  call setBufferSizes(1_intType, 1_intType, .true., .false., .false.)

  ! Loop to create the coarse grid levels.

  do level=2,nLevels

     ! Create the coarse grid blocks, its communication pattern, the
     ! coarse grid level 0 cooling parameters and check the
     ! communication buffer sizes.

     call createCoarseBlocks(level)
     call determineCommPattern(level)
     call coarseLevel0CoolingParameters(level)
     call setBufferSizes(level, 1_intType, .true., .false., .false.)

  enddo

  ! Synchronize the processors, just to be sure.

  call mpi_barrier(SUmb_comm_world, ierr)

  ! Allocate memory for the nonblocking point to point communication.

  allocate(sendBuffer(sendBufferSize), &
       recvBuffer(recvBufferSize), stat=ierr)
  if(ierr /= 0)                     &
       call returnFail("preprocessing", &
       "Memory allocation failure for sendBuffer &
       &and recvBuffer")

  ! Determine the cell range for the subfaces and initialize the
  ! arrays for the boundary condition data.
  ! Done for all grid levels.

  call cellRangeSubface
  call initBcdata

  ! Nullify the wallFringe poiter as initialization
  nullify(wallFringes, localWallFringes)

  !  ! Loop over the number of levels and perform a lot of tasks.
  !  ! See the corresponding subroutine header, although the
  !  ! names are pretty self-explaining

  level = 1
  call xhalo(level)
  call slidingComm(level, .true.)

  spectral: do sps=1,nTimeIntervalsSpectral
     domains: do nn=1,nDom
        nVolNeg = 0
        nVolPos = 0

        ib = flowDoms(nn,level,sps)%ib
        jb = flowDoms(nn,level,sps)%jb
        kb = flowDoms(nn,level,sps)%kb

        ! Allocate space for the volumes
        allocate(flowDoms(nn,level,sps)%vol(0:ib,0:jb,0:kb))

        call setPointers(nn, level, sps)

        ! Compute the volumes. The hexahedron is split into 6 pyramids
        ! whose volumes are computed. The volume is positive for a
        ! right handed block.
        ! Initialize the volumes to zero. The reasons is that the second
        ! level halo's must be initialized to zero and for convenience
        ! all the volumes are set to zero.

        vol = zero

        do k=1,ke
           n = k -1
           do j=1,je
              m = j -1
              do i=1,ie
                 l = i -1

                 ! Compute the coordinates of the center of gravity.

                 xp = eighth*(x(i,j,k,1) + x(i,m,k,1) &
                      +         x(i,m,n,1) + x(i,j,n,1) &
                      +         x(l,j,k,1) + x(l,m,k,1) &
                      +         x(l,m,n,1) + x(l,j,n,1))
                 yp = eighth*(x(i,j,k,2) + x(i,m,k,2) &
                      +         x(i,m,n,2) + x(i,j,n,2) &
                      +         x(l,j,k,2) + x(l,m,k,2) &
                      +         x(l,m,n,2) + x(l,j,n,2))
                 zp = eighth*(x(i,j,k,3) + x(i,m,k,3) &
                      +         x(i,m,n,3) + x(i,j,n,3) &
                      +         x(l,j,k,3) + x(l,m,k,3) &
                      +         x(l,m,n,3) + x(l,j,n,3))

                 ! Compute the volumes of the 6 sub pyramids. The
                 ! arguments of volpym must be such that for a (regular)
                 ! right handed hexahedron all volumes are positive.

                 vp1 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                      x(i,j,n,1), x(i,j,n,2), x(i,j,n,3), &
                      x(i,m,n,1), x(i,m,n,2), x(i,m,n,3), &
                      x(i,m,k,1), x(i,m,k,2), x(i,m,k,3))

                 vp2 = volpym(x(l,j,k,1), x(l,j,k,2), x(l,j,k,3), &
                      x(l,m,k,1), x(l,m,k,2), x(l,m,k,3), &
                      x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                      x(l,j,n,1), x(l,j,n,2), x(l,j,n,3))

                 vp3 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                      x(l,j,k,1), x(l,j,k,2), x(l,j,k,3), &
                      x(l,j,n,1), x(l,j,n,2), x(l,j,n,3), &
                      x(i,j,n,1), x(i,j,n,2), x(i,j,n,3))

                 vp4 = volpym(x(i,m,k,1), x(i,m,k,2), x(i,m,k,3), &
                      x(i,m,n,1), x(i,m,n,2), x(i,m,n,3), &
                      x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                      x(l,m,k,1), x(l,m,k,2), x(l,m,k,3))

                 vp5 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                      x(i,m,k,1), x(i,m,k,2), x(i,m,k,3), &
                      x(l,m,k,1), x(l,m,k,2), x(l,m,k,3), &
                      x(l,j,k,1), x(l,j,k,2), x(l,j,k,3))

                 vp6 = volpym(x(i,j,n,1), x(i,j,n,2), x(i,j,n,3), &
                      x(l,j,n,1), x(l,j,n,2), x(l,j,n,3), &
                      x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                      x(i,m,n,1), x(i,m,n,2), x(i,m,n,3))

                 ! Set the volume to 1/6 of the sum of the volumes of
                 ! the pyramid. Remember that volpym computes 6 times
                 ! the volume. For overset purposes, we can just take
                 ! the ABS for the cases of left handed blocks.

                 vol(i,j,k) = sixth*(vp1 + vp2 + vp3 + vp4 + vp5 + vp6)

                 if(vol(i,j,k) < zero) then
                    nVolNeg = nVolNeg + 1
                 else
                    nVolPos = nVolPos + 1
                 endif
              end do
           end do
        end do
        if(nVolPos == 0) then       ! Left handed block.
           flowDoms(nn,level,sps)%rightHanded = .false.
        else                        ! Right handed (or bad) block.
           flowDoms(nn,level,sps)%rightHanded = .true.
        endif

     end do domains
  end do spectral

  call determineNcellGlobal(level)
  call setGlobalCellsAndNodes(level)
  call allocMemBCData
  call oversetComm(level, .true., .false.)
contains
  function volpym(xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd)
    !
    !        ****************************************************************
    !        *                                                              *
    !        * volpym computes 6 times the volume of a pyramid. Node p,     *
    !        * whose coordinates are set in the subroutine metric itself,   *
    !        * is the top node and a-b-c-d is the quadrilateral surface.    *
    !        * It is assumed that the cross product vCa * vDb points in     *
    !        * the direction of the top node. Here vCa is the diagonal      *
    !        * running from node c to node a and vDb the diagonal from      *
    !        * node d to node b.                                            *
    !        *                                                              *
    !        ****************************************************************
    !
    use precision
    implicit none
    !
    !        Function type.
    !
    real(kind=realType) :: volpym
    !
    !        Function arguments.
    !
    real(kind=realType), intent(in) :: xa, ya, za, xb, yb, zb
    real(kind=realType), intent(in) :: xc, yc, zc, xd, yd, zd
    !
    !        ****************************************************************
    !        *                                                              *
    !        * Begin execution                                              *
    !        *                                                              *
    !        ****************************************************************
    !
    volpym = (xp - fourth*(xa + xb  + xc + xd))              &
         * ((ya - yc)*(zb - zd) - (za - zc)*(yb - yd))   + &
         (yp - fourth*(ya + yb  + yc + yd))              &
         * ((za - zc)*(xb - xd) - (xa - xc)*(zb - zd))   + &
         (zp - fourth*(za + zb  + zc + zd))              &
         * ((xa - xc)*(yb - yd) - (ya - yc)*(xb - xd))

  end function volpym
end subroutine preprocessingCustomOverset
