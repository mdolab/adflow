module coarseUtils

contains



  subroutine createCoarseBlocks(level)
    !
    !       createCoarseBlocks creates the block data structure for the    
    !       given coarse grid from the 1 level finer grid. Only direct     
    !       info is created, like owned coordinates, block sizes and       
    !       subface info. Indirect info, like face normals, volumes, wall  
    !       distances, etc. Are created later on. That info can be created 
    !       independent of the finer grid.                                 
    !
    use constants
    use block
    use inputTimeSpectral
    use coarse1to1Subface
    use coarseningInfo
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j, k, ii, jj, kk, n1to1
    integer(kind=intType) :: nn, mm, iil, jjl, kkl, il, jl, kl
    integer(kind=intType) :: iBeg, jBeg, kBeg, iEnd, jEnd, kEnd
    integer(kind=intType) :: l1, L2, l3
    integer(kind=intType) :: levm1, donorOffset, fact

    integer(kind=intType), dimension(:), allocatable :: imap, iimap
    integer(kind=intType), dimension(:), allocatable :: jmap, jjmap
    integer(kind=intType), dimension(:), allocatable :: kmap, kkmap

    integer(kind=intType), dimension(:), pointer :: dfine

    logical, dimension(:), pointer :: iCo, jCo, kCo

    ! Store the finer grid level in levm1

    levm1 = level -1

    ! Determine the total number of 1 to 1 block faces on this
    ! processor for the fine grid level.

    nSubface1to1 = 0
    do nn=1,nDom
       nSubface1to1 = nSubface1to1 + flowDoms(nn,levm1,1)%n1to1
    enddo

    ! Allocate the memory for subface1to1 and coarseInfo.

    allocate(subface1to1(nSubface1to1), coarseInfo(nDom), stat=ierr)
    if(ierr /= 0)                            &
         call terminate("createCoarseBlocks", &
         "Memory allocation failure for subface1to1 &
         &and coarseInfo")

    ! Loop over the number of blocks on this processor.

    n1to1 = 0
    domains: do nn=1,nDom

       ! If levm1 is 1, i.e. the finest grid, set the coarsenings to
       ! regular, as this is a regular grid.

       if(levm1 == 1) then
          flowDoms(nn,levm1,1)%iCoarsened = regular
          flowDoms(nn,levm1,1)%jCoarsened = regular
          flowDoms(nn,levm1,1)%kCoarsened = regular
       endif

       ! Copy cgns block id, blockIsMoving and addGridVelocities.
       ! Although identical to the values on the finest grid level,
       ! it is copied anyway for consistency reasons.

       flowDoms(nn,level,1)%cgnsBlockID   = &
            flowDoms(nn,levm1,1)%cgnsBlockID
       flowDoms(nn,level,1)%blockIsMoving = &
            flowDoms(nn,levm1,1)%blockIsMoving
       flowDoms(nn,level,1)%addGridVelocities = &
            flowDoms(nn,levm1,1)%addGridVelocities

       ! Store the number of fine nodes a bit easier and allocate the
       ! memory for the logicals iCo, jCo and kCo and for coarseIs1to1.

       iil = flowDoms(nn,levm1,1)%il
       jjl = flowDoms(nn,levm1,1)%jl
       kkl = flowDoms(nn,levm1,1)%kl

       ii = flowDoms(nn,levm1,1)%n1to1
       allocate(flowDoms(nn,levm1,1)%iCo(iil), &
            flowDoms(nn,levm1,1)%jCo(jjl), &
            flowDoms(nn,levm1,1)%kCo(kkl), &
            coarseInfo(nn)%coarseIs1to1(ii), stat=ierr)
       if(ierr /= 0)                          &
            call terminate("createCoarseBlocks", &
            "Memory allocation failure for iCo, jCo, kCo &
            &and coarseIs1to1")

       ! Initialize iCo, jCo and kCo such that the block boundaries
       ! remain, but all other nodes can disappear. Set pointers to
       ! make the code more readable.

       iCo => flowDoms(nn,levm1,1)%iCo
       jCo => flowDoms(nn,levm1,1)%jCo
       kCo => flowDoms(nn,levm1,1)%kCo

       iCo = .false.; iCo(1) = .true.; iCo(iil) = .true.
       jCo = .false.; jCo(1) = .true.; jCo(jjl) = .true.
       kCo = .false.; kCo(1) = .true.; kCo(kkl) = .true.

       ! Loop over the subfaces to keep their boundaries. Also internal
       ! block boundaries are kept.

       do i=1,flowDoms(nn,levm1,1)%nSubface
          iCo(flowDoms(nn,levm1,1)%inBeg(i)) = .true.
          iCo(flowDoms(nn,levm1,1)%inEnd(i)) = .true.

          jCo(flowDoms(nn,levm1,1)%jnBeg(i)) = .true.
          jCo(flowDoms(nn,levm1,1)%jnEnd(i)) = .true.

          kCo(flowDoms(nn,levm1,1)%knBeg(i)) = .true.
          kCo(flowDoms(nn,levm1,1)%knEnd(i)) = .true.
       enddo

       ! Create the coarser grid in i-direction. In case the fine grid
       ! was already a nonregular coarsened block, start at the
       ! opposite boundary. For a regular grid simply start at the
       ! left boundary.

       if(flowDoms(nn,levm1,1)%iCoarsened == leftStarted) then

          flowDoms(nn,level,1)%iCoarsened = rightStarted

          do i=(iil-1),2,-1
             if(.not. iCo(i+1)) iCo(i) = .true.
          enddo

       else

          flowDoms(nn,level,1)%iCoarsened = leftStarted

          do i=2,(iil-1)
             if(.not. iCo(i-1)) iCo(i) = .true.
          enddo

       endif

       ! Create the coarser grid in j-direction. Same story as in
       ! i-direction.

       if(flowDoms(nn,levm1,1)%jCoarsened == leftStarted) then

          flowDoms(nn,level,1)%jCoarsened = rightStarted

          do j=(jjl-1),2,-1
             if(.not. jCo(j+1)) jCo(j) = .true.
          enddo

       else

          flowDoms(nn,level,1)%jCoarsened = leftStarted

          do j=2,(jjl-1)
             if(.not. jCo(j-1)) jCo(j) = .true.
          enddo

       endif

       ! Create the coarser grid in k-direction. Same story as in
       ! i- and j-direction.

       if(flowDoms(nn,levm1,1)%kCoarsened == leftStarted) then

          flowDoms(nn,level,1)%kCoarsened = rightStarted

          do k=(kkl-1),2,-1
             if(.not. kCo(k+1)) kCo(k) = .true.
          enddo

       else

          flowDoms(nn,level,1)%kCoarsened = leftStarted

          do k=2,(kkl-1)
             if(.not. kCo(k-1)) kCo(k) = .true.
          enddo

       endif

       ! Determine the number of points in each direction for the
       ! coarse grid.

       il = 0
       do i=1,iil
          if( iCo(i) ) il = il +1
       enddo

       jl = 0
       do j=1,jjl
          if( jCo(j) ) jl = jl +1
       enddo

       kl = 0
       do k=1,kkl
          if( kCo(k) ) kl = kl +1
       enddo

       ! Store the number of nodes and cells in the three directions
       ! and the dimensions for the halo based quantities.

       flowDoms(nn,level,1)%il = il
       flowDoms(nn,level,1)%jl = jl
       flowDoms(nn,level,1)%kl = kl

       flowDoms(nn,level,1)%nx = il - 1
       flowDoms(nn,level,1)%ny = jl - 1
       flowDoms(nn,level,1)%nz = kl - 1

       flowDoms(nn,level,1)%ie = il + 1
       flowDoms(nn,level,1)%je = jl + 1
       flowDoms(nn,level,1)%ke = kl + 1

       flowDoms(nn,level,1)%ib = il + 2
       flowDoms(nn,level,1)%jb = jl + 2
       flowDoms(nn,level,1)%kb = kl + 2

       ! If the coarsening was regular, i.e. if the number of fine
       ! grid cells is twice the number of coarse grid cells, reset
       ! the coarsening to regular.

       if(flowDoms(nn,levm1,1)%nx == 2*flowDoms(nn,level,1)%nx) &
            flowDoms(nn,level,1)%iCoarsened = regular
       if(flowDoms(nn,levm1,1)%ny == 2*flowDoms(nn,level,1)%ny) &
            flowDoms(nn,level,1)%jCoarsened = regular
       if(flowDoms(nn,levm1,1)%nz == 2*flowDoms(nn,level,1)%nz) &
            flowDoms(nn,level,1)%kCoarsened = regular
       !
       !         The variables, which control the restriction to and the      
       !         interpolation from the coarser grid level.                   
       !
       ! Allocate the memory.

       i = flowDoms(nn,level,1)%ie
       j = flowDoms(nn,level,1)%je
       k = flowDoms(nn,level,1)%ke

       allocate(flowDoms(nn,level,1)%mgIFine(1:i,2),     &
            flowDoms(nn,level,1)%mgJFine(1:j,2),     &
            flowDoms(nn,level,1)%mgKFine(1:k,2),     &
            flowDoms(nn,level,1)%mgIWeight(2:il),    &
            flowDoms(nn,level,1)%mgJWeight(2:jl),    &
            flowDoms(nn,level,1)%mgKWeight(2:kl),    &
            flowDoms(nn,levm1,1)%mgICoarse(2:iil,2), &
            flowDoms(nn,levm1,1)%mgJCoarse(2:jjl,2), &
            flowDoms(nn,levm1,1)%mgKCoarse(2:kkl,2), stat=ierr)
       if(ierr /= 0)                            &
            call terminate("createCoarseBlocks", &
            "Memory allocation failure for interpolation &
            &variables")

       ! Set the halo values to the halo indices of the fine level.

       flowDoms(nn,level,1)%mgIFine(1,1:2) = (/ 0, 1 /)
       flowDoms(nn,level,1)%mgJFine(1,1:2) = (/ 0, 1 /)
       flowDoms(nn,level,1)%mgKFine(1,1:2) = (/ 0, 1 /)

       flowDoms(nn,level,1)%mgIFine(i,1) = flowDoms(nn,levm1,1)%ie
       flowDoms(nn,level,1)%mgIFine(i,2) = flowDoms(nn,levm1,1)%ib
       flowDoms(nn,level,1)%mgJFine(j,1) = flowDoms(nn,levm1,1)%je
       flowDoms(nn,level,1)%mgJFine(j,2) = flowDoms(nn,levm1,1)%jb
       flowDoms(nn,level,1)%mgKFine(k,1) = flowDoms(nn,levm1,1)%ke
       flowDoms(nn,level,1)%mgKFine(k,2) = flowDoms(nn,levm1,1)%kb

       ! Determine the restriction variables in i-direction.

       ii = 2
       do i=2,iil
          if( iCo(i) ) then
             if( iCo(i-1) ) then
                flowDoms(nn,level,1)%mgIFine(ii,1) = i
                flowDoms(nn,level,1)%mgIFine(ii,2) = i
                flowDoms(nn,level,1)%mgIWeight(ii) = half
             else
                flowDoms(nn,level,1)%mgIFine(ii,1) = i-1
                flowDoms(nn,level,1)%mgIFine(ii,2) = i
                flowDoms(nn,level,1)%mgIWeight(ii) = one
             endif
             ii = ii+1
          endif
       enddo

       ! Determine the restriction variables in j-direction.

       jj = 2
       do j=2,jjl
          if( jCo(j) ) then
             if( jCo(j-1) ) then
                flowDoms(nn,level,1)%mgJFine(jj,1) = j
                flowDoms(nn,level,1)%mgJFine(jj,2) = j
                flowDoms(nn,level,1)%mgJWeight(jj) = half
             else
                flowDoms(nn,level,1)%mgJFine(jj,1) = j-1
                flowDoms(nn,level,1)%mgJFine(jj,2) = j
                flowDoms(nn,level,1)%mgJWeight(jj) = one
             endif
             jj = jj+1
          endif
       enddo

       ! Determine the restriction variables in k-direction.

       kk = 2
       do k=2,kkl
          if( kCo(k) ) then
             if( kCo(k-1) ) then
                flowDoms(nn,level,1)%mgKFine(kk,1) = k
                flowDoms(nn,level,1)%mgKFine(kk,2) = k
                flowDoms(nn,level,1)%mgKWeight(kk) = half
             else
                flowDoms(nn,level,1)%mgKFine(kk,1) = k-1
                flowDoms(nn,level,1)%mgKFine(kk,2) = k
                flowDoms(nn,level,1)%mgKWeight(kk) = one
             endif
             kk = kk+1
          endif
       enddo

       ! Determine the interpolation variables in i-direction.

       ii = 2
       do i=2,iil
          if( iCo(i) ) then
             if( iCo(i-1) ) then
                flowDoms(nn,levm1,1)%mgICoarse(i,1) = ii
                flowDoms(nn,levm1,1)%mgICoarse(i,2) = ii
             else
                flowDoms(nn,levm1,1)%mgICoarse(i,1) = ii
                flowDoms(nn,levm1,1)%mgICoarse(i,2) = ii+1
             endif
             ii = ii+1
          else
             flowDoms(nn,levm1,1)%mgICoarse(i,1) = ii
             flowDoms(nn,levm1,1)%mgICoarse(i,2) = ii-1
          endif
       enddo

       ! Determine the interpolation variables in j-direction.

       jj = 2
       do j=2,jjl
          if( jCo(j) ) then
             if( jCo(j-1) ) then
                flowDoms(nn,levm1,1)%mgJCoarse(j,1) = jj
                flowDoms(nn,levm1,1)%mgJCoarse(j,2) = jj
             else
                flowDoms(nn,levm1,1)%mgJCoarse(j,1) = jj
                flowDoms(nn,levm1,1)%mgJCoarse(j,2) = jj+1
             endif
             jj = jj+1
          else
             flowDoms(nn,levm1,1)%mgJCoarse(j,1) = jj
             flowDoms(nn,levm1,1)%mgJCoarse(j,2) = jj-1
          endif
       enddo

       ! Determine the interpolation variables in k-direction.

       kk = 2
       do k=2,kkl
          if( kCo(k) ) then
             if( kCo(k-1) ) then
                flowDoms(nn,levm1,1)%mgKCoarse(k,1) = kk
                flowDoms(nn,levm1,1)%mgKCoarse(k,2) = kk
             else
                flowDoms(nn,levm1,1)%mgKCoarse(k,1) = kk
                flowDoms(nn,levm1,1)%mgKCoarse(k,2) = kk+1
             endif
             kk = kk+1
          else
             flowDoms(nn,levm1,1)%mgKCoarse(k,1) = kk
             flowDoms(nn,levm1,1)%mgKCoarse(k,2) = kk-1
          endif
       enddo
       !
       !         The coordinate mapping from fine to coarse and coarse to     
       !         fine. These are needed to determine the coarse grid subface  
       !         info.                                                        
       !
       ! Allocate the memory.

       allocate(imap(iil), jmap(jjl), kmap(kkl), &
            iimap(il), jjmap(jl), kkmap(kl), stat=ierr)
       if(ierr /= 0)                          &
            call terminate("createCoarseBlocks", &
            "Memory allocation failure for imap, etc")

       ! Set the values of imap and iimap.

       ii = 1
       do i=1,iil
          if( iCo(i) ) then
             imap(i)   = ii
             iimap(ii) = i
             ii = ii+1
          endif
       enddo

       ! Set the values of jmap and jjmap.

       jj = 1
       do j=1,jjl
          if( jCo(j) ) then
             jmap(j)   = jj
             jjmap(jj) = j
             jj = jj+1
          endif
       enddo

       ! Set the values of kmap and kkmap.

       kk = 1
       do k=1,kkl
          if( kCo(k) ) then
             kmap(k)   = kk
             kkmap(kk) = k
             kk = kk+1
          endif
       enddo
       !
       !         The subface info. Except for the subface range all other     
       !         data can be copied. The range must be adapted and the donor  
       !         range is created later, because the coarsening info of the   
       !         donor block must be known.                                   
       !
       flowDoms(nn,level,1)%nSubface   = flowDoms(nn,levm1,1)%nSubface
       flowDoms(nn,level,1)%n1to1      = flowDoms(nn,levm1,1)%n1to1
       flowDoms(nn,level,1)%nBocos     = flowDoms(nn,levm1,1)%nBocos
       flowDoms(nn,level,1)%nViscBocos = flowDoms(nn,levm1,1)%nViscBocos

       ! Allocate the memory.

       mm = flowDoms(nn,level,1)%nSubface
       allocate(flowDoms(nn,level,1)%BCType(mm),      &
            flowDoms(nn,level,1)%BCFaceID(mm),    &
            flowDoms(nn,level,1)%cgnsSubface(mm), &
            flowDoms(nn,level,1)%neighBlock(mm),  &
            flowDoms(nn,level,1)%neighProc(mm),   &
            flowDoms(nn,level,1)%groupNum(mm),    &
            flowDoms(nn,level,1)%inBeg(mm),       &
            flowDoms(nn,level,1)%jnBeg(mm),       &
            flowDoms(nn,level,1)%knBeg(mm),       &
            flowDoms(nn,level,1)%inEnd(mm),       &
            flowDoms(nn,level,1)%jnEnd(mm),       &
            flowDoms(nn,level,1)%knEnd(mm),       &
            flowDoms(nn,level,1)%dinBeg(mm),      &
            flowDoms(nn,level,1)%djnBeg(mm),      &
            flowDoms(nn,level,1)%dknBeg(mm),      &
            flowDoms(nn,level,1)%dinEnd(mm),      &
            flowDoms(nn,level,1)%djnEnd(mm),      &
            flowDoms(nn,level,1)%dknEnd(mm),      &
            flowDoms(nn,level,1)%l1(mm),          &
            flowDoms(nn,level,1)%l2(mm),          &
            flowDoms(nn,level,1)%l3(mm),          &
            stat=ierr)
       if(ierr /= 0)                          &
            call terminate("createCoarseBlocks", &
            "Memory allocation failure for subface info")

       ! Loop over the subfaces.

       subfaces: do mm=1,flowDoms(nn,level,1)%nSubface

          ! Determine the range of the coarse subface.

          flowDoms(nn,level,1)%inBeg(mm) = &
               imap(flowDoms(nn,levm1,1)%inBeg(mm))
          flowDoms(nn,level,1)%jnBeg(mm) = &
               jmap(flowDoms(nn,levm1,1)%jnBeg(mm))
          flowDoms(nn,level,1)%knBeg(mm) = &
               kmap(flowDoms(nn,levm1,1)%knBeg(mm))

          flowDoms(nn,level,1)%inEnd(mm) = &
               imap(flowDoms(nn,levm1,1)%inEnd(mm))
          flowDoms(nn,level,1)%jnEnd(mm) = &
               jmap(flowDoms(nn,levm1,1)%jnEnd(mm))
          flowDoms(nn,level,1)%knEnd(mm) = &
               kmap(flowDoms(nn,levm1,1)%knEnd(mm))

          ! Copy the rest of the subface info.

          flowDoms(nn,level,1)%BCType(mm)     = &
               flowDoms(nn,levm1,1)%BCType(mm)
          flowDoms(nn,level,1)%BCFaceID(mm)   = &
               flowDoms(nn,levm1,1)%BCFaceID(mm)
          flowDoms(nn,level,1)%neighBlock(mm) = &
               flowDoms(nn,levm1,1)%neighBlock(mm)
          flowDoms(nn,level,1)%neighProc(mm) = &
               flowDoms(nn,levm1,1)%neighProc(mm)
          flowDoms(nn,level,1)%groupNum(mm)  = &
               flowDoms(nn,levm1,1)%groupNum(mm)
          flowDoms(nn,level,1)%cgnsSubface(mm) = &
               flowDoms(nn,levm1,1)%cgnsSubface(mm)

          flowDoms(nn,level,1)%l1(mm) = flowDoms(nn,levm1,1)%l1(mm)
          flowDoms(nn,level,1)%l2(mm) = flowDoms(nn,levm1,1)%l2(mm)
          flowDoms(nn,level,1)%l3(mm) = flowDoms(nn,levm1,1)%l3(mm)

          ! Create some info if this is a 1 to 1 subface. This is stored
          ! in the array subface1to1 and is needed to determine the donor
          ! info and to check if this is still a 1 to 1 subface on the
          ! coarse grid.

          subface_1to1: if(mm > flowDoms(nn,level,1)%nBocos .and. &
               mm <= (flowDoms(nn,level,1)%nBocos     &
               +  flowDoms(nn,level,1)%n1to1)) then

             ! Update the counter.

             n1to1 = n1to1 +1

             ! Store the range of the 1 to 1 subface.

             subface1to1(n1to1)%iBeg = flowDoms(nn,level,1)%inBeg(mm)
             subface1to1(n1to1)%jBeg = flowDoms(nn,level,1)%jnBeg(mm)
             subface1to1(n1to1)%kBeg = flowDoms(nn,level,1)%knBeg(mm)
             subface1to1(n1to1)%iEnd = flowDoms(nn,level,1)%inEnd(mm)
             subface1to1(n1to1)%jEnd = flowDoms(nn,level,1)%jnEnd(mm)
             subface1to1(n1to1)%kEnd = flowDoms(nn,level,1)%knEnd(mm)

             ! Store the new range of this subface a bit easier.
             ! Make sure that i1, etc contains the lowest index.

             iBeg = min(subface1to1(n1to1)%iBeg, subface1to1(n1to1)%iEnd)
             jBeg = min(subface1to1(n1to1)%jBeg, subface1to1(n1to1)%jEnd)
             kBeg = min(subface1to1(n1to1)%kBeg, subface1to1(n1to1)%kEnd)

             iEnd = max(subface1to1(n1to1)%iBeg, subface1to1(n1to1)%iEnd)
             jEnd = max(subface1to1(n1to1)%jBeg, subface1to1(n1to1)%jEnd)
             kEnd = max(subface1to1(n1to1)%kBeg, subface1to1(n1to1)%kEnd)

             ! And copy it back into subface1to1.

             subface1to1(n1to1)%iBeg = iBeg
             subface1to1(n1to1)%jBeg = jBeg
             subface1to1(n1to1)%kBeg = kBeg
             subface1to1(n1to1)%iEnd = iEnd
             subface1to1(n1to1)%jEnd = jEnd
             subface1to1(n1to1)%kEnd = kEnd

             ! Store the processor and block ID of the donor.

             subface1to1(n1to1)%neighProc  = &
                  flowDoms(nn,level,1)%neighProc(mm)
             subface1to1(n1to1)%neighBlock = &
                  flowDoms(nn,level,1)%neighBlock(mm)

             ! Store the shorthand for the transformation matrix a
             ! bit easier.

             l1 = flowDoms(nn,level,1)%l1(mm)
             L2 = flowDoms(nn,level,1)%l2(mm)
             l3 = flowDoms(nn,level,1)%l3(mm)

             ! Determine the fine grid donor indices of the subface.
             ! i-direction.

             fact = 1
             if(l1 < 0) fact = -1
             ii = iEnd - iBeg + 1

             l1 = abs(l1)
             select case(l1)
             case (1_intType)
                allocate(subface1to1(n1to1)%idfine(ii), stat=ierr)
                dfine => subface1to1(n1to1)%idfine
                subface1to1(n1to1)%ndi = ii
                donorOffset = flowDoms(nn,levm1,1)%dinBeg(mm)
             case (2_intType)
                allocate(subface1to1(n1to1)%jdfine(ii), stat=ierr)
                dfine => subface1to1(n1to1)%jdfine
                subface1to1(n1to1)%ndj = ii
                donorOffset = flowDoms(nn,levm1,1)%djnBeg(mm)
             case (3_intType)
                allocate(subface1to1(n1to1)%kdfine(ii), stat=ierr)
                dfine => subface1to1(n1to1)%kdfine
                subface1to1(n1to1)%ndk = ii
                donorOffset = flowDoms(nn,levm1,1)%dknBeg(mm)
             end select

             if(ierr /= 0)                          &
                  call terminate("createCoarseBlocks", &
                  "Memory allocation failure for idfine")

             ii = 1
             do i=iBeg,iEnd
                jj = fact*(iimap(i) - flowDoms(nn,levm1,1)%inBeg(mm))
                dfine(ii) = jj + donorOffset
                ii = ii+1
             enddo

             ! Determine the fine grid donor indices of the subface.
             ! j-direction.

             fact = 1
             if(l2 < 0) fact = -1
             ii = jEnd - jBeg + 1

             L2 = abs(l2)
             select case(l2)
             case (1_intType)
                allocate(subface1to1(n1to1)%idfine(ii), stat=ierr)
                dfine => subface1to1(n1to1)%idfine
                subface1to1(n1to1)%ndi = ii
                donorOffset = flowDoms(nn,levm1,1)%dinBeg(mm)
             case (2_intType)
                allocate(subface1to1(n1to1)%jdfine(ii), stat=ierr)
                dfine => subface1to1(n1to1)%jdfine
                subface1to1(n1to1)%ndj = ii
                donorOffset = flowDoms(nn,levm1,1)%djnBeg(mm)
             case (3_intType)
                allocate(subface1to1(n1to1)%kdfine(ii), stat=ierr)
                dfine => subface1to1(n1to1)%kdfine
                subface1to1(n1to1)%ndk = ii
                donorOffset = flowDoms(nn,levm1,1)%dknBeg(mm)
             end select

             if(ierr /= 0)                          &
                  call terminate("createCoarseBlocks", &
                  "Memory allocation failure for jdfine")

             jj = 1
             do j=jBeg,jEnd
                ii = fact*(jjmap(j) - flowDoms(nn,levm1,1)%jnBeg(mm))
                dfine(jj) = ii + donorOffset
                jj = jj+1
             enddo

             ! Determine the fine grid donor indices of the subface.
             ! k-direction.

             fact = 1
             if(l3 < 0) fact = -1
             ii = kEnd - kBeg + 1

             l3 = abs(l3)
             select case(l3)
             case (1_intType)
                allocate(subface1to1(n1to1)%idfine(ii), stat=ierr)
                dfine => subface1to1(n1to1)%idfine
                subface1to1(n1to1)%ndi = ii
                donorOffset = flowDoms(nn,levm1,1)%dinBeg(mm)
             case (2_intType)
                allocate(subface1to1(n1to1)%jdfine(ii), stat=ierr)
                dfine => subface1to1(n1to1)%jdfine
                subface1to1(n1to1)%ndj = ii
                donorOffset = flowDoms(nn,levm1,1)%djnBeg(mm)
             case (3_intType)
                allocate(subface1to1(n1to1)%kdfine(ii), stat=ierr)
                dfine => subface1to1(n1to1)%kdfine
                subface1to1(n1to1)%ndk = ii
                donorOffset = flowDoms(nn,levm1,1)%dknBeg(mm)
             end select

             if(ierr /= 0)                          &
                  call terminate("createCoarseBlocks", &
                  "Memory allocation failure for kdfine")

             kk = 1
             do k=kBeg,kEnd
                ii = fact*(kkmap(k) - flowDoms(nn,levm1,1)%knBeg(mm))
                dfine(kk) = ii + donorOffset
                kk = kk+1
             enddo

          endif subface_1to1

       enddo subfaces

       ! Release the local memory allocated inside this loop.

       deallocate(imap, jmap, kmap, iimap, jjmap, kkmap, stat=ierr)
       if(ierr /= 0)                          &
            call terminate("createCoarseBlocks", &
            "Deallocation error for imap, iimap, etc.")

       ! Allocate the memory for the coordinates of all time spectral
       ! solutions.

       ii = flowDoms(nn,level,1)%ie
       jj = flowDoms(nn,level,1)%je
       kk = flowDoms(nn,level,1)%ke

       do mm=1,nTimeIntervalsSpectral
          allocate(flowDoms(nn,level,mm)%x(0:ii,0:jj,0:kk,3), stat=ierr)
          if(ierr /= 0)                          &
               call terminate("createCoarseBlocks", &
               "Memory allocation failure for x")
       enddo

       ! Copy the scalars such that they are known for all time
       ! spectral solutions.

       do mm=2,nTimeIntervalsSpectral

          flowDoms(nn,level,mm)%nx = flowDoms(nn,level,1)%nx
          flowDoms(nn,level,mm)%ny = flowDoms(nn,level,1)%ny
          flowDoms(nn,level,mm)%nz = flowDoms(nn,level,1)%nz

          flowDoms(nn,level,mm)%il = flowDoms(nn,level,1)%il
          flowDoms(nn,level,mm)%jl = flowDoms(nn,level,1)%jl
          flowDoms(nn,level,mm)%kl = flowDoms(nn,level,1)%kl

          flowDoms(nn,level,mm)%ie = flowDoms(nn,level,1)%ie
          flowDoms(nn,level,mm)%je = flowDoms(nn,level,1)%je
          flowDoms(nn,level,mm)%ke = flowDoms(nn,level,1)%ke

          flowDoms(nn,level,mm)%ib = flowDoms(nn,level,1)%ib
          flowDoms(nn,level,mm)%jb = flowDoms(nn,level,1)%jb
          flowDoms(nn,level,mm)%kb = flowDoms(nn,level,1)%kb

          flowDoms(nn,level,mm)%nSubface   = flowDoms(nn,level,1)%nSubface
          flowDoms(nn,level,mm)%n1to1      = flowDoms(nn,level,1)%n1to1
          flowDoms(nn,level,mm)%nBocos     = flowDoms(nn,level,1)%nBocos
          flowDoms(nn,level,mm)%nViscBocos = flowDoms(nn,level,1)%nViscBocos

          flowDoms(nn,level,mm)%iCoarsened = flowDoms(nn,level,1)%iCoarsened
          flowDoms(nn,level,mm)%jCoarsened = flowDoms(nn,level,1)%jCoarsened
          flowDoms(nn,level,mm)%kCoarsened = flowDoms(nn,level,1)%kCoarsened

          flowDoms(nn,level,mm)%blockIsMoving = &
               flowDoms(nn,level,1)%blockIsMoving
          flowDoms(nn,level,mm)%cgnsBlockID   = &
               flowDoms(nn,level,1)%cgnsBlockID

          flowDoms(nn,level,mm)%addGridVelocities = &
               flowDoms(nn,level,1)%addGridVelocities

       enddo

    enddo domains

    ! Determine the owned coordinates.

    call coarseOwnedCoordinates(level)

    ! Determine the donor info for the internal block boundaries.

    call coarseDonorInfo(level)

    ! Remove possible coarse grid non-matching block boundaries from
    ! the 1 to 1 list.

    call checkCoarse1to1(level)

    ! Release the memory of the coarsening info.

    do nn=1,nDom
       deallocate(coarseInfo(nn)%coarseIs1to1, stat=ierr)
       if(ierr /= 0)                          &
            call terminate("createCoarseBlocks", &
            "Deallocation error for coarseIs1to1")
    enddo

    deallocate(coarseInfo, stat=ierr)
    if(ierr /= 0)                          &
         call terminate("createCoarseBlocks", &
         "Deallocation error for coarseInfo")

  end subroutine createCoarseBlocks

  !      ==================================================================

  subroutine coarseOwnedCoordinates(level)
    !
    !       coarseOwnedCoordinates determines from the coarsening info     
    !       the owned coordinates of the coarse grid. This is done in a    
    !       separate routine, because in unsteady moving mesh mode or for  
    !       deforming meshes only new coordinates need to be computed,     
    !       while the connectivity remains the same.                       
    !
    use block
    use inputTimeSpectral
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ii, jj, kk
    integer(kind=intType) :: nn, mm, il, jl, kl, levm1

    logical, dimension(:), pointer :: iCo, jCo, kCo

    ! Store the finer grid level in levm1

    levm1 = level -1

    ! Loop over the number of blocks on this processor

    domains: do nn=1,nDom

       ! Easier storage of some variables of this block.

       il = flowDoms(nn,levm1,1)%il
       jl = flowDoms(nn,levm1,1)%jl
       kl = flowDoms(nn,levm1,1)%kl

       iCo => flowDoms(nn,levm1,1)%iCo
       jCo => flowDoms(nn,levm1,1)%jCo
       kCo => flowDoms(nn,levm1,1)%kCo

       ! Loop over the fine grid lines in the three directions and
       ! determine which should be kept on the coarse grid.

       kk = 1
       do k=1,kl
          if( kCo(k) ) then
             jj = 1
             do j=1,jl
                if( jCo(j) ) then
                   ii = 1
                   do i=1,il
                      if( iCo(i) ) then

                         ! Loop over the spectral solutions and copy the
                         ! coordinates from the fine grid.

                         do mm=1,nTimeIntervalsSpectral
                            flowDoms(nn,level,mm)%x(ii,jj,kk,1) = &
                                 flowDoms(nn,levm1,mm)%x(i,j,k,1)
                            flowDoms(nn,level,mm)%x(ii,jj,kk,2) = &
                                 flowDoms(nn,levm1,mm)%x(i,j,k,2)
                            flowDoms(nn,level,mm)%x(ii,jj,kk,3) = &
                                 flowDoms(nn,levm1,mm)%x(i,j,k,3)
                         enddo

                         ii = ii + 1
                      endif
                   enddo
                   jj = jj + 1
                endif
             enddo
             kk = kk + 1
          endif
       enddo

    enddo domains

  end subroutine coarseOwnedCoordinates

  subroutine update1to1Coarse(level, subface)
    !
    !       update1to1Coarse determines whether or not the given 1 to 1    
    !       block boundary subface on the fine level is still a 1 to 1     
    !       subface on the coarser grid level.                             
    !
    use block
    use coarse1to1Subface
    use coarseningInfo
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType),        intent(in) :: level
    type(coarse1to1SubfaceType), intent(in) :: subface
    !
    !      Local variables.
    !
    integer(kind=intType) :: levm1, nn, mm, ii, jj, kk, ll
    integer(kind=intType) :: i, j, k
    integer(kind=intType) :: iBeg, jBeg, kBeg, iEnd, jEnd, kEnd

    logical :: subfaceFound, idir, jdir, kdir

    ! Store the finer grid level in levm1

    levm1 = level -1

    ! Store the local block id a bit easier and store the minimum and
    ! maximum index in the three coordinate directions for the fine
    ! grid.

    nn = subface%neighBlock

    mm = subface%ndi
    iBeg = min(subface%idfine(1), subface%idfine(mm))
    iEnd = max(subface%idfine(1), subface%idfine(mm))

    mm = subface%ndj
    jBeg = min(subface%jdfine(1), subface%jdfine(mm))
    jEnd = max(subface%jdfine(1), subface%jdfine(mm))

    mm = subface%ndk
    kBeg = min(subface%kdfine(1), subface%kdfine(mm))
    kEnd = max(subface%kdfine(1), subface%kdfine(mm))

    ! Find the subface id of the block, which corresponds to the
    ! given subface of the subroutine header.

    subfaceFound = .false.
    findSubface: do ii=1,flowDoms(nn,levm1,1)%n1to1

       mm = ii + flowDoms(nn,levm1,1)%nBocos

       ! Check if this subface is the correct one.

       kk = min(flowDoms(nn,levm1,1)%inBeg(mm), &
            flowDoms(nn,levm1,1)%inEnd(mm))
       ll = max(flowDoms(nn,levm1,1)%inBeg(mm), &
            flowDoms(nn,levm1,1)%inEnd(mm))

       if(kk == iBeg .and. ll == iEnd) then
          kk = min(flowDoms(nn,levm1,1)%jnBeg(mm), &
               flowDoms(nn,levm1,1)%jnEnd(mm))
          ll = max(flowDoms(nn,levm1,1)%jnBeg(mm), &
               flowDoms(nn,levm1,1)%jnEnd(mm))

          if(kk == jBeg .and. ll == jEnd) then
             kk = min(flowDoms(nn,levm1,1)%knBeg(mm), &
                  flowDoms(nn,levm1,1)%knEnd(mm))
             ll = max(flowDoms(nn,levm1,1)%knBeg(mm), &
                  flowDoms(nn,levm1,1)%knEnd(mm))

             if(kk == kBeg .and. ll == kEnd) subfaceFound = .true.
          endif
       endif

       ! Exit the loop if this is indeed the subface.

       if( subfaceFound ) exit
    enddo findSubface

    ! The subface must be found on the fine grid. Check this.

    if(.not. subfaceFound)               &
         call terminate("update1to1Coarse", &
         "Invalid fine grid 1 to 1 subface connectivity")

    ! Check the i-direction of the coarse grid subface to see if it
    ! is still a 1 to 1 subface. First check the number of grid lines.
    ! If these are identical, check each coarse grid line.

    idir = .true.
    ii = abs(flowDoms(nn,level,1)%inEnd(mm) &
         -     flowDoms(nn,level,1)%inBeg(mm)) + 1

    if(ii == subface%ndi) then
       do i=1,subface%ndi
          ii = subface%idfine(i)
          if(.not. flowDoms(nn,levm1,1)%ico(ii) ) idir = .false.
       enddo
    else
       idir = .false.
    endif

    ! Check the j-direction of the coarse grid subface to see if it
    ! is still a 1 to 1 subface. First check the number of grid lines.
    ! If these are identical, check each coarse grid line.

    jdir = .true.
    jj = abs(flowDoms(nn,level,1)%jnEnd(mm) &
         -     flowDoms(nn,level,1)%jnBeg(mm)) + 1

    if(jj == subface%ndj) then
       do j=1,subface%ndj
          jj = subface%jdfine(j)
          if(.not. flowDoms(nn,levm1,1)%jco(jj) ) jdir = .false.
       enddo
    else
       jdir = .false.
    endif

    ! Check the k-direction of the coarse grid subface to see if it
    ! is still a 1 to 1 subface. First check the number of grid lines.
    ! If these are identical, check each coarse grid line.

    kdir = .true.
    kk = abs(flowDoms(nn,level,1)%knEnd(mm) &
         -     flowDoms(nn,level,1)%knBeg(mm)) + 1

    if(kk == subface%ndk) then
       do k=1,subface%ndk
          kk = subface%kdfine(k)
          if(.not. flowDoms(nn,levm1,1)%kco(kk) ) kdir = .false.
       enddo
    else
       kdir = .false.
    endif

    ! Set coarseIs1to1 to .true. if the subface on the coarse grid
    ! is still a 1 to 1 subface. Otherwise set it to .false.

    ii = mm - flowDoms(nn,levm1,1)%nBocos
    if(idir .and. jdir .and. kdir) then
       coarseInfo(nn)%coarseIs1to1(ii) = .true.
    else
       coarseInfo(nn)%coarseIs1to1(ii) = .false.
    endif

    ! Store the donor range. It is possible that the lower and
    ! upper boundary must be reversed.

    if(flowDoms(nn,levm1,1)%dinBeg(mm) < &
         flowDoms(nn,levm1,1)%dinEnd(mm)) then
       flowDoms(nn,level,1)%dinBeg(mm) = subface%iBeg
       flowDoms(nn,level,1)%dinEnd(mm) = subface%iEnd
    else
       flowDoms(nn,level,1)%dinBeg(mm) = subface%iEnd
       flowDoms(nn,level,1)%dinEnd(mm) = subface%iBeg
    endif

    if(flowDoms(nn,levm1,1)%djnBeg(mm) < &
         flowDoms(nn,levm1,1)%djnEnd(mm)) then
       flowDoms(nn,level,1)%djnBeg(mm) = subface%jBeg
       flowDoms(nn,level,1)%djnEnd(mm) = subface%jEnd
    else
       flowDoms(nn,level,1)%djnBeg(mm) = subface%jEnd
       flowDoms(nn,level,1)%djnEnd(mm) = subface%jBeg
    endif

    if(flowDoms(nn,levm1,1)%dknBeg(mm) < &
         flowDoms(nn,levm1,1)%dknEnd(mm)) then
       flowDoms(nn,level,1)%dknBeg(mm) = subface%kBeg
       flowDoms(nn,level,1)%dknEnd(mm) = subface%kEnd
    else
       flowDoms(nn,level,1)%dknBeg(mm) = subface%kEnd
       flowDoms(nn,level,1)%dknEnd(mm) = subface%kBeg
    endif

  end subroutine update1to1Coarse
  subroutine coarseDonorInfo(level)
    !
    !       coarseDonorInfo creates the donor info for the internal        
    !       block boundaries on the given coarse grid level.               
    !
    use communication
    use coarse1to1Subface
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer :: proc, size, ierr

    integer, dimension(mpi_status_size) :: mpiStatus
    integer, dimension(nProc)           :: sizeMessage

    integer(kind=intType) :: nn, mm, ii, i
    integer(kind=intType) :: nMessageReceive, nMessageSend

    integer(kind=intType), dimension(0:nProc) :: size2proc
    integer(kind=intType), dimension(nProc)   :: tmp

    integer(kind=intType), dimension(:), allocatable :: sendBuf
    integer(kind=intType), dimension(:), allocatable, target :: recvBuf

    type(coarse1to1SubfaceType) :: tmpSubface

    ! Determine the number of integers sent to every processor.

    size2proc = 0
    do nn=1,nSubface1to1
       mm = subface1to1(nn)%neighProc + 1    ! Proc ID's start at 0.

       size2proc(mm) = size2proc(mm) + 13                          &
            + subface1to1(nn)%iEnd - subface1to1(nn)%iBeg &
            + subface1to1(nn)%jEnd - subface1to1(nn)%jBeg &
            + subface1to1(nn)%kEnd - subface1to1(nn)%kBeg
    enddo

    ! No message is sent to myself, so set the corresponding entry in
    ! size2proc to 0. The processor id's start at 0, which explains
    ! the +1.

    size2proc(myID+1) = 0

    ! Determine the number of messages i have to sent. Store in tmp
    ! a 0 or a 1, depending whether or not a message must be sent
    ! to the corresponding processor. This info is needed to determine
    ! the number of messages i will receive.

    nMessageSend = 0
    do nn=1,nProc
       if(size2proc(nn) > 0) then
          nMessageSend = nMessageSend +1
          tmp(nn) = 1
       else
          tmp(nn) = 0
       endif
    enddo

    ! Determine the number of messages i will receive.

    sizeMessage = 1
    call mpi_reduce_scatter(tmp, nMessageReceive, sizeMessage,      &
         adflow_integer, mpi_sum, ADflow_comm_world, &
         ierr)

    ! Put size2proc in cumulative storage format and store the
    ! starting entry in tmp, which is used as a counter.

    do nn=1,nProc
       size2proc(nn) = size2proc(nn) + size2proc(nn-1)
       tmp(nn)       = size2proc(nn-1)
    enddo

    ! Allocate the memory for the send buffer.

    allocate(sendBuf(size2proc(nProc)), stat=ierr)
    if(ierr /= 0)                       &
         call terminate("coarseDonorInfo", &
         "Memory allocation failure for sendBuf")

    ! Loop over the number of 1 to 1 subfaces to fill the send buffer.

    unOwnedSubfaces: do nn=1,nSubface1to1

       ! Only info that must be sent to other processors must be stored.

       if(subface1to1(nn)%neighProc /= myID) then

          ! Store the entry in tmp in mm.
          ! Note that the proc id's start at 0.

          mm = subface1to1(nn)%neighProc + 1

          ! Store the block id of the donor and the coarse grid range of
          ! the current subface in the send buffer. Update the counter
          ! tmp(mm) accordingly.

          tmp(mm) = tmp(mm) +1
          sendBuf(tmp(mm)) = subface1to1(nn)%neighBlock

          tmp(mm) = tmp(mm) +1
          sendBuf(tmp(mm)) = subface1to1(nn)%iBeg

          tmp(mm) = tmp(mm) +1
          sendBuf(tmp(mm)) = subface1to1(nn)%jBeg

          tmp(mm) = tmp(mm) +1
          sendBuf(tmp(mm)) = subface1to1(nn)%kBeg

          tmp(mm) = tmp(mm) +1
          sendBuf(tmp(mm)) = subface1to1(nn)%iEnd

          tmp(mm) = tmp(mm) +1
          sendBuf(tmp(mm)) = subface1to1(nn)%jEnd

          tmp(mm) = tmp(mm) +1
          sendBuf(tmp(mm)) = subface1to1(nn)%kEnd

          ! Store the number of points in the three coordinate
          ! directions of the coarse grid donor face.

          tmp(mm) = tmp(mm) +1
          sendBuf(tmp(mm)) = subface1to1(nn)%ndi

          tmp(mm) = tmp(mm) +1
          sendBuf(tmp(mm)) = subface1to1(nn)%ndj

          tmp(mm) = tmp(mm) +1
          sendBuf(tmp(mm)) = subface1to1(nn)%ndk

          ! Store the i-donor indices of the fine grid in sendBuf.

          do i=1,subface1to1(nn)%ndi
             tmp(mm) = tmp(mm) +1
             sendBuf(tmp(mm)) = subface1to1(nn)%idfine(i)
          enddo

          ! Store the j-donor indices of the fine grid in sendBuf.

          do i=1,subface1to1(nn)%ndj
             tmp(mm) = tmp(mm) +1
             sendBuf(tmp(mm)) = subface1to1(nn)%jdfine(i)
          enddo

          ! Store the k-donor indices of the fine grid in sendBuf.

          do i=1,subface1to1(nn)%ndk
             tmp(mm) = tmp(mm) +1
             sendBuf(tmp(mm)) = subface1to1(nn)%kdfine(i)
          enddo

       endif

    enddo unOwnedSubfaces

    ! Send the data i have to send.

    mm = 0
    sends: do nn=1,nProc

       if(size2proc(nn) > size2proc(nn-1)) then

          ! Update mm and store the processor id, the size of the
          ! message and the starting index in sendbuf.

          mm   = mm +1
          proc = nn -1
          size = size2proc(nn) - size2proc(nn-1)
          ii   = size2proc(nn-1) +1

          ! Send the message.

          call mpi_isend(sendBuf(ii), size, adflow_integer, proc, proc, &
               ADflow_comm_world, sendRequests(mm), ierr)

       endif

    enddo sends

    ! Loop over the number of 1 to 1 subfaces stored on this processor
    ! and update the locally stored info.

    do nn=1,nSubface1to1
       if(subface1to1(nn)%neighProc == myID) &
            call update1to1Coarse(level, subface1to1(nn))
    enddo

    ! Release the memory of subface1to1.

    do nn=1,nSubface1to1
       deallocate(subface1to1(nn)%idfine, subface1to1(nn)%jdfine, &
            subface1to1(nn)%kdfine, stat=ierr)
       if(ierr /= 0)                       &
            call terminate("coarseDonorInfo", &
            "Deallocation error for idfine, etc")
    enddo

    deallocate(subface1to1, stat=ierr)
    if(ierr /= 0)                       &
         call terminate("coarseDonorInfo", &
         "Deallocation error for subface1to1")

    ! Loop over the number of messages i must receive to determine
    ! info from externally stored neighboring subfaces.

    receives: do nn=1,nMessageReceive

       ! Block until a message arrives.

       call mpi_probe(mpi_any_source, myID, ADflow_comm_world, &
            mpiStatus, ierr)

       ! Find the source and size of the message.

       proc = mpiStatus(mpi_source)
       call mpi_get_count(mpiStatus, adflow_integer, size, ierr)

       ! Check in debug mode that the incoming message is of
       ! correct size.

       if( debug ) then
          if(size == mpi_undefined)             &
               call terminate("coarseDonorInfo", &
               "Unexpected size of message")
       endif

       ! Allocate the memory for the receive buffer.

       allocate(recvBuf(size), stat=ierr)
       if(ierr /= 0)                         &
            call terminate("coarseDonorInfo", &
            "Memory allocation failure for recvBuf")

       ! Receive the messsage. As it has already arrived a blocking
       ! receive can be used.

       call mpi_recv(recvBuf, size, adflow_integer, proc, myID, &
            ADflow_comm_world, mpiStatus, ierr)

       ! Loop to extract the 1 to 1 subface info from the buffer.

       ii = 1
       extractSubface: do
          ! Exit the loop if no more info is present.

          if(ii > size) exit

          ! Store the info of this subface in tmpSubface.

          tmpSubface%neighProc  = proc
          tmpSubface%neighBlock = recvBuf(ii); ii = ii +1

          tmpSubface%iBeg = recvBuf(ii); ii = ii +1
          tmpSubface%jBeg = recvBuf(ii); ii = ii +1
          tmpSubface%kBeg = recvBuf(ii); ii = ii +1

          tmpSubface%iEnd = recvBuf(ii); ii = ii +1
          tmpSubface%jEnd = recvBuf(ii); ii = ii +1
          tmpSubface%kEnd = recvBuf(ii); ii = ii +1

          tmpSubface%ndi = recvBuf(ii); ii = ii +1
          tmpSubface%ndj = recvBuf(ii); ii = ii +1
          tmpSubface%ndk = recvBuf(ii); ii = ii +1

          mm = ii + tmpSubface%ndi
          tmpSubface%idfine => recvBuf(ii:(mm-1))
          ii = mm

          mm = ii + tmpSubface%ndj
          tmpSubface%jdfine => recvBuf(ii:(mm-1))
          ii = mm

          mm = ii + tmpSubface%ndk
          tmpSubface%kdfine => recvBuf(ii:(mm-1))
          ii = mm

          ! Update the 1 to 1 subface info.

          call update1to1Coarse(level, tmpSubface)

       enddo extractSubface

       ! Release the memory of the receive buffer.

       deallocate(recvBuf, stat=ierr)
       if(ierr /= 0)                       &
            call terminate("coarseDonorInfo", &
            "Deallocation failure for recvBuf")

    enddo receives

    ! Complete the nonblocking sends.

    size = nMessageSend
    do nn=1,nMessageSend
       call mpi_waitany(size, sendRequests, proc, mpiStatus, ierr)
    enddo

    ! Release the memory of the send buffer.

    deallocate(sendBuf, stat=ierr)
    if(ierr /= 0)                       &
         call terminate("coarseDonorInfo", &
         "Deallocation failure for sendBuf")

    ! Synchronize the processors, because wild cards have been used.

    call mpi_barrier(ADflow_comm_world, ierr)

  end subroutine coarseDonorInfo

  subroutine checkCoarse1to1(level)
    !
    !       checkCoarse1to1 removes the nonmatching block boundaries       
    !       from the list of 1 to 1 matching ones. They are in there,      
    !       because they are 1 to 1 matching on the finer grids.           
    !
    use constants
    use cgnsGrid
    use block
    use inputTimeSpectral
    use coarseningInfo
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    character(len=maxStringLen) :: errorMessage

    integer(kind=intType) :: i, nn, mm, kk, ll

    ! Loop over the number of domains.

    domains: do nn=1,nDom

       ! Loop over the number of 1 to 1 subfaces. As this number can
       ! change during the loop, the control statement is done via
       ! an exit.

       i  = 1
       n1to1: do

          ! Exit the loop if the counter is larger than the number
          ! of 1 to 1 subfaces.

          if(i > flowDoms(nn,level,1)%n1to1) exit

          ! Add the offset to i to store the correct place in the
          ! subface info.

          mm = i + flowDoms(nn,level,1)%nBocos

          ! Test if this is still a 1 to 1 subface on the coarse grid.

          is1to1: if( coarseInfo(nn)%coarseIs1to1(i) ) then

             ! Subface is still a 1 to 1 subface on the coarse grid.
             ! Only the counter i must be updated.

             i = i + 1

          else is1to1

             ! Due to the coarsening the subface is not a 1 to 1 block
             ! boundary anymore. Swap the current entry with the entry
             ! kk.

             kk = flowDoms(nn,level,1)%nBocos + flowDoms(nn,level,1)%n1to1

             ! The range info.

             ll                             = flowDoms(nn,level,1)%inBeg(mm)
             flowDoms(nn,level,1)%inBeg(mm) = flowDoms(nn,level,1)%inBeg(kk)
             flowDoms(nn,level,1)%inBeg(kk) = ll

             ll                             = flowDoms(nn,level,1)%jnBeg(mm)
             flowDoms(nn,level,1)%jnBeg(mm) = flowDoms(nn,level,1)%jnBeg(kk)
             flowDoms(nn,level,1)%jnBeg(kk) = ll

             ll                             = flowDoms(nn,level,1)%knBeg(mm)
             flowDoms(nn,level,1)%knBeg(mm) = flowDoms(nn,level,1)%knBeg(kk)
             flowDoms(nn,level,1)%knBeg(kk) = ll

             ll                             = flowDoms(nn,level,1)%inEnd(mm)
             flowDoms(nn,level,1)%inEnd(mm) = flowDoms(nn,level,1)%inEnd(kk)
             flowDoms(nn,level,1)%inEnd(kk) = ll

             ll                             = flowDoms(nn,level,1)%jnEnd(mm)
             flowDoms(nn,level,1)%jnEnd(mm) = flowDoms(nn,level,1)%jnEnd(kk)
             flowDoms(nn,level,1)%jnEnd(kk) = ll

             ll                             = flowDoms(nn,level,1)%knEnd(mm)
             flowDoms(nn,level,1)%knEnd(mm) = flowDoms(nn,level,1)%knEnd(kk)
             flowDoms(nn,level,1)%knEnd(kk) = ll

             ! The donor info.

             ll                              = flowDoms(nn,level,1)%dinBeg(mm)
             flowDoms(nn,level,1)%dinBeg(mm) = flowDoms(nn,level,1)%dinBeg(kk)
             flowDoms(nn,level,1)%dinBeg(kk) = ll

             ll                              = flowDoms(nn,level,1)%djnBeg(mm)
             flowDoms(nn,level,1)%djnBeg(mm) = flowDoms(nn,level,1)%djnBeg(kk)
             flowDoms(nn,level,1)%djnBeg(kk) = ll

             ll                              = flowDoms(nn,level,1)%dknBeg(mm)
             flowDoms(nn,level,1)%dknBeg(mm) = flowDoms(nn,level,1)%dknBeg(kk)
             flowDoms(nn,level,1)%dknBeg(kk) = ll

             ll                              = flowDoms(nn,level,1)%dinEnd(mm)
             flowDoms(nn,level,1)%dinEnd(mm) = flowDoms(nn,level,1)%dinEnd(kk)
             flowDoms(nn,level,1)%dinEnd(kk) = ll

             ll                              = flowDoms(nn,level,1)%djnEnd(mm)
             flowDoms(nn,level,1)%djnEnd(mm) = flowDoms(nn,level,1)%djnEnd(kk)
             flowDoms(nn,level,1)%djnEnd(kk) = ll

             ll                              = flowDoms(nn,level,1)%dknEnd(mm)
             flowDoms(nn,level,1)%dknEnd(mm) = flowDoms(nn,level,1)%dknEnd(kk)
             flowDoms(nn,level,1)%dknEnd(kk) = ll

             ! The transformation matrix.

             ll                          = flowDoms(nn,level,1)%l1(mm)
             flowDoms(nn,level,1)%l1(mm) = flowDoms(nn,level,1)%l1(kk)
             flowDoms(nn,level,1)%l1(kk) = ll

             ll                          = flowDoms(nn,level,1)%l2(mm)
             flowDoms(nn,level,1)%l2(mm) = flowDoms(nn,level,1)%l2(kk)
             flowDoms(nn,level,1)%l2(kk) = ll

             ll                          = flowDoms(nn,level,1)%l3(mm)
             flowDoms(nn,level,1)%l3(mm) = flowDoms(nn,level,1)%l3(kk)
             flowDoms(nn,level,1)%l3(kk) = ll

             ! The rest of the subface info.

             ll = flowDoms(nn,level,1)%BCType(mm)
             flowDoms(nn,level,1)%BCType(mm) = &
                  flowDoms(nn,level,1)%BCType(kk)
             flowDoms(nn,level,1)%BCType(kk) = ll

             ll = flowDoms(nn,level,1)%BCFaceID(mm)
             flowDoms(nn,level,1)%BCFaceID(mm) = &
                  flowDoms(nn,level,1)%BCFaceID(kk)
             flowDoms(nn,level,1)%BCFaceID(kk) = ll

             ll = flowDoms(nn,level,1)%neighBlock(mm)
             flowDoms(nn,level,1)%neighBlock(mm) = &
                  flowDoms(nn,level,1)%neighBlock(kk)
             flowDoms(nn,level,1)%neighBlock(kk) = ll

             ll = flowDoms(nn,level,1)%neighProc(mm)
             flowDoms(nn,level,1)%neighProc(mm) = &
                  flowDoms(nn,level,1)%neighProc(kk)
             flowDoms(nn,level,1)%neighProc(kk) = ll

             ll = flowDoms(nn,level,1)%groupNum(mm)
             flowDoms(nn,level,1)%groupNum(mm) = &
                  flowDoms(nn,level,1)%groupNum(kk)
             flowDoms(nn,level,1)%groupNum(kk) = ll

             ll = flowDoms(nn,level,1)%cgnsSubface(mm)
             flowDoms(nn,level,1)%cgnsSubface(mm) = &
                  flowDoms(nn,level,1)%cgnsSubface(kk)
             flowDoms(nn,level,1)%cgnsSubface(kk) = ll

             ! Decrease the number of 1 to 1 block boundaries. Note the
             ! counter i should not be updated.

             flowDoms(nn,level,1)%n1to1 = flowDoms(nn,level,1)%n1to1 - 1

101          format("Non-matching block-to-block face on zone ",a, &
                  "...Support not implemented yet.")
             ll = flowDoms(nn,1,1)%cgnsBlockID
             write(errorMessage,101) cgnsDoms(ll)%zoneName
             call terminate("checkCoarse1to1", errorMessage)

          endif is1to1

       enddo n1to1

       ! Copy the number of internal subfaces for the rest of
       ! the spectral solutions.

       do i=2,nTimeIntervalsSpectral
          flowDoms(nn,level,i)%n1to1 = flowDoms(nn,level,1)%n1to1
       enddo

    enddo domains

  end subroutine checkCoarse1to1

end module coarseUtils
