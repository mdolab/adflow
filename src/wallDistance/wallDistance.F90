module wallDistance

  use constants, only : intType, realType
  use wallDistanceData
  implicit none
  save

#ifndef USE_TAPENADE
  ! nquadVisc:             Number of local quads on the viscous
  !                        bodies.
  ! nNodeVisc:             Number of local nodes on the viscous
  !                        bodies.
  ! nquadViscGlob:         Global number of viscous quads.
  ! connVisc(4,nquadVisc): Connectivity of the local viscous
  !                        quads.
  ! coorVisc(3,nNodeVisc): The coordinates of the local nodes.

  integer(kind=intType) :: nquadVisc, nNodeVisc
  integer(kind=intType) :: nquadViscGlob

  integer(kind=intType), dimension(:,:), allocatable :: connVisc
  real(kind=realType),   dimension(:,:), allocatable :: coorVisc

  ! rotMatrixSections(nSections,3,3): Rotation matrices needed
  !                                   for the alignment of the
  !                                   sections. The rotation
  !                                   matrix is a**n, where a is
  !                                   periodic transformation
  !                                   matrix; n is an integer.

  real(kind=realType), dimension(:,:,:), allocatable :: rotMatrixSections
#endif

contains

  subroutine updateWallDistancesQuickly(nn, level, sps)

    ! This is the actual update routine that uses xSurf. It is done on
    ! block-level-sps basis.  This is the used to update the wall
    ! distance. Most importantly, this routine is included in the
    ! reverse mode AD routines, but NOT the forward mode. Since it is
    ! done on a per-block basis, it is assumed that the required block
    ! pointers are already set. 

    use constants
    use blockPointers, only : nx, ny, nz, il, jl, kl, x, flowDoms, d2wall
    implicit none

    ! Subroutine arguments
    integer(kind=intType) :: nn, level, sps

    ! Local Variables
    integer(kind=intType) :: i, j, k, ii, ind(4)
    real(kind=realType) :: xp(3), xc(3), u, v

#ifdef TAPENADE_FAST
    !$AD II-LOOP
    do ii=0,nx*ny*nz-1
       i = mod(ii, nx) + 2
       j = mod(ii/nx, ny) + 2
       k = ii/(nx*ny) + 2
#else
       do k=2,kl
          do j=2,jl
             do i=2,il
#endif

                if (flowDoms(nn, level, sps)%surfNodeIndices(1, i, j, k) == 0) then 
                   ! This node is too far away and has no
                   ! association. Set the distance to a large constant. 
                   d2wall(i, j, k) = large
                   cycle
                end if

                ! Extract elemID and u-v position for the association of
                ! this cell:

                ind = flowDoms(nn,level,sps)%surfNodeIndices(:, i, j, k)
                u   = flowDoms(nn,level,sps)%uv(1,i,j,k)
                v    = flowDoms(nn,level,sps)%uv(2,i,j,k)

                ! Now we have the 4 corners, use bi-linear shape
                ! functions o to get target: (CCW ordering remember!)

                xp(:) = &
                     (one-u)*(one-v)*xSurf(3*(ind(1)-1)+1:3*ind(1)) + &
                     (    u)*(one-v)*xSurf(3*(ind(2)-1)+1:3*ind(2)) + &
                     (    u)*(    v)*xSurf(3*(ind(3)-1)+1:3*ind(3)) + &
                     (one-u)*(    v)*xSurf(3*(ind(4)-1)+1:3*ind(4)) 

                ! Get the cell center
                xc(1) = eighth*(x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
                     +         x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
                     +         x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                     +         x(i-1,j,  k,  1) + x(i,j,  k,  1))

                xc(2) = eighth*(x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
                     +         x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
                     +         x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                     +         x(i-1,j,  k,  2) + x(i,j,  k,  2))

                xc(3) = eighth*(x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
                     +         x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
                     +         x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                     +         x(i-1,j,  k,  3) + x(i,j,  k,  3))


                ! Now we have the two points...just take the norm of the
                ! distance between them

                d2wall(i,j,k) = sqrt(&
                     (xc(1)-xp(1))**2 + (xc(2)-xp(2))**2 + (xc(3)-xp(3))**2)
#ifdef TAPENADE_FAST
             end do
#else
          enddo
       enddo
    enddo
#endif

  end subroutine updateWallDistancesQuickly

  ! ----------------------------------------------------------------------
  !                                                                      |
  !                    No Tapenade Routine below this line               |
  !                                                                      |
  ! ----------------------------------------------------------------------

#ifndef USE_TAPENADE
  subroutine computeWallDistance(level, allocMem)
    !
    !       wallDistance computes the distances of the cell centers to     
    !       the nearest viscous wall. An adt type of method is used, which 
    !       guarantees to find the minimum distance to the wall. Possible  
    !       periodic transformations are taken into account, such that     
    !       also in case of a periodic problem the correct distance is     
    !       computed; the nearest wall point may lie in a periodic domain. 
    use constants
    use blockPointers, only : nDom
    use communication, only : sendBuffer, recvBuffer, myid, sumb_comm_world, &
         sendBufferSize, recvBufferSize
    use inputPhysics, only : equations, wallDistanceNeeded
    use inputTimeSpectral, only :nTimeIntervalsSpectral
    use inputDiscretization, only: useApproxWallDistance
    use utils, only : setPointers, EChk, terminate, &
         deallocateTempMemory, allocateTempMemory
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    logical, intent(in)               :: allocMem
    !
    !      Local variables.
    !
    integer :: ierr, i, j, k, nn, ii, l

    integer(kind=intType) :: sps, sps2, ll, nLevels
    logical :: tempLogical
    double precision :: t0
    character(len=3) :: integerString

    ! Check if the RANS equations are solved. If not, the wall
    ! distance is not needed and a return can be made.

    if(equations /= RANSEquations) return

    ! If the turbulence model is wall distance free just compute the
    ! normal spacing of the first cell and store this in the wall
    ! distance. It may be needed for the boundary conditions and
    ! the monitoring of the y+. Return afterwards.

    if(.not. wallDistanceNeeded) then 

       ! Loop over the number of spectral solutions, initialize the
       ! distance and compute the initial normal spacing.

       do sps=1,nTimeIntervalsSpectral
          call initWallDistance(level, sps, allocMem)    
          call computeNormalSpacing(level, sps)
       enddo

       ! And return.

       return
    endif

    ! Write a message to stdout to indicate that the wall distance
    ! computation starts for the given level.

    write(integerString,"(i2)") level
    integerString = adjustl(integerString)

    ! Store the start time.

    t0 = mpi_wtime()

    ! Release temporarily some memory such that the overall memory
    ! requirement is not dictated by this routine. What memory is
    ! released depends on allocMem. If this is .True. It means that
    ! this is the first time the wall distances are computed, i.e. in
    ! the preprocessing phase. Then only send and receive buffers are
    ! released. If allocMem is .False., this means that the wall
    ! distances are computed in a moving mesh computation and
    ! consequently some more memory can be released.

    if( allocMem ) then
       deallocate(sendBuffer, recvBuffer, stat=ierr)
       if(ierr /= 0)                    &
            call terminate("wallDistance", &
            "Deallocation error for communication buffers")
    else
       call deallocateTempMemory(.false.)
    endif

    ! There are two different searches we can do: the original code
    ! always works and it capable to dealing with rotating/periodic
    ! geometries. It uses constant memory and is slow. The
    ! alternative method uses memory that scales with the size of
    ! the surface grid per processor and only works for
    ! steady/unsteady simulations without periodic/rotating
    ! components. But it is fast. It is designed to be used for
    ! updating the wall distances between iterations of
    ! aerostructural solutions. 

    ! Normal, original wall distance calc. Cannot be used when
    ! overset is present due to possibility of overlapping walls. 
    if (.not. useApproxWallDistance) then 
       ! Loop over the number of spectral solutions.
       spectralLoop: do sps=1,nTimeIntervalsSpectral

          ! Initialize the wall distances.

          call initWallDistance(level, sps, allocMem)

          ! Build the viscous surface mesh.

          call viscousSurfaceMesh(level, sps)

          ! If there are no viscous faces, processor 0 prints a warning
          ! and the wall distances are not computed.

          if(nquadViscGlob == 0) then

             if(myID == 0) then
                print "(a)", "#"
                print "(a)", "#               Warning!!!!"
                print "(a)", "# No viscous boundary found. Wall &
                     &distances are set to infinity"
                print "(a)", "#"
             endif

          else
             ! Determine the wall distances for the owned cells.
             call determineDistance(level, sps)
          end if
       end do spectralLoop
    else ! The user wants to use approx wall distance calcs OR we
       ! have overset mesh. :

       if (updateWallAssociation(level)) then 

          ! Initialize the wall distance
          spectralLoop2: do sps=1,nTimeIntervalsSpectral
             call initWallDistance(level, sps, allocMem)
          end do spectralLoop2

          ! Destroy the PETSc wall distance data if necessary
          call destroyWallDistanceDataLevel(level)

          ! Do the associtaion. This allocates the data destroyed
          ! in the destroyWallDistanceData call

          do sps=1, nTimeIntervalsSpectral
             call determineWallAssociation(level, sps)
          end do

          updateWallAssociation(level) = .False.
       end if

       ! Update the xsurf vector from X
       call updateXSurf(level)

       ! Call the actual update routine, on each of the sps instances and blocks
       do sps=1, nTimeIntervalsSpectral

          ! Now extract the vector of the surface data we need
          call VecGetArrayF90(xSurfVec(level, sps), xSurf, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          do nn=1,nDom
             call setPointers(nn, level, sps)
             call updateWallDistancesQuickly(nn, level, sps)
          end do

          call VecRestoreArrayF90(xSurfVec(level, sps), xSurf, ierr)
          call EChk(ierr,__FILE__,__LINE__)

       end do
    end if

    ! Allocate the temporarily released memory again. For more info
    ! see the comments at the beginning of this routine.

    if( allocMem ) then
       allocate(sendBuffer(sendBufferSize), &
            recvBuffer(recvBufferSize), stat=ierr)
       if(ierr /= 0)                    &
            call terminate("wallDistance", &
            "Memory allocation failure for comm buffers")
    else
       call allocateTempMemory(.false.)
    endif

    ! Synchronize the processors.

    call mpi_barrier(SUmb_comm_world, ierr)

    ! Write a message to stdout with the amount of time it
    ! took to compute the distances.

    !        if(myID == 0) then
    !          print 102, trim(integerString)
    !          print 103, mpi_wtime() - t0
    !          print "(a)", "#"
    !        endif
102 format("# End wall distances level",1X,A)
103 format("# Wall clock time:",E12.5," sec.")

  end subroutine computeWallDistance


  subroutine computeNormalSpacing(level, sps)
    !
    !       computeNormalSpacing computes the normal spacing of the first  
    !       cell center from the viscous wall for the given multigrid      
    !       level and spectral solution. This routine is called for        
    !       turbulence models, which do not need the wall distance.        
    !       However, they do need info of the first normal spacing for the 
    !       monitoring of y+ and possibly for the boundary conditions.     
    !       This is computed in this routine.                              
    !
    use constants
    use blockPointers, only : x, d2wall, nViscBocos, BCFaceID, BCData, &
         nx, ny, nz, il, jl, kl, nDom
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm, i, j
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

    real(kind=realType) :: nnx, nny, nnz, vecx, vecy, vecz, dot

    real(kind=realType), dimension(:,:,:), pointer :: xFace, xInt
    real(kind=realType), dimension(:,:),   pointer :: dd2Wall

    ! Loop over the domains.

    domain: do nn=1, nDom

       ! Set the pointers for this block.

       call setPointers(nn, level, sps)

       ! Loop over the viscous subfaces of this block. Note that
       ! these are numbered first.

       bocos: do mm=1,nViscBocos

          ! Set the pointers for the plane on the surface, one plane
          ! into the computational domain and the wall distance.
          ! This depends on the block face on which the subface is
          ! located. Note that the starting index of d2Wall is 2 and
          ! therefore a pointer offset will be needed later on.

          select case (BCFaceID(mm))

          case (iMin)
             xFace   => x(1, 1:,1:,:); xInt => x(2, 1:,1:,:)
             dd2Wall => d2Wall(2, :,:)

          case (iMax)
             xFace   => x(il,1:,1:,:); xInt => x(nx,1:,1:,:)
             dd2Wall => d2Wall(il,:,:)

          case (jMin)
             xFace   => x(1:,1, 1:,:); xInt => x(1:,2, 1:,:)
             dd2Wall => d2Wall(:,2 ,:)

          case (jMax)
             xFace   => x(1:,jl,1:,:); xInt => x(1:,ny,1:,:)
             dd2Wall => d2Wall(:,jl,:)

          case (kMin)
             xFace   => x(1:,1:,1, :); xInt => x(1:,1:,2 ,:)
             dd2Wall => d2Wall(:,:,2 )

          case (kMax)
             xFace   => x(1:,1:,kl,:); xInt => x(1:,1:,nz,:)
             dd2Wall => d2Wall(:,:,kl)

          end select

          ! Store the face range of this subface a bit easier.

          jBeg = BCData(mm)%jnBeg+1; jEnd = BCData(mm)%jnEnd
          iBeg = BCData(mm)%inBeg+1; iEnd = BCData(mm)%inEnd

          ! Loop over the faces of the subfaces.

          do j=jBeg,jEnd
             do i=iBeg,iEnd

                ! Store the three components of the unit normal a
                ! bit easier.

                nnx = BCData(mm)%norm(i,j,1)
                nny = BCData(mm)%norm(i,j,2)
                nnz = BCData(mm)%norm(i,j,3)

                ! Compute the vector from centroid of the adjacent cell
                ! to the centroid of the face.

                vecx = eighth*(xFace(i-1,j-1,1) + xFace(i-1,j,1) &
                     +         xFace(i,  j-1,1) + xFace(i,  j,1) &
                     -          xInt(i-1,j-1,1) -  xInt(i-1,j,1) &
                     -          xInt(i,  j-1,1) -  xInt(i,  j,1))

                vecy = eighth*(xFace(i-1,j-1,2) + xFace(i-1,j,2) &
                     +         xFace(i,  j-1,2) + xFace(i,  j,2) &
                     -          xInt(i-1,j-1,2) -  xInt(i-1,j,2) &
                     -          xInt(i,  j-1,2) -  xInt(i,  j,2))

                vecz = eighth*(xFace(i-1,j-1,3) + xFace(i-1,j,3) &
                     +         xFace(i,  j-1,3) + xFace(i,  j,3) &
                     -          xInt(i-1,j-1,3) -  xInt(i-1,j,3) &
                     -          xInt(i,  j-1,3) -  xInt(i,  j,3))

                ! Compute the projection of this vector onto the normal
                ! vector of the face. For a decent mesh there will not be
                ! much of a difference between the projection and the
                ! original mesh, but it does not hurt to do it.

                dot = nnx*vecx + nny*vecy + nnz*vecz

                ! As (nnx,nny,nnz) is a unit vector the distance to the
                ! wall of the first cell center is given by the absolute
                ! value of dot. Due to the use of pointers and the fact
                ! that the original d2Wall array starts at 2 and offset
                ! of -1 must be used to store the data at the correct
                ! location.

                dd2Wall(i-1,j-1) = abs(dot)

             enddo
          enddo

       enddo bocos

    enddo domain

  end subroutine computeNormalSpacing

  subroutine initWallDistance(level, sps, allocMem)
    !
    !       initWallDistance allocates the memory for the wall distance,   
    !       if needed, and initializes the wall distance to a large value. 
    !
    use constants
    use blockPointers, only : nDom, flowDoms
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps
    logical, intent(in)               :: allocMem
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, il, jl, kl

    ! Loop over the domains.

    domain: do nn=1, nDom

       ! Allocate the memory for d2Wall, if desired.

       if( allocMem ) then

          il = flowDoms(nn,level,sps)%il
          jl = flowDoms(nn,level,sps)%jl
          kl = flowDoms(nn,level,sps)%kl

          allocate(flowDoms(nn,level,sps)%d2Wall(2:il,2:jl,2:kl), &
               stat=ierr)
          if(ierr /= 0)                          &
               call terminate("initWallDistance", &
               "Memory allocation failure for d2Wall")
       endif

       ! Initialize the wall distances to a large value.

       flowDoms(nn,level,sps)%d2Wall = large

    enddo domain

  end subroutine initWallDistance

  subroutine determineDistance(level, sps)
    !
    !       determineDistance determines the distance from the center      
    !       of the cell to the nearest viscous wall for owned cells.       
    !
    use constants
    use adtAPI, only :adtBuildSurfaceADT, adtMinDistanceSearch, adtDeallocateADTs
    use blockPointers, only : x, flowDoms, kl, jl, il, nDom, nx, ny, nz, &
         sectionID, d2Wall
    use communication, only : sumb_comm_world
    use section, only : sections
    use inputPhysics, only : wallOffset
    use utils, only : setPointers, terminate
    implicit none
    !
    !      Subroutine arguments
    !
    integer(kind=intType), intent(in) :: level, sps
    !
    !      Local parameter, which defines the name of the adt to be create
    !      and used here. Only needed because of the api of the adt library
    !
    character(len=10), parameter :: viscAdt = "ViscousADT"
    !
    !      Local variables.
    !
    integer :: ierr

    integer, dimension(:), allocatable :: procID

    integer(kind=intType) :: nCell, nCellPer, nTria

    integer(kind=intType), dimension(1,1) :: connTria
    real(kind=realType),   dimension(3,2) :: dummy

    integer(kind=intType), dimension(:), allocatable :: elementID

    real(kind=realType), dimension(:),   allocatable :: dist2
    real(kind=realType), dimension(:),   allocatable :: dist2per
    real(kind=realType), dimension(:,:), allocatable :: coor, uvw
    real(kind=realType), dimension(:,:), allocatable :: coorPer

    integer(kind=adtElementType), dimension(:), allocatable :: &
         elementType

    integer(kind=intType) :: nn, mm, ll, ii, jj, i, j, k

    real(kind=realType), dimension(3) :: xc

    ! Build the adt of the surface grid. As the api requires to
    ! specify both the quadrilateral and the triangular connectivity,
    ! some dummy variables must be passed.

    nTria    = 0
    connTria = 0
    dummy    = zero

    call adtBuildSurfaceADT(nTria,    nquadVisc, nNodeVisc,       &
         coorVisc, connTria,  connVisc,        &
         dummy,    .false.,   SUmb_comm_world, &
         viscAdt)

    ! Determine the number of cell centers for which the distance
    ! must be computed. Also determine how many of them are periodic.

    nCell    = 0
    nCellPer = 0

    do nn=1,nDom
       ll    = flowDoms(nn,level,sps)%nx * flowDoms(nn,level,sps)%ny &
            * flowDoms(nn,level,sps)%nz
       nCell = nCell + ll

       mm = flowDoms(nn,level,sps)%sectionID
       if( sections(mm)%periodic ) nCellPer = nCellPer + ll
    enddo

    ! Allocate the memory for the arrays needed by the ADT.

    allocate(coor(3,nCell), procID(nCell), elementType(nCell), &
         elementID(nCell), uvw(3,nCell), dist2(nCell),      &
         coorPer(3,nCellPer), dist2Per(nCellPer), stat=ierr)
    if(ierr /= 0)                         &
         call terminate("determineDistance", &
         "Memory allocation failure for the variables &
         &needed by the adt.")
    !
    !       Step 1: The search of the original coordinates; possibly a     
    !               rotational periodic transformation is applied to align 
    !               the sections.                                          
    !
    ! Loop over the domains to store the coordinates of the cell
    ! centers of the owned cells. Apply the transformation such that
    ! the sections are aligned.

    mm = 0
    domains: do nn=1,nDom

       ! Set the pointers for this block, store the section id a bit
       ! easier in ll and loop over the cell centers.

       call setPointers(nn, level, sps)
       ll = sectionID

       do k=2,kl
          do j=2,jl
             do i=2,il

                ! Compute the coordinates of the cell center relative
                ! to the rotation center of this section.

                xc(1) = eighth*(x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
                     +         x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
                     +         x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                     +         x(i-1,j,  k,  1) + x(i,j,  k,  1)) &
                     - sections(ll)%rotCenter(1)

                xc(2) = eighth*(x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
                     +         x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
                     +         x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                     +         x(i-1,j,  k,  2) + x(i,j,  k,  2)) &
                     - sections(ll)%rotCenter(2)

                xc(3) = eighth*(x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
                     +         x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
                     +         x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                     +         x(i-1,j,  k,  3) + x(i,j,  k,  3)) &
                     - sections(ll)%rotCenter(3)

                ! Apply the periodic transformation for this section to
                ! align it with other sections and store this coordinate
                ! in the appropriate place in coor.

                mm = mm + 1
                coor(1,mm) = rotMatrixSections(ll,1,1)*xc(1) &
                     + rotMatrixSections(ll,1,2)*xc(2) &
                     + rotMatrixSections(ll,1,3)*xc(3) &
                     + sections(ll)%rotCenter(1)

                coor(2,mm) = rotMatrixSections(ll,2,1)*xc(1) &
                     + rotMatrixSections(ll,2,2)*xc(2) &
                     + rotMatrixSections(ll,2,3)*xc(3) &
                     + sections(ll)%rotCenter(2)

                coor(3,mm) = rotMatrixSections(ll,3,1)*xc(1) &
                     + rotMatrixSections(ll,3,2)*xc(2) &
                     + rotMatrixSections(ll,3,3)*xc(3) &
                     + sections(ll)%rotCenter(3)

                ! Initialize the distance squared, because this is an
                ! inout argument in the call to adtMinDistanceSearch.

                dist2(mm) = d2Wall(i,j,k)

             enddo
          enddo
       enddo
    enddo domains

    ! Perform the search. As no no interpolations are required,
    ! some dummies are passed.

    call adtMinDistanceSearch(nCell,  coor,        viscAdt,      &
         procID, elementType, elementID,    &
         uvw,    dist2,       0_intType, &
         dummy,  dummy)
    !        if (myid == 0) then
    !           print *,'procID = '
    !           print *,procID
    !           print *,'elemID:',elementID
    !        end if

    !
    !       Step 2: For periodic sections the nearest wall may be in the   
    !               periodic part of the grid that is not stored.          
    !               Therefore apply the periodic transformation to the     
    !               node and compute the minimum distance for this         
    !               coordinate.                                            
    !
    ! Initialize the counters mm and ii. Mm is the counter for coor;
    ! ii is the counter for coorPer.

    mm = 0
    ii = 0

    ! Loop over the domains and find the periodic ones.

    domainsPer1: do nn=1,nDom
       jj    = flowDoms(nn,level,sps)%nx * flowDoms(nn,level,sps)%ny &
            * flowDoms(nn,level,sps)%nz

       ll = flowDoms(nn,level,sps)%sectionID

       ! Check if the section is periodic.

       if( sections(ll)%periodic ) then

          ! Loop over the corresponding entries in coor of this block
          ! and apply the periodic transformation. The transformed
          ! coordinates are stored in coorPer. Also initialize
          ! the wall distance squared to the value just computed.

          do i=1,jj
             mm = mm + 1
             ii = ii + 1

             xc(1) = coor(1,mm) - sections(ll)%rotCenter(1)
             xc(2) = coor(2,mm) - sections(ll)%rotCenter(2)
             xc(3) = coor(3,mm) - sections(ll)%rotCenter(3)

             coorPer(1,ii) = sections(ll)%rotMatrix(1,1)*xc(1) &
                  + sections(ll)%rotMatrix(1,2)*xc(2) &
                  + sections(ll)%rotMatrix(1,3)*xc(3) &
                  + sections(ll)%rotCenter(1)         &
                  + sections(ll)%translation(1)

             coorPer(2,ii) = sections(ll)%rotMatrix(2,1)*xc(1) &
                  + sections(ll)%rotMatrix(2,2)*xc(2) &
                  + sections(ll)%rotMatrix(2,3)*xc(3) &
                  + sections(ll)%rotCenter(2)         &
                  + sections(ll)%translation(2)

             coorPer(3,ii) = sections(ll)%rotMatrix(3,1)*xc(1) &
                  + sections(ll)%rotMatrix(3,2)*xc(2) &
                  + sections(ll)%rotMatrix(3,3)*xc(3) &
                  + sections(ll)%rotCenter(3)         &
                  + sections(ll)%translation(3)

             dist2Per(ii) = dist2(mm)
          enddo

       else
          ! Section is not periodic. Update the counter mm.

          mm = mm + jj
       endif

    enddo domainsPer1

    ! Perform the adt search of this set of periodic coordinates.
    ! As no no interpolations are required, some dummies are passed.

    call adtMinDistanceSearch(nCellPer, coorPer,     viscAdt,      &
         procID,   elementType, elementID,    &
         uvw,      dist2Per,    0_intType, &
         dummy,    dummy)
    !
    !       Step 3: Also apply the inverse periodic transformation.        
    !
    ! Initialize the counters mm and ii. Mm is the counter for coor;
    ! ii is the counter for coorPer.

    mm = 0
    ii = 0

    ! Loop over the domains and find the periodic ones.

    domainsPer2: do nn=1,nDom
       jj    = flowDoms(nn,level,sps)%nx * flowDoms(nn,level,sps)%ny &
            * flowDoms(nn,level,sps)%nz

       ll = flowDoms(nn,level,sps)%sectionID

       ! Check if the section is periodic.

       if( sections(ll)%periodic ) then

          ! Loop over the corresponding entries in coor of this block
          ! and apply the inverse periodic transformation. Again the
          ! transformed coordinates are stored in coorPer. Note that
          ! the inverse of the rotation matrix is the transpose and
          ! that the translation vector should be multiplied by the
          ! inverse of the rotation matrix. Note that the wall distance
          ! has already been initialized in the previous periodic
          ! search.

          do i=1,jj
             mm = mm + 1
             ii = ii + 1

             xc(1) = coor(1,mm) - sections(ll)%rotCenter(1) &
                  -              sections(ll)%translation(1)
             xc(2) = coor(2,mm) - sections(ll)%rotCenter(2) &
                  -              sections(ll)%translation(2)
             xc(3) = coor(3,mm) - sections(ll)%rotCenter(3) &
                  -              sections(ll)%translation(3)

             coorPer(1,ii) = sections(ll)%rotMatrix(1,1)*xc(1) &
                  + sections(ll)%rotMatrix(2,1)*xc(2) &
                  + sections(ll)%rotMatrix(3,1)*xc(3) &
                  + sections(ll)%rotCenter(1)

             coorPer(2,ii) = sections(ll)%rotMatrix(1,2)*xc(1) &
                  + sections(ll)%rotMatrix(2,2)*xc(2) &
                  + sections(ll)%rotMatrix(3,2)*xc(3) &
                  + sections(ll)%rotCenter(2)

             coorPer(3,ii) = sections(ll)%rotMatrix(1,3)*xc(1) &
                  + sections(ll)%rotMatrix(2,3)*xc(2) &
                  + sections(ll)%rotMatrix(3,3)*xc(3) &
                  + sections(ll)%rotCenter(3)
          enddo

       else
          ! Section is not periodic. Update the counter mm.

          mm = mm + jj
       endif

    enddo domainsPer2

    ! Perform the adt search of this set of periodic coordinates.
    ! As no no interpolations are required, some dummies are passed.

    call adtMinDistanceSearch(nCellPer, coorPer,     viscAdt,      &
         procID,   elementType, elementID,    &
         uvw,      dist2Per,    0_intType, &
         dummy,    dummy)
    !
    !       Step 4: Store the minimum distance in the block type.          
    !
    mm = 0
    ii = 0

    domainsStore: do nn=1,nDom

       ! Set the pointers for this block and store the section id a
       ! bit easier in ll.

       call setPointers(nn, level, sps)
       ll = sectionID

       ! Check if the section is periodic. If so the distance is set
       ! to the minimum of the value computed in step 1 and the
       ! periodic values. Note that mm should not be updated in this
       ! loop, because it is updated in the loop over the cell centers
       ! of this block. Instead the counter j is used.

       if( sections(ll)%periodic ) then

          jj = nx*ny*nz

          j = mm
          do i=1,jj
             j  = j  + 1
             ii = ii + 1

             dist2(j) = min(dist2(j), dist2Per(ii))
          enddo
       endif

       ! Loop over the cell centers of the block to store the wall
       ! distance. Note that dist2 stores the distance squared and
       ! thus a square root must be taken.
       ! Add a possible offset for the debugging of wall functions.

       do k=2,kl
          do j=2,jl
             do i=2,il
                mm = mm + 1
                d2Wall(i,j,k) = sqrt(dist2(mm)) + wallOffset
             enddo
          enddo
       enddo

    enddo domainsStore

    ! Release the memory of the ADT and the arrays of the module
    ! viscSurface.

    call adtDeallocateADTs(viscAdt)

    deallocate(connVisc, coorVisc, rotMatrixSections, stat=ierr)
    if(ierr /= 0)                         &
         call terminate("determineDistance", &
         "Deallocation error for the arrays &
         &of viscSurface")

    ! Release the variables needed by the ADT.

    deallocate(coor, procID, elementType, elementID, uvw, dist2, &
         coorPer, dist2Per, stat=ierr)
    if(ierr /= 0)                          &
         call terminate("determineDistance", &
         "Deallocation failure for the variables &
         &needed by the adt.")

  end subroutine determineDistance


  subroutine localViscousSurfaceMesh(multSections, level, sps)
    !
    !       localViscousSurfaceMesh stores the local viscous surface       
    !       mesh (with possible periodic extensions in conn and coor.      
    !
    use constants
    use blockPointers, only : BCData, x, il, jl, kl, BCFaceID, sectionID, &
         flowDoms, nBocos, nDom, BCType
    use communication, only : myid, sumb_comm_world
    use section, only : sections, nSections
    use utils, only : setPointers, terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps
    integer(kind=intType), dimension(*), intent(in) :: multSections
    !
    !      Local variables.
    !
    integer :: size, ierr

    integer(kind=intType) :: nn, mm, i, j, k
    integer(kind=intType) :: np, nq, npOld, np1, nqOld, mp, mq
    integer(kind=intType) :: nq1, nq2, nq3, nq4, sec, row, col
    integer(kind=intType) :: iBeg, jBeg, iEnd, jEnd

    integer(kind=intType), dimension(3) :: ind

    real(kind=realType) :: length, dot, xx, yy, zz, r1, r2, aaa, bbb
    real(kind=realType) :: theta, cosTheta, sinTheta

    real(kind=realType), dimension(3,3) :: a
    real(kind=realType), dimension(nSections) :: thetaNMin, &
         thetaNMax, &
         thetaPMin, &
         thetaPMax, tmp
    real(kind=realType), dimension(nSections,3) :: rad1, rad2, axis

    real(kind=realType), dimension(:,:,:), pointer :: xface

    !       Determine the unit vectors of the local coordinate system      
    !       aligned with the rotation axis of the possible rotational      
    !       periodic section.                                              
    !
    do nn=1,nSections
       if(sections(nn)%nSlices == 1) cycle

       ! Section is rotational periodic. First determine the rotation
       ! axis. This is the eigenvector which corresponds to the
       ! eigenvalue 1 of the transformation matrix.

       ! Store rot - i in a and initialize ind.

       ind(1) = 1; ind(2) = 2; ind(3) = 3

       a(1,1) = sections(nn)%rotMatrix(1,1) - one
       a(1,2) = sections(nn)%rotMatrix(1,2)
       a(1,3) = sections(nn)%rotMatrix(1,3)

       a(2,1) = sections(nn)%rotMatrix(2,1)
       a(2,2) = sections(nn)%rotMatrix(2,2) - one
       a(2,3) = sections(nn)%rotMatrix(2,3)

       a(3,1) = sections(nn)%rotMatrix(3,1)
       a(3,2) = sections(nn)%rotMatrix(3,2)
       a(3,3) = sections(nn)%rotMatrix(3,3) - one

       ! Loop over the two times that Gaussian elimination must be
       ! applied.

       loopGauss: do k=1,2

          ! Find the largest value in the sub-matrix.

          aaa = abs(a(k,k)); row = k; col = k
          do j=k,3
             do i=k,3
                bbb = abs(a(i,j))
                if(bbb > aaa) then
                   aaa  = bbb
                   row = i
                   col = j
                endif
             enddo
          enddo

          ! Swap the rows k and row.

          do j=1,3
             aaa       = a(k,j)
             a(k,j)   = a(row,j)
             a(row,j) = aaa
          enddo

          ! Swap the colums k and col; also swap ind(k) and ind(col).

          i        = ind(k)
          ind(k)   = ind(col)
          ind(col) = i
          do i=1,3
             aaa       = a(i,k)
             a(i,k)   = a(i,col)
             a(i,col) = aaa
          enddo

          ! Perform gaussian eliMination, because now it's sure that
          ! the element (k,k) is non-zero.

          aaa = one/a(k,k)
          do i=(k+1),3
             bbb = a(i,k)*aaa
             do j=k,3
                a(i,j) = a(i,j) - bbb*a(k,j)
             enddo
          enddo

       enddo loopGauss

       ! Due to the full pivoting it is now guaranteed that the elements
       ! a(ind(1),ind(1)) and a(ind(2),ind(2)) are nonzero and
       ! a(ind(3),ind(3)) == zero. Remember that the rotation matrix
       ! only has 1 eigenvalue of one. Set axis(ind(3)) to one and
       ! determine the other two elements of the eigen vector.

       axis(nn,ind(3)) =   one
       axis(nn,ind(2)) = -(a(2,3)*axis(nn,ind(3)))/a(2,2)
       axis(nn,ind(1)) = -(a(1,3)*axis(nn,ind(3)) &
            +   a(1,2)*axis(nn,ind(2)))/a(1,1)

       ! Create a unit vector.

       length = one/sqrt(axis(nn,1)**2 + axis(nn,2)**2 + axis(nn,3)**2)
       axis(nn,1) = axis(nn,1)*length
       axis(nn,2) = axis(nn,2)*length
       axis(nn,3) = axis(nn,3)*length

       ! Make sure that the largest component of this vector is
       ! positive, such that a unique definition of the rotation
       ! axis is obtained. Use dot and length as a temporary
       ! storage.

       dot = axis(nn,1); length = abs(dot)
       if(abs(axis(nn,2)) > length) then
          dot = axis(nn,2); length = abs(dot)
       endif
       if(abs(axis(nn,3)) > length) then
          dot = axis(nn,3); length = abs(dot)
       endif

       if(dot < zero) then
          axis(nn,1) = -axis(nn,1)
          axis(nn,2) = -axis(nn,2)
          axis(nn,3) = -axis(nn,3)
       endif

       ! Determine the two vectors which determine the plane normal
       ! to the axis of rotation.

       ! Initial guess of rad1. First try the y-axis. If not good
       ! enough try the z-axis.

       if(abs(axis(nn,2)) < 0.707107_realType) then
          rad1(nn,1) = zero
          rad1(nn,2) = one
          rad1(nn,3) = zero
       else
          rad1(nn,1) = zero
          rad1(nn,2) = zero
          rad1(nn,3) = one
       endif

       ! Make sure that rad1 is normal to axis. Create a unit
       ! vector again.

       dot = rad1(nn,1)*axis(nn,1) + rad1(nn,2)*axis(nn,2) &
            + rad1(nn,3)*axis(nn,3)
       rad1(nn,1) = rad1(nn,1) - dot*axis(nn,1)
       rad1(nn,2) = rad1(nn,2) - dot*axis(nn,2)
       rad1(nn,3) = rad1(nn,3) - dot*axis(nn,3)

       length = one/(rad1(nn,1)**2 + rad1(nn,2)**2 + rad1(nn,3)**2)
       rad1(nn,1) = rad1(nn,1)*length
       rad1(nn,2) = rad1(nn,2)*length
       rad1(nn,3) = rad1(nn,3)*length

       ! Create the second vector which spans the radIal plane. This
       ! must be normal to both axis and rad1, i.e. the cross-product.

       rad2(nn,1) = axis(nn,2)*rad1(nn,3) - axis(nn,3)*rad1(nn,2)
       rad2(nn,2) = axis(nn,3)*rad1(nn,1) - axis(nn,1)*rad1(nn,3)
       rad2(nn,3) = axis(nn,1)*rad1(nn,2) - axis(nn,2)*rad1(nn,1)

    enddo

    ! Initialize the values of thetaNMin, etc.

    thetaNMin =  zero
    thetaNMax = -pi
    thetaPMin =  pi
    thetaPMax =  zero
    !
    !       Determine the local values of thetaNMin, etc. for the          
    !       different sections.                                            
    !
    do nn=1,nDom

       ! Store the section id of the block a bit easier. Continue with
       ! the next block if this section consist of only one slice.

       sec = flowDoms(nn,level,sps)%sectionID
       if(sections(sec)%nSlices == 1) cycle

       ! Set the pointers for this block.

       call setPointers(nn, level,sps)

       ! Initialize nq1, nq2, nq3 and nq4 to 0. These integers store
       ! the number of nodes in the first, second, third and fourth
       ! quadrant respectivily.

       nq1 = 0; nq2 = 0; nq3 = 0; nq4 = 0

       ! Loop over the nodes of this block.

       do k=1,kl
          do j=1,jl
             do i=1,il

                ! Determine the coordinates relative to the
                ! center of rotation.

                xx = x(i,j,k,1) - sections(sec)%rotCenter(1)
                yy = x(i,j,k,2) - sections(sec)%rotCenter(2)
                zz = x(i,j,k,3) - sections(sec)%rotCenter(3)

                ! Determine the radIal components in the local
                ! cylindrical coordinate system of the section.

                r1 = xx*rad1(sec,1) + yy*rad1(sec,2) + zz*rad1(sec,3)
                r2 = xx*rad2(sec,1) + yy*rad2(sec,2) + zz*rad2(sec,3)

                ! Determine the angle if r1 or r2 is nonzero.

                if((abs(r1) >= eps) .or. (abs(r2) >= eps)) then

                   theta = atan2(r2,r1)

                   ! Update the minimum and maximum angle for this
                   ! section, depending on the sign of theta.

                   if(theta >= zero) then
                      thetaPMin(sec) = min(thetaPMin(sec),theta)
                      thetaPMax(sec) = max(thetaPMax(sec),theta)
                   endif

                   if(theta <= zero) then
                      thetaNMin(sec) = min(thetaNMin(sec),theta)
                      thetaNMax(sec) = max(thetaNMax(sec),theta)
                   endif

                   ! Determine the quadrant in which this node is located
                   ! and update the corresponding counter.

                   if(theta <= -half*pi) then
                      nq3 = nq3 + 1
                   else if(theta <= zero) then
                      nq4 = nq4 + 1
                   else if(theta <= half*pi) then
                      nq1 = nq1 + 1
                   else
                      nq2 = nq2 + 1
                   endif

                endif

             enddo
          enddo
       enddo

       ! Modify the minimum and maximum angles if nodes are present
       ! in multiple quadrants.

       if(nq1 > 0 .and. nq4 > 0) then

          ! Nodes in both the 1st and 4th quadrant. Update the
          ! corresponding minimum and maximum angle.

          thetaNMax(sec) = zero
          thetaPMin(sec) = zero

       endif

       if(nq2 > 0 .and. nq3 > 0) then

          ! Nodes in both the 2nd and 3rd quadrant. Update the
          ! corresponding minimum and maximum angle.

          thetaNMin(sec) = -pi
          thetaPMax(sec) =  pi

       endif

    enddo

    ! Determine the minimum of the minimum angles and the maximum of
    ! the maximum angles for all sections.

    size = nSections
    call mpi_allreduce(thetaNMax, tmp, size, sumb_real, mpi_max, &
         SUmb_comm_world, ierr)
    thetaNMax = tmp

    call mpi_allreduce(thetaPMax, tmp, size, sumb_real, mpi_max, &
         SUmb_comm_world, ierr)
    thetaPMax = tmp

    call mpi_allreduce(thetaNMin, tmp, size, sumb_real, mpi_min, &
         SUmb_comm_world, ierr)
    thetaNMin = tmp

    call mpi_allreduce(thetaPMin, tmp, size, sumb_real, mpi_min, &
         SUmb_comm_world, ierr)
    thetaPMin = tmp

    ! Allocate the memory for rotMatrixSections, the rotation
    ! matrices of the sections needed for the alignment.

    allocate(rotMatrixSections(nSections,3,3), stat=ierr)
    if(ierr /= 0)                                  &
         call terminate("localViscousSurfaceMesh", &
         "Memory allocation failure for &
         &rotMatrixSections")
    !
    !       Determine the rotation matrix for each section, which aligns   
    !       the rotational periodic sections with other sections.          
    !
    do nn=1,nSections

       ! Test if a rotation is actually needed.

       testRot: if(sections(nn)%nSlices == 1 .or. &
            thetaPMin(nn) == zero) then

          ! Section consist out of 1 slice or the slice crosses the
          ! line theta == 0. For both cases the rotation matrix is
          ! the identity matrix.

          rotMatrixSections(nn,1,1) = one
          rotMatrixSections(nn,1,2) = zero
          rotMatrixSections(nn,1,3) = zero

          rotMatrixSections(nn,2,1) = zero
          rotMatrixSections(nn,2,2) = one
          rotMatrixSections(nn,2,3) = zero

          rotMatrixSections(nn,3,1) = zero
          rotMatrixSections(nn,3,2) = zero
          rotMatrixSections(nn,3,3) = one

       else testRot

          ! Section consist out of multiple slices and the current slice
          ! does not cross the line theta == 0. The rotation matrix
          ! for alignment must be computed.

          theta = two*pi/sections(nn)%nSlices

          ! Determine the number of rotations needed to align the mesh.

          if(thetaNMin(nn) < zero) then

             ! The section lies (at least partially) in the third and
             ! fourth quadrant. Determine the number of rotations for
             ! alignment; this is a positive number.

             mm = -thetaNMax(nn)/theta + 1

          else

             ! The section lies (completely) in the first and second
             ! quadrant. The number of rotations will be a negative
             ! number now.

             mm = -thetaPMin(nn)/theta - 1

          endif

          ! Compute the rotation angle in the local cylindrical frame
          ! and its sine and cosine.

          theta = mm*theta
          cosTheta = cos(theta)
          sinTheta = sin(theta)

          ! Apply the transformation to obtain the matrix in the
          ! original cartesian frame.

          rotMatrixSections(nn,1,1) = axis(nn,1)*axis(nn,1)             &
               + cosTheta*(rad1(nn,1)*rad1(nn,1) + rad2(nn,1)*rad2(nn,1))
          rotMatrixSections(nn,1,2) = axis(nn,1)*axis(nn,2)             &
               + cosTheta*(rad1(nn,1)*rad1(nn,2) + rad2(nn,1)*rad2(nn,2)) &
               + sinTheta*(rad1(nn,2)*rad2(nn,1) - rad1(nn,1)*rad2(nn,2))
          rotMatrixSections(nn,1,3) = axis(nn,1)*axis(nn,3)             &
               + cosTheta*(rad1(nn,1)*rad1(nn,3) + rad2(nn,1)*rad2(nn,3)) &
               + sinTheta*(rad1(nn,3)*rad2(nn,1) - rad1(nn,1)*rad2(nn,3))

          rotMatrixSections(nn,2,1) = axis(nn,1)*axis(nn,2)             &
               + cosTheta*(rad1(nn,1)*rad1(nn,2) + rad2(nn,1)*rad2(nn,2)) &
               - sinTheta*(rad1(nn,2)*rad2(nn,1) - rad1(nn,1)*rad2(nn,2))
          rotMatrixSections(nn,2,2) = axis(nn,2)*axis(nn,2)             &
               + cosTheta*(rad1(nn,2)*rad1(nn,2) + rad2(nn,2)*rad2(nn,2))
          rotMatrixSections(nn,2,3) = axis(nn,2)*axis(nn,3)             &
               + cosTheta*(rad1(nn,2)*rad1(nn,3) + rad2(nn,2)*rad2(nn,3)) &
               + sinTheta*(rad1(nn,3)*rad2(nn,2) - rad1(nn,2)*rad2(nn,3))

          rotMatrixSections(nn,3,1) = axis(nn,1)*axis(nn,3)             &
               + cosTheta*(rad1(nn,1)*rad1(nn,3) + rad2(nn,1)*rad2(nn,3)) &
               - sinTheta*(rad1(nn,3)*rad2(nn,1) - rad1(nn,1)*rad2(nn,3))
          rotMatrixSections(nn,3,2) = axis(nn,2)*axis(nn,3)             &
               + cosTheta*(rad1(nn,2)*rad1(nn,3) + rad2(nn,2)*rad2(nn,3)) &
               - sinTheta*(rad1(nn,3)*rad2(nn,2) - rad1(nn,2)*rad2(nn,3))
          rotMatrixSections(nn,3,3) = axis(nn,3)*axis(nn,3)             &
               + cosTheta*(rad1(nn,3)*rad1(nn,3) + rad2(nn,3)*rad2(nn,3))

       endif testRot

    enddo
    !
    !       Determine the local viscous surface grid.                      
    !
    np = 0
    nq = 0

    loopDomains: do nn=1,nDom

       ! Set the pointers for this block and store the section id
       ! a bit easier.

       call setPointers(nn, level, sps)
       sec = sectionID

       ! Loop over the subfaces of this block and test if this is
       ! a viscous subface.

       loopBocos: do mm=1,nBocos
          testViscous: if(BCType(mm) == NSWallAdiabatic .or. &
               BCType(mm) == NSWallIsothermal) then

             ! Viscous subface. Set the pointer for the coordinates of
             ! the face.

             select case (BCFaceID(mm))

             case (iMin)
                xface => x(1,1:,1:,:)

             case (iMax)
                xface => x(il,1:,1:,:)

             case (jMin)
                xface => x(1:,1,1:,:)

             case (jMax)
                xface => x(1:,jl,1:,:)

             case (kMin)
                xface => x(1:,1:,1,:)

             case (kMax)
                xface => x(1:,1:,kl,:)

             end select

             ! Store the nodal range of this subface a bit easier.

             jBeg = BCData(mm)%jnBeg; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd

             ! Store the old value of the number of points stored and
             ! determine the new coordinates.

             npOld = np
             do j=jBeg,jEnd
                do i=iBeg,iEnd

                   ! Determine the coordinates relative to the rotation
                   ! center of this section.

                   xx = xface(i,j,1) - sections(sec)%rotCenter(1)
                   yy = xface(i,j,2) - sections(sec)%rotCenter(2)
                   zz = xface(i,j,3) - sections(sec)%rotCenter(3)

                   ! Update the counter and determine the surface mesh
                   ! coordinates.

                   np = np + 1
                   coorVisc(1,np) = rotMatrixSections(sec,1,1)*xx &
                        + rotMatrixSections(sec,1,2)*yy &
                        + rotMatrixSections(sec,1,3)*zz &
                        + sections(sec)%rotCenter(1)

                   coorVisc(2,np) = rotMatrixSections(sec,2,1)*xx &
                        + rotMatrixSections(sec,2,2)*yy &
                        + rotMatrixSections(sec,2,3)*zz &
                        + sections(sec)%rotCenter(2)

                   coorVisc(3,np) = rotMatrixSections(sec,3,1)*xx &
                        + rotMatrixSections(sec,3,2)*yy &
                        + rotMatrixSections(sec,3,3)*zz &
                        + sections(sec)%rotCenter(3)
                enddo
             enddo

             ! Determine and store the connectivity of this subface.

             np1 = iEnd - iBeg + 1
             nqOld = nq

             do j=(jBeg+1),jEnd
                do i=(iBeg+1),iEnd

                   ! Update the counter nq and determine the 4 indices
                   ! of the surface quad.

                   nq = nq + 1

                   connVisc(1,nq) = npOld + (j-jBeg-1)*np1 + i - iBeg
                   connVisc(2,nq) = connVisc(1,nq) + 1
                   connVisc(3,nq) = connVisc(2,nq) + np1
                   connVisc(4,nq) = connVisc(3,nq) - 1

                enddo
             enddo

             ! Loop over the number of times the subface must be stored.
             ! This happens when the rotational periodicity differs from
             ! section to section. Note that this loop starts at k == 2.

             loopMultiplicity: do k=2,multSections(sec)

                ! Store the current number of nodes and quads
                ! in mp and mq respectivily.

                mp = np
                mq = nq

                ! Loop over the of points on this subface.

                do i=(npOld+1),mp

                   ! Determine the coordinates relative to the center
                   ! of rotation.

                   np = np + 1

                   xx = coorVisc(1,i) - sections(sec)%rotCenter(1)
                   yy = coorVisc(2,i) - sections(sec)%rotCenter(2)
                   zz = coorVisc(3,i) - sections(sec)%rotCenter(3)

                   ! Update the counter np and determine the new
                   ! coordinates after the transformation.

                   coorVisc(1,np) = sections(sec)%rotMatrix(1,1)*xx &
                        + sections(sec)%rotMatrix(1,2)*yy &
                        + sections(sec)%rotMatrix(1,3)*zz &
                        + sections(sec)%rotCenter(1)      &
                        + sections(sec)%translation(1)

                   coorVisc(2,np) = sections(sec)%rotMatrix(2,1)*xx &
                        + sections(sec)%rotMatrix(2,2)*yy &
                        + sections(sec)%rotMatrix(2,3)*zz &
                        + sections(sec)%rotCenter(2)      &
                        + sections(sec)%translation(2)

                   coorVisc(3,np) = sections(sec)%rotMatrix(3,1)*xx &
                        + sections(sec)%rotMatrix(3,2)*yy &
                        + sections(sec)%rotMatrix(3,3)*zz &
                        + sections(sec)%rotCenter(3)      &
                        + sections(sec)%translation(3)
                enddo

                ! Store the number of nodes in this subface in j
                ! and determine the connectivity of this rotated part.

                j = np - mp
                do i=(nqOld+1),mq

                   ! Update the counter nq and set the new connectivity,
                   ! which is the old connectivity plus an offset.

                   nq = nq + 1
                   connVisc(1,nq) = connVisc(1,i) + j
                   connVisc(2,nq) = connVisc(2,i) + j
                   connVisc(3,nq) = connVisc(3,i) + j
                   connVisc(4,nq) = connVisc(4,i) + j

                enddo

                ! Copy the values of mp and mq into npOld and nqOld for
                ! the next multiple of the slice. Idem for np in mp, etc.

                npOld = mp
                nqOld = mq

                mp = np
                mq = nq

             enddo loopMultiplicity

          endif testViscous
       enddo loopBocos
    enddo loopDomains

  end subroutine localViscousSurfaceMesh

  subroutine updateWallDistanceAllLevels
    !
    !       updateWallDistanceAllLevels updates the wall distances for     
    !       the cell centers on all grid levels. This routine is typically 
    !       called when grid parts have been moved, either due to a        
    !       physical motion of some parts or due to deformation.           
    !
    use constants
    use block, only : flowDoms
    use inputPhysics, only : equations
    use iteration, only : groundLevel
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: nLevels, nn

    ! Return immediately if the rans equations are not solved.

    if(equations /= RANSEquations) return

    ! Loop over the grid levels and call wallDistance.

    nLevels = ubound(flowDoms,2)
    do nn=groundLevel,nLevels
       call computeWallDistance(nn, .false.)
    enddo

  end subroutine updateWallDistanceAllLevels

  subroutine viscousSurfaceMesh(level, sps)
    !
    !       viscousSurfaceMesh determines and stores the entire viscous    
    !       surface possibly extended by periodic parts.                   
    !
    use constants
    use block, only : flowDoms, nDom
    use communication, only : sumb_comm_world, myid
    use section, only : nsections, sections
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, mm, ii
    integer(kind=intType) :: ni, nj, nk

    integer(kind=intType), dimension(nSections) :: multSections

    ! Determine the minimum number of slices present in a certain
    ! section. For time accurate computations the number of slices
    ! is identical for all sections, but for steady flow using the
    ! mixing plane assumption this is not necessarily the case.

    mm = sections(1)%nSlices
    do nn=2,nSections
       mm = min(mm,sections(nn)%nSlices)
    enddo

    ! Determine the multiplicity of every section needed in the
    ! surface mesh, such that every part covers an angle which is
    ! at least equal to the angle of the largest section. Again note
    ! that this multiplicity is 1 for all sections if a time
    ! accurate computation is performed.

    do nn=1,nSections
       multSections(nn) = sections(nn)%nSlices/mm
       if(sections(nn)%nSlices > mm*multSections(nn)) &
            multSections(nn) = multSections(nn) + 1
    enddo

    ! Determine the local number of viscous nodes and quads.
    ! Note that these numbers are identical for all spectral
    ! solutions and thus it is okay to take the 1st one.

    nNodeVisc = 0
    nquadVisc = 0

    do nn=1,nDom
       do mm=1,flowDoms(nn,level,1)%nBocos
          if(flowDoms(nn,level,1)%BCType(mm) == NSWallAdiabatic .or. &
               flowDoms(nn,level,1)%BCType(mm) == NSWallIsothermal) then

             ! Determine the number of nodes of the subface in the
             ! three directions.

             ni = flowDoms(nn,level,1)%inEnd(mm) &
                  - flowDoms(nn,level,1)%inBeg(mm)
             nj = flowDoms(nn,level,1)%jnEnd(mm) &
                  - flowDoms(nn,level,1)%jnBeg(mm)
             nk = flowDoms(nn,level,1)%knEnd(mm) &
                  - flowDoms(nn,level,1)%knBeg(mm)

             ! Determine the multiplication factor, because of the
             ! possible multiple sections.

             ii = flowDoms(nn,level,1)%sectionId
             ii = multSections(ii)

             ! Update the number of nodes and quads. Take the
             ! multiplicity into account.

             nNodeVisc = nNodeVisc + ii*(ni+1)*(nj+1)*(nk+1)
             nquadVisc = nquadVisc + ii*max(ni,1_intType) &
                  *                max(nj,1_intType) &
                  *                max(nk,1_intType)
          endif
       enddo
    enddo

    ! Determine the global number of elements on the viscous
    ! surfaces. Return if there are no viscous quads present.

    call mpi_allreduce(nQuadVisc, nquadViscGlob, 1, sumb_integer, &
         mpi_sum, SUmb_comm_world, ierr)

    if(nquadViscGlob == 0) return

    ! Allocate the memory for the local connectivity and coordinates.

    allocate(connVisc(4,nquadVisc), coorVisc(3,nNodeVisc), &
         stat=ierr)
    if(ierr /= 0)                            &
         call terminate("viscousSurfaceMesh", &
         "Memory allocation failure for connVisc &
         &and coorVisc.")

    ! Determine the local viscous surface mesh, possibly rotated
    ! to align the other sections.

    call localViscousSurfaceMesh(multSections, level, sps)

  end subroutine viscousSurfaceMesh

  subroutine determineWallAssociation(level, sps)

    ! This routine will determine the closest surface point for every
    ! field cell. Special treatment is required for overlapping surfaces. 

    use constants
    use adtAPI, only : minDistanceSearch
    use adtData, only : adtBBOXTargetType
    use adtLocalSearch, only :  mindistancetreesearchsinglepoint
    use adtUtils, only : stack
    use adtBuild, only : buildSerialQuad, destroyserialquad
    use blockPointers
    use communication
    use inputphysics
    use inputTimeSpectral
    use overset, only : oversetPresent, oversetWall, nClusters, clusters, cumDomProc
    use inputOverset
    use adjointVars
    use utils, only : setPointers, EChk
    use sorting, only : unique
    implicit none

    ! Input Variables
    integer(kind=intType), intent(in) :: level, sps

    ! Local Variables
    integer(kind=intType) :: i, j, k, l, ii, jj, kk, nn, mm, iNode, iCell, c
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ni, nj, nUnique, cellID, cellID2
    integer(kind=intType) :: ierr, iDim

    ! Data for local surface
    integer(kind=intType) :: nNodes, nCells
    logical :: gridHasOverset

    ! Overset Walls for storing the surface ADT's
    type(oversetWall), dimension(:), allocatable, target :: walls
    type(oversetWall), target :: fullWall
    integer(kind=intType), dimension(:),  allocatable :: link, indicesToGet

    ! Data for the ADT
    integer(kind=intType) :: intInfo(3), intInfo2(3)
    real(kind=realType) :: coor(4), uvw(5), uvw2(5)
    real(kind=realType), dimension(3, 2) :: dummy
    real(kind=realType), parameter :: tol=1e-12
    integer(kind=intType), dimension(:), pointer :: frontLeaves, frontLeavesNew, BBint
    type(adtBBoxTargetType), dimension(:), pointer :: BB
    real(kind=realType), dimension(3) :: xp

    ! The first thing we do is gather all the surface nodes to
    ! each processor such that every processor can make it's own copy of
    ! the complete surface mesh to use to search. Note that this
    ! procedure *DOES NOT SCALE IN MEMORY*...ie eventually the surface
    ! mesh will become too large to store on a single processor,
    ! although this will probably not happen until the sizes get up in
    ! the hundreds of millions of cells. 

    allocate(walls(nClusters))
    call buildClusterWalls(level, sps, .False., walls)

    if (oversetPresent) then 
       ! Finally build up a "full wall" that is made up of all the cluster
       ! walls. 

       nNodes = 0
       nCells = 0
       do i=1, nClusters
          nNodes = nNodes+ walls(i)%nNodes
          nCells = nCells + walls(i)%nCells
       end do

       allocate(fullWall%x(3, nNodes))
       allocate(fullWall%conn(4, nCells))
       allocate(fullWall%ind(nNodes))

       nNodes = 0
       nCells = 0
       ii = 0
       do i=1, nClusters

          ! Add in the nodes/elements from this cluster

          do j=1, walls(i)%nNodes
             nNodes = nNodes + 1
             fullWall%x(:, nNodes) = walls(i)%x(:, j)
             fullWall%ind(nNodes) = walls(i)%ind(j)
          end do

          do j=1, walls(i)%nCells
             nCells = nCells + 1
             fullWall%conn(:, nCells) = walls(i)%conn(:, j) + ii
          end do

          ! Increment the node offset
          ii = ii + walls(i)%nNodes
       end do

       ! Finish the setup of the full wall.
       fullWall%nCells = nCells
       fullWall%nNodes = nNodes
       call buildSerialQuad(nCells, nNodes, fullWall%x, fullWall%conn, fullWall%ADT)
    end if

    ! Allocate the (pointer) memory that may be resized as necessary for
    ! the singlePoint search routine. 
    allocate(stack(100), BB(20), BBint(20), frontLeaves(25), frontLeavesNew(25))

    ! We need to store the 4 global node indices defining the quad that
    ! each point has the closest point wrt. We also ned to store the uv
    ! values. This allows us to recompute the exact surface point, after
    ! the rquired nodes are fetched from (a possibly) remote proc. 

    do nn=1,nDom
       call setPointers(nn, level, sps)

       ! Check if elemID and uv are allocated yet.
       if (.not. associated(flowDoms(nn,level,sps)%surfNodeIndices)) then
          allocate(flowDoms(nn,level,sps)%surfNodeIndices(4, 2:il, 2:jl, 2:kl))
          allocate(flowDoms(nn,level,sps)%uv(2, 2:il, 2:jl, 2:kl))
       end if

       ! Set the cluster for this block
       c = clusters(cumDomProc(myid) + nn)

       do k=2, kl
          do j=2, jl
             do i=2, il

                ! Compute the coordinates of the cell center 
                coor(1) = eighth*(x(i-1,j-1,k-1,1) + x(i,j-1,k-1,1)  &
                     +         x(i-1,j,  k-1,1) + x(i,j,  k-1,1)  &
                     +         x(i-1,j-1,k,  1) + x(i,j-1,k,  1)  &
                     +         x(i-1,j,  k,  1) + x(i,j,  k,  1))

                coor(2) = eighth*(x(i-1,j-1,k-1,2) + x(i,j-1,k-1,2)  &
                     +         x(i-1,j,  k-1,2) + x(i,j,  k-1,2)  &
                     +         x(i-1,j-1,k,  2) + x(i,j-1,k,  2)  &
                     +         x(i-1,j,  k,  2) + x(i,j,  k,  2))

                coor(3) = eighth*(x(i-1,j-1,k-1,3) + x(i,j-1,k-1,3)  &
                     +         x(i-1,j,  k-1,3) + x(i,j,  k-1,3)  &
                     +         x(i-1,j-1,k,  3) + x(i,j-1,k,  3)  &
                     +         x(i-1,j,  k,  3) + x(i,j,  k,  3))

                if (.not. oversetPresent) then 
                   ! No overset present. Simply search our own wall,
                   ! walls(c), (the only one we have) up to the wall
                   ! cutoff.
                   coor(4) = wallDistCutoff**2
                   intInfo(3) = 0 ! Must be initialized since the search
                   ! may not find closer point.
                   call minDistancetreeSearchSinglePoint(walls(c)%ADT, coor, intInfo, &
                        uvw, dummy, 0, BB, frontLeaves, frontLeavesNew)

                   cellID = intInfo(3)
                   if (cellID > 0) then 
                      do kk=1,4
                         flowDoms(nn, level, sps)%surfNodeIndices(kk, i, j, k) = &
                              walls(c)%ind(walls(c)%conn(kk, cellID))
                      end do
                      flowDoms(nn, level, sps)%uv(:, i, j, k) = uvw(1:2)
                   else
                      ! Just set dummy values. These will never be used. 
                      flowDoms(nn, level, sps)%surfNodeIndices(:, i, j, k) = 0
                      flowDoms(nn, level, sps)%uv(:, i, j, k) = 0
                   end if

                   ! We are done with this point. 
                   cycle
                end if

                ! This is now the overset (possibly) overlapping surface
                ! mesh case. It is somewhat more complex since we use
                ! the same searches to flag cells that are inside the
                ! body.

                coor(4) = wallDistCutoff**2
                intInfo(3) = 0
                call minDistancetreeSearchSinglePoint(fullWall%ADT, coor, &
                     intInfo, uvw, dummy, 0, BB, frontLeaves, frontLeavesNew)
                cellID = intInfo(3)

                if (cellID > 0) then
                   ! We found the cell:

                   ! If the cell is outside of near-wall distance or our
                   ! cluster doesn't have any owned cells. Just accept it. 
                   if (uvw(4) > nearWallDist**2 .or. walls(c)%nCells == 0) then 

                      do kk=1,4
                         flowDoms(nn, level, sps)%surfNodeIndices(kk, i, j, k) = &
                              fullWall%ind(fullWall%conn(kk, cellID))
                      end do
                      flowDoms(nn, level, sps)%uv(:, i, j, k) = uvw(1:2)

                   else

                      ! This point is *closer* than the nearWallDist AND
                      ! it has a wall. Search on our own wall.

                      coor(4) = large
                      call minDistancetreeSearchSinglePoint(walls(c)%ADT, coor, &
                           intInfo2, uvw2, dummy, 0, BB, frontLeaves, frontLeavesNew)
                      cellID2 = intInfo2(3)

                      if (uvw2(4) < nearWallDist**2) then 
                         ! Both are close to the wall. Accept the one
                         ! from our own wall unconditionally.
                         do kk=1,4
                            flowDoms(nn, level, sps)%surfNodeIndices(kk, i, j, k) = &
                                 walls(c)%ind(walls(c)%conn(kk, cellID2))
                         end do
                         flowDoms(nn, level, sps)%uv(:, i, j, k) = uvw2(1:2)
                      else
                         ! The full wall distance is better. Take that. 

                         do kk=1,4
                            flowDoms(nn, level, sps)%surfNodeIndices(kk, i, j, k) = &
                                 fullWall%ind(fullWall%conn(kk, cellID))
                         end do
                         flowDoms(nn, level, sps)%uv(:, i, j, k) = uvw(1:2)

                      end if
                   end if
                else

                   ! What happend here is a cell is outside the
                   ! wallDistCutoff. We don't care about wall distance
                   ! info here so just set dummy info.

                   flowDoms(nn, level, sps)%surfNodeIndices(:, i, j, k) = 0
                   flowDoms(nn, level, sps)%uv(:, i, j, k) = 0

                end if
             end do
          end do
       end do
    end do

    ! Now determine all the node indices this processor needs to get. 
    mm = 0
    allocate(indicesToGet(nCellsLocal(level)*4), link(nCellsLocal(level)*4))
    do nn=1, nDom
       call setPointers(nn, level, sps)
       do k=2, kl
          do j=2, jl
             do i=2, il
                do kk=1,4
                   mm = mm + 1
                   indicesToGet(mm) = flowDoms(nn, level, sps)%surfNodeIndices(kk, i, j, k)
                end do
             end do
          end do
       end do
    end do

    ! This unique-ifies the indices. 
    call unique(indicesToGet, 4*nCellsLocal(level), nUnique, link)

    ! we need to update the stored indices to use the ordering of the nodes we will receive. 
    mm = 0
    do nn=1, nDom
       call setPointers(nn, level, sps)
       do k=2, kl
          do j=2, jl
             do i=2, il
                do kk=1,4
                   mm = mm + 1
                   flowDoms(nn, level, sps)%surfNodeIndices(kk, i, j, k) = link(mm)
                end do
             end do
          end do
       end do
    end do
    deallocate(link)

    ! Now create the index set for the nodes we need to get. We have to
    ! expand "indices to get" to include the DOF. Use link for this
    ! temporary array operation.

    allocate(link(nUnique*3))
    do i=1, nUnique
       link((i-1)*3+1) = indicesToGet(i)*3
       link((i-1)*3+2) = indicesToGet(i)*3+1
       link((i-1)*3+3) = indicesToGet(i)*3+2
    end do

    call ISCreateGeneral(sumb_comm_world, nUnique*3, link, PETSC_COPY_VALUES, IS1, ierr)
    call EChk(ierr,__FILE__,__LINE__)
    deallocate(link)

    ! Create the volume vector the nodes will be scatter from. Note that
    ! this vector contains all the spectal instances. It is therefore
    ! only allocated on the first call with sps=1
    if (sps == 1) then 
       call VecCreateMPI(SUMB_COMM_WORLD, 3*nNodesLocal(level)*nTimeIntervalsSpectral, &
            PETSC_DETERMINE, xVolumeVec(level), ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end if

    ! This is the vector we will scatter the nodes into. 
    call VecCreateMPI(SUMB_COMM_WORLD, 3*nUnique, PETSC_DETERMINE, &
         xSurfVec(level, sps), ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecGetOwnershipRange(xSurfVec(level, sps), i, j, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call ISCreateStride(SUMB_COMM_WORLD, j-i, i, 1, IS2, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Create the actual final scatter context.
    call  VecScatterCreate(xVolumeVec(level), IS1, xSurfVec(level, sps), IS2, &
         wallScatter(level, sps), ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call ISDestroy(IS1, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call ISDestroy(IS2, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Deallocate all the remaining temporary data
    deallocate(stack, BB, frontLeaves, frontLeavesNew, BBint)

    do i=1, nClusters
       deallocate(walls(i)%x, walls(i)%conn, walls(i)%ind)
       call destroySerialQuad(walls(i)%ADT)
    end do
    deallocate(walls)

    if (oversetPresent) then 
       deallocate(fullWall%x, fullWall%conn, fullWall%ind)
       call destroySerialQuad(fullWall%ADT)
    end if

  end subroutine determineWallAssociation

  subroutine updateXSurf(level)

    use blockPointers
    use inputTimeSpectral
    use utils, only : EChk, setPointers
    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in) :: level

    ! Working Parameters
    integer(kind=intType) :: ii, i,j,k,l, nn, sps, ierr

    ! Fill up xVolumeVec 
    call VecGetArrayF90(xVolumeVec(level), xVolume, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ii = 0
    do nn=1,nDom
       do sps=1,nTimeIntervalsSpectral
          call setPointers(nn, level, sps)
          do k=1, kl
             do j=1, jl
                do i=1, il
                   do l= 1,3
                      ii = ii + 1
                      xVolume(ii) = X(i, j, k, l)
                   end do
                end do
             end do
          end do
       end do
    end do
    call vecRestoreArrayF90(xVolumeVec(level), xVolume, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Perform the scatter from the global x vector to xSurf. SPS loop since the xSurfVec is done by SPS instance.
    do sps=1, nTimeIntervalsSpectral
       call VecScatterBegin(wallScatter(level, sps), xVolumeVec(level), &
            xSurfVec(level, sps), INSERT_VALUES, SCATTER_FORWARD, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecScatterEnd(wallScatter(level, sps), xVolumeVec(level), &
            xSurfVec(level, sps), INSERT_VALUES, SCATTER_FORWARD, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end do

  end subroutine updateXSurf

  subroutine destroyWallDistanceData

    use constants
    use block, only : flowDoms
    implicit none

    ! Working
    integer(kind=intType) :: level, nLevels, l

    nLevels = ubound(flowDoms,2)
    do l=1,nLevels
       call destroyWallDistanceDataLevel(l)
    end do

    deallocate(xSurfVec, xVolumeVec, wallScatter)

  end subroutine destroyWallDistanceData

  subroutine destroyWallDistanceDataLevel(level)
    use constants
    use inputTimeSpectral, only :nTimeIntervalsspectral
    use utils, onlY : EChk
    use block, only : flowDoms
   
    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in) :: level

    ! Working
    integer(kind=intType) :: ierr, sps    

    ! Determine if we need to deallocate the PETSc data for
    ! this level
    if (wallDistanceDataAllocated(level)) then 
       call VecDestroy(xVolumeVec(level), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       do sps=1, nTimeIntervalsSpectral
          call VecDestroy(xSurfVec(level, sps), ierr)
          call EChk(ierr,__FILE__,__LINE__)
          
          call VecScatterDestroy(wallScatter(level, sps), ierr)
          call EChk(ierr,__FILE__,__LINE__)
       end do
       
       wallDistanceDataAllocated(level) = .False.
    end if
  end subroutine destroyWallDistanceDataLevel

#endif
end module wallDistance

