module userSurfaceIntegrations

  use constants
  use communication, only : commType, internalCommType

  type userSurfCommType
     ! Data required on each proc:

     ! nDonor: The number of donor points the proc will provide
     ! frac (3, nDonor) : The uvw coordinates of the interpolation point
     ! donorInfo(4, nDonor) : Donor information. 1 is the local block ID and 2-4 is the 
     !    starting i,j,k indices for the interpolation. 
     ! procSizes(0:nProc-1) : The number of donors on each proc
     ! procDisps(0:nProc) : Cumulative form of procSizes

     ! inv(nConn) : Array allocated only on root processor used to
     ! reorder the nodes or elements back to the original order. 

     integer(kind=intType) :: nDonor
     real(kind=realType), dimension(:,:), allocatable :: frac
     integer(kind=intType), dimension(:, :), allocatable :: donorInfo
     integer(kind=intTYpe), dimension(:), allocatable :: procSizes, procDisps
     integer(kind=intTYpe), dimension(:), allocatable :: inv
     logical, dimension(:), allocatable :: valid

  end type userSurfCommType

  type userIntSurf

     character(len=maxStringLen) :: famName
     integer(Kind=intType) :: famID
     real(kind=realType), dimension(:, :), allocatable :: pts
     integer(kind=intType), dimension(:, :), allocatable :: conn

     ! Two separate commes: One for the nodes (based on the primal
     ! mesh) and one for the variables (based on the dual mesh)
     type(userSurfCommType) :: nodeComm, faceComm

  end type userIntSurf

  integer(kind=intType), parameter :: nUserIntSurfsMax=25
  type(userIntSurf), dimension(nUserIntSurfsMax), target :: userIntSurfs
  integer(kind=intTYpe) :: nUserIntSurfs=0

contains


  subroutine addIntegrationSurface(pts, conn, famName, famID, nPts, nConn)
    ! Add a user-supplied integration surface. 

    use communication, only : myID
    use constants

    implicit none

    ! Input variables
    real(kind=realType), dimension(3, nPts), intent(in) :: pts
    integer(kind=intType), dimension(4, nConn), intent(in) :: conn
    integer(kind=intType), intent(in) :: nPts, nConn, famID
    character(len=*) :: famName
    type(userIntSurf), pointer :: surf

    ! Not really much to do here...we just have to save the data
    ! into the data structure untilly we actual have to do the
    ! search.
    nUserIntSurfs = nUserIntSurfs + 1
    if (nUserIntSurfs > nUserIntSurfsMax) then 
       print *,"Error: Exceeded the maximum number of user-supplied "&
            &"integration slices. Increase nUserIntSurfsMax"
       stop
    end if

    if (myid == 0) then 
       surf => userIntSurfs(nUserIntSurfs)

       allocate(surf%pts(3, nPts), surf%conn(4, nConn))
       surf%pts = pts
       surf%conn = conn
       surf%famName = famName
       surf%famID = famID
    end if

  end subroutine addIntegrationSurface

  subroutine buildVolumeADTs(oBlocks, useDual)

    ! This builds volume ADTs for the the owned blocks. It will build
    ! either the dual mesh or the primal mesh depening on the flag
    ! useDual. 

    use constants
    use overset, only : oversetBlock
    use blockPointers, only : nDom, x, ie, je, ke, il, jl, kl,  vol, ib, jb, kb
    use adtBuild, only : buildSerialHex
    use utils, only : setPointers, EChk
    implicit none

    ! Input/Output Parameters
    type(oversetBlock), dimension(:), target :: oBlocks
    logical :: useDual

    ! Working Parameters
    integer(kind=intType) :: nInterpol, nn, i, j, k, iii, jjj, kkk
    integer(kind=intType) :: ii, jj, kk, mm, nADT, nHexa, planeOffset
    type(oversetBlock), pointer :: oBlock

    nInterpol = 1 ! we get the ADT to compute the interpolated volume for us. 

    domainLoop: do nn=1, nDom

       call setPointers(nn, 1, 1)
       oBlock => oBlocks(nn)

       primalOrDual: if (useDual) then 

          ! Now setup the data for the ADT
          nHexa = il * jl * kl
          nADT = ie * je * ke
          oBlock%il = il
          oBlock%jl = jl
          oBlock%kl = kl

          allocate(oBlock%xADT(3, nADT), oBlock%hexaConn(8, nHexa), &
               oBlock%qualDonor(1, nADT))
          ! Fill up the xADT using cell centers (dual mesh)
          mm = 0
          do k=1, ke
             do j=1, je
                do i=1, ie
                   mm = mm + 1
                   oBlock%xADT(:, mm) = eighth*(&
                        x(i-1, j-1, k-1, :) + &
                        x(i  , j-1, k-1, :) + &
                        x(i-1, j  , k-1, :) + &
                        x(i  , j  , k-1, :) + &
                        x(i-1, j-1, k  , :) + &
                        x(i  , j-1, k  , :) + &
                        x(i-1, j  , k  , :) + &
                        x(i  , j  , k  , :))
                   oBlock%qualDonor(1, mm) = vol(i, j, k)
                end do
             end do
          end do

          mm = 0
          ! These are the 'elements' of the dual mesh.
          planeOffset = ie * je
          do k=2, ke
             do j=2, je
                do i=2, ie
                   mm = mm + 1
                   oBlock%hexaConn(1, mm) = (k-2)*planeOffset + (j-2)*ie + (i-2) + 1
                   oBlock%hexaConn(2, mm) = oBlock%hexaConn(1, mm) + 1 
                   oBlock%hexaConn(3, mm) = oBlock%hexaConn(2, mm) + ie
                   oBlock%hexaConn(4, mm) = oBlock%hexaConn(3, mm) - 1 

                   oBlock%hexaConn(5, mm) = oBlock%hexaConn(1, mm) + planeOffset
                   oBlock%hexaConn(6, mm) = oBlock%hexaConn(2, mm) + planeOffset
                   oBlock%hexaConn(7, mm) = oBlock%hexaConn(3, mm) + planeOffset
                   oBlock%hexaConn(8, mm) = oBlock%hexaConn(4, mm) + planeOffset
                end do
             end do
          end do
       else
          ! Note that we will be including the halo primal cells. This
          ! should slightly increase robusness for viscous off-wall
          ! spacing. This means the primal mesh has 1 MORE node/cell
          ! in each direction. 

          ! Now setup the data for the ADT
          nHexa = ie * je * ke
          nADT = ib * jb * kb
          oBlock%il = ie
          oBlock%jl = je
          oBlock%kl = ke

          allocate(oBlock%xADT(3, nADT), oBlock%hexaConn(8, nHexa), &
               oBlock%qualDonor(1, nADT))

          oBlock%qualDonor = zero
          ! Fill up the xADT using the primal nodes
          mm = 0
          do k=0, ke
             do j=0, je
                do i=0, ie
                   mm = mm + 1
                   oBlock%xADT(:, mm) = x(i, j, k, :)

                   ! Since we don't have all 8 volumes surrounding the
                   ! halo nodes, clip the volumes to be between 0 and ib etc.
                   do iii=0,1
                      do jjj=0,1
                         do kkk=0,1
                            ii = min(max(0, iii+i), ib)
                            jj = min(max(0, jjj+j), jb)
                            kk = min(max(0, kkk+k), kb)

                            oBlock%qualDonor(1, mm) = oBlock%qualDonor(1, mm) + &
                                 vol(ii, jj, kk)
                         end do
                      end do
                   end do

                   ! Dividing by 8 isn't strictly necessary but we'll
                   ! do it anyway.
                   oBlock%qualDonor(1, mm) = oBlock%qualDonor(1, mm) * eighth
                end do
             end do
          end do

          mm = 0
          ! These are the 'elements' of the dual mesh.
          planeOffset = ib * jb
          do k=1, ke
             do j=1, je
                do i=1, ie
                   mm = mm + 1
                   oBlock%hexaConn(1, mm) = (k-1)*planeOffset + (j-1)*ib + (i-1) + 1
                   oBlock%hexaConn(2, mm) = oBlock%hexaConn(1, mm) + 1 
                   oBlock%hexaConn(3, mm) = oBlock%hexaConn(2, mm) + ib
                   oBlock%hexaConn(4, mm) = oBlock%hexaConn(3, mm) - 1 

                   oBlock%hexaConn(5, mm) = oBlock%hexaConn(1, mm) + planeOffset
                   oBlock%hexaConn(6, mm) = oBlock%hexaConn(2, mm) + planeOffset
                   oBlock%hexaConn(7, mm) = oBlock%hexaConn(3, mm) + planeOffset
                   oBlock%hexaConn(8, mm) = oBlock%hexaConn(4, mm) + planeOffset
                end do
             end do
          end do
       end if primalOrDual

       ! Call the custom build routine -- Serial only, only Hexa volumes,
       ! we supply our own ADT Type

       call buildSerialHex(nHexa, nADT, oBlock%xADT, oBlock%hexaConn, oBlock%ADT)
    end do domainLoop

  end subroutine buildVolumeADTs

  subroutine performInterpolation(pts, oBlocks, useDual, comm)

    ! This routine performs the actual searches for the slices. It is
    ! generic in the sense that it will search an arbtitrary set of
    ! points on either the primal or dual meshes. The final required
    ! communication data is then written into the supplied comm. 

    use constants
    use block, only : interpPtType
    use communication, only : adflow_comm_world, myid, nProc
    use overset, only : oversetBlock
    use blockPointers, only : nDom, x, ie, je, ke, il, jl, kl, x, iBlank, vol
    use adtLocalSearch, only :  mindistancetreesearchsinglepoint, &
         containmenttreesearchsinglepoint
    use adtUtils, only : stack
    use adtData, only : adtBBoxTargetType
    use utils, only : setPointers, mynorm2, EChk
    use inputOverset, only : oversetProjTol
    use oversetUtilities, only : fracToWeights2, getCumulativeForm

    implicit none

    ! Input parameters:
    real(kind=realType), dimension(:, :), intent(in) :: pts
    type(userIntSurf) :: surf
    type(oversetBlock), dimension(:), target, intent(in) :: oBlocks
    logical, intent(in) :: useDual
    !type(oversetSurf), dimension(:), intent(in) :: oSurfs
    type(userSurfCommType) :: comm

    ! Working parameters
    type(oversetBlock), pointer :: oBlock
    type(interpPtType), dimension(:), allocatable :: surfFringes

    integer(Kind=intType) :: i, j, k, ii, jj, kk, iii, jjj, kkk, nn, mm 
    integer(kind=intType) :: iSurf, ierr, nInterpol, iProc
    integer(kind=intType) :: nHexa, nAdt, planeOffset, elemID, nPts
    real(kind=realType) :: xc(4), weight(8)
    integer mpiStatus(MPI_STATUS_SIZE) 

    real(kind=realType) :: uvw(5), uvw2(5), donorQual, xcheck(3)
    integer(kind=intType) :: intInfo(3), intInfo2(3)
    logical :: failed, invalid

    integer(kind=intType), dimension(:, :), allocatable :: donorInfo, intSend
    integer(kind=intType), dimension(:), allocatable :: procSizes
    real(kind=realType), dimension(:, :), allocatable :: donorFrac, realSend

    ! Variables we have to pass the ADT search routine
    integer(kind=intType), dimension(:), pointer :: BB
    type(adtBBoxTargetType), dimension(:), pointer :: BB2
    integer(kind=intType), dimension(:), pointer :: frontLeaves
    integer(kind=intType), dimension(:), pointer :: frontLeavesNew

    ! Data for the search
    allocate(BB(20), frontLeaves(25), frontLeavesNew(25), stack(100))

    nPts = size(pts, 2)

    ! Allocate donor information arrays
    allocate(donorFrac(4, nPts), donorInfo(5, nPts))
    donorInfo = -1
    donorFrac(4, :) = large
    nInterpol = 1
    domainSearch: do nn=1, nDom
       oBlock => oBlocks(nn)
       call setPointers(nn, 1, 1)

       ! Search the supplied pts one at a time
       elemLoop: do i=1, nPts

          xc(1:3) = pts(:, i)

          ! Call the standard tree search
          call containmentTreeSearchSinglePoint(oBlock%ADT, xc, intInfo, uvw, &
               oBlock%qualDonor, nInterpol, BB, frontLeaves, frontLeavesNew, failed)

          ! Make sure this point is not garbage.
          if (intInfo(1) >= 0) then 
             call fracToWeights2(uvw(1:3), weight)
             xcheck = zero
             do j=1,8
                xcheck = xcheck + weight(j)*oBlock%xADT(:, oBlock%hexaConn(j, intInfo(3)))
             end do

             if (mynorm2(xcheck - xc(1:3)) > oversetProjTol) then 
                failed = .True.
             end if
          end if

          if (intInfo(1) >= 0 .and. failed) then 
             ! we "found" a point but it is garbage. Do the failsafe search
             xc(4) = large
             call minDistanceTreeSearchSinglePoint(oBlock%ADT, xc, intInfo, uvw, &
                  oBlock%qualDonor, nInterpol, BB2, frontLeaves, frontLeavesNew)

             ! Check this one:
             call fracToWeights2(uvw(1:3), weight)
             xcheck = zero
             do j=1,8
                xcheck = xcheck + weight(j)*oBlock%xADT(:, oBlock%hexaConn(j, intInfo(3)))
             end do

             ! Since this is the last line of defence, relax the tolerance a bit
             if (mynorm2(xcheck - xc(1:3)) > 100*oversetProjTol) then 
                ! This fringe has not found a donor
                intInfo(1) = -1
             else
                ! This one has now passed.

                ! Important! uvw(4) is the distance squared for this search
                ! not
                uvw(4) = uvw(5)
             end if
          end if

          elemFound: if (intInfo(1) >= 0) then 

             ! Donor and block and index information for this donor. 
             donorQual = uvw(4)
             elemID = intInfo(3) - 1 ! Make it zero based

             ! The dual mesh needs an offset of 1 becuase it only used
             ! 1:ie values. This is not necessary for the primal. 
             if (useDual) then 
                ii = mod(elemID, oBlock%il) + 1
                jj = mod(elemID/oBlock%il, oBlock%jl) + 1
                kk = elemID/(oBlock%il*oBlock%jl) + 1
             else
                ii = mod(elemID, oBlock%il) 
                jj = mod(elemID/oBlock%il, oBlock%jl) 
                kk = elemID/(oBlock%il*oBlock%jl) 
             end if
             ! Rememebr donorFrac(4, i) is the current best quality 
             if ( donorQual < donorFrac(4, i)) then 

                invalid = .False.

                ! For the dual mesh search, we have to make sure the
                ! potential donors are valid. Such a check is not
                ! necessary for the primal search since all nodes are
                ! considered valid. 

                if (useDual) then 
                   ! Check if the point is invalid. We can do this
                   ! with i-blank array. We can only accept compute
                   ! cells (iblank=1) or interpolated
                   ! cells(iblank=-1).
                   do kkk=0,1
                      do jjj=0,1
                         do iii=0,1
                            if (.not. (iblank(ii+iii, jj+jjj, kk+kkk) == 1 .or. &
                                 iblank(ii+iii, jj+jjj, kk+kkk) == -1)) then 
                               invalid = .True.
                            end if
                         end do
                      end do
                   end do
                end if

                if (.not. invalid) then 

                   ! Set the quality of the donor to the one we
                   ! just found. Save the rest of the necessary
                   ! information. 
                   donorInfo(1, i) = myid
                   donorInfo(2, i) = nn
                   donorInfo(3, i) = ii
                   donorInfo(4, i) = jj
                   donorInfo(5, i) = kk
                   donorFrac(1:3, i) = uvw(1:3)
                   donorFrac(4, i) = donorQual
                end if
             end if
          end if elemFound
       end do elemLoop
    end do domainSearch

    ! Next count up the number of valid donors we've found and compact
    ! the info back to that length.
    if (myid /=0) then 
       j = 0
       do i=1, nPts
          if (donorInfo(1, i) /= -1) then 
             j = j + 1
          end if
       end do
       allocate(intSend(6, j), realSend(4, i))
       if (j > 0) then 
          j = 0
          do i=1, nPts
             if (donorInfo(1, i) /= -1) then 
                j = j + 1
                intSend(1:5, j) = donorInfo(:, i)
                intSend(6, j) = i
                realSend(:, j) = donorFrac(:, i)
             end if
          end do
       end if
    else
       ! On the root proc, use intSend and realSend as the receiver
       ! buffer. These can be at most nPts sized.
       allocate(intSend(6, nPts), realSend(4, nPts))
    end if

    ! Gather up the sizes (j) to the root processor so he know who to
    ! expect data from.
    allocate(procSizes(0:nProc-1))

    call mpi_gather(j, 1, adflow_integer, procSizes, 1, &
         adflow_integer, 0, adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Next all the procs need to send all the information back to the
    ! root processor where we will determine the proper donors for
    ! each of the cells

    ! All procs except root fire off their data.
    if (myid >= 1) then 
       if (j > 0) then 
          call mpi_send(intSend, j*6, adflow_integer, 0, myid, &
               adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          call mpi_send(realSend, j*4, adflow_real, 0, myid, &
               adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)
       end if
    end if

    ! And the root processor recieves it...
    if (myid == 0) then 
       do iProc=1, nProc-1
          ! Determine if this proc has sent anything:
          if (procSizes(iProc) /= 0) then 

             call MPI_recv(intSend, 6*nPts, adflow_integer, iProc, iProc,&
                  adflow_comm_world, mpiStatus, ierr)
             call EChk(ierr,__FILE__,__LINE__)

             call MPI_recv(realSend, 4*nPts, adflow_real, iProc, iProc,&
                  adflow_comm_world, mpiStatus, ierr)
             call EChk(ierr,__FILE__,__LINE__)

             ! Now process the data (intSend and realSend) that we
             ! just received. We don't need to check the status for
             ! the sizes becuase we already know the sizes from the
             ! initial gather we did. 

             do i=1, procSizes(iProc)
                ii = intSend(6, i)

                if (realSend(4, i) < donorFrac(4, ii)) then 
                   ! The incoming quality is better. Accept it. 
                   donorInfo(1:5, ii) = intSend(1:5, i)
                   donorFrac(:, ii) = realSend(:, i)
                end if
             end do
          end if
       end do

       ! To make this easier, convert the information we have to a
       ! 'fringeType' array so we can use the pre-existing sorting
       ! routine.
       allocate(surfFringes(nPts))
       do i=1, nPts
          surfFringes(i)%donorProc = donorInfo(1, i)
          surfFringes(i)%donorBlock = donorInfo(2, i)
          surfFringes(i)%dI = donorInfo(3, i)
          surfFringes(i)%dJ = donorInfo(4, i)
          surfFringes(i)%dK = donorInfo(5, i)
          surfFringes(i)%donorFrac = donorFrac(1:3, i)
          ! Use the myBlock attribute to keep track of the original
          ! index. When we sort the fringes, they will no longer be
          ! in the same order
          surfFringes(i)%myBlock = i
       end do

       ! Perform the actual sort. 
       call qsortInterpPtType(surfFringes, nPts)

       ! We will reuse-proc sizes to now mean the number of elements
       ! that the processor *actually* has to send. We will include
       ! the root proc itself in the calc becuase that will tell us
       ! the size of the internal comm structure. 

       procSizes = 0
       allocate(comm%valid(nPts))
       comm%valid = .True.
       do i=1, nPts
          if (surfFringes(i)%donorProc < 0) then 
             ! We dont have a donor. Flag this point as invalid
             comm%valid(i) = .False.
          end if

          ! Dump the points without donors on the root proc by making
          ! sure j is at least 0 for the root proc. These will just
          ! simply be ignored during the comm. 
          j = max(surfFringes(i)%donorProc, 0)
          procSizes(j) = procSizes(j) + 1
       end do
    end if

    ! Simply broadcast out the the proc sizes back to everyone so all
    ! processors know if they are to receive anything back. 
    call mpi_bcast(procSizes, nProc, adflow_Integer, 0, adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! We can now save some of the final comm data required on the
    ! root proc for this surf. 
    allocate(comm%procSizes(0:nProc-1), comm%procDisps(0:nProc))

    ! Copy over procSizes and generate the cumulative form of
    ! the size array, procDisps
    comm%procSizes = procSizes
    call getCumulativeForm(comm%procSizes, nProc, comm%procDisps)

    ! Record the elemInverse which is necessary to index into
    ! the original conn array. 
    if (myid == 0) then 
       allocate(comm%Inv(nPts))
       do i=1, nPts
          comm%Inv(i) = surfFringes(i)%myBlock
       end do
    end if

    ! Now we can send out the final donor information to the
    ! processors that must supply it. 
    comm%nDonor = procSizes(myID)
    allocate(comm%frac(3, comm%nDonor), comm%donorInfo(4, comm%nDonor))

    if (myid >= 1) then 
       if (comm%nDonor > 0) then 
          ! We are responible for at least 1 donor. We have to make
          ! use of the intSend and realSend buffers again (which
          ! are guaranteed to be big enough). The reason we can't
          ! dump the data in directlyis that intSend and realSend
          ! have a different leading index than we need on the
          ! final data structure. 

          call MPI_recv(intSend, 6*comm%nDonor, adflow_integer, 0, myid, &
               adflow_comm_world, mpiStatus, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          call MPI_recv(realSend, 4*comm%nDonor, adflow_real, 0, myID, &
               adflow_comm_world, mpiStatus, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Copy into final structure
          do i=1, comm%nDonor
             comm%donorInfo(:, i) = intSend(1:4, i)
             comm%frac(:, i) = realSend(1:3, i)
          end do
       end if
    else
       ! We are the root processor.
       if (comm%nDonor > 0) then 
          ! We need to copy out our donor info on the root proc if we have any
          do i= comm%procDisps(myID)+1, comm%procDisps(myID+1)
             comm%donorInfo(1, i) = surfFringes(i)%donorBlock
             comm%donorInfo(2, i) = surfFringes(i)%dI
             comm%donorInfo(3, i) = surfFringes(i)%dJ
             comm%donorInfo(4, i) = surfFringes(i)%dK
             comm%frac(1:3, i)    = surfFringes(i)%donorFrac
          end do
       end if

       ! Now loop over the rest of the procs and send out the info we
       ! need. We have to temporarily copy the data back out of
       ! fringes to the intSend and realSend arrays
       do iProc=1, nProc-1

          if (comm%procSizes(iProc) > 0) then 
             ! Have something to send here:
             j = 0
             do i=comm%procDisps(iProc)+1, comm%procDisps(iProc+1)
                j = j + 1

                intSend(1, j) = surfFringes(i)%donorBlock
                intSend(2, j) = surfFringes(i)%dI      
                intSend(3, j) = surfFringes(i)%dJ
                intSend(4, j) = surfFringes(i)%dK
                realSend(1:3, j) = surfFringes(i)%donorFrac
             end do

             call mpi_send(intSend, j*6, adflow_integer, iProc, iProc, &
                  adflow_comm_world, ierr)
             call EChk(ierr,__FILE__,__LINE__)

             call mpi_send(realSend, j*4, adflow_real, iProc, iProc, &
                  adflow_comm_world, ierr)
             call EChk(ierr,__FILE__,__LINE__)
          end if
       end do

       ! Deallocate data allocatd only on root proc
       deallocate(surfFringes)
    end if

    ! Nuke rest of allocated on all procs
    deallocate(intSend, realSend, procSizes, donorInfo, donorFrac)

  end subroutine performInterpolation

  subroutine interpolateIntegrationSurfaces

    ! This routine performs the actual searches for the slices. We
    ! reuse much of the same machinery as is used in the overset code.

    use constants
    use communication, only : adflow_comm_world, myid
    use overset, only : oversetBlock
    use blockPointers, only : nDom, ie, je, ke, il, jl, kl
    use adtBuild, only : buildSerialHex, destroySerialHex
    use utils, only : setPointers, EChk

    implicit none

    ! Working parameters
    type(oversetBlock), dimension(nDom), target :: oBlocks
    type(userIntSurf), pointer :: surf

    integer(Kind=intType) :: iSurf, ii, i, nn, nPts, ierr
    real(kind=realType), dimension(:, :), allocatable :: pts
    logical :: useDual
    if (nUserIntSurfs == 0) then 
       return
    end if
    primalDualLoop: do ii=1, 2
       if (ii==1) then 
          useDual = .True.
       else
          useDual = .False.
       end if

       call buildVolumeADTs(oBlocks, useDual)

       masterLoop: do iSurf=1, nUserIntSurfs

          surf => userIntSurfs(iSurf)

          if (myid == 0) then 

             if (ii==1) then 
                nPts = size(surf%conn, 2)
                allocate(pts(3, nPts))
                ! Use the dual (cell center values)
                elemLoop: do i=1, nPts

                   ! Compute the center of the cell. 
                   pts(:, i) = fourth*( & 
                        surf%pts(:, surf%conn(1, i)) + &
                        surf%pts(:, surf%conn(2, i)) + &
                        surf%pts(:, surf%conn(3, i)) + &
                        surf%pts(:, surf%conn(4, i)))
                end do elemLoop
             else if (ii==2) then 
                nPts = size(surf%pts, 2)
                allocate(pts(3, nPts))

                ! Use the primal nodal values
                nodeLoop: do i=1, nPts
                   ! Compute the center of the cell. 
                   pts(:, i) = surf%pts(:, i)
                end do nodeLoop
             end if
          end if

          ! Send the number of points back to all procs:
          call mpi_bcast(nPts, 1, adflow_integer, 0, adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! All other procs except the root allocate space and receive
          ! the pt array. 
          if (myid /= 0) then 
             allocate(pts(3, nPts))
          end if

          call mpi_bcast(pts, 3*nPts, adflow_real, 0, adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Call the actual interpolation routine
          if (ii==1) then 
             call performInterpolation(pts, oBlocks, .True., surf%faceComm)
          else
             call performInterpolation(pts, oBlocks, .False., surf%nodeComm)
          end if

          deallocate(pts)
       end do masterLoop

       ! Destroy the ADT Data and allocated values
       do nn=1, nDom
          call destroySerialHex(oBlocks(nn)%ADT)
          deallocate(oBlocks(nn)%xADT, oBlocks(nn)%hexaConn, oBlocks(nn)%qualDonor)
       end do
    end do primalDualLoop
  end subroutine interpolateIntegrationSurfaces

  subroutine commUserIntegrationSurfaceVars(recvBuffer, varStart, varEnd, comm)

    use constants
    use block, onlY : flowDoms, nDom
    use communication, only : myid, adflow_comm_world
    use utils, only : EChk
    use oversetUtilities, only :fracToWeights

    implicit none

    ! Input/Output
    real(kind=realType), dimension(:) :: recvBuffer
    integer(kind=intType), intent(in) :: varStart, varEnd
    type(userSurfCommType) :: comm

    ! Working
    real(kind=realType), dimension(:), allocatable :: sendBuffer
    integer(Kind=intType) :: d1, i1, j1, k1, jj, k, nvar, i, ierr
    real(kind=realType), dimension(8) :: weight

    ! The number of variables we are transferring:
    nVar = varEnd - varStart + 1

    ! We assume that the pointers to the realCommVars have already been set. 

    allocate(sendBuffer(nVar*comm%nDonor))

    ! First generate the interpolated data necessary
    jj = 0
    donorLoop: do i=1, comm%nDonor
       ! Convert the frac to weights
       call fracToWeights(comm%frac(:, i), weight)

       ! Block and indices for easier reading. The +1 is due to the
       ! pointer offset on realCommVars.

       d1 = comm%donorInfo(1, i) ! Block Index
       i1 = comm%donorInfo(2, i)+1 ! donor I index
       j1 = comm%donorInfo(3, i)+1 ! donor J index
       k1 = comm%donorInfo(4, i)+1 ! donor K index

       ! We are interpolating nVar variables
       do k=varStart, varEnd
          jj = jj + 1
          if (d1 > 0) then ! Is this pt valid?
             sendBuffer(jj) = &
                  weight(1)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1+1,k1+1)
          end if
       end do
    end do donorLoop

    ! Now we can do an mpi_gatherv to the root proc:
    call mpi_gatherv(sendBuffer, nVar*comm%nDonor, adflow_real, recvBuffer, &
         nVar*comm%procSizes, nVar*comm%procDisps, adflow_real, 0, adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)
    deallocate(sendBuffer)

  end subroutine commUserIntegrationSurfaceVars

  subroutine integrateUserSurfaces(localValues, famList, sps)

    use constants
    use block, onlY : flowDoms, nDom
    use flowVarRefState, only : pRef, rhoRef, pRef, timeRef, LRef, TRef
    use communication, only : myid, adflow_comm_world
    use utils, only : EChk, mynorm2
    use flowUtils, only : computePtot, computeTtot
    use sorting, only : bsearchIntegers
    implicit none

    ! Input Parameters
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues
    integer(kind=intType), dimension(:), intent(in) :: famList
    integer(kind=intType), intent(in) :: sps

    ! Working parameters
    integer(kind=intType) :: iSurf, i, j, k, jj, ierr, nn, iDim
    real(kind=realType), dimension(:), allocatable :: recvBuffer1, recvBuffer2
    type(userIntSurf), pointer :: surf
    real(kind=realType) :: mReDim, rho, vx, vy, vz, PP, massFLowRate, vn
    real(Kind=realType) :: pTot, Ttot, mass_Ptot, mass_Ttot, mass_ps, massFlowRateLocal
    real(kind=realType), dimension(3) :: v1, v2, x1, x2, x3, x4, n
    logical :: valid
    logical, dimension(:), allocatable :: ptValid
    ! Set the pointers for the required communication variables: rho,
    ! vx, vy, vz, P and the x-coordinates

    domainLoop:do nn=1, nDom
       flowDoms(nn, 1, sps)%realCommVars(1)%var => flowDoms(nn, 1, sps)%w(:, :, :, 1)
       flowDoms(nn, 1, sps)%realCommVars(2)%var => flowDoms(nn, 1, sps)%w(:, :, :, 2)
       flowDoms(nn, 1, sps)%realCommVars(3)%var => flowDoms(nn, 1, sps)%w(:, :, :, 3)
       flowDoms(nn, 1, sps)%realCommVars(4)%var => flowDoms(nn, 1, sps)%w(:, :, :, 4)
       flowDoms(nn, 1, sps)%realCommVars(5)%var => flowDoms(nn, 1, sps)%P(:, :, :)
       flowDoms(nn, 1, sps)%realCommVars(6)%var => flowDoms(nn, 1, sps)%x(:, :, :, 1)
       flowDoms(nn, 1, sps)%realCommVars(7)%var => flowDoms(nn, 1, sps)%x(:, :, :, 2)
       flowDoms(nn, 1, sps)%realCommVars(8)%var => flowDoms(nn, 1, sps)%x(:, :, :, 3)
    end do domainLoop

    masterLoop: do iSurf=1, nUserIntSurfs

       ! Pointer for easier reading
       surf => userIntSurfs(iSurf)

       ! Do we need to include this surface?
       famInclude: if (bsearchIntegers(surf%famID, famList) > 0) then

          ! Communicate the face values and the nodal values
          if (myid == 0) then 
             allocate(recvBuffer1(5*size(surf%conn, 2)))
             allocate(recvBuffer2(3*size(surf%pts, 2)))
          else
             allocate(recvBuffer1(0))
             allocate(recvBuffer2(0))
          end if

          call commUserIntegrationSurfaceVars(recvBuffer1, 1, 5, surf%faceComm)
          call commUserIntegrationSurfaceVars(recvBuffer2, 6, 8, surf%nodeComm)

          ! *Finally* we can do the actual integrations
          if (myid == 0)  then 
             allocate(ptValid(size(surf%pts, 2)))

             ! Before we do the integration, update the pts array from
             ! the values in recvBuffer2
             do i=1, size(surf%pts, 2)
                j = surf%nodeComm%inv(i)
                surf%pts(:, j) = recvBuffer2(3*i-2:3*i)
                ptValid(j) = surf%nodeComm%valid(i)
             end do

             massFlowRate = zero
             mass_Ptot = zero
             mass_Ttot = zero
             mass_Ps = zero

             mReDim = sqrt(pRef*rhoRef)

             elemLoop: do i=1, size(surf%conn, 2)

                ! First check if this cell is valid to integrate. Both
                ! the cell center *and* all nodes must be valid

                ! j is now the element index into the original conn array. 
                j = surf%faceComm%Inv(i)

                valid = surf%faceComm%valid(i) .and. &
                     ptValid(surf%conn(1, j)) .and. ptValid(surf%conn(2, j)) .and. &
                     ptValid(surf%conn(3, j)) .and. ptValid(surf%conn(4, j))
                if (valid) then 
                   ! Extract out the interpolated quantities
                   rho = recvBuffer1(5*i-4)
                   vx  = recvBuffer1(5*i-3)
                   vy  = recvBuffer1(5*i-2)
                   vz  = recvBuffer1(5*i-1)
                   PP  = recvBuffer1(5*i  )

                   call computePtot(rho, vx, vy, vz, pp, Ptot)
                   call computeTtot(rho, vx, vy, vz, pp, Ttot)

                   ! Get the coordinates of the quad
                   x1 = surf%pts(:, surf%conn(1, j))
                   x2 = surf%pts(:, surf%conn(2, j))
                   x3 = surf%pts(:, surf%conn(3, j))
                   x4 = surf%pts(:, surf%conn(4, j))

                   ! Diagonal vectors
                   v1 = x3 - x1
                   v2 = x4 - x2

                   ! Take the cross product of the two diagaonal vectors.
                   n(1) = half*(v1(2)*v2(3) - v1(3)*v2(2))
                   n(2) = half*(v1(3)*v2(1) - v1(1)*v2(3))
                   n(3) = half*(v1(1)*v2(2) - v1(2)*v2(1))

                   ! Normal face velocity
                   vn = vx*n(1) + vy*n(2) + vz*n(3)

                   ! Finally the mass flow rate
                   massFlowRateLocal = rho*vn*mReDim
                   massFlowRate = massFlowRate + massFlowRateLocal
                   mass_Ptot = mass_pTot + Ptot * massFlowRateLocal * Pref
                   mass_Ttot = mass_Ttot + Ttot * massFlowRateLocal * Tref
                   mass_Ps = mass_Ps + pp*massFlowRateLocal*Pref
                end if

             end do elemLoop
             deallocate(ptValid)
          end if

          ! Accumulate the final values
          localValues(iMassFlow) = localValues(iMassFlow) + massFlowRate
          localValues(iMassPtot) = localValues(iMassPtot) + mass_Ptot
          localValues(iMassTtot) = localValues(iMassTtot) + mass_Ttot
          localValues(iMassPs)   = localValues(iMassPs)   + mass_Ps

          deallocate(recvBuffer1, recvBuffer2)


       end if famInclude
    end do masterLoop

  end subroutine integrateUserSurfaces

  subroutine qsortInterpPtType(arr, nn)

    use constants
    use block ! Cannot use-only becuase of <= operator
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn

    type(interpPtType), dimension(*), intent(inout) :: arr
    !
    !      Local variables.
    !
    integer(kind=intType), parameter :: m = 7

    integer(kind=intType) :: nStack
    integer(kind=intType) :: i, j, k, r, l, jStack, ii

    integer :: ierr

    type(interpPtType) :: a, tmp

    integer(kind=intType), allocatable, dimension(:) :: stack
    integer(kind=intType), allocatable, dimension(:) :: tmpStack

    ! Allocate the memory for stack.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                         &
         call terminate("qsortinterpPt", &
         "Memory allocation failure for stack")

    ! Initialize the variables that control the sorting.

    jStack = 0
    l      = 1
    r      = nn

    ! Start of the algorithm

    do

       ! Check for the size of the subarray.

       if((r-l) < m) then

          ! Perform insertion sort

          do j=l+1,r
             a = arr(j)
             do i=(j-1),l,-1
                if(arr(i) <= a) exit
                arr(i+1) = arr(i)
             enddo
             arr(i+1) = a
          enddo

          ! In case there are no more elements on the stack, exit from
          ! the outermost do-loop. Algorithm has finished.

          if(jStack == 0) exit

          ! Pop stack and begin a new round of partitioning.

          r = stack(jStack)
          l = stack(jStack-1)
          jStack = jStack - 2

       else

          ! Subarray is larger than the threshold for a linear sort.
          ! Choose median of left, center and right elements as
          ! partitioning element a.
          ! Also rearrange so that (l) <= (l+1) <= (r).

          k = (l+r)/2
          tmp      = arr(k)             ! Swap the elements
          arr(k)   = arr(l+1)           ! k and l+1.
          arr(l+1) = tmp

          if(arr(r) < arr(l)) then
             tmp    = arr(l)             ! Swap the elements
             arr(l) = arr(r)             ! r and l.
             arr(r) = tmp
          endif

          if(arr(r) < arr(l+1)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(r)           ! r and l+1.
             arr(r)   = tmp
          endif

          if(arr(l+1) < arr(l)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(l)           ! l and l+1.
             arr(l)   = tmp
          endif

          ! Initialize the pointers for partitioning.

          i = l+1
          j = r
          a = arr(l+1)

          ! The innermost loop

          do

             ! Scan up to find element >= a.
             do
                i = i+1
                if(a <= arr(i)) exit
             enddo

             ! Scan down to find element <= a.
             do
                j = j-1
                if(arr(j) <= a) exit
             enddo

             ! Exit the loop in case the pointers i and j crossed.

             if(j < i) exit

             ! Swap the element i and j.

             tmp    = arr(i)
             arr(i) = arr(j)
             arr(j) = tmp
          enddo

          ! Swap the entries j and l+1. Remember that a equals
          ! arr(l+1).

          arr(l+1) = arr(j)
          arr(j)   = a

          ! Push pointers to larger subarray on stack,
          ! process smaller subarray immediately.

          jStack = jStack + 2
          if(jStack > nStack) then

             ! Storage of the stack is too small. Reallocate.

             allocate(tmpStack(nStack), stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortinterpPt", &
                  "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortinterpPt", &
                  "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortinterpPt", &
                  "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortinterpPt", &
                  "Deallocation error for tmpStack")
          endif

          if((r-i+1) >= (j-l)) then
             stack(jStack)   = r
             r               = j-1
             stack(jStack-1) = j
          else
             stack(jStack)   = j-1
             stack(jStack-1) = l
             l               = j
          endif

       endif
    enddo

    ! Release the memory of stack.

    deallocate(stack, stat=ierr)
    if(ierr /= 0)                         &
         call terminate("qsortinterpPt", &
         "Deallocation error for stack")

    ! Check in debug mode whether the array is really sorted.

    if( debug ) then
       do i=1,(nn-1)
          if(arr(i+1) < arr(i))                 &
               call terminate("qsortinterpPt", &
               "Array is not sorted correctly")
       enddo
    endif

  end subroutine qsortInterpPtType
end module userSurfaceIntegrations
