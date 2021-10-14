module tecplotIO

  use constants, only : realType, intType, maxStringLen, maxCGNSNameLen
  implicit none
  save

  type slice

     ! nNodes : Number of nodes for this slice
     ! ind(2, nNodes) : Indices of the nodes in the global node list on either
     !                  side of node
     ! w(2, nNodes) : Weights used to multiply the two global nodes defined in
     !                ind to get compute nodal values (positions, forces etc)
     ! pL, vL, pD, vD : Pressure and viscous lift, pressure and viscous drag
     ! CLp, CLv, CDp, CDv : Coefficients of pressure and viscous lift and drag
     ! chord: chord of section

     character(len=maxStringLen) :: sliceName
     integer(kind=intType) :: sps
     integer(kind=intType), dimension(:,:), allocatable :: ind, conn
     real(kind=realType), dimension(:, :), allocatable :: w, vars
     integer(kind=intType) :: nNodes
     real(kind=realType) :: pL, vL, pD, vD, CLp, CLv, CDp, CDv
     real(kind=realType) :: chord, twist, thickness
     real(kind=realType), dimension(3) :: le, te
     real(kind=realType), dimension(3) :: pt, dir
     integer(kind=intType), allocatable, dimension(:) :: famList
  end type slice

  type liftDist
     ! nSegments: Number of nodes to use for distribution
     ! dir_ind: Index of direction..1 for x, 2 for y, 3 for z
     ! dir: Slice direction
     ! distName: Name of lift distribution
     ! slices: The list of slices this distribution will use
     ! delta: The current delta spacing for the distribution
     ! slicePoints: The list of points where the slices are taken
     character(len=maxStringLen) :: distName
     integer(kind=intType) :: nSegments, dir_ind
     integer(kind=intType), dimension(:), allocatable :: famList
     real(kind=realType) :: dir(3)
     real(kind=realType) :: delta
     real(kind=realType), dimension(:,:), allocatable :: slicePts
  end type liftDist

  logical :: liftDistInitialized = .False.
  integer(kind=intType) :: msCon1(16, 5), msCon2(4, 2)

   ! Data for the user supplied slices:
  integer(kind=intType), parameter :: nSliceMax=1000
  integer(kind=intType) :: nParaSlices=0
  integer(kind=intType) :: nAbsSlices=0
  type(slice), dimension(:, :), allocatable :: paraSlices, absSlices

  ! Data for the user supplied lift distributions
  integer(kind=intType), parameter :: nLiftDistMax=100
  integer(kind=intType) :: nLiftDists=0
  type(liftDist), dimension(nLiftDistMax), target :: liftDists

  ! Tecplot Variable names of the data in the lift distribution data file:
  character(len=maxCGNSNameLen), dimension(:), allocatable :: liftDistName
  integer(kind=intType), parameter :: nLiftDistVar=18

contains
  subroutine addParaSlice(sliceName, pt, direction, famList, n)
    !
    !       This subroutine is intended to be called from python.
    !       This routine will add a parametric slice to the list of user
    !       supplied slices.
    use constants
    use communication
    use surfaceFamilies
    use surfaceUtils
    use inputTimeSpectral
    implicit none

    ! Input parameters
    character*(*), intent(in) :: sliceName
    real(kind=realType), dimension(3), intent(in) :: pt, direction
    integer(kind=intType), intent(in) :: famList(n), n

    ! Working
    integer(kind=intType) :: sps, sizeNode, sizeCell
    integer(kind=intType), dimension(:), pointer :: wallList
    real(kind=realType), dimension(:, :), allocatable :: pts
    integer(kind=intType), dimension(:, :), allocatable :: conn
    integer(kind=intType), dimension(:), allocatable :: elemFam, cgnsBlockID

    if (.not. allocated(paraSlices)) then
       allocate(paraSlices(nSliceMax, nTimeIntervalsSpectral))
    end if

    ! We have to add a slice for each spectral instance.
    do sps=1, nTimeIntervalsSpectral
       nParaSlices = nParaSlices + 1

       if (nParaSlices > nSliceMax) then
          print *,'Error: Exceeded the maximum number of slices. Increase nSliceMax'
          stop
       end if

       ! Slices are created on walls and walls only. Retrieve the
       ! points, connectivity and familyID of all the walls.
       wallList =>  BCFamGroups(iBCGroupWalls)%famList
       call getSurfaceSize(sizeNode, sizeCell, wallList, size(wallList), .True.)
       allocate(pts(3, sizeNode), conn(4, sizeCell), elemFam(sizeCell), cgnsBlockID(sizeCell))
       call getSurfaceConnectivity(conn, cgnsBlockID, sizeCell, wallList, size(wallList), .True.)
       call getSurfacePoints(pts, sizeNode, sps, wallList, size(wallList), .True.)
       call getSurfaceFamily(elemFam, sizeCell, wallList, size(wallList), .True.)

       ! Create actual slice
       call createSlice(pts, conn, elemFam, paraSlices(nParaSlices, sps), pt, direction, &
            sliceName, famList)

       ! Clean up memory.
       deallocate(pts, conn, elemFam)
    end do

  end subroutine addParaSlice

  subroutine addAbsSlice(sliceName, pt, direction, famList, n)
    !
    !       This subroutine is intended to be called from python.
    !       This routine will add an absolute slice to the list of user
    !       supplied slices.
    use constants
    use communication
    use surfaceFamilies
    use surfaceUtils
    use inputTimeSpectral
    implicit none

    ! Input parameters
    character*(*), intent(in) :: sliceName
    real(kind=realType), dimension(3), intent(in) :: pt, direction
    integer(kind=intType), intent(in) :: famList(n), n

    ! Working
    integer(kind=intType) :: sps, sizeNode, sizeCell
    integer(kind=intType), dimension(:), pointer :: wallList
    real(kind=realType), dimension(:, :), allocatable :: pts
    integer(kind=intType), dimension(:, :), allocatable :: conn
    integer(kind=intType), dimension(:), allocatable :: elemFam, cgnsBlockID

    if (.not. allocated(absSlices)) then
       allocate(absSlices(nSliceMax, nTimeIntervalsSpectral))
    end if

    ! do sps=1, nTimeIntervalsSpectral
    do sps=1, 1
       nAbsSlices = nAbsSlices + 1

       if (nAbsSlices > nSliceMax) then
          print *,'Error: Exceeded the maximum number of slices. Increase nSliceMax'
          stop
       end if

       wallList =>  BCFamGroups(iBCGroupWalls)%famList
       call getSurfaceSize(sizeNode, sizeCell, wallList, size(wallList), .True.)
       allocate(pts(3, sizeNode), conn(4, sizeCell), elemFam(sizeCell), cgnsBlockID(sizeCell))
       call getSurfaceConnectivity(conn, cgnsBlockID, sizeCell, wallList, size(wallList), .True.)
       call getSurfacePoints(pts, sizeNode, sps, wallList, size(wallList), .True.)
       call getSurfaceFamily(elemFam, sizeCell, wallList, size(wallList), .True.)
       call createSlice(pts, conn, elemFam, absSlices(nAbsSlices, sps), pt, direction, &
            sliceName, famList)

       ! Clean up memory.
       deallocate(pts, conn, elemFam)
    end do

  end subroutine addAbsSlice

  subroutine addLiftDistribution(nSegments, dir_vec, dir_ind, distName, famList, n)
    !
    !       This subroutine is intended to be called from python.
    !       This routine will add the description of a lift distribution

    use constants
    use communication
    use surfaceFamilies
    implicit none

    ! Input parameters
    character*(*), intent(in) :: distName
    integer(kind=intType), intent(in) :: nSegments
    real(kind=realType), dimension(3) :: dir_vec
    integer(kind=intType), intent(in) :: dir_ind
    integer(kind=intType), intent(in) :: famList(n), n

    nLiftDists = nLiftDists + 1
    if (nLiftDists > nLiftDistMax) then
       print *,'Error: Exceeded the maximum number of lift distributions. &
            &Increase nLiftDistMax'
       stop
    end if

    liftDists(nLIftDists)%nSegments = nSegments
    liftDists(nLiftDists)%dir = dir_vec
    liftDists(nLiftDists)%distName = distName
    liftDists(nLIftDists)%dir_ind = dir_ind

    allocate(liftDists(nLiftDists)%famList(n))
    liftDists(nLiftDists)%famList(:) = famList

  end subroutine addLiftDistribution

  subroutine writeTecplot(sliceFile, writeSlices, liftFile, writeLift, &
       surfFile, writeSurf, famList, nFamList)
    !
    !       This is the master routine for writing tecplot data from adflow.
    !       This routine will write the slice, lift and surface files
    !       depending on the flags writeSlics, writeLift and writeSurface.
    !       The reason for the combined routine is that we can safely only
    !       perform the nodal averaging once which is required for all
    !       three output files.
    use constants
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use surfaceFamilies, only : BCFamGroups, BCFamExchange
    use surfaceUtils, only : getSurfaceSize
    use oversetData, only : zipperMeshes
    use outputMod, only : numberOfSurfSolVariables
    implicit none

    ! Input Params
    character*(*), intent(in) :: sliceFile, liftFile, surfFile
    logical, intent(in) :: writeSlices, writeLift, writeSurf
    integer(kind=intType), intent(in), dimension(nFamList) :: famList
    integer(kind=intType), intent(in) :: nFamList
    real(kind=realType), dimension(:, :, :), allocatable :: nodalValues
    ! Working
    integer(kind=intType) :: sps, nSolVar, sizeNOde, sizeCell
    integer(kind=intType), dimension(:), pointer :: wallLIST

    ! Determine the number of surface variables we have
    call numberOfSurfSolVariables(nSolVar)
    wallList =>  BCFamGroups(iBCGroupWalls)%famList
    call getSurfaceSize(sizeNode, sizeCell, wallList, size(wallList), .True.)

    ! Allocate and compute the wall-based surface data for hte slices
    ! and lift distributions.
    allocate(nodalValues(max(sizeNode,1), nSolVar+6+3, nTimeIntervalsSpectral))

    do sps=1, nTimeIntervalsSpectral
       call computeSurfaceOutputNodalData(BCFamExchange(iBCGroupWalls, sps), &
            zipperMeshes(iBCGroupWalls), .True., nodalValues(:, :, sps))
    end do

    if (writeSlices) then
       call writeSlicesFile(sliceFile, nodalValues)
    end if

    if (writeLift) then
       call writeLiftDistributionFile(liftFile, nodalValues)
    end if

    deallocate(nodalValues)

    if (writeSurf) then
       call writeTecplotSurfaceFile(surfFile, famList)
    end if

  end subroutine writeTecplot

  subroutine writeSlicesFile(fileName, nodalValues)
    !
    ! This subroutine is intended to be called from python.  This
    ! routine will write the user defined slics to an to the (ascii)
    ! tecplot file fileName. ASCII files are used for simplicity since
    ! very little information is actually written.
    use constants
    use communication
    use outputMod
    use inputTimeSpectral
    use inputPhysics
    use inputIteration
    use inputIO
    use surfaceFamilies
    use surfaceUtils
    use utils, only : EChk
    implicit none

    ! Input Params
    character*(*), intent(in) :: fileName
    real(kind=realType), intent(inout), dimension(:, :, :) :: nodalValues

    ! Working parameters
    integer(kind=intType) :: file, i, sps, nSolVar, ierr
    character(len=maxStringLen) :: fname
    character(len=7) :: intString
    character(len=maxCGNSNameLen), dimension(:), allocatable :: solNames
    integer(kind=intType), allocatable, dimension(:) :: famList
    type(slice) :: globalSlice
    integer(kind=intType) :: sizeNode, sizeCell
    integer(kind=intType), dimension(:), pointer :: wallList
    real(kind=realType), dimension(:, :), allocatable :: pts
    integer(kind=intType), dimension(:, :), allocatable :: conn
    integer(kind=intType), dimension(:), allocatable :: elemFam, cgnsBlockID

    ! Only write if we actually have lift distributions
    testwriteSlices: if(nParaSlices + nAbsSlices > 0) then

       if(myID == 0 .and. printIterations) then
          print "(a)", "#"
          print "(a)", "# Writing slices file(s) ..."
       endif

      !  do sps=1,nTimeIntervalsSpectral
       do sps=1,1

          ! If it is time spectral we need to agument the filename
          if (equationMode == timeSpectral) then
             write(intString,"(i7)") sps
             intString = adjustl(intString)
             fname = trim(fileName)//"Spectral"//trim(intString)
          else
             fname = fileName
          end if

          file = 11
          ! Open file on root proc:
          if (myid == 0) then
             open(unit=file, file=trim(fname))

             ! Write Header Information
             write (file,*) "Title = ""ADflow Slice Data"""
             write (file,"(a)", advance="no") "Variables = "
             write(file,"(a)",advance="no") " ""CoordinateX"" "
             write(file,"(a)",advance="no") " ""CoordinateY"" "
             write(file,"(a)",advance="no") " ""CoordinateZ"" "
             write(file,"(a)",advance="no") " ""XoC"" "
             write(file,"(a)",advance="no") " ""YoC"" "
             write(file,"(a)",advance="no") " ""ZoC"" "

             ! Number of additional variables on slices
             call numberOfSurfSolVariables(nSolVar)
             allocate(solNames(nSolVar))
             call surfSolNames(solNames)

             ! Write the rest of the variables
             do i=1,nSolVar
                write(file,"(a,a,a)",advance="no") """",trim(solNames(i)),""" "
             end do

             write(file,"(1x)")
             deallocate(solNames)
          end if
          call mpi_bcast(nSolVar, 1, adflow_integer, 0, adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)


          ! Slices are created on walls and walls only. Retrieve the
          ! points, connectivity and familyID of all the walls.
          wallList =>  BCFamGroups(iBCGroupWalls)%famList
          call getSurfaceSize(sizeNode, sizeCell, wallList, size(wallList), .True.)
          allocate(pts(3, sizeNode), conn(4, sizeCell), elemFam(sizeCell), cgnsBlockID(sizeCell))
          call getSurfaceConnectivity(conn, cgnsBlockID, sizeCell, wallList, size(wallList), .True.)
          call getSurfacePoints(pts, sizeNode, sps, wallList, size(wallList), .True.)
          call getSurfaceFamily(elemFam, sizeCell, wallList, size(wallList), .True.)

          ! Integration is performed in parallel
          do i=1, nParaSlices
             call integrateSlice(paraSlices(i, sps), globalSlice, &
                  nodalValues(:, :, sps), nSolVar, .True.)
             if (myid == 0) then
                call writeSlice(globalSlice, file, nSolVar)
             end if
             call destroySlice(globalSlice)
          end do

          do i=1, nAbsSlices
             ! 'Destroy' the slice...just dealloc the allocated data.
             ! before we do, save the family list
             allocate(famList(size(absSlices(i, sps)%famList)))
             famList = absSlices(i, sps)%famList
             call destroySlice(absSlices(i, sps))

             ! Make new one in the same location
             call createSlice(pts, conn, elemFam, absSlices(i, sps), &
                  absSlices(i, sps)%pt, absSlices(i, sps)%dir, &
                  absSlices(i, sps)%sliceName, famList)

             call integrateSlice(absSlices(i, sps), globalSlice, &
                  nodalValues(:, :, sps), nSolVar, .True.)
             if (myid == 0) then
                call writeSlice(globalSlice, file, nSolVar)
             end if
             call destroySlice(globalSlice)
             deallocate(famList)
          end do

          !Close file on root proc
          if (myid == 0) then
             close(file)
          end if
       end do

       if(myID == 0 .and. printIterations) then
          print "(a)", "# Slices file(s) written"
          print "(a)", "#"
       endif
    end if testwriteSlices
  end subroutine writeSlicesFile

  subroutine writeLiftDistributionFile(fileName, nodalValues)
    !
    !
    !       This subroutine is intended to be called from python.
    !       This routine will write the added lift distributions
    !       to the (ascii) tecplot file fileName. ASCII files are
    !       used for siplicity since very little informatin is actually
    !       written.
    use constants
    use communication
    use outputMod
    use inputPhysics
    use inputTimeSpectral
    use inputIteration
    use surfaceFamilies
    implicit none

    ! Input Params
    character*(*), intent(in) :: fileName
    real(kind=realType), dimension(:, :, :), allocatable :: nodalValues

    ! Working parameters
    integer(kind=intType) :: file, sps
    character(len=maxStringLen) :: fname
    character(len=7) :: intString

    ! Only write if we actually have lift distributions
    testwriteLiftDists: if(nLiftDists > 0) then

       if(myID == 0 .and. printIterations) then
          print "(a)", "#"
          print "(a)", "# Writing lift distribution file(s) ..."
       endif

       do sps=1,nTimeIntervalsSpectral

          ! If it is time spectral we need to agument the filename
          if (equationMode == timeSpectral) then
             write(intString,"(i7)") sps
             intString = adjustl(intString)
             fname = trim(fileName)//"Spectral"//trim(intString)
          else
             fname = fileName
          end if

          file = 11
          ! Open file on root proc:
          if (myid == 0) then
             open(unit=file, file=trim(fname))
          end if

          call writeLiftDistributions(sps, file, nodalValues(:, :, sps))

          ! Close file on root proc
          if (myid == 0) then
             close(file)
          end if
       end do

       if(myID == 0 .and. printIterations) then
          print "(a)", "# Lift distribution file(s) written"
          print "(a)", "#"
       endif

    end if testwriteLiftDists
  end subroutine writeLiftDistributionFile

  subroutine writeLiftDistributions(sps, fileID, nodalValues)
    !
    !       This subroutine writes the liftdistribution for the specified
    !       spectral instance. It is assumed that the required file handles
    !       are already open and can be written to
    use constants
    use communication
    use outputMod
    use su_cgns
    use cgnsNames
    use surfaceFamilies, only : BCFamGroups, BCFamExchange
    use surfaceUtils
    use utils, only : EChk
    use sorting, only : famInList
    implicit none

    ! Input parameters
    integer(kind=intType), intent(in) :: sps, fileID
    real(kind=realType), dimension(:, :), intent(in) :: nodalValues

    real(kind=realType), dimension(3) :: xmin, xmax, xmin_local, xmax_local
    real(kind=realType), parameter :: tol=1e-8
    type(liftDist), pointer :: d
    integer(kind=intType) :: i, j, ii, jj, iDist, ierr
    real(kind=realType), dimension(:,:), allocatable :: values
    character(len=maxCGNSNameLen), dimension(:), allocatable :: liftDistNames
    real(kind=realType) :: dmin, dmax, sumL, sumD, span, delta, xCur(3)
    type(slice) :: localSlice, globalSlice
    integer(kind=intType) :: sizeNode, sizeCell
    integer(kind=intType), dimension(:), pointer :: wallList
    real(kind=realType), dimension(:, :), allocatable :: pts
    integer(kind=intType), dimension(:, :), allocatable :: conn
    integer(kind=intType), dimension(:), allocatable :: elemFam, cgnsBlockID

    ! Slices are created on walls and walls only. Retrieve the
    ! points, connectivity and familyID of all the walls.
    wallList =>  BCFamGroups(iBCGroupWalls)%famList
    call getSurfaceSize(sizeNode, sizeCell, wallList, size(wallList), .True.)
    allocate(pts(3, sizeNode), conn(4, sizeCell), elemFam(sizeCell), cgnsBlockID(sizeCell))
    call getSurfaceConnectivity(conn, cgnsBlockID, sizeCell, wallList, size(wallList), .True.)
    call getSurfacePoints(pts, sizeNode, sps, wallList, size(wallList), .True.)
    call getSurfaceFamily(elemFam, sizeCell, wallList, size(wallList), .True.)

    do iDist=1,nLiftDists

       d => liftDists(iDist)
       xmin_local = huge(real(zero))
       xmax_local = -huge(real(zero))

       ! Get the bounding box for the entire geometry we have been slicing.
       elemLoop: do i=1, size(conn, 2)
          if (famInList(elemFam(i), d%FamList)) then

             ! Extract each of the 4 nodes on this quad:
             do jj=1,4
                xCur = pts(:, conn(jj, i))
                ! Check the max/min on each index
                do ii=1,3
                   xmin_local(ii) = min(xmin_local(ii) , xCur(ii))
                   xmax_local(ii) = max(xmax_local(ii) , xCur(ii))
                end do
             end do
          end if
       end do elemLoop

       ! Globalize all min/max values.
       call mpi_allreduce(xmin_local, xmin, 3, adflow_real, MPI_MIN, &
            adflow_comm_world, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call mpi_allreduce(xmax_local, xmax, 3, adflow_real, MPI_MAX, &
            adflow_comm_world, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       d%delta = (xMax(d%dir_ind) - xMin(d%dir_ind))/dble((d%nSegments - 1))
       allocate(d%slicePts(3, d%nSegments))

       ! Zero out all segments
       d%slicePts = zero

       ! These are the variable names for the lift distribution:
       allocate(liftDistNames(nLiftDistVar))
       ! Set the names here
       liftDistNames(1) = "eta"
       liftDistNames(2) = cgnsCoorX
       liftDistNames(3) = cgnsCoorY
       liftDistNames(4) = cgnsCoorZ
       liftDistNames(5) = "Lift"
       liftDistNames(6) = "Drag"
       liftDistNames(7) = "Normalized Lift"
       liftDistNames(8) = "Normalized Drag"
       liftDistNames(9) = "CL"
       liftDistNames(10) = "CD"
       liftDistNames(11) = "CLp"
       liftDistNames(12) = "CDp"
       liftDistNames(13) = "CLv"
       liftDistNames(14) = "CDv"
       liftDistNames(15) = "Ellptical"
       liftDistNames(16) = "thickness"
       liftDistNames(17) = "twist"
       liftDistNames(18) = "chord"

       ! Only write header info for first distribution only
       if (myid == 0) then
          if (iDist == 1) then
             write (fileID,*) "Title= ""ADflow Lift Distribution Data"""
             write (fileID,"(a)", advance="no") "Variables = "
             do i=1,nLIftDistVar
                write(fileID,"(a,a,a)",advance="no") """",trim(liftDistNames(i)),""" "
             end do
             write(fileID,"(1x)")
          end if

          write (fileID,"(a,a,a)") "Zone T= """,trim(d%distName),""""
          write (fileID,*) "I= ",d%nSegments
          write (fileID,*) "DATAPACKING=BLOCK"
15        format (E14.6)
       end if

       allocate(values(d%nSegments, nLiftDistVar))
       values = zero

       do i=1,d%nSegments
          if (i==1) then
             d%slicePts(d%dir_ind, i) = xMin(d%dir_ind) + tol
          else if (i == d%nSegments) then
             d%slicePts(d%dir_ind, i) = xMin(d%dir_ind) + (i-1)*d%delta - tol
          else
             d%slicePts(d%dir_ind, i) = xMin(d%dir_ind) + (i-1)*d%delta
          end if
       end do

       ! Scaled Eta values
       dmin = minVal(d%slicePts(d%dir_ind, :))
       dmax = maxval(d%slicePts(d%dir_ind, :))

       values(:, 1) = (d%slicePts(d%dir_ind, :) - dmin)/(dmax-dmin)
       ! Coordinate Varaibles
       if (d%dir_ind == 1) then! X slices
          values(:, 2) = d%slicePts(1, :)
       else if (d%dir_ind == 2) then ! Y slices
          values(:, 3) = d%slicePts(2, :)
       else if (d%dir_ind == 3) then ! Z slices
          values(:, 4) = d%slicePts(3, :)
       end if

       do i=1, d%nSegments
          ! Make new one in the same location
          call createSlice(pts, conn, elemFam, localSlice, d%slicePts(:, i), &
               d%dir, "does_not_matter", d%famList)
          call integrateSlice(localSlice, globalSlice, nodalValues, 0, .False.)

          ! Total lift and drag
          values(i, 5)  = globalSlice%pL + globalSlice%vL
          values(i, 6)  = globalSlice%pD + globalSlice%vD

          ! Total CL and CD
          values(i, 9)  = globalSlice%CLp + globalSlice%CLv
          values(i, 10) = globalSlice%CDp + globalSlice%CDv

          ! Pressure lift and drag coefficients
          values(i, 11) = globalSlice%CLp
          values(i, 12) = globalSlice%CDp

          ! Viscous lift and drag coefficients
          values(i, 13) = globalSlice%CLv
          values(i, 14) = globalSlice%CDv

          ! t/c, twist, chord
          values(i, 16) = globalSlice%thickness
          values(i, 17) = globalSlice%twist
          values(i, 18) = globalSlice%chord

          call destroySlice(localSlice)
          call destroySlice(globalSlice)

       end do

       ! Sum up the lift and drag values from the slices:
       sumL = zero
       sumD = zero
       do i=1,d%nSegments-1
          sumL = sumL + half*(values(i, 5) + values(i+1, 5))
          sumD = sumD + half*(values(i, 6) + values(i+1, 6))
       end do

       ! Now compute the normalized lift, drag and elliptical since
       ! we know the sum. Note delta is non-dimensional!
       delta = values(2, 1) - values(1, 1)

       ! This is the "nondimensional" span...it basically takes into account if you have
       ! a wing not at the sym plane
       span = maxval(values(:, 1)) - minval(values(:, 1))
       dmin = minval(values(:, 1))
       sumL = sumL * delta
       sumD = sumD * delta

       do i=1,d%nSegments

          ! Normalized Lift
          values(i, 7) = values(i, 5) / sumL
          values(i, 8) = values(i, 6) / sumD

          ! elliptical value
          values(i, 15) = four/pi/span*sqrt(one-(values(i, 1)-dmin)**2/span**2)
       end do

       ! Write all variables in block format
       if (myid == 0) then
          do j=1,nLiftDistVar
             do i=1,d%nSegments
                write(fileID,15) values(i, j)
             end do
          end do
       end if

       ! Deallocate slice list and point list
       deallocate(d%slicePts)

       ! Destroy temp variables
       deallocate(liftDistNames, values)
    end do

  end subroutine writeLiftDistributions

 subroutine writeTecplotSurfaceFile(fileName, famList)
    use constants
    use communication, only : myid, adflow_comm_world, nProc
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputPhysics, only : equationMode
    use inputIteration, only : printIterations
    use inputIO, only :precisionsurfgrid, precisionsurfsol
    use outputMod, only : surfSolNames, numberOfSurfSolVariables
    use surfaceFamilies, onlY : BCFamExchange, famNames, familyExchange
    use utils, only : EChk, setPointers, setBCPointers
    use BCPointers, only : xx
    use sorting, only : famInList
    use extraOutput, only : surfWriteBlank
    use oversetData, only : zipperMesh, zipperMeshes
    use surfaceUtils
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input Params
    character*(*), intent(in) :: fileName
    integer(kind=intType), intent(in), dimension(:) :: famList

    ! Working parameters
    integer(kind=intType) :: i, j, nn, mm, fileID, iVar, ii, ierr, iSize
    integer(Kind=intType) :: nSolVar, iBeg, iEnd, jBeg, jEnd, sps, sizeNode, sizeCell
    integer(kind=intType) :: iBCGroup, iFam, iProc, nCells, nNodes, nCellsToWrite, iZone, lastZoneSharing
    character(len=maxStringLen) :: fname
    character(len=7) :: intString
    integer(kind=intType), dimension(:), allocatable :: nodeSizes, nodeDisps
    integer(kind=intType), dimension(:), allocatable :: cellSizes, cellDisps
    character(len=maxCGNSNameLen), dimension(:), allocatable :: solNames
    real(kind=realType), dimension(:, :), allocatable :: nodalValues
    integer(kind=intType), dimension(:, :), allocatable :: conn, localConn
    real(kind=realType), dimension(:, :), allocatable :: vars
    integer(kind=intType), dimension(:), allocatable :: mask, elemFam, localElemFam, cgnsBlockID
    logical :: blankSave, BCGroupNeeded, dataWritten
    type(zipperMesh), pointer :: zipper
    type(familyExchange), pointer :: exch
    if(myID == 0 .and. printIterations) then
       print "(a)", "#"
       print "(a)", "# Writing tecplot surface file(s) ..."
    endif

    ! Number of surface variables. Note that we *explictly*
    ! remove the potential for writing the surface blanks as
    ! these are not necessary for the tecplot IO as we write the
    ! zipper mesh directly. We must save and restore the
    ! variable in case the CGNS otuput still wants to write it.
    blankSave = surfWriteBlank
    surfWriteBlank = .False.
    call numberOfSurfSolVariables(nSolVar)
    allocate(solNames(nSolVar))
    call surfSolNames(solNames)

    spectralLoop: do sps=1,nTimeIntervalsSpectral

       ! If it is time spectral we need to agument the filename
       if (equationMode == timeSpectral) then
          write(intString,"(i7)") sps
          intString = adjustl(intString)
          fname = trim(fileName)//"Spectral"//trim(intString)
       else
          fname = fileName
       end if

       fileID = 11
       ! Open file on root proc:

       if (myid == 0) then
          open(unit=fileID, file=trim(fname), form='UNFORMATTED', access='stream', status='replace')

          ! Tecplot magic number
          write(fileID) "#!TDV112"

          ! Integer value of 1 (4 bytes)
          call writeInteger(1)

          ! Integer for FileType: 0 = Full, 1= Grid, 2 = Solution
          call writeInteger(0)

          ! Write the title of the file
          call writeString("ADflow Surface Solution Data")

          ! Write the number of variable names
          call writeInteger(3 + nSolVar)

          ! Write the variable names
          call writeString("CoordinateX")
          call writeString("CoordinateY")
          call writeString("CoordinateZ")

          ! Write the rest of the variables
          do i=1, nSolVar
             call writeString(trim(solNames(i)))
          end do

          deallocate(solNames)

       end if

       ! First pass through to generate and write header information
       masterBCLoop1: do iBCGroup=1, nFamExchange

          ! Pointers for easier reading
          exch => BCFamExchange(iBCGroup, sps)
          zipper => zipperMeshes(iBCGroup)

          ! First thing we do is figure out if we actually need to do
          ! anything with this BCgroup at all. If none the requested
          ! families are in this BCExcahnge we don't have to do
          ! anything.

          BCGroupNeeded = .False.
          do i=1,size(famList)
             if (famInLIst(famList(i), exch%famList)) then
                BCGroupNeeded = .True.
             end if
          end do

          ! Keep going if we don't need this.
          if (.not. BCGroupNeeded) then
             cycle
          end if

          ! Get the sizes of this BCGroup
          call getSurfaceSize(sizeNode, sizeCell, exch%famList, size(exch%famList), .True.)
          call mpi_reduce(sizeNode, nNodes, 1, adflow_integer, MPI_SUM, 0, adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Now gather up the cell family info.
          allocate(cellDisps(0:nProc), cellSizes(nProc))

          call mpi_gather(sizeCell, 1, adflow_integer, &
               cellSizes, 1, adflow_integer, 0, adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          allocate(localElemFam(sizeCell))
          call getSurfaceFamily(localElemFam, sizeCell, exch%famList, size(exch%famList), .True.)

          if (myid == 0) then
             cellDisps(0) = 0
             do iProc=1, nProc
                cellDisps(iProc) = cellDisps(iProc-1) + cellSizes(iProc)
             end do
             nCells = sum(cellSizes)
             allocate(elemFam(nCells))
          end if

          call mpi_gatherv(localElemFam, &
               size(localElemFam), adflow_integer, elemFam, &
               cellSizes, cellDisps, adflow_integer, 0, adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Deallocate temp memory
          deallocate(localElemFam, cellSizes, cellDisps)

          rootProc: if (myid == 0 .and. nCells > 0) then
             do iFam=1, size(exch%famList)

                ! Check if we have to write this one:
                famInclude: if (famInList(exch%famList(iFam), famList)) then
                   nCellsToWrite = 0
                   do i=1, nCells
                      ! Check if this elem is to be included
                      if (elemFam(i) == exch%famList(iFam)) then
                         nCellsToWrite = nCellsToWrite + 1
                      end if
                   end do

                   if (nCellsToWrite > 0) then
                      call writeFloat(zoneMarker) ! Zone Marker
                      call writeString(trim(famNames(exch%famList(iFam)))) ! Zone Name
                      call writeInteger(-1) ! Parent Zone (-1 for None)
                      call writeInteger(-1) ! Strand ID (-2 for tecplot assignment)
                      call writeDouble(zero) ! Solution Time
                      call writeInteger(-1)! Zone Color (Not used anymore) (-1)
                      call writeInteger(3) ! Zone Type (3 for FEQuadrilateral)
                      call writeInteger(0)! Data Packing (0 for block)
                      call writeInteger(0)! Specify Var Location (0=don't specify, all at nodes)
                      call writeInteger(0) ! Are raw 1-to-1 face neighbours supplied (0 for false)
                      call writeInteger(nNodes) ! Number of nodes in FE Zone
                      call writeInteger(nCellsToWrite) ! Number of elements in FE Zone
                      call writeInteger(0) ! ICellDim, jCellDim, kCellDim (for future use, set to 0)
                      call writeInteger(0)
                      call writeInteger(0)
                      call writeInteger(0) ! Aux data specified (0 for no)
                   end if
                end if famInclude
             end do
          end if rootProc
          if (myid == 0) then
             deallocate(elemFam)
          end if
       end do masterBCLoop1

       if (myid == 0) then
          call writeFloat(dataSectionMarker) ! Eohmarker to mark difference between header and data section
       end if

       ! Now do everything again but for real.
       masterBCLoop: do iBCGroup=1,nFamExchange

          ! Pointers for easier reading
          exch => BCFamExchange(iBCGroup, sps)
          zipper => zipperMeshes(iBCGroup)

          ! First thing we do is figure out if we actually need to do
          ! anything with this BCgroup at all. If none the requested
          ! families are in this BCExcahnge we don't have to do
          ! anything.

          BCGroupNeeded = .False.
          do i=1,size(famList)
             if (famInList(famList(i), exch%famList)) then
                BCGroupNeeded = .True.
             end if
          end do

          ! Keep going if we don't need this.
          if (.not. BCGroupNeeded) then
             cycle
          end if

          ! Get the sizes of this BCGroup
          call getSurfaceSize(sizeNode, sizeCell, exch%famList, size(exch%famList), .True.)
          allocate(nodalValues(max(sizeNode,1), nSolVar+3+6))
          ! Compute the nodal data

          call computeSurfaceOutputNodalData(BCFamExchange(iBCGroup, sps), &
               zipperMeshes(iBCGroup), .False. , nodalValues(:, :))

          ! Gather up the number of nodes to be set to the root proc:
          allocate(nodeSizes(nProc), nodeDisps(0:nProc))
          nodeSizes = 0
          nodeDisps = 0

          call mpi_allgather(sizeNode, 1, adflow_integer, nodeSizes, 1, adflow_integer, &
               adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)
          nodeDisps(0) = 0
          do iProc=1, nProc
             nodeDisps(iProc) = nodeDisps(iProc-1) + nodeSizes(iProc)
          end do

          iSize = 3 + 6 + nSolVar
          if (myid == 0) then
             nNodes = sum(nodeSizes)
          else
             nNodes = 1
          end if

          ! Only root proc actually has any space allocated
          allocate(vars(nNodes, iSIze))

          ! Gather values to the root proc.
          do i=1, iSize
             call mpi_gatherv(nodalValues(:, i), sizeNode, &
                  adflow_real, vars(:, i), nodeSizes, nodeDisps, adflow_real, 0, adflow_comm_world, ierr)
             call EChk(ierr,__FILE__,__LINE__)
          end do
          deallocate(nodalValues)

          ! Now gather up the connectivity
          allocate(cellDisps(0:nProc), cellSizes(nProc))

          call mpi_gather(sizeCell, 1, adflow_integer, &
               cellSizes, 1, adflow_integer, 0, adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)


          if (allocated(cgnsBlockID)) then
             deallocate(cgnsBlockID)
          end if

          allocate(localConn(4, sizeCell), localElemFam(sizeCell), cgnsBlockID(sizeCell))
          call getSurfaceConnectivity(localConn, cgnsBlockID, sizeCell, exch%famList, size(exch%famList), .True.)
          call getSurfaceFamily(localElemFam, sizeCell, exch%famList, size(exch%famList), .True.)

          if (myid == 0) then
             cellDisps(0) = 0
             do iProc=1, nProc
                cellDisps(iProc) = cellDisps(iProc-1) + cellSizes(iProc)
             end do
             nCells = sum(cellSizes)
             allocate(conn(4, nCells))
             allocate(elemFam(nCells))
          end if

          ! We offset the conn array by nodeDisps(iProc) which
          ! automagically adjusts the connectivity to account for the
          ! number of nodes from different processors

          call mpi_gatherv(localConn+nodeDisps(myid), &
               4*size(localConn, 2), adflow_integer, conn, &
               cellSizes*4, cellDisps*4, adflow_integer, 0, adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          call mpi_gatherv(localElemFam, &
               size(localElemFam), adflow_integer, elemFam, &
               cellSizes, cellDisps, adflow_integer, 0, adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Local values are finished
          deallocate(localConn, localElemFam)
          iZone = 0
          rootProc2: if (myid == 0 .and. nCells > 0) then

             ! Need zero based for binary output.
             conn = conn - 1

             allocate(mask(nCells))
             dataWritten=.False.
             do iFam=1, size(exch%famList)

                ! Check if we have to write this one:
                famInclude2: if (famInList(exch%famList(iFam), famList)) then

                   ! Create a temporary mask
                   mask = 0
                   nCellsToWrite = 0
                   do i=1, nCells
                      ! Check if this elem is to be included
                      if (elemFam(i) == exch%famList(iFam)) then
                         mask(i) = 1
                         nCellsToWrite = nCellsToWrite + 1
                      end if
                   end do

                   actualWrite2 : if (nCellsToWrite > 0) then

                      ! Write Zone Data
                      call writeFloat(zoneMarker)

                      ! Data Type for each variable (1 for float, 2 for double)
                      do i=1, 3
                         if (precisionSurfGrid == precisionSingle) then
                            call writeInteger(1)
                         else if (precisionSurfGrid == precisionDouble) then
                            call writeInteger(2)
                         end if
                      end do

                      do i=1, nSolVar
                         if (precisionSurfSol == precisionSingle) then
                            call writeInteger(1)
                         else if (precisionSurfSol == precisionDouble) then
                            call writeInteger(2)
                         end if
                      end do

                      call writeInteger(0) ! Has passive variables (0 for no)

                      if (.not. dataWritten) then

                         ! Save the current 0-based zone so we know which to share with.
                         lastZoneSharing = iZone

                         call writeInteger(0) ! Has variable sharing (0 for no)
                         call writeInteger(-1)  ! Zone based zone number to share connectivity (-1 for no sharing)

                         ! Min/Max Value for coordinates
                         do j=1,3
                            call writeDouble(minval(vars(:, j)))
                            call writeDouble(maxval(vars(:, j)))
                         end do

                         ! Min/Max Value for solution variables
                         do j=1,nSolVar
                            call writeDouble(minval(vars(:, j+9)))
                            call writeDouble(maxval(vars(:, j+9)))
                         end do

                         ! Dump the coordinates
                         do j=1,3
                            if (precisionSurfGrid == precisionSingle) then
                               call writeFloats(vars(1:nNodes, j))
                            else if (precisionSurfSol == precisionDouble) then
                               call writeDoubles(vars(1:nNodes, j))
                            end if
                         end do

                         ! Dump the solution variables
                         do j=1,nSolVar
                            if (precisionSurfSol == precisionSingle) then
                               call writeFloats(vars(1:nNodes, j+9))
                            else if (precisionSurfSol == precisionDouble) then
                               call writeDoubles(vars(1:nNodes, j+9))
                            end if
                         end do

                      else
                         ! This will be sharing data from another zone.
                         call writeInteger(1) ! Has variable sharing (0 for no)

                         ! Write out the zone number for sharing the data.
                         do j=1,3+nSolVar
                            call writeInteger(lastZoneSharing)
                         end do

                         call writeInteger(-1)  ! Zone based zone number to share connectivity (-1 for no sharing
                      end if

                      ! Dump the connectivity
                      j = 0
                      do i=1, nCells
                         ! Check if this elem is to be included
                         if (mask(i) == 1) then
                            call writeIntegers(conn(:, i))
                            j = j + 1
                         end if
                      end do

                      iZone = iZone + 1
                   end if actualWrite2
                end if famInclude2
             end do
             deallocate(mask)
          end if rootProc2
          if (myid == 0) then
             deallocate(conn, elemFam)
          end if
          deallocate(cellSizes, cellDisps, nodeSizes, nodeDisps, vars)
       end do masterBCLoop

       if (myid == 0) then
          close(fileID)
       end if

    end do spectralLoop

    ! Restore the modified option
    surfWriteBlank = blankSave
    if(myID == 0 .and. printIterations) then
       print "(a)", "# Tecplot surface file(s) written"
       print "(a)", "#"
    endif

    contains

      subroutine writeFloat(adflowRealVal)
        use iso_fortran_env, only : real32
        implicit none
        real(kind=realType) :: adflowRealVal
        real(kind=real32) :: float
        float = adflowRealval
        write(fileID) float
      end subroutine writeFloat

      subroutine writeDouble(adflowRealVal)
        use iso_fortran_env, only : real64
        implicit none
        real(kind=realType) :: adflowRealVal
        real(kind=real64) :: dble
        dble = adFlowRealVal
        write(fileID) dble
      end subroutine writeDouble

      subroutine writeFloats(adflowRealVals)
        use iso_fortran_env, only : real32
        implicit none
        real(kind=realType) :: adflowRealVals(:)
        real(kind=real32) :: floats(size(adflowRealVals))
        integer :: i
        floats = adflowRealvals
        write(fileID) floats

      end subroutine writeFloats

      subroutine writeDoubles(adflowRealVals)
        use iso_fortran_env, only : real64
        implicit none
        real(kind=realType) :: adflowRealVals(:)
        real(kind=real64) :: dbles(size(adflowrealvals))
        integer :: i
        dbles = adflowrealvals
        write(fileID) dbles

      end subroutine writeDoubles

      subroutine writeInteger(adflowIntegerVal)
        use iso_fortran_env, only : int32
        implicit none
        integer(kind=intType) :: adflowIntegerVal
        integer(kind=int32) :: int

        int = adflowIntegerVal
        write(fileID) int
      end subroutine writeInteger

      subroutine writeIntegers(adflowIntegerVals)
        use iso_fortran_env, only : int32
        implicit none
        integer(kind=intType) :: adflowIntegerVals(:), i
        integer(kind=int32) :: ints(size(adflowintegervals))
        ints = adflowintegervals
        write(fileID) ints

      end subroutine writeIntegers

      subroutine writeString(str)

        implicit none

        character*(*)::  str
        integer(kind=intType) :: i

        do i=1,len(str)
           write(fileID) iachar(str(i:i))
        end do
        write(fileID) 0

      end subroutine writeString

    end subroutine writeTecplotSurfaceFile

  subroutine initializeLiftDistributionData
    use constants
    use communication
    use blockPointers
    use inputPhysics
    use outputMod
    use inputTimeSpectral
    use surfacefamilies
    implicit none

    ! Working Variables
    integer(kind=intType) :: nPts, nCells, sps

    if (liftDistInitialized) then
       return
    else

       ! Data for the marching squares method: Which edges are cut by
       ! the contour. We haven't dealt with the case of an ambiguous
       ! contour which is case 6 and 11. Since most of the slices we are
       ! doing are planes this won't matter.
       msCon1( 1 , : ) = (/ 0 , 0 , 0 , 0 , 0 /)
       msCon1( 2 , : ) = (/ 1 , 4 , 0 , 0 , 0 /)
       msCon1( 3 , : ) = (/ 1 , 2 , 0 , 0 , 0 /)
       msCon1( 4 , : ) = (/ 4 , 2 , 0 , 0 , 0 /)
       msCon1( 5 , : ) = (/ 2 , 3 , 0 , 0 , 0 /)
       msCon1( 6 , : ) = (/ 2 , 3 , 0 , 0 , 0 /) ! Should be 2, 3, 1, 4
       msCon1( 7 , : ) = (/ 1 , 3 , 0 , 0 , 0 /)
       msCon1( 8 , : ) = (/ 4 , 3 , 0 , 0 , 0 /)
       msCon1( 9 , : ) = (/ 4 , 3 , 0 , 0 , 0 /)
       msCon1( 10, : ) = (/ 1 , 3 , 0 , 0 , 0 /)
       msCon1( 11 , : ) = (/ 2 , 3 , 0 , 0 , 0 /) ! Should be 2, 3, 1, 4
       msCon1( 12 , : ) = (/ 2 , 3 , 0 , 0 , 0 /)
       msCon1( 13 , : ) = (/ 4 , 2 , 0 , 0 , 0 /)
       msCon1( 14 , : ) = (/ 1 , 2 , 0 , 0 , 0 /)
       msCon1( 15 , : ) = (/ 1 , 4 , 0 , 0 , 0 /)
       msCon1( 16 , : ) = (/ 0 , 0 , 0 , 0 , 0 /)

       msCon2(1, :) = (/1, 2/)
       msCon2(2, :) = (/2, 3/)
       msCon2(3, :) = (/3, 4/)
       msCon2(4, :) = (/4, 1/)
       liftDistInitialized = .True.
    end if

  end subroutine initializeLiftDistributionData

  subroutine computeSurfaceOutputNodalData(exch, zipper, includeTractions, nodalValues)
    !
    !       This purpose of this subroutine is to compute all nodal values
    !
    use constants
    use communication
    use inputPhysics
    use blockPointers
    use surfaceFamilies, only : BCFamGroups, familyExchange
    use outputMod, only : storeSurfSolInBuffer, numberOfSurfSolVariables, &
         surfSolNames
    use surfaceUtils
    use utils, only : setPointers, EChk
    use sorting, only : famInList
    use oversetData, only : zipperMesh
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none
    ! Input Param
    type(familyExchange) :: exch
    logical :: includeTractions
    real(kind=realType), dimension(:, :), intent(inout) :: nodalValues
    type(zipperMesh) :: zipper
    ! Working params
    integer(kind=intType) :: i, j, ii, jj, kk, nn, mm, iSol, ierr, nPts, nCells
    integer(kind=intType) :: nFields, nSolVar, iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
    integer(kind=intType) :: sizeNode, sizeCell, iDim
    integer(kind=intType), dimension(3,2) :: cellRangeCGNS
    character(len=maxCGNSNameLen), dimension(:), allocatable :: solNames
    real(kind=realType), dimension(:), allocatable :: buffer
    real(kind=realType), dimension(:), pointer :: weightPtr, localPtr
    real(kind=realType), dimension(:, :), allocatable :: tmp
    logical :: viscousSubFace

    nodalValues = zero

    call numberOfSurfSolVariables(nSolVar)
    allocate(solNames(nSolVar))
    call surfSolNames(solNames)

    ! The tractions have a sepcial routine so call that first before we
    ! mess with the family information.

    if (includeTractions) then

       ! And put the traction where it needs to go if necessary.  Note that
       ! this computation will only work corectly when the wallExchange
       ! which include all wallgroups is used.

       call computeNodalTractions(exch%sps)
       ii = 0
       do nn=1, nDom
          call setPointers(nn, 1_intType, exch%sps)
          do mm=1, nBocos
             iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
             jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd

             if (famInList(BCdata(mm)%famID, exch%famList)) then
                do j=jBeg, jEnd
                   do i=iBeg, iEnd
                      ii = ii + 1
                      nodalValues(ii, 4) = BCData(mm)%Tp(i, j, 1)
                      nodalValues(ii, 5) = BCData(mm)%Tp(i, j, 2)
                      nodalValues(ii, 6) = BCData(mm)%Tp(i, j, 3)

                      nodalValues(ii, 7) = BCData(mm)%Tv(i, j, 1)
                      nodalValues(ii, 8) = BCData(mm)%Tv(i, j, 2)
                      nodalValues(ii, 9) = BCData(mm)%Tv(i, j, 3)
                   end do
                end do
             end if
          end do
       end do

       ! Not quite dont yet with the nodal tractions; we need to send
       ! the nodal tractions that the duplices the *zipper* mesh needs
       ! to the root proc. We have a special scatter for this.

       if (zipper%allocated) then

          ! Loop over the 6 tractions
          do iDim=1,6

             ! Copy the values into localPtr
             call VecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
             call EChk(ierr,__FILE__,__LINE__)

             do i=1,exch%nNodes
                localPtr(i) = nodalValues(i, iDim+3)
             end do

             call VecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
             call EChk(ierr,__FILE__,__LINE__)

             ! Now use the *zipper* scatter
             call VecScatterBegin(zipper%scatter, exch%nodeValLocal,&
                  zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
             call EChk(ierr,__FILE__,__LINE__)

             call VecScatterEnd(zipper%scatter, exch%nodeValLocal,&
                  zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
             call EChk(ierr,__FILE__,__LINE__)

             ! Copy the zipper values out on the root proc
             if (myid == 0) then
                call VecGetArrayF90(zipper%localVal, localPtr, ierr)
                call EChk(ierr,__FILE__,__LINE__)

                do i=1,size(localPtr)
                   nodalValues(exch%nNodes+i, iDim+3) = localPtr(i)
                end do

                call VecRestoreArrayF90(zipper%localVal, localPtr, ierr)
                call EChk(ierr,__FILE__,__LINE__)
             end if
          end do
       end if
    end if

    ! Get the current set of surface points for the family we just set.
    call getSurfaceSize(sizeNode, sizeCell, exch%famList, size(exch%famlist), .True.)
    allocate(tmp(3, sizeNode))
    call getSurfacePoints(tmp, sizeNode, exch%sps, exch%famList, size(exch%famList), .True.)

    do i=1, sizeNode
       nodalValues(i, 1:3) = tmp(1:3, i)
    end do
    deallocate(tmp)
    ! For the remainder of the variables, use arithematic averaging.
    call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    localPtr = zero
    ! ii is the index into the pointer array
    ii = 0
    do nn=1, nDom
       call setPointers(nn, 1_intType, exch%sps)
       do mm=1, nBocos
          iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
          jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
          ni = iEnd - iBeg + 1
          nj = jEnd - jBeg + 1
          if (famInList(BCdata(mm)%famID, exch%FamList)) then
             do j=0,nj-2
                do i=0,ni-2

                   ! Scatter 1 to each node.
                   ind(1) = ii + (j  )*ni + i + 1
                   ind(2) = ii + (j  )*ni + i + 2
                   ind(3) = ii + (j+1)*ni + i + 2
                   ind(4) = ii + (j+1)*ni + i + 1
                   do jj=1,4
                      localPtr(ind(jj)) = localPtr(ind(jj)) + one
                   end do
                end do
             end do
             ii = ii + ni*nj
          end if
       end do
    end do

    call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Globalize the area
    call vecSet(exch%sumGlobal, zero, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
         exch%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
         exch%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Now compute the inverse of the weighting so that we can multiply
    ! instead of dividing.

    call vecGetArrayF90(exch%sumGlobal, localPtr, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    localPtr = one/localPtr

    call vecRestoreArrayF90(exch%sumGlobal, localPtr, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    varLoop: do iSol=1, nSolVar

       ! Extract the poitner to the local array
       call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! ii is the continous running element pointer
       ii = 0
       localPtr = zero
       ! Do each variable separately.
       domainLoop: do nn=1, nDom
          call setPointers(nn, 1, exch%sps)

          bocoLoop: do mm=1, nBocos
             if (famInList(BCdata(mm)%famID, exch%famList)) then

                ! storeSurfSolInBuffer needs to know if the subface is
                ! viscous or not.
                viscousSubFace = .False.
                if (BCType(mm) == NSWallAdiabatic .or. &
                     BCType(mm) == NSWallIsoThermal) then
                   viscousSubFace = .true.
                end if

                ! Determine the cell range *in the original cgns
                ! ordering*. The reason for this is the storeSurfSolInBufer
                ! is normally used for writing CGNS and thus it is that
                ! ording that is important. However, that isn't too hard to
                ! deal with.

                jBeg = BCData(mm)%jnBeg + 1
                jEnd = BCData(mm)%jnEnd
                iBeg = BCData(mm)%inBeg + 1
                iEnd = BCData(mm)%inEnd

                ! Dummy value of 1 for the face values not set.
                cellRangeCGNS = 1
                select case(BCFaceID(mm))
                case (iMin, iMax)
                   ! I range meaningless
                   cellRangeCGNS(2, 1) = iBeg + jBegOr - 1
                   cellRangeCGNS(2, 2) = iEnd + jBegOr - 1

                   cellRangeCGNS(3, 1) = jBeg + kBegOr - 1
                   cellRangeCGNS(3, 2) = jEnd + kBegOr - 1

                case (jMin, jMax)
                   ! J range meaningless
                   cellRangeCGNS(1, 1) = iBeg + iBegOr - 1
                   cellRangeCGNS(1, 2) = iEnd + iBegOr - 1

                   cellRangeCGNS(3, 1) = jBeg + kBegOr - 1
                   cellRangeCGNS(3, 2) = jEnd + kBegOr - 1

                case (kMin, kMax)
                   ! J range meaningless
                   cellRangeCGNS(1, 1) = iBeg + iBegOr - 1
                   cellRangeCGNS(1, 2) = iEnd + iBegOr - 1

                   cellRangeCGNS(2, 1) = jBeg + jBegOr - 1
                   cellRangeCGNS(2, 2) = jEnd + jBegOr - 1
                end select

                ! Allocate enough space for the 1D buffer
                allocate(buffer((iEnd-iBeg+1)*(jEnd-jBeg+1)))
                kk = 0
                call storeSurfsolInBuffer(exch%sps, buffer, kk, nn, BCfaceID(mm), &
                     cellRangeCGNS, solNames(iSol), viscousSubFace, .False.)

                ! Now since the storeSurfSol just put things in a flat
                ! array and are face based, here we take the 1D face data
                ! and scatter to the nodes.

                iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
                jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
                ni = iEnd - iBeg + 1
                nj = jEnd - jBeg + 1
                jj = 0
                do j=0,nj-2
                   do i=0,ni-2
                      jj = jj + 1
                      ! Scatter value to each node
                      ind(1) = ii + (j  )*ni + i + 1
                      ind(2) = ii + (j  )*ni + i + 2
                      ind(3) = ii + (j+1)*ni + i + 2
                      ind(4) = ii + (j+1)*ni + i + 1
                      do kk=1,4
                         localPtr(ind(kk)) = localPtr(ind(kk)) + buffer(jj)
                      end do
                   end do
                end do

                ii = ii + ni*nj
                deallocate(buffer)
             end if
          end do bocoLoop
       end do domainLoop

       ! Return our pointer
       call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Globalize the function value
       call vecSet(exch%nodeValGlobal, zero, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
            exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
            exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Now divide by the weighting sum. We can do this with a
       ! vecpointwisemult since we already divided by the weight.
       call vecPointwiseMult(exch%nodeValGlobal, exch%nodeValGlobal, &
            exch%sumGlobal, ierr)
       call EChk(ierr,__FILE__,__LINE__)


       ! Push back to the local values
       call VecScatterBegin(exch%scatter, exch%nodeValGlobal, &
            exch%nodeValLocal, INSERT_VALUES, SCATTER_REVERSE, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecScatterEnd(exch%scatter, exch%nodeValGlobal, &
            exch%nodeValLocal, INSERT_VALUES, SCATTER_REVERSE, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Copy the values into nodalValues
       call VecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       do i=1,size(localPtr)
          nodalValues(i, iSol+9) = localPtr(i)
       end do

       call VecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       if (zipper%allocated) then

          ! Now use the *zipper* scatter
          call VecScatterBegin(zipper%scatter, exch%nodeValLocal,&
               zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          call VecScatterEnd(zipper%scatter, exch%nodeValLocal,&
               zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Copy the zipper values out on the root proc
          if (myid == 0) then
             call VecGetArrayF90(zipper%localVal, localPtr, ierr)
             call EChk(ierr,__FILE__,__LINE__)

             do i=1,size(localPtr)
                nodalValues(exch%nNodes+i, 9+iSol) = localPtr(i)
             end do

             call VecRestoreArrayF90(zipper%localVal, localPtr, ierr)
             call EChk(ierr,__FILE__,__LINE__)
          end if
       end if
    end do varLoop
    deallocate(solNames)

  end subroutine computeSurfaceOutputNodalData

  subroutine createSlice(pts, conn, elemFam, slc, pt, dir, sliceName, famList)
    !
    !       This subroutine creates a slice on a plane defined by pt and
    !       and dir. It only uses the families specified in the famList.
    !       sps define which specral instance to use.
    !
    use constants
    use utils, only : reallocatereal2, reallocateinteger2, pointReduce
    use sorting, only : famInList
    implicit none

    ! Input param
    real(kind=realType), dimension(:, :), intent(in) :: pts
    integer(kind=intType), dimension(:,:), intent(in) :: conn
    integer(kind=intType), dimension(:), intent(in) :: elemFam
    type(slice), intent(inout) :: slc
    real(kind=realType), dimension(3), intent(in) :: pt, dir
    character*(*), intent(in) :: sliceName
    integer(kind=intType), dimension(:), intent(in) :: famList

    ! Working param
    integer(kind=intType) :: i, nMax, nUnique, oldInd, newInd
    integer(kind=intType) :: patchIndices(4), indexSquare, jj, kk, icon, iCoor, num1, num2
    real(kind=realType) :: f(4), d, ovrdnom, tol
    logical :: logic1, foundFam
    real(kind=realType), dimension(:, :), pointer :: tmpWeight, dummy, tmpNodes
    integer(kind=intType), dimension(:, :), pointer :: tmpInd
    integer(kind=intType), dimension(:), allocatable :: link
    real(kind=realType), dimension(:), allocatable :: fc

    ! Allocate the family list this slice is to use:
    slc%sliceName = sliceName
    ! Set the info for the slice:
    slc%pt = pt
    slc%dir = dir
    slc%nNodes = 0
    allocate(slc%famList(size(famList)))
    slc%famList = famList
    ! First step is to compute the 'function' value that will be used
    ! for the contour.

    ! Equation of plane: ax + by + cz + d = 0
    d = -pt(1)*dir(1) - pt(2)*dir(2) - pt(3)*dir(3)
    ovrdnom = one/sqrt(dir(1)**2 + dir(2)**2 + dir(3)**2)

    ! Compute the distance function on all possible surfaces on this
    ! processor.
    allocate(fc(size(pts,2)))
    do i=1, size(pts, 2)
       ! Now compute the signed distance
       fc(i) = (dir(1)*pts(1, i) + dir(2)*pts(2, i) + dir(3)*pts(3, i) + d)*ovrdnom
    end do

    ! Estimate size of slice by the sqrt of the number of nodes in the
    ! mesh. Exact size doesn't matter as we realloc if necessary.
    nMax = int(sqrt(dble(size(pts, 2))))
    allocate(tmpWeight(2, nMax), tmpInd(2, nMax), tmpNOdes(3, nMax))

    iCoor = 0
    oldInd = 1

    ! Loop over all elements
    elemLoop: do i=1, size(conn, 2)

       famInclude: if (famInList(elemFam(i), famList)) then

          ! Extract the indices and function values at each corner
          do jj=1,4
             patchIndices(jj) = conn(jj, i)
             f(jj) = fc(patchIndices(jj))
          end do

          ! Based on the values at each corner, determine which
          ! type contour we have
          indexSquare = 1

          if (f(1) .lt. zero) indexsquare = indexsquare + 1
          if (f(2) .lt. zero) indexsquare = indexsquare + 2
          if (f(3) .lt. zero) indexsquare = indexsquare + 4
          if (f(4) .lt. zero) indexsquare = indexsquare + 8

          logic1 = .true.

          kk = 1
          do while (logic1)
             ! This is the edge
             icon = mscon1(indexSquare, kk)

             if (icon == 0) then
                logic1=.false.
             else

                ! num1, num2 are node indices
                num1 = mscon2(icon,1)
                num2 = mscon2(icon,2)

                iCoor = iCoor + 1
                if (iCoor > nMax) then
                   ! Need to reallocate the arrays. Make it double the size
                   call reallocateReal2(tmpWeight, 2, 2*nMax, 2, nMax, .true.)
                   call reallocateReal2(tmpNodes, 3, 2*nMax, 3, nMax, .true.)
                   call reallocateInteger2(tmpInd, 2, 2*nMax, 2, nMax, .true.)
                   nMax = nMax * 2
                end if

                ! Weight factors
                tmpWeight(2, iCoor) = (zero - f(num1))/(f(num2) - f(num1))
                tmpWeight(1, iCoor) = one - tmpWeight(2, icoor)

                ! Store the weight factors
                tmpInd(:, iCoor) = (/patchIndices(num1), patchIndices(num2)/)

                ! Store the physical nodes so we know how to reduce
                tmpNodes(:, iCoor) = &
                     tmpWeight(1, iCoor)*pts(:, tmpInd(1, iCoor)) + &
                     tmpWeight(2, iCoor)*pts(:, tmpInd(2, iCoor))
                kk = kk + 1
             end if
          end do
       end if famInclude
    end do  ElemLoop

    ! To save disk space, we can compact out the doubly defined nodes that
    ! were created during the slicing process. Then we can allocate the
    ! final weight array and index array to be the exact correct
    ! length

    allocate(dummy(3, iCoor), link(iCoor))
    tol=1e-12
    call pointReduce(tmpNodes, iCoor, tol, dummy, link, nUnique)
    allocate(slc%w(2, nUnique), slc%ind(2, nUnique), slc%conn(2, iCoor/2))
    slc%nNodes = nUnique

    ! Modify the data accordingly
    do i=1, iCoor
       slc%w(:, link(i)) = tmpWeight(:, i)
       slc%ind(:, link(i)) = tmpInd(:, i)
    end do

    ! The connectivity is actually link reshaped since the original conn
    ! would have been (1,2), (3,4), (5,6) etc.
    do i=1, iCoor/2
       slc%conn(1, i) = link(2*i-1)
       slc%conn(2, i) = link(2*i)
    end do

    deallocate(tmpNodes, tmpWeight, tmpInd, dummy, link, fc)
  end subroutine createSlice

  subroutine destroySlice(slc)
    !
    !       This subroutine destroys a slice created by the createSlice
    !       routine
    !
    use constants
    implicit none

    ! Input param
    type(slice), intent(inout) :: slc

    ! Deallocate weights and indices if they are already allocated
    if (allocated(slc%w)) then
       deallocate(slc%w)
    end if

    if (allocated(slc%ind)) then
       deallocate(slc%ind)
    end if

    if (allocated(slc%conn)) then
       deallocate(slc%conn)
    end if

    if (allocated(slc%famList)) then
       deallocate(slc%famList)
    end if

    if (allocated(slc%vars)) then
       deallocate(slc%vars)
    end if

  end subroutine destroySlice

  subroutine integrateSlice(lSlc, gSlc, nodalValues, nFields, doConnectivity)
    !
    !       This subroutine integrates the forces on slice slc and computes
    !       the integrated quantities for lift, drag, cl and cd with
    !       contributions from both the pressure and visoucs forces.
    !       It optionally interpolates solution variables as well.
    !
    use constants
    use inputPhysics
    use flowVarRefState
    use communication
    use utils, only : EChk
    implicit none

    ! Input Variables
    type(slice) :: lSlc, gSlc
    integer(kind=intType), intent(in) :: nFields
    logical, intent(in) :: doConnectivity
    real(kind=realType), dimension(:, :), intent(in) :: nodalValues
    ! Working variables
    integer(kind=intType) :: i, j, i1, i2
    real(kind=realType), dimension(3) :: x1, x2, pT1, pT2, vT1, vT2, pF, vF
    real(kind=realType) :: len, dmax, dmin, dist, fact, M(3,3), tmp(4)
    real(kind=realType) :: r(3), r_new(3), hyp, te(3), le(3), theta, w1, w2
    integer(kind=intType) :: bestPair(2), dir_ind, iProc, ierr, iSize
    real(kind=realtype), dimension(:,:), allocatable :: tempCoords
    real(kind=realtype), dimension(:,:), allocatable :: localVals
    integer(kind=intType), dimension(:), allocatable :: sliceNodeSizes, sliceCellSizes
    integer(kind=intType), dimension(:), allocatable :: nodeDisps, cellDisps

    ! Copy the info related to slice
    gSlc%sliceName = trim(lSlc%sliceName)
    gSlc%pt = lSlc%pt
    gSlc%dir = lSlc%dir

    pF = zero
    vF = zero
    iSize = 3 + 6 + nFields

    allocate(localVals(iSize,lSlc%nNodes))

    ! Compute all the local variables we need.

    ! Interpolate the required values
    do i=1, lSlc%nNodes
       i1 = lSlc%ind(1, i)
       i2 = lSlc%ind(2, i)
       w1 = lSlc%w(1, i)
       w2 = lSlc%w(2, i)
       localVals(1:iSize, i) = w1*nodalValues(i1, 1:iSize) + w2*nodalValues(i2, 1:iSize)
    end do

    do i=1, size(lSlc%conn, 2)

       ! extract nodes:
       i1 = lslc%conn(1, i)
       i2 = lslc%conn(2, i)

       x1 = localVals(1:3, i1)
       x2 = localVals(1:3, i2)

       ! extract pressure tractions
       pT1 = localVals(4:6, i1)
       pT2 = localVals(4:6, i2)

       ! extract viscous tractions
       vT1 = localVals(7:9, i1)
       vT2 = localVals(7:9, i2)

       ! Length of this segment
       len = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2 + (x1(3)-x2(3))**2)

       ! Integrate the pressure and viscous forces separately
       pF = pF + half*(pT1 + pT2)*len
       vF = vF + half*(vT1 + vT2)*len

    end do

    ! That is as far as we can go in parallel. We now have to gather up
    ! pL, pD, vL vD as well as the nodes to the root proc.

    ! Gather up the number of nodes to be set to the root proc:
    allocate(sliceNodeSizes(nProc), nodeDisps(0:nProc))
    sliceNodeSizes = 0
    nodeDisps = 0
    call mpi_allgather(lSlc%nNodes,1, adflow_integer, sliceNodeSizes, 1, adflow_integer, &
         adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)
    nodeDisps(0) = 0
    do iProc=1, nProc
       nodeDisps(iProc) = nodeDisps(iProc-1) + sliceNodeSizes(iProc)*iSize
    end do

    if (myid == 0) then
       gSlc%nNodes = sum(sliceNodeSizes)
       allocate(gSlc%vars(iSize, gSlc%nNodes))

    end if

    call mpi_gatherv(localVals, iSize*lSlc%nNodes, adflow_real, gSlc%vars, sliceNodeSizes*iSize, &
         nodeDisps, adflow_real, 0, adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! We may also need to gather the connectivity if the slice will have
    ! to be written to a file.
    if (doConnectivity) then
       i = size(lslc%conn, 2)
       allocate(cellDisps(0:nProc), sliceCellSizes(nProc))

       call mpi_gather(i, 1, adflow_integer, sliceCellSizes, 1, adflow_integer, &
            0, adflow_comm_world, ierr)
       call EChk(ierr,__FILE__,__LINE__)


       if (myid == 0) then
          cellDisps(0) = 0
          do iProc=1, nProc
             cellDisps(iProc) = cellDisps(iProc-1) + sliceCellSizes(iProc)*2
          end do
          allocate(gSlc%conn(2, sum(sliceCellSizes)))
       end if

       ! We offset the conn array by nodeDisps(iProc) which
       ! automagically adjust the connectivity to account for the
       ! number of nodes from different processors

       call mpi_gatherv(lSlc%conn+nodeDisps(myid)/iSize, 2*size(lSlc%conn, 2), adflow_integer, gSlc%conn, &
            sliceCellSizes*2, cellDisps, adflow_integer, 0, adflow_comm_world, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Not quite finished yet since we will have gathered nodes from
       ! multiple procs we have to adjust the connectivity

       deallocate(sliceCellSizes, cellDisps)
    end if

    deallocate(sliceNodeSizes, nodeDisps)
    ! Set the local values we can in the slice
    lSlc%pL = liftDirection(1)*pF(1) + liftDirection(2)*pF(2) + liftDirection(3)*pF(3)
    lSlc%pD = dragDirection(1)*pF(1) + dragDirection(2)*pF(2) + dragDirection(3)*pF(3)
    lSlc%vL = liftDirection(1)*vF(1) + liftDirection(2)*vF(2) + liftDirection(3)*vF(3)
    lSlc%vD = dragDirection(1)*vF(1) + dragDirection(2)*vF(2) + dragDirection(3)*vF(3)

    ! Reduce the lift/drag values
    call mpi_reduce((/lSlc%pL, lSlc%pD, lSlc%vL, lSlc%vD/), tmp, 4, adflow_real, MPI_SUM, &
         0, adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    if (myid == 0) then
       gSlc%pL = tmp(1)
       gSlc%pD = tmp(2)
       gSlc%vL = tmp(3)
       gSlc%vD = tmp(4)

       ! Compute the chord as the max length between any two nodes...this
       ! is n^2, so should be changed in the future

       dmax = zero
       bestPair = (/1, 1/)
       do i=1,size(gSlc%vars, 2)
          ! extract node:
          x1 = gSlc%vars(1:3, i)

          do j=i+1,size(gSlc%vars, 2)
             ! extract node:
             x2 = gSlc%vars(1:3, j)

             dist = sqrt((x1(1)-x2(1))**2 +(x1(2)-x2(2))**2 + (x1(3)-x2(3))**2)

             if (dist > dmax) then
                dmax = dist
                bestPair = (/i, j/)
             end if
          end  do
       end do

       ! Set chord, protected from zero
       gSlc%chord = max(dmax, 1e-12)

       ! Compute factor to get coefficient
       fact = two/(gammaInf*pInf*MachCoef*MachCoef*pRef)

       ! Take dmax as chord and compute coefficients
       gSlc%CLp = gSlc%pL / gSlc%chord * fact
       gSlc%CDp = gSlc%pD / gSlc%chord * fact
       gSlc%CLv = gSlc%vL / gSlc%chord * fact
       gSlc%CDv = gSlc%vD / gSlc%chord * fact

       ! Default values
       gSlc%twist = zero
       gSlc%thickness = zero

       if (gSlc%nNodes == 0) then
          return
       end if

       ! Lastly we need the twist and the twist and the thickness
       i1 = bestPair(1)
       i2 = bestPair(2)
       x1 = gSlc%vars(1:3, i1)
       x2 = gSlc%vars(1:3, i2)

       if (x1(1) > x2(1)) then
          te = x1
          le = x2
       else
          te = x2
          le = x1
       end if

       ! Save the leading and trailing edges so we can do scaled output
       ! later
       gSlc%le = le
       gSlc%te = te

       ! Finally we need to get the thickness. For this, compute temporary
       ! section nodes and rotate them by the twist values we just computed
       ! and take the max and min

       ! Back out what is the main index of the slice, x, y or z based on
       ! the direction. Not the best approach, but that's ok
       dir_ind = maxloc(abs(gSlc%dir),1)

       ! Length of hyptoneuse is the same
       hyp = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2 + (x1(3)-x2(3))**2)

       if (dir_ind == 1) then
          ! Xslice...we don't how what to do here..could be y or z. Don't
          ! do anything.
          gSlc%twist = zero
       else if (dir_ind == 2) then
          ! Yslice
          theta = asin((le(3)-te(3))/hyp)
          gSlc%twist = theta*180.0/pi
       else
          ! Zslice
          theta = asin((le(2)-te(2))/hyp)
          gSlc%twist = theta*180.0/pi
       end if

       if (dir_ind == 1) then
          M(1,1) = one; M(1,2) = zero; M(1, 3) = zero;
          M(2,1) = zero; M(2,2) = one; M(2, 3) = zero;
          M(3,1) = zero; M(3,2) = zero; M(3,3) = one;
       else if(dir_ind == 2) then
          ! Y-rotation matrix
          M(1,1) = cos(-theta); M(1,2) = zero; M(1, 3) = sin(-theta);
          M(2,1) = zero; M(2,2) = one; M(2, 3) = zero;
          M(3,1) = -sin(-theta); M(3,2) = zero; M(3,3) = cos(-theta);
       else
          ! Z rotation Matrix
          M(1,1) = cos(theta); M(1,2) = -sin(theta); M(1,3) = zero;
          M(2,1) = sin(theta); M(2,2) = cos(theta); M(2,3) = zero;
          M(3,1) = zero; M(3,2) = zero; M(3,3) = one;
       end if

       allocate(tempCoords(3, size(gSlc%vars, 2)))
       do i=1, size(gSlc%vars, 2)
          ! extract node:
          r = gSlc%vars(1:3, i) - te
          r_new = matmul(M, r)
          tempCoords(:, i) = r_new + te
       end do

       ! Now get the max and the min and divide by the chord for t/c
       if (dir_ind == 1) then
          gSlc%thickness = 0 ! Again, don't know what to do here
       else if(dir_ind == 2) then
          dmax = maxval(tempCoords(3, :))
          dmin = minval(tempCoords(3, :))
          gSlc%thickness = (dmax-dmin)/hyp
       else if(dir_ind == 3) then
          dmax = maxval(tempCoords(2, :))
          dmin = minval(tempCoords(2, :))
          gSlc%thickness = (dmax-dmin)/hyp
       end if
       deallocate(tempCoords)
    end if
  end subroutine integrateSlice

  subroutine writeSlice(slc, fileID, nFields)
    ! Write the data in slice 'slc' to openfile ID fileID
    use constants
    use inputIO
    implicit none

    ! Input Parameters
    type(slice), intent(in) :: slc
    integer(kind=intType), intent(in) :: fileID, nFields

    ! Working Variables
    integer(kind=intType) :: i, j
    real(kind=realType) :: tmp, tx, ty, tz
    write (fileID,"(a,a,a)") "Zone T= """,trim(slc%sliceName),""""

    ! IF we have nodes actually write:
    if (slc%nNodes > 0) then
       write (fileID,*) "Nodes = ", slc%nNodes, " Elements= ",  size(slc%conn, 2), " ZONETYPE=FELINESEG"
       write (fileID,*) "DATAPACKING=POINT"
13     format (E14.6)

       do i=1,slc%nNodes
          ! Write the coordinates
          do j=1,3
             write(fileID,13, advance='no') slc%vars(j, i)
          end do

          ! Write the scaled coordiantes with the LE at (0,0,0)
          do j=1,3
             tmp = slc%vars(j, i)
             write(fileID,13, advance='no') (tmp - slc%le(j))/slc%chord
          end do

          ! Write field data. Starts at 9 (after 3 coordindates and the 6 tractions)
          do j=1,nFields
             write(fileID,13, advance='no') slc%vars(9+j, i)
          end do

          write(fileID,"(1x)")
       end do

15     format(I5, I5)
       do i=1, size(slc%conn, 2)
          write(fileID, 15) slc%conn(1, i), slc%conn(2, i)
       end do
    else ! Write dummy data so the number of zones are the same

       write (fileID,*) "Nodes = ", 2, " Elements= ", 1, " ZONETYPE=FELINESEG"
       write (fileID,*) "DATAPACKING=POINT"
       do i=1,2
          do j=1,6
             write(fileID,13, advance='no') zero
          end do

          do j=1,nFields
             write(fileID,13, advance='no') zero
          end do

          write(fileID,"(1x)")
       end do
       write(fileID, 15) 1, 2
    end if
  end subroutine writeSlice

end module tecplotIO
