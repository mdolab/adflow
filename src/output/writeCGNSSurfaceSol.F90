!
!      ******************************************************************
!      *                                                                *
!      * File:          writeCgnsSurfaceSol.F90                         *
!      * Author:        Edwin van der Weide , Gaetan K. W. Kenway       *
!      * Starting date: 05-15-2003                                      *
!      * Last modified: 09-29-2015                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine writeCGNSSurfaceSol
  !
  !      ******************************************************************
  !      *                                                                *
  !      * writeCGNSSurfaceSol and its subroutines write the surface      *
  !      * solution file(s). The unknowns are stored in the center of the *
  !      * surface quadrilaterals when nodalOuput is False and stored at  *
  !      * nodes when nodalOutput is True.                                *
  !      *                                                                *
  !      ******************************************************************
  !
#ifdef USE_NO_CGNS
  call terminate("writeCGNSSurfaceSol", &
       "Routine should not be called if no cgns support &
       &is selected.")
#else
  use cgnsGrid
  use communication
  use su_cgns
  use outputMod
  use inputIteration
  use block
  use blockPointers
  use cgnsNames
  use extraOutput
  use inputIO
  implicit none
  !
  !      Local parameter, the cell dimension.
  !
  integer, parameter :: celldim = 2
  !
  !      Local variables.
  !
  integer :: cgnsInd, ierr

  integer(kind=intType) :: nn, mm, ll
  integer(kind=intType) :: nSolVar, nZonesWritten
  character(len=maxStringLen) :: errorMessage
  integer(kind=intType) :: iSurf, nisoSurfVar
  character(len=maxCGNSNameLen), dimension(:), allocatable :: &
       solNames, isoSurfSolNames
  character(len=maxCGNSNameLen) :: contourName
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Determine the number and names of the solution files.

  call surfSolFileNamesWrite

  ! Return immediately if no surface solution files must
  ! be written.

  if(nSurfSolToWrite == 0) return

  ! If we have nodal output computed the weights of all boundary
  ! conditions once to save doing the Newton Rapson search for each variable
  if (nodalOutput) then 
     call computeSurfaceNodalWeights
  end if

  ! Write a message that the solution file(s) are being written.
  ! Of course only processor 0 does this.

  if(myID == 0 .and. printIterations) then
     print "(a)", "#"
     print "(a,a)", "# Writing surface solution file(s): ", trim(surfSolFileNames(1))
  endif

  ! Allocate the memory for the fileIDs and the bases.

  allocate(fileIDs(nSurfSolToWrite), cgnsBases(nSurfSolToWrite), &
       stat=ierr)
  if(ierr /= 0)                           &
       call terminate("writeCGNSSurfaceSol", &
       "Memory allocation failure for fileIDs &
       &and cgnsBases")

  ! Open the cgns file(s) and write the header. This is only done
  ! by processor 0.

  testRootProc: if(myID == 0) then

     ! Loop over the number of surface solution files to write.

     solLoop: do nn=1,nSurfSolToWrite

        ! Open the cgns file for writing and check if it went okay.
        ! Store the file index for later purposes.
        call cg_open_f(surfSolFileNames(nn), mode_write, cgnsInd, &
             ierr)
        if(ierr /= CG_OK) then
           write(errorMessage,101) trim(surfSolFileNames(nn))
101        format("File",1X,A,1X,"could not be opened by cgns for &
                &writing")

           call terminate("writeCGNSSurfaceSol", errorMessage)
        endif

        fileIDs(nn) = cgnsInd

        ! Create the base.

        call cg_base_write_f(cgnsInd, "BaseSurfaceSol", celldim, &
             cgnsPhysdim, cgnsBases(nn), ierr)
        if(ierr /= CG_OK)                      &
             call terminate("writeCGNSSurfaceSol", &
             "Something wrong when calling &
             &cg_base_write_f")

        ! Write the header in the cgns file.
        call writeCGNSHeader(cgnsInd, cgnsBases(nn))

     enddo solLoop

  endif testRootProc

  ! Determine the number of variables to be written to the surface
  ! solution file as well as the cgns names.

  call numberOfSurfSolVariables(nSolVar)

  allocate(solNames(nSolVar), stat=ierr)
  if(ierr /= 0)                           &
       call terminate("writeCGNSSurfaceSol", &
       "Memory allocation failure for solNames")

  call surfSolNames(solNames)

  ! Loop over the number of cgns blocks and its boundary subfaces
  ! and write the cell centered surface solution of the subface.

  nZonesWritten = 0
  zoneLoop: do nn=1,cgnsNDom

     ! Determine the number of blocks on this processor that belong
     ! to this cgns block.

     mm = nblocksCGNSblock(nn) - nblocksCGNSblock(nn-1)

     ! Loop over the number of boundary subfaces of the original
     ! cgns block and write the cell centered surface solution to
     ! the cgns surface file.

     do ll=1,cgnsDoms(nn)%nBocos

        ! Only write the solution to file if this is a true subface.
        if( cgnsDoms(nn)%bocoInfo(ll)%actualFace )                 &
             call writeSurfsolCGNSZone(nn, mm, ll, nSolVar, solNames, &
             nZonesWritten, .false.)
     enddo

     ! Loop over the number of internal block boundaries of the
     ! original grid and write the periodic boundaries.

     do ll=1,cgnsDoms(nn)%n1to1

        ! Only periodic boundaries are written; check for this.

        if( cgnsDoms(nn)%conn1to1(ll)%periodic )                   &
             call writeSurfsolCGNSZone(nn, mm, ll, nSolVar, solNames, &
             nZonesWritten, .true.)
     enddo
  enddo zoneLoop


  ! Check if isosurface will be written. These will be written to
  ! a new base

  testIsoSurafce: if (nIsoSurface > 0)  then

     allocate (cgnsIsoSurfBases(nSurfSolToWrite), stat=ierr)
     testRootProc2: if (myID == 0) then

        ! Loop over the number of surface solution files

        solLoop2: do nn=1,nSurfSolToWrite

           ! Create the new base

           cgnsInd = fileIDs(nn)
           call cg_base_write_f(cgnsInd, "IsoSurfaces", celldim, &
                cgnsPhysDim, cgnsIsoSurfBases(nn), ierr)
           if (ierr /= CG_OK) &
                call terminate("WriteCGNSSurfaceSol", &
                "Something wrong when calling cg_base_write_f for &
                isoSurface")
        end do solLoop2
     end if testRootProc2

     ! Determine the number of variables to be written to the
     ! isosurface itself well as the cgns names. 

     call numberOfIsoSurfVariables(nIsoSurfVar)

     if (nIsoSurfVar > 0) then
        allocate(isoSurfSolNames(nIsoSurfVar), stat=ierr)
        if(ierr /= 0)                           &
             call terminate("writeCGNSSurfaceSol", &
             "Memory allocation failure for isoNames")
        call isoSurfNames(isoSurfSolNames)
     end if

     solLoop3: do ll=1,nSurfSolToWrite ! Numer of spectral instances!
        ! Allocate fn and fc for each domain:
        do nn=1,nDom
           call setPointers(nn, 1, ll)
           allocate(flowDoms(nn, 1, ll)%fn(il, jl, kl))
           allocate(flowDoms(nn, 1, ll)%fc(1:ie, 1:je, 1:ke))
        end do

        ! Finally loop over the required isoSurfaces
        do iSurf=1,nIsoSurface
           call computeIsoVariable(isoSurfaceNames(iSurf), ll, isoValues(iSurf))

11         format(A,A,A,F7.4)
           write(contourName, 11), "Contour ", trim(isoSurfaceNames(iSurf)), "=", isoValues(iSurf)
           call writeIsoSurface(contourName, ll, nIsoSurfVar, isoSurfSolNames)
        end do

        ! deAllocate fn and fc for each domain:
        do nn=1,nDom
           deallocate(flowDoms(nn, 1, ll)%fn, flowDoms(nn, 1, ll)%fc)
        end do
     end do solLoop3

     ! Free memory for bases
     deallocate(cgnsIsoSurfBases, stat=ierr)
     if (nIsoSurfVar > 0) then
        deallocate(isoSurfSolNames)
     end if
  end if testIsoSurafce

  ! Close the cgns file(s). Only processor 0 does this.

  if(myID == 0) then
     do nn=1,nSurfSolToWrite
        call cg_close_f(fileIDs(nn), ierr)
        if(ierr /= CG_OK)                      &
             call terminate("writeCGNSSurfaceSol", &
             "Something wrong when calling cg_close_f")
     enddo
  end if

  ! Deallocate the memory of solNames, fileIDs and cgnsBases.

  deallocate(solNames, fileIDs, cgnsBases, stat=ierr)
  if(ierr /= 0)                           &
       call terminate("writeCGNSSurfaceSol", &
       "Deallocation error for solNames, fileIDs &
       &and cgnsBases")

  ! If we have nodal output computed, deallocate the nodal weights
  if (nodalOutput) then 
     call deallocateSurfaceNodalWeights
  end if

  ! Wait until all processors (especially processor 0) reach
  ! this point.

  call mpi_barrier(SUmb_comm_world, ierr)

  ! Write a message that the solution file(s) have been written.
  ! Of course only processor 0 does this.

  if(myID == 0 .and. printIterations) then
     print "(a)", "# Surface solution file(s) written"
     print "(a)", "#"
  endif

#endif

contains

  subroutine computeSurfaceNodalWeights

    use bcTypes
    use blockPointers
    use inputTimeSpectral
    use overset
    integer(kind=intType) :: nn, sps, mm, i, j, k, iSize, jSize
    real(kind=realType), dimension(:,:,:),   pointer :: xx
    real(kind=realType), dimension(3) :: n1, n2, n3, n4, x1, x21, x41, x3142, xf
    real(kind=realType), dimension(3) :: an, bn, vt, b, a, vf
    real(kind=realType) :: u, v, du, dv, uv, uold, vvold, val, vn, invLen
    real(kind=realType), dimension(3) :: norm

    integer(kind=intType), parameter :: iterMax   = 15
    real(kind=realType),   parameter :: adteps    = 1.e-25_realType
    real(kind=realType),   parameter :: thresConv = 1.e-10_realType
    logical, dimension(6) :: boundaryOnFace

    ! We do each spectal instance of each owned blocks
    donaminLoop: do nn=1,nDom
       spectralLoop: do sps=1,nTimeIntervalsSpectral
          call setPointers(nn, 1_intType, sps)
          boundaryOnFace = .False.
          bocoLoop: do mm=1,nBocos
             ! Flag this face as having a boundary
             boundaryOnFace(BCFaceID(mm)) = .True.
          end  do bocoLoop

          ! Nullify all the weight pointers. This also serves as
          ! poitner initialization. 
          do mm=1,6
             nullify(flowDoms(nn, 1_intType, sps)%nodalWeights(mm)%weight)
          end do

          ! Generically loop from 1 to 6 (6 total faces) and only
          ! compute the nodal weights if that face had a BC on it

          faceIDLoop: do mm=1, 6

             boundaryNeeded: if (boundaryOnFace(mm)) then 

                ! Extract the pointer to the nodes on the surface
                select case (mm)
                case (iMin)
                   xx => x(1, :, :, :)
                   iSize = jl
                   jSize = kl
                case (iMax)
                   xx => x(il, :, :, :)
                   iSize = jl
                   jSize = kl
                case (jMin)
                   xx => x(:, 1, :, :)
                   iSize = il
                   jSize = kl
                case (jMax)
                   xx => x(:, jl, :, :)
                   iSize = il
                   jSize = kl
                case (kMin)
                   xx => x(:, :, 1, :)
                   iSize = il
                   jSize = jl
                case (kMax)
                   xx => x(:, :, kl, :)
                   iSize = il
                   jSize = jl
                end select

                ! Allocate space for the nodal weighting 
                allocate(flowDoms(nn, 1_intType, sps)%nodalWeights(mm)%weight(4, iSize, jSize))

                ! Loop over the nodal range
                do j=1, jSize
                   do i=1, iSize

                      ! Extract out the dual mesh surrounding
                      ! this node. Everything is shifted by 1 due to pointer offset

                      n1 = fourth*(xx(i  , j  , :) + xx(i+1, j, :) + xx(i+1, j+1, :) + xx(i  , j+1, :))
                      n2 = fourth*(xx(i+1, j  , :) + xx(i+2, j, :) + xx(i+2, j+1, :) + xx(i+1, j+1, :))
                      n3 = fourth*(xx(i  , j+1, :) + xx(i+1, j+1, :) + xx(i+1, j+2, :) + xx(i  , j+2, :))
                      n4 = fourth*(xx(i+1, j+1, :) + xx(i+2, j+1, :) + xx(i+2, j+2, :) + xx(i+1, j+2, :))

                      ! We now need to find the weights for point xx(i, j, :) in
                      ! the quad defined n1, through n4 (cyclic). This will tell
                      ! us the 4 weights each of the 4 cells. 

                      x1(1) = n1(1)
                      x1(2) = n1(2)
                      x1(3) = n1(3)

                      x21(1) = n2(1) - x1(1)
                      x21(2) = n2(2) - x1(2)
                      x21(3) = n2(3) - x1(3)

                      x41(1) = n3(1) - x1(1)
                      x41(2) = n3(2) - x1(2)
                      x41(3) = n3(3) - x1(3)

                      x3142(1) = n4(1) - x1(1) - x21(1) - x41(1)
                      x3142(2) = n4(2) - x1(2) - x21(2) - x41(2)
                      x3142(3) = n4(3) - x1(3) - x21(3) - x41(3)

                      ! Initialize u and v to 0.5 and determine the
                      ! corresponding coordinates on the face, which is the
                      ! centroid.

                      u  = half
                      v  = half
                      uv = u*v

                      xf(1) = x1(1) + u*x21(1) + v*x41(1) + uv*x3142(1)
                      xf(2) = x1(2) + u*x21(2) + v*x41(2) + uv*x3142(2)
                      xf(3) = x1(3) + u*x21(3) + v*x41(3) + uv*x3142(3)

                      ! Newton loop to determine the point on the surface,
                      ! which minimizes the distance to the given coordinate.

                      NewtonQuads: do ll=1,iterMax

                         ! Store the current values of u and v for a stop
                         ! criterion later on.

                         uold = u
                         vvold = v

                         ! Determine the vector vf from xf to given coordinate.

                         vf(1) = xx(i+1, j+1, 1) - xf(1)
                         vf(2) = xx(i+1, j+1, 2) - xf(2)
                         vf(3) = xx(i+1, j+1, 3) - xf(3)

                         ! Determine the tangent vectors in u- and v-direction.
                         ! Store these in a and b respectively.

                         a(1) = x21(1) + v*x3142(1)
                         a(2) = x21(2) + v*x3142(2)
                         a(3) = x21(3) + v*x3142(3)

                         b(1) = x41(1) + u*x3142(1)
                         b(2) = x41(2) + u*x3142(2)
                         b(3) = x41(3) + u*x3142(3)

                         ! Determine the normal vector of the face by taking the
                         ! cross product of a and b. Afterwards this vector will
                         ! be scaled to a unit vector.

                         norm(1) = a(2)*b(3) - a(3)*b(2)
                         norm(2) = a(3)*b(1) - a(1)*b(3)
                         norm(3) = a(1)*b(2) - a(2)*b(1)

                         invLen = one/max(adtEps, sqrt(norm(1)*norm(1) &
                              +                        norm(2)*norm(2) &
                              +                        norm(3)*norm(3)))

                         norm(1) = norm(1)*invLen
                         norm(2) = norm(2)*invLen
                         norm(3) = norm(3)*invLen

                         ! Determine the projection of the vector vf onto
                         ! the face.

                         vn = vf(1)*norm(1) + vf(2)*norm(2) + vf(3)*norm(3)
                         vt(1) = vf(1) - vn*norm(1)
                         vt(2) = vf(2) - vn*norm(2)
                         vt(3) = vf(3) - vn*norm(3)

                         ! The vector vt points from the current point on the
                         ! face to the new point. However this new point lies on
                         ! the plane determined by the vectors a and b, but not
                         ! necessarily on the face itself. The new point on the
                         ! face is obtained by projecting the point in the a-b
                         ! plane onto the face. this can be done by determining
                         ! the coefficients du and dv, such that vt = du*a + dv*b.
                         ! To solve du and dv the vectors normal to a and b
                         ! inside the plane ab are needed.

                         an(1) = a(2)*norm(3) - a(3)*norm(2)
                         an(2) = a(3)*norm(1) - a(1)*norm(3)
                         an(3) = a(1)*norm(2) - a(2)*norm(1)

                         bn(1) = b(2)*norm(3) - b(3)*norm(2)
                         bn(2) = b(3)*norm(1) - b(1)*norm(3)
                         bn(3) = b(1)*norm(2) - b(2)*norm(1)

                         ! Solve du and dv. the clipping of vn should not be
                         ! active, as this would mean that the vectors a and b
                         ! are parallel. This corresponds to a quad degenerated
                         ! to a line, which should not occur in the surface mesh.

                         vn = a(1)*bn(1) + a(2)*bn(2) + a(3)*bn(3)
                         vn = sign(max(eps,abs(vn)),vn)
                         du = (vt(1)*bn(1) + vt(2)*bn(2) + vt(3)*bn(3))/vn

                         vn = b(1)*an(1) + b(2)*an(2) + b(3)*an(3)
                         vn = sign(max(eps,abs(vn)),vn)
                         dv = (vt(1)*an(1) + vt(2)*an(2) + vt(3)*an(3))/vn

                         ! Determine the new parameter values uu and vv. These
                         ! are limited to 0 <= (uu,vv) <= 1.

                         u = u + du; u = min(one,max(zero,u))
                         v = v + dv; v = min(one,max(zero,v))

                         ! Determine the final values of the corrections.

                         du = abs(u-uold)
                         dv = abs(v-vvold)

                         ! Determine the new coordinates of the point xf.

                         uv  = u*v
                         xf(1) = x1(1) + u*x21(1) + v*x41(1) + uv*x3142(1)
                         xf(2) = x1(2) + u*x21(2) + v*x41(2) + uv*x3142(2)
                         xf(3) = x1(3) + u*x21(3) + v*x41(3) + uv*x3142(3)

                         ! Exit the loop if the update of the parametric
                         ! weights is below the threshold

                         val = sqrt(du*du + dv*dv)
                         if(val <= thresConv) exit NewtonQuads

                      enddo NewtonQuads
             
                      ! Finally store the weights
                      flowDoms(nn, 1_intType, sps)%nodalWeights(mm)%weight(1, i, j) = (one - u)*(one -v)
                      flowDoms(nn, 1_intType, sps)%nodalWeights(mm)%weight(2, i, j) = (      u)*(one -v)
                      flowDoms(nn, 1_intType, sps)%nodalWeights(mm)%weight(3, i, j) = (one - u)*(     v)
                      flowDoms(nn, 1_intType, sps)%nodalWeights(mm)%weight(4, i, j) = (      u)*(     v)
                   end do
                end do
             end if boundaryNeeded
          end do faceIDLoop
       end do spectralLoop
    end do donaminLoop
  end subroutine computeSurfaceNodalWeights

  subroutine deallocateSurfaceNodalWeights

    use inputTimeSpectral
    integer(kind=intType) :: nn, sps, mm

    donaminLoop: do nn=1,nDom
       spectralLoop: do sps=1,nTimeIntervalsSpectral
          do mm=1,6
             if (associated(flowDoms(nn, 1_intType, sps)%nodalWeights(mm)%weight)) then
                deallocate(flowDoms(nn, 1_intType, sps)%nodalWeights(mm)%weight)
                nullify(flowDoms(nn, 1_intType, sps)%nodalWeights(mm)%weight)
             end if
          end do
       end do spectralLoop
    end do donaminLoop
  end subroutine deallocateSurfaceNodalWeights

end subroutine writeCGNSSurfaceSol
