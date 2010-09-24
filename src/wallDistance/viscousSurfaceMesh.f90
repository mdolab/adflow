!
!      ******************************************************************
!      *                                                                *
!      * File:          viscousSurfaceMesh.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-12-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine viscousSurfaceMesh(level, sps)
!
!      ******************************************************************
!      *                                                                *
!      * viscousSurfaceMesh determines and stores the entire viscous    *
!      * surface possibly extended by periodic parts.                   *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use block
       use communication
       use section
       use viscSurface
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
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
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

       call mpi_allreduce(nQuadVisc, nquadViscGlob, 1, adt_integer, &
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
