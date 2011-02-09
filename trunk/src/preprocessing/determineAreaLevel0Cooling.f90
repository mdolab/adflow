!
!      ******************************************************************
!      *                                                                *
!      * File:          determineAreaLevel0Cooling.f90                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-30-2005                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineAreaLevel0Cooling(level)
!
!      ******************************************************************
!      *                                                                *
!      * determineAreaLevel0Cooling computes the axial area of the      *
!      * injection planes of the level 0 cooling model for the given    *
!      * multigrid. If only a slice is simulated the area is scaled up  *
!      * to the full wheel, such that no discrepencies occur when       *
!      * sections are present with different periodic angles.           *
!      *                                                                *
!      * This model has been developed by Pratt and Whitney and should  *
!      * not be given to third parties. This implementation assumes     *
!      * the x-direction is the axial direction.                        *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use block
       use communication
       use coolingModelLevel0
       use section
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer :: size, ierr

       integer(kind=intType) :: nn, mm, i, j
       integer(kind=intType) :: nSlices
       integer(kind=intType) :: blockID, indNorm, indexDir
       integer(kind=intType) :: jcBeg, jcEnd, icBeg, icEnd

       real(kind=realType), dimension(nPlanesLevel0CoolingModel) :: localArea
       real(kind=realType), dimension(nPlanesLevel0CoolingModel) :: area

       real(kind=realType), dimension(:,:,:), pointer :: ss
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if no cooling planes are present.

       if(nPlanesLevel0CoolingModel == 0) return

       ! Loop over the number of cooling planes to determine the local
       ! contributions to the area.

       planeLoop: do nn=1,nPlanesLevel0CoolingModel

         ! Initialize the local area and the number of slices to zero.

         localArea(nn) = zero
         nSlices       = 0

         ! Loop over the local subfaces of the cooling plane.

         subfaceLoop: do mm=1,level0Cooling(nn,level)%nSubfaces

           ! Store the needed subface parameters a bit easier.

           blockID  = level0Cooling(nn,level)%blockID(mm)
           indNorm  = level0Cooling(nn,level)%indNorm(mm)
           indexDir = level0Cooling(nn,level)%indexDir(mm)
           jcBeg    = level0Cooling(nn,level)%jcBeg(mm)
           jcEnd    = level0Cooling(nn,level)%jcEnd(mm)
           icBeg    = level0Cooling(nn,level)%icBeg(mm)
           icEnd    = level0Cooling(nn,level)%icEnd(mm)

           ! Set the pointer for the normals, depending on the
           ! index direction.

           select case (indexDir)
             case (iMin, iMax)
               ss => flowDoms(blockID,level,1)%sI(indNorm,:,:,:)
             case (jMin, jMax)
               ss => flowDoms(blockID,level,1)%sJ(:,indNorm,:,:)
             case (kMin, kMax)
               ss => flowDoms(blockID,level,1)%sK(:,:,indNorm,:)
           end select

           ! Loop over faces of the current subface and update the
           ! local axial area.

           do j=jcBeg,jcEnd
             do i=icBeg,icEnd
               localArea(nn) = localArea(nn) + ss(i,j,1)
             enddo
           enddo

           ! Determine the number of slices for this section.
           ! Throw an error if a different periodicity is encountered.

           i = flowDoms(blockID,level,1)%sectionID
           j = sections(i)%nSlices

           if(nSlices == 0) then
             nSlices = j
           else if(nSlices /= j) then
             call terminate("determineAreaLevel0Cooling", &
                            "Different periodicity encountered for &
                            &cooling plane.")
           endif

         enddo subfaceLoop

         ! Multiply the local axial area by the number of slices to
         ! obtain the value for the entire wheel.

         localArea(nn) = localArea(nn)*nSlices

       enddo planeLoop

       ! Do an allreduce to obtain the total area.

       size = nPlanesLevel0CoolingModel
       call mpi_allreduce(localArea, area, size, sumb_real, &
                          mpi_sum, SUmb_comm_world, ierr)

       ! Copy the area's into the appropriate places of level0Cooling.

       do nn=1,nPlanesLevel0CoolingModel
         level0Cooling(nn,level)%area = area(nn)
       enddo

       end subroutine determineAreaLevel0Cooling
