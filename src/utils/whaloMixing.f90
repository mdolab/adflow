!
!      ******************************************************************
!      *                                                                *
!      * File:          whaloMixing.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-26-2005                                      *
!      * Last modified: 09-16-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine whaloMixing(level, startIn, endIn, commPressure,   &
                              commVarGamma, commLamVis, commEddyVis, &
                              nLayers)
!
!      ******************************************************************
!      *                                                                *
!      * whaloMixing determines the halo values for the cells adjacent  *
!      * to sliding interfaces for the given communication pattern.     *
!      * The mixing plane approximation is just a azimuthal averaging   *
!      * of the variables just that a steady flow can be computed for   *
!      * the inherently unsteady problem of rotor/stator interaction.   *
!      * It is possible to exchange a range of variables and not the    *
!      * entire set, e.g. only the flow variables or only the turbulent *
!      * variables. This is controlled by the arguments start, end,     *
!      * commPressure and commViscous. The exchange takes place for     *
!      * the given grid level.                                          *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use commMixing
       use constants
       use inputPhysics
       use inputTimeSpectral
       use interfaceGroups
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, startIn, endIn
       integer(kind=intType), intent(in) :: nLayers

       logical, intent(in) :: commPressure, commVarGamma
       logical, intent(in) :: commLamVis, commEddyVis
!
!      Local variables.
!
       integer :: comm, size, ierr

       integer(kind=intType) :: nVar, start, end
       integer(kind=intType) :: sps, nn, mm, ll
       integer(kind=intType) :: nInter, nDonor, nHalo

       integer(kind=intType), dimension(:),   pointer :: bd, bh
       integer(kind=intType), dimension(:),   pointer :: indListD
       integer(kind=intType), dimension(:),   pointer :: nintD
       integer(kind=intType), dimension(:,:), pointer :: indListH
       integer(kind=intType), dimension(:,:), pointer :: indD, indH

       real(kind=realType), dimension(:),     pointer :: wd
       real(kind=realType), dimension(:,:),   pointer :: wh
       real(kind=realType), dimension(:,:,:), pointer :: matD, matH

       logical :: commVel
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if either not a steady problem is solved or
       ! if no sliding mesh interfaces are present.

       if(equationMode /= steady .or. nInterfaceGroups == 0) return

       ! Set the logical whether or not the velocity components are
       ! communicated. If they have to be communicated, possibly correct
       ! the start and end indices.

       commVel = .false.
       if(startIn <= ivz .and. endIn >= ivx) commVel = .true.

       start = startIn
       end   = endIn

       if( commVel ) then
         if(start > ivx) start = ivx
         if(end   < ivz) end   = ivz
       endif

       ! Determine the number of variables per cell to be set.

       nVar = max(0_intType,(end - start + 1))
       if( commPressure ) nVar = nVar + 1
       if( commVarGamma ) nVar = nVar + 1
       if( commLamVis )   nVar = nVar + 1
       if( commEddyVis )  nVar = nVar + 1

       if(nVar == 0) return

       ! Loop over the number of spectral solutions. For a steady state
       ! computation this will be one, but it is added for consistency.

       spectralLoop: do sps=1,nTimeIntervalsSpectral

         ! Loop over the number of interface groups.

         interfaceLoop: do nn=1,nInterfaceGroups

           ! If i do not contribute to the current "color" continue
           ! to the next color.

           if(.not. myInterfaces(nn)%procContributes) cycle

           ! Store the communicator a bit easier.

           comm = myInterfaces(nn)%commSlide

           ! Loop over the two sides of the interface.

           sideLoop: do ll=1,2

             ! Some abbreviations to make the code more readable.

             nInter = commPatternMixing(level,nn,ll)%nInter
             nDonor = commPatternMixing(level,nn,ll)%nDonor
             nHalo  = commPatternMixing(level,nn,ll)%nHalo

             wd       => commPatternMixing(level,nn,ll)%weightDonor
             matD     => commPatternMixing(level,nn,ll)%rotMatDonor
             bd       => commPatternMixing(level,nn,ll)%blockDonor
             nintD    => commPatternMixing(level,nn,ll)%nIntervalsDonor
             indListD => commPatternMixing(level,nn,ll)%indListDonor

             wh       => commPatternMixing(level,nn,ll)%weightHalo
             matH     => commPatternMixing(level,nn,ll)%rotMatHalo
             bh       => commPatternMixing(level,nn,ll)%blockHalo
             indListH => commPatternMixing(level,nn,ll)%indListHalo

             ! Loop over the number of halo layers to be set. This is
             ! either 1 or 2.

             layerLoop: do mm=1,nLayers

               ! Set the pointers for the donor and halo indices to make
               ! the code more readable.

               indD => commPatternMixing(level,nn,ll)%indD(:,:,mm)
               indH => commPatternMixing(level,nn,ll)%indH(:,:,mm)

               ! Determine the local contribution to interpolation
               ! list. This data is stored in sendBuffer.

               call localPartMixingPlane(sendBuffer)

               ! Determine the sum of the local contribution for the
               ! processor group of this sliding interface. Store the
               ! global sum in recvBuffer.

               size = nVar*nInter
               call mpi_allreduce(sendBuffer, recvBuffer, size, &
                                  sumb_real, mpi_sum, comm, ierr)

               ! Determine the local mixing plane halo cells.

               call setMixingPlaneHalos(recvBuffer)

             enddo layerLoop
           enddo sideLoop
         enddo interfaceLoop
       enddo spectralLoop

       !=================================================================

       contains

         !===============================================================

         subroutine localPartMixingPlane(buffer)
!
!        ****************************************************************
!        *                                                              *
!        * LocalPartMixingPlane stores the local contribution to the    *
!        * currently active mixing plane into buffer.                   *
!        *                                                              *
!        ****************************************************************
!
         use block
         implicit none
!
!        Subroutine arguments.
!
         real(kind=realType), dimension(nVar,*), intent(out) :: buffer
!
!        Local variables.
!
         integer(kind=intType) :: i, j, k, b, ii, jj, nn, mm
         integer(kind=intType) :: iax, irad, itheta, iof, iof2, start2

         real(kind=realType) :: vx, vy, vz, vax, vrad, vtheta
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Initialize the buffer to zero; it will contain the sum of the
         ! local values.

         do i=1,nInter
           do j=1,nVar
             buffer(j,i) = zero
           enddo
         enddo

         ! First treat the velocities, because these require a
         ! transformation from the global cartesian to the local
         ! cylindrical frame.

         testCommVel: if( commVel ) then

           ! Determine the indices in the buffer to store the three
           ! cylindrical velocity components.

           iax    = ivx - start + 1
           irad   = iax + 1
           itheta = irad + 1

           ! Loop over the number of donors.

           velDonorLoop: do nn=1,nDonor

             ! Store the indices of the donor a bit easier.

             b = bd(nn)
             i = indD(nn,1)
             j = indD(nn,2)
             k = indD(nn,3)

             ! Store the 3 cartesian velocity components a bit easier.

             vx = flowDoms(b,level,sps)%w(i,j,k,ivx)
             vy = flowDoms(b,level,sps)%w(i,j,k,ivy)
             vz = flowDoms(b,level,sps)%w(i,j,k,ivz)

             ! Transform the velocity to the components of the local
             ! cylindrical coordinate system

             vax    = matD(nn,1,1)*vx + matD(nn,1,2)*vy &
                    + matD(nn,1,3)*vz
             vrad   = matD(nn,2,1)*vx + matD(nn,2,2)*vy &
                    + matD(nn,2,3)*vz
             vtheta = matD(nn,3,1)*vx + matD(nn,3,2)*vy &
                    + matD(nn,3,3)*vz

             ! Loop over the number of intervals to which this
             ! donor contributes.

             do jj=(nintD(nn-1)+1),nintD(nn)

               ! Store the index in the buffer a bit easier and
               ! update the velocity entries in the buffer.

               ii = indListD(jj)

               buffer(iax,ii)    = buffer(iax,ii)    + wd(jj)*vax
               buffer(irad,ii)   = buffer(irad,ii)   + wd(jj)*vrad
               buffer(itheta,ii) = buffer(itheta,ii) + wd(jj)*vtheta
             enddo

           enddo velDonorLoop
         endif testCommVel

         ! Determine the offset such that the index start will be stored
         ! at the first position in the buffer. Also set start2, needed
         ! to store the variables with an index larger than the velocity.
         ! Iof2 is the offset to store the other variables than the
         ! working ones at the correct location in buffer

         iof    = 1 - start
         iof2   = max(0_intType,(end - start + 1))
         start2 = max(start,ivz+1)

         ! Store the other variables in the buffer.

         donorLoop: do nn=1,nDonor

           ! Store the indices of the donor a bit easier.

           b = bd(nn)
           i = indD(nn,1)
           j = indD(nn,2)
           k = indD(nn,3)

           ! Loop over the number of intervals to which this
           ! donor contributes.

           contribLoop: do jj=(nintD(nn-1)+1),nintD(nn)

             ! Store the index in the list a bit easier.

             ii = indListD(jj)

             ! The working variables. Be careful not to overwrite the
             ! already computed velocity components.

             do mm=start,(ivx-1)
               buffer(iof+mm,ii) = buffer(iof+mm,ii) &
                                 + wd(jj)*flowDoms(b,level,sps)%w(i,j,k,mm)
             enddo

             do mm=start2,end
               buffer(iof+mm,ii) = buffer(iof+mm,ii) &
                                 + wd(jj)*flowDoms(b,level,sps)%w(i,j,k,mm)
             enddo

             ! The other variables. Note that for gamma and rlv the level
             ! is 1, because these variables are only allocated on the
             ! finest grid; they are not considered true mg variables in
             ! the sense that they depend on other variables.

             mm = iof2

             if( commPressure ) then
               mm = mm + 1
               buffer(mm,ii) = buffer(mm,ii) &
                             + wd(jj)*flowDoms(b,level,sps)%p(i,j,k)
             endif

             if( commVarGamma ) then
               mm = mm + 1
               buffer(mm,ii) = buffer(mm,ii) &
                             + wd(jj)*flowDoms(b,1,sps)%gamma(i,j,k)
             endif

             if( commLamVis ) then
               mm = mm + 1
               buffer(mm,ii) = buffer(mm,ii) &
                             + wd(jj)*flowDoms(b,1,sps)%rlv(i,j,k)
             endif

             if( commEddyVis ) then
               mm = mm + 1
               buffer(mm,ii) = buffer(mm,ii) &
                             + wd(jj)*flowDoms(b,level,sps)%rev(i,j,k)
             endif

           enddo contribLoop
         enddo donorLoop

         end subroutine localPartMixingPlane

         !===============================================================

         subroutine setMixingPlaneHalos(buffer)
!
!        ****************************************************************
!        *                                                              *
!        * SetMixingPlaneHalos set the values in the halo cells         *
!        * adjacent to the active mixing plane. The variables are       *
!        * interpolated from buffer.                                    *
!        *                                                              *
!        ****************************************************************
!
         use block
         implicit none
!
!        Subroutine arguments.
!
         real(kind=realType), dimension(nVar,*), intent(in) :: buffer
!
!        Local variables.
!
         integer(kind=intType) :: i, j, k, b, nn, mm, i1, i2
         integer(kind=intType) :: iof, iof2

         real(kind=realType) :: w1, w2, vax, vrad, vtheta, vx, vy, vz
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Determine the offset iof such that the index start will get
         ! its data from index 1 in the buffer. Also set the offset iof2
         ! for the "other" variables.

         iof  = 1 - start
         iof2 = max(0_intType,(end - start + 1))

         ! Loop over the number of halos to be set.

         haloLoop: do nn=1,nHalo

           ! Store the indices of the halo as well as the indices and
           ! the weights in the list a bit easier.

           b = bh(nn)
           i = indH(nn,1)
           j = indH(nn,2)
           k = indH(nn,3)

           i1 = indListH(nn,1)
           i2 = indListH(nn,2)

           w1 = wh(nn,1)
           w2 = wh(nn,2)

           ! Loop over the working variables to be set. Note that the
           ! velocity components of the local cylindrical coordinate
           ! system are interpolated. These are transformed to the
           ! cartesian frame later on.

           do mm=start,end
             flowDoms(b,level,sps)%w(i,j,k,mm) = w1*buffer(iof+mm,i1) &
                                               + w2*buffer(iof+mm,i2)
           enddo

           ! The other variables. Again level = 1 for gamma and rlv,
           ! see the explanation earlier.

           mm = iof2

           if( commPressure ) then
             mm = mm + 1
             flowDoms(b,level,sps)%p(i,j,k) = w1*buffer(mm,i1) &
                                            + w2*buffer(mm,i2)
           endif

           if( commVarGamma ) then
             mm = mm + 1
             flowDoms(b,1,sps)%gamma(i,j,k) = w1*buffer(mm,i1) &
                                            + w2*buffer(mm,i2)
           endif

           if( commLamVis ) then
             mm = mm + 1
             flowDoms(b,1,sps)%rlv(i,j,k) = w1*buffer(mm,i1) &
                                          + w2*buffer(mm,i2)
           endif

           if( commEddyVis ) then
             mm = mm + 1
             flowDoms(b,level,sps)%rev(i,j,k) = w1*buffer(mm,i1) &
                                              + w2*buffer(mm,i2)
           endif

         enddo haloLoop

         ! Transform the velocities to the cartesian components if
         ! velocities must be communicated.

         testCommVel: if( commVel ) then

           ! Loop over the number of halo cells to be set.

           velHaloLoop: do nn=1,nHalo

             ! Store the indices of the halo a bit easier.

             b = bh(nn)
             i = indH(nn,1)
             j = indH(nn,2)
             k = indH(nn,3)

             ! Compute the cartesian components of the velocity.

             vax    = flowDoms(b,level,sps)%w(i,j,k,ivx)
             vrad   = flowDoms(b,level,sps)%w(i,j,k,ivy)
             vtheta = flowDoms(b,level,sps)%w(i,j,k,ivz)

             vx = matH(nn,1,1)*vax    + matH(nn,1,2)*vrad &
                + matH(nn,1,3)*vtheta
             vy = matH(nn,2,1)*vax    + matH(nn,2,2)*vrad &
                + matH(nn,2,3)*vtheta
             vz = matH(nn,3,1)*vax    + matH(nn,3,2)*vrad &
                + matH(nn,3,3)*vtheta

             flowDoms(b,level,sps)%w(i,j,k,ivx) = vx
             flowDoms(b,level,sps)%w(i,j,k,ivy) = vy
             flowDoms(b,level,sps)%w(i,j,k,ivz) = vz

           enddo velHaloLoop
         endif testCommVel

         end subroutine setMixingPlaneHalos

       end subroutine whaloMixing
