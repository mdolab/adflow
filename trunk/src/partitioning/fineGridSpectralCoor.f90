!
!      ******************************************************************
!      *                                                                *
!      * File:          fineGridSpectralCoor.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-23-2004                                      *
!      * Last modified: 10-06-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine fineGridSpectralCoor
!
!      ******************************************************************
!      *                                                                *
!      * fineGridSpectralCoor computes the coordinates of all but       *
!      * the first spectral solution from the known coordinates of the  *
!      * first time instance.                                           *
!      *                                                                *
!      ******************************************************************
!
       use block
       use inputPhysics
       use inputTimeSpectral
       use IOModule
       use iteration
       use monitor
       use section
       use partitionMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, ll, i, j, k
       integer(kind=intType) :: il, jl, kl

       real(kind=realType), dimension(nSections) :: dt, t
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! This routine is only used for the spectral solutions. Return
       ! immediately if a different mode is solved.

       if(equationMode /= timeSpectral) return

       ! Also return immediately if all coordinates were already read
       ! from the grid file.

       if(.not. interpolSpectral) return
!
!      ******************************************************************
!      *                                                                *
!      * Step 1. Perform a rigid body motion of the coordinates of the  *
!      *         1st time instance to the other instances.              *
!      *                                                                *
!      ******************************************************************
!
       ! Set currentLevel to 1, such that updateCoorFineMesh
       ! updates the coordinates of the correct level.

       currentLevel = 1

       ! Determine the delta t for every section. Remember it is possible
       ! that every section has a different periodic time.

       do nn=1,nSections
         dt(nn) = sections(nn)%timePeriod &
                / real(nTimeIntervalsSpectral,realType)
       enddo

       ! Loop over the number of spectral solutions, starting at 2,
       ! because the first is already known.

       timeUnsteady = zero
       spectralLoop: do ll=2,nTimeIntervalsSpectral

         ! Set the owned coordinates to the coordinates of the first
         ! solution such that updateCoorFineMesh can modify them.

         do nn=1,nDom

           do k=1,flowDoms(nn,1,1)%kl
             do j=1,flowDoms(nn,1,1)%jl
               do i=1,flowDoms(nn,1,1)%il
                 flowDoms(nn,1,ll)%x(i,j,k,1) = flowDoms(nn,1,1)%x(i,j,k,1)
                 flowDoms(nn,1,ll)%x(i,j,k,2) = flowDoms(nn,1,1)%x(i,j,k,2)
                 flowDoms(nn,1,ll)%x(i,j,k,3) = flowDoms(nn,1,1)%x(i,j,k,3)
               enddo
             enddo
           enddo

         enddo

         ! Compute the corresponding times for this spectral solution
         ! and call updateCoorFineMesh to determine the coordinates.

         do nn=1,nSections
           t(nn) = (ll-1)*dt(nn)
         enddo

         call updateCoorFineMesh(t, ll)

       enddo spectralLoop

       ! Return if only one grid has been read.

       if(nGridsRead == 1) return
!
!      ******************************************************************
!      *                                                                *
!      * Step 2. Multiple grids have been read, but the number is not   *
!      *         equal to the number of time instances used in the      *
!      *         computation. As multiple grids have been read this     *
!      *         means that a time spectral computation on a deforming  *
!      *         mesh is performed. Therefore the deformations,         *
!      *         relative to the first grid, must be interpolated.      *
!      *                                                                *
!      ******************************************************************
!
       ! First allocate the memory of IOVar(..,1)%w.
       ! This will serve as temporary storage for the coordinates of
       ! the 1st spectral instance, because those will be used in the
       ! call to updateCoorFineMesh. In this way this routine can be
       ! used without modification.

       do nn=1,nDom

         kl = flowDoms(nn,1,1)%kl
         jl = flowDoms(nn,1,1)%jl
         il = flowDoms(nn,1,1)%il

         allocate(IOVar(nn,1)%w(il,jl,kl,3), stat=ierr)
         if(ierr /= 0)                            &
           call terminate("fineGridSpectralCoor", &
                          "Memory allocation failure for &
                          &IOVar(nn,1)%w")

         do k=1,kl
           do j=1,jl
             do i=1,il
               IOVar(nn,1)%w(i,j,k,1) = flowDoms(nn,1,1)%x(i,j,k,1)
               IOVar(nn,1)%w(i,j,k,2) = flowDoms(nn,1,1)%x(i,j,k,2)
               IOVar(nn,1)%w(i,j,k,3) = flowDoms(nn,1,1)%x(i,j,k,3)
             enddo
           enddo
         enddo

       enddo

       ! Determine the delta t for every section. Remember it is possible
       ! that every section has a different periodic time.

       do nn=1,nSections
         dt(nn) = sections(nn)%timePeriod &
                / real(nGridsRead,realType)
       enddo

       ! Loop over the number of spectral solutions read, starting at
       ! 2, and determine the displacements relative to the rigid body
       ! motion of the grid of the 1st time instance.

       timeUnsteady = zero
       spectralLoopRead: do ll=2,nGridsRead

         ! Compute the corresponding times for this spectral solution
         ! and call updateCoorFineMesh to determine the coordinates.

         do nn=1,nSections
           t(nn) = (ll-1)*dt(nn)
         enddo

         call updateCoorFineMesh(t, 1_intType)

         ! Determine the relative displacements for this time instance
         ! and initialize flowDoms(nn,1,1) for the next round.

         do nn=1,nDom

           do k=1,flowDoms(nn,1,1)%kl
             do j=1,flowDoms(nn,1,1)%jl
               do i=1,flowDoms(nn,1,1)%il
                 IOVar(nn,ll)%w(i,j,k,1) = IOVar(nn,ll)%w(i,j,k,1)     &
                                         - flowDoms(nn,1,1)%x(i,j,k,1)
                 IOVar(nn,ll)%w(i,j,k,2) = IOVar(nn,ll)%w(i,j,k,2)     &
                                         - flowDoms(nn,1,1)%x(i,j,k,2)
                 IOVar(nn,ll)%w(i,j,k,3) = IOVar(nn,ll)%w(i,j,k,3)     &
                                         - flowDoms(nn,1,1)%x(i,j,k,3)

                 flowDoms(nn,1,1)%x(i,j,k,1) = IOVar(nn,1)%w(i,j,k,1)
                 flowDoms(nn,1,1)%x(i,j,k,2) = IOVar(nn,1)%w(i,j,k,2)
                 flowDoms(nn,1,1)%x(i,j,k,3) = IOVar(nn,1)%w(i,j,k,3)
               enddo
             enddo
           enddo

         enddo

       enddo spectralLoopRead

       ! The coordinates of IOVar now contain the relative
       ! displacements compared to the rigid body motion of the first
       ! time instance, except for the first time instance.
       ! Set these to zero.

       do nn=1,nDom
         IOVar(nn,1)%w = zero
       enddo

       ! Interpolate the displacements and add them to the currently
       ! stored coordinates.

       call terminate("fineGridSpectralCoor", &
                      "Arti should do this interpolation stuff")

       ! Release the memory of the variable w in IOVar.

       do nn=1,nDom
         do ll=1,nGridsRead
           deallocate(IOVar(nn,ll)%w, stat=ierr)
           if(ierr /= 0)                            &
             call terminate("fineGridSpectralCoor", &
                            "Deallocation failure for IOVar%w")
         enddo
       enddo

       end subroutine fineGridSpectralCoor
