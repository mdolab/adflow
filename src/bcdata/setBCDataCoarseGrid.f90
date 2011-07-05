!
!      ******************************************************************
!      *                                                                *
!      * File:          setBCDataCoarseGrid.f90                         *
!      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
!      * Starting date: 07-07-2004                                      *
!      * Last modified: 08-09-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setBCDataCoarseGrid
!
!      ******************************************************************
!      *                                                                *
!      * setBCDataCoarseGrid determines the boundary condition info     *
!      * on the coarse grid from the known info on the fine grid. It    *
!      * will be stored in the BCData arrays of flowDoms.               *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use flowVarRefState
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, l, sps
       integer(kind=intType) :: iBeg, jBeg, iEnd, jEnd, iiMax, jjMax
       integer(kind=intType) :: nLevels, level, levm1

       integer(kind=intType), dimension(:,:), pointer :: iFine, jFine

       real(kind=realType) :: var

       real(kind=realType), dimension(3) :: dir
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of grid levels.

       nLevels = ubound(flowDoms,2)

       ! Loop over the coarser grid levels. It is assumed that the
       ! bc data of groundLevel is set correctly.

       coarseLevelLoop: do level=(groundLevel+1),nLevels

         ! Store the fine grid level a bit easier.

         levm1 = level - 1

         ! Loop over the number of spectral solutions and local blocks.

         spectralLoop: do sps=1,nTimeIntervalsSpectral
           domainsLoop: do i=1,nDom

             ! Set the pointers to the coarse block.

             call setPointers(i, level, sps)

             ! Loop over the boundary subfaces and interpolate the
             ! prescribed boundary data for this grid level.

             bocoLoop: do j=1,nBocos

               ! Determine the block face on which the subface is
               ! located and set some multigrid variables accordingly.

               select case (BCFaceID(j))

                 case (iMin,iMax)
                   iiMax = jl; jjMax = kl
                   iFine => mgJFine; jFine => mgKFine

                 case (jMin,jMax)
                   iiMax = il; jjMax = kl
                   iFine => mgIFine; jFine => mgKFine

                 case (kMin,kMax)
                   iiMax = il; jjMax = jl
                   iFine => mgIFine; jFine => mgJFine

               end select

               ! Abbreviate the size of the subface a bit easier.

               iBeg = BCData(j)%icBeg; iEnd = BCData(j)%icEnd
               jBeg = BCData(j)%jcBeg; jEnd = BCData(j)%jcEnd

               ! Copy the subsonic boundary conditions treatment.

               BCData(j)%subsonicInletTreatment = &
                 flowDoms(i,levm1,sps)%BCData(j)%subsonicInletTreatment

               ! Interpolate the data for the possible prescribed boundary
               ! data.

               call interpolateBcData(BCData(j)%TNS_Wall, &
                   flowDoms(i,levm1,sps)%BCData(j)%TNS_Wall)

               call interpolateBcData(BCData(j)%ptInlet, &
                   flowDoms(i,levm1,sps)%BCData(j)%ptInlet)
               call interpolateBcData(BCData(j)%ttInlet, &
                   flowDoms(i,levm1,sps)%BCData(j)%ttInlet)
               call interpolateBcData(BCData(j)%flowXdirInlet, &
                   flowDoms(i,levm1,sps)%BCData(j)%flowXdirInlet)
               call interpolateBcData(BCData(j)%flowYdirInlet, &
                   flowDoms(i,levm1,sps)%BCData(j)%flowYdirInlet)
               call interpolateBcData(BCData(j)%flowZdirInlet, &
                   flowDoms(i,levm1,sps)%BCData(j)%flowZdirInlet)

               call interpolateBCVecData(BCData(j)%turbInlet, &
                       flowDoms(i,levm1,sps)%BCData(j)%turbInlet, &
                       nt1, nt2)

               call interpolateBcData(BCData(j)%rho,  &
                   flowDoms(i,levm1,sps)%BCData(j)%rho)
               call interpolateBcData(BCData(j)%velx, &
                   flowDoms(i,levm1,sps)%BCData(j)%velx)
               call interpolateBcData(BCData(j)%vely, &
                   flowDoms(i,levm1,sps)%BCData(j)%vely)
               call interpolateBcData(BCData(j)%velz, &
                   flowDoms(i,levm1,sps)%BCData(j)%velz)
               call interpolateBcData(BCData(j)%ps,   &
                   flowDoms(i,levm1,sps)%BCData(j)%ps)

               ! Some additional variables should be computed/corrected
               ! for some boundary conditions. Determine the type of
               ! boundary condition.

               if((BCType(j) == SubsonicInflow .and.                         &
                   BCData(j)%subsonicInletTreatment == totalConditions) .or. &
                  BCType(j) == DomainInterfaceTotal) then

                 ! Total conditions are specified for subsonic inflow
                 ! or domain interfaces.

                 ! Compute the total enthalpy and make
                 ! sure that the unit vector is a unit vector.

                 ! Loop over the faces of the subface.

                 do l=jBeg,jEnd
                   do k=iBeg,iEnd

                     ! Compute the total enthalpy.

                     call computeHtot(BCData(j)%ttInlet(k,l), &
                                      BCData(j)%htInlet(k,l))

                     ! Flow direction.

                     dir(1) = BCData(j)%flowXdirInlet(k,l)
                     dir(2) = BCData(j)%flowYdirInlet(k,l)
                     dir(3) = BCData(j)%flowZdirInlet(k,l)

                     var = one/max(eps,sqrt(dir(1)**2 + dir(2)**2 &
                         +                  dir(3)**2))

                     BCData(j)%flowXdirInlet(k,l) = var*dir(1)
                     BCData(j)%flowYdirInlet(k,l) = var*dir(2)
                     BCData(j)%flowZdirInlet(k,l) = var*dir(3)

                   enddo
                 enddo

               endif

             enddo bocoLoop
           enddo domainsLoop
         enddo spectralLoop
       enddo coarseLevelLoop

       !=================================================================

       contains

         !===============================================================

         subroutine interpolateBcData(varCoarse, varFine)
!
!        ****************************************************************
!        *                                                              *
!        * InterpolateBcData interpolates the given data array from   *
!        * the fine to the coarse grid. Of course only if the fine      *
!        * array is associated with some data.                          *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         real(kind=realType), dimension(:,:), pointer :: varCoarse
         real(kind=realType), dimension(:,:), pointer :: varFine
!
!        Local variables.
!
         integer(kind=intType) :: i, j, if1, if2, jf1, jf2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Check if varFine is associated to data. If not return.

         if(.not. associated(varFine)) return

         ! Loop over the faces of the given subface.
         ! First the j-direction.

         do j=jBeg,jEnd

           ! Determine the two children in this direction. Take care of
           ! the halo's, as this info is only available for owned cells.

           if(j < 2) then
             jf1 = 1; jf2 = 1
           else if(j > jjMax) then
             jf1 = jFine(jjMax,2) +1; jf2 = jf1
           else
             jf1 = jFine(j,1); jf2 = jFine(j,2)
           endif

           ! Loop in the i-direction.

           do i=iBeg,iEnd

             ! Determine the two children in this direction.
             ! Same story as in j-direction.

             if(i < 2) then
               if1 = 1; if2 = 1
             else if(i > iiMax) then
               if1 = iFine(iiMax,2) +1; if2 = if1
             else
               if1 = iFine(i,1); if2 = iFine(i,2)
             endif

             ! Compute the coarse grid data as the average of the
             ! 4 fine grid values.

             varCoarse(i,j) = fourth*(varFine(if1,jf1) &
                            +         varFine(if2,jf1) &
                            +         varFine(if1,jf2) &
                            +         varFine(if2,jf2))
           enddo
         enddo

         end subroutine interpolateBcData

         !===============================================================

         subroutine interpolateBCVecData(varCoarse, varFine, &
                                         nstart, nend)
!
!        ****************************************************************
!        *                                                              *
!        * interpolateBCVecData interpolates the given data array       *
!        * from the fine to the coarse grid. Of course only if the fine *
!        * array is associated with some data.                          *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer(kind=intType), intent(in) :: nstart, nend

         real(kind=realType), dimension(:,:,:), pointer :: varCoarse
         real(kind=realType), dimension(:,:,:), pointer :: varFine
!
!        Local variables.
!
         integer(kind=intType) :: nn, i, j, if1, if2, jf1, jf2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Check if varFine is associated to data. If not return.

         if(.not. associated(varFine)) return

         ! Loop over the faces of the given subface.
         ! First the j-direction.

         do j=jBeg,jEnd

           ! Determine the two children in this direction. Take care of
           ! the halo's, as this info is only available for owned cells.

           if(j < 2) then
             jf1 = 1; jf2 = 1
           else if(j > jjMax) then
             jf1 = jFine(jjMax,2) +1; jf2 = jf1
           else
             jf1 = jFine(j,1); jf2 = jFine(j,2)
           endif

           ! Loop in the i-direction.

           do i=iBeg,iEnd

             ! Determine the two children in this direction.
             ! Same story as in j-direction.

             if(i < 2) then
               if1 = 1; if2 = 1
             else if(i > iiMax) then
               if1 = iFine(iiMax,2) +1; if2 = if1
             else
               if1 = iFine(i,1); if2 = iFine(i,2)
             endif

             ! Compute the coarse grid data as the average of the
             ! 4 fine grid values.

             do nn=nstart,nend
               varCoarse(i,j,nn) = fourth*(varFine(if1,jf1,nn) &
                                 +         varFine(if2,jf1,nn) &
                                 +         varFine(if1,jf2,nn) &
                                 +         varFine(if2,jf2,nn))
             enddo
           enddo
         enddo

         end subroutine interpolateBCVecData

       end subroutine setBCDataCoarseGrid
