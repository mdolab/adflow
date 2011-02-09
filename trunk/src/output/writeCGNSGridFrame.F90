!
!      ******************************************************************
!      *                                                                *
!      * File:          writeCGNSGridFrame.F90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-21-2004                                      *
!      * Last modified: 10-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeCGNSGridFrame(cgnsZone, ind)
!
!      ******************************************************************
!      *                                                                *
!      * writeCGNSGridFrame writes the framework for the grid file      *
!      * gridNames(ind) using the information stored in the module      *
!      * cgnsGrid. Basically all information but the coordinates is     *
!      * written by this routine.                                       *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use su_cgns
       use outputMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)  :: ind
       integer, dimension(*), intent(out) :: cgnsZone

#ifdef USE_NO_CGNS
       call terminate("writeCGNSGridFrame", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
!
!      Local variables.
!
       integer :: ierr, ii, jj, cgnsInd, cgnsBase
       integer, dimension(9)   :: sizes
       integer, dimension(3,2) :: zoneRange, donorRange
       integer, dimension(3)   :: transform
       integer, dimension(:,:), allocatable :: donorData

       integer(kind=intType) :: nn, mm, ll, i, j, k
       integer(kind=intType) :: s1, s2, s3

       real, dimension(3) :: rotCenter, rotRate, translation

       real(kind=realType) :: LRefInv

       character(len=maxStringLen)   :: errorMessage

       type(cgnsBcDatasetType), pointer, dimension(:) :: dataSet
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Open the CGNS file for writing and check if it went okay.
       ! Store the file index afterwards.

       call cg_open_f(gridFileNames(ind), mode_write, cgnsInd, ierr)
       if(ierr /= all_ok) then
         write(errorMessage,*) "File ", trim(gridfileNames(ind)), &
                               " could not be opened by cgns &
                               &for writing"
         call terminate("writeCGNSGridFrame", errorMessage)
       endif

       fileIDs(ind) = cgnsInd

       ! Create the base. Copy the cell and physical dimensions and
       ! store the base ID for this index.

       call cg_base_write_f(cgnsInd, cgnsBaseName, cgnsCelldim, &
                            cgnsPhysdim, cgnsBase, ierr)
       if(ierr /= all_ok)                     &
         call terminate("writeCGNSGridFrame", &
                        "Something wrong when calling cg_base_write_f")

       cgnsBases(ind) = cgnsBase
!
!      ******************************************************************
!      *                                                                *
!      * Write the family info.                                         *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of families.

       familyLoop: do nn=1,cgnsNfamilies

         ! Create the family node.

         call cg_family_write_f(cgnsInd, cgnsBase, &
                                cgnsFamilies(nn)%familyName, ii, ierr)
         if(ierr /= all_ok)                     &
           call terminate("writeCGNSGridFrame", &
                          "Something wrong when calling &
                          &cg_family_write_f")

         ! Write the family BC, if this is present.

         if(cgnsFamilies(nn)%BCTypeCGNS /= Null) then

           call cg_fambc_write_f(cgnsInd, cgnsBase, ii,   &
                                 cgnsFamilies(nn)%bcName, &
                                 cgnsFamilies(nn)%BCTypeCGNS, jj, ierr)
           if(ierr /= all_ok)                     &
             call terminate("writeCGNSGridFrame", &
                            "Something wrong when calling &
                            &cg_fambc_write_f")

           ! If the boundary condition is UserDefined add the
           ! description what type of user defined BC it is.

           if(cgnsFamilies(nn)%BCTypeCGNS == UserDefined) then

             ! Ultimately you would like to create the
             ! UserDefinedData_t as a subnode of the family boundary
             ! condition node. However, at the moment CGNS does not
             ! allow this and therefore it is put one level higher.
             ! As only 1 boundary condition per family is allowed,
             ! this is not really problem.

             call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                            "Family_t", ii, "end")
             if(ierr /= all_ok)                     &
               call terminate("writeCGNSGridFrame", &
                              "Something wrong when calling cg_goto_f")

             call cg_user_data_write_f(cgnsFamilies(nn)%userDefinedName, &
                                       ierr)
             if(ierr /= all_ok)                     &
               call terminate("writeCGNSGridFrame", &
                              "Something wrong when calling &
                              &cg_user_data_write_f")
           endif

         endif
       enddo familyLoop
!
!      ******************************************************************
!      *                                                                *
!      * Write all the zone info, except the coordinates.               *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of zones in the original grid.

       zoneLoop: do nn=1,cgnsNDom

         ! Store the inverse of the scaling factor to meters.

         LRefInv = one/cgnsDoms(nn)%LRef

         ! Store the dimensions of the zone in sizes and create the zone.

         sizes(1) = cgnsDoms(nn)%il
         sizes(2) = cgnsDoms(nn)%jl
         sizes(3) = cgnsDoms(nn)%kl
         sizes(4) = cgnsDoms(nn)%nx
         sizes(5) = cgnsDoms(nn)%ny
         sizes(6) = cgnsDoms(nn)%nz
         sizes(7) = 0
         sizes(8) = 0
         sizes(9) = 0

         call cg_zone_write_f(cgnsInd, cgnsBase,                   &
                              cgnsDoms(nn)%zoneName, sizes,        &
                              cgnsDoms(nn)%zoneType, cgnsZone(nn), &
                              ierr)
         if(ierr /= all_ok)                     &
           call terminate("writeCGNSGridFrame", &
                          "Something wrong when calling &
                          &cg_zone_write_f")

         ! Go to the current zone. Needed when family and/or rotating
         ! frame info must be written.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", &
                        cgnsZone(nn), "end")
         if(ierr /= all_ok) &
           call terminate("writeCGNSGridFrame", &
                          "Something wrong when calling cg_goto_f")

         ! Check if the zone belongs to a family. If so, write the
         ! family name.

         mm = cgnsDoms(nn)%familyID
         if(mm > 0) then

           call cg_famname_write_f(cgnsFamilies(mm)%familyName, ierr)
           if(ierr /= all_ok) &
             call terminate("writeCGNSGridFrame", &
                            "Something wrong when calling &
                            &cg_famname_write_f")
         endif

         ! Write the rotating frame info, if the zone is rotating.

         if( cgnsDoms(nn)%rotatingFrameSpecified ) then

           ! Convert the rotation rate to degrees per second and store
           ! it in a single precision array.

           rotRate(1) = cgnsDoms(nn)%rotRate(1)*180.0_realType/pi
           rotRate(2) = cgnsDoms(nn)%rotRate(2)*180.0_realType/pi
           rotRate(3) = cgnsDoms(nn)%rotRate(3)*180.0_realType/pi

           ! Convert the rotation center to the original units
           ! and also in single precision.

           rotCenter(1) = LRefInv*cgnsDoms(nn)%rotCenter(1)
           rotCenter(2) = LRefInv*cgnsDoms(nn)%rotCenter(2)
           rotCenter(3) = LRefInv*cgnsDoms(nn)%rotCenter(3)

           ! Write the rotation rate and rotation center.

           call cg_rotating_write_f(rotRate, rotCenter, ierr)
           if(ierr /= all_ok) &
             call terminate("writeCGNSGridFrame", &
                            "Something wrong when calling &
                            &cg_rotating_write_f")

           ! Write the units of the rotation rate of the
           ! rotating frame.

           call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t",        &
                          cgnsZone(nn), "RotatingCoordinates_t", 1, &
                          "DataArray_t", 2, "end")
           if(ierr /= all_ok)                     &
             call terminate("writeCGNSGridFrame", &
                            "Something wrong when calling cg_goto_f")

           call cg_dataclass_write_f(Dimensional, ierr)
           if(ierr /= all_ok)                     &
             call terminate("writeCGNSGridFrame", &
                            "Something wrong when calling &
                            &cg_dataclass_write_f")

           call cg_units_write_f(Null, Null, Second, Null, &
                                 Degree, ierr)
           if(ierr /= all_ok)                     &
             call terminate("writeCGNSGridFrame", &
                            "Something wrong when calling &
                            &cg_units_write_f")
         endif

         ! Loop over all 1 to 1 connectivities of the block and
         ! write the data.

         loop1to1: do mm=1,cgnsDoms(nn)%n1to1

           ! Store the range of the subface in zoneRange and
           ! the range of the donor in donorRange.

           zoneRange(1,1) = cgnsDoms(nn)%conn1to1(mm)%iBeg
           zoneRange(2,1) = cgnsDoms(nn)%conn1to1(mm)%jBeg
           zoneRange(3,1) = cgnsDoms(nn)%conn1to1(mm)%kBeg

           zoneRange(1,2) = cgnsDoms(nn)%conn1to1(mm)%iEnd
           zoneRange(2,2) = cgnsDoms(nn)%conn1to1(mm)%jEnd
           zoneRange(3,2) = cgnsDoms(nn)%conn1to1(mm)%kEnd

           donorRange(1,1) = cgnsDoms(nn)%conn1to1(mm)%diBeg
           donorRange(2,1) = cgnsDoms(nn)%conn1to1(mm)%djBeg
           donorRange(3,1) = cgnsDoms(nn)%conn1to1(mm)%dkBeg

           donorRange(1,2) = cgnsDoms(nn)%conn1to1(mm)%diEnd
           donorRange(2,2) = cgnsDoms(nn)%conn1to1(mm)%djEnd
           donorRange(3,2) = cgnsDoms(nn)%conn1to1(mm)%dkEnd

           ! Check whether the subface is periodic or not.

           periodicTest: if( cgnsDoms(nn)%conn1to1(mm)%periodic ) then

             ! Subface is periodic. Due to the current limitations in
             ! cgns it is not possible to write this info as a 1 to 1
             ! subface and the general connectivity must be used.

             ! First allocate the memory for donorData.

             ll = (abs(donorRange(3,2) - donorRange(3,1)) + 1) &
                * (abs(donorRange(2,2) - donorRange(2,1)) + 1) &
                * (abs(donorRange(1,2) - donorRange(1,1)) + 1)

             allocate(donorData(3,ll), stat=ierr)
             if(ierr /= 0) &
               call terminate("writeCGNSGridFrame", &
                              "Memory allocation failure for &
                              &donorData")

             ! Determine the step for the three directions of
             ! donorData and fill the array.

             s1 = 1; s2 = 1; s3 = 1
             if(donorRange(1,2) < donorRange(1,1)) s1 = -1
             if(donorRange(2,2) < donorRange(2,1)) s2 = -1
             if(donorRange(3,2) < donorRange(3,1)) s3 = -1

             ll = 0
             do k=donorRange(3,1),donorRange(3,2),s3
               do j=donorRange(2,1),donorRange(2,2),s2
                 do i=donorRange(1,1),donorRange(1,2),s1
                   ll = ll+1

                   donorData(1,ll) = i
                   donorData(2,ll) = j
                   donorData(3,ll) = k
                 enddo
               enddo
             enddo

             ! Write the general connectivity.

             ii = ll
             call cg_conn_write_f(cgnsInd, cgnsBase, cgnsZone(nn),       &
                                  cgnsDoms(nn)%conn1to1(mm)%connectName, &
                                  Vertex, Abutting1to1, PointRange, 2,   &
                                  zoneRange,                             &
                                  cgnsDoms(nn)%conn1to1(mm)%donorName,   &
                                  Structured, PointListDonor, Integer,   &
                                  ii, donorData, jj, ierr)
             if(ierr /= all_ok)                     &
               call terminate("writeCGNSGridFrame", &
                              "Something wrong when calling &
                              &cg_conn_write_f")

             ! Deallocate the memory of donorData again.

             deallocate(donorData, stat=ierr)
             if(ierr /= 0)                          &
               call terminate("writeCGNSGridFrame", &
                              "Deallocation failure for donorData")

             ! Write the periodic info. First transform the rotation
             ! center and translation vector to the original coordinates,
             ! the angles to degrees and store everything in single
             ! precision arrays.

             rotCenter  = cgnsDoms(nn)%conn1to1(mm)%rotationCenter &
                        * LRefInv
             translation = cgnsDoms(nn)%conn1to1(mm)%translation &
                        * LRefInv
             rotRate    = cgnsDoms(nn)%conn1to1(mm)%rotationAngles &
                        * 180.0_realType/pi

             call cg_conn_periodic_write_f(cgnsInd, cgnsBase,           &
                                           cgnsZone(nn), jj, rotCenter, &
                                           rotRate, translation, ierr)
             if(ierr /= all_ok)                     &
               call terminate("writeCGNSGridFrame", &
                              "Something wrong when calling &
                              &cg_conn_periodic_write_f")

             ! Write the units of the periodic rotation.

             call cg_goto_f(cgnsInd, cgnsBase, ierr,         &
                            "Zone_t", cgnsZone(nn),          &
                            "ZoneGridConnectivity_t", 1,     &
                            "GridConnectivity_t", jj,        &
                            "GridConnectivityProperty_t", 1, &
                            "Periodic_t", 1, "DataArray_t", 2, "end")
             if(ierr /= all_ok)                     &
               call terminate("writeCGNSGridFrame", &
                              "Something wrong when calling cg_goto_f")

             call cg_dataclass_write_f(Dimensional, ierr)
             if(ierr /= all_ok)                     &
               call terminate("writeCGNSGridFrame", &
                              "Something wrong when calling &
                              &cg_dataclass_write_f")

             call cg_units_write_f(Null, Null, Null, Null, &
                                   Degree, ierr)
             if(ierr /= all_ok)                     &
               call terminate("writeCGNSGridFrame", &
                              "Something wrong when calling &
                              &cg_units_write_f")

           else periodicTest

             ! Normal 1 to 1 subface. Set the elements for the
             ! abbreviation of the transformation matrix.

             transform(1) = cgnsDoms(nn)%conn1to1(mm)%l1
             transform(2) = cgnsDoms(nn)%conn1to1(mm)%l2
             transform(3) = cgnsDoms(nn)%conn1to1(mm)%l3

             ! Write the connectivity.

             call cg_1to1_write_f(cgnsInd, cgnsBase, cgnsZone(nn),       &
                                  cgnsDoms(nn)%conn1to1(mm)%connectName, &
                                  cgnsDoms(nn)%conn1to1(mm)%donorName,   &
                                  zoneRange, donorRange, transform,      &
                                  ii, ierr)
             if(ierr /= all_ok)                     &
               call terminate("writeCGNSGridFrame", &
                              "Something wrong when calling &
                              &cg_1to1_write_f")

           endif periodicTest

         enddo loop1to1

         ! Loop over the boundary subfaces and write the data.

         loopBocos: do mm=1,cgnsDoms(nn)%nBocos

           ! Check if this is an actual face. If not, continue with
           ! the next face.

           if(.not. cgnsDoms(nn)%bocoInfo(mm)%actualFace) cycle

           ! Store the range of the subface in zoneRange.

           zoneRange(1,1) = cgnsDoms(nn)%bocoInfo(mm)%iBeg
           zoneRange(2,1) = cgnsDoms(nn)%bocoInfo(mm)%jBeg
           zoneRange(3,1) = cgnsDoms(nn)%bocoInfo(mm)%kBeg

           zoneRange(1,2) = cgnsDoms(nn)%bocoInfo(mm)%iEnd
           zoneRange(2,2) = cgnsDoms(nn)%bocoInfo(mm)%jEnd
           zoneRange(3,2) = cgnsDoms(nn)%bocoInfo(mm)%kEnd

           ! Write the boundary condition. As the preprocessing
           ! overwrites the BCType for a family specified BC, the
           ! boundary condition is constructed first and stored in jj.

           jj = cgnsDoms(nn)%bocoInfo(mm)%BCTypeCGNS
           ll = cgnsDoms(nn)%bocoInfo(mm)%familyID
           if(ll > 0) jj = FamilySpecified

           call cg_boco_write_f(cgnsInd, cgnsBase, cgnsZone(nn),    &
                                cgnsDoms(nn)%bocoInfo(mm)%bocoName, &
                                jj, PointRange, 2, zoneRange, ii, ierr)
           if(ierr /= all_ok)                     &
             call terminate("writeCGNSGridFrame", &
                            "Something wrong when calling &
                            &cg_boco_write_f")

           ! Write the family name if the boundary condition is
           ! specified per family.

           if(ll > 0) then

             ! Go to the current boundary condition and write
             ! the appropriate family name.

             call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                            "Zone_t", cgnsZone(nn),  &
                            "ZoneBC_t", 1, "BC_t", ii, "end")
             if(ierr /= all_ok)                     &
               call terminate("writeCGNSGridFrame", &
                              "Something wrong when calling cg_goto_f")

             call cg_famname_write_f(cgnsFamilies(ll)%familyName, ierr)
             if(ierr /= all_ok)                     &
               call terminate("writeCGNSGridFrame", &
                              "Something wrong when calling &
                              &cg_famname_write_f")
           endif

           ! If the boundary condition is UserDefined, write the
           ! description of what type of user defined BC.

           if(jj == UserDefined) then

             ! Go to the current boundary condition and write
             ! the appropriate data.

             call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                            "Zone_t", cgnsZone(nn),  &
                            "ZoneBC_t", 1, "BC_t", ii, "end")
             if(ierr /= all_ok)                     &
               call terminate("writeCGNSGridFrame", &
                              "Something wrong when calling cg_goto_f")

             call cg_user_data_write_f(    &
                         cgnsDoms(nn)%bocoInfo(mm)%userDefinedName, ierr)
             if(ierr /= all_ok)                     &
               call terminate("writeCGNSGridFrame", &
                              "Something wrong when calling &
                              &cg_user_data_write_f")
           endif

           ! If this boundary condition has allocated memory for data
           ! sets, write them.

           if( cgnsDoms(nn)%bocoInfo(mm)%dataSetAllocated ) then

             ! Set the pointer for the data sets to make the code
             ! more readable.

             dataSet => cgnsDoms(nn)%bocoInfo(mm)%dataSet

             ! Loop over the number of data sets for this boundary face.

             do ll=1,cgnsDoms(nn)%bocoInfo(mm)%nDataSet

               ! Create the bc dataset node.

               call cg_dataset_write_f(cgnsInd, cgnsBase,       &
                                       cgnsZone(nn), ii,        &
                                       dataSet(ll)%datasetName, &
                                       dataSet(ll)%BCType, jj, ierr)
               if(ierr /= all_ok)                     &
                 call terminate("writeCGNSGridFrame", &
                                "Something wrong when calling &
                                &cg_dataset_write_f")

               ! Write the Dirichlet and Neumann boundary condition
               ! data sets if present.

               call writeBcdataArrays(dataSet(ll)%ndirichletArrays, &
                                      dataSet(ll)%dirichletArrays,  &
                                      Dirichlet)

               call writeBcdataArrays(dataSet(ll)%nneumannArrays, &
                                      dataSet(ll)%neumannArrays,  &
                                      Neumann)
             enddo
           endif

         enddo loopBocos
       enddo zoneLoop

       !=================================================================

       contains

         !===============================================================

         subroutine writeBcdataArrays(narr, arr, DirNeu)
!
!        ****************************************************************
!        *                                                              *
!        * writeBcdataArrays writes the given bc data set arrays,       *
!        * either of the dirichlet or neumann type, to the correct      *
!        * position in the CGNS file.                                   *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         integer, intent(in) :: DirNeu

         integer(kind=intType), intent(in) :: narr
         type(cgnsBcdataArray), pointer, dimension(:) :: arr
!
!        Local variables.
!
         integer :: ierr
         integer :: realTypeCGNS

         integer(kind=intType) :: i, j, kk

         real(kind=cgnsRealType), dimension(:), allocatable :: tmp
!
!        Function definition.
!
         integer :: setCGNSRealType
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution.                                             *
!        *                                                              *
!        ****************************************************************
!
         ! Return immediately if narr == 0, i.e. if there is nothing
         ! to write.

         if(narr == 0) return

         ! Set the cgns real type.

         realTypeCGNS = setCGNSRealType()

         ! Create the BCData node.

         call cg_bcdata_write_f(cgnsInd, cgnsBase, cgnsZone(nn), &
                                  ii, jj, DirNeu, ierr)
         if(ierr /= all_ok)                    &
           call terminate("writeBcdataArrays", &
                          "Something wrong when calling &
                          &cg_bcdata_write_f")

         ! Loop over the number of data arrays.

         loopDataArrays: do kk=1,narr

           ! Go to the main node of the data arrays.

           call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t",       &
                          cgnsZone(nn), "ZoneBC_t", 1, "BC_t", ii, &
                          "BCDataSet_t", jj, "BCData_t", DirNeu, "end")
           if(ierr /= all_ok)                    &
             call terminate("writeBcdataArrays", &
                            "Something wrong when calling cg_goto_f")

           ! Determine the total size of the prescribed data,
           ! allocate the memory for tmp and copy the data into it.

           j = arr(kk)%dataDim(1)
           do i=2,arr(kk)%nDimensions
             j = j*arr(kk)%dataDim(i)
           enddo

           allocate(tmp(j), stat=ierr)
           if(ierr /= 0)                         &
             call terminate("writeBcdataArrays", &
                            "Memory allocation failure for tmp")

           tmp = arr(kk)%dataArr

           ! Write the data array and release the memory of tmp
           ! afterwards.

           call cg_array_write_f(arr(kk)%arrayName,   realTypeCGNS,    &
                                 arr(kk)%nDimensions, arr(kk)%dataDim, &
                                 tmp, ierr)
           if(ierr /= all_ok)                    &
             call terminate("writeBcdataArrays", &
                            "Something wrong when calling &
                            &cg_array_write_f")

           deallocate(tmp, stat=ierr)
           if(ierr /= 0)                         &
             call terminate("writeBcdataArrays", &
                            "Deallocation failure for tmp")

           ! Write the dimensional info for this array.

           call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t",       &
                          cgnsZone(nn), "ZoneBC_t", 1, "BC_t", ii, &
                          "BCDataSet_t", jj, "BCData_t", DirNeu,   &
                          "DataArray_t", kk, "end")
           if(ierr /= all_ok)                    &
             call terminate("writeBcdataArrays", &
                            "Something wrong when calling cg_goto_f")

           call cg_dataclass_write_f(Dimensional, ierr)
           if(ierr /= all_ok)                    &
             call terminate("writeBcdataArrays", &
                            "Something wrong when calling &
                            &cg_dataclass_write_f")

           call cg_units_write_f(arr(kk)%mass,  arr(kk)%len,  &
                                 arr(kk)%time,  arr(kk)%temp, &
                                 arr(kk)%angle, ierr)
           if(ierr /= all_ok)                    &
             call terminate("writeBcdataArrays", &
                            "Something wrong when calling &
                            &cg_units_write_f")

         enddo loopDataArrays

         end subroutine writeBcdataArrays

#endif

       end subroutine writeCGNSGridFrame
