!
!      ******************************************************************
!      *                                                                *
!      * File:          readFamilyInfo.F90                              *
!      * Author:        Edwin van der Weide, Steve Repsher,             *
!      *                Seonghyeon Hahn                                 *
!      * Starting date: 08-21-2003                                      *
!      * Last modified: 10-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readFamilyInfo(cgnsInd, cgnsBase)
!
!      ******************************************************************
!      *                                                                *
!      * readFamilyInfo determines the number of families in the        *
!      * given base of the cgns grid and determines their possible      *
!      * boundary condition, including some user defined ones.          *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use cgnsGrid
       use communication
       use iteration
       use su_cgns
       implicit none
!
!      Subroutine arguments
!
       integer, intent(in) :: cgnsInd, cgnsBase
!
!      Local variables.
!
       integer :: nn, bc, nFamBC, nGeo, nUserData, ierr

       character(len=maxStringLen) :: errorMessage
!
!      Function definition.
!
       integer(kind=intType) :: internalBC
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS

       call terminate("readFamilyInfo", &
                      "Routine should not be called if no cgns support &
                      &is selected.")

#else
       ! Determine the number of families in the given base.

       call cg_nfamilies_f(cgnsInd, cgnsBase, nn, ierr)
       if(ierr /= all_ok)                 &
         call terminate("readFamilyInfo", &
                        "Something wrong when calling cg_nfamilies_f")
       cgnsNFamilies = nn

       ! Allocate the memory for cgnsFamilies.

       allocate(cgnsFamilies(nn), stat=ierr)
       if(ierr /= 0)                      &
         call terminate("readFamilyInfo", &
                        "Memory allocation failure for cgnsFamilies")

       ! Loop over the number of families and read the info.

       nFam: do nn=1,cgnsNFamilies

         ! Initialize slidingID to 0 to indicate that this family does
         ! not belong to a sliding mesh interface. Idem for the
         ! bleedRegionID.

         cgnsFamilies(nn)%slidingID   = 0
         cgnsFamilies(nn)%bleedRegionID = 0

         ! Initialize the logical to monitor the mass flow to .false.

         cgnsFamilies(nn)%monitorMassFlow = .false.

         ! Nullify the pointer for the prescribed boundary data.

         nullify(cgnsFamilies(nn)%dataSet)

         ! Read the family name and the number of boundary conditions
         ! specified.

         call cg_family_read_f(cgnsInd, cgnsBase, nn,       &
                               cgnsFamilies(nn)%familyName, &
                               nFamBC, nGeo, ierr)
         if(ierr /= all_ok)               &
         call terminate("readFamilyInfo", &
                        "Something wrong when calling cg_family_read_f")

         ! Determine the boundary condition for this family, if specified.

         select case (nFamBC)

           case (0)
             cgnsFamilies(nn)%BCTypeCGNS = Null
             cgnsFamilies(nn)%BCType     = BCNull
             cgnsFamilies(nn)%bcName     = ""

           !=============================================================

           case (1)
             bc = 1
             call cg_fambc_read_f(cgnsInd, cgnsBase, nn, bc, &
                                  cgnsFamilies(nn)%bcName,   &
                                  cgnsFamilies(nn)%BCTypeCGNS, ierr)
             if(ierr /= all_ok)                 &
               call terminate("readFamilyInfo", &
                              "Something wrong when calling &
                              &cg_fambc_read_f")

             ! If this is a user defined boundary condition it must
             ! contain more information to determine the internally
             ! used BC.

             testUserDefined: if(cgnsFamilies(nn)%BCTypeCGNS == &
                                 UserDefined) then

               ! Move to the family and determine the number of
               ! user defined data nodes.

               call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                              "Family_t", nn, "end")
               if(ierr /= all_ok)                 &
                 call terminate("readFamilyInfo", &
                                "Something wrong when calling cg_goto_f")

               call cg_nuser_data_f(nUserData, ierr)
               if(ierr /= all_ok)                 &
                 call terminate("readFamilyInfo", &
                                "Something wrong when calling &
                                &cg_nuser_data_f")

               ! nUserData should be 1. Check this.

               if(nUserData /= 1) then
                 write(errorMessage,101) trim(cgnsFamilies(nn)%familyName)
                 if(myID == 0) &
                   call terminate("readFamilyInfo", errorMessage)
                 call mpi_barrier(SUmb_comm_world, ierr)
               endif

               ! Read the name of the user defined data node.

               call cg_user_data_read_f(nUserData,                        &
                                        cgnsFamilies(nn)%userDefinedName, &
                                        ierr)
               if(ierr /= all_ok)                 &
                 call terminate("readFamilyInfo", &
                                "Something wrong when calling &
                                &cg_user_data_read_f")

             else testUserDefined

               ! Set the user defined name to an empty string.

               cgnsFamilies(nn)%userDefinedName = ""

             endif testUserDefined

             ! Determine the internal BC type from the CGNS type and
             ! possibly the user defined name.

             cgnsFamilies(nn)%BCType = &
                            internalBC(cgnsFamilies(nn)%BCTypeCGNS, &
                                       cgnsFamilies(nn)%userDefinedName)

           !=============================================================

           case default
             write(errorMessage,201) trim(cgnsFamilies(nn)%familyName)
             if(myID == 0) &
               call terminate("readFamilyInfo", errorMessage)
             call mpi_barrier(SUmb_comm_world, ierr)

         end select

       enddo nFam

       ! Format statements.

 101   format("Family",1x,a,": Need 1 UserDefinedData_t node for &
              &user defined boundary condition")
 102   format("Family",1x,a,": Unknown user defined boundary &
              &condition",1x,a)
 201   format("Family",1x,a,": More than 1 boundary condition specified")

#endif

       end subroutine readFamilyInfo
