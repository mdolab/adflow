!
!      ******************************************************************
!      *                                                                *
!      * File:          setFamilyInfoFaces.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-18-2007                                      *
!      * Last modified: 10-30-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setFamilyInfoFaces(level)
!
!      ******************************************************************
!      *                                                                *
!      * setFamilyInfoFaces sets the values of the family parameters    *
!      * for faces on the given multigrid level. The default values for *
!      * indFamily is 0, which means that the mass flow through that    *
!      * face does not contribute to the mass flow that must be         *
!      * monitored. For sliding mesh interfaces both sides of the       *
!      * interface are monitored and the value of indFamily corresponds *
!      * to one of the two entries in the local monitoring arrays. The  *
!      * values of factFamily are such that the mass flow entering the  *
!      * block is defined positive.                                     *
!      * Note that only the 1st spectral solution is treated, because   *
!      * this informations is the same for all of them.                 *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use inputTimeSpectral
       use monitor
       use section
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm, i, j, k, ii, nSlices
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

       integer(kind=intType), dimension(cgnsNFamilies) :: orToMassFam
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the offset ii for the values of orToMassFam. If the mass
       ! flow through sliding mesh interfaces must be monitored this
       ! offset if 2*cgnsNSliding. This means that in the arrays to store
       ! the mass flow the sliding mesh interfaces are stored first,
       ! followed by the families.

       if( monMassSliding ) then
         ii = 2*cgnsNSliding
       else
         ii = 0
       endif

       ! Determine the number of families for which the mass flow must
       ! be monitored and set the entries of orToMassFam accordingly,
       ! i.e. the offset ii is included.

       mm = ii
       do nn=1,cgnsNFamilies
         if(cgnsFamilies(nn)%monitorMassflow            .and. &
            cgnsFamilies(nn)%BCType /= MassBleedInflow  .and. &
            cgnsFamilies(nn)%BCType /= MassBleedOutflow .and. &
            cgnsFamilies(nn)%BCType /= SlidingInterface) then
           mm = mm + 1
           orToMassFam(nn) = mm
         else
           orToMassFam(nn) = 0
         endif
       enddo

       ! Set monMassFamilies to .true. if the mass flow of at least one
       ! family must be monitored. Otherwise set it to .false.

       if(mm > ii) then
         monMassFamilies = .true.
       else
         monMassFamilies = .false.
       endif

       ! If this is the first level, allocate the memory for
       ! massFlowFamilyInv and massFlowFamilyDiss.

       if(level == 1) then
         nn = nTimeIntervalsSpectral
         allocate(massFlowFamilyInv(0:mm,nn), &
                  massFlowFamilyDiss(0:mm,nn), stat=ierr)
         if(ierr /= 0) &
           call terminate("setFamilyInfoFaces", &
                          "Memory allocation failure for &
                           &massFlowFamilyInv and massFlowFamilyDiss")
       endif

       ! Loop over the number of domains.

       domains: do nn=1,nDom

         ! Allocate the memory for indFamily and factFamily.

         il = flowDoms(nn,level,1)%il
         jl = flowDoms(nn,level,1)%jl
         kl = flowDoms(nn,level,1)%kl

         allocate(flowDoms(nn,level,1)%indFamilyI (1:il,2:jl,2:kl), &
                  flowDoms(nn,level,1)%indFamilyJ (2:il,1:jl,2:kl), &
                  flowDoms(nn,level,1)%indFamilyK (2:il,2:jl,1:kl), &
                  flowDoms(nn,level,1)%factFamilyI(1:il,2:jl,2:kl), &
                  flowDoms(nn,level,1)%factFamilyJ(2:il,1:jl,2:kl), &
                  flowDoms(nn,level,1)%factFamilyK(2:il,2:jl,1:kl), &
                  stat=ierr)
         if(ierr /= 0)                          &
           call terminate("setFamilyInfoFaces", &
                          "Memory allocation failure for indFamily &
                          &and factFamily")

         ! Set the pointers for this domain.

         call setPointers(nn, level, 1_intType)

         ! Determine the number of slices for this block to make
         ! the full wheel.

         nSlices = sections(sectionID)%nSlices

         ! Initialize the values of indFamily and factFamily.

         indFamilyI = 0_intType
         indFamilyJ = 0_intType
         indFamilyK = 0_intType

         factFamilyI = 0_intType
         factFamilyJ = 0_intType
         factFamilyK = 0_intType

         ! Loop over the boundary conditions.

         boco: do mm=1,nBocos

           ! Test for the boundary condition.

           select case (BCType(mm))
             case (SlidingInterface)

               ! Sliding mesh boundary.
               ! If the mass flow through sliding interfaces must be monitored,
               ! determine the index in the arrays massFlowFamilyInv and
               ! massFlowFamilyDiss where to store the contribution of this
               ! subface. If the sliding mesh mass flows are not monitored,
               ! set the index to 0.

               if( monMassSliding ) then
                 ii = 2*abs(groupNum(mm))
                 if(groupNum(mm) < 0) ii = ii - 1
               else
                 ii = 0
               endif

             case (MassBleedInflow, MassBleedOutflow)

               ! Inflow or outflow bleed. These boundary conditions are
               ! handled separetely and need not be monitored.

               ii = 0

             case default

               ! Subface is an ordinary boundary condition. Determine the
               ! family ID and set the index ii in the arrays massFlowFamilyInv
               ! and massFlowFamilyDiss accordingly.

               if(groupNum(mm) > 0) then
                 ii = orToMassFam(groupNum(mm))
               else
                 ii = 0
               endif

           end select

           ! Set the owned cell range for the faces on this subface.
           ! As icBeg, etc. may contain halo cells, inBeg, etc. is
           ! used.

           iBeg = min(inBeg(mm), inEnd(mm)) +1
           iEnd = max(inBeg(mm), inEnd(mm))

           jBeg = min(jnBeg(mm), jnEnd(mm)) +1
           jEnd = max(jnBeg(mm), jnEnd(mm))

           kBeg = min(knBeg(mm), knEnd(mm)) +1
           kEnd = max(knBeg(mm), knEnd(mm))

           ! Determine the block this subface is located on and set
           ! the corresponding values of indFamily and factFamily.
           ! Note that factFamily is set to nSlices on min faces and to
           ! -nSlices on max faces, such that the mass flow entering the
           ! domain is defined positive and the mass flow of the entire
           ! wheel is monitored.

           select case( BCFaceID(mm) )
             case (iMin)
               do k=kBeg,kEnd
                 do j=jBeg,jEnd
                   indFamilyI (1,j,k) = ii
                   factFamilyI(1,j,k) = nSlices
                 enddo
               enddo

             !===========================================================

             case (iMax)
               do k=kBeg,kEnd
                 do j=jBeg,jEnd
                   indFamilyI (il,j,k) = ii
                   factFamilyI(il,j,k) = -nSlices
                 enddo
               enddo

             !===========================================================

             case (jMin)
               do k=kBeg,kEnd
                 do i=iBeg,iEnd
                   indFamilyJ (i,1,k) = ii
                   factFamilyJ(i,1,k) = nSlices
                 enddo
               enddo

             !===========================================================

             case (jMax)
               do k=kBeg,kEnd
                 do i=iBeg,iEnd
                   indFamilyJ (i,jl,k) = ii
                   factFamilyJ(i,jl,k) = -nSlices
                 enddo
               enddo

             !===========================================================

             case (kMin)
               do j=jBeg,jEnd
                 do i=iBeg,iEnd
                   indFamilyK (i,j,1) = ii
                   factFamilyK(i,j,1) = nSlices
                 enddo
               enddo

             !===========================================================

             case (kMax)
               do j=jBeg,jEnd
                 do i=iBeg,iEnd
                   indFamilyK (i,j,kl) = ii
                   factFamilyK(i,j,kl) = -nSlices
                 enddo
               enddo

           end select

         enddo boco
       enddo domains

       end subroutine setFamilyInfoFaces
