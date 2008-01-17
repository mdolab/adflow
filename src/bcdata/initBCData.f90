!
!      ******************************************************************
!      *                                                                *
!      * File:          initBCData.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-07-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initBCData
!
!      ******************************************************************
!      *                                                                *
!      * initBCData allocates and initializes the arrays BCData for     *
!      * all boundary subfaces on all grid levels for all spectral      *
!      * solutions.                                                     *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use inputTimeSpectral
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, sps
       integer(kind=intType) :: nLevels, level
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of grid levels.

       nLevels = ubound(flowDoms,2)

       ! Loop over the number of grid levels.

       levelLoop: do level=1,nLevels

         ! Loop over the number of spectral solutions and number of
         ! blocks stored on this processor.

         spectralLoop: do sps=1,nTimeIntervalsSpectral
           domainsLoop: do i=1,nDom

             ! Allocate the memory for the array of the boundary
             ! condition data.

             j = flowDoms(i,level,sps)%nBocos
             allocate(flowDoms(i,level,sps)%BCData(j), stat=ierr)
             if(ierr /= 0)                   &
               call terminate("initBCData", &
                              "Memory allocation failure for BCData")

             ! Set the pointers to make it more readable.

             call setPointers(i,level,sps)

             ! Copy the range of the subfaces in BCData and nullify its
             ! pointers.

             bocoLoop: do j=1,nBocos

               ! Determine the block face on which the subface is located
               ! and set the dimensions accordingly.

               select case (BCFaceID(j))

                 case (iMin,iMax)
                   BCData(j)%inBeg = jnBeg(j)
                   BCData(j)%inEnd = jnEnd(j)
                   BCData(j)%jnBeg = knBeg(j)
                   BCData(j)%jnEnd = knEnd(j)

                   BCData(j)%icbeg = jcbeg(j)
                   BCData(j)%icend = jcend(j)
                   BCData(j)%jcbeg = kcbeg(j)
                   BCData(j)%jcend = kcend(j)

                 case (jMin,jMax)
                   BCData(j)%inBeg = inBeg(j)
                   BCData(j)%inEnd = inEnd(j)
                   BCData(j)%jnBeg = knBeg(j)
                   BCData(j)%jnEnd = knEnd(j)

                   BCData(j)%icbeg = icbeg(j)
                   BCData(j)%icend = icend(j)
                   BCData(j)%jcbeg = kcbeg(j)
                   BCData(j)%jcend = kcend(j)

                 case (kMin,kMax)
                   BCData(j)%inBeg = inBeg(j)
                   BCData(j)%inEnd = inEnd(j)
                   BCData(j)%jnBeg = jnBeg(j)
                   BCData(j)%jnEnd = jnEnd(j)

                   BCData(j)%icbeg = icbeg(j)
                   BCData(j)%icend = icend(j)
                   BCData(j)%jcbeg = jcbeg(j)
                   BCData(j)%jcend = jcend(j)

               end select

               ! Initialize the boundary condition treatment for
               ! subsonic inlet to noSubInlet.

               BCData(j)%subsonicInletTreatment = noSubInlet

               ! Nullify the pointers of BCData.
               ! Some compilers require this.

               nullify(BCData(j)%norm)
               nullify(BCData(j)%rface)

               nullify(BCData(j)%uSlip)
               nullify(BCData(j)%TNS_Wall)

               nullify(BCData(j)%ptInlet)
               nullify(BCData(j)%ttInlet)
               nullify(BCData(j)%htInlet)
               nullify(BCData(j)%flowXdirInlet)
               nullify(BCData(j)%flowYdirInlet)
               nullify(BCData(j)%flowZdirInlet)

               nullify(BCData(j)%turbInlet)

               nullify(BCData(j)%rho)
               nullify(BCData(j)%velx)
               nullify(BCData(j)%vely)
               nullify(BCData(j)%velz)
               nullify(BCData(j)%ps)

             enddo bocoLoop
           enddo domainsLoop
         enddo spectralLoop
       enddo levelLoop

       end subroutine initBCData
