!
!      ******************************************************************
!      *                                                                *
!      * File:          setInletFreestreamTurb.f90                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-29-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setInletFreestreamTurb
!
!      ******************************************************************
!      *                                                                *
!      * setInletFreestreamTurb sets for all boundary subfaces          *
!      * stored in turbFreestreamSubfaces the turbulence variables to   *
!      * the free stream variables. This is done for all multigrid      *
!      * levels starting from groundLevel.                              *
!      *                                                                *
!      ******************************************************************
!
       use block
       use flowVarRefState
       use iteration
       use BCDataMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, i, j, l
       integer(kind=intType) :: mm, boco, sps
       integer(kind=intType) :: level, nLevels

       type(BCDataType), dimension(:), pointer :: BCData
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! If turbFreestreamSubfaces has not been allocated there are
       ! no such subfaces present and an immediate return should be made.

       if(.not. allocated(turbFreestreamSubfaces) ) return

       ! Determine the maximum number of multigrid levels present
       ! in flowDoms.

       nLevels = ubound(flowDoms,2)

       ! Loop over the number of flagged subfaces.

       subfaceLoop: do nn=1,nTurbFreestreamSubfaces

         ! Store the local block ID, the boundary subface and the
         ! spectral solution a bit easier.

         mm   = turbFreestreamSubfaces(nn,1)
         boco = turbFreestreamSubfaces(nn,2)
         sps  = turbFreestreamSubfaces(nn,3)

         ! Loop over the number multigrid levels, starting at the
         ! current groundLevel.

         do level=groundLevel,nLevels

           ! Set the pointer for BCData to make the code more readable.

           BCData => flowDoms(mm,level,sps)%BCData

           ! Loop over the range of this subface.

           do j=BCData(boco)%jcBeg,BCData(boco)%jcEnd
             do i=BCData(boco)%icBeg,BCData(boco)%icEnd

               ! Set the turbulence variables.

               do l=nt1,nt2
                 BCData(boco)%turbInlet(i,j,l) = wInf(l)
               enddo
             enddo
           enddo
         enddo

       enddo subfaceLoop

       ! Release the memory of turbFreestreamSubfaces.

       deallocate(turbFreestreamSubfaces, stat=ierr)
       if(ierr /= 0) &
         call terminate("setInletFreestreamTurb", &
                        "Deallocation failure for &
                        &turbFreestreamSubfaces")

       end subroutine setInletFreestreamTurb
