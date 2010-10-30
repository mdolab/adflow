!
!      ******************************************************************
!      *                                                                *
!      * File:          setSupersonicInletFreeStream.f90                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-13-2005                                      *
!      * Last modified: 07-11-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setSupersonicInletFreeStream
!
!      ******************************************************************
!      *                                                                *
!      * setSupersonicInletFreeStream sets for all boundary subfaces    *
!      * stored in freestreamSubfaces the primitive flow variables to   *
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

       integer(kind=intType) :: nn, i, j
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
       ! Return immediately if freestreamSubfaces has not been allocated,
       ! i.e. if no subfaces present for which the inflow data should
       ! be set to free stream values.

       if(.not. allocated(freestreamSubfaces)) return

       ! Determine the maximum number of multigrid levels present
       ! in flowDoms.

       nLevels = ubound(flowDoms,2)

       ! Loop over the number of flagged subfaces.

       subfaceLoop: do nn=1,nFreestreamSubfaces

         ! Store the local block ID, the boundary subface and the
         ! spectral solution a bit easier.

         mm   = freestreamSubfaces(nn,1)
         boco = freestreamSubfaces(nn,2)
         sps  = freestreamSubfaces(nn,3)

         ! Loop over the number multigrid levels, starting at the
         ! current groundLevel.

         do level=groundLevel,nLevels

           ! Set the pointer for BCData to make the code more readable.

           BCData => flowDoms(mm,level,sps)%BCData

           ! Loop over the range of this subface.

           do j=BCData(boco)%jcBeg,BCData(boco)%jcEnd
             do i=BCData(boco)%icBeg,BCData(boco)%icEnd

               ! Set the flow field variables.

               BCData(boco)%rho(i,j)  = wInf(irho)
               BCData(boco)%velx(i,j) = wInf(ivx)
               BCData(boco)%vely(i,j) = wInf(ivy)
               BCData(boco)%velz(i,j) = wInf(ivz)
               BCData(boco)%ps(i,j)   = pInfCorr

             enddo
           enddo
         enddo

       enddo subfaceLoop

       ! Release the memory of freestreamSubfaces.

       deallocate(freestreamSubfaces, stat=ierr)
       if(ierr /= 0)                                    &
         call terminate("setSupersonicInletFreeStream", &
                        "Deallocation failure for freestreamSubfaces")

       end subroutine setSupersonicInletFreeStream
