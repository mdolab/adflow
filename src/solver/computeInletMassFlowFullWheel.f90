!
!      ******************************************************************
!      *                                                                *
!      * File:          computeInletMassFlowFullWheel.f90               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-27-2005                                      *
!      * Last modified: 11-20-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine computeInletMassFlowFullWheel(sps, massFlowWheel, &
                                                computeGlobalValue)
!
!      ******************************************************************
!      *                                                                *
!      * computeInletMassFlowFullWheel computes the mass flow at the    *
!      * inlet for the given spectral solution. If only a slice is      *
!      * simulated the mass flow is scaled up to the full wheel, such   *
!      * that no discrepencies occur when sections are present with     *
!      * different periodic angles and the mass flow is used to compute *
!      * other quantities.                                              *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use communication
       use iteration
       use section
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)  :: sps
       real(kind=realType),   intent(out) :: massFlowWheel
       logical,               intent(in)  :: computeGlobalValue
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm, i, j
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
       integer(kind=intType) :: nSlices

       real(kind=realType) :: localMassFlow, fact
       real(kind=realType) :: rho, ux, uy, uz, un

       real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2, ss
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Compute the local contribution to the mass flow rate.

       localMassFlow = zero
       nSlices       = 0

       domainLoop: do nn=1,nDom

         ! Set the pointers for this block to make the code
         ! more readable.

         call setPointers(nn, currentLevel,sps)

         ! Loop over the local number of boundary subfaces and test for
         ! a subsonic inflow boundary.

         subfaceLoop: do mm=1,nBocos

           testInflow: if((BCType(mm) == SubsonicInflow)       .or. &
                          (BCType(mm) == DomainInterfaceTotal) .or. &
                          (BCType(mm) == DomainInterfaceAll)) then

             ! Subface is an inflow boundary. Set a couple of pointers
             ! depending on the block face to make the code more
             ! readable; fact is present, such that the dot product
             ! with the inward normal is taken to compute the
             ! mass flow.

             select case( BCFaceID(mm) )

               case(iMin)
                 ww1 => w(1,1:,1:,:)
                 ww2 => w(2,1:,1:,:)
                 ss  => sI(1,:,:,:)
                 fact = one

               !=========================================================

               case(iMax)
                 ww1 => w(ie,1:,1:,:)
                 ww2 => w(il,1:,1:,:)
                 ss  => sI(il,:,:,:)
                 fact = -one

               !=========================================================

               case(jMin)
                 ww1 => w(1:,1,1:,:)
                 ww2 => w(1:,2,1:,:)
                 ss  => sJ(:,1,:,:)
                 fact = one

               !=========================================================

               case(jMax)
                 ww1 => w(1:,je,1:,:)
                 ww2 => w(1:,jl,1:,:)
                 ss  => sJ(:,jl,:,:)
                 fact = -one

               !=========================================================

               case(kMin)
                 ww1 => w(1:,1:,1,:)
                 ww2 => w(1:,1:,2,:)
                 ss  => sK(:,:,1,:)
                 fact = one

               !=========================================================

               case(kMax)
                 ww1 => w(1:,1:,ke,:)
                 ww2 => w(1:,1:,kl,:)
                 ss  => sK(:,:,kl,:)
                 fact = -one

             end select

             ! Loop over the internal faces, i.e. without halo's,
             ! of the subface. As the cell range may contain the halo's,
             ! it is more convenient to use the nodal range. This range
             ! is guaranteed to contain the internals only. As there is
             ! one more node than face an offset of +1 is added to inBeg
             ! and jnBeg.

             jBeg = BCData(mm)%jnBeg + 1
             jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg + 1
             iEnd = BCData(mm)%inEnd

             do j=jBeg,jEnd
               do i=iBeg,iEnd

                 ! Compute the average quantities at the interface
                 ! and update the local massflow.

                 rho = half*(ww1(i,j,irho) + ww2(i,j,irho))
                 ux  = half*(ww1(i,j,ivx)  + ww2(i,j,ivx))
                 uy  = half*(ww1(i,j,ivy)  + ww2(i,j,ivy))
                 uz  = half*(ww1(i,j,ivz)  + ww2(i,j,ivz))

                 un = ux*ss(i,j,1) + uy*ss(i,j,2) + uz*ss(i,j,3)
                 localMassFlow = localMassFlow + fact*rho*un

               enddo
             enddo

             ! Determine the number of slices for this section.
             ! Throw an error if a different periodicity is encountered.

             j = sections(sectionID)%nSlices

             if(nSlices == 0) then
               nSlices = j
             else if(nSlices /= j) then
               call terminate("computeInletMassFlowFullWheelnSlices", &
                              "Different periodicity encountered for &
                              &inlet mass flow.")
             endif

           endif testInflow

         enddo subfaceLoop

       enddo domainLoop

       ! Multiply the local mass flow by the number of slices to
       ! obtain the local mass flow proportional to the entire wheel.

       localMassFlow = localMassFlow*nSlices

       ! Do an allreduce to obtain the total mass flow if the global
       ! value is to be computed. Otherwise simply copy the local
       ! value and the allreduce is performed by the calling routine.

       if( computeGlobalValue ) then
         call mpi_allreduce(localMassFlow, massFlowWheel, 1, sumb_real, &
                            mpi_sum, SUmb_comm_world, ierr)
       else
         massFlowWheel = localMassFlow
       endif

       end subroutine computeInletMassFlowFullWheel
