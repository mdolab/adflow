!
!      ******************************************************************
!      *                                                                *
!      * File:          bleedFlowParameters.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-15-2005                                      *
!      * Last modified: 10-31-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bleedFlowParameters(sps, useInternalOnly)
!
!      ******************************************************************
!      *                                                                *
!      * The routine bleedFlowParameters determines the global          *
!      * parameters for the given spectral solution on the active MG    *
!      * level such that the bleed flow boundary conditions can be      *
!      * applied. If no bleed flow regions are present an immediate     *
!      * return is made.                                                *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use bleedFlows
       use blockPointers
       use communication
       use flowVarRefState
       use iteration
       use section
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: sps
       logical,               intent(in) :: useInternalOnly
!
!      Local variables.
!
       integer :: nSize, ierr

       integer(kind=intType) :: nBleeds, nSlices, nn, mm, i, j
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

       real(kind=realType) :: fluxSubface, fact, sF, vn1, vn2
       real(kind=realType), dimension(nInflowBleeds+nOutflowBleeds) :: &
                                            massFluxLocal, massFluxGlobal

       real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2, ss
       real(kind=realType), dimension(:,:),   pointer :: sFace
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if no bleed flow regions are present.

       if(nInflowBleeds == 0 .and. nOutflowBleeds == 0) return

       ! Initialize the local contribution to the mass fluxes to 0.

       nBleeds = nInflowBleeds + nOutflowBleeds
       do nn=1,nBleeds
         massFluxLocal(nn) = zero
       enddo

       ! Compute the local contribution to the mass fluxes of both
       ! the inflow and the outflow bleeds.

       domainLoop: do nn=1,nDom

         ! Set the pointers for this block to make the code
         ! more readable.

         call setPointers(nn, currentLevel,sps)

         ! Loop over the boundary subfaces.

         subfaceLoop: do mm=1,nBocos

           ! Check for a bleed boundary condition.

           testBleed: if(BCType(mm) == MassBleedInflow .or. &
                         BCType(mm) == MassBleedOutflow) then

             ! Set the pointers for the normal vector and solution in the
             ! inside of the block and in the halos. This depends on the
             ! block face on which the subface is located. The variable
             ! fact is present, such that the dot product with the inward
             ! normal is computed to obtain the velocity flux. For an
             ! outflow this will be reversed later on.

             select case( BCFaceID(mm) )
               case (iMin)
                 ss  => sI(1,:,:,:);  fact = one
                 ww2 => w(2,1:,1:,:); ww1 => w(1,1:,1:,:)
                 if( addGridVelocities ) sFace => sFaceI(1,:,:)

               !=========================================================

               case (iMax)
                 ss  => sI(il,:,:,:);  fact = -one
                 ww2 => w(il,1:,1:,:); ww1 => w(ie,1:,1:,:)
                 if( addGridVelocities ) sFace => sFaceI(il,:,:)

               !=========================================================

               case (jMin)
                 ss  => sJ(:,1,:,:);  fact = one
                 ww2 => w(1:,2,1:,:); ww1 => w(1:,1,1:,:)
                 if( addGridVelocities ) sFace => sFaceJ(:,1,:)

               !=========================================================

               case (jMax)
                 ss  => sJ(:,jl,:,:);  fact = -one
                 ww2 => w(1:,jl,1:,:); ww1 => w(1:,je,1:,:)
                 if( addGridVelocities ) sFace => sFaceJ(:,jl,:)

               !=========================================================

               case (kMin)
                 ss  => sK(:,:,1,:);  fact = one
                 ww2 => w(1:,1:,2,:); ww1 => w(1:,1:,1,:)
                 if( addGridVelocities ) sFace => sFaceK(:,:,1)

               !=========================================================

               case (kMax)
                 ss  => sK(:,:,kl,:);  fact = -one
                 ww2 => w(1:,1:,kl,:); ww1 => w(1:,1:,ke,:)
                 if( addGridVelocities ) sFace => sFaceK(:,:,kl)
             end select

             ! Reset the pointer for ww1 if only internal values must
             ! be used.

             if( useInternalOnly ) ww1 => ww2

             ! Reverse fact for an outflow boundary, such that the
             ! mass flow is positive when leaving the domain.

             if(BCType(mm) == MassBleedOutflow) fact = -fact

             ! Initialize sF to zero. This value will be used if the
             ! block is not moving.

             sF = zero

             ! Loop over the internal faces, i.e. without halo's, of the
             ! subface, and compute the contribution to velocity flux.
             ! As the cell range may contain the halo's, it is more
             ! convenient to use the nodal range. This range is
             ! guaranteed to contain the internals only. As there is
             ! one more node than face an offset of +1 is added to inBeg
             ! and jnBeg.

             jBeg = BCData(mm)%jnBeg + 1
             jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg + 1
             iEnd = BCData(mm)%inEnd

             fluxSubface = zero

             do j=jBeg,jEnd
               do i=iBeg,iEnd

                 ! Set the dot product of the grid velocity and the
                 ! normal for a moving face.

                 if( addGridVelocities ) sF = sFace(i,j)

                 ! Compute the normal velocity to the left and right
                 ! of the face.

                 vn1 = ww1(i,j,ivx)*ss(i,j,1) + ww1(i,j,ivy)*ss(i,j,2) &
                     + ww1(i,j,ivz)*ss(i,j,3) - sF
                 vn2 = ww2(i,j,ivx)*ss(i,j,1) + ww2(i,j,ivy)*ss(i,j,2) &
                     + ww2(i,j,ivz)*ss(i,j,3) - sF

                 ! Update the value of fluxSubface by computing
                 ! the massflow through the subface.

                 fluxSubface = fluxSubface             &
                             + half*(ww1(i,j,irho)*vn1 &
                             +       ww2(i,j,irho)*vn2)
               enddo
             enddo

             ! Determine the entry in massFluxLocal (the outflow regions
             ! are stored after the inflow regions) as well as the number
             ! of slices of the blade row for this subface. Update the
             ! corresponding value of massFluxLocal. Multiply by nSlices
             ! to obtain the value of the entire wheel.

             nSlices = sections(sectionID)%nSlices
             i       = groupNum(mm)
             if(BCType(mm) == MassBleedOutflow) i = i + nInflowBleeds

             massFluxLocal(i) = massFluxLocal(i) &
                              + nSlices*fluxSubface*fact

           endif testBleed

         enddo subfaceLoop
       enddo domainLoop

       ! Perform an allreduce to obtain the values of the global mass
       ! fluxes for the bleed regions.

       nSize = nBleeds
       call mpi_allreduce(massFluxLocal, massFluxGlobal, nSize, &
                          sumb_real, mpi_sum, SUmb_comm_world, ierr)

       ! The currently stored value is the nondimensional mass flow.
       ! As the prescribed value is dimensional, the dimensional value
       ! is needed. Therefore a multiplication by the reference values
       ! is needed.

       fact = sqrt(pRef*rhoRef)

       do nn=1,nInflowBleeds
         inflowBleeds(nn)%curMassFlux = fact*massFluxGlobal(nn)
       enddo

       do nn=1,nOutflowBleeds
         outflowBleeds(nn)%curMassFlux = fact &
                                       * massFluxGlobal(nn+nInflowBleeds)
       enddo

       end subroutine bleedFlowParameters
