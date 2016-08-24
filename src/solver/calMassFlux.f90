!
!      ******************************************************************
!      *                                                                *
!      * File:          calMassFlux.f90                                 *
!      * Author:        Eran Arad                                       *
!      * Starting date: 30-01-2006                                      *
!      * Last modified: 25-01-2007                                               *
!      *                                                                *
!      ******************************************************************
!
       subroutine calMassFlux(massFluxG)
!
!      ******************************************************************
!      *                                                                *
!      * Calculate mass flux on        subsonic outflow boundary        *
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
       use utils, only : setPointers
       implicit none
!
!      Subroutine arguments
!
       real(kind=realType) :: massFluxG

!      Local variables.
!
       integer  :: ierr
       integer(kind=intType) :: nn, mm, iBeg, iEnd, jBeg, jEnd, j, i
       real(kind=realType) :: massFluxL, sF, vn1, vn2, fluxSubface, fact

       real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2, ss
       real(kind=realType), dimension(:,:),   pointer :: sFace

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if no 

       if(nOutflowBleeds + nOutflowSubsonic + nInflowSubsonic  == 0) return

 ! Initialize the local contribution to the mass fluxes to 0.

       massFluxL = zero

      ! Compute the local contribution to the mass flux

       domains: do nn=1,nDom

         ! Set the pointers for this block to make the code
         ! more readable.

          call setPointers(nn, currentLevel,1)


       ! Loop over the boundary condition subfaces of this block.

          bocos: do mm=1,nBocos

         ! Check for  outfloww/inflow  boundary condition.

             outflowSubsonic: if(BCType(mm) == SubsonicOutflow  .or. &
                  BCType(mm) == MassBleedOutflow .or. BCType(mm) == SubsonicInflow) then

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

             ! Initialize sF to zero. This value will be used if the
             ! block is not moving.

                sF = zero

             ! Loop over the internal faces, i.e. without halo's, of the
             ! subface, and compute the contribution to mass flux.
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

                massFluxL = massFluxL + fluxSubface*fact


             endif outflowSubsonic

          enddo bocos
       enddo domains


       ! Perform an allreduce to obtain the values of the global mass
       ! fluxes for the bleed region as well as the global mass flux


       massFluxG = zero
       call mpi_reduce(massFluxL, massFluxG, 1, sumb_real, &
            mpi_sum, 0, SUmb_comm_world, ierr)


      ! Compute the dimensional mass flow

       massFluxG =  sqrt(pRef*rhoRef)*massFluxG ! rho*ur*md, ur=sqrt(p/rho)


       end subroutine calMassFlux
