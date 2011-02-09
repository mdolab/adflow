!
!      ******************************************************************
!      *                                                                *
!      * File:          initres.f90                                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-18-2003                                      *
!      * Last modified: 06-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initres(varStart, varEnd)
!
!      ******************************************************************
!      *                                                                *
!      * initres initializes the given range of the residual. Either to *
!      * zero, steady computation, or to an unsteady term for the time  *
!      * spectral and unsteady modes. For the coarser grid levels the   *
!      * residual forcing term is taken into account.                   *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use flowVarRefState
       use inputIteration
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: varStart, varEnd
!
!      Local variables.
!
       integer(kind=intType) :: sps, nn, mm, ll, ii, jj, i, j, k, l, m
       real(kind=realType)   :: oneOverDt, tmp

       real(kind=realType), dimension(:,:,:,:), pointer :: ww, wsp, wsp1
       real(kind=realType), dimension(:,:,:),   pointer :: volsp
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately of no variables are in the range.

       if(varEnd < varStart) return

       ! Loop over the number of spectral solutions and domains.

       spectralLoop: do sps=1,nTimeIntervalsSpectral

         domains: do nn=1,nDom

           ! Set the pointers to this block.

           call setPointers(nn, currentLevel, sps)

           ! Determine the equation mode and act accordingly.

           select case (equationMode)
             case (steady)

               ! Steady state computation.
               ! Determine the currently active multigrid level.

               steadyLevelTest: if(currentLevel == groundLevel) then

                 ! Ground level of the multigrid cycle. Initialize the
                 ! owned residuals to zero.

                 do l=varStart,varEnd
                   do k=2,kl
                     do j=2,jl
                       do i=2,il
                         dw(i,j,k,l) = zero
                       enddo
                     enddo
                   enddo
                 enddo

               else steadyLevelTest

                 ! Coarse grid level. Initialize the owned cells to the
                 ! residual forcing terms.

                 do l=varStart,varEnd
                   do k=2,kl
                     do j=2,jl
                       do i=2,il
                         dw(i,j,k,l) = wr(i,j,k,l)
                       enddo
                     enddo
                   enddo
                 enddo

               endif steadyLevelTest

             !===========================================================

             case (unsteady)

               ! Unsteady computation.
               ! A further distinction must be made.

               select case(timeIntegrationScheme)

                 case (explicitRK)

                   ! We are always on the finest grid.
                   ! Initialize the residual to zero.

                   do l=varStart,varEnd
                     do k=2,kl
                       do j=2,jl
                         do i=2,il
                           dw(i,j,k,l) = zero
                         enddo
                       enddo
                     enddo
                   enddo

                 !=======================================================

                 case (implicitRK)
                   call terminate("initRes", &
                                  "Implicit RK not implemented yet")

                 !=======================================================

                 case (BDF)

                   ! Store the inverse of the physical nonDimensional
                   ! time step a bit easier.

                   oneOverDt = timeRef/deltaT

                   ! Store the pointer for the variable to be used to compute
                   ! the unsteady source term. For a runge-kutta smoother this
                   ! is the solution of the zeroth runge-kutta stage. As for
                   ! rkStage == 0 this variable is not yet set w is used.
                   ! For other smoothers w is to be used as well.

                   if(smoother == RungeKutta .and. rkStage > 0) then
                     ww => wn
                   else
                     ww => w
                   endif

                   ! Determine the currently active multigrid level.

                   unsteadyLevelTest: if(currentLevel == groundLevel) then

                     ! Ground level of the multigrid cycle. Initialize the
                     ! owned cells to the unsteady source term. First the
                     ! term for the current time level. Note that in w the
                     ! velocities are stored and not the momentum variables.
                     ! Therefore the if-statement is present to correct this.

                     do l=varStart,varEnd

                       if(l == ivx .or. l == ivy .or. l == ivz) then

                         ! Momentum variables.

                         do k=2,kl
                           do j=2,jl
                             do i=2,il
                               dw(i,j,k,l) = coefTime(0)*vol(i,j,k) &
                                           * ww(i,j,k,l)*ww(i,j,k,irho)
                             enddo
                           enddo
                         enddo

                       else

                         ! Non-momentum variables, for which the variable
                         ! to be solved is stored; for the flow equations this
                         ! is the conservative variable, for the turbulent
                         ! equations the primitive variable.

                         do k=2,kl
                           do j=2,jl
                             do i=2,il
                               dw(i,j,k,l) = coefTime(0)*vol(i,j,k) &
                                           * ww(i,j,k,l)
                             enddo
                           enddo
                         enddo

                       endif

                     enddo

                     ! The terms from the older time levels. Here the
                     ! conservative variables are stored. In case of a
                     ! deforming mesh, also the old volumes must be taken.

                     deformingTest: if( deforming_Grid ) then

                       ! Mesh is deforming and thus the volumes can change.
                       ! Use the old volumes as well.

                       do m=1,nOldLevels
                         do l=varStart,varEnd
                           do k=2,kl
                             do j=2,jl
                               do i=2,il
                                 dw(i,j,k,l) = dw(i,j,k,l)                 &
                                             + coefTime(m)*volOld(m,i,j,k) &
                                             * wOld(m,i,j,k,l)
                               enddo
                             enddo
                           enddo
                         enddo
                       enddo

                     else deformingTest

                       ! Rigid mesh. The volumes remain constant.

                       do m=1,nOldLevels
                         do l=varStart,varEnd
                           do k=2,kl
                             do j=2,jl
                               do i=2,il
                                 dw(i,j,k,l) = dw(i,j,k,l)            &
                                             + coefTime(m)*vol(i,j,k) &
                                             * wOld(m,i,j,k,l)
                               enddo
                             enddo
                           enddo
                         enddo
                       enddo

                     endif deformingTest

                     ! Multiply the time derivative by the inverse of the
                     ! time step to obtain the true time derivative.
                     ! This is done after the summation has been done, because
                     ! otherwise you run into finite accuracy problems for
                     ! very small time steps.

                     do l=varStart,varEnd
                       do k=2,kl
                         do j=2,jl
                           do i=2,il
                             dw(i,j,k,l) = oneOverDt*dw(i,j,k,l)
                          enddo
                         enddo
                       enddo
                     enddo

                   else unsteadyLevelTest

                     ! Coarse grid level. Initialize the owned cells to the
                     ! residual forcing term plus a correction for the
                     ! multigrid treatment of the time derivative term.
                     ! As the velocities are stored instead of the momentum,
                     ! these terms must be multiplied by the density.

                     tmp = oneOverDt*coefTime(0)

                     do l=varStart,varEnd

                       if(l == ivx .or. l == ivy .or. l == ivz) then

                         ! Momentum variables.

                         do k=2,kl
                           do j=2,jl
                             do i=2,il
                               dw(i,j,k,l) = tmp*vol(i,j,k)               &
                                           * (ww(i,j,k,l)*ww(i,j,k,irho)  &
                                           -  w1(i,j,k,l)*w1(i,j,k,irho))
                               dw(i,j,k,l) = dw(i,j,k,l) + wr(i,j,k,l)
                             enddo
                           enddo
                         enddo

                       else

                         ! Non-momentum variables.

                         do k=2,kl
                           do j=2,jl
                             do i=2,il
                               dw(i,j,k,l) = tmp*vol(i,j,k)             &
                                           * (ww(i,j,k,l) - w1(i,j,k,l))
                               dw(i,j,k,l) =  dw(i,j,k,l) + wr(i,j,k,l)
                             enddo
                           enddo
                         enddo

                       endif

                     enddo

                   endif unsteadyLevelTest

               end select

             !===========================================================

             case (timeSpectral)

               ! Time spectral computation. The time derivative of the
               ! current solution is given by a linear combination of
               ! all other solutions, i.e. a matrix vector product.

               ! First store the section to which this block belongs
               ! in jj.

               jj = sectionID

               ! Determine the currently active multigrid level.

               spectralLevelTest: if(currentLevel == groundLevel) then

                 ! Finest multigrid level. The residual must be
                 ! initialized to the time derivative.

                 ! Initialize it to zero.

                 do l=varStart,varEnd
                   do k=2,kl
                     do j=2,jl
                       do i=2,il
                         dw(i,j,k,l) = zero
                       enddo
                     enddo
                   enddo
                 enddo

                 ! Loop over the number of terms which contribute
                 ! to the time derivative.

                 timeLoopFine: do mm=1,nTimeIntervalsSpectral

                   ! Store the pointer for the variable to be used to
                   ! compute the unsteady source term and the volume.
                   ! Also store in ii the offset needed for vector
                   ! quantities.

                   wsp   => flowDoms(nn,currentLevel,mm)%w
                   volsp => flowDoms(nn,currentLevel,mm)%vol
                   ii    =  3*(mm-1)

                   ! Loop over the number of variables to be set.

                   varLoopFine: do l=varStart,varEnd

                     ! Test for a momentum variable.

                     if(l == ivx .or. l == ivy .or. l == ivz) then

                       ! Momentum variable. A special treatment is
                       ! needed because it is a vector and the velocities
                       ! are stored instead of the momentum. Set the
                       ! coefficient ll, which defines the row of the
                       ! matrix used later on.

                       if(l == ivx) ll = 3*sps - 2
                       if(l == ivy) ll = 3*sps - 1
                       if(l == ivz) ll = 3*sps

                       ! Loop over the owned cell centers to add the
                       ! contribution from wsp.

                       do k=2,kl
                         do j=2,jl
                           do i=2,il

                             ! Store the matrix vector product with the
                             ! velocity in tmp.

                             tmp = dvector(jj,ll,ii+1)*wsp(i,j,k,ivx) &
                                 + dvector(jj,ll,ii+2)*wsp(i,j,k,ivy) &
                                 + dvector(jj,ll,ii+3)*wsp(i,j,k,ivz)

                             ! Update the residual. Note the
                             ! multiplication with the density to obtain
                             ! the correct time derivative for the
                             ! momentum variable.

                             dw(i,j,k,l) = dw(i,j,k,l) &
                                         + tmp*volsp(i,j,k)*wsp(i,j,k,irho)

                           enddo
                         enddo
                       enddo

                     else

                       ! Scalar variable.  Loop over the owned cells to
                       ! add the contribution of wsp to the time
                       ! derivative.

                       do k=2,kl
                         do j=2,jl
                           do i=2,il
                             dw(i,j,k,l) = dw(i,j,k,l)        &
                                         + dscalar(jj,sps,mm) &
                                         * volsp(i,j,k)*wsp(i,j,k,l)
                           enddo
                         enddo
                       enddo

                     endif

                   enddo varLoopFine

                 enddo timeLoopFine

               else spectralLevelTest

                 ! Coarse grid level. Initialize the owned cells to the
                 ! residual forcing term plus a correction for the
                 ! multigrid treatment of the time derivative term.

                 ! Initialization to the residual forcing term.

                 do l=varStart,varEnd
                   do k=2,kl
                     do j=2,jl
                       do i=2,il
                         dw(i,j,k,l) = wr(i,j,k,l)
                       enddo
                     enddo
                   enddo
                 enddo

                 ! Loop over the number of terms which contribute
                 ! to the time derivative.

                 timeLoopCoarse: do mm=1,nTimeIntervalsSpectral

                   ! Store the pointer for the variable to be used to
                   ! compute the unsteady source term and the pointers
                   ! for wsp1, the solution when entering this MG level
                   ! and for the volume.
                   ! Furthermore store in ii the offset needed for
                   ! vector quantities.

                   wsp   => flowDoms(nn,currentLevel,mm)%w
                   wsp1  => flowDoms(nn,currentLevel,mm)%w1
                   volsp => flowDoms(nn,currentLevel,mm)%vol
                   ii    =  3*(mm-1)

                   ! Loop over the number of variables to be set.

                   varLoopCoarse: do l=varStart,varEnd

                     ! Test for a momentum variable.

                     if(l == ivx .or. l == ivy .or. l == ivz) then

                       ! Momentum variable. A special treatment is
                       ! needed because it is a vector and the velocities
                       ! are stored instead of the momentum. Set the
                       ! coefficient ll, which defines the row of the
                       ! matrix used later on.

                       if(l == ivx) ll = 3*sps - 2
                       if(l == ivy) ll = 3*sps - 1
                       if(l == ivz) ll = 3*sps

                       ! Add the contribution of wps to the correction
                       ! of the time derivative. The difference between
                       ! the current time derivative and the one when
                       ! entering this grid level must be added, because
                       ! the residual forcing term only takes the spatial
                       ! part of the coarse grid residual into account.

                       do k=2,kl
                         do j=2,jl
                           do i=2,il

                             ! Store the matrix vector product with the
                             ! momentum in tmp.

                             tmp = dvector(jj,ll,ii+1)                &
                                 * (wsp( i,j,k,irho)*wsp( i,j,k,ivx)  &
                                 -  wsp1(i,j,k,irho)*wsp1(i,j,k,ivx)) &
                                 +  dvector(jj,ll,ii+2)               &
                                 * (wsp( i,j,k,irho)*wsp( i,j,k,ivy)  &
                                 -  wsp1(i,j,k,irho)*wsp1(i,j,k,ivy)) &
                                 + dvector(jj,ll,ii+3)                &
                                 * (wsp( i,j,k,irho)*wsp( i,j,k,ivz)  &
                                 -  wsp1(i,j,k,irho)*wsp1(i,j,k,ivz))

                             ! Add tmp to the residual. Multiply it by
                             ! the volume to obtain the finite volume
                             ! formulation of the  derivative of the
                             ! momentum.

                             dw(i,j,k,l) = dw(i,j,k,l) + tmp*volsp(i,j,k)

                           enddo
                         enddo
                       enddo

                     else

                       ! Scalar variable. Loop over the owned cells
                       ! to add the contribution of wsp to the correction
                       ! of the time derivative.

                       do k=2,kl
                         do j=2,jl
                           do i=2,il
                             dw(i,j,k,l) = dw(i,j,k,l)        &
                                         + dscalar(jj,sps,mm) &
                                         * volsp(i,j,k)       &
                                         * (wsp(i,j,k,l) - wsp1(i,j,k,l))
                           enddo
                         enddo
                       enddo

                     endif

                   enddo varLoopCoarse

                 enddo timeLoopCoarse

               endif spectralLevelTest

           end select

           ! Set the residual in the halo cells to zero. This is just
           ! to avoid possible problems. Their values do not matter.

           do l=varStart,varEnd
             do k=0,kb
               do j=0,jb
                 dw(0,j,k,l)  = zero
                 dw(1,j,k,l)  = zero
                 dw(ie,j,k,l) = zero
                 dw(ib,j,k,l) = zero
               enddo
             enddo

             do k=0,kb
               do i=2,il
                 dw(i,0,k,l)  = zero
                 dw(i,1,k,l)  = zero
                 dw(i,je,k,l) = zero
                 dw(i,jb,k,l) = zero
               enddo
             enddo

             do j=2,jl
               do i=2,il
                 dw(i,j,0,l)  = zero
                 dw(i,j,1,l)  = zero
                 dw(i,j,ke,l) = zero
                 dw(i,j,kb,l) = zero
               enddo
             enddo
           enddo

         enddo domains

       enddo spectralLoop

       end subroutine initres
