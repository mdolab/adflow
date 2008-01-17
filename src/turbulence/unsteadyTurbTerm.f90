!
!      ******************************************************************
!      *                                                                *
!      * File:          unsteadyTurbTerm.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-09-2004                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine unsteadyTurbTerm(mAdv, nAdv, offset, qq)
!
!      ******************************************************************
!      *                                                                *
!      * unsteadyTurbTerm discretizes the time derivative of the        *
!      * turbulence transport equations and add it to the residual.     *
!      * As the time derivative is the same for all turbulence models,  *
!      * this generic routine can be used; both the discretization of   *
!      * the time derivative and its contribution to the central        *
!      * jacobian are computed by this routine.                         *
!      *                                                                *
!      * Only nAdv equations are treated, while the actual system has   *
!      * size mAdv. The reason is that some equations for some          *
!      * turbulence equations do not have a time derivative, e.g. the   *
!      * f equation in the v2-f model. The argument offset indicates    *
!      * the offset in the w vector where this subsystem starts. As a   *
!      * consequence it is assumed that the indices of the current      *
!      * subsystem are contiguous, e.g. if a 2*2 system is solved the   *
!      * Last index in w is offset+1 and offset+2 respectively.         *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use flowVarRefState
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use iteration
       use section
       use turbMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: mAdv, nAdv, offset

       real(kind=realType), dimension(2:il,2:jl,2:kl,mAdv,mAdv), &
                                                      intent(inout) :: qq
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ii, jj, nn

       real(kind=realType) :: oneOverDt, tmp
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the equation mode.

       select case (equationMode)

       ! case (steady)
         case (steady, timeSpectral)

           ! Steady computation. No time derivative present.

           return

         !===============================================================

         case (unsteady)

           ! The time deritvative term depends on the integration
           ! scheme used.

           select case (timeIntegrationScheme)

             case (BDF)

               ! Backward difference formula is used as time
               ! integration scheme.

               ! Store the inverse of the physical nonDimensional
               ! time step a bit easier.

               oneOverDt = timeRef/deltaT

               ! Loop over the number of turbulent transport equations.

               nAdvLoopUnsteady: do ii=1,nAdv

                 ! Store the index of the current turbulent variable in jj.

                 jj = ii + offset

                 ! Loop over the owned cells of this block to compute the
                 ! time derivative.

                 do k=2,kl
                   do j=2,jl
                     do i=2,il

                       ! Initialize tmp to the value of the current
                       ! level multiplied by the corresponding coefficient
                       ! in the time integration scheme.

                       tmp = coefTime(0)*w(i,j,k,jj)

                       ! Loop over the old time levels and add the
                       ! corresponding contribution to tmp.

                       do nn=1,noldLevels
                         tmp = tmp + coefTime(nn)*wold(nn,i,j,k,jj)
                       enddo

                       ! Update the residual. Note that in the turbulent
                       ! routines the residual is defined with an opposite
                       ! sign compared to the residual of the flow equations.
                       ! Therefore the time derivative must be substracted
                       ! from dvt.

                       dvt(i,j,k,ii) = dvt(i,j,k,ii) - oneOverDt*tmp

                       ! Update the central jacobian.

                       qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) &
                                       + coefTime(0)*oneOverDt
                     enddo
                   enddo
                 enddo

               enddo nAdvLoopUnsteady

             !===========================================================

             case (explicitRK)

               ! Explicit time integration scheme. The time derivative
               ! is handled differently.

               return

             !===========================================================

             case (implicitRK)

               call terminate("unsteadyTurbTerm", &
                              "Implicit RK not implemented yet")

           end select

         !===============================================================

       ! case (timeSpectral)
         case default

           ! Time spectral method.

           ! Loop over the number of turbulent transport equations.

           nAdvLoopSpectral: do ii=1,nAdv

             ! Store the index of the current turbulent variable in jj.

             jj = ii + offset

             ! The time derivative has been computed earlier in
             ! unsteadyTurbSpectral and stored in entry jj of dw.
             ! Substract this value for all owned cells. It must be
             ! substracted, because in the turbulent routines the
             ! residual is defined with an opposite sign compared to
             ! the residual of the flow equations.
             ! Also add a term to the diagonal matrix, which corresponds
             ! to to the contribution of the highest frequency. This is
             ! equivalent to an explicit treatment of the time derivative
             ! and may need to be changed.

             tmp = nTimeIntervalsSpectral*pi*timeRef &
                 / sections(sectionID)%timePeriod

             do k=2,kl
               do j=2,jl
                 do i=2,il
                   dvt(i,j,k,ii)   = dvt(i,j,k,ii)   - dw(i,j,k,jj)
                   qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) + tmp
                 enddo
               enddo
             enddo

           enddo nAdvLoopSpectral

       end select

       end subroutine unsteadyTurbTerm
