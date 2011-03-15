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
subroutine initresNKPC(varStart, varEnd,wAdj,volAdj,dwAdj,nn,level,sps)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * initres initializes the given range of the residual. Either to *
  !      * zero, steady computation, or to an unsteady term for the time  *
  !      * spectral and unsteady modes. For the coarser grid levels the   *
  !      * residual forcing term is taken into account.                   *
  !      *                                                                *
  !      * This is a local routine, so assume that pointers are already   *
  !      * set.                                                           *
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

  real(kind=realType),dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral):: wAdj
  real(kind=realType),dimension(nw,nTimeIntervalsSpectral),  intent(inout) :: dwAdj
  real(kind=realType),dimension(nTimeIntervalsSpectral):: volAdj
  integer(kind=intType), intent(in)::nn,level,sps
  !
  !      Local variables.
  !
  integer(kind=intType) ::  mm, ll, ii, jj, i, j, k, l, m
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw)::  wspAdj
  real(kind=realType)  :: volspAdj
  !unsteady and timespectral variables
  real(kind=realType)   :: oneOverDt, tmp
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !

  ! Return immediately of no variables are in the range.

  if(varEnd < varStart) return


  ! Determine the equation mode and act accordingly.
  !print *,'equation Mode',equationMode,'ref',steady,timespectral,unsteady
  !switch to if statments. this particular case setup doesn't work
  !with tapenade. The steady case dissappears and Tapenade doesn't
  !know how to handle the empty case....
  !           select case (equationMode)
  !             case (steady)
  if(equationMode==steady)then
     ! Steady state computation.
     ! Determine the currently active multigrid level.
     steadyLevelTest: if(currentLevel == groundLevel) then

        do l=varStart,varEnd
           dwAdj(l,:) = zero
        enddo
     else steadyLevelTest
        call terminate("initResAD", &
             "ADjoint does not function on coarse grid level")

     endif steadyLevelTest

     !===========================================================

  elseif(equationMode==Unsteady)then
     !case (unsteady)

     ! Unsteady computation.
     ! A further distinction must be made.

     select case(timeIntegrationScheme)

     case (explicitRK)

        do l=varStart,varEnd
           dwAdj(l,:) = zero
        enddo

        !=======================================================

     case (implicitRK)
        call terminate("initRes", &
             "Implicit RK not implemented yet")

        !=======================================================

     case (BDF)
        !call terminate("initRes", &
        !     "BDF ADjoint not yet implemented")
        do l=varStart,varEnd
           dwAdj(l,:) = zero
        enddo
     end select

     !===========================================================
  elseif(equationMode==timespectral)then

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
           dwAdj(l,sps) = zero
        enddo
        ! Loop over the number of terms which contribute
        ! to the time derivative.

        timeLoopFine: do mm=1,nTimeIntervalsSpectral

           ! Store the pointer for the variable to be used to
           ! compute the unsteady source term and the volume.
           ! Also store in ii the offset needed for vector
           ! quantities.

           wspAdj   = wAdj(:,:,:,:,mm)
           volspAdj = volAdj(mm)!(:,:,:,mm)
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

                 ! Store the matrix vector product with the
                 ! velocity in tmp.

                 tmp = dvector(jj,ll,ii+1)*wspAdj(0,0,0,ivx) &
                      + dvector(jj,ll,ii+2)*wspAdj(0,0,0,ivy) &
                      + dvector(jj,ll,ii+3)*wspAdj(0,0,0,ivz)

                 ! Update the residual. Note the
                 ! multiplication with the density to obtain
                 ! the correct time derivative for the
                 ! momentum variable.


                 dwAdj(l,sps) = dwAdj(l,sps) &
                      + tmp*volspAdj*wspAdj(0,0,0,irho)


              else

                 ! Scalar variable.  Loop over the owned cells to
                 ! add the contribution of wsp to the time
                 ! derivative.

                 dwAdj(l,sps) = dwAdj(l,sps)        &
                      + dscalar(jj,sps,mm) &
                      * volspAdj*wspAdj(0,0,0,l)
              endif

           enddo varLoopFine

        enddo timeLoopFine

     else spectralLevelTest

        call terminate("initRes", &
             "Coarse levels not supported in ADjoint...")
     endif spectralLevelTest
  else
     call terminate("initResAdj", &
          "Not a valid equation Mode...")
  endif

end subroutine initresNKPC
