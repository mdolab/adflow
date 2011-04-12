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
subroutine initres_block_TS(varStart, varEnd,nn,sps,mm)
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
  use forwardAdjointVars
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

  real(kind=realType), dimension(:,:,:,:), pointer :: ww, wsp1, wsp
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


  ! Determine the equation mode and act accordingly.

  select case (equationMode)
  case (timeSpectral)

     ! Time spectral computation. The time derivative of the
     ! current solution is given by a linear combination of
     ! all other solutions, i.e. a matrix vector product.

     ! First store the section to which this block belongs
     ! in jj.

     jj = sectionID

     ! Determine the currently active multigrid level.

     ! Finest multigrid level. The residual must be
     ! initialized to the time derivative.
     
     ! Initialize it to zero.
     
     ! Loop over the number of terms which contribute
     ! to the time derivative.
     
     
     ! Store the pointer for the variable to be used to
     ! compute the unsteady source term and the volume.
     ! Also store in ii the offset needed for vector
     ! quantities.
     
     !wsp   => flowDoms(nn,currentLevel,mm)%w
     !volsp => flowDoms(nn,currentLevel,mm)%vol
     wsp   => w_offTimeInstance
     volsp => vol_offTimeInstance
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
        
     end do varLoopFine
  end select
end subroutine initres_block_TS

