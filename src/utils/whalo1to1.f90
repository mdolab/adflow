!
!      ******************************************************************
!      *                                                                *
!      * File:          whalo1to1.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-07-2003                                      *
!      * Last modified: 09-16-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine whalo1to1(level, start, end, commPressure,       &
                            commVarGamma, commLamVis, commEddyVis, &
                            commPattern, internal)
!
!      ******************************************************************
!      *                                                                *
!      * whalo1to1 exchanges the 1 to 1 internal halo's for the cell    *
!      * centered variables for the given communication pattern. It     *
!      * is possible to send a range of variables and not the entire    *
!      * set, e.g. only the flow variables or only the turbulent        *
!      * variables. This is controlled by the arguments start, end,     *
!      * commPressure and commViscous. The exchange takes place for     *
!      * the given grid level.                                          *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       use inputTimeSpectral
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, start, end
       logical, intent(in) :: commPressure, commVarGamma
       logical, intent(in) :: commLamVis, commEddyVis

       type(commType), dimension(*), intent(in)         :: commPattern
       type(internalCommType), dimension(*), intent(in) :: internal

       integer(kind=intType) :: nVar, nn, k, sps

       logical :: correctPeriodic
!
!      Interfaces
!
       interface
         subroutine correctPeriodicVelocity(level, sp, nPeriodic, &
                                            periodicData)
         use block
         use communication
         use constants
         implicit none

         integer(kind=intType), intent(in) :: level, sp, nPeriodic
         type(periodicDataType), dimension(:), pointer :: periodicData

         end subroutine correctPeriodicVelocity
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the logical correctPeriodic. Only if a momentum variable
       ! is communicated it is needed to apply the periodic
       ! transformations.


        correctPeriodic = .false.
        if(start <= ivx .and. end >= ivz) correctPeriodic = .true.

       spectralModes: do sps=1,nTimeIntervalsSpectral

          ! Set the pointers for the required variables
          domainLoop:do nn=1, nDom
             nVar = 0

             do k=start, end
                nVar = nVar + 1 
                flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                     flowDoms(nn, level, sps)%w(:, :, :, k)
             end do

             if( commPressure )  then 
                nVar = nVar + 1
                flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                     flowDoms(nn, level, sps)%P(:, :, :)
             end if

             if( commVarGamma ) then 
                nVar = nVar + 1
                flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                     flowDoms(nn, 1, sps)%gamma(:, :, :)
             end if

             if( commLamVis ) then 
                nVar = nVar + 1
                flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                     flowDoms(nn, 1, sps)%rlv(:, :, :)
             end if

             if( commEddyVis ) then 
                nVar = nVar + 1
                flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                     flowDoms(nn, level, sps)%rev(:, :, :)
             end if

             if(nVar == 0) return
          end do domainLoop

          ! Run the generic exchange
          call wHalo1to1RealGeneric(nVar, level, sps, commPattern, internal)

          if (correctPeriodic) then 
             if ( internal(level)%nPeriodic > 0 )   then 
                call correctPeriodicVelocity(level, sps,             & 
                     internal(level)%nPeriodic, internal(level)%periodicData)
             end if

             if ( commPattern(level)%nPeriodic > 0 )   then 
                call correctPeriodicVelocity(level, sps,             & 
                     commPattern(level)%nPeriodic, commPattern(level)%periodicData)
             end if
          end if

       end do spectralModes

     end subroutine whalo1to1

!      ==================================================================

       subroutine correctPeriodicVelocity(level, sp, nPeriodic, &
                                          periodicData)
!
!      ******************************************************************
!      *                                                                *
!      * correctPeriodicVelocity applies the periodic transformation    *
!      * to the velocity of the cell halo's in periodicData.            *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       use constants
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sp, nPeriodic
       type(periodicDataType), dimension(:), pointer :: periodicData
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, ii, i, j, k
       real(kind=realType)   :: vx, vy, vz

       real(kind=realType), dimension(3,3) :: rotMatrix
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of periodic transformations.

       do nn=1,nPeriodic

         ! Store the rotation matrix a bit easier.

         rotMatrix = periodicData(nn)%rotMatrix

         ! Loop over the number of halo cells for this transformation.

         do ii=1,periodicData(nn)%nhalos

           ! Store the block and the indices a bit easier.

           mm = periodicData(nn)%block(ii)
           i  = periodicData(nn)%indices(ii,1)
           j  = periodicData(nn)%indices(ii,2)
           k  = periodicData(nn)%indices(ii,3)

           ! Store the original velocities in vx, vy, vz.

           vx = flowDoms(mm,level,sp)%w(i,j,k,ivx)
           vy = flowDoms(mm,level,sp)%w(i,j,k,ivy)
           vz = flowDoms(mm,level,sp)%w(i,j,k,ivz)

           ! Compute the new velocity vector.

           flowDoms(mm,level,sp)%w(i,j,k,ivx) = rotMatrix(1,1)*vx &
                                              + rotMatrix(1,2)*vy &
                                              + rotMatrix(1,3)*vz
           flowDoms(mm,level,sp)%w(i,j,k,ivy) = rotMatrix(2,1)*vx &
                                              + rotMatrix(2,2)*vy &
                                              + rotMatrix(2,3)*vz
           flowDoms(mm,level,sp)%w(i,j,k,ivz) = rotMatrix(3,1)*vx &
                                              + rotMatrix(3,2)*vy &
                                              + rotMatrix(3,3)*vz
         enddo

       enddo

       end subroutine correctPeriodicVelocity

