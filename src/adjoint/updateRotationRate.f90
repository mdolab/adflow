!
!     ******************************************************************
!     *                                                                *
!     * File:          updateRotationRate.f90                          *
!     * Authors:       C.A(Sandy) Mader                                *
!     * Starting date: 02-09-2010                                      *
!     * Last modified: 02-09-2010                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine updateRotationRate(rotCenter, rotRate, blocks, nblocks)

  use inputTimeSpectral 
  use section
  use inputPhysics  
  use inputMotion   
  use cgnsGrid      
  use monitor
  use iteration
  use blockpointers
  implicit none

  real(kind=realType),intent(in)::rotCenter(3), rotRate(3)
  integer(kind=intType), intent(in) :: nblocks
  integer(kind=intType), intent(in) :: blocks(nblocks)

  integer(kind=intType) ::mm,nnn,nn, level, sps, i
  real(kind=realType), dimension(nSections) :: t

  groundlevel = 1

  do nn=1,nblocks
     cgnsDoms(nn)%rotRate = rotRate
     cgnsDoms(nn)%rotCenter = rotCenter

     do i=1,cgnsDoms(nn)%nBocos
        cgnsDoms(nn)%bocoInfo(i)%rotRate = rotRate
     end do
  enddo

  do sps=1,nTimeIntervalsSpectral
     do level=1,ubound(flowDoms,2)
        do nn=1,nDom
           flowDoms(nn,level,sps)%blockIsMoving = .True.
           flowDoms(nn,level,sps)%addGridVelocities = .True.
        end do
     end do
  end do
  groundlevel = 1
  do mm=1,nTimeIntervalsSpectral

     ! Compute the time, which corresponds to this spectral solution.
     ! For steady and unsteady mode this is simply the restart time;
     ! for the spectral mode the periodic time must be taken into
     ! account, which can be different for every section.

     t = timeUnsteadyRestart

     if(equationMode == timeSpectral) then
        do nnn=1,nSections
           t(nnn) = t(nnn) + (mm-1)*sections(nnn)%timePeriod &
                /         real(nTimeIntervalsSpectral,realType)
        enddo
     endif

     call gridVelocitiesFineLevel(.false., t, mm)
     call gridVelocitiesCoarseLevels(mm)
     call normalVelocitiesAllLevels(mm)
     call slipVelocitiesFineLevel(.false., t, mm)
     call slipVelocitiesCoarseLevels(mm)

  enddo

end subroutine updateRotationRate
