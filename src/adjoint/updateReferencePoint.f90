!
!     ******************************************************************
!     *                                                                *
!     * File:          updateReferencePoint.f90                        *
!     * Authors:       C.A(Sandy) Mader                                *
!     * Starting date: 11-08-2010                                      *
!     * Last modified: 11-08-2010                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine updateReferencePoint
  !
  !     ******************************************************************
  !     *                                                                *
  !     *  reruns the initialization routines to update AOA and other    *
  !     *  flow variables after a design change                          *
  !     *                                                                *
  !     ******************************************************************
  !

  use inputTimeSpectral
  use section
  use inputPhysics 
  use inputMotion  
  use cgnsGrid     
  use monitor
  use iteration
  use blockpointers
  implicit none

  ! Working variables
  integer(kind=intType) ::mm,nnn,nn
  real(kind=realType), dimension(nSections) :: t

  groundlevel = 1
  do mm=1,nTimeIntervalsSpectral
     do nn=1,nDom
        ! Set the pointers for this block.
        call setPointers(nn, groundLevel, mm)
        !lref is outside
        cgnsDoms(nbkglobal)%rotCenter = rotPoint
     enddo
  enddo

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
end subroutine updateReferencePoint
