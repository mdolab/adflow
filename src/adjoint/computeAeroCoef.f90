!
!     ******************************************************************
!     *                                                                *
!     * File:          computeAeroCoef.f90                             *
!     * Author:        Andre C. Marta,C.A.(Sandy) Mader                *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 01-17-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine computeAeroCoef(CL,CD,CFx,CFy,CFz,CMx,CMy,CMz,level,sps)
!
!     ******************************************************************
!     *                                                                *
!     * Compute the aerodynamic coefficients from the force and moment *
!     * produced by the pressure and shear stresses on the body walls: *
!     *                                                                *
!     *   CL  - Lift coefficient                                       *
!     *   CD  - Drag coefficient                                       *
!     *   CFx - Force coefficient in x-direction                       *
!     *   CFy - Force coefficient in y-direction                       *
!     *   CFz - Force coefficient in z-direction                       *
!     *   CMx - Moment coefficient in x-direction                      *
!     *   CMy - Moment coefficient in y-direction                      *
!     *   CMz - Moment coefficient in z-direction                      *
!     *                                                                *
!     ******************************************************************
!
      use blockPointers  ! nDom
      use communication  ! my_ID SUmb_comm_world
      use inputPhysics   ! liftDirection, dragDirection
      use iteration      ! groundLevel
      !use monitor        ! monLoc, monGlob, nMonSum
      use precision      ! sumb_real
      implicit none
!
!     Subroutine arguments.
!
      real(kind=realType), intent(out)  :: CL, CD,CFx, CFy,CFz, CMx, CMy, CMz
      integer(kind=intType), intent(in) :: level, sps
!
!     Local variables.
!
      integer :: ierr, nCoefSum,nn
      real(kind=realType) :: yplusMax
      real(kind=realType), dimension(3) :: cFp, cFv, cMp, cMv
      real(kind=realType), dimension(8) :: coefLoc,coefGlob
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!

      groundLevel = level
      !Reset the local storage to zero
      coefLoc(:) = 0
      domains: do nn=1,nDom
         !print *,'domain',nn
         ! Set the pointers for this block.
         !print *,'setting pointers',nn, groundLevel, sps
         call setPointers(nn, groundLevel, sps)
         
         ! Compute the forces and moments for this block.
         !print *,'calling forces'
         call forcesAndMoments(cfp, cfv, cmp, cmv, yplusMax)
         !print *,'summing'
         coefLoc(1) =coefLoc(1) + (cFp(1) + cFv(1))*liftDirection(1) &
              + (cFp(2) + cFv(2))*liftDirection(2) &
              + (cFp(3) + cFv(3))*liftDirection(3)
         coefLoc(2) = coefLoc(2) + (cFp(1) + cFv(1))*dragDirection(1) &
              + (cFp(2) + cFv(2))*dragDirection(2) &
              + (cFp(3) + cFv(3))*dragDirection(3)
         coefLoc(3) =coefLoc(3) +  cFp(1) + cFv(1)
         coefLoc(4) =coefLoc(4) +  cFp(2) + cFv(2)
         coefLoc(5) = coefLoc(5) + cFp(3) + cFv(3)
         coefLoc(6) =coefLoc(6) +   cMp(1) + cMv(1)
         coefLoc(7) = coefLoc(7) +  cMp(2) + cMv(2)
         coefLoc(8) = coefLoc(8) +  cMp(3) + cMv(3)
         nCoefSum = 8

      enddo domains

      !print *,'loop finished...reducing'
      ! Determine the global sum of the summation monitoring
      ! variables. The sum is made known to all processors.

      call mpi_allreduce(coefLoc, coefGlob, nCoefSum, sumb_real, &
                         mpi_sum, SUmb_comm_world, ierr)

      ! Transfer the cost function values to output arguments.

      CL  = coefGlob(1)
      CD  = coefGlob(2)
      CFx = coefGlob(3)
      CFy = coefGlob(4)
      CFz = coefGlob(5)
      CMx = coefGlob(6)
      CMy = coefGlob(7)
      CMz = coefGlob(8)

      end subroutine computeAeroCoef
