!
!      ******************************************************************
!      *                                                                *
!      * File:          initMemMixingPlane.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-27-2005                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initMemMixingPlane(level)
!
!      ******************************************************************
!      *                                                                *
!      * initMemMixingPlane initializes the memory for the mixing       *
!      * plane communication/interpolation pattern for the given grid   *
!      * level.                                                         *
!      *                                                                *
!      ******************************************************************
!
       use commMixing
       use interfaceGroups
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of interface groups and number of sides
       ! of an interface.

       groupLoop: do nn=1,nInterfaceGroups
         sideLoop: do mm=1,2

           ! Initialize all the numbers to 0.

           commPatternMixing(level,nn,mm)%nInter = 0
           commPatternMixing(level,nn,mm)%nDonor = 0
           commPatternMixing(level,nn,mm)%nHalo  = 0

           ! Nullify the pointers.

           nullify(commPatternMixing(level,nn,mm)%blockDonor)
           nullify(commPatternMixing(level,nn,mm)%indD)
           nullify(commPatternMixing(level,nn,mm)%rotMatDonor)
           nullify(commPatternMixing(level,nn,mm)%nIntervalsDonor)
           nullify(commPatternMixing(level,nn,mm)%weightDonor)
           nullify(commPatternMixing(level,nn,mm)%indListDonor)

           nullify(commPatternMixing(level,nn,mm)%indListHalo)
           nullify(commPatternMixing(level,nn,mm)%blockHalo)
           nullify(commPatternMixing(level,nn,mm)%indH)
           nullify(commPatternMixing(level,nn,mm)%weightHalo)
           nullify(commPatternMixing(level,nn,mm)%rotMatHalo)

         enddo sideLoop
       enddo groupLoop

       end subroutine initMemMixingPlane
