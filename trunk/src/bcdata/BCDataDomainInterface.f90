!
!      ******************************************************************
!      *                                                                *
!      * File:          BCDataDomainInterface.f90                       *
!      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
!      * Starting date: 07-19-2005                                      *
!      * Last modified: 05-10-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine BCDataDomainInterface(boco)
!
!      ******************************************************************
!      *                                                                *
!      * BCDataDomainInterface initializes the boundary data for the    *
!      * domain interfaces to zero to avoid unexpected behavior later.  *
!      * The true values of the data is set by the coupler, although    *
!      * a better initial guess is made as soon as the field data has   *
!      * been initialized.                                              *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use couplerParam
       use inputPhysics
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType) :: boco
!
!      Local variables.
!
       real(kind=realType) :: asound, velMag
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       asound = sqrt(gammaConstant*pIni/rhoIni)
       velMag = asound*MachIni

       ! Initialize the allocated data to proper values.

       if( associated(BCData(boco)%rho)  ) BCData(boco)%rho  = rhoIni
       if( associated(BCData(boco)%velx) ) BCData(boco)%velx = velMag*velDirIni(1)
       if( associated(BCData(boco)%vely) ) BCData(boco)%vely = velMag*velDirIni(2)
       if( associated(BCData(boco)%velz) ) BCData(boco)%velz = velMag*velDirIni(3)
       if( associated(BCData(boco)%ps  ) ) BCData(boco)%ps   = pIni

       if( associated(BCData(boco)%flowXdirInlet) ) &
         BCData(boco)%flowXdirInlet = velDirIni(1)
       if( associated(BCData(boco)%flowYdirInlet) ) &
         BCData(boco)%flowYdirInlet = velDirIni(2)
       if( associated(BCData(boco)%flowZdirInlet) ) &
         BCData(boco)%flowZdirInlet = velDirIni(3)

       if( associated(BCData(boco)%ptInlet) ) &
         BCData(boco)%ptInlet = zero
       if( associated(BCData(boco)%ttInlet) ) &
         BCData(boco)%ttInlet = zero
       if( associated(BCData(boco)%htInlet) ) &
         BCData(boco)%htInlet = zero

       if( associated(BCData(boco)%turbInlet) ) &
         BCData(boco)%turbInlet = zero

       end subroutine BCDataDomainInterface
