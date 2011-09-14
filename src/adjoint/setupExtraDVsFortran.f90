!
!     ******************************************************************
!     *                                                                *
!     * File:          setupExtraDVsFortran.f90                        *
!     * Authors:       C.A.(Sandy) Mader                               *
!     * Starting date: 09-14-2011                                      *
!     * Last modified: 09-14-2011                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine setupExtraDVsFortran()

!     *****************************************
!     * This subroutine sets up what extra DVs are included in the adjoint
!     ****************************************

      use ADjointVars
      use communication   ! myID
      implicit none

!
!     Local variables.
!

!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
      if(myID == 0) then
         print *,'Setting up Extra Variables...'
      endif
      
      nDesignAOA  = 0
      nDesignSSA    = -1 
      nDesignMach   = -1  
      nDesignMachGrid = -1
      nDesignRotX    = -1 
      nDesignRotY   = -1  
      nDesignRotZ   = -1  
      nDesignRotCenX  = -1
      nDesignRotCenY  = -1 
      nDesignRotCenZ  = -1  
      nDesignPointRefX   = -1
      nDesignPointRefY  = -1  
      nDesignPointRefZ  = -1  
      nDesignLengthRef = -1
      nDesignSurfaceRef = -1
      
      nDesignExtra = 1

      allocate(dIda(nDesignExtra))

    end subroutine setupExtraDVsFortran
