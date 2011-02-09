!
!      ******************************************************************
!      *                                                                *
!      * File:          wallComponentsAssiignment.f90                   *
!      * Author:        Eran Arad                                       *
!      * Starting date: 08-15-2009                                      *
!      * Last modified:                                                 *
!      *                                                                *
!      ******************************************************************
!
  subroutine wallComponentsAssignment(nZone,ibc,familyName)
!      ******************************************************************
!      *                                                                *
!      * Scan through wall bundary conditions of any type and create    *
!      * a list of all components (families) for                        *
!      * Component Break Down analysis                                  *
!      * As default, each component (wall family ) is set to contribute *
!      * to monitoring of total forces and moments. This setting can    *
!      * be modifies by the paramfile                                   *
!      * This routine is part of eran-cbd modifications                 *
!      *                                                                *
!      ******************************************************************
!
    use cgnsGrid  
    implicit none
!
!------- routine input
!
     integer(kind=intType) :: nZone,ibc
     character  (len=*), intent(in)  :: familyName
!
!------- local variables
!
     integer(kind=intType)  :: iname, idWBCL
     logical :: nameFound     
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
    cgnsDoms(nZone)%bocoInfo(ibc)%wallBCName = trim(familyName)

    nameFound = .false.

    if (cgnsNWallSurfaces >= 1)then
       do iname=1,cgnsNWallSurfaces 
          if (trim(familyName) == trim(WallBCNames(iname)))then
             nameFound = .true.
             idWBCL     = iname
             exit
          end if
       end do
    end if
    newName: if(.not.nameFound)then
       cgnsNWallSurfaces = cgnsNWallSurfaces + 1
       if(cgnsNWallSurfaces > 100)call terminate('readBocos ERROR:',&
            'Number of wall-type surfaces in grid > 100 (allocated size)')

       idWBCL        = cgnsNWallSurfaces
       WallBCNames(cgnsNWallSurfaces) = trim(familyName)
    end if newName
    cgnsDoms(nZone)%bocoInfo(ibc)%idWBC = idWBCL
    cgnsDoms(nZone)%bocoInfo(ibc)%contributeToForce = .true.                  

    return
  end subroutine wallComponentsAssignment


