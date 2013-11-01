!
!      ******************************************************************
!      *                                                                *
!      * File:          setDiffSizes.f90                                *
!      * Author:        Peter Zhoujie Lyu	                        *
!      * Starting date: 03-27-2013                                      *
!      * Last modified: 03-27-2013                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine setDiffSizes
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This routine set the sizes for the pointers that will be       *
  !      * used in the forward debug mode and reverse mode AD.            *
  !      *                                                                *
  !      ******************************************************************
  !
  use diffSizes
  use blockPointers
  use constants
  use flowVarRefState
  use inputTimeSpectral
  use inputPhysics
  use costFunctions
  implicit none

  ! local variables
  integer(kind=intType) :: nLevels

  ! Compute nlevels
  nLevels = ubound(flowDoms, 2)

  ! Set the size for dynamic pointers to zero for debug purpose
  ! bcdata%norm
  ISIZE1OFDrfDrfbcdata_norm = 0
  ISIZE2OFDrfDrfbcdata_norm = 0
  ISIZE3OFDrfDrfbcdata_norm = 0

  ! bcdata%rface
  ISIZE1OFDrfDrfbcdata_rface = 0
  ISIZE2OFDrfDrfbcdata_rface = 0

  ! bcdata%m
  ISIZE1OFDrfDrfbcdata_m = 0
  ISIZE2OFDrfDrfbcdata_m = 0
  ISIZE3OFDrfDrfbcdata_m = 0

  ! bcdata%fp
  ISIZE1OFDrfDrfbcdata_fp = 0
  ISIZE2OFDrfDrfbcdata_fp = 0
  ISIZE3OFDrfDrfbcdata_fp = 0

  ! bcdata%fv
  ISIZE1OFDrfDrfbcdata_fv = 0
  ISIZE2OFDrfDrfbcdata_fv = 0
  ISIZE3OFDrfDrfbcdata_fv = 0

  ! bcdata%fp
  ISIZE1OFDrfDrfbcdata_oarea = 0
  ISIZE2OFDrfDrfbcdata_oarea = 0
  ISIZE3OFDrfDrfbcdata_oarea = 0

  ! viscsubface%tau
  ISIZE1OFDrfDrfviscsubface_tau = 0
  ISIZE2OFDrfDrfviscsubface_tau = 0
  ISIZE3OFDrfDrfviscsubface_tau = 0
  
  ! prod
  ISIZE1OFDrfprod = 0
  ISIZE2OFDrfprod = 0
  ISIZE3OFDrfprod = 0
 
  ! vort
  ISIZE1OFDrfvort = 0
  ISIZE2OFDrfvort = 0
  ISIZE3OFDrfvort = 0

  ! dvt
  ISIZE1OFDrfdvt = 0
  ISIZE2OFDrfdvt = 0
  ISIZE3OFDrfdvt = 0
  ISIZE4OFDrfdvt = 0

  ! vol
  ISIZE1OFDrfvol = 0
  ISIZE2OFDrfvol = 0
  ISIZE3OFDrfvol = 0

  ! rho,etot,u,v,w,p,k
  ISIZE1OFrho = 0
  ISIZE1OFetot = 0
  ISIZE1OFu = 0
  ISIZE1OFv = 0
  ISIZE1OFw = 0
  ISIZE1OFp = 0
  ISIZE1OFk = 0

  ! Du1, Du2, Du3
  ISIZE1OFDu1 = 0
  ISIZE1OFDu2 = 0
  ISIZE1OFDu3 = 0

  ! Left, Right, Flux
  ISIZE1OFLeft = 0
  ISIZE1OFRight = 0
  ISIZE1OFFlux = 0

  ! bcdata
  ISIZE1OFDrfbcdata = 0!nbocos

  ! s
  ISIZE1OFDrfs = 0!ie
  ISIZE2OFDrfs = je
  ISIZE3OFDrfs = ke
  ISIZE4OFDrfs = 3

  ! sfacei
  ISIZE3OFDrfsfaceI = 0!ie + 1
  ISIZE2OFDrfsfaceI = je
  ISIZE1OFDrfsfaceI = ke

  ! sfacej
  ISIZE3OFDrfsfaceJ = 0!ie
  ISIZE2OFDrfsfaceJ = je + 1
  ISIZE1OFDrfsfaceJ = ke

  ! sfacek
  ISIZE3OFDrfsfaceK = 0!ie 
  ISIZE2OFDrfsfaceK = je
  ISIZE1OFDrfsfaceK = ke + 1

  ! Define size for the pointers
  ! flowdoms
  ISIZE1OFDrfflowdoms = nDom
  ISIZE2OFDrfflowdoms = nLevels
  ISIZE3OFDrfflowdoms = nTimeIntervalsSpectral
  
  !viscSubface
  ISIZE1OFDrfviscsubface = nViscBocos

  ! x
  ISIZE4OFDrfx = 3
  ISIZE1OFDrfx = ie + 1
  ISIZE2OFDrfx = je + 1
  ISIZE3OFDrfx = ke + 1
  
  ! flowdoms_x
  ISIZE4OFDRFFLOWDOMS_X = 3
  ISIZE1OFDRFFLOWDOMS_X = ie + 1
  ISIZE2OFDRFFLOWDOMS_X = je + 1
  ISIZE3OFDRFFLOWDOMS_X = ke + 1

  if ( equations == RANSEquations ) then
     ! rev
     ISIZE1OFDrfrev = ib + 1
     ISIZE2OFDrfrev = jb + 1
     ISIZE3OFDrfrev = kb + 1
  else
     ! rev
     ISIZE1OFDrfrev = 0
     ISIZE2OFDrfrev = 0
     ISIZE3OFDrfrev = 0
  end if

  ! rlv
  ISIZE1OFDrfrlv = ib + 1
  ISIZE2OFDrfrlv = jb + 1
  ISIZE3OFDrfrlv = kb + 1

  ! w
  ISIZE4OFDrfw = nw
  ISIZE1OFDrfw = ib + 1
  ISIZE2OFDrfw = jb + 1
  ISIZE3OFDrfw = kb + 1

  ! flowdoms_x
  ISIZE4OFDRFFLOWDOMS_W = nw
  ISIZE1OFDRFFLOWDOMS_W = ib + 1
  ISIZE2OFDRFFLOWDOMS_W = jb + 1
  ISIZE3OFDRFFLOWDOMS_W = kb + 1

  ! flowdoms_dw
  ISIZE4OFDRFFLOWDOMS_dw = nw
  ISIZE1OFDRFFLOWDOMS_dw = ib + 1
  ISIZE2OFDRFFLOWDOMS_dw = jb + 1
  ISIZE3OFDRFFLOWDOMS_dw = kb + 1

  ! flowdoms_vol
  ISIZE1OFDRFFLOWDOMS_vol = ib + 1
  ISIZE2OFDRFFLOWDOMS_vol = jb + 1
  ISIZE3OFDRFFLOWDOMS_vol = kb + 1

  ! fw
  ISIZE4OFDrffw = nwf
  ISIZE1OFDrffw = ib + 1
  ISIZE2OFDrffw = jb + 1
  ISIZE3OFDrffw = kb + 1

  ! dw
  ISIZE4OFDrfdw = nw
  ISIZE1OFDrfdw = ib + 1
  ISIZE2OFDrfdw = jb + 1
  ISIZE3OFDrfdw = kb + 1

  ! p
  ISIZE1OFDrfp = ib + 1
  ISIZE2OFDrfp = jb + 1
  ISIZE3OFDrfp = kb + 1

  ! gamma
  ISIZE1OFDrfgamma = ib + 1
  ISIZE2OFDrfgamma = jb + 1
  ISIZE3OFDrfgamma = kb + 1

  ! radI
  ISIZE1OFDrfradI = ie
  ISIZE2OFDrfradI = je
  ISIZE3OFDrfradI = ke

  ! radJ
  ISIZE1OFDrfradJ = ie
  ISIZE2OFDrfradJ = je
  ISIZE3OFDrfradJ = ke
  
  ! radK
  ISIZE1OFDrfradK = ie
  ISIZE2OFDrfradK = je
  ISIZE3OFDrfradK = ke
  
  ! sI
  ISIZE1OFDRFsI = ie + 1
  ISIZE2OFDRFsI = je
  ISIZE3OFDRFsI = ke
  ISIZE4OFDRFsI = 3

  ! sJ
  ISIZE1OFDRFsJ = ie
  ISIZE2OFDRFsJ = je + 1
  ISIZE3OFDRFsJ = ke
  ISIZE4OFDRFsJ = 3

  ! sK
  ISIZE1OFDRFsK = ie
  ISIZE2OFDRFsK = je
  ISIZE3OFDRFsK = ke + 1
  ISIZE4OFDRFsK = 3
  
  !bmti1
  ISIZE1OFDrfbmti1 = je
  ISIZE2OFDrfbmti1 = ke
  ISIZE3OFDrfbmti1 = nt2 - nt1 + 1
  ISIZE4OFDrfbmti1 = nt2 - nt1 + 1 

  !bmti2
  ISIZE1OFDrfbmti2 = je
  ISIZE2OFDrfbmti2 = ke
  ISIZE3OFDrfbmti2 = nt2 - nt1 + 1
  ISIZE4OFDrfbmti2 = nt2 - nt1 + 1 

  !bmtj1
  ISIZE1OFDrfbmtj1 = ie
  ISIZE2OFDrfbmtj1 = ke
  ISIZE3OFDrfbmtj1 = nt2 - nt1 + 1
  ISIZE4OFDrfbmtj1 = nt2 - nt1 + 1

  !bmtj2
  ISIZE1OFDrfbmtj2 = ie
  ISIZE2OFDrfbmtj2 = ke
  ISIZE3OFDrfbmtj2 = nt2 - nt1 + 1
  ISIZE4OFDrfbmtj2 = nt2 - nt1 + 1

  !bmtk1
  ISIZE1OFDrfbmtk1 = ie
  ISIZE2OFDrfbmtk1 = je
  ISIZE3OFDrfbmtk1 = nt2 - nt1 + 1
  ISIZE4OFDrfbmtk1 = nt2 - nt1 + 1

  !bmtk2
  ISIZE1OFDrfbmtk2 = ie
  ISIZE2OFDrfbmtk2 = je
  ISIZE3OFDrfbmtk2 = nt2 - nt1 + 1
  ISIZE4OFDrfbmtk2 = nt2 - nt1 + 1

  !bvti1
  ISIZE1OFDrfbvti1 = je
  ISIZE2OFDrfbvti1 = ke
  ISIZE3OFDrfbvti1 = nt2 - nt1 + 1

  !bvti2
  ISIZE1OFDrfbvti2 = je
  ISIZE2OFDrfbvti2 = ke
  ISIZE3OFDrfbvti2 = nt2 - nt1 + 1
  !bvti1
  ISIZE1OFDrfbvti1 = je
  ISIZE2OFDrfbvti1 = ke
  ISIZE3OFDrfbvti1 = nt2 - nt1 + 1

  !bvti2
  ISIZE1OFDrfbvti2 = je
  ISIZE2OFDrfbvti2 = ke
  ISIZE3OFDrfbvti2 = nt2 - nt1 + 1

  !bvtk1
  ISIZE1OFDrfbvtk1 = ie
  ISIZE2OFDrfbvtk1 = je
  ISIZE3OFDrfbvtk1 = nt2 - nt1 + 1

  !bvtk2
  ISIZE1OFDrfbvtk2 = ie
  ISIZE2OFDrfbvtk2 = je
  ISIZE3OFDrfbvtk2 = nt2 - nt1 + 1

end subroutine setDiffSizes
