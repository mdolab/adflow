  module diffSizes
    use precision
    implicit none
    save

    ! These are the diff sizes reqruied for the forward mode AD
    integer(kind=intType), parameter :: ISIZE3ofviscsubface = 3
    integer(kind=intType) :: ISIZE1OFDrfbcdata

    integer(kind=intType) :: ISIZE1OFDrfviscsubface
    integer(kind=intType) :: ISIZE1OFDrfflowdoms
    integer(kind=intType) :: ISIZE2OFDrfflowdoms 
    integer(kind=intType) :: ISIZE3OFDrfflowdoms
    integer(kind=intType) :: ISIZE1OFDrfflowdoms_bcdata

    ! These are the diff sizes reqruied for the reverse mode AD
    integer(kind=intType) :: ISIZE3OFDrfrlv, ISIZE2OFDrfrlv, ISIZE1OFDrfrlv
    integer(kind=intType) :: ISIZE4OFDrfw, ISIZE3OFDrfw, ISIZE2OFDrfw, ISIZE1OFDrfw
    integer(kind=intType) :: ISIZE4OFDrffw, ISIZE3OFDrffw, ISIZE2OFDrffw, ISIZE1OFDrffw
    integer(kind=intType) :: ISIZE4OFDrfdw, ISIZE3OFDrfdw, ISIZE2OFDrfdw, ISIZE1OFDrfdw
    integer(kind=intType) :: ISIZE3OFDrfgamma, ISIZE2OFDrfgamma, ISIZE1OFDrfgamma
    integer(kind=intType) :: ISIZE3OFDrfp, ISIZE2OFDrfp, ISIZE1OFDrfp
    integer(kind=intType) :: ISIZE3OFDrfrev, ISIZE2OFDrfrev, ISIZE1OFDrfrev
    integer(kind=intType) :: ISIZE3OFDrfsfaceI, ISIZE2OFDrfsfaceI, ISIZE1OFDrfsfaceI
    integer(kind=intType) :: ISIZE3OFDrfsfaceJ, ISIZE2OFDrfsfaceJ, ISIZE1OFDrfsfaceJ
    integer(kind=intType) :: ISIZE3OFDrfsfaceK, ISIZE2OFDrfsfaceK, ISIZE1OFDrfsfaceK
    integer(kind=intType) :: ISIZE3OFDrfradI, ISIZE2OFDrfradI, ISIZE1OFDrfradI
    integer(kind=intType) :: ISIZE3OFDrfradJ, ISIZE2OFDrfradJ, ISIZE1OFDrfradJ
    integer(kind=intType) :: ISIZE3OFDrfradK, ISIZE2OFDrfradK, ISIZE1OFDrfradK

    integer(kind=intType) :: ISIZE4OFDRFsI, ISIZE3OFDRFsI, ISIZE2OFDRFsI, ISIZE1OFDRFsI
    integer(kind=intType) :: ISIZE4OFDRFsJ, ISIZE3OFDRFsJ, ISIZE2OFDRFsJ, ISIZE1OFDRFsJ
    integer(kind=intType) :: ISIZE4OFDRFsK, ISIZE3OFDRFsK, ISIZE2OFDRFsK, ISIZE1OFDRFsK

    integer(kind=intType) :: ISIZE4OFDRFFLOWDOMS_X, ISIZE3OFDRFFLOWDOMS_X
    integer(kind=intType) :: ISIZE2OFDRFFLOWDOMS_X, ISIZE1OFDRFFLOWDOMS_X

    integer(kind=intType) :: ISIZE4OFDRFFLOWDOMS_W, ISIZE3OFDRFFLOWDOMS_W
    integer(kind=intType) :: ISIZE2OFDRFFLOWDOMS_W, ISIZE1OFDRFFLOWDOMS_W

    integer(kind=intType) :: ISIZE4OFDrfbmtj2
    integer(kind=intType) :: ISIZE3OFDrfbmtj2, ISIZE3OFDrfbvtj2
    integer(kind=intType) :: ISIZE2OFDrfbmtj2, ISIZE2OFDrfbvtj2
    integer(kind=intType) :: ISIZE1OFDrfbmtj2, ISIZE1OFDrfbvtj2

    integer(kind=intType) :: ISIZE4OFDrfbmtj1
    integer(kind=intType) :: ISIZE3OFDrfbmtj1, ISIZE3OFDrfbvtj1
    integer(kind=intType) :: ISIZE2OFDrfbmtj1, ISIZE2OFDrfbvtj1
    integer(kind=intType) :: ISIZE1OFDrfbmtj1, ISIZE1OFDrfbvtj1

    integer(kind=intType) :: ISIZE4OFDrfbmti2
    integer(kind=intType) :: ISIZE3OFDrfbmti2, ISIZE3OFDrfbvti2
    integer(kind=intType) :: ISIZE2OFDrfbmti2, ISIZE2OFDrfbvti2
    integer(kind=intType) :: ISIZE1OFDrfbmti2, ISIZE1OFDrfbvti2

    integer(kind=intType) :: ISIZE4OFDrfbmti1
    integer(kind=intType) :: ISIZE3OFDrfbmti1, ISIZE3OFDrfbvti1
    integer(kind=intType) :: ISIZE2OFDrfbmti1, ISIZE2OFDrfbvti1
    integer(kind=intType) :: ISIZE1OFDrfbmti1, ISIZE1OFDrfbvti1

    integer(kind=intType) :: ISIZE4OFDrfbmtk2
    integer(kind=intType) :: ISIZE3OFDrfbmtk2, ISIZE3OFDrfbvtk2
    integer(kind=intType) :: ISIZE2OFDrfbmtk2, ISIZE2OFDrfbvtk2
    integer(kind=intType) :: ISIZE1OFDrfbmtk2, ISIZE1OFDrfbvtk2

    integer(kind=intType) :: ISIZE4OFDrfbmtk1
    integer(kind=intType) :: ISIZE3OFDrfbmtk1, ISIZE3OFDrfbvtk1
    integer(kind=intType) :: ISIZE2OFDrfbmtk1, ISIZE2OFDrfbvtk1
    integer(kind=intType) :: ISIZE1OFDrfbmtk1, ISIZE1OFDrfbvtk1

    ! These are the diff sizes reqruied for the forward mode debug
    integer(kind=intType) :: ISIZE1OFDrfDrfbcdata_m, ISIZE2OFDrfDrfbcdata_m, ISIZE3OFDrfDrfbcdata_m
    integer(kind=intType) :: ISIZE1OFDrfDrfbcdata_fp, ISIZE2OFDrfDrfbcdata_fp, ISIZE3OFDrfDrfbcdata_fp
    integer(kind=intType) :: ISIZE1OFDrfDrfbcdata_fv, ISIZE2OFDrfDrfbcdata_fv, ISIZE3OFDrfDrfbcdata_fv
    integer(kind=intType) :: ISIZE1OFDrfDrfbcdata_norm, ISIZE2OFDrfDrfbcdata_norm, ISIZE3OFDrfDrfbcdata_norm
    integer(kind=intType) :: ISIZE1OFDrfDrfbcdata_oarea, ISIZE2OFDrfDrfbcdata_oarea, ISIZE3OFDrfDrfbcdata_oarea
    integer(kind=intType) :: ISIZE1OFDrfDrfbcdata_sepsensor, ISIZE2OFDrfDrfbcdata_sepSensor

    integer(kind=intType) :: ISIZE1OFDrfDrfbcdata_rface, ISIZE2OFDrfDrfbcdata_rface, ISIZE3OFDrfDrfbcdata_rface
    integer(kind=intType) :: ISIZE1OFDrfDrfbcdata_uslip, ISIZE2OFDrfDrfbcdata_uslip, ISIZE3OFDrfDrfbcdata_uslip
    integer(kind=intType) :: ISIZE1OFDrfprod, ISIZE2OFDrfprod, ISIZE3OFDrfprod
    integer(kind=intType) :: ISIZE1OFDrfdvt, ISIZE2OFDrfdvt, ISIZE3OFDrfdvt, ISIZE4OFDrfdvt
    integer(kind=intType) :: ISIZE1OfDrfvort, ISIZE2OFDrfvort, ISIZE3OFDrfvort
    integer(kind=intType) :: ISIZE1OFDrfvol, ISIZE2OFDrfvol, ISIZE3OFDrfvol

    integer(kind=intType) :: ISIZE1OFDrfDrfviscsubface_tau, ISIZE2OFDrfDrfviscsubface_tau
    integer(kind=intType) :: ISIZE3OFDrfDrfviscsubface_tau

    integer(kind=intType) :: ISIZE1OFDRFFLOWDOMS_dw, ISIZE2OFDRFFLOWDOMS_dw
    integer(kind=intType) :: ISIZE3OFDRFFLOWDOMS_dw, ISIZE4OFDRFFLOWDOMS_dw

    integer(kind=intType) :: ISIZE1OFDRFFLOWDOMS_vol, ISIZE2OFDRFFLOWDOMS_vol
    integer(kind=intType) :: ISIZE3OFDRFFLOWDOMS_vol, ISIZE4OFDRFFLOWDOMS_vol

    integer(kind=intType) :: ISIZE4OFDRFX, ISIZE3OFDRFX
    integer(kind=intType) :: ISIZE2OFDRFX, ISIZE1OFDRFX

    integer(kind=intType) :: ISIZE4OFDRFS, ISIZE3OFDRFS
    integer(kind=intType) :: ISIZE2OFDRFS, ISIZE1OFDRFS

    integer(kind=intType) :: ISIZE1OFrho, ISIZE1OFetot, ISIZE1OFu, ISIZE1OFv
    integer(kind=intType) :: ISIZE1OFw, ISIZE1OFp, ISIZE1OFk

    integer(kind=intType) :: ISIZE1OFDu1, ISIZE1OFDu2, ISIZE1OFDu3
    integer(kind=intType) :: ISIZE1OFLeft, ISIZE1OFRight, ISIZE1OFFlux
  end module diffSizes
