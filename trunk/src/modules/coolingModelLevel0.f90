!
!      ******************************************************************
!      *                                                                *
!      * File:          coolingModelLevel0.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-27-2005                                      *
!      * Last modified: 04-29-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module coolingModelLevel0
!
!      ******************************************************************
!      *                                                                *
!      * Module which contains the derived data type as well as the     *
!      * corresponding array to store the information for the level 0   *
!      * cooling model used in the turbine.                             *
!      *                                                                *
!      * This model has been developed by Pratt and Whitney and should  *
!      * not be given to third parties. This implementation assumes     *
!      * the x-direction is the axial direction.                        *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived data type level0CoolingType.     *
!      *                                                                *
!      ******************************************************************
!
       type level0CoolingType

         ! nSubfaces: Local number of subfaces which are part of this
         !            cooling plane.
         ! mDotRatio: Relative amount of cooling mass flow added.
         ! dpLog:     Log of the total pressure loss coefficient.
         ! dTLog:     Log of the total temperature loss coefficient.
         ! area:      Reference area. Based on the entire wheel such that
         !            different periodic angles are treated correctly.

         ! blockID(:):  The local block ID of the subfaces.
         ! indexDir(:): The index direction of this subface, iMin, Imax,
         !              jMin, jMax, kMin or kMin. These names are a bit
         !              abused here, because it can also be an internal
         !              plane.
         ! indSol(:):   The cell index of the cell centered variables
         !              whose residuals must be modified.
         ! indNorm(:):  The face index for the normals.
         ! indX1(:):    The face index of the 1st adjacent plane to the
         !              cell center; normally identical to indNorm.
         ! indX2(:):    The face index of the 2nd adjacent plane.
         ! jcBeg(:):    Start index of the cell in j-direction of the
         !              subface.
         ! jcEnd(:):    Idem, but then for the end index.
         ! icBeg(:):    Idem, but then starting index in i-direction.
         ! icEnd(:):    Idem, but then for the end index.

         integer(kind=intType) :: nSubfaces

         integer(kind=intType), dimension(:), pointer :: blockID
         integer(kind=intType), dimension(:), pointer :: indexDir
         integer(kind=intType), dimension(:), pointer :: indSol
         integer(kind=intType), dimension(:), pointer :: indNorm
         integer(kind=intType), dimension(:), pointer :: indX1
         integer(kind=intType), dimension(:), pointer :: indX2
         integer(kind=intType), dimension(:), pointer :: jcBeg, jcEnd
         integer(kind=intType), dimension(:), pointer :: icBeg, icEnd

         real(kind=realType) :: mDotRatio, dpLog, dTLog
         real(kind=realType) :: area

       end type level0CoolingType
!
!      ******************************************************************
!      *                                                                *
!      * Variables stored in this module.                               *
!      *                                                                *
!      ******************************************************************
!
       ! nPlanesLevel0CoolingModel:  Number of injection planes used for
       !                             the level 0 cooling model.
       ! level0Cooling(:,nMGLevels): The array of the derived datatype
       !                             to store the info to implement the
       !                             cooling model. The first dimension
       !                             is nPlanesLevel0CoolingModel.

       integer(kind=intType) :: nPlanesLevel0CoolingModel

       type(level0CoolingType), dimension(:,:), allocatable :: &
                                                           level0Cooling

       end module coolingModelLevel0
