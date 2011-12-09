!     ******************************************************************
!     *                                                                *
!     * File:          matrixInfio.F90                                 *
!     * Author:        Gaetan K.W. Kenway                              *
!     * Starting date: 11-05-2010                                      *
!     * Last modified: 11-05-2010                                      *
!     *                                                                *
!     ******************************************************************

! This function takes the import PETSc Matrices, dRdwT, dRdwPreT,
! dRdx, dFdw and dFdx and outputs nicely formatted information about
! them. This is useful for debugging purposes. 

subroutine matrixInfo(pdRdwT,pdRdwPreT,pdRdx,pdRda,&
     pdFdw,pdFdx,pLocal,pSum,pMax)
#ifndef USE_NO_PETSC
  use ADjointPETSc
  use ADjointVars     ! nCellsLocal,nNodesLocal, nDesignExtra
  use communication   ! myID, nProc
  use inputTimeSpectral !nTimeIntervalsSpectral
  use flowVarRefState ! 
  use inputADjoint    !ApproxPC

  implicit None

  ! Logicals on wether or not to print these matrix
  logical , intent(in) :: pdRdwT,pdRdwPreT,pdRdx,pdFdx,pdFdw,pdRda
  logical , intent(in) :: pLocal,pSum,pMax
  integer(kind=intType) :: ierr
  if (pdRdwT) then

     if (pLocal) then 
        call MatGetInfo(dRdwT,MAT_LOCAL,localInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printLocal("dRdwT")
     end if

     if (pMax) then 
        call MatGetInfo(dRdwT,MAT_GLOBAL_MAX,maxInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printMax("dRdwT")
     end if

     if (pSum) then
        call MatGetInfo(dRdwT,MAT_GLOBAL_SUM,sumInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printSum("dRdwT")
     end if
  end if

  if (pdRdwPreT .and. ApproxPC) then

     if (pLocal) then 
        call MatGetInfo(dRdwPreT,MAT_LOCAL,localInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printLocal("dRdwPreT")
     end if

     if (pMax) then 
        call MatGetInfo(dRdwPreT,MAT_GLOBAL_MAX,maxInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printMax("dRdwPreT")
     end if

     if (pSum) then
        call MatGetInfo(dRdwPreT,MAT_GLOBAL_SUM,sumInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printSum("dRdwPreT")
     end if

  end if

  if (pdRdx) then

     if (pLocal) then 
        call MatGetInfo(dRdx,MAT_LOCAL,localInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printLocal("dRdx")
     end if

     if (pMax) then 
        call MatGetInfo(dRdx,MAT_GLOBAL_MAX,maxInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printMax("dRdx")
     end if

     if (pSum) then
        call MatGetInfo(dRdx,MAT_GLOBAL_SUM,sumInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printSum("dRdx")
     end if

  end if

  if (pdRda) then

     if (pLocal) then 
        call MatGetInfo(dRda,MAT_LOCAL,localInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printLocal("dRda")
     end if

     if (pMax) then 
        call MatGetInfo(dRda,MAT_GLOBAL_MAX,maxInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printMax("dRda")
     end if

     if (pSum) then
        call MatGetInfo(dRda,MAT_GLOBAL_SUM,sumInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printSum("dRda")
     end if

  end if


  if (pdFdx) then

     if (pLocal) then 
        call MatGetInfo(dFdx,MAT_LOCAL,localInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printLocal("dFdx")
     end if

     if (pMax) then 
        call MatGetInfo(dFdx,MAT_GLOBAL_MAX,maxInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printMax("dFdx")
     end if

     if (pSum) then
        call MatGetInfo(dFdx,MAT_GLOBAL_SUM,sumInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printSum("dFdx")
     end if

  end if

  if (pdFdw) then

     if (pLocal) then 
        call MatGetInfo(dFdw,MAT_LOCAL,localInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printLocal("dFdw")
     end if

     if (pMax) then 
        call MatGetInfo(dFdw,MAT_GLOBAL_MAX,maxInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printMax("dFdw")
     end if

     if (pSum) then
        call MatGetInfo(dFdw,MAT_GLOBAL_SUM,sumInfo,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call printSum("dFdw")
     end if

  end if

contains

  subroutine printLocal(matrixName)
    implicit none
    character*(*) :: matrixName
    ! Each Proc should print this info:

    print *,'================================================================='
900 format(A,A,A,I3)
901 format(A,I2,A,F20.10)
902 format(A,A)

    write(*,900) "Local Matrix Information for ",matrixName," on proc :",myid
    write(*,901) "[",myid,"] Block Size: ",localInfo(1)
    write(*,901) "[",myid,"] Non Zeros Allocated  : ",localInfo(2)
    write(*,901) "[",myid,"] Non Zeros Used       : ",localInfo(3)
    write(*,901) "[",myid,"] Non Zeros Unused     : ",localInfo(4)
    write(*,901) "[",myid,"] Memory (MBytes)       : ",localInfo(5)/1000000
    write(*,901) "[",myid,"] Number Mat Assemblies: ",localInfo(6)
    write(*,901) "[",myid,"] Mallocs              : ",localInfo(7)
    write(*,901) "[",myid,"] Fill Ratio Given for ILU   : ",localInfo(8)
    write(*,901) "[",myid,"] Fill Ratio Required for ILU: ",localInfo(9)
    write(*,901) "[",myid,"] Mallocs for Factor         : ",localInfo(10)
    print *,'================================================================='

  end subroutine printLocal

  subroutine printMax(matrixName)

    implicit none
    character*(*) :: matrixName

    if (myid == 0) then 
       print *,'================================================================='
       write(*,902) "Max Matrix Information for ",matrixName
       write(*,901) "[",myid,"] Block Size: ",maxInfo(1)
       write(*,901) "[",myid,"] Non Zeros Allocated  : ",maxInfo(2)
       write(*,901) "[",myid,"] Non Zeros Used       : ",maxInfo(3)
       write(*,901) "[",myid,"] Non Zeros Unused     : ",maxInfo(4)
       write(*,901) "[",myid,"] Memory (MBytes)       : ",maxInfo(5)/1000000
       write(*,901) "[",myid,"] Number Mat Assemblies: ",maxInfo(6)
       write(*,901) "[",myid,"] Mallocs              : ",maxInfo(7)
       write(*,901) "[",myid,"] Fill Ratio Given for ILU   : ",maxInfo(8)
       write(*,901) "[",myid,"] Fill Ratio Required for ILU: ",maxInfo(9)
       write(*,901) "[",myid,"] Mallocs for Factor         : ",maxInfo(10)
       print *,'================================================================='
    end if

900 format(A,A,A,I3)
901 format(A,I2,A,F20.10)
902 format(A,A)
  end subroutine printMax

  subroutine printSum(matrixName)

    implicit none
    character*(*) :: matrixName
    if (myid == 0) then 
       print *,'================================================================='
       write(*,902) "Summed Matrix Information for ",matrixName
       write(*,901) "[",myid,"] Block Size: ",sumInfo(1)
       write(*,901) "[",myid,"] Non Zeros Allocated  : ",sumInfo(2)
       write(*,901) "[",myid,"] Non Zeros Used       : ",sumInfo(3)
       write(*,901) "[",myid,"] Non Zeros Unused     : ",sumInfo(4)
       write(*,901) "[",myid,"] Memory (MBytes)       : ",sumInfo(5)/1000000
       write(*,901) "[",myid,"] Number Mat Assemblies: ",sumInfo(6)
       write(*,901) "[",myid,"] Mallocs              : ",sumInfo(7)
       write(*,901) "[",myid,"] Fill Ratio Given for ILU   : ",sumInfo(8)
       write(*,901) "[",myid,"] Fill Ratio Required for ILU: ",sumInfo(9)
       write(*,901) "[",myid,"] Mallocs for Factor         : ",sumInfo(10)
       print *,'================================================================='
    end if

900 format(A,A,A,I3)
901 format(A,I2,A,G25.18)
902 format(A,A)
  end subroutine printSum
#endif
end subroutine matrixInfo

