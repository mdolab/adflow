!
!      ******************************************************************
!      *                                                                *
!      * FILE:          warp.f90                                        *
!      * AUTHOR:        Juan J. Alonso                                  *
!      * STARTING DATE: 03-04-2004                                      *
!      * LAST MODIFIED: 04-23-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       SUBROUTINE WARP(NPTS,XYZ,BLOCKNUMS,IJK)
!
!      ******************************************************************
!      *                                                                *
!      * THIS IS THE MAIN DRIVER ROUTINE FOR THE MULTIBLOCK MESH        *
!      * PERTURBATION ALGORITHM BASED ON WARP-MB.  THIS ROUTINE         *
!      * RECEIVES THE SURFACE MESH PERTURBATIONS AS A NUMBER NPTS OF    *
!      * SURFACE POINTS, THEIR DISPLACED (X,Y,Z) COORDINATES, AND THEIR *
!      * ID GIVEN BY (BLOCK,I,J,K) AND CORRESPONDING TO THE BLOCK       *
!      * NUMBER OF THE ORIGINAL MESH.  THE SUBROUTINE EXPECTS THAT BOTH *
!      * A MASTER CORNER AND A MASTER EDGE LIST HAVE ALREADY BEEN       *
!      * CREATED AND GOES THROUGH EACH OF THE FOLLOWING STEPS IF        *
!      * NECESSARY:                                                     *
!      *                                                                *
!      * 1) CHECK TO MAKE SURE THAT MASTER CORNER AND EDGE LISTS HAVE   *
!      *    ALREADY BEEN CREATED.  IF THEY HAVE NOT, EXIT WITH A        *
!      *    MESSAGE STATING THAT THE MASTER CORNER AND EDGE LISTS NEED  *
!      *    TO BE CREATED.                                              *
!      * 2) GO THROUGH THE LIST OF SURFACE NODES TO BE MOVED AND        *
!      *    IDENTIFY WHICH MASTER CORNERS AND MASTER EDGES ARE BEING    *
!      *    EXPLICITLY PERTURBED.                                       *
!      * 3) STORE THE SURFACE DISPLACEMENTS OF THE DISPLACED SURFACE    *
!      *    NODES IN THE APPROPRIATE DATA STRUCTURES.                   *
!      * 4) TRANSFER EXPLICITLY PERTURBED CORNERS AND EDGES TO ALL OF   *
!      *    THE CORNERS AND EDGES THAT BELONG TO THE SAME MASTER CORNER *
!      *    AND EDGES.                                                  *
!      * 5) FLAG ALL EDGES THAT TOUCH AN EXPLICITLY PERTURBED CORNER    *
!      *    AS IMPLICITLY PERTURBED EDGES.                              *
!      * 6) RUN WARPBLK ON EVERY BLOCK IN THE MESH THAT HAS AT LEAST    *
!      *    ONE EXPLICITLY OR IMPLICITLY PERTURBED EDGE OR AN           *
!      *    EXPLICITLY PERTURBED CORNER.                                *  
!      * 7) STORE THE PERTURBED BLOCKS FOR LATER TRANSFER TO THE        *
!      *    RECEIVING PROGRAM (TYPICALLY TFLO2000).                     *
!      *                                                                *
!      ******************************************************************
!
       USE MESH_BLOCK 
       USE EDGES

       use complexify 
       IMPLICIT NONE

       INTERFACE
          SUBROUTINE FIND_EDGE_INDEX(BLOCKNUM,IJK,ME,E,EDGE_INDEX,ARR_SIZE)

            USE MESH_BLOCK
            USE EDGES

       use complexify 
            IMPLICIT NONE

            INTEGER(KIND=INT_TYPE) :: BLOCKNUM,EDGE_INDEX,ARR_SIZE
            TYPE(MASTER_EDGE), POINTER :: ME
            TYPE(EDGE), POINTER :: E
            INTEGER(KIND=INT_TYPE), DIMENSION(3) :: IJK

          END SUBROUTINE FIND_EDGE_INDEX

          SUBROUTINE FIND_MASTER_CORNER(BLOCKNUM,INDEX1,MCORNER)

            USE MESH_BLOCK
            USE CORNERS  

       use complexify 
            IMPLICIT NONE

            INTEGER(KIND=INT_TYPE) :: BLOCKNUM
            INTEGER(KIND=INT_TYPE), DIMENSION(3) :: INDEX1
            TYPE(MASTER_CORNER), POINTER :: MCORNER

          END SUBROUTINE FIND_MASTER_CORNER

          SUBROUTINE PERTURB_EDGE(E)

            USE MESH_BLOCK
            USE EDGES
            
       use complexify 
            IMPLICIT NONE
          
            TYPE(EDGE), POINTER :: E
          END SUBROUTINE PERTURB_EDGE

       END INTERFACE

!
!      SUBROUTINE ARGUMENTS.
!
       INTEGER(KIND= INT_TYPE), INTENT(IN) :: NPTS
       INTEGER(KIND= INT_TYPE), DIMENSION(NPTS), INTENT(IN) :: BLOCKNUMS
       INTEGER(KIND= INT_TYPE), DIMENSION(3,NPTS), INTENT(IN) :: IJK
       complex   (KIND=REAL_TYPE), DIMENSION(3,NPTS), INTENT(IN) :: XYZ
!
!      LOCAL VARIABLES.
!
       TYPE(CORNER),        POINTER :: C
       TYPE(MASTER_CORNER), POINTER :: MC
       TYPE(EDGE),          POINTER :: E
       TYPE(MASTER_EDGE),   POINTER :: ME

       INTEGER(KIND=INT_TYPE) :: N,BLOCK_NUM,EDGE_INDEX,ARR_SIZE,COUNT
       INTEGER(KIND=INT_TYPE) :: IMAX, JMAX, KMAX, IP3D, I, J, K, L
       INTEGER(KIND=INT_TYPE), DIMENSION(3) :: IJK_NUM
       INTEGER(KIND=INT_TYPE), DIMENSION( 6,NMESHBLOCKS) :: IFACEPTB
       INTEGER(KIND=INT_TYPE), DIMENSION(12,NMESHBLOCKS) :: IEDGEPTB

       complex(KIND=REAL_TYPE), ALLOCATABLE, DIMENSION(:,:,:,:) :: XYZNEW,XYZ0,S0
       complex(KIND=REAL_TYPE), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: DFACEI
       complex(KIND=REAL_TYPE), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: DFACEJ
       complex(KIND=REAL_TYPE), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: DFACEK

       LOGICAL :: CORNER_POINT,EDGE_POINT,FACE_POINT
!
!      FUNCTION DEFINITION.
!
       LOGICAL :: IS_CORNER,ON_EDGE,ON_FACE

!
!      ******************************************************************
!      *                                                                *
!      * INITIALIZE ALL NECESSARY DATA FOR MASTER LISTS AND BLOCK       *
!      *  PERTURBATIONS                                                 *
!      *                                                                *
!      ******************************************************************
!
       NULLIFY(C)
       NULLIFY(MC)
       NULLIFY(E)
       NULLIFY(ME)

       IJK_NUM  = 0
       IFACEPTB = 0
       IEDGEPTB = 0
!       
!      ******************************************************************
!      *                                                                *
!      * STEP 1: CHECK FOR THE EXISTENCE OF THE MASTER CORNER AND EDGE  *
!      *         LISTS AND EXIT IF THEY ARE NOT PRESENT.                *
!      *                                                                *
!      ******************************************************************
!
!      MASTER CORNER LIST
!
       IF (NCORNERS .ceq. 0) THEN
         WRITE(*,*)
         WRITE(*,*)'ERROR!!!'
         WRITE(*,*)'MASTER CORNER LIST HAS NOT YET BEEN INITIALIZED'
         WRITE(*,*)
         WRITE(*,*)'EXITING...'
         WRITE(*,*)
         STOP 'WARP: NO MASTER CORNER LIST FOR THIS MESH'
       END IF
!
!      MASTER EDGE LIST
!
       IF (NEDGES .ceq. 0) THEN
         WRITE(*,*)
         WRITE(*,*)'ERROR!!!'
         WRITE(*,*)'MASTER EDGE LIST HAS NOT YET BEEN INITIALIZED'
         WRITE(*,*)
         WRITE(*,*)'EXITING...'
         WRITE(*,*)
         STOP 'WARP: NO MASTER EDGE LIST FOR THIS MESH'
       END IF
!
!      ******************************************************************
!      *                                                                *
!      * STEP 2: IDENTIFY MASTER CORNERS AND EDGES THAT ARE BEING       *
!      *         EXPLICITLY PERTURBED.                                  *
!      *                                                                *
!      * STEP 3: STORE SURFACE NODE DISPLACEMENTS IN THE BLOCKS DATA    *
!      *         STRUCTURE FOR THE GIVEN BLOCK                          *
!      *                                                                *
!      ******************************************************************

!
!      INITIALIZE EXP_PERTURBED FLAGS IN BOTH MASTER CORNER AND EDGE
!      LISTS TO .FALSE.
!
       MC => FIRST_MC
       DO
         IF (ASSOCIATED(MC)) THEN
           MC % EXP_PERTURBED = .FALSE.
           MC => MC % NEXT
         ELSE
           EXIT
         END IF
       END DO

       ME => FIRST_MEDG
       DO
         IF (ASSOCIATED(ME)) THEN
           ME % EXP_PERTURBED = .FALSE.
           ME => ME % NEXT
         ELSE
           EXIT
         END IF
       END DO

!
!      FLAG MASTER CORNERS AND EDGES THAT ARE BEING EXPLICITLY PERTURBED
!
       DO N=1,NPTS

          BLOCK_NUM    = BLOCKNUMS(N)
          IJK_NUM(1)   = IJK(1,N)
          IJK_NUM(2)   = IJK(2,N)
          IJK_NUM(3)   = IJK(3,N)

          CORNER_POINT = IS_CORNER(BLOCK_NUM,IJK_NUM)
          EDGE_POINT   = ON_EDGE  (BLOCK_NUM,IJK_NUM)
          FACE_POINT   = ON_FACE  (BLOCK_NUM,IJK_NUM)

          ! GO THROUGH EVERY EXPLICITLY PERTURBED CORNER AND RECORD ITS DXYZ
          ! VALUES, NOTING HOW MANY CORNERS ASSOCIATED WITH EACH MASTER CORNER
          ! HAVE BEEN EXPLICITLY PERTURBED.  IN THIS LOOP, THE DXYZ VALUES ARE
          ! JUST ADDED TO THE DXYZ STORED IN THE MC.  THESE DXYZ VALUES ARE
          ! AVERAGED ONCE ALL OF THE DXYZ VALUES HAVE BEEN RECORDED (IN CASE
          ! THERE ARE DIFFERENT PERTURBATIONS FOR CORNERS ASSOCIATED WITH
          ! THE SAME MASTER CORNER: AN AVERAGE PERTURBATION WILL RESULT).  

          IF (CORNER_POINT) THEN
             
             CALL FIND_MASTER_CORNER(BLOCK_NUM,IJK_NUM,MC)
             
             IF (.NOT. MC % EXP_PERTURBED) THEN
                MC % EXP_PERTURBED = .TRUE. 
                MC % DXYZ = 0.0
                MC % EXP_COUNT = 0
             END IF

             MC % DXYZ(1) = MC % DXYZ(1) +(XYZ(1,N) -MESHBLOCKS(BLOCK_NUM) % &
                            X0(IJK_NUM(1),IJK_NUM(2),IJK_NUM(3),1))
             MC % DXYZ(2) = MC % DXYZ(2) +(XYZ(2,N) -MESHBLOCKS(BLOCK_NUM) % &
                            X0(IJK_NUM(1),IJK_NUM(2),IJK_NUM(3),2))
             MC % DXYZ(3) = MC % DXYZ(3) +(XYZ(3,N) -MESHBLOCKS(BLOCK_NUM) % &
                            X0(IJK_NUM(1),IJK_NUM(2),IJK_NUM(3),3))
             MC % EXP_COUNT = MC % EXP_COUNT + 1

          END IF
          
          ! SEE IF THE POINT LIES ON AN EDGE, AND IF SO, RECORD THE DXYZ
          ! VALUES FOR THE NON-CORNER POINTS ON THAT EDGE IN THE APPROPRIATE
          ! ARRAY.  THE DXYZ VALUES FOR EACH INDIVIDUAL EDGE WILL BE AVERAGED
          ! FOR THE MASTER EDGE LATER IN THE SUBROUTINE.

          IF (EDGE_POINT) THEN

             CALL FIND_EDGE_INDEX(BLOCK_NUM,IJK_NUM,ME,E,EDGE_INDEX,ARR_SIZE)

             IF ((ARR_SIZE .cne. 0) .AND. ASSOCIATED(E)) THEN
                IF (.NOT. ASSOCIATED(ME % DXYZ)) ALLOCATE(ME%DXYZ(3,ARR_SIZE))
                IF (.NOT. ASSOCIATED(E % DXYZ)) THEN
                   ALLOCATE(E % DXYZ(3,ARR_SIZE))
                   E % DXYZ = 0.0
                END IF
                E % DXYZ(1,EDGE_INDEX) = XYZ(1,N) - MESHBLOCKS(BLOCK_NUM) % &
                     X0(IJK_NUM(1),IJK_NUM(2),IJK_NUM(3),1)
                E % DXYZ(2,EDGE_INDEX) = XYZ(2,N) - MESHBLOCKS(BLOCK_NUM) % &
                     X0(IJK_NUM(1),IJK_NUM(2),IJK_NUM(3),2)
                E % DXYZ(3,EDGE_INDEX) = XYZ(3,N) - MESHBLOCKS(BLOCK_NUM) % &
                     X0(IJK_NUM(1),IJK_NUM(2),IJK_NUM(3),3)
                E % EXP_PERTURBED = .TRUE.
                ME % EXP_PERTURBED = .TRUE.

             END IF
          END IF
         
          ! IF THE POINT IS NOT A CORNER POINT OR AN EDGE POINT, 
          ! DETERMINE WHICH FACE THE PERTURBED POINT IS ON, AND MARK THIS FACE
          ! AS BEING EXPLICITLY PERTURBED, IN ADDITION TO EXPLICITLY PETURBING
          ! THE X VALUES IN BLOCKS          
          IF ((.NOT. CORNER_POINT) .AND.                                     &
              (.NOT. EDGE_POINT)   .AND. FACE_POINT) THEN
             MESHBLOCKS(BLOCK_NUM) % X(IJK_NUM(1),IJK_NUM(2),IJK_NUM(3),1) = XYZ(1,N)
             MESHBLOCKS(BLOCK_NUM) % X(IJK_NUM(1),IJK_NUM(2),IJK_NUM(3),2) = XYZ(2,N)
             MESHBLOCKS(BLOCK_NUM) % X(IJK_NUM(1),IJK_NUM(2),IJK_NUM(3),3) = XYZ(3,N)
             IF (IJK_NUM(1) .ceq. 1) THEN
                IFACEPTB(1,BLOCK_NUM) = 2
             ELSEIF (IJK_NUM(1) .ceq. MESHBLOCKS(BLOCK_NUM) % IL) THEN
                IFACEPTB(2,BLOCK_NUM) = 2
             ELSEIF (IJK_NUM(2) .ceq. 1) THEN 
                IFACEPTB(3,BLOCK_NUM) = 2
             ELSEIF (IJK_NUM(2) .ceq. MESHBLOCKS(BLOCK_NUM) % JL) THEN
                IFACEPTB(4,BLOCK_NUM) = 2
             ELSEIF (IJK_NUM(3) .ceq. 1) THEN
                IFACEPTB(5,BLOCK_NUM) = 2
             ELSEIF (IJK_NUM(3) .ceq. MESHBLOCKS(BLOCK_NUM) % KL) THEN
                IFACEPTB(6,BLOCK_NUM) = 2
             END IF
          END IF
          
          IF (.NOT. FACE_POINT) THEN
             WRITE(*,*) 'WARNING - ILLEGAL PERTURBATION!  PERTURBATION PASSED TO '
             WRITE(*,*) 'WARP IS IN THE INTERIOR OF THE BLOCK, RATHER THAN ON A'
             WRITE(*,*) 'BLOCK FACE AS IT SHOULD BE - WARP WILL LIKELY FAIL'
             STOP 'ILLEGAL BLOCK PERTURBATION OF A INTERIOR POINT IN A BLOCK'
          END IF
       END DO

       ! AVERAGE THE DXYZ VALUES FOR EACH MASTER CORNER
       MC => FIRST_MC
       DO WHILE(ASSOCIATED(MC))
          IF (MC % EXP_PERTURBED) THEN
             IF (MC % EXP_COUNT .ceq. 0) then
                WRITE(*,*) 'WARNING - THE EXP_PERTURBED FLAG FOR THIS MASTER CORNER'
                WRITE(*,*) 'IS MARKED AS TRUE, BUT THE NUMBER OF EXPLICITLY'
                WRITE(*,*) 'PERTURBED CORNERS ASSOCIATED WITH THIS MASTER CORNER'
                WRITE(*,*) 'IS ZERO!'  
                STOP
             END IF
             MC % DXYZ = MC % DXYZ / REAL(MC % EXP_COUNT)
          END IF
          MC => MC % NEXT
       END DO

       ! AVERAGE THE DXYZ VALUES FOR EACH MASTER EDGE
       ME => FIRST_MEDG
       DO WHILE(ASSOCIATED(ME))
          IF (ME % EXP_PERTURBED) THEN
             COUNT     = 0
             ME % DXYZ = 0.0
             E => ME % FIRST_EDGE
             DO WHILE(ASSOCIATED(E))
                IF (E % EXP_PERTURBED) THEN
                   COUNT = COUNT + 1
                   IF (E % DIRECTION .ceq. 1) THEN
                      ME % DXYZ = ME % DXYZ + E % DXYZ
                   ELSE
                      ! IF DIRECTION IS NOT 1, THEN READ THE DXYZ VALUES IN BACKWARDS
                      DO N=1,SIZE(ME % DXYZ,2)
                         ME%DXYZ(:,N) = ME%DXYZ(:,N) +E%DXYZ(:,SIZE(E % DXYZ,2)+1-N)
                      END DO
                   END IF
                END IF
                E => E % NEXT
             END DO
             ME % DXYZ = ME % DXYZ / REAL(COUNT)
          END IF
          ME => ME % NEXT
       END DO

!
!      ******************************************************************
!      *                                                                *
!      * STEP 4: TRANSFER THE DISPLACEMENTS OF ALL EXPLICITLY PERTURBED *
!      *         CORNERS AND EDGES TO ALL THE MEMBERS OF THE            *
!      *         CORRESPONDING MASTER CORNER AND EDGE.                  *
!      *                                                                *
!      ******************************************************************

!
!      FIND EXPLICITLY PERTURBED MASTER CORNERS AND TRANSFER THEIR
!      DISPLACEMENTS TO ALL CORNERS ASSOCIATED WITH THAT MASTER CORNER
!
       MC => FIRST_MC
       DO WHILE((ASSOCIATED(MC)))
          IF (MC % EXP_PERTURBED) THEN
             ! LOOP THROUGH ALL EXPLICITLY PERTURBED MASTER CORNERS AND
             ! PROPAGATE THE AVERAGE DXYZ VALUES FOR THAT MASTER CORNER
             ! TO EVERY CORNER ASSOCIATED WITH THAT MASTER CORNER
             C => MC % FIRST_CORNER
             DO WHILE(ASSOCIATED(C))
                MESHBLOCKS(C % BLOCK) % X (C % IJK(1),C % IJK(2),C % IJK(3),1) = &
                MESHBLOCKS(C % BLOCK) % X0(C % IJK(1),C % IJK(2),C % IJK(3),1) + &
                     MC % DXYZ(1)
                MESHBLOCKS(C % BLOCK) % X (C % IJK(1),C % IJK(2),C % IJK(3),2) = &
                MESHBLOCKS(C % BLOCK) % X0(C % IJK(1),C % IJK(2),C % IJK(3),2) + &
                     MC % DXYZ(2)
                MESHBLOCKS(C % BLOCK) % X (C % IJK(1),C % IJK(2),C % IJK(3),3) = &
                MESHBLOCKS(C % BLOCK) % X0(C % IJK(1),C % IJK(2),C % IJK(3),3) + &
                     MC % DXYZ(3)
                C => C % NEXT
             END DO              
          END IF
          MC => MC % NEXT
       END DO
!
!      DO THE SAME FOR EXPLICITLY PERTURBED MASTER EDGES
!
       ME => FIRST_MEDG
       DO WHILE(ASSOCIATED(ME)) 
          IF (ME % EXP_PERTURBED) THEN
             E => ME % FIRST_EDGE
             DO WHILE(ASSOCIATED(E))
                CALL PERTURB_EDGE(E)
                E => E % NEXT
             END DO
          END IF
          ME => ME % NEXT
       END DO

!
!      ******************************************************************
!      *                                                                *
!      * STEPS 5 & 6: FLAG ALL IMPLICITLY PERTURBED EDGES AND FACES AND *
!      * THEN RUN WARPBLK IN ALL BLOCKS THAT NEED TO BE PERTURBED       *
!      *                                                                *
!      ******************************************************************
!
       
       ! FLAG ALL FACES AND EDGES IN THE MESH THAT ARE IMPLICITLY PERTURBED
       DO N=1,NMESHBLOCKS
          CALL FLAG_IMPLICITS(N, IFACEPTB(1,N), IEDGEPTB(1,N))
!          CALL FLAG_IMPLICITS(N, IFACEPTB(:,N), IEDGEPTB(:,N))
       END DO
       
       ! LOOP THROUGH ALL BLOCKS AND CALL WARPBLK WHERE APPROPRIATE
       DO N=1,NMESHBLOCKS

          ! SAVE NEW AND INITIAL XYZ VALUES TO BE PASSED TO WARPBLK
          IMAX = MESHBLOCKS(N) % IL
          JMAX = MESHBLOCKS(N) % JL
          KMAX = MESHBLOCKS(N) % KL

          ALLOCATE(XYZ0(3,0:IMAX+1,0:JMAX+1,0:KMAX+1),XYZNEW(3,0:IMAX+1,0:JMAX+1,0:KMAX+1))
          xyz0 = 0
          xyznew = 0

          XYZ0(1,1:IMAX,1:JMAX,1:KMAX) =  MESHBLOCKS(N) % X0(1:IMAX,1:JMAX,1:KMAX,1)
          XYZ0(2,1:IMAX,1:JMAX,1:KMAX) = MESHBLOCKS(N) % X0(1:IMAX,1:JMAX,1:KMAX,2)
          XYZ0(3,1:IMAX,1:JMAX,1:KMAX) = MESHBLOCKS(N) % X0(1:IMAX,1:JMAX,1:KMAX,3)
          XYZNEW(1,1:IMAX,1:JMAX,1:KMAX) = MESHBLOCKS(N) % X(1:IMAX,1:JMAX,1:KMAX,1)
          XYZNEW(2,1:IMAX,1:JMAX,1:KMAX) = MESHBLOCKS(N) % X(1:IMAX,1:JMAX,1:KMAX,2)
          XYZNEW(3,1:IMAX,1:JMAX,1:KMAX) = MESHBLOCKS(N) % X(1:IMAX,1:JMAX,1:KMAX,3)
          ! SAVE DFACE VALUES TO BE PASSED TO WARPBLK
!          ALLOCATE(DFACEI(3,0:JMAX+1,0:KMAX+1,2,4),DFACEJ(3,0:IMAX+1,0:KMAX+1,2,4),&
!               DFACEK(3,0:IMAX+1,0:JMAX+1,2,4))
          ALLOCATE(DFACEI(3,1:JMAX,1:KMAX,2,4),DFACEJ(3,1:IMAX,1:KMAX,2,4),&
               DFACEK(3,1:IMAX,1:JMAX,2,4))

          DFACEI = 0.0
          DFACEJ = 0.0
          DFACEK = 0.0
          DFACEI(1,1:JMAX,1:KMAX,1,1) = MESHBLOCKS(N) % X(1,1:JMAX,1:KMAX,1)-&
               MESHBLOCKS(N) % X0(1,1:JMAX,1:KMAX,1)
          DFACEI(2,1:JMAX,1:KMAX,1,1) = MESHBLOCKS(N) % X(1,1:JMAX,1:KMAX,2)-&
               MESHBLOCKS(N) % X0(1,1:JMAX,1:KMAX,2)
          DFACEI(3,1:JMAX,1:KMAX,1,1) = MESHBLOCKS(N) % X(1,1:JMAX,1:KMAX,3)-&
               MESHBLOCKS(N) % X0(1,1:JMAX,1:KMAX,3)
          DFACEI(1,1:JMAX,1:KMAX,2,1) = MESHBLOCKS(N) % X(IMAX,1:JMAX,1:KMAX,1)-&
               MESHBLOCKS(N) % X0(IMAX,1:JMAX,1:KMAX,1)
          DFACEI(2,1:JMAX,1:KMAX,2,1) = MESHBLOCKS(N) % X(IMAX,1:JMAX,1:KMAX,2)-&
               MESHBLOCKS(N) % X0(IMAX,1:JMAX,1:KMAX,2)
          DFACEI(3,1:JMAX,1:KMAX,2,1) = MESHBLOCKS(N) % X(IMAX,1:JMAX,1:KMAX,3)-&
               MESHBLOCKS(N) % X0(IMAX,1:JMAX,1:KMAX,3)
          DFACEJ(1,1:IMAX,1:KMAX,1,1) = MESHBLOCKS(N) % X(1:IMAX,1,1:KMAX,1)-&
               MESHBLOCKS(N) % X0(1:IMAX,1,1:KMAX,1)
          DFACEJ(2,1:IMAX,1:KMAX,1,1) = MESHBLOCKS(N) % X(1:IMAX,1,1:KMAX,2)-&
               MESHBLOCKS(N) % X0(1:IMAX,1,1:KMAX,2)
          DFACEJ(3,1:IMAX,1:KMAX,1,1) = MESHBLOCKS(N) % X(1:IMAX,1,1:KMAX,3)-&
               MESHBLOCKS(N) % X0(1:IMAX,1,1:KMAX,3)
          DFACEJ(1,1:IMAX,1:KMAX,2,1) = MESHBLOCKS(N) % X(1:IMAX,JMAX,1:KMAX,1)-&
               MESHBLOCKS(N) % X0(1:IMAX,JMAX,1:KMAX,1)
          DFACEJ(2,1:IMAX,1:KMAX,2,1) = MESHBLOCKS(N) % X(1:IMAX,JMAX,1:KMAX,2)-&
               MESHBLOCKS(N) % X0(1:IMAX,JMAX,1:KMAX,2)
          DFACEJ(3,1:IMAX,1:KMAX,2,1) = MESHBLOCKS(N) % X(1:IMAX,JMAX,1:KMAX,3)-&
               MESHBLOCKS(N) % X0(1:IMAX,JMAX,1:KMAX,3)
          DFACEK(1,1:IMAX,1:JMAX,1,1) = MESHBLOCKS(N) % X(1:IMAX,1:JMAX,1,1)-&
               MESHBLOCKS(N) % X0(1:IMAX,1:JMAX,1,1)
          DFACEK(2,1:IMAX,1:JMAX,1,1) = MESHBLOCKS(N) % X(1:IMAX,1:JMAX,1,2)-&
               MESHBLOCKS(N) % X0(1:IMAX,1:JMAX,1,2)
          DFACEK(3,1:IMAX,1:JMAX,1,1) = MESHBLOCKS(N) % X(1:IMAX,1:JMAX,1,3)-&
               MESHBLOCKS(N) % X0(1:IMAX,1:JMAX,1,3)
          DFACEK(1,1:IMAX,1:JMAX,2,1) = MESHBLOCKS(N) % X(1:IMAX,1:JMAX,KMAX,1)-&
               MESHBLOCKS(N) % X0(1:IMAX,1:JMAX,KMAX,1)
          DFACEK(2,1:IMAX,1:JMAX,2,1) = MESHBLOCKS(N) % X(1:IMAX,1:JMAX,KMAX,2)-&
               MESHBLOCKS(N) % X0(1:IMAX,1:JMAX,KMAX,2)
          DFACEK(3,1:IMAX,1:JMAX,2,1) = MESHBLOCKS(N) % X(1:IMAX,1:JMAX,KMAX,3)-&
               MESHBLOCKS(N) % X0(1:IMAX,1:JMAX,KMAX,3)
          ALLOCATE(S0(3,0:IMAX+1,0:JMAX+1,0:KMAX+1))
          s0 = 0.0
          ! CALL WARPBLK FOR THE CURRENT BLOCK IF ANY OF THE FACES OR EDGE IN THAT 
          ! BLOCK ARE IMPLICITLY OR EXPLICITLY PERTURBED
          IF (MAXVAL(IFACEPTB(1:6,N)) >= 1 .OR. MAXVAL(IEDGEPTB(1:12,N)) >= 1) THEN
            
             CALL WARPBLK(IFACEPTB(1:6,N),IEDGEPTB(1:12,N),-3,0,1,IMAX,JMAX,&
                  KMAX,XYZ0,S0,DFACEI,DFACEJ,DFACEK,XYZNEW)

!            ASSIGN THESE NEW XYZ VALUES TO THE MESH ITSELF
             DO I=1,IMAX
                DO J=1,JMAX
                   DO K=1,KMAX
                      MESHBLOCKS(N) % X(I,J,K,1) = XYZNEW(1,I,J,K)
                      MESHBLOCKS(N) % X(I,J,K,2) = XYZNEW(2,I,J,K)
                      MESHBLOCKS(N) % X(I,J,K,3) = XYZNEW(3,I,J,K)
                   END DO
                END DO
             END DO
          END IF
          DEALLOCATE(XYZNEW,XYZ0,S0,DFACEI,DFACEJ,DFACEK)
       END DO
       END SUBROUTINE WARP


       LOGICAL FUNCTION IS_CORNER(BLOCK_NUM, IJK)
!
!      ******************************************************************
!      *                                                                *
!      * DETERMINES WHETHER OR NOT THE POINT ASSOCIATED WITH THE GIVEN  *
!      * BLOCK NUMBER AND I,J,K VALUES IS A CORNER OF THAT BLOCK        *
!      *                                                                *
!      ******************************************************************
!         
       USE MESH_BLOCK  

       use complexify 
       IMPLICIT NONE
!
!      LOCAL VARIABLES.
!       
       INTEGER(KIND=INT_TYPE) :: BLOCK_NUM, IMAX, JMAX, KMAX
       INTEGER(KIND=INT_TYPE), DIMENSION(3) :: IJK
!
!      ******************************************************************
!      *                                                                *
!      * BEGIN EXECUTION                                                *
!      *                                                                *
!      ******************************************************************
!      
       IS_CORNER = .FALSE.

       IMAX = MESHBLOCKS(BLOCK_NUM) % IL
       JMAX = MESHBLOCKS(BLOCK_NUM) % JL
       KMAX = MESHBLOCKS(BLOCK_NUM) % KL
       
       ! CHECKS THE EIGHT POSSIBILITIES TO SEE IF THE POINT IS A CORNER
       IF ((IJK(1) .ceq. 1) .AND. (IJK(2) .ceq. 1) .AND. (IJK(3) .ceq. 1)) THEN 
          IS_CORNER = .TRUE.
       ELSEIF ((IJK(1) .ceq. IMAX) .AND. (IJK(2) .ceq. 1) .AND. (IJK(3) .ceq. 1)) THEN 
          IS_CORNER = .TRUE.
       ELSEIF ((IJK(1) .ceq. 1) .AND. (IJK(2) .ceq. JMAX) .AND. (IJK(3) .ceq. 1)) THEN 
          IS_CORNER = .TRUE.
       ELSEIF ((IJK(1) .ceq. IMAX) .AND. (IJK(2) .ceq. JMAX) .AND. (IJK(3) .ceq. 1)) THEN 
          IS_CORNER = .TRUE.
       ELSEIF ((IJK(1) .ceq. 1) .AND. (IJK(2) .ceq. 1) .AND. (IJK(3) .ceq. KMAX)) THEN 
          IS_CORNER = .TRUE.
       ELSEIF ((IJK(1) .ceq. IMAX) .AND. (IJK(2) .ceq. 1) .AND. (IJK(3) .ceq. KMAX)) THEN 
          IS_CORNER = .TRUE.
       ELSEIF ((IJK(1) .ceq. 1) .AND. (IJK(2) .ceq. JMAX) .AND. (IJK(3) .ceq. KMAX)) THEN 
          IS_CORNER = .TRUE.
       ELSEIF ((IJK(1) .ceq. IMAX) .AND. (IJK(2) .ceq. JMAX) .AND. (IJK(3) .ceq. KMAX)) THEN 
          IS_CORNER = .TRUE.
       END IF
          
       END FUNCTION IS_CORNER

       LOGICAL FUNCTION ON_FACE(BLOCK_NUM, IJK)
!
!      ******************************************************************
!      *                                                                *
!      * DETERMINES WHETHER OR NOT THE POINT ASSOCIATED WITH THE GIVEN  *
!      * BLOCK NUMBER AND I,J,K VALUES IS ON A FACE OF THAT BLOCK       *
!      *                                                                *
!      ******************************************************************
!         
       USE MESH_BLOCK  

       use complexify 
       IMPLICIT NONE
!
!      LOCAL VARIABLES.
!       
       INTEGER(KIND=INT_TYPE) :: BLOCK_NUM, IMAX, JMAX, KMAX
       INTEGER(KIND=INT_TYPE), DIMENSION(3) :: IJK
!
!      ******************************************************************
!      *                                                                *
!      * BEGIN EXECUTION                                                *
!      *                                                                *
!      ******************************************************************
!      
       ON_FACE = .FALSE.

       IMAX = MESHBLOCKS(BLOCK_NUM) % IL
       JMAX = MESHBLOCKS(BLOCK_NUM) % JL
       KMAX = MESHBLOCKS(BLOCK_NUM) % KL

       IF (((IJK(1) .ceq. 1)) .OR. ((IJK(1) .ceq. IMAX)) .OR.                     &
           ((IJK(2) .ceq. 1)) .OR. ((IJK(2) .ceq. JMAX)) .OR.                     &
           ((IJK(3) .ceq. 1)) .OR. ((IJK(3) .ceq. KMAX))) THEN
          ON_FACE = .TRUE.
       END IF
          
       END FUNCTION ON_FACE


       LOGICAL FUNCTION ON_EDGE(BLOCK_NUM, IJK)
!
!      ******************************************************************
!      *                                                                *
!      * DETERMINES WHETHER OR NOT THE POINT ASSOCIATED WITH THE GIVEN  *
!      * BLOCK NUMBER AND I,J,K VALUES IS ON AN EDGE OF THAT BLOCK      *
!      *                                                                *
!      ******************************************************************
!         
       USE MESH_BLOCK  

       use complexify 
       IMPLICIT NONE
!
!      LOCAL VARIABLES.
!       
       INTEGER(KIND=INT_TYPE) :: BLOCK_NUM, IMAX, JMAX, KMAX
       INTEGER(KIND=INT_TYPE), DIMENSION(3) :: IJK
!
!      ******************************************************************
!      *                                                                *
!      * BEGIN EXECUTION                                                *
!      *                                                                *
!      ******************************************************************
!      
       ON_EDGE = .FALSE.

       IMAX = MESHBLOCKS(BLOCK_NUM) % IL
       JMAX = MESHBLOCKS(BLOCK_NUM) % JL
       KMAX = MESHBLOCKS(BLOCK_NUM) % KL

       IF ((IJK(1) .ceq. 1) .AND. (IJK(2) .ceq. 1)) THEN 
          ON_EDGE = .TRUE.
       ELSEIF ((IJK(1) .ceq. IMAX) .AND. (IJK(2) .ceq. 1)) THEN 
          ON_EDGE = .TRUE.
       ELSEIF ((IJK(1) .ceq. 1) .AND. (IJK(2) .ceq. JMAX)) THEN 
          ON_EDGE = .TRUE.
       ELSEIF ((IJK(1) .ceq. IMAX) .AND. (IJK(2) .ceq. JMAX)) THEN 
          ON_EDGE = .TRUE.
       ELSEIF ((IJK(1) .ceq. 1) .AND. (IJK(3) .ceq. 1)) THEN 
          ON_EDGE = .TRUE.
       ELSEIF ((IJK(1) .ceq. IMAX) .AND. (IJK(3) .ceq. 1)) THEN 
          ON_EDGE = .TRUE.
       ELSEIF ((IJK(1) .ceq. 1) .AND. (IJK(3) .ceq. KMAX)) THEN 
          ON_EDGE = .TRUE.
       ELSEIF ((IJK(1) .ceq. IMAX) .AND. (IJK(3) .ceq. KMAX)) THEN 
          ON_EDGE = .TRUE.
       ELSEIF ((IJK(2) .ceq. 1) .AND. (IJK(3) .ceq. 1)) THEN 
          ON_EDGE = .TRUE.
       ELSEIF ((IJK(2) .ceq. JMAX) .AND. (IJK(3) .ceq. 1)) THEN 
          ON_EDGE = .TRUE.
       ELSEIF ((IJK(2) .ceq. 1) .AND. (IJK(3) .ceq. KMAX)) THEN 
          ON_EDGE = .TRUE.
       ELSEIF ((IJK(2) .ceq. JMAX) .AND. (IJK(3) .ceq. KMAX)) THEN 
          ON_EDGE = .TRUE.   
       END IF
          
       END FUNCTION ON_EDGE

       SUBROUTINE FIND_EDGE_INDEX(BLOCKNUM,IJK,ME,E,EDGE_INDEX,ARR_SIZE)
!
!      ******************************************************************
!      *                                                                *
!      * RETURNS THE EDGE INDEX OF THE POINT ASSOCIATED WITH THE        *
!      * GIVEN BLOCKNUM AND IJK VALUES AS WELL AS A POINTER TO THE      *
!      * MASTER EDGE, ME, AND A POINTER TO THE ACTUAL EDGE, E,          *
!      * ASSOCIATED WITH THIS POINT.                                    *
!      *                                                                *
!      * THE SUBROUTINE ALSO RETURNS THE VALUE ARR_SIZE WHICH IS THE    *
!      * SIZE OF THE DXYZ ARRAY NEEDED TO STORE THE DXYZ VALUES FOR     *
!      * MASTER EDGE ME.  EDGE_INDEX REFERS TO THE ORDER NUMBER OF THE  *
!      * GIVEN POINT ALONG THE ASSOCIATED EDGE, ASSUMING DIRECTION = 1. *
!      *                                                                *
!      * THIS SUBROUTINE ONLY USES NON-CORNER POINTS (CORNER            *
!      * PERTURBATIONS ARE HANDLED ELSEWHERE). AN INDEX NUMBER OF 1     *
!      * CORRESPONDS TO THE FIRST POINT AFTER THE CORNER ALONG A GIVEN  *
!      * EDGE. FOR EXAMPLE, ON AN EDGE WHERE THE I INDEX CHANGES, THE   *
!      * EDGE_INDEX CAN VARY FROM 1 TO IL-2                             *
!      *                                                                *
!      ******************************************************************
!        
       USE MESH_BLOCK
       USE EDGES

       use complexify 
       IMPLICIT NONE
       
       INTERFACE

          SUBROUTINE FIND_MASTER_CORNER(BLOCKNUM, INDEX1, MCORNER)
            USE MESH_BLOCK
            USE CORNERS  
            
       use complexify 
            IMPLICIT NONE

            INTEGER(KIND=INT_TYPE) :: BLOCKNUM
            INTEGER(KIND=INT_TYPE), DIMENSION(3) :: INDEX1
            TYPE(MASTER_CORNER), POINTER :: MCORNER

          END SUBROUTINE FIND_MASTER_CORNER
          
          SUBROUTINE FIND_MASTER_EDGE(MC1,MC2,MEDG,EXIST_MEDG)

            USE MESH_BLOCK
            USE CORNERS
            USE EDGES

       use complexify 
            IMPLICIT NONE

            TYPE(MASTER_CORNER), POINTER :: MC1, MC2
            TYPE(MASTER_EDGE), POINTER :: MEDG
            LOGICAL :: EXIST_MEDG
          END SUBROUTINE FIND_MASTER_EDGE

       END INTERFACE
!
!      SUBROUTINE ARGUMENTS
!
       INTEGER(KIND=INT_TYPE) :: BLOCKNUM, EDGE_INDEX, ARR_SIZE
       TYPE(MASTER_EDGE), POINTER :: ME
       TYPE(EDGE), POINTER :: E
       INTEGER(KIND=INT_TYPE), DIMENSION(3) :: IJK
!
!      LOCAL VARIABLES.
!       
       TYPE(MASTER_CORNER), POINTER :: MC1, MC2
       INTEGER(KIND=INT_TYPE) :: IMAX, JMAX, KMAX, INDEX_NUM
       INTEGER(KIND=INT_TYPE), DIMENSION(3) :: INDEX1, INDEX2
       LOGICAL :: ME_EXISTS
       TYPE(EDGE), POINTER :: ETEMP
!
!      ******************************************************************
!      *                                                                *
!      * BEGIN EXECUTION                                                *
!      *                                                                *
!      ******************************************************************
!
       IMAX   = MESHBLOCKS(BLOCKNUM) % IL
       JMAX   = MESHBLOCKS(BLOCKNUM) % JL
       KMAX   = MESHBLOCKS(BLOCKNUM) % KL

       INDEX1 = IJK
       INDEX2 = IJK
       
       ! DETERMINE WHICH INDEX ON THIS EDGE IS THE CHANGING INDEX, AND RECORD
       ! THE NUMBER OF THIS INDEX IN INDEX_NUM.  ALSO, IJK VALUES ARE GENERATED
       ! FOR THE TWO CORNERS CORRESPONDING TO THIS EDGE

       IF (IJK(1) < IMAX .AND. IJK(1) > 1) THEN
          INDEX_NUM  = 1
          EDGE_INDEX = IJK(1)-1
          INDEX1(1)  = 1
          INDEX2(1)  = IMAX
          ARR_SIZE   = IMAX-2
       ELSEIF (IJK(2) < JMAX .AND. IJK(2) > 1) THEN
          INDEX_NUM  = 2
          EDGE_INDEX = IJK(2)-1
          INDEX1(2)  = 1
          INDEX2(2)  = JMAX
          ARR_SIZE   = JMAX-2
       ELSEIF (IJK(3) < KMAX .AND. IJK(3) > 1) THEN
          INDEX_NUM  = 3
          EDGE_INDEX = IJK(3)-1
          INDEX1(3)  = 1
          INDEX2(3)  = KMAX
          ARR_SIZE   = KMAX-2
       ELSE
          ! ARR_SIZE OF ZERO SIGNALS THAT THE POINT IS A CORNER AND ITS
          ! PERTURBATION WILL NOT BE STORED IN THE MASTER EDGE DXYZ ARRAY
          ARR_SIZE = 0
          RETURN
       END IF
       
       ! FIND THE MASTER EDGE ASSOCIATED WITH THIS POINT SO THAT WE CAN 
       ! DETERMINE THE ORIENTATION OF THE MASTER CORNERS (TO BE USED LATER
       ! BY CHECKING THE DIRECTION FLAGS OF EACH EDGE)

       CALL FIND_MASTER_CORNER(BLOCKNUM, INDEX1, MC1)
       CALL FIND_MASTER_CORNER(BLOCKNUM, INDEX2, MC2)

       CALL FIND_MASTER_EDGE(MC1,MC2,ME,ME_EXISTS)

       IF (.NOT. ME_EXISTS) THEN
          CALL FIND_MASTER_EDGE(MC2,MC1,ME,ME_EXISTS)
       END IF

       ! QUICK SANITY CHECK IN CASE THE POINT BEING PASSED TO THIS ROUTINE
       ! CANNOT BE FOUND IN THE MASTER EDGE LIST

       IF (.NOT. ME_EXISTS) THEN
          WRITE(*,*) 'EDGE POINT NOT ON ANY MASTER EDGE'
          WRITE(*,*) 'SOURCE CODE PROBLEM IN FIND_EDGE_INDEX'
          STOP 'FIND_EDGE_INDEX: SOURCE CODE PROBLEM'
       END IF
 
       ! FIND THE ACTUAL EGDE ASSOCIATED WITH THIS POINT
       NULLIFY(E)
       ETEMP => ME % FIRST_EDGE

       DO WHILE(ASSOCIATED(ETEMP))
          IF (((ETEMP % BLOCK1 .ceq. BLOCKNUM) .AND.                         &
               (SUM(ABS(INDEX1 - ETEMP % IJK1)) .ceq. 0)  .AND.              &
               (ETEMP % BLOCK2 .ceq. BLOCKNUM) .AND.                         &
               (SUM(ABS(INDEX2 - ETEMP % IJK2)) .ceq. 0)) .OR.               &
              ((ETEMP % BLOCK2 .ceq. BLOCKNUM) .AND.                         &
               (SUM(ABS(INDEX2 - ETEMP % IJK1)) .ceq. 0)  .AND.              &
               (ETEMP % BLOCK1 .ceq. BLOCKNUM) .AND.                         &
               (SUM(ABS(INDEX1 - ETEMP % IJK2)) .ceq. 0))) THEN
             E => ETEMP
             EXIT
          END IF
          ETEMP => ETEMP % NEXT
       END DO
       IF (.NOT. ASSOCIATED(E)) THEN
          WRITE(*,*) 'ERROR IN FIND_EDGE_INDEX- EDGE NOT FOUND'
          RETURN
       END IF

       ! IF THE EDGE NUMBERING FOR THE ACTUAL EDGE BEING EXAMINED STARTS AT 
       ! THE MAX FOR THAT INDEX AND GOES TO 1, CHANGE EDGE_INDEX TO REFLECT
       ! THIS SITUATION.
       IF (E % IJK1(INDEX_NUM) > EDGE_INDEX) THEN
          EDGE_INDEX = (E % IJK1(INDEX_NUM) - 1) - EDGE_INDEX
       END IF
       
       END SUBROUTINE FIND_EDGE_INDEX

       SUBROUTINE PERTURB_EDGE(E)
!
!      ******************************************************************
!      *                                                                * 
!      * CHANGES THE X, Y, AND Z VALUES OF THE PASSED EDGE BY ADDING    *
!      * THE DXYZ VALUES FROM THE CORRESPONDING MASTER EDGE.  THE       *
!      * ROUTINE CHECKS THE EDGE TO DETERMINE THE DIRECTION AND ASSIGNS *
!      * THE DXYZ VALUES IN THE CORRECT ORDER.                          *
!      *                                                                *
!      ******************************************************************
!        
       USE MESH_BLOCK
       USE EDGES

       use complexify 
       IMPLICIT NONE
!
!      SUBROUTINE ARGUMENTS
!
       TYPE(EDGE), POINTER :: E
!
!      LOCAL VARIABLES.
!       
       INTEGER(KIND=INT_TYPE) :: IMAX, JMAX, KMAX, INDEX_NUM, I, N
!
!      ******************************************************************
!      *                                                                *
!      * BEGIN EXECUTION                                                *
!      *                                                                *
!      ******************************************************************
!  
      ! print *,'in perturb edge'
       IMAX = MESHBLOCKS(E % BLOCK1) % IL
       JMAX = MESHBLOCKS(E % BLOCK1) % JL
       KMAX = MESHBLOCKS(E % BLOCK1) % KL

       ! DETERMINE THE CHANGING INDEX FOR THE CURRENT EDGE
       IF (E % IJK1(1) .cne. E % IJK2(1)) THEN
          INDEX_NUM = 1
       ELSE IF (E % IJK1(2) .cne. E % IJK2(2)) THEN
          INDEX_NUM = 2
       ELSE 
          INDEX_NUM = 3
       END IF

       ! LOOP THROUGH THE ELEMENTS OF THE DXYZ ARRAY IN THE MASTER EDGE 
       ! ASSOCIATED WITH THE CURRENT EDGE AND ADD THESE VALUES TO THE 
       ! INITIAL XYZ VALUES (X0)
       DO I=1,SIZE(E % MASTER % DXYZ,2)
          
          ! N IS THE NUMBER ADDED TO THE CHANGING INDEX OF THE FIRST CORNER
          ! FOR THIS EDGE.  N DEPENDS ON THE DIRECTION OF THE EDGE AS WELL AS
          ! THE NUMBERING SCHEME FOR THE EDGE (IF THE CHANGING INDEX FOR THIS
          ! EDGE BEGINS AT ITS MAX VALUE AND GOES TO 1, THIS WILL AFFECT THE
          ! VALUE OF N)
          IF ((E % DIRECTION .ceq. 1) .AND. (E % IJK1(INDEX_NUM) .ceq. 1)) THEN
             N = I
          ELSE IF ((E % DIRECTION .ceq. 1) .AND. (E % IJK1(INDEX_NUM) .cne. 1)) THEN
             N = -I
          ELSE IF ((E % DIRECTION .ceq. -1) .AND. (E % IJK1(INDEX_NUM) .ceq. 1)) THEN
             N = E % IJK2(INDEX_NUM) -1 -I
          ELSE IF ((E % DIRECTION .ceq. -1) .AND. (E % IJK1(INDEX_NUM) .cne. 1)) THEN
             N = I + 1 -E % IJK1(INDEX_NUM)
          ELSE 
             WRITE(*,*) 'WEIRD EDGE FOUND. STOPPING...'
          END IF
           
          IF (INDEX_NUM .ceq. 1) THEN
             MESHBLOCKS(E % BLOCK1) % X(E % IJK1(1)+N,E % IJK1(2),E % IJK1(3),1) = &
                  MESHBLOCKS(E % BLOCK1) % X0(E % IJK1(1)+N,E % IJK1(2),E % IJK1(3),1) &
                  + E % MASTER % DXYZ(1,I)
             MESHBLOCKS(E % BLOCK1) % X(E % IJK1(1)+N,E % IJK1(2),E % IJK1(3),2) = &
                  MESHBLOCKS(E % BLOCK1) % X0(E % IJK1(1)+N,E % IJK1(2),E % IJK1(3),2) &
                  + E % MASTER % DXYZ(2,I)
             MESHBLOCKS(E % BLOCK1) % X(E % IJK1(1)+N,E % IJK1(2),E % IJK1(3),3) = &
                  MESHBLOCKS(E % BLOCK1) % X0(E % IJK1(1)+N,E % IJK1(2),E % IJK1(3),3) &
                  + E % MASTER % DXYZ(3,I)
          ELSEIF (INDEX_NUM .ceq. 2) THEN
             MESHBLOCKS(E % BLOCK1) % X(E % IJK1(1),E % IJK1(2)+N,E % IJK1(3),1) = &
                  MESHBLOCKS(E % BLOCK1) % X0(E % IJK1(1),E % IJK1(2)+N,E % IJK1(3),1) &
                  + E % MASTER % DXYZ(1,I)
             MESHBLOCKS(E % BLOCK1) % X(E % IJK1(1),E % IJK1(2)+N,E % IJK1(3),2) = &
                  MESHBLOCKS(E % BLOCK1) % X0(E % IJK1(1),E % IJK1(2)+N,E % IJK1(3),2) &
                  + E % MASTER % DXYZ(2,I)
             MESHBLOCKS(E % BLOCK1) % X(E % IJK1(1),E % IJK1(2)+N,E % IJK1(3),3) = &
                  MESHBLOCKS(E % BLOCK1) % X0(E % IJK1(1),E % IJK1(2)+N,E % IJK1(3),3) &
                  + E % MASTER % DXYZ(3,I)
          ELSE
             MESHBLOCKS(E % BLOCK1) % X(E % IJK1(1),E % IJK1(2),E % IJK1(3)+N,1) = &
                  MESHBLOCKS(E % BLOCK1) % X0(E % IJK1(1),E % IJK1(2),E % IJK1(3)+N,1) &
                  + E % MASTER % DXYZ(1,I)
             MESHBLOCKS(E % BLOCK1) % X(E % IJK1(1),E % IJK1(2),E % IJK1(3)+N,2) = &
                  MESHBLOCKS(E % BLOCK1) % X0(E % IJK1(1),E % IJK1(2),E % IJK1(3)+N,2) &
                  + E % MASTER % DXYZ(2,I)
             MESHBLOCKS(E % BLOCK1) % X(E % IJK1(1),E % IJK1(2),E % IJK1(3)+N,3) = &
                  MESHBLOCKS(E % BLOCK1) % X0(E % IJK1(1),E % IJK1(2),E % IJK1(3)+N,3) &
                  + E % MASTER % DXYZ(3,I)
          END IF
       END DO
              
       END SUBROUTINE PERTURB_EDGE
