MODULE FUSIONUTILS
  ! subroutines to simulate fusion events between particles
  ! VERSION BY: SAURABH MOGRE 
  ! DATE: MAY 2020
  
  CONTAINS
  SUBROUTINE FUSEPARTICLES(PGROUP,PARAMP)
    ! subroutine to implement fusion of particles
    USE PARTICLEUTILS, ONLY: PARTICLEGROUP,PARTICLE,PARTICLEPARAMS,&
    REMOVEPARTICLE
    
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    TYPE(PARTICLEPARAMS), POINTER :: PARAMP
    TYPE(PARTICLE), POINTER :: PARTP1,PARTP2,PARTP
    
    INTEGER :: PC1,PC2,NPROC,TOTNPART,DELIND,KEEPIND,NCHECK,REMPART
    LOGICAL :: CHECKFUSE
    
    TOTNPART = SUM(PGROUP%NPART)
    PC1 = 1
    NPROC = 0
    DO WHILE(NPROC.LT.(TOTNPART-1).AND.SUM(PGROUP%NPART).GT.0)
      NPROC = NPROC + 1
      PARTP1=>PGROUP%PLIST(PC1)
      
      ! run through other particles
      NCHECK = 0
      PC2 = PC1+1
      REMPART = SUM(PGROUP%NPART)-PC1
      DO WHILE(NCHECK.LT.REMPART.AND.SUM(PGROUP%NPART).GT.0) 
        NCHECK = NCHECK+1
        PARTP2=>PGROUP%PLIST(PC2)
        
        ! check if fusion occurs between a given particle pair
        CALL CHECKFUSION(CHECKFUSE,PARTP1,PARTP2,PARAMP)
        
        ! if fusion has occured, update particle states
        IF(CHECKFUSE) THEN
          ! decide which particle to keep
          ! default: particle with larger PTYPE  keeps its identity
          IF(PARTP1%PTYPE.GT.PARTP2%PTYPE) THEN
            DELIND = PC2
            KEEPIND = PC1
          ELSE 
            DELIND = PC1
            KEEPIND = PC2
          END IF
          
          ! update properties of the fused particle
          PARTP=>PGROUP%PLIST(KEEPIND)
          PARTP%NFUSE = PARTP%NFUSE+1
          
          ! update protein conversion rate caused by acquiring material during fusion
          CALL UPDATECONVRATE(PARTP,PARAMP)
          
          ! if particle switches state on fusion, update the state
          IF(PARAMP%SWITCHONFUSE(PARTP%PTYPE)) THEN
            PARTP%CURSTATE = 1
            PARTP%CURDIR = -1
          END IF
          
          ! remove unused particle
          CALL REMOVEPARTICLE(PGROUP,DELIND)
          IF(DELIND.EQ.PC1) THEN
            PC1 = PC1-1
            EXIT
          ELSE
            NPROC = NPROC+1
            CYCLE
          END IF
        END IF
        
        PC2 = PC2 + 1
      END DO
      PC1 = PC1+1
    END DO
    
  END SUBROUTINE FUSEPARTICLES
  
  SUBROUTINE CHECKFUSION(CHECKFUSE,PARTP1,PARTP2,PARAMP)
    ! subroutine to check whether fusion occurs
    ! the checkfuse flag returns the result of the fusion process
    USE PARTICLEUTILS, ONLY: PARTICLE,PARTICLEPARAMS
    USE MT19937, ONLY: GRND
    TYPE(PARTICLE),POINTER :: PARTP1,PARTP2
    TYPE(PARTICLEPARAMS),POINTER :: PARAMP
    DOUBLE PRECISION :: POS1,POS2,PREVPOS1,PREVPOS2,VEC1,VEC2
    INTEGER :: PTYPE1,PTYPE2
    DOUBLE PRECISION :: PF
    LOGICAL,INTENT(OUT) :: CHECKFUSE
    LOGICAL :: AREINTERACTING
    
    CHECKFUSE = .FALSE.
    AREINTERACTING = .FALSE.
    
    
    ! end check if particles have reached maximum number of fusion events allowed
    IF((PARTP1%NFUSE.GE.PARAMP%MAXNFUSE(PARTP1%PTYPE)).OR.&
    (PARTP2%NFUSE.GE.PARAMP%MAXNFUSE(PARTP2%PTYPE))) RETURN
    
    ! get probability of fusion
    IF(PARAMP%USEKFUSE) THEN
      PF = 1-EXP(-PARAMP%KFUSE*PARAMP%DELT)
      ! end check if fusion probability is zero
      IF(PF.LE.0D0) RETURN
    ELSE
      PF = PARAMP%PFUSE(PTYPE1,PTYPE2)
      ! end check if fusion probability is zero
      IF(PF.LE.0D0) RETURN
      
      POS1 = PARTP1%POS
      PREVPOS1 = PARTP1%PREVPOS
      PTYPE1 = PARTP1%PTYPE
      POS2 = PARTP2%POS
      PREVPOS2 = PARTP2%PREVPOS
      PTYPE2 = PARTP2%PTYPE
      
      ! check relative particle positions
      VEC1 = PREVPOS1-PREVPOS2
      VEC2 = POS1-POS2
      IF(VEC1*VEC2.GT.0D0) RETURN
    END IF
        
    ! set flag if fusion is successful
    IF(GRND().LE.PF) CHECKFUSE = .TRUE.
    
    ! set error if fusion is unsuccessful with pf=1
    IF(.NOT.CHECKFUSE.AND.PF.EQ.1D0) THEN
      PRINT*, 'UNSUCCESSFUL FUSION BETWEEN ID1,ID2,TYPE1,TYPE2:',&
      PARTP1%ID,PARTP2%ID,PARTP1%PTYPE,PARTP2%PTYPE
      STOP 1                                    
    END IF
    
    ! update interaction
    ! CALL UPDATEINTERACTION(PARTP1,PARTP2,.NOT.CHECKFUSE)
    
  END SUBROUTINE CHECKFUSION
  
  SUBROUTINE CHECKINTERACTION(PARTP1,PARTP2,STATUS)
    ! checks whether particles are currently interacting
    USE PARTICLEUTILS,ONLY: PARTICLE
    IMPLICIT NONE
    
    TYPE(PARTICLE),POINTER :: PARTP1,PARTP2
    LOGICAL,INTENT(OUT) :: STATUS
    INTEGER :: PC
    LOGICAL :: INT12,INT21
    
    STATUS = .FALSE.
    INT12 = .FALSE.
    INT21 = .FALSE.
    
    PC = 1
    DO WHILE(PC.LE.MIN(PARTP1%NINT,SIZE(PARTP1%INTERACTING)))
      IF(PARTP1%INTERACTING(PC).EQ.PARTP2%ID) THEN
        INT12 = .TRUE.
        EXIT
      END IF
      PC = PC+1
    END DO
    
    PC = 1
    DO WHILE(PC.LE.MIN(PARTP2%NINT,SIZE(PARTP2%INTERACTING)))
      IF(PARTP2%INTERACTING(PC).EQ.PARTP1%ID) THEN
        INT21 = .TRUE.
        EXIT
      END IF
      PC = PC+1
    END DO
    
    
    IF(INT12.AND.INT21) THEN
      STATUS = .TRUE.
    ELSE IF(XOR(INT12,INT21)) THEN
      PRINT*, 'ERROR IN CHECKINTERACTION: check status'
      STOP 1
    END IF
    
  END SUBROUTINE CHECKINTERACTION
  
  SUBROUTINE UPDATEINTERACTION(PARTP1,PARTP2,STATUS)
    ! updates the interaction between particles to that indicated by STATUS
    USE PARTICLEUTILS,ONLY: PARTICLE
    IMPLICIT NONE
    
    TYPE(PARTICLE),POINTER :: PARTP1,PARTP2
    LOGICAL,INTENT(IN) :: STATUS
    INTEGER :: PC
    LOGICAL :: INTFLAG
    
    ! PRINT*, 'TESTX: UPDATING INTERACTION, STATUS, NINT1, NINT2:', STATUS,PARTP1%NINT,PARTP2%NINT
    
    IF(STATUS) THEN
      ! record new interaction
      IF(PARTP1%NINT.LT.SIZE(PARTP1%INTERACTING)) THEN
        PARTP1%INTERACTING(PARTP1%NINT+1) = PARTP2%ID
        PARTP1%NINT = PARTP1%NINT+1
      ELSE
        PRINT*, 'ERROR IN UPDATEINTERACTION: MAXIMUM NUMBER OF INTERACTIONS REACHED'
        PRINT*, 'POS, TYPE: ',PARTP1%POS,PARTP1%PTYPE
        PRINT*, 'NINT: ',PARTP1%NINT
        DO PC = 1,PARTP1%NINT
          PRINT*, PARTP1%INTERACTING(PC)
        END DO
        STOP 1
      END IF
      
      IF(PARTP2%NINT.LT.SIZE(PARTP2%INTERACTING)) THEN
        PARTP2%INTERACTING(PARTP2%NINT+1) = PARTP1%ID
        PARTP2%NINT = PARTP2%NINT+1
      ELSE
        PRINT*, 'ERROR IN UPDATEINTERACTION: MAXIMUM NUMBER OF INTERACTIONS REACHED'
        PRINT*, 'POS, TYPE: ',PARTP2%POS,PARTP2%PTYPE
        PRINT*, 'NINT: ',PARTP2%NINT
        DO PC = 1,PARTP2%NINT
          PRINT*, PARTP2%INTERACTING(PC)
        END DO
        STOP 1
      END IF
    ELSE
      ! end interaction
      PC = 1
      INTFLAG = .FALSE.
      DO WHILE(PC.LE.MIN(PARTP1%NINT,SIZE(PARTP1%INTERACTING)).AND.(.NOT.INTFLAG))
        IF(PARTP1%INTERACTING(PC).EQ.PARTP2%ID) THEN
          PARTP1%INTERACTING(PC:(PARTP1%NINT-1)) = PARTP1%INTERACTING((PC+1):PARTP1%NINT)
          PARTP1%INTERACTING(PARTP1%NINT) = 0
          PARTP1%NINT = MAX(PARTP1%NINT-1,0)
          INTFLAG = .TRUE.
          EXIT
        END IF
        PC = PC+1
      END DO
      
      PC = 1
      INTFLAG = .FALSE.
      DO WHILE(PC.LE.MIN(PARTP2%NINT,SIZE(PARTP2%INTERACTING)).AND.(.NOT.INTFLAG))
        IF(PARTP2%INTERACTING(PC).EQ.PARTP1%ID) THEN
          PARTP2%INTERACTING(PC:(PARTP2%NINT-1)) = PARTP2%INTERACTING((PC+1):PARTP2%NINT)
          PARTP2%INTERACTING(PARTP2%NINT) = 0
          PARTP2%NINT = MAX(PARTP2%NINT-1,0)
          INTFLAG = .TRUE.
          EXIT
        END IF
        PC = PC+1
      END DO
      
    END IF
    
  END SUBROUTINE UPDATEINTERACTION
  
  SUBROUTINE UPDATEACIDIFICATION(PARTP,PARAMP)
    ! update acidification of particle
    USE KEYS,ONLY: DELT
    USE PARTICLEUTILS,ONLY: PARTICLE,PARTICLEPARAMS
    IMPLICIT NONE
    
    TYPE(PARTICLE),POINTER :: PARTP
    TYPE(PARTICLEPARAMS),POINTER :: PARAMP
    
    IF(.NOT.(PARAMP%USEDILUTION)) THEN
      PARTP%ACIDITY = MIN(PARTP%ACIDITY+PARTP%ACIDIFICATIONRATE*DELT,PARAMP%STARTACID(1));
    END IF
    
  END SUBROUTINE UPDATEACIDIFICATION  
  
  SUBROUTINE UPDATECONVRATE(PARTP,PARAMP)
    ! update acidification rate/acidity and
    ! protein conversion rate for a particle, dependent on acidity
    USE PARTICLEUTILS,ONLY: PARTICLE,PARTICLEPARAMS
    IMPLICIT NONE
    
    TYPE(PARTICLE),POINTER :: PARTP
    TYPE(PARTICLEPARAMS),POINTER :: PARAMP
    
    IF(PARAMP%USEDILUTION) THEN
      ! volumetric dilution of acidity upon fusion
      PARTP%ACIDITY = PARAMP%STARTACID(1)*PARTP%NFUSE/&
      &(PARTP%NFUSE+PARAMP%PRAD(2)**3/PARAMP%PRAD(1)**3)
    ELSE
      ! fusion increases the number of pumps, thus increasing acidification rate
      ! acidification rate for each pump is proportional to the difference between
      ! the current acidification and the maximal lysosome acidification
      PARTP%ACIDIFICATIONRATE  = PARTP%NFUSE*PARAMP%ACIDIFICATIONRATE(PARTP%PTYPE)*&
      & (PARAMP%STARTACID(1)-PARTP%ACIDITY)
      ! PARTP%ACIDIFICATIONRATE  = PARAMP%ACIDIFICATIONRATE(PARTP%PTYPE)
    END IF
    
    IF(PARAMP%KM(PARTP%PTYPE).GT.0D0) THEN
      PARTP%PROTCONV = PARAMP%KMAX(:,:,PARTP%PTYPE)*PARTP%ACIDITY/&
      (PARTP%ACIDITY+PARAMP%KM(PARTP%PTYPE))
    ELSE
      PARTP%PROTCONV = PARAMP%KMAX(:,:,PARTP%PTYPE)
    END IF
    
  END SUBROUTINE UPDATECONVRATE  
  
END MODULE FUSIONUTILS
