MODULE PARTDYNAMICS
  ! subroutines to simulate transport and dynamics of 1d particles
  ! VERSION BY: SAURABH MOGRE 
  ! DATE: MAY 2020
  
  USE MT19937, ONLY: GRND,RNORM,GRND
  
  IMPLICIT NONE
  
  CONTAINS
  
  SUBROUTINE TRANSPORTSIM(GROUPLIST,PARAMP,NTRIALS,DOFUSION)
    !simulates the transport of phagosomes and lysosomes in a 1D domain
    USE PARTICLEUTILS, ONLY: PARTICLE,PARTICLEGROUP,PARTICLEPARAMS,&
    OUTPUTPROT,OUTPUTSNAPSHOT,REMOVEPARTICLE,&
    WRITENTOUCH
    USE KEYS, ONLY: MAXNPROT,PRINTEVERY,SNAPFILE,SNAPEVERY,OUTFILE,OUTEVERY,DELT,&
    NSTEPS,CALCNTOUCH,NTOUCHFILE,NSTATE
    USE FUSIONUTILS, ONLY: FUSEPARTICLES,UPDATEACIDIFICATION
    
    INTEGER, INTENT(IN) :: NTRIALS
    TYPE(PARTICLEGROUP), TARGET :: GROUPLIST(NTRIALS)
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    TYPE(PARTICLE), POINTER :: PLIST(:)
    TYPE(PARTICLE), POINTER :: PARTP
    TYPE(PARTICLEPARAMS), POINTER :: PARAMP
    
    LOGICAL, INTENT(IN) :: DOFUSION
    
    DOUBLE PRECISION, ALLOCATABLE :: INFO(:)
    INTEGER :: TC,PC,TRC,NPROC,TOTNPART,SURVCT,SC
    LOGICAL :: APPENDTOFILE
    DOUBLE PRECISION :: CELLPROT(MAXNPROT),AXPROT(MAXNPROT),TOTAXPROT(MAXNPROT)
    DOUBLE PRECISION :: TOTRATE,TMP,RATEVEC(NSTATE)
    
    IF(PRINTEVERY.NE.0) THEN
      PRINT'(A5,I1,A4,I8)', 'Step ',0,' of ',NSTEPS
    END IF
    
    ! number of lysosomes that reach the end of the domain without fusing
    SURVCT = 0
    
    ! output initial snapshot
    APPENDTOFILE = .FALSE.
    IF(SNAPEVERY.GT.0) THEN
      ALLOCATE(INFO(3))
      INFO = (/0D0,DFLOAT(NTRIALS),PARAMP%DOMLEN/)
      PRINT*, 'SAVING INITIAL SNAPSHOT'
      CALL OUTPUTSNAPSHOT(GROUPLIST,SNAPFILE,INFO,NTRIALS,APPENDTOFILE)
    END IF
    
    ! output initial protein content
    IF(OUTEVERY.GT.0) THEN
      CELLPROT = 0D0
      AXPROT = 0D0
      TOTAXPROT = 0D0
      CALL OUTPUTPROT(0D0,CELLPROT,TOTAXPROT,OUTFILE,MAXNPROT,NTRIALS,APPENDTOFILE)
    END IF
    
    ! output initial number of contacts
    ALLOCATE(PLIST(PARAMP%MAXNPART))
    PARTP=>PLIST(1)
    IF(CALCNTOUCH) CALL WRITENTOUCH(PARTP,NTOUCHFILE,APPENDTOFILE,NTRIALS)
    
    APPENDTOFILE = .TRUE.
    DO TC = 1,NSTEPS
      
      ! print information on simulation status
      IF((PRINTEVERY.NE.0)) THEN
        IF(MOD(TC,PRINTEVERY).EQ.0) THEN
          PRINT'(A5,I8,A4,I8)', 'Step ',TC,' of ',NSTEPS
        END IF
      END IF
      
      PARAMP%CURTIME = TC*DELT
      
      ! initialize arrays to record protein content
      TOTAXPROT = 0D0
      
      ! go through all particle groups
      DO TRC = 1,NTRIALS
        PGROUP=>GROUPLIST(TRC)
        
        ! check for fusion
        IF(DOFUSION) THEN
          CALL FUSEPARTICLES(PGROUP,PARAMP)
        END IF
        
        ! move particles, change direction, and delete out of bounds
        PC = 1
        NPROC = 0
        TOTNPART = SUM(PGROUP%NPART)
        DO WHILE(NPROC.LT.TOTNPART.AND.SUM(PGROUP%NPART).GT.0)
          PARTP=>PGROUP%PLIST(PC)
          NPROC = NPROC + 1
          
          ! if particle is in the tip region, check if it exits
          IF(PARTP%ISINTIP) THEN
            IF(PARTP%KEXIT.LT.0D0.OR.GRND().LE.1-EXP(-PARTP%KEXIT*DELT)) THEN
              PARTP%ISINTIP = .FALSE.
              PARTP%CURDIR = -1
              PGROUP%NTIP(PARTP%PTYPE) = MAX(PGROUP%NTIP(PARTP%PTYPE)-1,0) 
              IF(PARTP%PTYPE.EQ.1) PARTP%CURSTATE = 2
              PARTP%POS = 1D0
            ELSE IF(PARTP%KEXIT.EQ.0D0) THEN
              CALL REMOVEPARTICLE(PGROUP,PC)
              CYCLE
            ELSE
              PC=PC+1
              CYCLE
            END IF
          END IF
          
          ! move particle
          PARTP%PREVPOS = PARTP%POS
          IF(PARTP%CURDIR.EQ.0) THEN
            PARTP%VEL = 0
          ELSE
            PARTP%VEL = PARAMP%VEL(PARTP%PTYPE,(PARTP%CURDIR+3)/2)
            PARTP%POS = PARTP%POS+PARTP%CURDIR*PARTP%VEL*DELT
          END IF
          
          ! check if particle changes direction in its current state
          IF(PARTP%CURDIR.EQ.0) THEN ! stationary particle
            IF(GRND().LE.1-EXP(-PARTP%KWALK(PARTP%CURSTATE)*DELT)) THEN !start a walk
              IF(GRND().LE.PARTP%PRETRO(PARTP%CURSTATE)) THEN !choose direction
                PARTP%CURDIR = -1
              ELSE
                PARTP%CURDIR = 1
              END IF
            END IF
          ELSE 
            !PRINT*, 'TESTX, CURSTATE, CURDIR, KHALT, KWALK, VEL: ',&
            !& PARTP%CURSTATE,PARTP%CURDIR, PARTP%KHALT(PARTP%CURSTATE),&
            !& PARTP%KHALT(PARTP%CURSTATE),PARTP%VEL
            IF(GRND().LE.1-EXP(-PARTP%KHALT(PARTP%CURSTATE)*DELT)) THEN ! stop walking
              PARTP%CURDIR = 0
            END IF
          END IF
          
          ! check if particle changes state
          TOTRATE = SUM(PARTP%KSWITCH(PARTP%CURSTATE,:))
          IF(TOTRATE.GT.0D0) THEN
            RATEVEC = PARTP%KSWITCH(PARTP%CURSTATE,:)/TOTRATE
            IF(GRND().LE.1-EXP(-TOTRATE*DELT)) THEN
              TMP = GRND()
              DO SC = 1,NSTATE
                IF(SC.EQ.PARTP%CURSTATE) CYCLE
                IF(TMP.LE.SUM(RATEVEC(1:SC))) THEN
                  PARTP%CURSTATE = SC
                  PARTP%CURDIR = PARAMP%STARTDIR(PARTP%CURSTATE)
                  PARTP%TSWITCH = PARAMP%CURTIME-PARTP%TMAKE
                  EXIT
                END IF
              END DO
            END IF
          END IF
          
          ! check if fused particle has matured
          IF(PARTP%CANFUSE.AND.PARTP%NFUSE.GT.0) THEN
            IF(GRND().LE.(1-EXP(-PARTP%KMAT*DELT))) THEN
              ! PRINT*,'TESTX: PARTICLE HAS MATURED, KMAT:',PARTP%KMAT,' NFUSE:',PARTP%NFUSE
              PARTP%CANFUSE = .FALSE.
            END IF
          END IF
          
          ! update time for which degradation has occured
          IF(PARTP%NFUSE.GT.0) THEN
            PARTP%TDEG = PARTP%TDEG+DELT
          END IF
          
          
          ! protein conversion
          ! ---------------------------------------
          ! TODO: change to matrix exponentiation?
          ! ---------------------------------------
          IF(PARTP%NPROT.GT.0) THEN
            PARTP%PROT = PARTP%PROT+MATMUL(PARTP%PROTCONV,PARTP%PROT)*DELT
          END IF
          
          ! update acidification level of particle
          CALL UPDATEACIDIFICATION(PARTP,PARAMP)

          ! deal with particles that move out of bounds depending on boundary conditions
          IF(PARTP%POS.LT.0D0) THEN
            IF(PARAMP%REFBOUND(1)) THEN
              PARTP%POS = -PARTP%POS
              PARTP%CURDIR = -PARTP%CURDIR
            ELSE
              ! record total protein if output is requested
              IF(PARTP%NPROT.GT.0.AND.OUTEVERY.GT.0) THEN
                CELLPROT(1:PARTP%NPROT) = CELLPROT(1:PARTP%NPROT)+PARTP%PROT
              END IF
              
              ! record a touch if required
              IF(CALCNTOUCH) CALL WRITENTOUCH(PARTP,NTOUCHFILE,APPENDTOFILE)
              
              ! remove particle from simulation
              CALL REMOVEPARTICLE(PGROUP,PC)
              CYCLE
            END IF
          END IF
          
          IF(PARTP%POS.GT.PARAMP%DOMLEN) THEN
            IF(PARAMP%REFBOUND(2)) THEN
              PARTP%POS = 2*PARAMP%DOMLEN-PARTP%POS
              PARTP%CURDIR = -PARTP%CURDIR
              ! PRINT*, 'PARTICLE REVERSED DIRECTION'
            ELSE
              ! particle is now in the tip region
              PARTP%ISINTIP = .TRUE.
              PARTP%POS = 1D0
              PGROUP%NTIP(PARTP%PTYPE) = PGROUP%NTIP(PARTP%PTYPE)+1
              PC = PC+1
              CYCLE
            END IF
          END IF
          
          PC = PC+1
        END DO
        
        ! measure protein content in the axon
        IF(OUTEVERY.GT.0) THEN
          CALL GETAXPROT(AXPROT,PGROUP,MAXNPROT)
          TOTAXPROT(1:MAXNPROT) = TOTAXPROT(1:MAXNPROT)+AXPROT
        END IF
        
        ! particle production
        CALL PRODUCEPARTICLES(PGROUP,PARAMP)      
        
      END DO
      
      ! output snapshot to file
      IF(SNAPEVERY.GT.0) THEN
        IF(MOD(TC,SNAPEVERY).EQ.0) THEN
          INFO(1) = TC*DELT
          CALL OUTPUTSNAPSHOT(GROUPLIST,SNAPFILE,INFO,NTRIALS,APPENDTOFILE)
        END IF
      END IF
      
      ! output protein flux to file
      IF(OUTEVERY.GT.0) THEN
        IF(MOD(TC,OUTEVERY).EQ.0) THEN
          CALL OUTPUTPROT(TC*DELT,CELLPROT,TOTAXPROT,OUTFILE,MAXNPROT,NTRIALS,APPENDTOFILE)
        END IF
      END IF
      
      
    END DO
    
    PRINT*, 'DYNAMICS DONE'
    
  END SUBROUTINE TRANSPORTSIM
  
  SUBROUTINE PRODUCEPARTICLES(PGROUP,PARAMP)
    ! subroutine to produce new particles
    USE PARTICLEUTILS, ONLY: PARTICLEGROUP,PARTICLEPARAMS,PARTICLE,MAKEPARTICLE
    USE GENUTIL,ONLY: GETPOISSONEVENTS
    USE KEYS, ONLY: DELT
    IMPLICIT NONE
    
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    TYPE(PARTICLEPARAMS), POINTER :: PARAMP
    
    TYPE(PARTICLE),POINTER :: PARTP
    
    INTEGER :: NC,PID,NMAKE,PC
    
    NMAKE = 0
    DO NC = 1,PARAMP%NTYPE ! go through each type of particle
      
      ! obtain number of particles produced
      CALL GETPOISSONEVENTS(PARAMP%KPROD(NC)*DELT,NMAKE)
      
      IF(NMAKE.GT.0) THEN
        DO PC = 1,NMAKE
          IF(SUM(PGROUP%NPART).GE.PGROUP%MAXNPART) THEN
            PRINT*, 'ERROR: too many particles added.',&
            PGROUP%NPART, PGROUP%MAXNPART
            STOP 1
          END IF
          
          PARTP=>PGROUP%PLIST(SUM(PGROUP%NPART)+1)
          PID = PGROUP%LASTID+1
          CALL MAKEPARTICLE(PARTP,PID,NC,PARAMP,CANFUSE=.TRUE.)
          PGROUP%NPART(NC)  = PGROUP%NPART(NC) + 1
          PGROUP%LASTID = PID  
          ! IF(PARTP%PTYPE.EQ.2) THEN
          !   PRINT*, 'AV PRODUCED AT: ', PARTP%POS
          ! END IF
        END DO
      END IF
    END DO
    
  END SUBROUTINE PRODUCEPARTICLES
  
  SUBROUTINE GETAXPROT(AXPROT,PGROUP,MAXNPROT)
    ! subroutine to measure amount of protein in the axon
    USE PARTICLEUTILS, ONLY: PARTICLEGROUP,PARTICLE
    IMPLICIT NONE
    DOUBLE PRECISION :: AXPROT(MAXNPROT)
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    TYPE(PARTICLE), POINTER :: PARTP
    INTEGER, INTENT(IN) :: MAXNPROT
    
    INTEGER :: PC
    
    AXPROT = 0D0
    DO PC = 1,SUM(PGROUP%NPART)
      PARTP=>PGROUP%PLIST(PC)
      IF(PARTP%NPROT.EQ.0) CYCLE
      AXPROT(1:PARTP%NPROT) = AXPROT(1:PARTP%NPROT) + PARTP%PROT
    END DO
    
  END SUBROUTINE GETAXPROT
  
  
END MODULE PARTDYNAMICS
