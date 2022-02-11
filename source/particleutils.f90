MODULE PARTICLEUTILS
  ! utilities for defining parameters and states of individual particles,
  ! particle sets, and domain information

  TYPE PARTICLE
  ! object containing information about individual particles

    ! Particle ID
    INTEGER :: ID

    ! type of particle
    ! total number of particle types depends on system to be simulated
    INTEGER :: PTYPE

    ! particle radius
    DOUBLE PRECISION :: PRAD

    ! current particle position
    DOUBLE PRECISION :: POS

    ! time at which particle was made
    DOUBLE PRECISION :: TMAKE

    ! previous particle position
    DOUBLE PRECISION :: PREVPOS

    ! current particle direction
    ! default: +1 = moving forward; -1 = moving backward; 0 = stationary
    INTEGER :: CURDIR

    ! particle state
    ! can be used for particles that switch between states
    ! 1 = particles do not switch direction
    ! 2 = particles can stop and switch direction (bidirectional)
    ! 3 = particles are stationary
    INTEGER :: CURSTATE
    
    ! walk velocity in anterograde and retrograde direction respectively
    DOUBLE PRECISION :: VEL

    ! diffusion coefficient for diffusive particles
    DOUBLE PRECISION :: D

    ! rates to stop walking in each state 
    DOUBLE PRECISION, POINTER :: KHALT(:)

    ! rate to start walking in each state
    DOUBLE PRECISION, POINTER :: KWALK(:)

    ! is the particle located within the cell tip?
    LOGICAL :: ISINTIP

    ! rate to exit cell tip
    DOUBLE PRECISION :: KEXIT

    ! probability of moving retrograde while starting walk
    DOUBLE PRECISION, POINTER :: PRETRO(:)

    ! rate to switch states
    DOUBLE PRECISION, POINTER :: KSWITCH(:,:)

    ! rate to mature
    DOUBLE PRECISION :: KMAT

    ! acidity of the particle
    DOUBLE PRECISION :: ACIDITY

    ! rate of acidification increase
    ! accounts for vATPase pumps that import H+ ions to increase acidification with time
    ! fusion with lysosomes increases the number of available pumps
    DOUBLE PRECISION :: ACIDIFICATIONRATE

    ! can the particle fuse?
    LOGICAL :: CANFUSE

    ! list of particles currently interacting with this particle
    INTEGER,POINTER :: INTERACTING(:)

    ! current number of interactions
    INTEGER :: NINT

    ! time for which degradation has occured
    ! organelle starts degrading material after fusion
    DOUBLE PRECISION :: TDEG    

    ! location of first fusion
    DOUBLE PRECISION :: XFUSE

    ! time of first fusion
    DOUBLE PRECISION :: TFUSE

    ! time of first switch
    DOUBLE PRECISION :: TSWITCH

    ! number of fusion events
    INTEGER :: NFUSE

    ! number of fusion events needed to saturate decay rate
    INTEGER :: MAXNFUSE

    ! Number of different types of possible fusion events
    INTEGER :: NFUSETYPE

    ! Number of contacts
    INTEGER :: NTOUCH
    
    ! probability of fusion per contact
    DOUBLE PRECISION,POINTER :: PFUSE(:)

    ! fusion rate
    DOUBLE PRECISION,POINTER :: KFUSE(:)

    ! number of different protein types
    INTEGER :: NPROT
    
    ! array of protein content
    ! contains amounts of different protein types    
    DOUBLE PRECISION, POINTER :: PROT(:)

    ! conversion rates between protein types
    DOUBLE PRECISION, POINTER :: PROTCONV(:,:)

    ! have arrays been allocated?
    LOGICAL :: ARRAYSET = .FALSE.

  END TYPE PARTICLE

  TYPE PARTICLEGROUP
  ! object containing information about a group of mixed particle types

    ! number of active particles in the group
    INTEGER,POINTER :: NPART(:)

    ! maximum number of particles allowed
    INTEGER :: MAXNPART

    ! ID of the last added particle
    INTEGER :: LASTID

    ! pointer containing particles within the group
    TYPE(PARTICLE), POINTER :: PLIST(:)

    ! number of particles in the tip
    INTEGER, POINTER :: NTIP(:)
    
    ! arrays have been allocated
    LOGICAL :: ARRAYSET = .FALSE.      
  END TYPE PARTICLEGROUP


  TYPE PARTICLEPARAMS
  ! object containing the various parameters associated with each type of particle

    DOUBLE PRECISION :: DELT ! simulation timestep

    DOUBLE PRECISION :: CURTIME ! current time

    INTEGER :: NTYPE ! how many types of particles?

    INTEGER :: MAXNPART ! maximum number of particles allowed

    ! velocity and diffusivity for each type
    DOUBLE PRECISION, POINTER  :: VEL(:,:), D(:)

    ! production rate
    ! and reversal rate
    DOUBLE PRECISION, POINTER :: KPROD(:), KREV(:)

    ! maturation rate
    DOUBLE PRECISION, POINTER :: KMAT(:)

    ! rate to switch states
    DOUBLE PRECISION, POINTER :: KSWITCH(:,:,:)

    ! halting rate
    DOUBLE PRECISION, POINTER :: KHALT(:,:)

    ! walking rate
    DOUBLE PRECISION, POINTER :: KWALK(:,:)

    ! retrograde probability
    DOUBLE PRECISION, POINTER :: PRETRO(:,:)

    ! rate to exit tip region
    DOUBLE PRECISION, POINTER :: KEXIT(:)

    ! starting position of each particle type
    DOUBLE PRECISION, POINTER :: STARTPOS(:)

    ! start particle in tip?
    LOGICAL, POINTER :: STARTINTIP(:)

    ! distance over which particles spawn
    ! (decay length scale for exponetial process)
    DOUBLE PRECISION, POINTER :: SPAWNLENGTH(:)

    ! starting direction for each particle state
    INTEGER, POINTER :: STARTDIR(:)

    ! starting state for each particle type
    INTEGER, POINTER :: STARTSTATE(:)

    ! initial protein content of each particle type
    DOUBLE PRECISION, POINTER :: STARTPROT(:,:)

    ! initial acidity of each particle type
    DOUBLE PRECISION, POINTER :: STARTACID(:)

    ! switch particle state after fusing
    LOGICAL, POINTER :: SWITCHONFUSE(:)

    ! rate of acidification increase for each particle type
    DOUBLE PRECISION, POINTER :: ACIDIFICATIONRATE(:)

    ! particle radius
    DOUBLE PRECISION , POINTER :: PRAD(:)
    
    ! number of proteins for each type
    INTEGER, POINTER :: NPROT(:)

    ! maximum conversion rates btwn proteins for each type
    DOUBLE PRECISION, POINTER :: KMAX(:,:,:)

    ! rate of fusion and/or prob of fusion
    INTEGER, POINTER :: MAXNFUSE(:)
    DOUBLE PRECISION, POINTER :: PFUSE(:,:)
    DOUBLE PRECISION, POINTER :: KFUSE(:,:)
    DOUBLE PRECISION, POINTER :: KM(:) !michaelis-menten constant for protein degradation
    LOGICAL :: USEKFUSE
    LOGICAL :: USEINTERACTIONS

    ! use dilution to increase acidification
    LOGICAL :: USEDILUTION

    ! length of the domain
    DOUBLE PRECISION :: DOMLEN

    ! reflecting boundaries?
    ! status of the boundary at x=0 and x=DOMLEN
    ! default is [.FALSE.,.FALSE.]
    LOGICAL :: REFBOUND(2)

    ! arrays allocated?
    LOGICAL :: ARRAYSET = .FALSE.

  END TYPE PARTICLEPARAMS

CONTAINS

  SUBROUTINE SETUPPARTICLEGROUP(PGROUP,NPART,PLIST)
  ! set up a cluster of particles
    USE KEYS, ONLY: MAXNPART,NTYPE
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    TYPE(PARTICLE), INTENT(IN), POINTER :: PLIST(:)
    INTEGER, INTENT(IN) :: NPART(NTYPE)
    INTEGER :: TOTNPART

    TOTNPART = SUM(NPART)

    IF(PGROUP%ARRAYSET) THEN
      PRINT*, 'ERROR IN SETUPPARTICLEGROUP: trying to overwrite'
      STOP 1
    END IF

    ALLOCATE(PGROUP%NPART(NTYPE))
    PGROUP%NPART = NPART

    ALLOCATE(PGROUP%PLIST(MAXNPART))
    PGROUP%PLIST(1:TOTNPART) = PLIST(1:TOTNPART)
    PGROUP%MAXNPART = MAXNPART

    ALLOCATE(PGROUP%NTIP(NTYPE))
    PGROUP%NTIP = 0

    PGROUP%LASTID = TOTNPART
    
    PGROUP%ARRAYSET = .TRUE.
  END SUBROUTINE SETUPPARTICLEGROUP

  SUBROUTINE SETUPPARTICLEPARAMS(PARAMP)
  ! set up and initiate object to track particle parameters
    USE KEYS, ONLY : NSTATE,NTYPE,VEL,DIFFCONST,KMAT,KPROD,KHALT,KWALK,KSWITCH,&
         &PRETRO,KEXIT,STARTPOS,STARTINTIP,STARTDIR,&
         & PRAD,NPROT,PFUSE,KFUSE,KMAX,MAXNPROT,STARTPROT,MAXNPART,DOMLEN,REFBOUND,&
         MAXNFUSE,KMVALS,STARTACID,ACIDIFICATIONRATE,STARTSTATE,SWITCHONFUSE,USEDILUTION,&
         SPAWNLENGTH,USEKFUSE,USEINTERACTIONS,DELT

    
    IMPLICIT NONE
    TYPE(PARTICLEPARAMS), POINTER :: PARAMP
    
    PARAMP%NTYPE = NTYPE
    ALLOCATE(PARAMP%VEL(NTYPE,2),PARAMP%D(NTYPE),PARAMP%KPROD(NTYPE),&
         & PARAMP%KSWITCH(NSTATE,NSTATE,NTYPE),PARAMP%KMAT(NTYPE),&
         PARAMP%STARTPOS(NTYPE),PARAMP%STARTINTIP(NTYPE))
    ALLOCATE(PARAMP%KHALT(NTYPE,NSTATE),PARAMP%KWALK(NTYPE,NSTATE))
    ALLOCATE(PARAMP%PRETRO(NTYPE,NSTATE),PARAMP%KEXIT(NTYPE))
    ALLOCATE(PARAMP%STARTDIR(NSTATE), PARAMP%PRAD(NTYPE), PARAMP%NPROT(NTYPE), &
         & PARAMP%STARTSTATE(NTYPE))
    ALLOCATE(PARAMP%KMAX(MAXNPROT,MAXNPROT,NTYPE),PARAMP%KM(NTYPE))
    ALLOCATE(PARAMP%PFUSE(NTYPE,NTYPE),&
              PARAMP%MAXNFUSE(NTYPE),PARAMP%KFUSE(NTYPE,NTYPE))
    ALLOCATE(PARAMP%STARTPROT(NTYPE,MAXNPROT))
    ALLOCATE(PARAMP%STARTACID(NTYPE),PARAMP%ACIDIFICATIONRATE(NTYPE))
    ALLOCATE(PARAMP%SWITCHONFUSE(NTYPE))
    ALLOCATE(PARAMP%SPAWNLENGTH(NTYPE))

    PARAMP%DELT = DELT
    PARAMP%CURTIME = 0D0

    PARAMP%MAXNPART = MAXNPART
    PARAMP%VEL = VEL(1:NTYPE,1:2)
    PARAMP%D = DIFFCONST(1:NTYPE)
    PARAMP%KPROD = KPROD(1:NTYPE)

    PARAMP%KSWITCH = KSWITCH(1:NSTATE,1:NSTATE,1:NTYPE)
    PARAMP%KHALT = KHALT(1:NTYPE,1:NSTATE)
    PARAMP%KWALK = KWALK(1:NTYPE,1:NSTATE)
    PARAMP%PRETRO = PRETRO(1:NTYPE,1:NSTATE)
    PARAMP%KEXIT = KEXIT(1:NTYPE)

    PARAMP%KMAT = KMAT(1:NTYPE)
    PARAMP%STARTPOS = STARTPOS(1:NTYPE)
    PARAMP%STARTINTIP = STARTINTIP(1:NTYPE)
    PARAMP%SPAWNLENGTH = SPAWNLENGTH(1:NTYPE)
    PARAMP%STARTDIR = STARTDIR(1:NSTATE)
    PARAMP%STARTSTATE = STARTSTATE(1:NTYPE)
    PARAMP%PRAD = PRAD(1:NTYPE)
    PARAMP%NPROT = NPROT(1:NTYPE)

    PARAMP%KMAX = KMAX(:,:,1:NTYPE)
    PARAMP%KM = KMVALS(1:NTYPE)
    PARAMP%STARTACID = STARTACID(1:NTYPE)
    PARAMP%ACIDIFICATIONRATE = ACIDIFICATIONRATE(1:NTYPE)
    
    PARAMP%PFUSE = PFUSE(1:NTYPE,1:NTYPE)
    PARAMP%KFUSE = KFUSE(1:NTYPE,1:NTYPE)
    PARAMP%MAXNFUSE = MAXNFUSE(1:NTYPE)
    PARAMP%USEKFUSE = USEKFUSE
    PARAMP%USEINTERACTIONS = USEINTERACTIONS

    PARAMP%STARTPROT = STARTPROT(1:NTYPE,:)

    PARAMP%SWITCHONFUSE = SWITCHONFUSE(1:NTYPE)

    PARAMP%USEDILUTION = USEDILUTION

    PARAMP%DOMLEN = DOMLEN
    PARAMP%REFBOUND = REFBOUND

    PARAMP%ARRAYSET = .TRUE.    
  END SUBROUTINE SETUPPARTICLEPARAMS
  
  SUBROUTINE MAKEPARTICLE(PARTP,PARTID,PTYPE,PARAMP,CANFUSE)
  ! make a new particle and assign parameters
    USE MT19937, ONLY: GRND
    USE KEYS,ONLY: MAXNPART,NSTATE
    IMPLICIT NONE

    TYPE(PARTICLE), POINTER :: PARTP
    INTEGER, INTENT(IN) :: PTYPE, PARTID
    TYPE(PARTICLEPARAMS), POINTER :: PARAMP
    LOGICAL, INTENT(IN), OPTIONAL :: CANFUSE
    DOUBLE PRECISION :: ELL,EXPRND
    
    INTEGER :: NTYPE, NPROT
    
    IF(PARTP%ARRAYSET) THEN
      PRINT*, 'ERROR IN MAKEPARTICLE: trying to overwrite particle'
      PRINT*,'Existing PID, new PID:', PARTP%ID,PARTID
      STOP 1
    END IF
    
    PARTP%ID = PARTID
    PARTP%PTYPE = PTYPE
    PARTP%PRAD = PARAMP%PRAD(PTYPE)
    PARTP%TMAKE = PARAMP%CURTIME
    
    IF(PARAMP%SPAWNLENGTH(PTYPE).EQ.0D0) THEN
    ! if spawnlength is zero, particle spawns at designated location
      PARTP%POS = PARAMP%STARTPOS(PTYPE)
    ELSEIF(PARAMP%SPAWNLENGTH(PTYPE).GT.0D0) THEN
    ! if spawnlength is positive, particle spawns with exponentially decaying position
      ELL = PARAMP%SPAWNLENGTH(PTYPE)
      EXPRND = -ELL*LOG(1-(1-EXP(-PARAMP%DOMLEN/ELL))*GRND())
      IF(PARAMP%STARTPOS(PTYPE).EQ.PARAMP%DOMLEN) THEN
        PARTP%POS = MAX(MIN(PARAMP%DOMLEN-EXPRND,PARAMP%DOMLEN),0D0)
      ELSE
        PARTP%POS = MAX(MIN(EXPRND,PARAMP%DOMLEN),0D0)
      END IF
    ELSE
    ! if spawnlength is negative, particle spawns uniformly over spawnlength
      IF(PARAMP%STARTPOS(PTYPE).EQ.PARAMP%DOMLEN) THEN
        PARTP%POS = PARAMP%DOMLEN-GRND()*ABS(PARAMP%SPAWNLENGTH(PTYPE))
      ELSE
        PARTP%POS = GRND()*ABS(PARAMP%SPAWNLENGTH(PTYPE))
      END IF
    END IF

    PARTP%ISINTIP = PARAMP%STARTINTIP(PTYPE)

    PARTP%PREVPOS = PARTP%POS

    PARTP%CURSTATE = PARAMP%STARTSTATE(PTYPE)
    PARTP%CURDIR = PARAMP%STARTDIR(PARTP%CURSTATE)

    IF(PARTP%CURDIR.EQ.0) THEN
      PARTP%VEL = 0
    ELSE
      PARTP%VEL = PARAMP%VEL(PTYPE,(PARTP%CURDIR+3)/2)
    END IF
    PARTP%D = PARAMP%D(PTYPE)

    ALLOCATE(PARTP%KSWITCH(NSTATE,NSTATE))    
    PARTP%KSWITCH = PARAMP%KSWITCH(:,:,PTYPE)
    ALLOCATE(PARTP%KHALT(NSTATE))
    PARTP%KHALT = PARAMP%KHALT(PTYPE,:)
    ALLOCATE(PARTP%KWALK(NSTATE))
    PARTP%KWALK = PARAMP%KWALK(PTYPE,:)
    PARTP%KEXIT = PARAMP%KEXIT(PTYPE)
    ALLOCATE(PARTP%PRETRO(NSTATE))
    PARTP%PRETRO = PARAMP%PRETRO(PTYPE,:)

    PARTP%KMAT = PARAMP%KMAT(PTYPE)

    IF(PRESENT(CANFUSE)) THEN
      PARTP%CANFUSE = CANFUSE
    ELSE
      PARTP%CANFUSE = .TRUE.
    END IF
    PARTP%NFUSE = 0
    PARTP%NTOUCH = 0
    ALLOCATE(PARTP%INTERACTING(MAXNPART))
    PARTP%INTERACTING = 0
    PARTP%NINT = 0

    PARTP%TDEG = 0D0
    PARTP%XFUSE = -1D0
    PARTP%TFUSE = -1D0
    PARTP%TSWITCH = -1D0

    NTYPE = PARAMP%NTYPE
    PARTP%NFUSETYPE = NTYPE

    ! allocate either rates or probs of fusion
    ALLOCATE(PARTP%PFUSE(NTYPE))
    PARTP%PFUSE = PARAMP%PFUSE(:,PTYPE)
    PARTP%MAXNFUSE = PARAMP%MAXNFUSE(PTYPE)
    ALLOCATE(PARTP%KFUSE(NTYPE))
    PARTP%KFUSE = PARAMP%KFUSE(:,PTYPE)

    NPROT = PARAMP%NPROT(PTYPE)
    PARTP%NPROT = NPROT
    IF(NPROT.GT.0) THEN
      ALLOCATE(PARTP%PROT(NPROT))
      PARTP%PROT = PARAMP%STARTPROT(PTYPE,1:PARTP%NPROT)
      ALLOCATE(PARTP%PROTCONV(NPROT,NPROT))
      IF(PARAMP%KM(PARTP%PTYPE).GT.0D0) THEN
        PARTP%PROTCONV = 0D0
      ELSE
        PARTP%PROTCONV = PARAMP%KMAX(:,:,PTYPE)
      END IF
    END IF

    ! allocate initial acidification
    PARTP%ACIDITY = PARAMP%STARTACID(PTYPE)

    ! allocate acidification rate
    PARTP%ACIDIFICATIONRATE = 0D0
    ! PARTP%ACIDIFICATIONRATE = PARAMP%ACIDIFICATIONRATE(PTYPE)

    PARTP%ARRAYSET= .TRUE.

  END SUBROUTINE MAKEPARTICLE

  SUBROUTINE REMOVEPARTICLE(PGROUP,IND,DELETEPART)
  ! remove particle from the group
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    INTEGER, INTENT(IN) :: IND
    TYPE(PARTICLE), POINTER :: PARTP
    TYPE(PARTICLE),TARGET :: PTMP
    INTEGER :: NPARTOLD,PTYPE,NPTYPE,PID
    LOGICAL, INTENT(IN), OPTIONAL :: DELETEPART
    LOGICAL :: DODELETE = .TRUE.

    IF(PRESENT(DELETEPART)) DODELETE = DELETEPART

    PTMP = PGROUP%PLIST(IND)

    NPARTOLD = SUM(PGROUP%NPART)
    PGROUP%PLIST(IND:(NPARTOLD-1)) = PGROUP%PLIST((IND+1):NPARTOLD)
    PGROUP%PLIST(NPARTOLD) = PTMP
    
    PARTP=>PGROUP%PLIST(NPARTOLD)
    PTYPE = PARTP%PTYPE
    NPTYPE = PGROUP%NPART(PTYPE)
    PID = PARTP%ID

    IF(DODELETE) CALL DELETEPARTICLE(PARTP)

    PARTP%ARRAYSET = .FALSE.
 
    PGROUP%NPART(PTYPE) = NPTYPE-1    
   
  END SUBROUTINE REMOVEPARTICLE
  
  SUBROUTINE DELETEPARTICLE(PARTP)
  ! delete particle and deallocate arrays
    IMPLICIT NONE
    TYPE(PARTICLE), POINTER :: PARTP
    INTEGER :: PID

    PID = PARTP%ID

    PARTP%ID = 0
    PARTP%TMAKE = -1D0
    PARTP%PTYPE = 0
    PARTP%PRAD = 0
    PARTP%POS = 0D0
    PARTP%PREVPOS = 0D0
    PARTP%CURDIR = 0
    PARTP%CURSTATE = 0
    PARTP%VEL = 0D0
    PARTP%D = 0D0
    PARTP%KSWITCH = 0D0
    PARTP%KHALT = 0D0
    PARTP%KWALK = 0D0
    PARTP%KEXIT = 0D0
    PARTP%PRETRO = 0D0
    PARTP%KMAT = 0D0

    PARTP%CANFUSE = .FALSE.
    DEALLOCATE(PARTP%INTERACTING)
    PARTP%NINT = 0
    PARTP%NFUSE = 0
    PARTP%NTOUCH = 0
    PARTP%NFUSETYPE = 0
    IF(SIZE(PARTP%PFUSE).GT.0) THEN
      DEALLOCATE(PARTP%PFUSE)
    END IF
    IF(SIZE(PARTP%KFUSE).GT.0) THEN
      DEALLOCATE(PARTP%KFUSE)
    END IF
    
    IF(PARTP%NPROT.GT.0) THEN
      PARTP%NPROT = 0
      PARTP%PROTCONV = 0D0
      IF(SIZE(PARTP%PROTCONV).GT.1) DEALLOCATE(PARTP%PROTCONV)
      IF(SIZE(PARTP%PROT).GT.1) DEALLOCATE(PARTP%PROT)
    END IF
    PARTP%ACIDITY = 0D0
    PARTP%ACIDIFICATIONRATE = 0D0

    PARTP%ARRAYSET = .FALSE.
  END SUBROUTINE DELETEPARTICLE

  SUBROUTINE WRITENTOUCH(PARTP,NTOUCHFILE,APPEND,NTRIALS)
  ! output number of touches
    IMPLICIT NONE

    TYPE(PARTICLE), POINTER :: PARTP
    CHARACTER(LEN=*), INTENT(IN) :: NTOUCHFILE
    LOGICAL, INTENT(IN) :: APPEND
    INTEGER, INTENT(IN), OPTIONAL :: NTRIALS

    IF (APPEND) THEN
      OPEN(UNIT=22,FILE=NTOUCHFILE,STATUS='UNKNOWN',POSITION='APPEND')
      WRITE(22,'(I8)') PARTP%NTOUCH
    ELSE
      OPEN(UNIT=22,FILE=NTOUCHFILE,STATUS='UNKNOWN')
      IF(PRESENT(NTRIALS)) WRITE(22,'(I8)') NTRIALS
    END IF

    CLOSE(22)

  END SUBROUTINE WRITENTOUCH
  
  SUBROUTINE OUTPUTSNAPSHOT(GROUPLIST,FILENAME,INFO,NTRIALS,APPEND)
  ! Output a snapshot of current configuration

    IMPLICIT NONE

    TYPE(PARTICLEGROUP), TARGET :: GROUPLIST(NTRIALS)
    TYPE(PARTICLEGROUP), POINTER :: PGROUP
    TYPE(PARTICLE), POINTER :: PLIST(:)
    INTEGER, INTENT(IN) :: NTRIALS
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    DOUBLE PRECISION, INTENT(IN) :: INFO(:)
    LOGICAL, INTENT(IN) :: APPEND

    INTEGER :: TRC,PC,NINFO,INTIP
    CHARACTER(LEN=2) :: X1
    CHARACTER(LEN=16) :: FMT1,FMT2

    IF (APPEND) THEN
      OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN',POSITION='APPEND')
    ELSE
      OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN')
    END IF

    ! Write general info for this snapshot
    NINFO = SIZE(INFO)
    WRITE (X1,'(I2)') NINFO
    FMT1 = '('//TRIM(X1)//'ES15.6E2)'
    WRITE(99,FMT1) INFO

    DO TRC = 1,NTRIALS
      PGROUP=>GROUPLIST(TRC)
      PLIST=>PGROUP%PLIST
      WRITE(99,'(2I8)') TRC,SUM(PGROUP%NPART)
      DO PC = 1,SUM(PGROUP%NPART)
        INTIP = 0
        IF(PGROUP%PLIST(PC)%ISINTIP) INTIP = 1
        WRITE(99,'(7I12)') PC,PLIST(PC)%ID,PLIST(PC)%PTYPE,PLIST(PC)%NPROT,&
                          PLIST(PC)%NFUSE,PLIST(PC)%CURSTATE,INTIP
        WRITE (X1,'(I2)') PLIST(PC)%NPROT+6        
        FMT2 = '('//TRIM(X1)//'ES15.6E2)'
        ! PRINT*, 'TESTX1:PC,TYPE,NPROT,FMT', PC, PLIST(PC)%PTYPE, PLIST(PC)%NPROT, fmt2
        IF(PLIST(PC)%NPROT.GT.0) THEN
          ! print*, 'TESTX2:', PLIST(PC)%POS,PLIST(PC)%PROT,PLIST(PC)%ACIDITY
          WRITE(99,FMT2) PLIST(PC)%POS,PLIST(PC)%PROT,PLIST(PC)%ACIDITY,PLIST(PC)%TDEG,&
                &PLIST(PC)%XFUSE,PLIST(PC)%TFUSE,PLIST(PC)%TSWITCH
        ELSE
          WRITE(99,FMT2) PLIST(PC)%POS,PLIST(PC)%TDEG,PLIST(PC)%XFUSE&
                &,PLIST(PC)%TFUSE,PLIST(PC)%TSWITCH
        END IF
      END DO
    END DO
    WRITE(99,*)
    CLOSE(99)

  END SUBROUTINE OUTPUTSNAPSHOT
  
  SUBROUTINE OUTPUTPROT(CURTIME,CELLPROT,AXPROT,FILENAME,NPROT,NTRIALS,APPEND)
  ! Output protein content in the cell

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: CURTIME
    DOUBLE PRECISION, INTENT(IN) :: CELLPROT(NPROT),AXPROT(NPROT)
    INTEGER, INTENT(IN) :: NPROT,NTRIALS

    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    LOGICAL, INTENT(IN) :: APPEND

    CHARACTER(LEN=2) :: X1

    IF (APPEND) THEN
      OPEN(UNIT=44,FILE=FILENAME,STATUS='UNKNOWN',POSITION='APPEND')
    ELSE
      OPEN(UNIT=44,FILE=FILENAME,STATUS='UNKNOWN')
    END IF

    IF(CURTIME.EQ.0D0) THEN
      WRITE(44,'(2I8)') NPROT, NTRIALS
    ELSE
      WRITE(X1,'(I2)') 2*NPROT+1
      WRITE(44,'('//TRIM(X1)//'ES15.6E2)') CURTIME,CELLPROT,AXPROT
    END IF
    CLOSE(44)
  
  END SUBROUTINE OUTPUTPROT
  
  SUBROUTINE CLEANUPPARTICLEGROUP(PGROUP)
  ! clean up particlegroup object and deallocate arrays
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PGROUP

    PGROUP%NPART = 0
    DEALLOCATE(PGROUP%PLIST)
    PGROUP%LASTID = 0
    PGROUP%ARRAYSET = .FALSE.
  END SUBROUTINE CLEANUPPARTICLEGROUP

  SUBROUTINE CLEANUPPARTICLEPARAMS(PARAMP)
  ! deallocate object arrays
    IMPLICIT NONE
    TYPE(PARTICLEPARAMS), POINTER :: PARAMP

    DEALLOCATE(PARAMP%VEL,PARAMP%D,PARAMP%KPROD,PARAMP%KSWITCH,PARAMP%KMAT,&
              & PARAMP%STARTPOS,PARAMP%STARTINTIP,PARAMP%KHALT,PARAMP%KWALK,&
              & PARAMP%PRETRO,PARAMP%KEXIT,PARAMP%STARTDIR,PARAMP%PRAD,&
              & PARAMP%NPROT,PARAMP%STARTSTATE,PARAMP%KMAX,&
              & PARAMP%KM,PARAMP%PFUSE,PARAMP%MAXNFUSE,PARAMP%KFUSE,PARAMP%STARTPROT,&
              & PARAMP%STARTACID,PARAMP%ACIDIFICATIONRATE,PARAMP%SWITCHONFUSE,&
              & PARAMP%SPAWNLENGTH)
         
    PARAMP%ARRAYSET = .FALSE.
    
  END SUBROUTINE CLEANUPPARTICLEPARAMS
  
END MODULE PARTICLEUTILS