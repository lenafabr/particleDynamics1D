SUBROUTINE READKEY
! this subroutine reads in keywords from a parameter file
! it sets the various global variables defined in KEYS module
! name of the parameter file is param.* where * is a keyword argument
! if no keyword argument is supplied, the default is just a file called param
! The EXTRAPARAMFILES keyword will allow extra parameter files to be 
! read in as well

  USE KEYS
  USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
  USE GENUTIL

  IMPLICIT NONE

  ! ---- stuff for inputing the parameter file in free format --------
  CHARACTER*100 :: ARG ! command line argument
  INTEGER :: NUMARG ! number of command line arguments
  INTEGER :: NITEMS ! number of items on the line in the parameter file
  INTEGER :: PF ! input file unit
  LOGICAL :: FILEEND=.FALSE. ! done reading file?
  CHARACTER*100 :: WORD ! keyword
  ! -------------- for reading multiple parameter files --------  
  INTEGER, PARAMETER :: MAXNFILES = 10
  CHARACTER*100 :: PARAMFILES(MAXNFILES)
  INTEGER :: NPARAMFILES, NPARAMREAD
  ! ------ for initializing random number generator
  INTEGER :: TIMEVAL(8), SEED
  ! ---------------- temporary variables ---------------
  INTEGER :: I,J,K
  CHARACTER*100 :: DUMSTR
  LOGICAL :: LDUM

  ! ------------------------
  ! set variable defaults
  ! ------------------------

  ! general program control
  ACTION = 'NONE'
  RNGSEED = 0
  VERBOSE = .FALSE.

  ! simulation info
  NTRIALS = 1 ! number of trials

  ! particle properties
  NPART = 1 ! number of particles
  NTYPE = 1 ! number of particle types
  VEL = 1D0 ! velocity of particles
  PRAD = 1D-2 ! particle radius (different radii can be assigned to different particles)
  
  KSWITCH = 0D0 ! rate of changing the state

  KHALT = 0D0 ! rate of stopping
  KWALK = 0D0 ! rate of starting a walk
  PRETRO = 0D0 ! probability of walking in the retrograde direction

  KMAT = 0D0 ! rate of maturation
  
  SWITCHONFUSE = .FALSE. ! switch states after fusion

  DIFFCONST = 0D0 ! diffusion coefficient

  ! protein dynamics
  NPROT = 1
  KMAX = 0D0
  KMVALS = 0D0
  KMAX_A = 0D0
  KMVALS_A = 0D0
  STARTACID = 0D0
  ACIDIFICATIONRATE = 0D0
  USEDILUTION = .TRUE.

  ! fusion details
  PFUSE = 0D0 ! matrix of fusion probabilities
  KFUSE = 0D0 ! matrix of fusion rates
  MAXNFUSE = HUGE(1) ! number of fusion events needed for saturation
  DOFUSION = .FALSE.
  USEKFUSE = .FALSE.
  USEINTERACTIONS = .FALSE.
  
  ! particle production
  KPROD = 0D0
  STARTPOS = -1D0 ! starting position
  STARTINTIP = .FALSE. ! start particles in tip?
  STARTPROT = 0D0 ! amount of protein to start with
  STARTSTATE = 1 ! initial state of particle
  STARTDIR = 0 ! intiial direction of particle
  SPAWNLENGTH = 0D0 ! distance over which particles spawn

  ! domain info
  DOMLEN = 1D0 ! domain length
  REFBOUND = (/.FALSE.,.FALSE./)

  ! input/output  
  SNAPFILE = '*.snap.out' ! snapshot file
  OUTFILE = '*.out' ! output file
  NFUSEFILE = '*.fuse.out' ! file to write spatial distribution of fusion events
  NTOUCHFILE = '*.touch.out' ! file to write number of touches
  DUMPSNAPSHOTS = .FALSE. ! periodically dump snapshots
  APPENDSNAPSHOTS = .FALSE. ! append snapshots to file rather than replacing
  CALCNTOUCH = .FALSE. ! calculate number of touches
  SNAPEVERY = 0 ! how often to dump snapshots
  OUTEVERY = 0 ! how often to write to output file
  FUSEDISTEVERY = 0 !how often to write distribution of fusion events
  
  ! dynamics
  DELT = 1D-4 ! time step
  NSTEPS = INT(1D6) ! number of brownian steps to run
  PRINTEVERY = 0 ! how often to print output
    
  ! -------------------------
  ! Read in all parameter files, starting with the ones specified on command line
  ! --------------------------

  PF = 55 ! i/o unit number to be used for parameter files

  ! get input parameter files from command line
  NPARAMFILES = 0
  NUMARG = COMMAND_ARGUMENT_COUNT()  
  IF (NUMARG==0) THEN
    NPARAMFILES = 1
    PARAMFILES(1) = 'param'
    ARG = ''
  ELSE
    DO I = 1,NUMARG
      CALL GETARG(I, ARG)
      NPARAMFILES = NPARAMFILES + 1
      WRITE(DUMSTR,'(A)') 'param.' //TRIM(ADJUSTL(ARG))
      PARAMFILES(NPARAMFILES) = DUMSTR
    END DO
    ! reset arg to its original value
    IF (NUMARG.GT.1) CALL GETARG(1,ARG)
  ENDIF

  NPARAMREAD = 0 ! keep track of how many files have been read
  DO WHILE (NPARAMREAD.LT.NPARAMFILES)
    NPARAMREAD = NPARAMREAD + 1

    PRINT*, 'Reading parameter file: ', PARAMFILES(NPARAMREAD)
    INQUIRE(FILE=PARAMFILES(NPARAMREAD),EXIST=LDUM)
    IF (.NOT.LDUM) THEN
      PRINT*, 'ERROR in READKEY: Parameter file ', TRIM(ADJUSTL(PARAMFILES(NPARAMREAD))), ' does not exist.'
      STOP 1
    ENDIF
    OPEN(UNIT=PF, FILE=PARAMFILES(NPARAMREAD), STATUS='OLD')

    ! read in the keywords one line at a time
    DO 
      CALL READLINE(PF,FILEEND,NITEMS)
      IF (FILEEND.and.nitems.eq.0) EXIT

      ! skip empty lines
      IF (NITEMS.EQ.0) CYCLE

      ! Read in the keyword for this line
      CALL READA(WORD,CASESET=1)

      ! Skip any empty lines or any comment lines
      IF (WORD(1:1).EQ.'#') CYCLE

      SELECT CASE(WORD) ! pick which keyword
        CASE('ACTION')
          CALL READA(ACTION, CASESET=1)
        CASE('ACIDIFICATIONRATE')
          CALL READF(ACIDIFICATIONRATE(1))
          ACIDIFICATIONRATE = ACIDIFICATIONRATE(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(ACIDIFICATIONRATE(I))
          END DO
        CASE('DIFFCONST')
          CALL READF(DIFFCONST(1))
          DIFFCONST = DIFFCONST(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(DIFFCONST(I))
          END DO
        CASE('DOMLEN')
          CALL READF(DOMLEN)
        CASE('DELT')
          CALL READF(DELT)
        CASE('DOFUSION')
          CALL READO(DOFUSION)
        CASE('FUSEDIST')
          IF (NITEMS.GT.1) CALL READI(FUSEDISTEVERY)
          IF (NITEMS.GT.2) CALL READA(NFUSEFILE)
          ! IF (NITEMS.GT.3) CALL READO(APPENDSNAPSHOTS)
        CASE('KEXIT')   
          CALL READF(KEXIT(1))
          KEXIT = KEXIT(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(KEXIT(I))
          END DO          
        CASE('KFUSE')
          CALL READI(I)
          DO J = 1,MIN(MAXNTYPE,NITEMS-2)
            CALL READF(KFUSE(I,J))
            KFUSE(J,I) = KFUSE(I,J)
          END DO
        CASE('KHALT')
          CALL READI(I)
          DO J = 1,MIN(NSTATE,NITEMS-2)
            CALL READF(KHALT(I,J))
          END DO
        CASE('KMAT')
          CALL READF(KMAT(1))
          KMAT = KMAT(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(KMAT(I))
          END DO           
        CASE('KMAX')
          CALL READI(I)
          CALL READI(J)
          DO K = 1,MIN(MAXNTYPE,NITEMS-3)
            CALL READF(KMAX(I,J,K))
          END DO
        CASE('KMAX_A')   
          CALL READF(KMAX_A(1))
          KMAX_A = KMAX_A(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(KMAX_A(I))
          END DO
        CASE('KMVALS')
          CALL READF(KMVALS(1))
          KMVALS = KMVALS(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(KMVALS(I))
          END DO
        CASE('KMVALS_A')
          KMVALS_A = KMVALS_A(1)
          DO I = 1,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(KMVALS_A(I))
          END DO
        CASE('KPROD')
          CALL READF(KPROD(1))
          KPROD = KPROD(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(KPROD(I))
          END DO
        CASE('KSWITCH')
          CALL READI(I)
          CALL READI(J)
          DO K = 1,MIN(MAXNTYPE,NITEMS-3)
            CALL READF(KSWITCH(I,J,K))
          END DO
        CASE('KWALK')
          CALL READI(I)
          DO J = 1,MIN(NSTATE,NITEMS-2)
            CALL READF(KWALK(I,J))
          END DO  
        CASE('MAXNFUSE')
          CALL READI(MAXNFUSE(1))
          MAXNFUSE = MAXNFUSE(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READI(MAXNFUSE(I))
          END DO
        CASE('NPART')
          CALL READI(NPART(1))
          NPART = NPART(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READI(NPART(I))
          END DO
        CASE('NPROT')
          CALL READI(NPROT(1))
          NPROT = NPROT(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READI(NPROT(I))
          END DO
        CASE('NSTEPS')
          CALL READI(NSTEPS)
        CASE('NTOUCHFILE')
          CALCNTOUCH = .TRUE.
          CALL READA(NTOUCHFILE)
        CASE('NTRIALS')
          CALL READI(NTRIALS)
        CASE('NTYPE')
          CALL READI(NTYPE)
          IF(NTYPE.GT.MAXNTYPE) THEN
            PRINT*, 'Maximum number of particle types allowed is ',MAXNTYPE
            NTYPE = MAXNTYPE
          END IF
        CASE('OUTPUT')
          IF (NITEMS.GT.1) CALL READI(OUTEVERY)
          IF (NITEMS.GT.2) CALL READA(OUTFILE)
        CASE('OUTFILE')
          CALL READA(OUTFILE)
        CASE('PFUSE')
          CALL READI(I)
          DO J = 1,MIN(MAXNTYPE,NITEMS-2)
            CALL READF(PFUSE(I,J))
            PFUSE(J,I) = PFUSE(I,J)
          END DO
        CASE('PRAD')   
          CALL READF(PRAD(1))
          PRAD = PRAD(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(PRAD(I))
          END DO
        CASE('PRETRO')
          CALL READI(I)
          DO J = 1,MIN(NSTATE,NITEMS-2)
            CALL READF(PRETRO(I,J))
          END DO
        CASE('PRINTEVERY')
          CALL READI(PRINTEVERY)
          IF(PRINTEVERY.LT.0) PRINTEVERY = NSTEPS/10
        CASE('REFBOUND')
          DO I = 1,NITEMS-1
            CALL READO(REFBOUND(I))
          END DO
        CASE('RNGSEED')
          CALL READI(RNGSEED)
        CASE('SNAPFILE')
          CALL READA(SNAPFILE)  
        CASE('SNAPSHOTS')
          DUMPSNAPSHOTS = .TRUE.
          IF (NITEMS.GT.1) THEN
            CALL READI(SNAPEVERY)
            IF(SNAPEVERY.LT.0) SNAPEVERY = NSTEPS 
          END IF
          IF (NITEMS.GT.2) CALL READA(SNAPFILE)
          IF (NITEMS.GT.3) CALL READO(APPENDSNAPSHOTS)
        CASE('SPAWNLENGTH')
          CALL READF(SPAWNLENGTH(1))
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(SPAWNLENGTH(I))
          END DO
        CASE('STARTACID')
          CALL READF(STARTACID(1))
          STARTACID = STARTACID(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(STARTACID(I))
          END DO
        CASE('STARTDIR')
          DO I = 1,MIN(NSTATE,NITEMS-1)
            CALL READI(STARTDIR(I))
          END DO
        CASE('STARTINTIP')
          DO I = 1,MIN(MAXNTYPE,NITEMS-1)
            CALL READO(STARTINTIP(I))
          END DO          
        CASE('STARTPOS')
          CALL READF(STARTPOS(1))
          STARTPOS = STARTPOS(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READF(STARTPOS(I))
          END DO
        CASE('STARTPROT')
          CALL READI(I)
          DO J = 1,MIN(MAXNPROT,NITEMS-2)
            CALL READF(STARTPROT(I,J))
          END DO
        CASE('STARTSTATE')
          CALL READI(STARTSTATE(1))
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READI(STARTSTATE(I))
          END DO
        CASE('SWITCHONFUSE')
          CALL READO(SWITCHONFUSE(1))
          SWITCHONFUSE = SWITCHONFUSE(1)
          DO I = 2,MIN(MAXNTYPE,NITEMS-1)
            CALL READO(SWITCHONFUSE(I))
          END DO
        CASE('USEINTERACTIONS')
          CALL READO(USEINTERACTIONS)
        CASE('USEKFUSE')
          CALL READO(USEKFUSE)
        CASE('USEDILUTION')
          CALL READO(USEDILUTION)
        CASE('VEL')
          CALL READI(I)
          DO J = 1,MIN(2,NITEMS-2)
            CALL READF(VEL(I,J))
          END DO
        CASE('VERBOSE')
          CALL READO(VERBOSE)
        CASE DEFAULT
          print*, 'ERROR: unidentified keyword ', TRIM(WORD), " Will ignore."
      END SELECT
    ENDDO
    CLOSE(PF)
  ENDDO

  ! -----------------
  ! check validity of some values, raise errors or adjust as necessary
  ! -----------------  

  ! ----------- fix file names -----------
  CALL REPLACESUBSTR(SNAPFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(OUTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(NFUSEFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(NTOUCHFILE,'*',TRIM(ADJUSTL(ARG)))
  ! ---------------------------

  ! Initiate random number generator 
  IF (RNGSEED.EQ.0) THEN
    ! use the current time of day in milliseconds
    CALL DATE_AND_TIME(VALUES=TIMEVAL)
    SEED = TIMEVAL(5)*3600*1000 + TIMEVAL(6)*60*1000 + TIMEVAL(7)*1000 + TIMEVAL(8)
  ELSEIF (RNGSEED.EQ.-1) THEN
    ! use the last 5 characters in the command-line argument
    SEED = STRING2NUM(TRIM(ADJUSTL(ARG)))    
  ELSEIF (RNGSEED.EQ.-2) THEN
    ! use the last 4 characters in the command-line argument 
    ! and additionally the millisecond time 
    CALL DATE_AND_TIME(VALUES=TIMEVAL)
    SEED = STRING2NUM(TRIM(ADJUSTL(ARG)),TIMEVAL(8))
  ELSE
    ! use this seed directly
    SEED = RNGSEED
  ENDIF

  print*, 'Initiating Mersenne twister random number generator with seed:', SEED
  CALL SGRND(SEED)

  print*, '------------Parameter values : -------------------'
  print*, 'ACTION: ', TRIM(ADJUSTL(ACTION))
  IF (DUMPSNAPSHOTS) THEN
    PRINT*, 'Dumping snapshot every', SNAPEVERY,'steps. In file:', TRIM(ADJUSTL(SNAPFILE))
  ENDIF
  
  print*, 'number of particle types:', NTYPE
  print*, 'velocities:', VEL

  IF(DOMLEN.GT.0D0) PRINT*, 'domain length: ',DOMLEN
  
  PRINT*, '----------------------------------------------------'
	
END SUBROUTINE READKEY
