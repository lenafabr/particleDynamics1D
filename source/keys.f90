MODULE KEYS
  ! keyword parameters that are globally used in many different places in the code
  IMPLICIT NONE

  ! ----------------------------------------------------------------------------
  ! general program control
  ! ----------------------------------------------------------------------------
  CHARACTER*100 :: ACTION
  INTEGER :: RNGSEED
  LOGICAL :: VERBOSE
  
  ! ----------------------------------------------------------------------------
  ! simulation info
  ! ----------------------------------------------------------------------------
  INTEGER:: NTRIALS

  ! ----------------------------------------------------------------------------
  ! particle properties
  ! ----------------------------------------------------------------------------
  INTEGER, PARAMETER :: MAXNPART = 5000 ! maximum number of particles allowed
  INTEGER, PARAMETER :: MAXNTYPE = 2 ! maximum number of particle types
  INTEGER, PARAMETER :: NSTATE = 4 ! number of particle states
                                   ! 1 = bidirectional with stops
                                   ! 2 = fully retrograde with stops
                                   ! 3 = fully anterograde with stops
                                   ! 4 = fully stationary
  INTEGER :: NTYPE ! number of particle types
  INTEGER :: NPART(MAXNTYPE) ! number of active particles in a group
  DOUBLE PRECISION :: VEL(MAXNTYPE,2) ! velocity of each particle type in each direction
                                      ! VEL(1) in the retrograde direction
                                      ! VEL(2) in the anterograde direction
  DOUBLE PRECISION :: PRAD(MAXNTYPE) ! radius of each particle type
  
  DOUBLE PRECISION :: KMAT(MAXNTYPE) ! maturation rates

  DOUBLE PRECISION :: KSWITCH(NSTATE,NSTATE,MAXNTYPE) ! rate to switch between states
  LOGICAL :: SWITCHONFUSE(MAXNTYPE) ! switch state on fusing?

  DOUBLE PRECISION :: DIFFCONST(MAXNTYPE) ! diffusion coefficient for diffusive particles
  INTEGER :: STARTDIR(NSTATE) ! initial particle direction in each state

  INTEGER :: STARTSTATE(MAXNTYPE) ! initial state of particle

  DOUBLE PRECISION :: KHALT(MAXNTYPE,NSTATE) ! stopping rates for particles in various states

  DOUBLE PRECISION :: KWALK(MAXNTYPE,NSTATE) ! walking rates for particles in various states

  DOUBLE PRECISION :: PRETRO(MAXNTYPE,NSTATE) ! probability of engaging in retrograde motion

  DOUBLE PRECISION :: KEXIT(MAXNTYPE) ! rate of exiting tip region

  ! ----------------------------------------------------------------------------
  ! protein and acidification dynamics
  ! ----------------------------------------------------------------------------
  INTEGER, PARAMETER :: MAXNPROT = 2 ! maximum number of protein types
  INTEGER :: NPROT(MAXNTYPE) ! number of protein types in each particle type
  DOUBLE PRECISION :: KMAX(MAXNPROT,MAXNPROT,MAXNTYPE) ! conversion rate between proteins
  DOUBLE PRECISION :: KMVALS(MAXNTYPE) !MM constant for protein degradation
  DOUBLE PRECISION :: KMAX_A(MAXNTYPE) ! maximum acidification rate
  DOUBLE PRECISION :: KMVALS_A(MAXNTYPE) !MM constant for protein degradation
  DOUBLE PRECISION :: STARTACID(MAXNTYPE) ! initial acidity in each particle type
  DOUBLE PRECISION :: ACIDIFICATIONRATE(MAXNTYPE) ! acidification rate in each particle type
  LOGICAL :: USEDILUTION ! use dilution to increase acidification?

  ! ----------------------------------------------------------------------------
  ! fusion parameters
  ! ----------------------------------------------------------------------------
  DOUBLE PRECISION :: PFUSE(MAXNTYPE,MAXNTYPE) !matrix for fusion probabilities
  DOUBLE PRECISION :: KFUSE(MAXNTYPE,MAXNTYPE) !matrix for fusion rates
  LOGICAL :: DOFUSION ! implement fusion?
  INTEGER :: MAXNFUSE(MAXNTYPE)
  LOGICAL :: USEKFUSE ! use rates instead of fusion probability?
  LOGICAL :: USEINTERACTIONS ! use explicit interactions to determine fusion?

  ! ----------------------------------------------------------------------------
  ! particle production
  ! ----------------------------------------------------------------------------
  DOUBLE PRECISION :: KPROD(MAXNTYPE) ! rate of production
  DOUBLE PRECISION :: STARTPOS(MAXNTYPE) ! starting position
  LOGICAL :: STARTINTIP(MAXNTYPE) ! start particles in tip?
  DOUBLE PRECISION :: SPAWNLENGTH(MAXNTYPE) ! length over which particles spawn
                                            ! positive values: exponential
                                            ! zero (default): start at startpos
                                            ! negative values: uniform over spawnlength
  DOUBLE PRECISION :: STARTPROT(MAXNTYPE,MAXNPROT) ! amount of protein to start with
 
  ! ----------------------------------------------------------------------------
  ! domain info
  ! ----------------------------------------------------------------------------
  DOUBLE PRECISION :: DOMLEN ! length of domain
  LOGICAL :: REFBOUND(2) ! reflecting boundary status
 
  ! ----------------------------------------------------------------------------
  ! input/output
  ! ----------------------------------------------------------------------------
  CHARACTER*200 :: SNAPFILE, OUTFILE, NFUSEFILE,NTOUCHFILE
  LOGICAL :: DUMPSNAPSHOTS, APPENDSNAPSHOTS,CALCNTOUCH
  INTEGER :: SNAPEVERY, OUTEVERY, FUSEDISTEVERY
  
  ! ----------------------------------------------------------------------------
  ! dynamics
  ! ----------------------------------------------------------------------------
  DOUBLE PRECISION :: DELT
  INTEGER :: NSTEPS, PRINTEVERY
	
END MODULE KEYS