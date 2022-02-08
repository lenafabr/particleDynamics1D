PROGRAM MAIN
  ! subroutines for testing parts of the code
  USE KEYS
  USE GENUTIL
  USE PARTICLEUTILS, ONLY: PARTICLEPARAMS,PARTICLE,PARTICLEGROUP,DOMAIN,&
                          &SETUPDOMAIN,SETUPPARTICLEPARAMS,MAKEPARTICLE,&
                          SETUPPARTICLEGROUP,CLEANUPPARTICLEGROUP,&
                          CLEANUPPARTICLEPARAMS,CLEANUPDOMAIN,&
                          OUTPUTSNAPSHOT, REMOVEPARTICLE

  ! CALL TESTPARTICLES
  CALL TESTLOC
CONTAINS
  ! SUBROUTINE TESTPARTICLES
  ! ! test making / deleting particles

  !   IMPLICIT NONE
  !   TYPE(PARTICLEPARAMS), TARGET :: PARAMS
  !   TYPE(PARTICLEPARAMS), POINTER :: PARAMP
  !   TYPE(PARTICLE), TARGET :: PARTLIST(MAXNPART)
  !   TYPE(PARTICLE), POINTER :: PLIST(:),PARTP
  !   TYPE(PARTICLEGROUP), TARGET :: GROUP
  !   TYPE(PARTICLEGROUP), POINTER :: PGROUP
  !   TYPE(DOMAIN), TARGET :: DOM
  !   TYPE(DOMAIN), POINTER :: DOMP
  !   INTEGER :: PTYPE,PID,PC,CC
  !   DOUBLE PRECISION :: INFO(1)
    
    
  !   CALL READKEY
  !   ! initialize particle arrays
  !   !  ALLOCATE(GROUPLIST(NTRIALS))
  !   ! ALLOCATE(PLIST(MAXNPART))

  !   ! set up simulation domain
  !   DOMP=>DOM
  !   CALL SETUPDOMAIN(DOMP,DOMLEN,NREG,REGPOS)
  !   ! set up parameters
  !   PARAMP=>PARAMS
  !   CALL SETUPPARTICLEPARAMS(PARAMP)
  !   ALLOCATE(PLIST(PARAMP%MAXNPART))
  !   ! PLIST=>PARTLIST

  !   PID = 1
  !   DO PC = 1,NTYPE
  !     DO CC = 1,NPART(PC)
  !       PARTP=>PLIST(PID)
  !       CALL MAKEPARTICLE(PARTP,PID,PC,PARAMP)
  !       PID = PID+1
  !     END DO
  !   END DO
  !   PGROUP=>GROUP
  !   ! PRINT*, 'TESTX, NPART: ', NPART
  !   CALL SETUPPARTICLEGROUP(PGROUP,NPART,PLIST)

  !   ! PRINT*, 'TESTX, PGROUP%NPART: ', PGROUP%NPART
    
  !   ! CALL MAKEPARTICLE(PARTP,PARTID,PTYPE,PARAMP,.TRUE.)
  !   ! CALL SETUPPARTICLEGROUP(PGROUP,1,PLIST)

  !   INFO = 1
  !   CALL OUTPUTSNAPSHOT(PGROUP,SNAPFILE,INFO,.FALSE.)

  !   INFO = 2
  !   PARTP=>PGROUP%PLIST(3)
  !   PTYPE = 2
  !   CALL MAKEPARTICLE(PARTP,PID,PTYPE,PARAMP)
  !   PGROUP%NPART(PTYPE) = PGROUP%NPART(PTYPE)+1
  !   PID = PID+1
  !   CALL OUTPUTSNAPSHOT(PGROUP,SNAPFILE,INFO,.TRUE.)
    
  !   INFO = 3
  !   ! PARTP=>PGROUP%PLIST(2)
  !   CALL REMOVEPARTICLE(PGROUP,1)
  !   CALL OUTPUTSNAPSHOT(PGROUP,SNAPFILE,INFO,.TRUE.)

  !   INFO = 4
  !   PARTP=>PGROUP%PLIST(SUM(PGROUP%NPART)+1)
  !   PTYPE = 1
  !   CALL MAKEPARTICLE(PARTP,PID,PTYPE,PARAMP)
  !   PGROUP%NPART(PTYPE) = PGROUP%NPART(PTYPE)+1
  !   PID = PID+1
  !   CALL OUTPUTSNAPSHOT(PGROUP,SNAPFILE,INFO,.TRUE.)
    
  !   CALL CLEANUPPARTICLEGROUP(PGROUP)
  !   CALL CLEANUPPARTICLEPARAMS(PARAMP)
  !   CALL CLEANUPDOMAIN(DOMP)
  
  ! END SUBROUTINE TESTPARTICLES

  SUBROUTINE TESTLOC
  ! subroutine to search for index
    USE GENUTIL,ONLY: LINSPACE
    IMPLICIT NONE
    DOUBLE PRECISION :: BINEDGES(11), POS
    INTEGER :: IND

    BINEDGES = LINSPACE(0D0,1D0,11)
    POS = 0.25
    IND = MAXLOC(BINEDGES,DIM=1,MASK=(POS.GE.BINEDGES))
    PRINT*,BINEDGES
    PRINT*,IND

  END SUBROUTINE TESTLOC
END PROGRAM MAIN
