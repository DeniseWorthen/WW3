!> @file
!> @brief Contains module WMWAVEMD.
!> 
!> @author H. L. Tolman @date 22-Mar-2021

#include "w3macros.h"
!/ ------------------------------------------------------------------- /
!>
!> @brief Running the multi-grid version of WAVEWATCH III up to a
!>  given ending time for each grid.
!> 
!> @author H. L. Tolman @date 22-Mar-2021
!>
MODULE WMWAVEMD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |           H. L. Tolman            |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         22-Mar-2021 |
  !/                  +-----------------------------------+
  !/
  !/    13-Jun-2005 : Origination.                        ( version 3.07 )
  !/    30-Jan-2006 : Add static nesting.                 ( version 3.08 )
  !/    25-May-2006 : Add overlapping grids.              ( version 3.09 )
  !/    09-Aug-2006 : Unified point output added.         ( version 3.10 )
  !/    22-Dec-2006 : Final algorith changes for tests.   ( version 3.10 )
  !/    25-Jan-2007 : Tweaking algorithm.                 ( version 3.10 )
  !/    02-Feb-2007 : Replacing MPI_BCAST with WMBCST.    ( version 3.10 )
  !/    07-Feb-2007 : Reintroduce pre-fetching.           ( version 3.10 )
  !/    10-May-2007 : Removing / streamlining WMBCST.     ( version 3.11 )
  !/    21-Jun-2007 : Dedicated output processes.         ( version 3.11 )
  !/    20-Sep-2007 : Fix reset of GRSTAT in 0.b          ( version 3.13 )
  !/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
  !/    20-Aug-2010 : Fix MAPSTA/MAPST2 bug.            ( version 3.14.6 )
  !/    12-Mar-2012 : Use MPI_COMM_NULL for checks.       ( version 3.14 )
  !/    28-Jan-2014 : Add memory hwm to profiling.        ( version 5.00 )
  !/    22-Mar-2021 : Support for air density input       ( version 7.13 )
  !/
  !/    Copyright 2009-2014 National Weather Service (NWS),
  !/       National Oceanic and Atmospheric Administration.  All rights
  !/       reserved.  WAVEWATCH III is a trademark of the NWS. 
  !/       No unauthorized use without permission.
  !/
  !  1. Purpose :
  !
  !     Running the multi-grid version of WAVEWATCH III up to a given
  !     ending time for each grid.
  !
  !  2. Variables and types :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !     ----------------------------------------------------------------
  !
  !  3. Subroutines and functions :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !      WMWAVE    Subr. Public   Wave model initialization.
  !      WMPRNT    Subr. Public   Print action table to log file.
  !      WMBCST    Subr. Public   Non-blocking MPI broadcast.
  !      WMWOUT    Subr. Public   Non-blocking MPI broadcast.
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines and functions used :
  !
  !     See subroutine documentation.
  !
  !  5. Remarks :
  !
  !  6. Switches :
  !
  !     See subroutine documentation.
  !
  !  7. Source code :
  !
  !/ ------------------------------------------------------------------- /
  use wav_shr_flags
  !
  PUBLIC
  !/
CONTAINS
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Run multi-grid version of WAVEWATCH III.
  !>
  !> @param[in] TEND Ending time for calculations for each grid.
  !> @author H. L. Tolman @date 22-Mar-2021
  !>        

  SUBROUTINE WMWAVE ( TEND )
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           H. L. Tolman            |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         22-Mar-2021 |
    !/                  +-----------------------------------+
    !/
    !/    13-Jun-2005 : Origination.                        ( version 3.07 )
    !/    30-Jan-2006 : Add static nesting.                 ( version 3.08 )
    !/    25-May-2006 : Add overlapping grids.              ( version 3.09 )
    !/    09-Aug-2006 : Unified point output added.         ( version 3.10 )
    !/    22-Dec-2006 : Final algorith changes for tests.   ( version 3.10 )
    !/    25-Jan-2007 : Tweaking algorithm.                 ( version 3.10 )
    !/    02-Feb-2007 : Replacing MPI_BCAST with WMBCST.    ( version 3.10 )
    !/    07-Feb-2007 : Reintroduce pre-fetching.           ( version 3.10 )
    !/    10-May-2007 : Removing / streamlining WMBCST.     ( version 3.11 )
    !/    21-Jun-2007 : Dedicated output processes.         ( version 3.11 )
    !/    20-Sep-2007 : Fix reset of GRSTAT in 0.b          ( version 3.13 )
    !/    20-Aug-2010 : Fix MAPSTA/MAPST2 bug sec. 9.a.   ( version 3.14.6 )
    !/    12-Mar-2012 : Use MPI_COMM_NULL for checks.       ( version 3.14 )
    !/    28-Jan-2014 : Add memory hwm to profiling.        ( version 5.00 )
    !/    22-Mar-2021 : Support for air density input       ( version 7.13 )
    !/
    !  1. Purpose :
    !
    !     Run multi-grid version of WAVEWATCH III.
    !
    !  2. Method :
    !
    !     See manual.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       TEND    I.A.   I   Ending time for calculations for each grid.
    !     ----------------------------------------------------------------
    !
    !     Local variables
    !     ----------------------------------------------------------------
    !       J       Int.   Group counter.
    !       JJ      Int.   Grid in group counter.
    !       I       Int.   Grid counter.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      W3SETG    Subr  W3GDATMD Point to grid/model.
    !      W3SETW    Subr  W3WDATMD Point to grid/model.
    !      W3SETA    Subr  W3ADATMD Point to grid/model.
    !      W3SETO    Subr  W3ODATMD Point to grid/model.
    !      W3IOPE    Subr  W3IOPOMD Extracting point output.
    !      W3WAVE    Subr  W3WAVEMD Actual ave model routine.
    !      STRACE    Subr  W3SERVMD Subroutine tracing.
    !      EXTCDE    Subr    Id.    Abort program with exit code.
    !      WWTIME    Subr    Id.    System time in readable format.
    !      PRTIME    Subr    Id.    Profiling routine ( !/MPRF )
    !      STME21    Subr  W3TIMEMD Print date and time readable.
    !      DSEC21    Func    Id.    Difference between times.             
    !      TICK21    Subr    Id.    Advance time.
    !      WMSETM    Subr  WMMDATMD Point to grid/model.
    !      WMUPDT    Subr  WMUPDTMD Update input fields at driver level.
    !      WMIOBG    Subr  WMINIOMD Gather staged boundary data.
    !      WMIOBS    Subr    Id.    Stage boundary data.
    !      WMIOBF    Subr    Id.    Finalize WMIOBS.            ( !/MPI )
    !      WMIOHS    Subr    Id.    Stage high-to-low data.
    !      WMIOHG    Subr    Id.    Gather staged high-to-low data.
    !      WMIOHF    Subr    Id.    Finalize WMIOHS.            ( !/MPI )
    !      WMIOES    Subr    Id.    Stage same-rank data.
    !      WMIOEG    Subr    Id.    Gather staged same-rank data.
    !      WMIOEF    Subr    Id.    Finalize WMIOES.            ( !/MPI )
    !      WMIOPO    Subr  WMIOPOMD Unified point output.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      W3MLTI    Prog.   N/A    Multi-grid model driver.
    !      ....                     Any coupled model.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !     See formats 1000 and following, or escape locations 2000 and
    !     following.
    !
    !  7. Remarks :
    !
    !     - If no action is taken in the endless loop, an error is 
    !       assumed (code 2099). This should never take place in the
    !       default driver, but may be a problem in a coupled model.
    !     - If output is requested for the initial model time, TSYNC
    !       is set to TIME instead of (-1,0) for GRSTAT = 0. In this case
    !       input is updated, after which GRSTAT is set to 6 instead
    !       of 1. This assures that restarts do not impact model results
    !       by spurious double reconciliations.
    !
    !  8. Structure :
    !
    !     See source code and manual.
    !
    !  9. Switches :
    !
    !       !/S     Enable subroutine tracing.
    !       !/T     Enable test output.
    !       !/MPIT  Enable test output (use with !/MPI only).
    !       !/MPRF  Profiling.
    !
    !       !/SHRD, !/DIST, !/MPI
    !               Shared / distributed program model.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
    USE CONSTANTS
    !/
    USE W3ADATMD, ONLY: W3SETA
    !
    USE W3GDATMD, ONLY: W3SETG
    USE W3GDATMD, ONLY: DTMAX, NX, NY, MAPSTA, MAPST2
    !
    USE W3IDATMD, ONLY: INFLAGS1
    USE W3IOPOMD, ONLY: W3IOPE
    !
    USE W3ODATMD, ONLY: W3SETO, NOTYPE
    USE W3ODATMD, ONLY: FLOUT, TONEXT, DTOUT, TOLAST, IAPROC
    USE W3ODATMD, ONLY: NAPPNT, NOPTS, UNIPTS
    USE W3ODATMD, ONLY: NRQPO, IRQPO1 ! W3_MPI
    !
    USE W3SERVMD, ONLY: EXTCDE, WWTIME
    USE W3SERVMD, ONLY: STRACE ! W3_S
    !
    USE W3TIMEMD, ONLY: PRTIME ! W3_MPRF
    USE W3TIMEMD, ONLY: DSEC21, STME21, TICK21
    !
    USE W3WDATMD, ONLY: W3SETW
    !
    USE W3WAVEMD, ONLY: W3WAVE
    USE W3WDATMD, ONLY: TIME, VA
    !
    USE WMINIOMD, ONLY: WMIOBG, WMIOBS, WMIOBF, WMIOHG, WMIOHS
    USE WMINIOMD, ONLY: WMIOHF, WMIOEG, WMIOES, WMIOEF
    !
    USE WMIOPOMD, ONLY: WMIOPO
    !
    USE WMUPDTMD, ONLY: WMUPDT
    !
    USE WMMDATMD, ONLY: WMSETM
    USE WMMDATMD, ONLY: MDSO, MDSS, MDST, MDSE, IMPROC
    USE WMMDATMD, ONLY: NMPROC, NMPSCR, NMPERR, NMPTST, NMPLOG
    USE WMMDATMD, ONLY: STIME, ETIME, NMV, TMV, AMV, DMV
    USE WMMDATMD, ONLY: NRGRD, NRGRP, GRANK, INGRP, GRDHGH, GRDEQL
    USE WMMDATMD, ONLY: GRDLOW, TSYNC, TMAX, TOUTP, TDATA, GRSTAT
    USE WMMDATMD, ONLY: FLLSTL, FLLSTI, FLLSTR, DTRES, FLGHG1
    USE WMMDATMD, ONLY: FLGHG2, MAPMSK
    USE WMMDATMD, ONLY: MPI_COMM_MWAVE, MPI_COMM_GRD, MPI_COMM_BCT, CROOT, FBCAST ! W3_MPI
    USE WMMDATMD, ONLY: MDSP ! W3_MPRF
    use WMMDATMD, only: mpi_comm_grd_not_null
    !/
    IMPLICIT NONE
    !
#ifdef W3_MPI
    INCLUDE "mpif.h"
#endif
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    INTEGER, INTENT(IN)     :: TEND(2,NRGRD)
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
    INTEGER              :: J, JJ, I, JO, TPRNT(2), TAUX(2)
    INTEGER              :: II, JJJ, IX, IY, UPNEXT(2), UPLAST(2)
    INTEGER              :: DUMMY2(35)=0
    INTEGER              :: ILOOP                           ! W3_T
    INTEGER              :: IENT = 0                        ! W3_S
    INTEGER              :: IERR_MPI, NMPSCS                ! W3_MPI
    REAL                 :: DTTST, DTMAXI
    REAL                 :: PRFT0, PRFTN, PRFTS             ! W3_MPRF
    REAL(KIND=8)         :: get_memory
    CHARACTER(LEN=8)     :: WTIME
    CHARACTER(LEN=23)    :: MTIME
    LOGICAL              :: DONE, TSTAMP, FLAGOK, DO_UPT, FLG_O1, FLG_O2
    LOGICAL              :: FLAG                            ! W3_MPI
    INTEGER, ALLOCATABLE :: STATUS(:,:)
    LOGICAL, ALLOCATABLE :: FLSYNC(:), GRSYNC(:), TMSYNC(:), FLEQOK(:)
    LOGICAL, ALLOCATABLE :: PREGTB(:), PREGTH(:), PREGTE(:) ! W3_MPI
    !/
    !/ ------------------------------------------------------------------- /
    !
    !
    ! Set mpi_comm_grd_not_null
    !
    !
    if (w3_s_flag) then
       CALL STRACE (IENT, 'WMWAVE')
    end if
    !
    if (w3_mprf_flag) then
       CALL PRTIME ( PRFT0 )
    end if
    !              
    if (w3_o10_flag) then
       IF ( MDSS.NE.MDSO .AND. NMPSCR.EQ.IMPROC ) WRITE (MDSS,900)
    end if
    !
    ! 0.  Initializations ------------------------------------------------ /
    ! 0.a Initial testing
    !     Test GRSTAT
    !
    DO I=1, NRGRD
       IF ( ( GRSTAT(I).LT.0 .OR. GRSTAT(I).GT.7 ) .AND. GRSTAT(I).NE.99 ) GOTO 2000
       !
       !     Consistency of times for grids
       !
       IF ( TSYNC(1,I) .NE. -1 ) THEN
          DTTST  = DSEC21 ( TSYNC(:,I), TEND(:,I) )
          IF ( DTTST .LT. 0. ) GOTO 2001
       END IF
    END DO
    !
    !     Consistency of times within groups, set global TSYNC(:,0)
    !
    DO J=1, NRGRP
       DO JJ=2, INGRP(J,0)
          IF ( DSEC21(TSYNC(:,INGRP(J,1)),TSYNC(:,INGRP(J,JJ))).NE.0. .OR. &
               DSEC21(TEND(:,INGRP(J,1)),TEND(:,INGRP(J,JJ))).NE.0. ) then
             GOTO 2002
          end IF
       END DO
       IF ( GRANK(INGRP(J,1)).EQ.1 .AND. TSYNC(1,0).EQ.-1 ) then
          TSYNC(:,0) = TSYNC(:,INGRP(J,1))
       end IF
    END DO
    !
    !     Check if FLSYNC initialized
    !
    IF ( .NOT. ALLOCATED(FLSYNC) ) THEN
       ALLOCATE ( FLSYNC(NRGRD), GRSYNC(NRGRP), TMSYNC(NRGRD), FLEQOK(NRGRD) )
       if (w3_mpi_flag) then
          ALLOCATE ( PREGTB(NRGRD), PREGTH(NRGRD), PREGTE(NRGRD) )
       end if
       FLSYNC = .FALSE.
       GRSYNC = .FALSE.
       TMSYNC = .TRUE.
       FLEQOK = .FALSE.
       if (w3_mpi_flag) then
          PREGTB = .FALSE.
          PREGTH = .FALSE.
          PREGTE = .FALSE.
       end if
    END IF
    !
    ! 0.b Reset GRSTAT as needed
    !
    DO I=1, NRGRD
       CALL W3SETW ( I, MDSE, MDST )
       DTTST  = DSEC21 ( TIME, TEND(:,I) )
       IF ( GRSTAT(I).EQ.99 .AND. DTTST.GT.0. ) GRSTAT(I) = 0
    END DO
    !
    ! 0.c Other initializations
    !
    DTRES  = 0.
    if (w3_mpi_flag) then
       NMPSCS = NMPSCR
    end if
    !
    IF ( UNIPTS ) THEN
       CALL W3SETO ( 0, MDSE, MDST )
       UPNEXT = TONEXT(:,2)
       UPLAST = TOLAST(:,2)
       DO_UPT = .TRUE.
    ELSE
       UPNEXT(1) = -1
       UPNEXT(2) =  0
       DO_UPT = .FALSE.
    END IF
    !
    ! 0.d Output
    !
    IF ( MDSS.NE.MDSO .AND. NMPSCR.EQ.IMPROC ) THEN
       CALL WMPRNT ( MDSO, NRGRD, TSYNC(:,0), GRSTAT )
       CALL STME21 ( TSYNC(:,0), MTIME )
       CALL WWTIME ( WTIME )
       WRITE (MDSS,901) MTIME, WTIME, MINVAL(GRSTAT), MAXVAL(GRSTAT)
       TPRNT  = TSYNC(:,0)
       TSTAMP = .TRUE.
    ENDIF
    !
    if (w3_mpi_flag) then
       CALL MPI_BARRIER (MPI_COMM_MWAVE,IERR_MPI)
    end if
    !
#ifdef W3_MPRF
    CALL PRTIME ( PRFTN )
    WRITE (MDSP,990) PRFT0, PRFTN, get_memory()
#endif
    !
    ! 1.  Setting up loop structure -------------------------------------- /
    !
    if (w3_t_flag) then
       ILOOP  = 0
    end if
    !
    LOOP_OUTER: DO
       !
       IF ( MDSS.NE.MDSO .AND. NMPSCR.EQ.IMPROC .AND. DSEC21(TPRNT,TSYNC(:,0)).NE.0. ) THEN
          IF ( .NOT. TSTAMP ) WRITE (MDSS,*) 
          CALL WMPRNT ( MDSO, NRGRD, TSYNC(:,0), GRSTAT )
          CALL STME21 ( TSYNC(:,0), MTIME )
          CALL WWTIME ( WTIME )
          WRITE (MDSS,901) MTIME, WTIME, MINVAL(GRSTAT), MAXVAL(GRSTAT)
          !
          TPRNT  = TSYNC(:,0)
          TSTAMP = .TRUE.
       ENDIF
       !
       if (w3_t_flag) then
          ILOOP  = ILOOP + 1
          WRITE (MDST,9000) ILOOP, TSYNC(:,0)
          DO I=1, NRGRD
             CALL W3SETW ( I, MDSE, MDST )
             WRITE (MDST,9001) I, GRSTAT(I), TIME, TSYNC(:,I), TEND(:,I)
          END DO
          IF ( ILOOP .EQ. -1 ) CALL EXTCDE ( 508 )
       end if
       !
       DONE   = .FALSE.
       TPRNT  = TSYNC(:,0)
       !
       LOOP_J: DO J=1, NRGRP
          !
          if (w3_mpi_flag) then
             GRSYNC(J) = .FALSE.
             DO JJ=1, INGRP(J,0)
                I = INGRP(J,JJ)
                CALL WMSETM ( I, MDSE, MDST )
                GRSYNC(J) = GRSYNC(J) .OR. FBCAST
             END DO
          end if
          !
          LOOP_JJ: DO JJ=1, INGRP(J,0)
             I      = INGRP(J,JJ)
             CALL WMSETM ( I, MDSE, MDST )
             !
             if (w3_mpi_flag) then
                IF ( GRSTAT(I).EQ.0 ) TMSYNC(I) = .NOT. FBCAST
                IF ( FBCAST ) THEN
                   NMPSCR = CROOT
                ELSE
                   NMPSCR = NMPSCS
                END IF
             end if
             !
             ! 2.  Update input fields -------------------------------------------- /
             !     ( GRSTAT = 0 )
             !
             ! 2.a Check TDATA and finish step if data is still OK
             !
             IF ( (w3_shrd_flag .and. (GRSTAT(I) .EQ. 0 )) .or. &
                  (w3_mpi_flag  .and. ( GRSTAT(I).EQ.0 .AND. .NOT.FLSYNC(I) )) ) then
                !
                if (w3_t_flag) then
                   WRITE (MDST,9002) I, GRSTAT(I), ' '
                end if
                !
                IF ( TDATA(1,I) .EQ. -1 ) THEN
                   DTTST  = 0.  
                ELSE 
                   CALL W3SETW ( I, MDSE, MDST )
                   DTTST  = DSEC21 ( TIME , TDATA(:,I) ) 
                END IF
                if (w3_t_flag) then
                   WRITE (MDST,9020) DTTST
                end if
                !
                IF ( DTTST .GT. 0. ) THEN
                   GRSTAT(I) = 1
                   if (w3_t_flag) then
                      WRITE (MDST,9003) I, GRSTAT(I)
                   end if
                   DONE      = .TRUE.
                END IF
                !
             END IF
             !
             ! 2.b Update input and TDATA
             !
             IF ( (w3_shrd_flag .and. (GRSTAT(I) .EQ. 0 )) .or. &
                  (w3_mpi_flag  .and. (GRSTAT(I).EQ.0 .AND. .NOT.FLSYNC(I) .AND. mpi_comm_grd_not_null))) then
                !
                if (w3_mprf_flag) then
                   CALL PRTIME ( PRFT0 )
                end if
                IF ( DTTST .LE. 0 ) THEN
                   IF ( MDSS.NE.MDSO .AND. NMPSCR.EQ.IMPROC ) WRITE (MDSS,*)
                   TSTAMP = .FALSE.
                   CALL WMUPDT ( I, TDATA(:,I) )
                   if (w3_t_flag) then
                      WRITE (MDST,9021) TIME, TDATA(:,I), TEND(:,I)
                   end if
                END IF
                !
                ! 2.c Finish up if possible ( !/SHRD or .NOT. FBCAST or no update )
                !
                if (w3_shrd_flag) then
                   GRSTAT(I) = 1
                   DONE      = .TRUE.
                end if
                !
                IF ( .NOT. GRSYNC(J) ) THEN
                   if (w3_mpit_flag) then
                      WRITE (MDST,9902) I, GRSTAT(I),'NO SYNC FOR TDATA NEEDED'
                   end if
                   GRSTAT(I) = 1
                   DONE      = .TRUE.
                END IF
                !
#ifdef W3_MPRF
                CALL PRTIME ( PRFTN )
                WRITE (MDSP,991) PRFT0, PRFTN, get_memory(),'ST00', I
#endif
                !
             END IF
             !
             ! 2.d Synchronize in parts ( !/MPI )
             !
             if (w3_mpi_flag) then
                IF ( GRSTAT(I).EQ.0 .AND. GRSYNC(J) ) THEN
                   DONE   = .TRUE.
                   !
                   IF ( FLSYNC(I) ) THEN
                      if (w3_mpit_flag) then
                         WRITE (MDST,9902) I, GRSTAT(I), 'SYNCING TDATA'
                      end if
                      if (w3_mprf_flag) then
                         IF (FLSYNC(I)) CALL PRTIME ( PRFT0 )
                      end if
                      IF ( FBCAST ) CALL WMBCST ( TDATA(1,I), 2, I, NRGRD, 1 )
#ifdef W3_MPRF
                      IF (FLSYNC(I)) CALL PRTIME ( PRFTN )
                      IF (FLSYNC(I)) WRITE (MDSP,991) PRFT0, PRFTN, get_memory(), 'BCST',I
#endif
                      if (w3_mpit_flag) then
                         WRITE (MDST,9902) I, GRSTAT(I), 'SYNCING DONE'
                      end if
                      GRSTAT(I) = 1
                      FLSYNC(I) = .FALSE.
                      IF ( GRSYNC(J) ) CYCLE LOOP_JJ
                   ELSE
                      if (w3_mpit_flag) then
                         WRITE (MDST,9902) I, GRSTAT(I), 'CYCLE BEFORE SYNCING TDATA'
                      end if
                      FLSYNC(I) = .TRUE.
                      CYCLE LOOP_JJ
                   END IF
                   !
                END IF
             end if
             !
             ! 3.  Update data from lower ranked grids ---------------------------- /
             !     ( GRSTAT = 1 )
             !
             ! 3.a Skip for initial output only
             !
             IF ( GRSTAT(I) .EQ. 1 .AND. TSYNC(1,I) .NE. -1 ) THEN
                if (w3_t_flag) then
                   WRITE (MDST,9002) I, GRSTAT(I), 'FIRST PART'
                end if
                CALL W3SETW ( I, MDSE, MDST )
                DTTST  = DSEC21 ( TIME, TSYNC(:,I) )
                IF ( DTTST .EQ. 0. ) THEN
                   GRSTAT(I) = 7
                   if (w3_t_flag) then
                      WRITE (MDST,9003) I, GRSTAT(I)
                   end if
                   DONE      = .TRUE.
                END IF
             END IF
             !
             ! 3.b Normal processing
             !

             IF ( GRSTAT(I) .EQ. 1 ) THEN
                if (w3_t_flag) then
                   WRITE (MDST,9002) I, GRSTAT(I), 'SECOND PART'
                end if
                if (w3_mprf_flag) then
                   CALL PRTIME ( PRFT0 )
                end if
                !
                ! 3.b.1 Test if data is there
                !
                FLAGOK = .TRUE.
                CALL W3SETW ( I, MDSE, MDST )
                TAUX   = TIME
                DO JJJ=1, GRDLOW(I,0)
                   CALL W3SETW ( GRDLOW(I,JJJ), MDSE, MDST )
                   FLAGOK = FLAGOK .AND. DSEC21(TAUX,TIME).GT.0.       &
                        .AND. GRSTAT(GRDLOW(I,JJJ)).EQ.5
                END DO
                CALL W3SETW ( I, MDSE, MDST )
                !
                if (w3_t_flag) then
                   WRITE (MDST,9004) FLAGOK
                end if
                !
                ! 3.b.1 Get the data
                !
                if (w3_mpi_flag) then
                   IF ( .NOT.FLAGOK .AND. .NOT.PREGTB(I) ) THEN
                      IF ( mpi_comm_grd_not_null ) CALL WMIOBG (I,FLAG)
                      PREGTB(I) = .TRUE.
                   END IF
                end if
                !
                IF ( FLAGOK ) THEN
                   if (w3_shrd_flag) then
                      CALL WMIOBG ( I, FLAGOK )
                   end if
                   if (w3_mpi_flag) then
                      IF ( mpi_comm_grd_not_null ) CALL WMIOBG ( I )
                      PREGTB(I) = .FALSE.
                   end if
                   GRSTAT(I) = 2
                   DONE      = .TRUE.
                END IF
                !
#ifdef W3_MPRF
                CALL PRTIME ( PRFTN )
                WRITE (MDSP,991) PRFT0, PRFTN, get_memory(),'ST01', I
#endif
             END IF ! IF ( GRSTAT(I) .EQ. 1 )
             !
             ! 4.  Update model time step ----------------------------------------- /
             !     ( GRSTAT = 2 )
             !
             IF ( GRSTAT(I) .EQ. 2 ) THEN
                if (w3_t_flag) then
                   WRITE (MDST,9002) I, GRSTAT(I), ' '
                end if
                if (w3_mprf_flag) then
                   CALL PRTIME ( PRFT0 )
                end if
                !
                ! 4.a Check TMAX and update as necessary ( needs !/MPI synchronizaion )
                !
                CALL W3SETW ( I, MDSE, MDST )
                IF ( TMAX(1,I) .EQ. -1 ) THEN
                   TMAX(:,I) = TIME
                   DTTST  = 0.
                ELSE
                   DTTST  = DSEC21(TIME,TMAX(:,I))
                END IF
                !
                IF ( DTTST .LE. 0 ) THEN
                   CALL W3SETG ( I, MDSE, MDST )
                   DTMAXI = REAL(NINT(DTMAX+DTRES(I)+0.0001))
                   DTRES(I)= DTRES(I)+ DTMAX - DTMAXI
                   IF ( ABS(DTRES(I)) .LT. 0.001 ) DTRES(I) = 0.
                   TMAX(:,I) = TIME
                   CALL TICK21 ( TMAX(:,I), DTMAXI )
                   TAUX   = TMAX(:,I)
                   IF ( TDATA(1,I) .NE. -1 ) THEN
                      IF ( DSEC21(TDATA(:,I),TMAX(:,I)) .GT. 0 )    &
                           TMAX(:,I) = TDATA(:,I)
                   END IF
                   IF ( TOUTP(1,I) .NE. -1 ) THEN
                      IF ( DSEC21(TOUTP(:,I),TMAX(:,I)) .GT. 0 )    &
                           TMAX(:,I) = TOUTP(:,I)
                   END IF
                   IF ( UNIPTS ) THEN
                      IF ( DSEC21(UPNEXT,TMAX(:,I)) .GT. 0 )        &
                           TMAX(:,I) = UPNEXT(:)
                   END IF
                   if (w3_t_flag) then
                      WRITE (MDST,9040) TMAX(:,I), DTRES(I), TAUX,   &
                           TDATA(:,I), TOUTP(:,I), UPNEXT
                   end if
                   DONE   = .TRUE.
                   CYCLE LOOP_JJ
                ELSE
                   if (w3_t_flag) then
                      WRITE (MDST,9041) TMAX(:,I)
                   end if
                END IF ! IF ( DTTST .LE. 0 )
                !
                ! 4.b Lowest ranked grids, minimum of all TMAXes
                !
                if (w3_t_flag) then
                   WRITE (MDST,9042) GRANK(I)
                end if
                !
                IF ( GRANK(I) .EQ. 1 ) THEN
                   !
                   TAUX   = TMAX(:,I)
                   FLAGOK = .TRUE.
                   !
                   ! 4.b.1 Check if all grids have reached previous sync point
                   !
                   DO II=1, NRGRD
                      CALL W3SETW ( II, MDSE, MDST )
                      if ( w3_shrd_flag .and. ( TIME(1) .NE. -1 ) .or. &
                           w3_mpi_flag  .and. ( TIME(1) .NE. -1 .AND. mpi_comm_grd_not_null) ) then
                         IF ( DSEC21(TIME,TSYNC(:,0)) .NE. 0 ) THEN
                            FLAGOK = .FALSE.
                            EXIT
                         END IF
                      END IF
                   END DO
                   !
                   ! 4.b.2 Check availability of data
                   !
                   DO II=1, NRGRD
                      IF ( GRANK(II) .EQ. 1 ) THEN
                         IF ( TMAX(1,II) .EQ. -1 ) THEN
                            FLAGOK = .FALSE.
                            EXIT
                         ELSE
                            IF ( DSEC21 (TAUX,TMAX(:,II)) .LT. 0. ) TAUX   = TMAX(:,II)
                         END IF
                      END IF
                   END DO
                   !
                   CALL W3SETW ( I, MDSE, MDST )
                   FLAGOK = FLAGOK .AND. DSEC21(TIME,TAUX).GT.0.
                   !
                   ! 4.b.3 Update TSYNC for all grids
                   !
                   IF ( FLAGOK ) THEN
                      !
                      TSYNC(:,0) = TAUX
                      DO_UPT = .TRUE.
                      DO II=1, NRGRD
                         IF ( GRANK(II) .EQ. 1 ) THEN
                            TSYNC(:,II) = TAUX
                            IF ( GRSTAT(II) .EQ. 2 ) GRSTAT(II) = 3
                            if (w3_t_flag) then
                               IF ( GRSTAT(II) .EQ. 3 ) WRITE (MDST,9003) II, GRSTAT(II)
                            end if
                         END IF
                      END DO
                      DONE   = .TRUE.
#ifdef W3_MPRF
                      CALL PRTIME ( PRFTS )
                      WRITE (MDSP,992) PRFTS, PRFTS, get_memory(), 'TIME', TSYNC(:,0)
#endif
                      !
                      ! 4.b.4 Output
                      !
                      if (w3_t_flag) then
                         WRITE (MDST,9043) TSYNC(:,0)
                         WRITE (MDST,9045)
                         WRITE (MDST,9046) (II,TSYNC(:,II),II=0,NRGRD)
                      end if
                      !
                      ! 4.b.5 Skip computations so that all grids start processing
                      !       simultaneously.
                      !
#ifdef W3_MPRF
                      CALL PRTIME ( PRFTN )
                      WRITE (MDSP,991) PRFT0, PRFTN, get_memory(), 'ST02', I
#endif
                      if (w3_t_flag) then
                         IF ( INGRP(J,0) .GT. 1 ) WRITE (MDST,9006)
                      end if
                      IF ( INGRP(J,0) .GT. 1 ) GOTO 1111
                      !
                   END IF ! IF ( FLAGOK )
                   !
                   ! 4.c Other grids, logical from relations and TMAXes
                   !
                ELSE IF ( TSYNC(1,0) .NE. -1 ) THEN
                   !
                   TAUX   = TSYNC(:,0)
                   FLAGOK = .TRUE.
                   !
                   ! 4.c.1 Check availability of data within group
                   !       Time within group needs to be the same for load balancing.
                   !
                   DO JJJ=1, INGRP(J,0)
                      II     = INGRP(J,JJJ)
                      IF ( TMAX(1,II) .EQ. -1 ) THEN
                         FLAGOK = .FALSE.
                         EXIT
                      ELSE
                         IF ( DSEC21 (TAUX,TMAX(:,II)) .LT. 0. ) TAUX   = TMAX(:,II)
                      END IF
                   END DO
                   !
                   ! 4.c.2 Check with dependent lower rank grids ( TSYNC )
                   !
                   DO JJJ=1, GRDLOW(I,0)
                      II     = GRDLOW(I,JJJ)
                      IF ( TSYNC(1,II) .EQ. -1 ) THEN
                         FLAGOK = .FALSE.
                         EXIT
                      ELSE
                         IF ( DSEC21 (TAUX,TSYNC(:,II)) .LT. 0. ) TAUX   = TSYNC(:,II)
                      END IF
                   END DO
                   !
                   ! 4.c.3 Check with dependent higher rank grids ( TSYNC )
                   !       No check needed
                   !
                   ! 4.c.4 Final check against grid time
                   !
                   CALL W3SETW ( I, MDSE, MDST )
                   FLAGOK = FLAGOK .AND. DSEC21(TIME,TAUX).GT.0.
                   !
                   ! 4.c.5 Update TSYNC throughout group
                   !
                   IF ( FLAGOK ) THEN
                      !
                      DO JJJ=1, INGRP(J,0)
                         II     = INGRP(J,JJJ)
                         TSYNC(:,II) = TAUX
                         IF ( GRSTAT(II) .EQ. 2 ) GRSTAT(II) = 3
                         if (w3_t_flag) then
                            IF ( GRSTAT(II) .EQ. 3 ) WRITE (MDST,9003) II, GRSTAT(II)
                         end if
                      END DO
                      DONE   = .TRUE.
                      !
                      if (w3_t_flag) then
                         WRITE (MDST,9044) TSYNC(:,I), TAUX
                         WRITE (MDST,9045)
                         WRITE (MDST,9046) (II,TSYNC(:,II),II=0,NRGRD)
                      end if
                      !
                      ! 4.c.6 Skip computations so that all grids in group are advanced
                      !       simultaneously.
#ifdef W3_MPRF
                      CALL PRTIME ( PRFTN )
                      WRITE (MDSP,991) PRFT0, PRFTN, get_memory(),  'ST02', I
#endif
                      if (w3_t_flag) then
                         IF ( INGRP(J,0) .GT. 1 ) WRITE (MDST,9006)
                      end if
                      IF ( INGRP(J,0) .GT. 1 ) GOTO 1111
                      !
                   END IF ! IF ( FLAGOK )
                   !
                END IF ! IF ( GRANK(I) .EQ. 1 )
                !
             END IF ! IF ( GRSTAT(I) .EQ. 2 )
             !
             ! 5.  Run the wave model --------------------------------------------- /
             !     ( GRSTAT = 3 ) w3xdatmd data structures set in W3WAVE
             !
             ! 5.a Run model
             !
             if ( (w3_shrd_flag .and. ( GRSTAT(I) .EQ. 3 )) .or. &
                  (w3_mpi_flag .and.  ( GRSTAT(I) .EQ. 3 ) .and. mpi_comm_grd_not_null) ) then
                CALL W3SETW ( I, MDSE, MDST )
                TIME      = TSYNC(:,I)
                GRSTAT(I) = 4
                DONE      = .TRUE.
             ELSE IF ( w3_mpi_flag .and. GRSTAT(I).EQ.3 .and. mpi_comm_grd_not_null) then
                !
                if (w3_t_flag) then
                   WRITE (MDST,9002) I, GRSTAT(I), 'RUNNING MODEL'
                end if
                if (w3_mprf_flag) then
                   CALL PRTIME ( PRFT0 )
                end if
                !
                CALL WMSETM ( I, MDSE, MDST )
                CALL W3WAVE ( I, DUMMY2, TSYNC(:,I), .FALSE., .TRUE. )
                IF ( FLLSTL ) INFLAGS1(1) = .FALSE.
                IF ( FLLSTI ) INFLAGS1(4) = .FALSE.
                IF ( FLLSTR ) INFLAGS1(6) = .FALSE.
                !
                ! 5.b Stage data for grids with equal rank
                !
                if (w3_mpi_flag) then
                   CALL WMIOEF ( I )
                end if
                CALL WMIOES ( I )
                !
                ! 5.c Finish up 
                !
                GRSTAT(I) = 4
                DONE      = .TRUE.
                !
#ifdef W3_MPRF
                CALL PRTIME ( PRFTN )
                WRITE (MDSP,991) PRFT0, PRFTN, get_memory(),'ST03', I
#endif
                !
             END IF
             !
             ! 6.  Reconcile grids with same rank --------------------------------- /
             !     and stage data transfer to higher and lower ranked grids.
             !     ( GRSTAT = 4 )
             !
             IF ( GRSTAT(I) .EQ. 4 ) THEN
                if (w3_mprf_flag) then
                   CALL PRTIME ( PRFT0 )
                end if
                !
                ! 6.a Test against times and statuses of dependent grids
                !     Note: This is done per GROUP, not per local equal grid dependence
                !           Therefore, it is essential that sync times per group are
                !           equal (4.c.1) and that all equal grid dependences are a
                !           subset of groups (WMGEQL 5.d)
                !
                if (w3_t_flag) then
                   WRITE (MDST,9002) I, GRSTAT(I), 'FIRST PART'
                   WRITE (MDST,9005) FLEQOK(I)
                end if
                !
                ! 6.a.1 Check if sync point is reached
                !
                IF ( .NOT. FLEQOK(I) ) THEN
                   !
                   FLAGOK = .TRUE.
                   CALL W3SETW ( I, MDSE, MDST )
                   TAUX   = TIME
                   DO JJJ=1, INGRP(J,0)
                      CALL W3SETW ( INGRP(J,JJJ), MDSE, MDST )
                      FLAGOK = FLAGOK .AND. DSEC21(TAUX,TIME).EQ.0. .AND. GRSTAT(INGRP(J,JJJ)).EQ.4
                   END DO
                   CALL W3SETW ( I, MDSE, MDST )
                   if (w3_t_flag) then
                      WRITE (MDST,9004) FLAGOK
                   end if
                   !
                   ! 6.a.2 Point reached, set flag for all in group and cycle
                   !
                   IF ( FLAGOK ) THEN
                      DO JJJ=1, INGRP(J,0)
                         FLEQOK(INGRP(J,JJJ)) = .TRUE.
                         if (w3_t_flag) then
                            WRITE (MDST,9061) INGRP(J,JJJ), FLEQOK(INGRP(J,JJJ))
                         end if
                      END DO
                      DONE      = .TRUE.
#ifdef W3_MPRF
                      CALL PRTIME ( PRFTN )
                      WRITE (MDSP,991) PRFT0, PRFTN, get_memory(), 'ST04', I
#endif
                      !
                      if (w3_t_flag) then
                         IF ( INGRP(J,0) .GT. 1 ) WRITE (MDST,9006)
                      end if
                      IF ( INGRP(J,0) .GT. 1 ) GOTO 1111
                   END IF
                   !
                END IF ! IF ( .NOT. FLEQOK(I) )
                !
                ! 6.b Call gathering routine, reset FLEQOK and cycle
                !
                if (w3_mpi_flag) then
                   IF ( .NOT.FLEQOK(I) .AND. .NOT.PREGTE(I) ) THEN
                      IF ( mpi_comm_grd_not_null) CALL WMIOEG (I,FLAG)
                      PREGTE(I) = .TRUE.
                   END IF
                end if
                !
                IF ( FLEQOK(I) ) THEN
                   if (w3_shrd_flag) then
                      CALL WMIOEG ( I )
                   end if
                   if (w3_mpi_flag) then
                      IF ( mpi_comm_grd_not_null) CALL WMIOEG ( I )
                      PREGTE(I) = .FALSE.
                   end if
                   GRSTAT(I) = 5
                   FLEQOK(I) = .FALSE.
                   DONE      = .TRUE.
                END IF
                !
                ! 6.c Stage data
                !
                IF ( GRSTAT(I) .EQ. 5 ) THEN
                   !
                   if (w3_t_flag) then
                      WRITE (MDST,9002) I, GRSTAT(I)-1, 'SECOND PART'
                   end if
                   !
                   if (w3_shrd_flag) then
                      CALL WMIOBS ( I )
                   end if
                   if (w3_mpi_flag) then
                      IF ( mpi_comm_grd_not_null) then
                         CALL WMIOBF ( I )
                         CALL WMIOBS ( I )
                      END IF
                   end if
                   !
#ifdef W3_MPRF
                   CALL PRTIME ( PRFTN )
                   WRITE (MDSP,991) PRFT0, PRFTN, get_memory(), 'ST04', I
#endif
                   CYCLE LOOP_JJ
                   !
                END IF
                !
#ifdef W3_MPRF
                CALL PRTIME ( PRFTN )
                WRITE (MDSP,991) PRFT0, PRFTN, get_memory(),  'ST04', I
#endif
                !
             END IF
             !
             ! 7.  Reconcile with higher ranked grids ----------------------------- /
             !     ( GRSTAT = 5 )
             !
             !     This needs to be a little more complicated than with boundary
             !     data to assure proper logic in cases where data providing
             !     data does not get data back (e.g., as for the boundary grid
             !     in mww3_test_04)
             !
             IF ( GRSTAT(I) .EQ. 5 ) THEN
                if (w3_mprf_flag) then
                   CALL PRTIME ( PRFT0 )
                end if
                if (w3_t_flag) then
                   WRITE (MDST,9002) I, GRSTAT(I), ' '
                end if
                !
                ! 7.a Test against times and statuses of dependent grids
                !
                IF ( GRDHGH(I,0) .EQ. 0 ) THEN
                   GRSTAT(I) = 6
                   DONE      = .TRUE.
                ELSE
                   !
                   FLAGOK = .TRUE.
                   CALL W3SETW ( I, MDSE, MDST )
                   TAUX   = TIME
                   DO JJJ=1, GRDHGH(I,0)
                      CALL W3SETW ( GRDHGH(I,JJJ), MDSE, MDST )
                      IF ( .NOT. ( DSEC21(TAUX,TIME).EQ.0. .AND.      &
                           ( GRSTAT(GRDHGH(I,JJJ)).GE.7  .OR. &
                           GRSTAT(GRDHGH(I,JJJ)).LE.2 ) ) ) &
                           FLAGOK = .FALSE.
                   END DO
                   CALL W3SETW ( I, MDSE, MDST )
                   !
                   if (w3_t_flag) then
                      WRITE (MDST,9004) FLAGOK
                   end if
                   !
                   ! 7.b Call gathering routine
                   !
                   if (w3_mpi_flag) then
                      IF ( .NOT.FLAGOK .AND. .NOT.PREGTH(I) ) THEN
                         IF ( mpi_comm_grd_not_null )  CALL WMIOHG ( I, FLAG )
                         PREGTH(I) = .TRUE.
                      END IF
                   end if

                   IF ( FLAGOK ) THEN
                      if (w3_shrd_flag) then
                         CALL WMIOHG ( I, FLAGOK )
                      end if
                      if (w3_mpi_flag) then
                         if (mpi_comm_grd_not_null) CALL WMIOHG ( I )
                         PREGTH(I) = .FALSE.
                      end if
                      GRSTAT(I) = 6
                      DONE      = .TRUE.
                   END IF
                   !
                END IF

                !
                ! 7.c Stage data
                !
                if (w3_shrd_flag) then
                   IF ( GRSTAT(I) .EQ. 6 ) CALL WMIOHS ( I )
                end if
                !
                if (w3_mpi_flag) then
                   IF ( GRSTAT(I) .EQ. 6 .AND. mpi_comm_grd_not_null) then
                      CALL WMIOHF ( I )
                      CALL WMIOHS ( I )
                   END IF
                end if
                !
                if (w3_t_flag) then
                   IF (GRSTAT(I).EQ.6) WRITE(MDST,9003) I, GRSTAT(I)
                end if
#ifdef W3_MPRF
                CALL PRTIME ( PRFTN )
                WRITE (MDSP,991) PRFT0, PRFTN, get_memory(), 'ST05', I
#endif
             END IF
             !
             ! 8.  Perform data assimmilation ------------------------------------- /
             !     ( GRSTAT = 6 ) Placeholder only .....
             !
             IF ( GRSTAT(I) .EQ. 6 ) THEN
                if (w3_mprf_flag) then
                   CALL PRTIME ( PRFT0 )
                end if
                if (w3_t_flag) then
                   WRITE (MDST,9002) I, GRSTAT(I), ' '
                end if
                GRSTAT(I) = 7
#ifdef W3_MPRF
                CALL PRTIME ( PRFTN )
                WRITE (MDSP,991) PRFT0, PRFTN, get_memory(), 'ST06', I
#endif
                DONE      = .TRUE.
             END IF
             !
             ! 9.  Perform output ------------------------------------------------- /
             !     ( GRSTAT = 7 ) w3xdatmd data structures set in W3WAVE
             !
             !
             ! 9.a Check times and finish step if no output to be made
             !
             if ( (w3_shrd_flag .and. ( GRSTAT(I) .EQ. 7 )) .or.  &
                  (w3_mpi_flag .and. ( GRSTAT(I).EQ.7 .AND. .NOT.FLSYNC(I))) ) then
                !
                if (w3_t_flag) then
                   WRITE (MDST,9002) I, GRSTAT(I), ' '
                end if
                !
                IF ( TOUTP(1,I) .EQ. -1 ) THEN
                   DTTST  = 1.  
                ELSE 
                   CALL W3SETW ( I, MDSE, MDST )
                   DTTST  = DSEC21 ( TIME , TOUTP(:,I) ) 
                END IF
                if (w3_t_flag) then
                   WRITE (MDST,9090) DTTST
                end if
                FLG_O1 = DTTST .EQ. 0.
                !
                IF ( UNIPTS ) THEN
                   CALL W3SETW ( I, MDSE, MDST )
                   DTTST  = DSEC21 ( TIME , UPNEXT )
                   FLG_O2 = DTTST .EQ. 0.
                ELSE
                   FLG_O2 = .FALSE.
                END IF
                !
                IF ( .NOT.FLG_O1 .AND. .NOT.FLG_O2 ) THEN
                   GRSTAT(I) = 8
                   if (w3_t_flag) then
                      WRITE (MDST,9003) I, GRSTAT(I)
                   end if
                   DONE      = .TRUE.
                END IF
                !
             END IF
             !
             ! 9.b Perform output
             !
             if ( (w3_shrd_flag .and. ( GRSTAT(I) .EQ. 7 )) .or. &
                  (w3_mpi_flag  .and. ( GRSTAT(I) .EQ. 7 ) .and. mpi_comm_grd_not_null) )  THEN
                !
                !!/MPRF CALL PRTIME ( PRFT0 )
                !!/MPRF CALL WMWOUT ( I, NRGRD, 3 )
                !!/MPRF CALL PRTIME ( PRFTN )
                !!/MPRF WRITE (MDSP,991) PRFT0, PRFTN, get_memory(), 'BCST',I
                !
                if (w3_mprf_flag) then
                   CALL PRTIME ( PRFT0 )
                end if
                !
                IF ( FLG_O1 ) THEN
                   CALL W3SETG ( I, MDSE, MDST )
                   CALL WMSETM ( I, MDSE, MDST )
                   !
                   IF ( FLGHG1 .AND. .NOT.FLGHG2 .AND. GRDHGH(I,0).GT.0 ) THEN
                      MAPST2 = MAPST2 - 8*MAPMSK
                      MAPSTA = ABS(MAPSTA)
                      DO IX=1, NX
                         DO IY=1, NY
                            IF ( MAPST2(IY,IX) .GT. 0 ) then
                               MAPSTA(IY,IX) = - MAPSTA(IY,IX)
                            end IF
                         END DO
                      END DO
                      !
                   END IF
                   !
                   CALL W3WAVE ( I, DUMMY2, TSYNC(:,I), .FALSE. )
                   !
                   IF ( FLGHG1 .AND. .NOT.FLGHG2 .AND. GRDHGH(I,0).GT.0 ) THEN
                      MAPST2 = MAPST2 + 8*MAPMSK
                      MAPSTA = ABS(MAPSTA)
                      DO IX=1, NX
                         DO IY=1, NY
                            IF ( MAPST2(IY,IX) .GT. 0 ) then
                               MAPSTA(IY,IX) = - MAPSTA(IY,IX)
                            end IF
                         end DO
                      END DO
                   END IF
                   !
                   IF ( FLLSTL ) INFLAGS1(1) = .FALSE.
                   IF ( FLLSTI ) INFLAGS1(4) = .FALSE.
                   IF ( FLLSTR ) INFLAGS1(6) = .FALSE.
                   !
                   ! 9.c Update TOUPT
                   !
                   TOUTP(1,I) = -1
                   TOUTP(2,I) =  0
                   !
                   DO JO=1, NOTYPE
                      IF ( .NOT.FLOUT(JO) ) CYCLE
                      IF ( TOUTP(1,I) .EQ. -1 ) THEN
                         TOUTP(:,I) = TONEXT(:,JO)
                      ELSE
                         DTTST = DSEC21 ( TOUTP(:,I) , TONEXT(:,JO) )
                         IF (DTTST.LT.0.) TOUTP(:,I) = TONEXT(:,JO)
                      ENDIF
                   END DO
                   ! CHECKPOINT
                   JO=8
                   IF ( .NOT.FLOUT(JO) ) CYCLE
                   IF ( TOUTP(1,I) .EQ. -1 ) THEN
                      TOUTP(:,I) = TONEXT(:,JO)
                   ELSE
                      DTTST = DSEC21 ( TOUTP(:,I) , TONEXT(:,JO) )
                      IF (DTTST.LT.0.) TOUTP(:,I) = TONEXT(:,JO)
                   ENDIF
                   ! END CHECKPOINT
                   !
                   if (w3_t_flag) then
                      WRITE (MDST,9091) TOUTP(:,I)
                   end if
                   !
                END IF !  IF ( FLG_O1 )

                !
                ! 9.d Process unified point output for selected grid
                !
                IF ( UNIPTS ) THEN
                   IF ( FLG_O2 ) THEN
                      CALL W3SETO ( I, MDSE, MDST )
                      !
#ifdef W3_MPI
                      IF ( NRQPO.NE.0 ) CALL MPI_STARTALL ( NRQPO, IRQPO1, IERR_MPI )
#endif
                      !
                      IF ( NOPTS.NE.0 .AND. IAPROC.EQ.NAPPNT ) THEN
                         CALL W3SETG ( I, MDSE, MDST )
                         CALL W3SETA ( I, MDSE, MDST )
                         CALL W3IOPE ( VA )
                      END IF
                      !
#ifdef W3_MPI
                      IF ( NRQPO .NE. 0 ) THEN
                         ALLOCATE ( STATUS(MPI_STATUS_SIZE,NRQPO) )
                         CALL MPI_WAITALL ( NRQPO, IRQPO1, STATUS, IERR_MPI )
                         DEALLOCATE ( STATUS )
                      END IF
#endif
                      !
                      if (w3_t_flag) then
                         WRITE (MDST,9092) NOPTS
                      end if
                      !
                   END IF
                   !
                END IF !  IF ( UNIPTS )
                !
#ifdef W3_MPRF
                CALL PRTIME ( PRFTN )
                WRITE (MDSP,991) PRFT0, PRFTN,get_memory(),'ST07', I
#endif
                !
                ! 9.e Update TOUPT outside communicator
                !
             ELSE IF ( w3_mpi_flag ) then
                !
                if ( FLG_O1 ) THEN
                   CALL W3SETO ( I, MDSE, MDST )
                   CALL W3SETW ( I, MDSE, MDST )
                   !
                   TIME       = TOUTP(:,I)
                   TOUTP(1,I) = -1
                   TOUTP(2,I) =  0
                   ! 
                   DO JO=1, NOTYPE
                      !
                      IF ( FLOUT(JO) ) THEN
                         DO
                            DTTST = DSEC21 ( TIME, TONEXT(:,JO) )
                            IF ( DTTST .LE. 0. ) THEN
                               CALL TICK21 ( TONEXT(:,JO), DTOUT(JO) )
                               DTTST = DSEC21 ( TONEXT(:,JO), TOLAST(:,JO) )
                               IF ( DTTST .LT. 0. ) THEN
                                  FLOUT(JO) = .FALSE.
                                  EXIT
                               END IF
                            ELSE
                               EXIT
                            END IF
                         END DO
                      END IF ! IF ( FLOUT(JO) )
                      !
                      IF ( .NOT.FLOUT(JO) ) CYCLE
                      IF ( TOUTP(1,I) .EQ. -1 ) THEN
                         TOUTP(:,I) = TONEXT(:,JO)
                      ELSE
                         DTTST = DSEC21 ( TOUTP(:,I) , TONEXT(:,JO) )
                         IF (DTTST.LT.0.) TOUTP(:,I) = TONEXT(:,JO)
                      ENDIF
                      !
                   END DO ! DO JO=1, NOTYPE
                   !
                   ! Checkpoint
                   !
                   JO=8
                   !
                   IF ( FLOUT(JO) ) THEN
                      DO
                         DTTST = DSEC21 ( TIME, TONEXT(:,JO) )
                         IF ( DTTST .LE. 0. ) THEN
                            CALL TICK21 ( TONEXT(:,JO), DTOUT(JO) )
                            DTTST = DSEC21 ( TONEXT(:,JO), TOLAST(:,JO) )
                            IF ( DTTST .LT. 0. ) THEN
                               FLOUT(JO) = .FALSE.
                               EXIT
                            END IF
                         ELSE
                            EXIT
                         END IF
                      END DO
                   END IF
                   !
                   IF ( .NOT.FLOUT(JO) ) CYCLE
                   IF ( TOUTP(1,I) .EQ. -1 ) THEN
                      TOUTP(:,I) = TONEXT(:,JO)
                   ELSE
                      DTTST = DSEC21 ( TOUTP(:,I) , TONEXT(:,JO) )
                      IF (DTTST.LT.0.) TOUTP(:,I) = TONEXT(:,JO)
                   ENDIF
                   !
                   ! End Checkpoint
                   if (w3_mpit_flag) then
                      WRITE (MDST,9991) TOUTP(:,I)
                   end if
                END IF ! if ( FLG_O1 )
                !
                ! 9.f Finish up
                !
                GRSTAT(I) = 8
                DONE      = .TRUE.
                !
             END IF ! if ( (w3_shrd_flag .and. ( GRSTAT(I) .EQ. 7 )) .or....
             !
             ! 10. Go to next time step ------------------------------------------- /
             !     ( GRSTAT = 8 ) ( 9 added for diagnostic output only ... )
             !                    ( Unified point output and synchronization added )
             !
             IF ( GRSTAT(I) .EQ. 8 ) THEN
                !
                if (w3_t_flag) then
                   WRITE (MDST,9002) I, GRSTAT(I), ' '
                end if
                !
                ! 10.a Processing unified point output
                !
                IF ( UNIPTS .AND. DO_UPT ) THEN
                   CALL W3SETW ( I, MDSE, MDST )
                   FLAGOK = DSEC21 ( TIME, UPNEXT ) .EQ. 0.
                   if (w3_t_flag) then
                      WRITE (MDST,9095) FLAGOK
                   end if
                ELSE
                   FLAGOK = .FALSE.
                END IF
                !
                IF ( FLAGOK ) THEN
                   !
                   DO II=1, NRGRD
                      CALL W3SETW ( II, MDSE, MDST )
                      FLAGOK = FLAGOK .AND. GRSTAT(II).EQ.8 .AND.     &
                           DSEC21(TIME,UPNEXT).EQ.0.
                   END DO
                   if (w3_t_flag) then
                      WRITE (MDST,9096) FLAGOK
                   end if
                   !
                   IF ( FLAGOK ) THEN
                      !
                      if (w3_mprf_flag) then
                         CALL PRTIME ( PRFT0 )
                      end if
                      CALL WMIOPO ( UPNEXT )
                      DO_UPT = .FALSE.
                      !
                      CALL W3SETO ( 0, MDSE, MDST )
                      CALL TICK21 ( UPNEXT, DTOUT(2) )
                      IF ( DSEC21(UPNEXT,UPLAST) .GE. 0. ) THEN
                         TONEXT(:,2) = UPNEXT
                      ELSE
                         UNIPTS = .FALSE.
                         UPNEXT(1) = -1
                         UPNEXT(2) =  0
                      END IF
                      !
                      DO II=1, NRGRD
                         CALL W3SETW ( II, MDSE, MDST )
                         DTTST  = DSEC21 ( TIME, TEND(:,II) )
                         IF ( DTTST .GT. 0. ) THEN
                            GRSTAT(II) = 9
                         ELSE IF ( DTTST .EQ. 0 ) THEN
                            GRSTAT(II) = 99
                         END IF
                         TSYNC(1,II) = -1
                         TSYNC(2,II) =  0
                         if (w3_t_flag) then
                            IF ( I .NE. II ) WRITE (MDST,9003) II, GRSTAT(II)
                         end if
                      END DO
                      !
                      DONE      = .TRUE.
#ifdef W3_MPRF
                      CALL PRTIME ( PRFTN )
                      WRITE (MDSP,991) PRFT0, PRFTN, get_memory(), 'UPTS',I
#endif
                   END IF
                   !
                ELSE
                   FLAGOK = .TRUE.
                END IF
                !
                ! 10.b Regular processing
                !
                IF ( FLAGOK ) THEN
                   CALL W3SETW ( I, MDSE, MDST )
                   DTTST  = DSEC21 ( TIME, TEND(:,I) )
                   IF ( DTTST .GT. 0. ) THEN
                      GRSTAT(I) = 9
                      DONE      = .TRUE.
                   ELSE IF ( DTTST .EQ. 0 ) THEN
                      GRSTAT(I) = 99
                      DONE      = .TRUE.
                   END IF
                   if (w3_t_flag) then
                      WRITE (MDST,9003) I, GRSTAT(I)
                   end if
                END IF
                !
                IF ( GRSTAT(I).EQ.9 .OR. GRSTAT(I).EQ.99 ) THEN
                   TSYNC(1,I) = -1
                   TSYNC(2,I) =  0
                END IF
                !
             END IF
             !
             ! ... End of loops started in 1. ------------------------------------- /
             !
          END DO LOOP_JJ
          !
1111      CONTINUE
          !
       END DO LOOP_J
       !
       if (w3_mpi_flag) then
          NMPSCR = NMPSCS
       end if
       IF ( MDSS.NE.MDSO .AND. NMPSCR.EQ.IMPROC ) then
          CALL WMPRNT ( MDSO, NRGRD, TSYNC(:,0), GRSTAT )
       end IF
       !
       DO I=1, NRGRD
          IF ( GRSTAT(I) .EQ. 9 ) GRSTAT(I) = 0
       END DO
       !
       IF ( .NOT. DONE ) GOTO 2099
       IF ( MINVAL(GRSTAT) .EQ. 99 ) EXIT LOOP_OUTER
    END DO LOOP_OUTER
    !
    !    End of routine -------------------------------------------------- /

    IF ( MDSS.NE.MDSO .AND. NMPSCR.EQ.IMPROC ) THEN
       CALL WWTIME ( WTIME )
       WRITE (MDSS,902) WTIME
    ENDIF
    !              
    if (w3_mpi_flag) then
       DO I=1, NRGRD
          CALL WMSETM ( I, MDSE, MDST )
          IF ( mpi_comm_grd_not_null) then
             CALL WMIOBF ( I )
             CALL WMIOHF ( I )
             CALL WMIOEF ( I )
          END IF
       END DO
    end if
    !
    if (w3_o10_flag) then
       IF ( MDSS.NE.MDSO .AND. NMPSCR.EQ.IMPROC ) WRITE (MDSS,999)
    end if
    !
    if (w3_t_flag) then
       WRITE (MDST,9100)
    end if
    !
    RETURN
    !
    ! Escape locations
    !
2000 CONTINUE
    IF ( IMPROC .EQ. NMPERR ) WRITE (MDSE,1000) I, GRSTAT(I)
    CALL EXTCDE ( 2000 )
    RETURN
    !
2001 CONTINUE
    IF ( IMPROC .EQ. NMPERR ) WRITE (MDSE,1001) I, TSYNC(:,I), TEND(:,I)
    CALL EXTCDE ( 2001 )
    RETURN
    !
2002 CONTINUE
    IF ( IMPROC .EQ. NMPERR ) WRITE (MDSE,1002) J, INGRP(J,1),      &
         INGRP(J,JJ), TSYNC(:,INGRP(J,1)), TSYNC(:,INGRP(J,JJ)),    &
         TEND(:,INGRP(J,1)), TEND(:,INGRP(J,JJ))
    CALL EXTCDE ( 2002 )
    RETURN
    !
2099 CONTINUE
    IF ( IMPROC .EQ. NMPERR ) WRITE (MDSE,1099)
    CALL EXTCDE ( 2099 )
    RETURN
    !
    ! Formats
    !
900 FORMAT (           ' ========== STARTING WAVE MODEL (WMWAVE) ==========', &
                       '============================'/)
901 FORMAT (           '  MWW3 calculating for ',A,' at ',A,'  status [',      &
         I2,           '-',I2,']')
902 FORMAT (           '  MWW3 reached the end of the computation loop at ',A)
990 FORMAT (1X,3F12.3, ' WMWAVE INIT') ! W3_MPRF
991 FORMAT (1X,3F12.3, ' WMWAVE ',A4,I6)
992 FORMAT (1X,3F12.3, ' WMWAVE ',A4,I9.8,I7.6)
999 FORMAT (/          ' ========== END OF WAVE MODEL (WMWAVE) ============', &
                       '============================'/)
    !
1000 FORMAT (/         ' *** WAVEWATCH III ERROR IN WMWAVE : *** '/           &
                       '     GRID',I3,' HAS ILLEGAL GRSTAT :',I8/)
    !
1001 FORMAT (/         ' *** WAVEWATCH III ERROR IN WMWAVE : *** '/           &
                       '     GRID',I3,' HAS ILLEGAL TSYNC / TEND '/           &
                       '     TSYNC :',I9.8,I7.6/                              &
                       '     TEND  :',I9.8,I7.6/)
    !
1002 FORMAT (/         ' *** WAVEWATCH III ERROR IN WMWAVE : *** '/           &
                       '     GROUP',I3,' HAS INCOMPATIBLE TIMES ',            &
                       'IN GRIDS ',I3,' AND ',I3/                             &
                       '     TSYNC :',I9.8,I7.6,1X,I9.8,I7.6/                 &
                       '     TEND  :',I9.8,I7.6,1X,I9.8,I7.6/)
    !
    ! Note: This 1099 error can occur when multi-grid time steps are not
    !       compatible.
1099 FORMAT (/         ' *** WAVEWATCH III ERROR IN WMWAVE : *** '/           &
                       '     ABORT FOR POSSIBLE ENDLESS LOOP '/)
    !
9000 FORMAT (          ' TEST WMWAVE : LOOP',I8,' ======================', & ! W3_T
                       '===== (',I9.8,I7.6,' ) =='/ &
                       '               GRID, GRSTAT, TIME, TSYNC, TEND')
9001 FORMAT (          '              ',I3,I3,3(I10.8,I7.6))
9002 FORMAT (          ' TEST WMWAVE : PROCESSING GRID',I3,                &
                       ' STATUS',I3,' ',A)
9902 FORMAT (          ' MPIT WMWAVE : PROCESSING GRID',I3,             & ! W3_MPIT
                       ' STATUS',I3,' ',A)
9003 FORMAT (          ' TEST WMWAVE : GRID',I3,' STATUS RESET TO',I3) ! W3_T
9004 FORMAT (          ' TEST WMWAVE : FLAGOK = ',L1)
9005 FORMAT (          ' TEST WMWAVE : FLEQOK = ',L1)
9006 FORMAT (          ' TEST WMWAVE : CYCLE GROUP')
    !
9020 FORMAT (          ' TEST WMWAVE : DTTST ',E10.3) ! W3_T
9021 FORMAT (          ' TEST WMWAVE : TIME  :',I10.8,I7.6/                &
                       '               TDATA :',I10.8,I7.6/                &
                       '               TEND  :',I10.8,I7.6)
    !
9040 FORMAT (          ' TEST WMWAVE : TMAX  :',I10.8,I7.6,F8.2/           & ! W3_T
                       '               DTMAX :',I10.8,I7.6/                &
                       '               TDATA :',I10.8,I7.6/                &
                       '               TOUTP :',I10.8,I7.6/                &
                       '               UPNEXT:',I10.8,I7.6)
9041 FORMAT (          ' TEST WMWAVE : TMAX  :',I10.8,I7.6)
9941 FORMAT (          ' MPIT WMWAVE : TMAX  :',I10.8,I7.6) ! W3_MPIT
9042 FORMAT (          ' TEST WMWAVE : GRANK :',I4,' FOR GRSTAT = 2') ! W3_T
9043 FORMAT (          ' TEST WMWAVE : GLOBAL TSYNC :',I10.8,I7.6)
9044 FORMAT (          ' TEST WMWAVE : LOCAL TSYNC :',I10.8,I7.6,          &
                       ' (',I8.8,I7.6,')')
9045 FORMAT (          ' TEST WMWAVE : GRID TSYNC')
9046 FORMAT (          '             ',I5,I10.8,I7.6)
    !
9061 FORMAT (          '               GRID',I4,',  FLEQOK = ',L1) ! W3_T
    !
9090 FORMAT (          ' TEST WMWAVE : DTTST ',E10.3) ! W3_T
9091 FORMAT (          ' TEST WMWAVE : NEXT TOUTP :',I10.8,I7.6)
9991 FORMAT (          ' MPIT WMWAVE : NEXT TOUTP :',I10.8,I7.6) ! W3_MPIT
9092 FORMAT (          ' TEST WMWAVE : UNIFIED POINT OUTPUT PREP DONE',I6) ! W3_T
    !
9095 FORMAT (          ' TEST WMWAVE : UNIFIED POINT OUTPUT, FLAGOK = ',L1) ! W3_T
9096 FORMAT (          '               ALL GRIDS, FLAGOK = ',L1)
    !
9100 FORMAT (          ' TEST WMWAVE : LOOP DONE  ======================', & ! W3_T
                       '==============================')
    !/
    !/ End of WMWAVE ----------------------------------------------------- /
    !/
  END SUBROUTINE WMWAVE
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Print out action table in the log file log.\ mww3.
  !>
  !> @param[in] MDSO
  !> @param[in] NRGRD Number of grids.
  !> @param[in] TSYNC Synchronization time.
  !> @param[in] GRSTAT Status array per grid.
  !>    
  !> @author H. L. Tolman  @date 22-Feb-2005
  !>     
  SUBROUTINE WMPRNT ( MDSO, NRGRD, TSYNC, GRSTAT )
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           H. L. Tolman            |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         22-Feb-2005 |
    !/                  +-----------------------------------+
    !/
    !/    22-Feb-2005 : Origination.                        ( version 3.07 )
    !/
    !  1. Purpose :
    !
    !     Print out action table in the log file log.ww3m
    !
    !  2. Method :
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       NRGRD   Int.   I   Number of grids.
    !       TSYN    I.A.   I   Synchronization time.
    !       GRSTAT  I.A.   I   Status array per grid.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD Subroutine tracing.
    !      STME21    Subr. W3TIMEMD Print date and time readable.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      WMWAVE    Subr. WMWAVEMD Multi-grid wave model routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !     None.
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    !     !/S     Enable subroutine tracing.
    !     !/T     Test output.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
    USE W3SERVMD, ONLY: STRACE ! W3_S
    USE W3TIMEMD, ONLY: STME21
    !/
    IMPLICIT NONE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    INTEGER, INTENT(IN)     :: MDSO, NRGRD, TSYNC(2), GRSTAT(NRGRD)
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
    INTEGER, PARAMETER      :: IW = 15
    INTEGER                 :: I, I0, IN
    INTEGER, SAVE           :: IENT = 0 ! W3_S
    INTEGER, SAVE           :: IDLAST(2)
    LOGICAL, SAVE           :: FIRST = .TRUE.
    CHARACTER(LEN=23)       :: IDTIME
    CHARACTER(LEN=3)        :: STR(IW), LNE(IW)
    !/
    !/ ------------------------------------------------------------------- /
    !/
    if (w3_s_flag) then
       CALL STRACE (IENT, 'WMPRNT')
    end if
    !
    DO I=1, IW
       LNE(I) = '---'
    END DO
    !
    IF ( FIRST ) THEN
       WRITE (MDSO,900) NRGRD, LNE, '-+'
       FIRST  = .FALSE.
       IDLAST(1) = -1
       IDLAST(2) =  0
    ELSE
       BACKSPACE (MDSO)
    END IF
    !
    CALL STME21 ( TSYNC, IDTIME )
    !
    DO I=1, MIN(IW,NRGRD)
       WRITE (STR(I),'(I3)') GRSTAT(I)
    END DO
    DO I=1+MIN(IW,NRGRD), IW
       STR(I) = '   '
    END DO
    !
    IF ( IDLAST(1).EQ.TSYNC(1) .AND. IDLAST(2).EQ.TSYNC(2) ) THEN
       if (w3_o11_flag) then
          WRITE (MDSO,903) STR, ' |'
       end if
    ELSE IF ( IDLAST(1) .EQ. TSYNC(1) ) THEN
       WRITE (MDSO,902) IDTIME(12:19), STR, ' |'
    ELSE
       WRITE (MDSO,901) IDTIME(01:19), STR, ' |'
    END IF
    IDLAST = TSYNC
    !
    IF ( NRGRD .GT. IW ) THEN
       I0     = 1
       IN     = IW
       DO
          I0     = I0 + IW
          IN     = IN + IW
          DO I=I0, MIN(IN,NRGRD)
             WRITE (STR(I-I0+1),'(I3)') GRSTAT(I)
          END DO
          DO I=1+MIN(IN,NRGRD), IN
             STR(I-I0+1) = '   '
          END DO
          WRITE (MDSO,903) STR, ' |'
          IF ( IN .GE. NRGRD ) EXIT
       END DO
    END IF
    !
    WRITE (MDSO,904) LNE, '-+'
    !
    RETURN
    !
    ! Formats
    !
900 FORMAT (1X,' Time (sync rank 1)  | Status for',I3,' grids'/     &
         1X,'---------------------+',16A)
901 FORMAT (2X,A19,' |',16A)
902 FORMAT (2X,11X,A8,' |',16A)
903 FORMAT (21X,' |',16A)
904 FORMAT (1X,'---------------------+',16A)
    !/
    !/ End of WMPRNT ----------------------------------------------------- /
    !/
  END SUBROUTINE WMPRNT
  !/ ------------------------------------------------------------------- /
  !>     
  !> @brief Non-blocking broadcast for integer arrays.
  !>
  !> @details Non-blocking broadcast, initially for times only, 
  !>  but made for any integer array. Sending data from first process
  !>  in the model communicator to all processes that are in the overall 
  !>  communicator but not in the model communicator.
  !> 
  !>  Standard send and receives using defined communicator. Send
  !>  form first processor in communicator.    
  !>
  !> @param[inout] DATA Data to be send/received.
  !> @param[in] NR Size of array.
  !> @param[in] IMOD Model number.
  !> @param[in] NMOD Number of models.
  !> @param[in] ID ID number, used with NMOD for ITAG.
  !> @author H. L. Tolman  @date 02-Feb-2007
  !>      
  SUBROUTINE WMBCST ( DATA, NR, IMOD, NMOD, ID )
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           H. L. Tolman            |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         02-Feb-2007 !
    !/                  +-----------------------------------+
    !/
    !/    02-Feb-2007 : Origination.                        ( version 3.10 )
    !/
    !  1. Purpose :
    !
    !     Non-blocking broadcast, initially for times only, but made for
    !     any integer array. Sending data from first process in the 
    !     model cummunicator to all processes that are in the overall 
    !     communicator but not in the model communicator.
    !
    !  2. Method :
    !
    !     Standard send and receives using defined communicator. Send
    !     form first processor in communicator.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       DATA    I.A.  I/O  Data to be send/received.
    !       NR      Int.   I   Size of array.
    !       IMOD    Int.   I   Model number.
    !       NMOD    Int.   I   Number of models.
    !       ID      Int.   I   ID number, used with NMOD for ITAG.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Sur.  W3SERVMD Subroutine tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      WMWAVE    Subr. WMWAVEMD Multi-grid wave model routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !  9. Switches :
    !
    !     !/S    Enable subroutine tracing.
    !     !/MPIT Enable test output
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
    !
    USE WMMDATMD, ONLY: MDST, MTAGB, IMPROC, NMPROC, ALLPRC ! W3_MPI 
    USE WMMDATMD, ONLY: CROOT, MPI_COMM_MWAVE               ! W3_MPI
    use WMMDATMD, only: mpi_comm_grd_not_null               ! W3_MPI
    USE W3SERVMD, ONLY: STRACE                              ! W3_S
    !
    IMPLICIT NONE
    !
#ifdef W3_MPI
    INCLUDE "mpif.h"
#endif
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    INTEGER, INTENT(IN)     :: NR, IMOD, NMOD, ID
    INTEGER, INTENT(INOUT)  :: DATA(NR)
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
#ifdef W3_MPI
    INTEGER :: STATUS(MPI_STATUS_SIZE)
#endif
    INTEGER :: ITAG, IP, IERR_MPI ! W3_MPI
    INTEGER :: IENT = 0           ! W3_S
    !
    ! -------------------------------------------------------------------- /
    ! 0.  Initializations
    !
    if (w3_s_flag) then
       CALL STRACE (IENT, 'WMBCST')
    end if

    if (w3_mpi_flag) then
       ITAG   = MTAGB + IMOD + ID*NMOD
    end if
    !
    !!/MPIT      WRITE (MDST,9000) IMOD, NMOD, ID, ITAG-MTAGB
    !!/MPIT      WRITE (MDST,9001) IMPROC, NMPROC
    !!/MPIT      WRITE (MDST,9002) ALLPRC(:,IMOD)
    !
    ! -------------------------------------------------------------------- /
#ifdef W3_MPI
    if ( ALLPRC(IMPROC,IMOD) .EQ. 1 ) THEN
       ! 1.  Processor to send data from
       !
       DO IP=1, NMPROC
          IF ( ALLPRC(IP,IMOD) .EQ. 0 ) THEN
             !!/MPIT WRITE (MDST,9010) ALLPRC(IMPROC,IMOD), IP
             CALL MPI_SEND ( DATA, NR, MPI_INTEGER, IP-1, ITAG, MPI_COMM_MWAVE, IERR_MPI )
          end IF
       END DO
       !
       ! -------------------------------------------------------------------- /
       ! 2.  Processor to receive data at
       !
    else if ( ALLPRC(IMPROC,IMOD) .EQ. 0 ) THEN
       !!/MPIT WRITE (MDST,9020) ALLPRC(IMPROC,IMOD), CROOT
       CALL MPI_RECV ( DATA, NR, MPI_INTEGER, CROOT-1, ITAG, MPI_COMM_MWAVE, STATUS, IERR_MPI )
       !
       ! -------------------------------------------------------------------- /
       ! 3.  Processor with no action
       !
       !!/MPIT ELSE
       !!/MPIT   WRITE (MDST,9030) ALLPRC(IMPROC,IMOD)
    end if
#endif
    !
    RETURN
    !
    ! Formats
    !
9000 FORMAT ( ' TEST WMBCST : INPUTS :',4I4) ! W3_MPIT
9001 FORMAT ( ' TEST WMBCST : IMPROC, NMPROC:',2I5,'  ALLPRC :')
9002 FORMAT (14X,13I5)
    !
9010 FORMAT ( ' TEST WMBCST : IAPROC =',I5,'    SENDING TO ',I5) ! W3_MPIT
    !
9020 FORMAT ( ' TEST WMBCST : IAPROC =',I5,                    & ! W3_MPIT
         '    RECEIVING FROM ',I5)
    !
9030 FORMAT ( ' TEST WMBCST : IAPROC =',I5,'    NO ACTION') ! W3_MPIT
    !/
    !/ End of WMBCST ----------------------------------------------------- /
    !/
  END SUBROUTINE WMBCST
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Non-blocking broadcast using dummy parameter to have output.
  !>
  !> @details Processes wait for computations on first node to be finished.
  !>  Needed for profiling purposes only.
  !>
  !>  Standard send and recieves using defined communicator. Send
  !>  form first processor in communicator.
  !>
  !> @param[in] IMOD Model number.
  !> @param[in] NMOD Number of models.
  !> @param[in] ID ID number, used with NMOD for ITAG.
  !>      
  !> @author H. L. Tolman  @date 21-Jun-2007
  !>      
  SUBROUTINE WMWOUT ( IMOD, NMOD, ID )
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           H. L. Tolman            |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         21-Jun-2007 !
    !/                  +-----------------------------------+
    !/
    !/    21-Jun-2007 : Origination.                        ( version 3.11 )
    !/
    !  1. Purpose :
    !
    !     Non-blocking broadcast using dummy parameter to have output!
    !     processes wait for computations on first node to be finished.
    !     Neede for profiling purposes only.
    !
    !  2. Method :
    !
    !     Standard send and recieves using defined communicator. Send
    !     form first processor in communicator.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       IMOD    Int.   I   Model number.
    !       NMOD    Int.   I   Number of models.
    !       ID      Int.   I   ID number, used with NMOD for ITAG.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Sur.  W3SERVMD Subroutine tracing.
    !      W3SETO    Subr. W3ODATMD Point to data structure
    !      W3SETA    Subr. W3ADATMD Point to data structure
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      WMWAVE    Subr. WMWAVEMD Multi-grid wave model routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !  9. Switches :
    !
    !     !/S    Enable subroutine tracing.
    !     !/MPIT Enable test output
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
    !
    USE W3ADATMD, ONLY: W3SETA
    USE W3ADATMD, ONLY: MPI_COMM_WAVE
    USE W3ODATMD, ONLY: W3SETO                 ! W3_MPI
    USE W3ODATMD, ONLY: IAPROC, NAPROC, NTPROC ! W3_MPI
    USE WMMDATMD, ONLY: MDST, MDSE, MTAGB
    use WMMDATMD, only: mpi_comm_grd_not_null  ! W3_MPI
    USE W3SERVMD, ONLY: STRACE ! W3_S
    !
    IMPLICIT NONE
    !
#ifdef W3_MPI
    INCLUDE "mpif.h"
#endif
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    INTEGER, INTENT(IN)     :: IMOD, NMOD, ID
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
#ifdef W3_MPI
    INTEGER :: STATUS(MPI_STATUS_SIZE)
#endif
    INTEGER :: ITAG, IP, IERR_MPI ! W3_MPI
    INTEGER :: IENT = 0 ! W3_S
    REAL    :: DUMMY = 999. ! W3_MPI
    !/
    if (w3_s_flag) then
       CALL STRACE (IENT, 'WMWOUT')
    end if
    !
    ! -------------------------------------------------------------------- /
    ! 0.  Initializations
    !
#ifdef W3_MPI
    CALL W3SETO ( IMOD, MDSE, MDST )
    CALL W3SETA ( IMOD, MDSE, MDST )
    ITAG   = MTAGB + IMOD + ID*NMOD
    !
    !!/MPIT      WRITE (MDST,9000) IMOD, NMOD, ID, ITAG-MTAGB
    !!/MPIT      WRITE (MDST,9001) IAPROC, NAPROC, NTPROC
    IF ( IAPROC .LT. 1 ) THEN
       !!/MPIT          WRITE (MDST,9002)
       RETURN
    END IF
    !
    ! -------------------------------------------------------------------- /
    ! 1.  Processor to send data from
    !
    IF ( IAPROC .EQ. 1 ) THEN
       DO IP=NAPROC+1, NTPROC
          !!/MPIT   WRITE (MDST,9010) IAPROC, IP
          CALL MPI_SEND ( DUMMY, 1, MPI_INTEGER, IP-1, ITAG, MPI_COMM_WAVE, IERR_MPI )
       END DO
       !
       ! -------------------------------------------------------------------- /
       ! 2.  Processor to receive data at
       !
    ELSE IF ( IAPROC .GT. NAPROC ) THEN
       !!/MPIT  WRITE (MDST,9020) IAPROC, 1
       CALL MPI_RECV ( DUMMY, 1, MPI_INTEGER, 0, ITAG, MPI_COMM_WAVE, STATUS, IERR_MPI )
    ELSE
       !
       ! -------------------------------------------------------------------- /
       ! 3.  Processor with no action
       !
       !!/MPIT        ELSE
       !!/MPIT          WRITE (MDST,9030) IAPROC
    END IF
#endif
    !
    RETURN
    !
    ! Formats
    !
9000 FORMAT ( ' TEST WMWOUT : INPUTS :',4I4) ! W3_MPIT
9001 FORMAT ( ' TEST WMWOUT : IAPROC, NAPROC, NTPROC :',3I5)
9002 FORMAT ( ' TEST WMWOUT : NOT IN COMMUNICATOR')
    !
9010 FORMAT ( ' TEST WMWOUT : IAPROC =',I5,'    SENDING TO ',I5) ! W3_MPIT
9020 FORMAT ( ' TEST WMWOUT : IAPROC =',I5,'    RECEIVING FROM ',I5) ! W3_MPIT
9030 FORMAT ( ' TEST WMWOUT : IAPROC =',I5,'    NO ACTION') ! W3_MPIT
    !/
    !/ End of WMWOUT ----------------------------------------------------- /
    !/
  END SUBROUTINE WMWOUT
  !/
  !/ End of module WMWAVEMD -------------------------------------------- /
  !/
END MODULE WMWAVEMD
