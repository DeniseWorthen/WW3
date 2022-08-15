!> @file
!> @brief Contains module WMFINLMD.
!> 
!> @author H. L. Tolman @date 04-Feb-2014
!> 

#include "w3macros.h"
!/ ------------------------------------------------------------------- /
!>
!> @brief Finalization of the multi-grid wave model.
!> 
!> @author H. L. Tolman @date 04-Feb-2014
!>
MODULE WMFINLMD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |           H. L. Tolman            |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         04-Feb-2014 |
  !/                  +-----------------------------------+
  !/
  !/    06-May-2005 : Origination.                        ( version 3.07 )
  !/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
  !/    03-Sep-2012 : Output of initilization time.       ( version 4.10 )
  !/    28-Jan-2014 : Add memory hwm to profiling.        ( version 5.00 )
  !/    04-Feb-2014 : Switched clock to DATE_AND_TIME     ( version 4.18 )
  !/                  (A. Chawla and Mark Szyszka) 
  !/
  !/    Copyright 2009-2014 National Weather Service (NWS),
  !/       National Oceanic and Atmospheric Administration.  All rights
  !/       reserved.  WAVEWATCH III is a trademark of the NWS. 
  !/       No unauthorized use without permission.
  !/
  !  1. Purpose :
  !
  !     Finalization of the multi-grid wave model.
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
  !      WMFINL    Subr. Public   Wave model initialization.
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
  !> @brief Initialize multi-grid version of WAVEWATCH III.
  !>
  !> @author H. L. Tolman  @date 28-Jan-2014
  !>        
  SUBROUTINE WMFINL
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           H. L. Tolman            |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         28-Jan-2014 |
    !/                  +-----------------------------------+
    !/
    !/    06-May-2005 : Origination.                        ( version 3.07 )
    !/    03-Sep-2012 : Output of initilization time.       ( version 4.10 )
    !/    28-Jan-2014 : Add memory hwm to profiling.        ( version 5.00 )
    !/
    !  1. Purpose :
    !
    !     Initialize multi-grid version of WAVEWATCH III.
    !
    !  2. Method :
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      PRTIME    Subr. W3SERVMD Profiling routine ( !/MPRF )
    !      MPI_BARRIER
    !                Subr.          Standard MPI routines.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      WW3_MULTI Prog.   N/A    Multi-grid model driver.
    !      ....                     Any coupled model.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    !       !/MPI   MPI routines.
    !
    !       !/O10   Enable output identifying start and end of routine
    !      
    !       !/S     Enable subroutine tracing.
    !       !/T     Enable test output
    !       !/MPRF  Profiling.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
    ! use w3getmem ; fake use statement for make_makefile.sh
    !
    USE W3TIMEMD, ONLY: TDIFF
    USE W3TIMEMD, ONLY: PRTIME ! W3_MPRF
    USE W3SERVMD, ONLY: STRACE ! W3_S
    USE WMMDATMD, ONLY: MDSS, MDSO, NMPSCR, NMPLOG, IMPROC
    USE WMMDATMD, ONLY: CLKDT1, CLKDT2, CLKDT3, CLKFIN
    USE WMMDATMD, ONLY: MDSP ! W3_MPRF
    USE WMMDATMD, ONLY: MPI_COMM_MWAVE ! W3_MPI
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
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
    INTEGER      :: IERR_MPI ! W3_MPI
    REAL         :: PRFT0, PRFTN ! W3_MPRF
#ifdef W3_MPRF
    REAL(KIND=8) :: get_memory
#endif
    INTEGER      :: IENT = 0 ! W3_S
    !/
    !/ ------------------------------------------------------------------- /
    ! 1.  Identification at start
    !
    if (w3_s_flag) then
       CALL STRACE (IENT, 'WMFINL')
    end if
    if (w3_mprf_flag) then
       CALL PRTIME ( PRFT0 )
    end if
    !
    if (w3_o10_flag) then
       IF ( MDSS.NE.MDSO .AND. NMPSCR.EQ.IMPROC ) WRITE (MDSS,900)
    end if
    !
    !/ ------------------------------------------------------------------- /
    ! 2.  Finalization
    !
    if (w3_mpi_flag) then
       CALL MPI_BARRIER ( MPI_COMM_MWAVE, IERR_MPI )
    end if
    !
    IF ( MDSS.NE.MDSO .AND. NMPSCR.EQ.IMPROC ) WRITE (MDSS,920) CLKFIN
    IF ( NMPLOG.EQ.IMPROC ) WRITE (MDSO,920) CLKFIN

    CALL DATE_AND_TIME ( VALUES=CLKDT3 )

    CLKFIN = TDIFF ( CLKDT1,CLKDT3 )
    IF ( MDSS.NE.MDSO .AND. NMPSCR.EQ.IMPROC ) WRITE (MDSS,921) CLKFIN
    IF ( NMPLOG.EQ.IMPROC ) WRITE (MDSO,921) CLKFIN
    !
    !/ ------------------------------------------------------------------- /
    ! 3.  Identification at end
    !
    if (w3_o10_flag) then
       IF ( MDSS.NE.MDSO .AND. NMPSCR.EQ.IMPROC ) WRITE (MDSS,999)
    end if
    !
#ifdef W3_MPRF
    CALL PRTIME ( PRFTN )
    WRITE (MDSP,990) PRFT0, PRFTN, get_memory()
#endif
    !
    RETURN
    !
    ! Formats
    !
900 FORMAT ( ' ========== STARTING MWW3 FINALIZATION (WMFINL) ===', &
         '============================' )
920 FORMAT (/'  Initialization time :',F10.2,' s')
921 FORMAT ( '  Elapsed time        :',F10.2,' s')

    !
990 FORMAT (1X,3F12.3,' WMFINL') ! W3_MPRF
    !
999 FORMAT (/' ========== END OF MWW3 INITIALIZATION (WMFINL) ===', &
         '============================'/)
    !/
    !/ End of WMFINL ----------------------------------------------------- /
    !/
  END SUBROUTINE WMFINL
  !/
  !/ End of module WMFINLMD -------------------------------------------- /
  !/
END MODULE WMFINLMD
