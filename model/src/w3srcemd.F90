!> @file
!> @brief Source term integration routine.
!>
!> @author H. L. Tolman
!> @author F. Ardhuin
!> @date   22-Mar-2021
!>

#include "w3macros.h"
!/ ------------------------------------------------------------------- /

!>
!> @brief Source term integration routine.
!>
!> @author H. L. Tolman
!> @author F. Ardhuin
!> @date   22-Mar-2021
!>
!> @copyright Copyright 2009-2022 National Weather Service (NWS),
!>       National Oceanic and Atmospheric Administration.  All rights
!>       reserved.  WAVEWATCH III is a trademark of the NWS.
!>       No unauthorized use without permission.
!>
MODULE W3SRCEMD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |           H. L. Tolman            |
  !/                  |            F. Ardhuin             |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         22-Mar-2021 |
  !/                  +-----------------------------------+
  !/
  !/    For updates see subroutine.
  !/
  !  1. Purpose :
  !
  !     Source term integration routine.
  !
  !  2. Variables and types :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !      OFFSET    R.P.  Private  Offset in time integration scheme.
  !                               0.5 in original WAM, now 1.0
  !     ----------------------------------------------------------------
  !
  !  3. Subroutines and functions :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !      W3SRCE    Subr. Public   Calculate and integrate source terms.
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines and functions used :
  !
  !     See corresponding documentation of W3SRCE.
  !
  !  5. Remarks :
  !
  !  6. Switches :
  !
  !       See section 9 of W3SRCE.
  !
  !  7. Source code :
  !
  !/ ------------------------------------------------------------------- /
  !/
  REAL, PARAMETER, PRIVATE:: OFFSET = 1.
  !/
CONTAINS
  !/ ------------------------------------------------------------------- /

  !>
  !> @brief Calculate and integrate source terms for a single grid point.
  !>
  !> @verbatim
  !>     Physics  : see manual and corresponding subroutines.
  !>
  !>     Numerics :
  !>
  !>     Dynamic-implicit integration of the source terms based on
  !>     WW-II (Tolman 1992). The dynamic time step is calculated
  !>     given a maximum allowed change of spectral densities for
  !>     frequencies / wavenumbers below the usual cut-off.
  !>     The maximum change is given by the minimum of a parametric
  !>     and a relative change. The parametric change relates to a
  !>     PM type equilibrium range
  !>
  !>                                -1  (2pi)**4       1
  !>       dN(k)     =  Xp alpha  pi   ---------- ------------
  !>            max                       g**2     k**3 sigma
  !>
  !>                              1                                     .
  !>                 =  FACP ------------                              (1)
  !>                          k**3 sigma                                .
  !>
  !>     where
  !>           alpha = 0.62e-4                       (set in W3GRID)
  !>           Xp      fraction of PM shape          (read in W3GRID)
  !>           FACP    combined factor               (set in W3GRID)
  !>
  !>     The maximum relative change is given as
  !>
  !>                           /            +-                  -+ \    .
  !>       dN(k)     =  Xr max | N(k) , max | Nx , Xfilt N(k)    | |   (2)
  !>            max            \            +-               max-+ /    .
  !>
  !>     where
  !>           Xr      fraction of relative change   (read in W3GRID)
  !>           Xfilt   filter level                  (read in W3GRID)
  !>           Nx      Maximum parametric change (1)
  !>                   for largest wavenumber.
  !> @endverbatim
  !>
  !> @param[in]    srce_call
  !> @param[in]    IT
  !> @param[in]    ISEA
  !> @param[in]    JSEA
  !> @param[in]    IX         Discrete grid point counters.
  !> @param[in]    IY         Discrete grid point counters.
  !> @param[in]    IMOD       Model number.
  !> @param[in]    SPECOLD
  !> @param[inout] SPEC       Spectrum (action) in 1-D form.
  !> @param[out]   VSIO
  !> @param[out]   VDIO
  !> @param[out]   SHAVEIO
  !> @param[inout] ALPHA      Nondimensional 1-D spectrum corresponding
  !>                          to above full spectra (Phillip's const.).
  !> @param[inout] WN1        Discrete wavenumbers.
  !> @param[inout] CG1        Id. group velocities.
  !> @param[in]    CLATSL
  !> @param[in]    D_INP      Depth, compared to DMIN to get DEPTH.
  !> @param[in]    U10ABS     Wind speed at reference height.
  !> @param[in]    U10DIR     Id. wind direction.
  !> @param[inout] TAUA       Magnitude of total atmospheric stress.
  !> @param[inout] TAUADIR    Direction of atmospheric stress.
  !> @param[in]    AS         Air-sea temperature difference.
  !> @param[inout] USTAR      Friction velocity.
  !> @param[inout] USTDIR     Idem, direction.
  !> @param[in]    CX         Current velocity component.
  !> @param[in]    CY         Current velocity component.
  !> @param[in]    ICE        Sea ice concentration.
  !> @param[in]    ICEH       Sea ice thickness.
  !> @param[inout] ICEF       Sea ice floe diameter.
  !> @param[in]    ICEDMAX    Sea ice maximum floe diameter
  !> @param[in]    REFLEC     Reflection coefficients.
  !> @param[in]    REFLED     Reflection direction.
  !> @param[in]    DELX       Grid cell size in X direction.
  !> @param[in]    DELY       Grid cell size in Y direction.
  !> @param[in]    DELA       Grid cell area.
  !> @param[in]    TRNX       Grid transparency in X.
  !> @param[in]    TRNY       Grid transparency in Y.
  !> @param[in]    BERG       Iceberg damping coefficient.
  !> @param[inout] FPI        Peak-input frequency.
  !> @param[out]   DTDYN      Average dynamic time step.
  !> @param[out]   FCUT       Cut-off frequency for tail.
  !> @param[in]    DTG        Global time step.
  !> @param[inout] TAUWX
  !> @param[inout] TAUWY
  !> @param[inout] TAUOX
  !> @param[inout] TAUWIX
  !> @param[inout] TAUWIY
  !> @param[inout] TAUWNX
  !> @param[inout] TAUWNY
  !> @param[inout] PHIAW
  !> @param[inout] CHARN
  !> @param[inout] TWS
  !> @param[inout] PHIOC
  !> @param[inout] WHITECAP   Whitecap statistics.
  !> @param[in]    D50        Sand grain size.
  !> @param[in]    PSIC       Critical shields.
  !> @param[inout] BEDFORM    Bedform parameters.
  !> @param[inout] PHIBBL     Energy flux to BBL.
  !> @param[inout] TAUBBL     Momentum flux to BBL.
  !> @param[inout] TAUICE     Momentum flux to sea ice.
  !> @param[inout] PHICE      Energy flux to sea ice.
  !> @param[inout] TAUOCX     Total ocean momentum component.
  !> @param[inout] TAUOCY     Total ocean momentum component.
  !> @param[inout] WNMEAN     Mean wave number.
  !> @param[in]    DAIR       Air density.
  !> @param[in]    COEF
  !>
  !> @author H. L. Tolman
  !> @author F. Ardhuin
  !> @author A. Roland
  !> @author M. Dutour Sikiric
  !> @date   22-Mar-2021
  !>
  SUBROUTINE W3SRCE ( srce_call, IT, ISEA, JSEA, IX, IY, IMOD,          &
       SPECOLD, SPEC, VSIO, VDIO, SHAVEIO,         &
       ALPHA, WN1, CG1, CLATSL,                    &
       D_INP, U10ABS, U10DIR,                      &
#ifdef W3_FLX5
       TAUA, TAUADIR,                              &
#endif
       AS, USTAR, USTDIR,                          &
       CX, CY,  ICE, ICEH, ICEF, ICEDMAX,          &
       REFLEC, REFLED, DELX, DELY, DELA, TRNX,     &
       TRNY, BERG, FPI, DTDYN, FCUT, DTG, TAUWX,   &
       TAUWY, TAUOX, TAUOY, TAUWIX, TAUWIY, TAUWNX,&
       TAUWNY, PHIAW, CHARN, TWS, PHIOC, WHITECAP, &
       D50, PSIC, BEDFORM , PHIBBL, TAUBBL, TAUICE,&
       PHICE, TAUOCX, TAUOCY, WNMEAN, DAIR, COEF)

    !/ ------------------------------------------------------------------- /
    USE CONSTANTS, ONLY: DWAT, srce_imp_post, srce_imp_pre,         &
         srce_direct, GRAV, TPI, TPIINV, LPDLIB
#ifdef W3_T
    USE CONSTANTS, ONLY: RADE
#endif
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, TH, DMIN, DTMAX,       &
         DTMIN, FACTI1, FACTI2, FACSD, FACHFA, FACP, &
         XFC, XFLT, XREL, XFT, FXFM, FXPM, DDEN,     &
         FTE, FTF, FHMAX, ECOS, ESIN, IICEDISP,      &
         ICESCALES, IICESMOOTH
    USE W3GDATMD, ONLY: FSSOURCE, optionCall
    USE W3GDATMD, ONLY: B_JGS_NLEVEL, B_JGS_SOURCE_NONLINEAR
#ifdef W3_REF1
    USE W3GDATMD, ONLY: IOBP, IOBPD, IOBDP, GTYPE, UNGTYPE, REFPARS
#endif
    USE W3WDATMD, ONLY: TIME
    USE W3ODATMD, ONLY: NDSE, NDST, IAPROC
    USE W3IDATMD, ONLY: INFLAGS2, ICEP2
    USE W3DISPMD
#ifdef W3_NNT
    USE W3ODATMD, ONLY: IAPROC, SCREEN, FNMPRE
#endif
#ifdef W3_FLD1
    USE W3FLD1MD, ONLY: W3FLD1
    USE W3GDATMD, ONLY: AALPHA
#endif
#ifdef W3_FLD2
    USE W3FLD2MD, ONLY: W3FLD2
    USE W3GDATMD, ONLY: AALPHA
#endif
#ifdef W3_FLX1
    USE W3FLX1MD
#endif
#ifdef W3_FLX2
    USE W3FLX2MD
#endif
#ifdef W3_FLX3
    USE W3FLX3MD
#endif
#ifdef W3_FLX4
    USE W3FLX4MD
#endif
#ifdef W3_FLX5
    USE W3FLX5MD
#endif
#ifdef W3_LN1
    USE W3SLN1MD
#endif
#ifdef W3_ST0
    USE W3SRC0MD
#endif
#ifdef W3_ST1
    USE W3SRC1MD
#endif
#ifdef W3_ST2
    USE W3SRC2MD
    USE W3GDATMD, ONLY : ZWIND
#endif
#ifdef W3_ST3
    USE W3SRC3MD
    USE W3GDATMD, ONLY : ZZWND, FFXFM, FFXPM
#endif
#ifdef W3_ST4
    USE W3SRC4MD, ONLY : W3SPR4, W3SIN4, W3SDS4
    USE W3GDATMD, ONLY : ZZWND, FFXFM, FFXPM, FFXFA
#endif
#ifdef W3_ST6
    USE W3SRC6MD
    USE W3SWLDMD, ONLY : W3SWL6
    USE W3GDATMD, ONLY : SWL6S6
#endif
#ifdef W3_NL1
    USE W3SNL1MD
#endif
#ifdef W3_NL2
    USE W3SNL2MD
#endif
#ifdef W3_NL3
    USE W3SNL3MD
#endif
#ifdef W3_NL4
    USE W3SNL4MD
#endif
#ifdef W3_NL5
    USE W3SNL5MD
    USE W3TIMEMD, ONLY: TICK21
#endif
#ifdef W3_NLS
    USE W3SNLSMD
#endif
#ifdef W3_BT1
    USE W3SBT1MD
#endif
#ifdef W3_BT4
    USE W3SBT4MD
#endif
#ifdef W3_BT8
    USE W3SBT8MD
#endif
#ifdef W3_BT9
    USE W3SBT9MD
#endif
#ifdef W3_IC1
    USE W3SIC1MD
#endif
#ifdef W3_IC2
    USE W3SIC2MD
#endif
#ifdef W3_IC3
    USE W3SIC3MD
#endif
#ifdef W3_IC4
    USE W3SIC4MD
#endif
#ifdef W3_IC5
    USE W3SIC5MD
#endif
#ifdef W3_IS1
    USE W3SIS1MD
#endif
#ifdef W3_IS2
    USE W3SIS2MD
    USE W3GDATMD, ONLY : IS2PARS
#endif
#ifdef W3_DB1
    USE W3SDB1MD
#endif
#ifdef W3_TR1
    USE W3STR1MD
#endif
#ifdef W3_BS1
    USE W3SBS1MD
#endif
#ifdef W3_REF1
    USE W3REF1MD
#endif
#ifdef W3_IG1
    USE W3GDATMD, ONLY : IGPARS
#endif
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
#ifdef W3_NNT
    USE W3SERVMD, ONLY: EXTCDE
#endif
#ifdef W3_UOST
    USE W3UOSTMD, ONLY : UOST_SRCTRMCOMPUTE
#endif
#ifdef W3_PDLIB
    USE PDLIB_W3PROFSMD, ONLY : B_JAC, ASPAR_JAC, ASPAR_DIAG_SOURCES, ASPAR_DIAG_ALL
    USE yowNodepool,    ONLY: PDLIB_CCON, NPA, PDLIB_I_DIAG, PDLIB_JA, PDLIB_IA_P, PDLIB_SI
    USE W3GDATMD, ONLY: IOBP_LOC, IOBPD_LOC, IOBPA_LOC, IOBDP_LOC
    USE W3WDATMD, ONLY: VA
    USE W3PARALL, ONLY: ONESIXTH, ZERO, THR, IMEM, LSLOC
#endif
    !debug
    use w3odatmd     , only : naproc, iaproc
    !/
    IMPLICIT NONE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    INTEGER, INTENT(IN)     :: srce_call, IT, ISEA, JSEA, IX, IY, IMOD
    REAL, intent(in)        :: SPECOLD(NSPEC), CLATSL
    REAL, INTENT(OUT)       :: VSIO(NSPEC), VDIO(NSPEC)
    LOGICAL, INTENT(OUT)    :: SHAVEIO
    REAL, INTENT(IN)        :: D_INP, U10ABS,     &
         U10DIR, AS, CX, CY, DTG, D50,PSIC,   &
         ICE, ICEH
#ifdef W3_FLX5
    REAL, INTENT(IN)        :: TAUA, TAUADIR
#endif
    INTEGER, INTENT(IN)     :: REFLED(6)
    REAL, INTENT(IN)        :: REFLEC(4), DELX, DELY, DELA,         &
         TRNX, TRNY, BERG, ICEDMAX, DAIR
    REAL, INTENT(INOUT)     :: WN1(NK), CG1(NK), &
         SPEC(NSPEC), ALPHA(NK), USTAR,       &
         USTDIR, FPI, TAUOX, TAUOY,           &
         TAUWX, TAUWY, PHIAW, PHIOC, PHICE,   &
         CHARN, TWS, BEDFORM(3), PHIBBL,      &
         TAUBBL(2), TAUICE(2), WHITECAP(4),   &
         TAUWIX, TAUWIY, TAUWNX, TAUWNY,      &
         ICEF, TAUOCX, TAUOCY, WNMEAN
    REAL, INTENT(OUT)       :: DTDYN, FCUT
    REAL, INTENT(IN)        :: COEF
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
    INTEGER                 :: IK, ITH, IS, IS0, NSTEPS,  NKH, NKH1,&
         IKS1, IS1, NSPECH, IDT, IERR, NKD, ISP
    INTEGER                 :: IOBPIP, IOBPDIP, IOBDPIP
#ifdef W3_S
    INTEGER, SAVE           :: IENT = 0
#endif
#ifdef W3_NNT
    INTEGER, SAVE           :: NDSD = 89, NDSD2 = 88, J
#endif
#ifdef W3_NL5
    INTEGER                 :: QI5TSTART(2)
    REAL                    :: QR5KURT
    INTEGER, PARAMETER      :: NL5_SELECT = 1
    REAL, PARAMETER         :: NL5_OFFSET = 0.  ! explicit dyn.
#endif
    REAL                    :: DTTOT, FHIGH, DT, AFILT, DAMAX, AFAC,&
         HDT, ZWND, FP, DEPTH, TAUSCX, TAUSCY, FHIGI
    ! Scaling factor for SIN, SDS, SNL
    REAL                    :: ICESCALELN, ICESCALEIN, ICESCALENL, ICESCALEDS
    REAL                    :: EMEAN, FMEAN, AMAX, CD, Z0, SCAT,    &
         SMOOTH_ICEDISP
    REAL                    :: WN_R(NK), CG_ICE(NK),ALPHA_LIU(NK), ICECOEF2,&
         R(NK)
    DOUBLE PRECISION        :: ATT, ISO
#ifdef W3_ST1
    REAL                    :: FH1, FH2
#endif
#ifdef W3_ST2
    REAL                    :: FHTRAN, DFH, FACDIA, FACPAR
#endif
#ifdef W3_ST3
    REAL                    :: FMEANS, FH1, FH2
#endif
#ifdef W3_ST4
    REAL                    :: FMEANS, FH1, FH2, FAGE, DLWMEAN
#endif
    REAL                    :: QCERR  = 0.     !/XNL2 and !/NNT
#ifdef W3_SEED
    REAL                    :: UC, SLEV
#endif
#ifdef W3_MLIM
    REAL                    :: HM, EM
#endif
#ifdef W3_NNT
    REAL                    :: FACNN
#endif
#ifdef W3_T
    REAL                    :: DTRAW
#endif
    REAL                    :: EBAND, DIFF, EFINISH, HSTOT, PHINL,       &
         FMEAN1, FMEANWS, MWXINIT, MWYINIT,        &
         FACTOR, FACTOR2, DRAT, TAUWAX, TAUWAY,    &
         MWXFINISH, MWYFINISH, A1BAND, B1BAND,     &
         COSI(2)
    REAL                    :: SPECINIT(NSPEC), SPEC2(NSPEC), FRLOCAL, JAC2
    REAL                    :: DAM (NSPEC), DAM2(NSPEC), WN2 (NSPEC),            &
         VSLN(NSPEC),                         &
         VSIN(NSPEC), VDIN(NSPEC),            &
         VSNL(NSPEC), VDNL(NSPEC),            &
         VSDS(NSPEC), VDDS(NSPEC),            &
#ifdef W3_ST6
         VSWL(NSPEC), VDWL(NSPEC),            &
#endif
         VSBT(NSPEC), VDBT(NSPEC),            &
#ifdef W3_IC1
         VSIC(NSPEC), VDIC(NSPEC),            &
#endif
#ifdef W3_IC2
         VSIC(NSPEC), VDIC(NSPEC),            &
#endif
#ifdef W3_IC3
         VSIC(NSPEC), VDIC(NSPEC),            &
#endif
#ifdef W3_IC4
         VSIC(NSPEC), VDIC(NSPEC),            &
#endif
#ifdef W3_IC5
         VSIC(NSPEC), VDIC(NSPEC),            &
#endif
#ifdef W3_DB1
         VSDB(NSPEC), VDDB(NSPEC),            &
#endif
#ifdef W3_TR1
         VSTR(NSPEC), VDTR(NSPEC),            &
#endif
#ifdef W3_BS1
         VSBS(NSPEC), VDBS(NSPEC),            &
#endif
#ifdef W3_REF1
         VREF(NSPEC),                         &
#endif
#ifdef W3_IS1
         VSIR(NSPEC), VDIR(NSPEC),            &
#endif
#ifdef W3_IS2
         VSIR(NSPEC), VDIR(NSPEC),VDIR2(NSPEC), &
#endif
#ifdef W3_UOST
         VSUO(NSPEC), VDUO(NSPEC),            &
#endif
         VS(NSPEC), VD(NSPEC), EB(NK)
#ifdef W3_ST3
    LOGICAL                 :: LLWS(NSPEC)
#endif
#ifdef W3_ST4
    LOGICAL                 :: LLWS(NSPEC)
    REAL                    :: BRLAMBDA(NSPEC)
#endif
#ifdef W3_IS2
    DOUBLE PRECISION        :: SCATSPEC(NTH)
#endif
    REAL                    :: FOUT(NK,NTH), SOUT(NK,NTH), DOUT(NK,NTH)
    REAL, SAVE              :: TAUNUX, TAUNUY
#ifdef W3_OMPG
    !$omp threadprivate( TAUNUX, TAUNUY)
#endif
    LOGICAL, SAVE           :: FLTEST = .FALSE., FLAGNN = .TRUE.
#ifdef W3_OMPG
    !$omp threadprivate( FLTEST, FLAGNN )
#endif
    LOGICAL                 :: SHAVE
    LOGICAL                 :: LBREAK
    LOGICAL, SAVE           :: FIRST = .TRUE.
#ifdef W3_OMPG
    !$omp threadprivate( FIRST )
#endif
    LOGICAL                 :: PrintDeltaSmDA
    REAL                    :: eInc1, eInc2, eVS, eVD, JAC
    REAL                    :: DeltaSRC(NSPEC)
    REAL, PARAMETER         :: DTMINTOT = 0.01
    LOGICAL                 :: LNEWLIMITER = .FALSE.
#ifdef W3_PDLIB
    REAL                 :: PreVS, FAK, DVS, SIDT, FAKS, MAXDAC
#endif

#ifdef W3_NNT
    CHARACTER(LEN=17), SAVE :: FNAME = 'test_data_nnn.ww3'
#endif
    logical :: onde1, onde2
    onde1 = .false.
    onde2 = .false.
    if(naproc .eq. 10 .and. iaproc .eq. 2)onde2 = .true.
    if(naproc .eq.  5 .and. iaproc .eq. 1)onde1 = .true.
    !/
    !/ ------------------------------------------------------------------- /
    !/
#ifdef W3_S
    CALL STRACE (IENT, 'W3SRCE')
#endif
    !
#ifdef W3_T
    FLTEST = .TRUE.
#endif
    !
    VDIO   = 0.
    VSIO   = 0.
    DEPTH  = MAX ( DMIN , D_INP )

    IKS1 = 1
    ICESCALELN = MAX(0.,MIN(1.,1.-ICE*ICESCALES(1)))
    ICESCALEIN = MAX(0.,MIN(1.,1.-ICE*ICESCALES(2)))
    ICESCALENL = MAX(0.,MIN(1.,1.-ICE*ICESCALES(3)))
    ICESCALEDS = MAX(0.,MIN(1.,1.-ICE*ICESCALES(4)))
#ifdef W3_IG1
    !
    ! Does not integrate source terms for IG band if IGPARS(12) = 0.
    !
    IF (NINT(IGPARS(12)).EQ.0) IKS1 = NINT(IGPARS(5))
#endif
    IS1=(IKS1-1)*NTH+1
    !
#ifdef W3_LN0
    VSLN = 0.
#endif
#ifdef W3_LN1
    VSLN = 0.
#endif
#ifdef W3_SEED
    VSLN = 0.
#endif
#ifdef W3_ST0
    VSIN = 0.
    VDIN = 0.
#endif
#ifdef W3_ST3
    VSIN = 0.
    VDIN = 0.
#endif
#ifdef W3_ST4
    VSIN = 0.
    VDIN = 0.
#endif

#ifdef W3_NL0
    VSNL = 0.
    VDNL = 0.
#endif
#ifdef W3_NL1
    VSNL = 0.
    VDNL = 0.
#endif
#ifdef W3_TR1
    VSTR = 0.
    VDTR = 0.
#endif
#ifdef W3_ST0
    VSDS = 0.
    VDDS = 0.
#endif
#ifdef W3_ST4
    VSDS = 0.
    VDDS = 0.
#endif
    VSBT = 0.
    VDBT = 0.
#ifdef W3_DB1
    VSDB = 0.
    VDDB = 0.
#endif
#ifdef W3_IC1
    VSIC = 0.
    VDIC = 0.
#endif
#ifdef W3_IC2
    VSIC = 0.
    VDIC = 0.
#endif
#ifdef W3_IC3
    VSIC = 0.
    VDIC = 0.
#endif
#ifdef W3_IC4
    VSIC = 0.
    VDIC = 0.
#endif
#ifdef W3_UOST
    VSUO = 0.
    VDUO = 0.
#endif
#ifdef W3_IC5
    VSIC = 0.
    VDIC = 0.
#endif
    !
#ifdef W3_IS1
    VSIR = 0.
    VDIR = 0.
#endif
#ifdef W3_IS2
    VSIR = 0.
    VDIR = 0.
    VDIR2= 0.
#endif
    !
#ifdef W3_ST6
    VSWL = 0.
    VDWL = 0.
#endif
    !
#ifdef W3_ST0
    ZWND   = 10.
#endif
#ifdef W3_ST1
    ZWND   = 10.
#endif
#ifdef W3_ST2
    ZWND   = ZWIND
#endif
#ifdef W3_ST4
    ZWND   = ZZWND
#endif
#ifdef W3_ST6
    ZWND   = 10.
#endif
    !
    DRAT  = DAIR / DWAT
#ifdef W3_T
    WRITE (NDST,9000)
    WRITE (NDST,9001) DEPTH, U10ABS, U10DIR*RADE
#endif
    !
    ! 1.  Preparations --------------------------------------------------- *
    !
    ! 1.a Set maximum change and wavenumber arrays.
    !
    !XP     = 0.15
    !FACP   = XP / PI * 0.62E-3 * TPI**4 / GRAV**2
    !
    DO IK=1, NK
      DAM(1+(IK-1)*NTH) = FACP / ( SIG(IK) * WN1(IK)**3 )
      WN2(1+(IK-1)*NTH) = WN1(IK)
    END DO
    !
    DO IK=1, NK
      IS0    = (IK-1)*NTH
      DO ITH=2, NTH
        DAM(ITH+IS0) = DAM(1+IS0)
        WN2(ITH+IS0) = WN2(1+IS0)
      END DO
    END DO
    !
    ! 1.b Prepare dynamic time stepping
    !
    DTDYN  = 0.
    DTTOT  = 0.
    NSTEPS = 0
    PHIAW  = 0.
    CHARN  = 0.
    TWS    = 0.
    PHINL  = 0.
    PHIBBL = 0.
    TAUWIX = 0.
    TAUWIY = 0.
    TAUWNX = 0.
    TAUWNY = 0.
    TAUWAX = 0.
    TAUWAY = 0.
    TAUSCX = 0.
    TAUSCY = 0.
    TAUBBL = 0.
    TAUICE = 0.
    PHICE  = 0.
    TAUOCX = 0.
    TAUOCY = 0.
    WNMEAN = 0.

    !
    ! TIME is updated in W3WAVEMD prior to the call of W3SCRE, we should
    ! move 'TIME' one time step backward (QL)
#ifdef W3_ST4
    DLWMEAN= 0.
    BRLAMBDA(:)=0.
    WHITECAP(:)=0.
#endif
    !
    ! 1.c Set mean parameters
    !
    if(onde1 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2a ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
    if(onde2 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2a ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
#ifdef W3_ST4
    TAUWX=0.
    TAUWY=0.
    IF ( IT .eq. 0 ) THEN
      LLWS(:) = .TRUE.
      USTAR=0.
      USTDIR=0.
    ELSE
      ! tauwx,tauwy are intent(in) for spr4
      CALL W3SPR4 (SPEC, CG1, WN1, EMEAN, FMEAN, FMEAN1, WNMEAN, &
           AMAX, U10ABS, U10DIR,                           &
           USTAR, USTDIR,                                  &
           TAUWX, TAUWY, CD, Z0, CHARN, LLWS, FMEANWS, DLWMEAN)
      if(onde1 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2b_1 ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
      if(onde2 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2b_1 ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
      ! TAUWX, TAUWY, TAUWNX, TAUWNY are intent(out) for sin4
      CALL W3SIN4 ( SPEC, CG1, WN2, U10ABS, USTAR, DRAT, AS,       &
           U10DIR, Z0, CD, TAUWX, TAUWY, TAUWAX, TAUWAY,       &
           VSIN, VDIN, LLWS, IX, IY, BRLAMBDA )
      if(onde1 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2b_2 ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
      if(onde2 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2b_2 ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
    END IF
    if(onde1 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2b_3 ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
    if(onde2 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2b_3 ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
    CALL W3SPR4 (SPEC, CG1, WN1, EMEAN, FMEAN, FMEAN1, WNMEAN, &
         AMAX, U10ABS, U10DIR,                         &
         USTAR, USTDIR,                                &
         TAUWX, TAUWY, CD, Z0, CHARN, LLWS, FMEANWS, DLWMEAN)
    TWS = 1./FMEANWS
#endif
    if(onde1 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2b ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
    if(onde2 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2b ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
    !
    ! 1.c2 Stores the initial data
    !
    SPECINIT = SPEC
    !
    ! 1.d Stresses
    !
    !
    ! 1.e Prepare cut-off beyond which the tail is imposed with a power law
    !
#ifdef W3_ST4
    FAGE   = 0.
    FHIGH  = MAX( (FFXFM + FAGE ) * MAX(FMEAN1,FMEANWS), FFXPM / USTAR)
    FHIGI  = FFXFA * FMEAN1
#endif
    !
    ! ... Branch point dynamic integration - - - - - - - - - - - - - - - -
    !
    DO
      !
      NSTEPS = NSTEPS + 1
      !
      !
      ! 2.  Calculate source terms ----------------------------------------- *
      !
      ! 2.a Input.
      !
#ifdef W3_LN1
      CALL W3SLN1 (       WN1, FHIGH, USTAR, U10DIR , VSLN       )
#endif
      !
#ifdef W3_ST4
      CALL W3SIN4 ( SPEC, CG1, WN2, U10ABS, USTAR, DRAT, AS,       &
           U10DIR, Z0, CD, TAUWX, TAUWY, TAUWAX, TAUWAY,       &
           VSIN, VDIN, LLWS, IX, IY, BRLAMBDA )
#endif
      if(onde1 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2c ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
      if(onde2 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2c ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
      !
      ! 2.b Nonlinear interactions.
      !
#ifdef W3_NL1
      CALL W3SNL1 ( SPEC, CG1, WNMEAN*DEPTH,        VSNL, VDNL )
#endif
      !
      !
      ! 2.c Dissipation... except for ST4
      ! 2.c1   as in source term package
      !
#ifdef W3_ST4
      CALL W3SDS4 ( SPEC, WN1, CG1, USTAR, USTDIR, DEPTH, DAIR, VSDS,   &
           VDDS, IX, IY, BRLAMBDA, WHITECAP, DLWMEAN )
#endif

      !
#ifdef W3_PDLIB
      IF (.NOT. FSSOURCE .or. LSLOC) THEN
#endif
#ifdef W3_DB1
        CALL W3SDB1 ( IX, SPEC, DEPTH, EMEAN, FMEAN, WNMEAN, CG1,       &
             LBREAK, VSDB, VDDB )
#endif
#ifdef W3_PDLIB
      ENDIF
#endif
      !
      ! 2.c2   optional dissipation parameterisations
      !
      !
      ! 2.d Bottom interactions.
      !
#ifdef W3_BT1
      CALL W3SBT1 ( SPEC, CG1, WN1, DEPTH,            VSBT, VDBT )
#endif
      !
      ! 2.e Unresolved Obstacles Source Term
      !
#ifdef W3_UOST
      ! UNRESOLVED OBSTACLES
      CALL UOST_SRCTRMCOMPUTE(IX, IY, SPEC, CG1, DT,            &
           U10ABS, U10DIR, VSUO, VDUO)
#endif
      !
      ! 3.  Set frequency cut-off ------------------------------------------ *
      !
      NKH    = MIN ( NK , INT(FACTI2+FACTI1*LOG(MAX(1.E-7,FHIGH))) )
      NKH1   = MIN ( NK , NKH+1 )
      NSPECH = NKH1*NTH
      !
      ! 4.  Summation of source terms and diagonal term and time step ------ *
      !
      DT     = MIN ( DTG-DTTOT , DTMAX )
      AFILT  = MAX ( DAM(NSPEC) , XFLT*AMAX )
      !
      !     For input and dissipation calculate the fraction of the ice-free
      !     surface. In the presence of ice, the effective water surface
      !     is reduce to a fraction of the cell size free from ice, and so is
      !     input :
      !             SIN = (1-ICE)**ISCALEIN*SIN and SDS=(1-ICE)**ISCALEDS*SDS ------------------ *
      !     INFLAGS2(4) is true if ice concentration was ever read during
      !             this simulation
      IF ( INFLAGS2(4) ) THEN
        VSNL(1:NSPECH) = ICESCALENL * VSNL(1:NSPECH)
        VDNL(1:NSPECH) = ICESCALENL * VDNL(1:NSPECH)
        VSLN(1:NSPECH) = ICESCALELN * VSLN(1:NSPECH)
        VSIN(1:NSPECH) = ICESCALEIN * VSIN(1:NSPECH)
        VDIN(1:NSPECH) = ICESCALEIN * VDIN(1:NSPECH)
        VSDS(1:NSPECH) = ICESCALEDS * VSDS(1:NSPECH)
        VDDS(1:NSPECH) = ICESCALEDS * VDDS(1:NSPECH)
      END IF
      !
      VS = 0
      VD = 0
      DO IS=IS1, NSPECH
        VS(IS) = VSLN(IS) + VSIN(IS) + VSNL(IS)  &
             + VSDS(IS) + VSBT(IS)
#ifdef W3_ST6
        VS(IS) = VS(IS) + VSWL(IS)
#endif
#ifdef W3_TR1
        VS(IS) = VS(IS) + VSTR(IS)
#endif
#ifdef W3_BS1
        VS(IS) = VS(IS) + VSBS(IS)
#endif
#ifdef W3_UOST
        VS(IS) = VS(IS) + VSUO(IS)
#endif
        VD(IS) =  VDIN(IS) + VDNL(IS)  &
             + VDDS(IS) + VDBT(IS)
#ifdef W3_ST6
        VD(IS) = VD(IS) + VDWL(IS)
#endif
#ifdef W3_TR1
        VD(IS) = VD(IS) + VDTR(IS)
#endif
#ifdef W3_BS1
        VD(IS) = VD(IS) + VDBS(IS)
#endif
#ifdef W3_UOST
        VD(IS) = VD(IS) + VDUO(IS)
#endif
        DAMAX  = MIN ( DAM(IS) , MAX ( XREL*SPECINIT(IS) , AFILT ) )
        AFAC   = 1. / MAX( 1.E-10 , ABS(VS(IS)/DAMAX) )
          DT     = MIN ( DT , AFAC / ( MAX ( 1.E-10,                  &
               1. + OFFSET*AFAC*MIN(0.,VD(IS)) ) ) )
      END DO  ! end of loop on IS
      !
      DT     = MAX ( 0.5, DT ) ! The hardcoded min. dt is a problem for certain cases e.g. laborotary scale problems.
      !
      DTDYN  = DTDYN + DT
      IDT     = 1 + INT ( 0.99*(DTG-DTTOT)/DT ) ! number of iterations
      DT      = (DTG-DTTOT)/REAL(IDT)           ! actualy time step
      SHAVE   = DT.LT.DTMIN .AND. DT.LT.DTG-DTTOT   ! limiter check ...
      SHAVEIO = SHAVE
      DT      = MAX ( DT , MIN (DTMIN,DTG-DTTOT) ) ! override dt with input time step or last time step if it is bigger ... anyway the limiter is on!
      !
      IF (srce_call .eq. srce_imp_post) DT = DTG  ! for implicit part
        HDT    = OFFSET * DT
      DTTOT  = DTTOT + DT


#ifdef W3_PDLIB
      IF (srce_call .eq. srce_imp_pre) THEN
        IF (LSLOC) THEN
          IF (IMEM == 1) THEN
            SIDT  = PDLIB_SI(JSEA) * DTG
            DO IK = 1, NK
              JAC = CLATSL/CG1(IK)
              DO ITH = 1, NTH
                ISP = ITH + (IK-1)*NTH
                VD(ISP) = MIN(0., VD(ISP))
                IF (LNEWLIMITER) THEN
                  MAXDAC = MAX(DAM(ISP),DAM2(ISP))
                ELSE
                  MAXDAC = DAM(ISP)
                ENDIF
                FAKS   = DTG / MAX ( 1. , (1.-DTG*VD(ISP)))
                DVS    = VS(ISP) * FAKS
                DVS    = SIGN(MIN(MAXDAC,ABS(DVS)),DVS)
                PreVS  = DVS / FAKS
                eVS    = PreVS / CG1(IK) * CLATSL
                eVD    = MIN(0.,VD(ISP))
                B_JAC(ISP,JSEA)                   = B_JAC(ISP,JSEA) + SIDT * (eVS - eVD*SPEC(ISP)*JAC)
                ASPAR_JAC(ISP,PDLIB_I_DIAG(JSEA)) = ASPAR_JAC(ISP,PDLIB_I_DIAG(JSEA)) - SIDT * eVD
#ifdef W3_DB1
                eVS = VSDB(ISP) * JAC
                eVD = MIN(0.,VDDB(ISP))
                IF (eVS .gt. 0.) THEN
                  evS = 2*evS
                  evD = -evD
                ELSE
                  evS = -evS
                  evD = 2*evD
                ENDIF
#endif
                B_JAC(ISP,JSEA)                   = B_JAC(ISP,JSEA) + SIDT * eVS
                ASPAR_JAC(ISP,PDLIB_I_DIAG(JSEA)) = ASPAR_JAC(ISP,PDLIB_I_DIAG(JSEA)) - SIDT * eVD

                B_JAC(ISP,JSEA)                   = B_JAC(ISP,JSEA) + SIDT * eVS
                ASPAR_JAC(ISP,PDLIB_I_DIAG(JSEA)) = ASPAR_JAC(ISP,PDLIB_I_DIAG(JSEA)) - SIDT * eVD
              END DO
            END DO

          ELSEIF (IMEM == 2) THEN

            SIDT   = PDLIB_SI(JSEA) * DTG
            DO IK=1,NK
              JAC = CLATSL/CG1(IK)
              DO ITH=1,NTH
                ISP=ITH + (IK-1)*NTH
                VD(ISP) = MIN(0., VD(ISP))
                IF (LNEWLIMITER) THEN
                  MAXDAC    = MAX(DAM(ISP),DAM2(ISP))
                ELSE
                  MAXDAC    = DAM(ISP)
                ENDIF
                FAKS      = DTG / MAX ( 1. , (1.-DTG*VD(ISP)))
                DVS       = VS(ISP) * FAKS
                DVS       = SIGN(MIN(MAXDAC,ABS(DVS)),DVS)
                PreVS     = DVS / FAKS
                eVS = PreVS / CG1(IK) * CLATSL
                eVD = VD(ISP)
#ifdef W3_DB1
                eVS = eVS + DBLE(VSDB(ISP)) * JAC
                eVD = evD + MIN(0.,DBLE(VDDB(ISP)))
#endif
                B_JAC(ISP,JSEA)          = B_JAC(ISP,JSEA) + SIDT * (eVS - eVD*VA(ISP,JSEA))
                ASPAR_DIAG_ALL(ISP,JSEA) = ASPAR_DIAG_ALL(ISP,JSEA) - SIDT * eVD
              END DO
            END DO
          ENDIF
        ENDIF

        IF (.not. LSLOC) THEN
          IF (optionCall .eq. 1) THEN
            CALL SIGN_VSD_PATANKAR_WW3(SPEC,VS,VD)
          ELSE IF (optionCall .eq. 2) THEN
            CALL SIGN_VSD_SEMI_IMPLICIT_WW3(SPEC,VS,VD)
          ELSE IF (optionCall .eq. 3) THEN
            CALL SIGN_VSD_SEMI_IMPLICIT_WW3(SPEC,VS,VD)
          ENDIF
          VSIO = VS
          VDIO = VD
        ENDIF

        RETURN ! return everything is done for the implicit ...

      END IF ! srce_imp_pre
      if(onde1 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2d ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
      if(onde2 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2d ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY

#endif W3_PDLIB
      !
      !
      ! 5.  Increment spectrum --------------------------------------------- *
      !
      IF (srce_call .eq. srce_direct) THEN
        IF ( SHAVE ) THEN
          DO IS=IS1, NSPECH
            eInc1 = VS(IS) * DT / MAX ( 1. , (1.-HDT*VD(IS)))
            eInc2 = SIGN ( MIN (DAM(IS),ABS(eInc1)) , eInc1 )
            SPEC(IS) = MAX ( 0. , SPEC(IS)+eInc2 )
          END DO
        ELSE
          !
          DO IS=IS1, NSPECH
            eInc1 = VS(IS) * DT / MAX ( 1. , (1.-HDT*VD(IS)))
            SPEC(IS) = MAX ( 0. , SPEC(IS)+eInc1 )
          END DO
        END IF
        !
#ifdef W3_DB1
        DO IS=IS1, NSPECH
          eInc1 = VSDB(IS) * DT / MAX ( 1. , (1.-HDT*VDDB(IS)))
          SPEC(IS) = MAX ( 0. , SPEC(IS)+eInc1 )
        END DO
#endif

      END IF


      !
      ! 5.b  Computes
      !              atmos->wave flux PHIAW-------------------------------- *
      !              wave ->BBL  flux PHIBBL------------------------------- *
      !              wave ->ice  flux PHICE ------------------------------- *
      !
      WHITECAP(3)=0.
      HSTOT=0.
      DO IK=IKS1, NK
        FACTOR = DDEN(IK)/CG1(IK)                    !Jacobian to get energy in band
        FACTOR2= FACTOR*GRAV*WN1(IK)/SIG(IK)         ! coefficient to get momentum

        ! Wave direction is "direction to"
        ! therefore there is a PLUS sign for the stress
        DO ITH=1, NTH
          IS   = (IK-1)*NTH + ITH
          COSI(1)=ECOS(IS)
          COSI(2)=ESIN(IS)
          PHIAW = PHIAW + (VSIN(IS))* DT * FACTOR                    &
               / MAX ( 1. , (1.-HDT*VDIN(IS))) ! semi-implict integration scheme

          PHIBBL= PHIBBL- (VSBT(IS))* DT * FACTOR                    &
               / MAX ( 1. , (1.-HDT*VDBT(IS))) ! semi-implict integration scheme
          PHINL = PHINL + VSNL(IS)* DT * FACTOR                      &
               / MAX ( 1. , (1.-HDT*VDNL(IS))) ! semi-implict integration scheme
          IF (VSIN(IS).GT.0.) WHITECAP(3) = WHITECAP(3) + SPEC(IS)  * FACTOR
          HSTOT = HSTOT + SPEC(IS) * FACTOR
        END DO
      END DO
      WHITECAP(3)=4.*SQRT(WHITECAP(3))
      HSTOT=4.*SQRT(HSTOT)
      TAUWIX= TAUWIX+ TAUWX * DRAT *DT
      TAUWIY= TAUWIY+ TAUWY * DRAT *DT
      TAUWNX= TAUWNX+ TAUWAX * DRAT *DT
      TAUWNY= TAUWNY+ TAUWAY * DRAT *DT
      if(onde1 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2e ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
      if(onde2 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2e ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY

      ! MISSING: TAIL TO BE ADDED ?
      !
      !
      ! 6.  Add tail ------------------------------------------------------- *
      !   a Mean parameters
      !
      !
#ifdef W3_ST4
      CALL W3SPR4 (SPEC, CG1, WN1, EMEAN, FMEAN, FMEAN1, WNMEAN,&
           AMAX, U10ABS, U10DIR,                          &
           USTAR, USTDIR,                                 &
           TAUWX, TAUWY, CD, Z0, CHARN, LLWS, FMEANWS, DLWMEAN)

      ! Introduces a Long & Resio (JGR2007) type dependance on wave age
      FAGE   = FFXFA*TANH(0.3*U10ABS*FMEANWS*TPI/GRAV)
      FH1    = (FFXFM+FAGE) * FMEAN1

      FH2    = FFXPM / USTAR
      FHIGH  = MIN ( SIG(NK) , MAX ( FH1 , FH2 ) )
      NKH    = MAX ( 2 , MIN ( NKH1 ,                           &
           INT ( FACTI2 + FACTI1*LOG(MAX(1.E-7,FHIGH)) ) ) )
#endif
      !
      !
      ! 6.b Limiter for shallow water or Miche style criterion
      !     Last time step ONLY !
      !     uses true depth (D_INP) instead of limited depth
      !
#ifdef W3_MLIM
      IF ( DTTOT .GE. 0.9999*DTG ) THEN
        HM     = FHMAX *TANH(WNMEAN*MAX(0.,D_INP)) / MAX(1.E-4,WNMEAN )
        EM     = HM * HM / 16.
        IF ( EMEAN.GT.EM .AND. EMEAN.GT.1.E-30 ) THEN
          SPEC   = SPEC / EMEAN * EM
          EMEAN  = EM
        END IF
      END IF
#endif
      if(onde1 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2f ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
      if(onde2 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2f ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY

      !
      ! 6.c Seeding of spectrum
      !     alpha = 0.005 , 0.5 in eq., 0.25 for directional distribution
      !
#ifdef W3_SEED
      ! DO IK=MIN(NK,NKH), NK
      !   UC     = FACSD * GRAV / SIG(IK)
      !   SLEV   = MIN ( 1. , MAX ( 0. , U10ABS/UC-1. ) ) * 6.25E-4 / WN1(IK)**3 / SIG(IK)
      !   IF (INFLAGS2(4)) SLEV=SLEV*(1-ICE)
      !   DO ITH=1, NTH
      !     SPEC(ITH+(IK-1)*NTH) = MAX ( SPEC(ITH+(IK-1)*NTH) , SLEV * MAX ( 0. , COS(U10DIR-TH(ITH)) )**2 )
      !   END DO
      ! END DO
#endif
      !
      ! 6.d Add tail
      !
      DO IK=NKH+1, NK
        DO ITH=1, NTH
          SPEC(ITH+(IK-1)*NTH) = SPEC(ITH+(IK-2)*NTH) * FACHFA + 0.
        END DO
      END DO
      !
      ! 6.e  Update wave-supported stress----------------------------------- *
      !
#ifdef W3_ST4
      CALL W3SIN4 ( SPEC, CG1, WN2, U10ABS, USTAR, DRAT, AS,      &
           U10DIR, Z0, CD, TAUWX, TAUWY, TAUWAX, TAUWAY, &
           VSIN, VDIN, LLWS, IX, IY, BRLAMBDA )
#endif
      if(onde1 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2g ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
      if(onde2 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2g ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY

      !
      ! 7.  Check if integration complete ---------------------------------- *
      !
      IF (srce_call .eq. srce_imp_post) THEN
        EXIT
      ENDIF
      IF ( DTTOT .GE. 0.9999*DTG ) THEN
        !            IF (IX == DEBUG_NODE) WRITE(*,*) 'DTTOT, DTG', DTTOT, DTG
        EXIT
      ENDIF
    END DO ! INTEGRATIN LOOP
    !
    ! ... End point dynamic integration - - - - - - - - - - - - - - - - - -
    !
    ! 8.  Save integration data ------------------------------------------ *
    !
    DTDYN  = DTDYN / REAL(MAX(1,NSTEPS))
    FCUT   = FHIGH * TPIINV
    !
    GOTO 888
    !
    ! Error escape locations
    !
#ifdef W3_NNT
800 CONTINUE
    WRITE (NDSE,8000) FNAME, IERR
    CALL EXTCDE (1)
#endif
    !
#ifdef W3_NNT
801 CONTINUE
    WRITE (NDSE,8001) IERR
    CALL EXTCDE (2)
#endif
    !
888 CONTINUE
    if(onde1 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2h ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
    if(onde2 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2h ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY

    !
    ! 9.a  Computes PHIOC------------------------------------------ *
    !     The wave to ocean flux is the difference between initial energy
    !     and final energy, plus wind input plus the SNL flux to high freq.,
    !     minus the energy lost to the bottom boundary layer (BBL)
    !
    EFINISH  = 0.
    MWXFINISH  = 0.
    MWYFINISH  = 0.
    DO IK=1, NK
      EBAND = 0.
      A1BAND = 0.
      B1BAND = 0.
      DO ITH=1, NTH
        DIFF = SPECINIT(ITH+(IK-1)*NTH)-SPEC(ITH+(IK-1)*NTH)
        EBAND = EBAND + DIFF
        A1BAND = A1BAND + DIFF*ECOS(ITH)
        B1BAND = B1BAND + DIFF*ESIN(ITH)
      END DO
      EFINISH  = EFINISH  + EBAND * DDEN(IK) / CG1(IK)
      MWXFINISH  = MWXFINISH  + A1BAND * DDEN(IK) / CG1(IK)        &
           * WN1(IK)/SIG(IK)
      MWYFINISH  = MWYFINISH  + B1BAND * DDEN(IK) / CG1(IK)        &
           * WN1(IK)/SIG(IK)
    END DO
    !
    ! Transformation in momentum flux in m^2 / s^2
    !
    TAUOX=(GRAV*MWXFINISH+TAUWIX-TAUBBL(1))/DTG
    TAUOY=(GRAV*MWYFINISH+TAUWIY-TAUBBL(2))/DTG
    TAUWIX=TAUWIX/DTG
    TAUWIY=TAUWIY/DTG
    TAUWNX=TAUWNX/DTG
    TAUWNY=TAUWNY/DTG
    TAUBBL(:)=TAUBBL(:)/DTG
    TAUOCX=DAIR*COEF*COEF*USTAR*USTAR*COS(USTDIR) + DWAT*(TAUOX-TAUWIX)
    TAUOCY=DAIR*COEF*COEF*USTAR*USTAR*SIN(USTDIR) + DWAT*(TAUOY-TAUWIY)
    !
    ! Transformation in wave energy flux in W/m^2=kg / s^3
    !
    PHIOC =DWAT*GRAV*(EFINISH+PHIAW-PHIBBL)/DTG
    PHIAW =DWAT*GRAV*PHIAW /DTG
    PHINL =DWAT*GRAV*PHINL /DTG
    PHIBBL=DWAT*GRAV*PHIBBL/DTG
    !
    ! 10.1  Adds ice scattering and dissipation: implicit integration---------------- *
    !     INFLAGS2(4) is true if ice concentration was ever read during
    !             this simulation
    !
    if(onde1 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2i ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
    if(onde2 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2i ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY

    IF ( INFLAGS2(4).AND.ICE.GT.0 ) THEN

      IF (IICEDISP) THEN
        ICECOEF2 = 1E-6
        CALL LIU_FORWARD_DISPERSION (ICEH,ICECOEF2,DEPTH, &
             SIG,WN_R,CG_ICE,ALPHA_LIU)
        !
        IF (IICESMOOTH) THEN
        END IF
      ELSE
        WN_R=WN1
        CG_ICE=CG1
      END IF
      !
      R(:)=1 ! In case IC2 is defined but not IS2
      !
      SPEC2 = SPEC
      !
      TAUICE(:) = 0.
      PHICE = 0.
      DO IK=1,NK
        IS = 1+(IK-1)*NTH
        !
        ! First part of ice term integration: dissipation part
        !
        ATT=1.
        SPEC(1+(IK-1)*NTH:NTH+(IK-1)*NTH) = ATT*SPEC2(1+(IK-1)*NTH:NTH+(IK-1)*NTH)
        !
        ! Second part of ice term integration: scattering including re-distribution in directions
        !
        !
        ! 10.2  Fluxes of energy and momentum due to ice effects
        !
        FACTOR = DDEN(IK)/CG1(IK)                    !Jacobian to get energy in band
        FACTOR2= FACTOR*GRAV*WN1(IK)/SIG(IK)         ! coefficient to get momentum
        DO ITH = 1,NTH
          IS = ITH+(IK-1)*NTH
          PHICE = PHICE + (SPEC(IS)-SPEC2(IS)) * FACTOR
          COSI(1)=ECOS(IS)
          COSI(2)=ESIN(IS)
          TAUICE(:) = TAUICE(:) - (SPEC(IS)-SPEC2(IS))*FACTOR2*COSI(:)
        END DO
      END DO
      PHICE =-1.*DWAT*GRAV*PHICE /DTG
      TAUICE(:)=TAUICE(:)/DTG
    ELSE
    END IF
    if(onde1 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2j ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
    if(onde2 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2j ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
    !
    !
    ! - - - - - - - - - - - - - - - - - - - - - -
    ! 11. Sea state dependent stress routine calls
    ! - - - - - - - - - - - - - - - - - - - - - -
    !Note the Sea-state dependent stress calculations are primarily for high-wind
    !conditions (>10 m/s).  It is not recommended to use these at lower wind
    !in their current state.
    !

    ! FLD1/2 requires the calculation of FPI:
#ifdef W3_FLD2
    CALL CALC_FPI(SPEC, CG1, FPI, VSIN )

    IF (U10ABS.GT.10. .and. HSTOT.gt.0.5) then
      CALL W3FLD2 ( SPEC,min(FPI/TPI,2.0),COEF*U10ABS*COS(U10DIR),        &
           COEF*U10ABS*Sin(U10DIR), ZWND, DEPTH, 0.0, &
           DAIR, USTAR, USTDIR, Z0,TAUNUX,TAUNUY,CHARN)
    ELSE
      CHARN = AALPHA
    ENDIF
#endif
    if(onde1 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2k ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
    if(onde2 .and. ix .eq. 8442)print '(a,2i12,8g15.7)','DEBUGX2k ',time,TAUWIX, TAUWIY, TAUWNX, TAUWNY, TAUWX, TAUWY, TAUWAX, TAUWAY
    !
    ! 12. includes shoreline reflection --------------------------------------------- *
    !

    FIRST  = .FALSE.

    IF (IT.EQ.0) SPEC = SPECINIT

    SPEC = MAX(0., SPEC)
    !
    RETURN
    !
    ! Formats
    !
#ifdef W3_NNT
8000 FORMAT (/' *** ERROR W3SRCE : ERROR IN OPENING FILE ',A,' ***'/ &
         '                    IOSTAT = ',I10/)
8001 FORMAT (/' *** ERROR W3SRCE : ERROR IN WRITING TO FILE ***'/    &
         '                    IOSTAT = ',I10/)
#endif
    !
#ifdef W3_T
9000 FORMAT (' TEST W3SRCE : COUNTERS   : NO LONGER AVAILABLE')
9001 FORMAT (' TEST W3SRCE : DEPTH      :',F8.1/                     &
         '               WIND SPEED :',F8.1/                     &
         '               WIND DIR   :',F8.1)
#endif
#ifdef W3_ST1
9004 FORMAT (' TEST W3SRCE : FHIGH (3X) : ',3F8.4/                   &
         ' ------------- NEW DYNAMIC INTEGRATION LOOP',          &
         ' ------------- ')
#endif
#ifdef W3_ST2
9005 FORMAT (' TEST W3SRCE : FHIGH      : ',F8.4/                    &
         ' ------------- NEW DYNAMIC INTEGRATION LOOP',          &
         ' ------------- ')
#endif
#ifdef W3_ST3
9006 FORMAT (' TEST W3SRCE : FHIGH (3X) : ',3F8.4/                   &
         ' ------------- NEW DYNAMIC INTEGRATION LOOP',          &
         ' ------------- ')
#endif
#ifdef W3_ST4
9006 FORMAT (' TEST W3SRCE : FHIGH (3X) : ',3F8.4/                   &
         ' ------------- NEW DYNAMIC INTEGRATION LOOP',          &
         ' ------------- ')
#endif
    !
#ifdef W3_T
9020 FORMAT (' TEST W3SRCE : NSTEP : ',I4,'    DTTOT :',F6.1)
9021 FORMAT (' TEST W3SRCE : NKH (3X)   : ',2I3,I6)
#endif
    !
#ifdef W3_T
9040 FORMAT (' TEST W3SRCE : DTRAW, DT, SHAVE :',2F6.1,2X,L1)
#endif
    !
#ifdef W3_ST1
9060 FORMAT (' TEST W3SRCE : FHIGH (3X) : ',3F8.4/                   &
         '               NKH        : ',I3)
#endif
#ifdef W3_ST2
9061 FORMAT (' TEST W3SRCE : FHIGH (2X) : ',2F8.4/                   &
         '               NKH        : ',I3)
#endif
#ifdef W3_ST3
9062 FORMAT (' TEST W3SRCE : FHIGH (3X) : ',3F8.4/                   &
         '               NKH        : ',I3)
#endif
#ifdef W3_ST4
9062 FORMAT (' TEST W3SRCE : FHIGH (3X) : ',3F8.4/                   &
         '               NKH        : ',I3)
#endif
#ifdef W3_ST6
9063 FORMAT (' TEST W3SRCE : FHIGH      : ',F8.4/                    &
         '               NKH        : ',I3)
#endif
    !/
    !/ End of W3SRCE ----------------------------------------------------- /
    !/
  END SUBROUTINE W3SRCE
  !/ ------------------------------------------------------------------- /

  !>
  !> @brief Calculate equivalent peak frequency.
  !>
  !> @details Tolman and Chalikov (1996), equivalent peak frequency from source.
  !>
  !> @param[in]  A    Action density spectrum (1-D).
  !> @param[in]  CG   Group velocities for k-axis of spectrum.
  !> @param[out]  FPI  Input 'peak' frequency.
  !> @param[in] S    Source term (1-D version).
  !>
  !> @author Jessica Meixner
  !> @date   06-Jun-2018
  !>
  SUBROUTINE CALC_FPI( A, CG, FPI, S )
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |          Jessica Meixner          |
    !/                  |                                   |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         06-Jun-2018 |
    !/                  +-----------------------------------+
    !/
    !/    06-Jul-2016 : Origination                         ( version 5.12 )
    !/    06-Jul-2016 : Add SUBROUTINE SIGN_VSD_SEMI_IMPLICIT_WW3
    !/                  Add optional DEBUGSRC/PDLIB           ( version 6.04 )
    !/
    !  1. Purpose :
    !
    !     Calculate equivalent peak frequency
    !
    !  2. Method :
    !
    !     Tolman and Chalikov (1996), equivalent peak frequency from source

    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       A       R.A.  I   Action density spectrum (1-D).
    !       CG      R.A.  I   Group velocities for k-axis of spectrum.
    !       FPI     R.A.  O   Input 'peak' frequency.
    !       S       R.A.  I   Source term (1-D version).
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD Subroutine tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      W3SRCE Subr.
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
    !       !/S      Enable subroutine tracing.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
    USE CONSTANTS
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, XFR, DDEN, SIG,FTE, FTTR
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
    !
    IMPLICIT NONE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    REAL, INTENT(IN)        :: A(NSPEC), CG(NK), S(NSPEC)
    REAL, INTENT(OUT)       :: FPI
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
    INTEGER                 :: IS, IK
#ifdef W3_S
    INTEGER, SAVE           :: IENT = 0
#endif
    REAL                    ::  M0, M1, SIN1A(NK)
    !/
    !/ ------------------------------------------------------------------- /
    !/
#ifdef W3_S
    CALL STRACE (IENT, 'CALC_FPI')
#endif
    !
    !     Calculate FPI: equivalent peak frequncy from wind source term
    !     input
    !
    DO IK=1, NK
      SIN1A(IK) = 0.
      DO IS=(IK-1)*NTH+1, IK*NTH
        SIN1A(IK) = SIN1A(IK) + MAX ( 0. , S(IS) )
      END DO
    END DO
    !
    M0     = 0.
    M1     = 0.
    DO IK=1, NK
      SIN1A(IK) = SIN1A(IK) * DDEN(IK) / ( CG(IK) * SIG(IK)**3 )
      M0        = M0 + SIN1A(IK)
      M1        = M1 + SIN1A(IK)/SIG(IK)
    END DO
    !
    SIN1A(NK) = SIN1A(NK) / DDEN(NK)
    M0        = M0 + SIN1A(NK) * FTE
    M1        = M1 + SIN1A(NK) * FTTR
    IF ( M1 .LT. 1E-20 ) THEN
      FPI    = XFR * SIG(NK)
    ELSE
      FPI    = M0 / M1
    END IF

  END SUBROUTINE CALC_FPI
  !/ ------------------------------------------------------------------- /!

  !>
  !> @brief Put source term in matrix same as done always.
  !>
  !> @param[in]    SPEC
  !> @param[inout] VS
  !> @param[inout] VD
  !>
  !> @author Aron Roland
  !> @author Mathieu Dutour-Sikiric
  !> @date   01-Jun-2018
  !>
  SUBROUTINE SIGN_VSD_SEMI_IMPLICIT_WW3(SPEC, VS, VD)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |                                   |
    !/                  | Aron Roland (BGS IT&E GmbH)       |
    !/                  | Mathieu Dutour-Sikiric (IRB)      |
    !/                  |                                   |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :        01-June-2018 |
    !/                  +-----------------------------------+
    !/
    !/    01-June-2018 : Origination.                        ( version 6.04 )
    !/
    !  1. Purpose : Put source term in matrix same as done always
    !  2. Method :
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
    !      STRACE    Subr. W3SERVMD Subroutine tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !  7. Remarks
    !  8. Structure :
    !  9. Switches :
    !
    !     !/S  Enable subroutine tracing.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
    !
    USE W3GDATMD, only : NTH, NK, NSPEC
    IMPLICIT NONE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local PARAMETERs
    !/
#ifdef W3_S
    INTEGER, SAVE           :: IENT = 0
#endif
    !/
    !/ ------------------------------------------------------------------- /
    !/

    INTEGER             :: ISP, ITH, IK, IS
    REAL, INTENT(IN)    :: SPEC(NSPEC)
    REAL, INTENT(INOUT) :: VS(NSPEC), VD(NSPEC)
#ifdef W3_S
    CALL STRACE (IENT, 'SIGN_VSD_SEMI_IMPLICIT_WW3')
#endif
    DO IS=1,NSPEC
      VD(IS) = MIN(0., VD(IS))
    END DO
  END SUBROUTINE SIGN_VSD_SEMI_IMPLICIT_WW3
  !/ ------------------------------------------------------------------- /

  !>
  !> @brief Put source term in matrix Patankar style (experimental).
  !>
  !> @param[in]    SPEC
  !> @param[inout] VS
  !> @param[inout] VD
  !>
  !> @author Aron Roland
  !> @author Mathieu Dutour-Sikiric
  !> @date   01-Jun-2018
  !>
  SUBROUTINE SIGN_VSD_PATANKAR_WW3(SPEC, VS, VD)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |                                   |
    !/                  | Aron Roland (BGS IT&E GmbH)       |
    !/                  | Mathieu Dutour-Sikiric (IRB)      |
    !/                  |                                   |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :        01-June-2018 |
    !/                  +-----------------------------------+
    !/
    !/    01-June-2018 : Origination.                        ( version 6.04 )
    !/
    !  1. Purpose : Put source term in matrix Patankar style (experimental)
    !  2. Method :
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
    !      STRACE    Subr. W3SERVMD Subroutine tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !  7. Remarks
    !  8. Structure :
    !  9. Switches :
    !
    !     !/S  Enable subroutine tracing.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
    !

    USE W3GDATMD, only : NTH, NK, NSPEC
    IMPLICIT NONE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local PARAMETERs
    !/
#ifdef W3_S
    INTEGER, SAVE           :: IENT = 0
#endif
    !/
    !/ ------------------------------------------------------------------- /
    !/
    INTEGER             :: ISP, ITH, IK, IS
    REAL, INTENT(IN)    :: SPEC(NSPEC)
    REAL, INTENT(INOUT) :: VS(NSPEC), VD(NSPEC)
#ifdef W3_S
    CALL STRACE (IENT, 'SIGN_VSD_PATANKAR_WW3')
#endif
    DO IS=1,NSPEC
      VD(IS) = MIN(0., VD(IS))
      VS(IS) = MAX(0., VS(IS))
    END DO
  END SUBROUTINE SIGN_VSD_PATANKAR_WW3
  !/
  !/ End of module W3SRCEMD -------------------------------------------- /
  !/
END MODULE W3SRCEMD
