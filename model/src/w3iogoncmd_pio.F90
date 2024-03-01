!> @file w3iogoncmd
!!
!> @brief Write gridded model output as netCDF
!!
!> @author mvertens@ucar.edu, Denise.Worthen@noaa.gov
!> @date 01-05-2022

module w3iogoncmd_pio

  use w3parall          , only : init_get_isea
  use w3gdatmd          , only : xgrd, ygrd
  use w3gdatmd          , only : nk, nx, ny, mapsf, mapsta, nsea
  use w3odatmd          , only : noswll, undef
  use w3odatmd          , only : nds, iaproc, napout
  use w3adatmd          , only : mpi_comm_wave
  use wav_import_export , only : nseal_cpl
  use pio
  use netcdf

  implicit none

  private

  interface wav_initdecomp
     module procedure wav_initdecomp_2d
     module procedure wav_initdecomp_3d
  end interface wav_initdecomp

  public :: w3iogonc_pio
  public :: wav_initdecomp

  ! used/reused in module
  integer             :: isea, jsea, ix, iy, ierr
  integer             :: len_s, len_m, len_p, len_k
  character(len=1024) :: fname

  real, allocatable, target :: var3ds(:,:)
  real, allocatable, target :: var3dm(:,:)
  real, allocatable, target :: var3dp(:,:)
  real, allocatable, target :: var3dk(:,:)

  ! output variable for (nx,ny,nz) fields
  real, pointer :: var3d(:,:)

  type(iosystem_desc_t) :: wav_pio_subsystem
  type(file_desc_t)     :: pioid
  type(io_desc_t)       :: iodesc2d    !2d only
  type(io_desc_t)       :: iodesc3ds   !s-axis variables
  type(io_desc_t)       :: iodesc3dm   !m-axis variables
  type(io_desc_t)       :: iodesc3dp   !p-axis variables
  type(io_desc_t)       :: iodesc3dk   !k-axis variables
  type(var_desc_t)      :: vardesc

  integer :: pio_iotype

  !===============================================================================
contains
  !===============================================================================

  subroutine w3iogonc_pio ( )

    use w3odatmd   , only : fnmpre, naproc, iaproc
    use w3gdatmd   , only : filext, trigp, ntri, ungtype, gtype
    use w3servmd   , only : extcde
    use w3wdatmd   , only : w3setw, w3dimw, time, wlv, ice, icef, iceh, berg, ust, ustdir, asf, rhoair
    use w3gdatmd   , only : e3df, p2msf, us3df, usspf, w3setg
    use w3odatmd   , only : nogrp, ngrpp, idout, ndst, ndse,  noswll, w3seto
    use w3adatmd   , only : w3seta, w3dima
    use w3adatmd   , only : ainit, dw, ua, ud, as, cx, cy, wn, taua, tauadir
    use w3adatmd   , only : hs, wlm, t02, t0m1, t01, fp0, thm, ths, thp0, wbt, wnmean
    use w3adatmd   , only : dtdyn
    use w3adatmd   , only : fcut, aba, abd, uba, ubd, sxx, syy, sxy
    use w3adatmd   , only : phs, ptp, plp, pdir, psi, pws, pwst, pnr
    use w3adatmd   , only : pthp0, pqp, ppe, pgw, psw, ptm1, pt1, pt2
    use w3adatmd   , only : pep, usero, tauox, tauoy, tauwix, tauwiy
    use w3adatmd   , only : phiaw, phioc, tusx, tusy, prms, tpms
    use w3adatmd   , only : ussx, ussy, mssx, mssy, mssd, mscx, mscy
    use w3adatmd   , only : mscd, qp, tauwnx, tauwny, charn, tws, bhd
    use w3adatmd   , only : phibbl, taubbl, whitecap, bedforms, cge, ef
    use w3adatmd   , only : cflxymax, cflthmax, cflkmax, p2sms, us3d
    use w3adatmd   , only : th1m, sth1m, th2m, sth2m, hsig, phice, tauice
    use w3adatmd   , only : stmaxe, stmaxd, hmaxe, hcmaxe, hmaxd, hcmaxd, ussp, tauocx, tauocy
#ifdef W3_CESMCOUPLED
    use w3adatmd   , only : langmt
#endif
    use wav_grdout , only : varatts, outvars

    use w3timemd   , only : set_user_timestring
    use w3odatmd   , only : time_origin, calendar_name, elapsed_secs
    use w3odatmd   , only : use_user_histname, user_histfname
    !TODO: use unstr_mesh from wav_shr_mod; currently fails due to CI
    !use wav_shr_mod      , only : unstr_mesh

    !integer, intent(in)   :: timen(2)

    ! local variables
    integer             :: igrd
    integer    ,target  :: dimid3(3)
    integer    ,target  :: dimid4(4)
    integer    ,pointer :: dimid(:)
    character(len=12)   :: vname
    character(len=16)   :: user_timestring    !YYYY-MM-DD-SSSSS

    integer :: n, xtid, ytid, xeid, ztid, stid, mtid, ptid, ktid, timid
    logical :: s_axis = .false., m_axis = .false., p_axis = .false., k_axis = .false.

    ! pio
    integer :: nprocs
    integer :: istride
    integer :: basetask
    integer :: numiotasks
    integer :: rearranger
    integer :: my_task
    integer :: master_task
    integer :: status

    !-------------------------------------------------------------------------------

    igrd   = 1
    !TODO: are these needed?
    call w3seto ( igrd, ndse, ndst )
    call w3setg ( igrd, ndse, ndst )
    call w3seta ( igrd, ndse, ndst )  ! sets pointers into wadats in w3adatmd
    call w3setw ( igrd, ndse, ndst )  ! sets pointers into wdatas in w3wdatmd

    ! TODO: for now, hardwire  the io system
    pioid%fh = -1
    pio_iotype = PIO_IOTYPE_PNETCDF
    nprocs = naproc
    my_task = iaproc - 1
    master_task = 0
    istride = 4
    basetask = 1
    !numiotasks = max((nprocs-basetask)/istride,1)
    numiotasks = 2
    rearranger = PIO_REARR_BOX
    print '(a,6i8)','SETUP ',nprocs,iaproc,my_task,numiotasks, nseal_cpl, size(hs)

    call pio_init(my_task, MPI_COMM_WAVE, numiotasks, master_task, istride, rearranger, &
         wav_pio_subsystem, base=basetask)
    call pio_seterrorhandling(wav_pio_subsystem, PIO_RETURN_ERROR)

    ! -------------------------------------------------------------
    ! create the netcdf file
    ! -------------------------------------------------------------

    if (use_user_histname) then
      if (len_trim(user_histfname) == 0 ) then
        call extcde (60, msg="user history filename requested but not provided")
      end if
      call set_user_timestring(time,user_timestring)
      fname = trim(user_histfname)//trim(user_timestring)//'.nc'
    else
      write(fname,'(a,i8.8,a1,i6.6,a)')trim(fnmpre),time(1),'.',time(2),'.out_grd.'//trim(filext)//'.nc'
    end if

    len_s = noswll + 1                  ! 0:noswll
    len_m = p2msf(3)-p2msf(2) + 1       ! ?
    len_p = usspf(2)                    ! partitions
    len_k = e3df(3,1) - e3df(2,1) + 1   ! frequencies

    ! define the dimensions required for the requested gridded fields
    do n = 1,size(outvars)
      if (outvars(n)%validout) then
        if(trim(outvars(n)%dims) == 's')s_axis = .true.
        if(trim(outvars(n)%dims) == 'm')m_axis = .true.
        if(trim(outvars(n)%dims) == 'p')p_axis = .true.
        if(trim(outvars(n)%dims) == 'k')k_axis = .true.
      end if
    end do

    ! allocate arrays if needed
    if (s_axis) allocate(var3ds(1:nseal_cpl,len_s))
    if (m_axis) allocate(var3dm(1:nseal_cpl,len_m))
    if (p_axis) allocate(var3dp(1:nseal_cpl,len_p))
    if (k_axis) allocate(var3dk(1:nseal_cpl,len_k))

    ! create the netcdf file
    ierr = pio_createfile(wav_pio_subsystem, pioid, pio_iotype, trim(fname), pio_clobber)
    call handle_err(ierr, 'pio_create')

    ierr = pio_def_dim(pioid, 'nx', nx, xtid)
    ierr = pio_def_dim(pioid, 'ny', ny, ytid)
    ierr = pio_def_dim(pioid, 'time', PIO_UNLIMITED, timid)

    if (s_axis) ierr = pio_def_dim(pioid, 'noswll', len_s, stid)
    if (m_axis) ierr = pio_def_dim(pioid, 'nm'    , len_m, mtid)
    if (p_axis) ierr = pio_def_dim(pioid, 'np'    , len_p, ptid)
    if (k_axis) ierr = pio_def_dim(pioid, 'freq'  , len_k, ktid)
    if (gtype .eq. ungtype) then
      ierr = pio_def_dim(pioid, 'ne'  , ntri, xeid)
      ierr = pio_def_dim(pioid, 'nn'  ,    3, ztid)
    END IF

    ! define the time variable
    ierr = pio_def_var(pioid, 'time', PIO_DOUBLE, (/timid/), vardesc)
    call handle_err(ierr,'def_timevar')
    ierr = pio_put_att(pioid, vardesc, 'units', trim(time_origin))
    call handle_err(ierr,'def_time_units')
    ierr = pio_put_att(pioid, vardesc, 'calendar', trim(calendar_name))
    call handle_err(ierr,'def_time_calendar')

    ! define the spatial axis variables (lat,lon)
    ierr = pio_def_var(pioid, 'lon', PIO_DOUBLE, (/xtid,ytid/), vardesc)
    call handle_err(ierr,'def_lonvar')
    ierr = pio_put_att(pioid, vardesc, 'units', 'degrees_east')
    ierr = pio_def_var(pioid, 'lat', PIO_DOUBLE, (/xtid,ytid/), vardesc)
    call handle_err(ierr,'def_latvar')
    ierr = pio_put_att(pioid, vardesc, 'units', 'degrees_north')

    ! add mapsta
    ierr = pio_def_var(pioid, 'mapsta', PIO_INT, (/xtid, ytid, timid/), vardesc)
    call handle_err(ierr, 'def_mapsta')
    ierr = pio_put_att(pioid, vardesc, 'units', 'unitless')
    ierr = pio_put_att(pioid, vardesc, 'long_name', 'map status')

    if (gtype .eq. ungtype) then
      ierr = pio_def_var(pioid, 'nconn', PIO_INT, (/ztid,xeid/), vardesc)
      call handle_err(ierr,'def_nodeconnections')
      ierr = pio_put_att(pioid, vardesc, 'units', 'unitless')
      ierr = pio_put_att(pioid, vardesc, 'long_name', 'node connectivity')
    end if

    ! define the variables
    dimid3(1:2) = (/xtid, ytid/)
    dimid4(1:2) = (/xtid, ytid/)
    do n = 1,size(outvars)
      if (trim(outvars(n)%dims) == 's') then
        dimid4(3:4) = (/stid, timid/)
        dimid => dimid4
      else if (trim(outvars(n)%dims) == 'm') then
        dimid4(3:4) = (/mtid, timid/)
        dimid => dimid4
      else if (trim(outvars(n)%dims) == 'p') then
        dimid4(3:4) = (/ptid, timid/)
        dimid => dimid4
      else if (trim(outvars(n)%dims) == 'k') then
        dimid4(3:4) = (/ktid, timid/)
        dimid => dimid4
      else
        dimid3(3) = timid
        dimid => dimid3
      end if

      ierr = pio_def_var(pioid, trim(outvars(n)%var_name), PIO_REAL, dimid, vardesc)
      call handle_err(ierr, 'define variable '//trim((outvars(n)%var_name)))
      ierr = pio_put_att(pioid, vardesc, 'units'     , trim(outvars(n)%unit_name))
      ierr = pio_put_att(pioid, vardesc, 'long_name' , trim(outvars(n)%long_name))
      ierr = pio_put_att(pioid, vardesc, '_FillValue', undef)
    end do
    ! end variable definitions
    ierr = pio_enddef(pioid)
    call handle_err(ierr, 'end variable definition')

    call wav_initdecomp(iodesc=iodesc2d)
    if (s_axis)call wav_initdecomp(nz=len_s, iodesc=iodesc3ds)
    if (m_axis)call wav_initdecomp(nz=len_m, iodesc=iodesc3dm)
    if (p_axis)call wav_initdecomp(nz=len_p, iodesc=iodesc3dp)
    if (k_axis)call wav_initdecomp(nz=len_k, iodesc=iodesc3dk)

    ! write the time and spatial axis values (lat,lon,time)
    ierr = pio_inq_varid(pioid,  'lat', vardesc)
    call handle_err(ierr, 'inquire variable lat ')
    ierr = pio_put_var(pioid, vardesc, transpose(ygrd))
    call handle_err(ierr, 'put lat')

    ierr = pio_inq_varid(pioid,  'lon', vardesc)
    call handle_err(ierr, 'inquire variable lon ')
    ierr = pio_put_var(pioid, vardesc, transpose(xgrd))
    call handle_err(ierr, 'put lon')

    ierr = pio_inq_varid(pioid,  'time', vardesc)
    call handle_err(ierr, 'inquire variable time ')
    ierr = pio_put_var(pioid, vardesc, (/1/), real(elapsed_secs,8))
    call handle_err(ierr, 'put time')

    if (gtype .eq. ungtype) then
      ierr = pio_inq_varid(pioid,  'nconn', vardesc)
      call handle_err(ierr, 'inquire variable nconn ')
      ierr = pio_put_var(pioid, vardesc, trigp)
      call handle_err(ierr, 'put trigp')
    end if

    !maps
    ierr = pio_inq_varid(pioid,  'mapsta', vardesc)
    call handle_err(ierr, 'inquire variable mapsta ')
    ierr = pio_put_var(pioid, vardesc, transpose(mapsta))
    call handle_err(ierr, 'put mapsta')

    ! write the requested variables
    do n = 1,size(outvars)
       vname = trim(outvars(n)%var_name)
      if (trim(outvars(n)%dims) == 's') then
        var3d => var3ds
        ! Group 4
        if(vname .eq.      'PHS') call write_var3d(iodesc3ds, vname, phs      (1:nseal_cpl,0:noswll) )
        if(vname .eq.      'PTP') call write_var3d(iodesc3ds, vname, ptp      (1:nseal_cpl,0:noswll) )
        if(vname .eq.      'PLP') call write_var3d(iodesc3ds, vname, plp      (1:nseal_cpl,0:noswll) )
        if(vname .eq.     'PDIR') call write_var3d(iodesc3ds, vname, pdir     (1:nseal_cpl,0:noswll) )
        if(vname .eq.      'PSI') call write_var3d(iodesc3ds, vname, psi      (1:nseal_cpl,0:noswll) )
        if(vname .eq.      'PWS') call write_var3d(iodesc3ds, vname, pws      (1:nseal_cpl,0:noswll) )
        if(vname .eq.      'PDP') call write_var3d(iodesc3ds, vname, pthp0    (1:nseal_cpl,0:noswll) )
        if(vname .eq.      'PQP') call write_var3d(iodesc3ds, vname, pqp      (1:nseal_cpl,0:noswll) )
        if(vname .eq.      'PPE') call write_var3d(iodesc3ds, vname, ppe      (1:nseal_cpl,0:noswll) )
        if(vname .eq.      'PGW') call write_var3d(iodesc3ds, vname, pgw      (1:nseal_cpl,0:noswll) )
        if(vname .eq.      'PSW') call write_var3d(iodesc3ds, vname, psw      (1:nseal_cpl,0:noswll) )
        if(vname .eq.     'PTM1') call write_var3d(iodesc3ds, vname, ptm1     (1:nseal_cpl,0:noswll) )
        if(vname .eq.      'PT1') call write_var3d(iodesc3ds, vname, pt1      (1:nseal_cpl,0:noswll) )
        if(vname .eq.      'PT2') call write_var3d(iodesc3ds, vname, pt2      (1:nseal_cpl,0:noswll) )
        if(vname .eq.      'PEP') call write_var3d(iodesc3ds, vname, pep      (1:nseal_cpl,0:noswll) )

      else if (trim(outvars(n)%dims) == 'm') then                          ! m axis
        var3d => var3dm
        ! Group 6
        if (vname .eq.   'P2SMS') call write_var3d(iodesc3dm, vname, p2sms    (1:nseal_cpl,p2msf(2):p2msf(3)) )

      else if (trim(outvars(n)%dims) == 'p') then                          ! partition axis
        var3d => var3dp
        ! Group 6
        if (vname .eq.   'USSPX') call write_var3d(iodesc3dp, vname, ussp     (1:nseal_cpl,   1:usspf(2)) )
        if (vname .eq.   'USSPY') call write_var3d(iodesc3dp, vname, ussp     (1:nseal_cpl,nk+1:nk+usspf(2)) )

      else if (trim(outvars(n)%dims) == 'k') then                           ! freq axis
        var3d => var3dk
        ! Group 3
        if(vname .eq.       'EF') call write_var3d(iodesc3dk, vname, ef       (1:nseal_cpl,e3df(2,1):e3df(3,1)) )
        if(vname .eq.     'TH1M') call write_var3d(iodesc3dk, vname, ef       (1:nseal_cpl,e3df(2,2):e3df(3,2)) )
        if(vname .eq.    'STH1M') call write_var3d(iodesc3dk, vname, ef       (1:nseal_cpl,e3df(2,3):e3df(3,3)) )
        if(vname .eq.     'TH2M') call write_var3d(iodesc3dk, vname, ef       (1:nseal_cpl,e3df(2,4):e3df(3,4)) )
        if(vname .eq.    'STH2M') call write_var3d(iodesc3dk, vname, ef       (1:nseal_cpl,e3df(2,5):e3df(3,5)) )
        !TODO: wn has reversed indices (1:nk, 1:nseal_cpl)
        ! Group 6
        if (vname .eq.   'US3DX') call write_var3d(iodesc3dk, vname, us3d     (1:nseal_cpl,   us3df(2):us3df(3)) )
        if (vname .eq.   'US3DY') call write_var3d(iodesc3dk, vname, us3d     (1:nseal_cpl,nk+us3df(2):nk+us3df(3)) )

      else
        ! Group 1
        if (vname .eq.      'DW') call write_var2d(vname, dw       (1:nseal_cpl), init0='false')
        if (vname .eq.      'CX') call write_var2d(vname, cx       (1:nseal_cpl), init0='false')
        if (vname .eq.      'CY') call write_var2d(vname, cy       (1:nseal_cpl), init0='false')
        if (vname .eq.     'UAX') call write_var2d(vname, ua       (1:nseal_cpl), dir=cos(ud(1:nseal_cpl)), init0='false')
        if (vname .eq.     'UAY') call write_var2d(vname, ua       (1:nseal_cpl), dir=sin(ud(1:nseal_cpl)), init0='false')
        if (vname .eq.      'AS') call write_var2d(vname, as       (1:nseal_cpl), init0='false')
        if (vname .eq.     'WLV') call write_var2d(vname, wlv      (1:nseal_cpl), init0='false')
        if (vname .eq.     'ICE') call write_var2d(vname, ice      (1:nseal_cpl), init0='false')
        if (vname .eq.    'BERG') call write_var2d(vname, berg     (1:nseal_cpl), init0='false')
        if (vname .eq.    'TAUX') call write_var2d(vname, taua     (1:nseal_cpl), dir=cos(tauadir(1:nseal_cpl)), init0='false')
        if (vname .eq.    'TAUY') call write_var2d(vname, taua     (1:nseal_cpl), dir=sin(tauadir(1:nseal_cpl)), init0='false')
        if (vname .eq.  'RHOAIR') call write_var2d(vname, rhoair   (1:nseal_cpl), init0='false')
        if (vname .eq.    'ICEH') call write_var2d(vname, iceh     (1:nseal_cpl), init0='false')
        if (vname .eq.    'ICEF') call write_var2d(vname, icef     (1:nseal_cpl), init0='false')

        ! Group 2
        if (vname .eq.      'HS') call write_var2d(vname, hs       (1:nseal_cpl) )
        if (vname .eq.     'WLM') call write_var2d(vname, wlm      (1:nseal_cpl) )
        if (vname .eq.     'T02') call write_var2d(vname, t02      (1:nseal_cpl) )
        if (vname .eq.    'T0M1') call write_var2d(vname, t0m1     (1:nseal_cpl) )
        if (vname .eq.     'T01') call write_var2d(vname, t01      (1:nseal_cpl) )
        if (vname .eq.     'FP0') call write_var2d(vname, fp0      (1:nseal_cpl) )
        if (vname .eq.     'THM') call write_var2d(vname, thm      (1:nseal_cpl) )
        if (vname .eq.     'THS') call write_var2d(vname, ths      (1:nseal_cpl) )
        if (vname .eq.    'THP0') call write_var2d(vname, thp0     (1:nseal_cpl) )
        if (vname .eq.    'HSIG') call write_var2d(vname, hsig     (1:nseal_cpl) )
        if (vname .eq.  'STMAXE') call write_var2d(vname, stmaxe   (1:nseal_cpl) )
        if (vname .eq.  'STMAXD') call write_var2d(vname, stmaxd   (1:nseal_cpl) )
        if (vname .eq.   'HMAXE') call write_var2d(vname, hmaxe    (1:nseal_cpl) )
        if (vname .eq.  'HCMAXE') call write_var2d(vname, hcmaxe   (1:nseal_cpl) )
        if (vname .eq.   'HMAXD') call write_var2d(vname, hmaxd    (1:nseal_cpl) )
        if (vname .eq.  'HCMAXD') call write_var2d(vname, hcmaxd   (1:nseal_cpl) )
        if (vname .eq.     'WBT') call write_var2d(vname, wbt      (1:nseal_cpl) )
        if (vname .eq.  'WNMEAN') call write_var2d(vname, wnmean   (1:nseal_cpl), init0='false')

        ! Group 4
        if(vname .eq.     'PWST') call write_var2d(vname, pwst     (1:nseal_cpl) )
        if(vname .eq.      'PNR') call write_var2d(vname, pnr      (1:nseal_cpl) )

        ! Group 5
        if (vname .eq.    'USTX') call write_var2d(vname, ust      (1:nseal_cpl)*asf(1:nseal_cpl), dir=cos(ustdir(1:nseal_cpl)), usemask='true')
        if (vname .eq.    'USTY') call write_var2d(vname, ust      (1:nseal_cpl)*asf(1:nseal_cpl), dir=sin(ustdir(1:nseal_cpl)), usemask='true')
        if (vname .eq.   'CHARN') call write_var2d(vname, charn    (1:nseal_cpl) )
        if (vname .eq.     'CGE') call write_var2d(vname, cge      (1:nseal_cpl) )
        if (vname .eq.   'PHIAW') call write_var2d(vname, phiaw    (1:nseal_cpl),   init2='true')
        if (vname .eq.  'TAUWIX') call write_var2d(vname, tauwix   (1:nseal_cpl),   init2='true')
        if (vname .eq.  'TAUWIY') call write_var2d(vname, tauwiy   (1:nseal_cpl),   init2='true')
        if (vname .eq.  'TAUWNX') call write_var2d(vname, tauwnx   (1:nseal_cpl),   init2='true')
        if (vname .eq.  'TAUWNY') call write_var2d(vname, tauwny   (1:nseal_cpl),   init2='true')
        if (vname .eq.     'WCC') call write_var2d(vname, whitecap (1:nseal_cpl,1), init2='true')
        if (vname .eq.     'WCF') call write_var2d(vname, whitecap (1:nseal_cpl,2), init2='true')
        if (vname .eq.     'WCH') call write_var2d(vname, whitecap (1:nseal_cpl,3), init2='true')
        if (vname .eq.     'WCM') call write_var2d(vname, whitecap (1:nseal_cpl,4), init2='true')
        if (vname .eq.     'TWS') call write_var2d(vname, tws      (1:nseal_cpl) )

        ! Group 6
        if (vname .eq.     'SXX') call write_var2d(vname, sxx      (1:nseal_cpl) )
        if (vname .eq.     'SYY') call write_var2d(vname, syy      (1:nseal_cpl) )
        if (vname .eq.     'SXY') call write_var2d(vname, sxy      (1:nseal_cpl) )
        if (vname .eq.   'TAUOX') call write_var2d(vname, tauox    (1:nseal_cpl),   init2='true')
        if (vname .eq.   'TAUOY') call write_var2d(vname, tauoy    (1:nseal_cpl),   init2='true')
        if (vname .eq.     'BHD') call write_var2d(vname, bhd      (1:nseal_cpl) )
        if (vname .eq.   'PHIOC') call write_var2d(vname, phioc    (1:nseal_cpl),   init2='true')
        if (vname .eq.    'TUSX') call write_var2d(vname, tusx     (1:nseal_cpl) )
        if (vname .eq.    'TUSY') call write_var2d(vname, tusy     (1:nseal_cpl) )
        if (vname .eq.    'USSX') call write_var2d(vname, ussx     (1:nseal_cpl) )
        if (vname .eq.    'USSY') call write_var2d(vname, ussy     (1:nseal_cpl) )
        if (vname .eq.    'PRMS') call write_var2d(vname, prms     (1:nseal_cpl) )
        if (vname .eq.    'TPMS') call write_var2d(vname, tpms     (1:nseal_cpl) )
        if (vname .eq. 'TAUICEX') call write_var2d(vname, tauice   (1:nseal_cpl,1) )
        if (vname .eq. 'TAUICEY') call write_var2d(vname, tauice   (1:nseal_cpl,2) )
        if (vname .eq.   'PHICE') call write_var2d(vname, phice    (1:nseal_cpl) )
        if (vname .eq.  'TAUOCX') call write_var2d(vname, tauocx   (1:nseal_cpl) )
        if (vname .eq.  'TAUOCY') call write_var2d(vname, tauocy   (1:nseal_cpl) )
#ifdef W3_CESMCOUPLED
        if (vname .eq.  'LANGMT') call write_var2d(vname, langmt   (1:nseal_cpl) )
#endif
        ! Group 7
        if (vname .eq.    'ABAX') call write_var2d(vname, aba      (1:nseal_cpl), dir=cos(abd(1:nseal_cpl)) )
        if (vname .eq.    'ABAY') call write_var2d(vname, aba      (1:nseal_cpl), dir=sin(abd(1:nseal_cpl)) )
        if (vname .eq.    'UBAX') call write_var2d(vname, uba      (1:nseal_cpl), dir=cos(ubd(1:nseal_cpl)) )
        if (vname .eq.    'UBAY') call write_var2d(vname, uba      (1:nseal_cpl), dir=sin(ubd(1:nseal_cpl)) )
        if (vname .eq.     'BED') call write_var2d(vname, bedforms (1:nseal_cpl,1), init2='true')
        if (vname .eq. 'RIPPLEX') call write_var2d(vname, bedforms (1:nseal_cpl,2), init2='true')
        if (vname .eq. 'RIPPLEY') call write_var2d(vname, bedforms (1:nseal_cpl,3), init2='true')
        if (vname .eq.  'PHIBBL') call write_var2d(vname, phibbl   (1:nseal_cpl),   init2='true')
        if (vname .eq. 'TAUBBLX') call write_var2d(vname, taubbl   (1:nseal_cpl,1), init2='true')
        if (vname .eq. 'TAUBBLY') call write_var2d(vname, taubbl   (1:nseal_cpl,2), init2='true')

        ! Group 8
        if (vname .eq.    'MSSX') call write_var2d(vname, mssx     (1:nseal_cpl) )
        if (vname .eq.    'MSSY') call write_var2d(vname, mssy     (1:nseal_cpl) )
        if (vname .eq.    'MSCX') call write_var2d(vname, mscx     (1:nseal_cpl) )
        if (vname .eq.    'MSCY') call write_var2d(vname, mscy     (1:nseal_cpl) )
        !TODO: remaining variables have inconsistency between shel_inp listing and iogo code

        ! Group 9
        if (vname .eq.   'DTDYN') call write_var2d(vname, dtdyn    (1:nseal_cpl) )
        if (vname .eq.    'FCUT') call write_var2d(vname, fcut     (1:nseal_cpl) )
        if (vname .eq.'CFLXYMAX') call write_var2d(vname, cflxymax (1:nseal_cpl) )
        if (vname .eq.'CFLTHMAX') call write_var2d(vname, cflthmax (1:nseal_cpl) )
        if (vname .eq. 'CFLKMAX') call write_var2d(vname, cflkmax  (1:nseal_cpl) )

        ! Group 10
      end if
    end do

    if (s_axis) deallocate(var3ds)
    if (m_axis) deallocate(var3dm)
    if (p_axis) deallocate(var3dp)
    if (k_axis) deallocate(var3dk)

    call pio_freedecomp(pioid,iodesc2d)
    if (s_axis) call pio_freedecomp(pioid, iodesc3ds)
    if (m_axis) call pio_freedecomp(pioid, iodesc3dm)
    if (p_axis) call pio_freedecomp(pioid, iodesc3dp)
    if (k_axis) call pio_freedecomp(pioid, iodesc3dk)

    call pio_closefile(pioid)

  end subroutine w3iogonc_pio

  !===============================================================================
  subroutine write_var2d(vname, var, dir, usemask, init0, init2)
    ! write (nseal) array as (nx,ny)
    ! if dir is present, write x or y component of (nsea) array as (nx,ny)
    ! if mask is present and true, use mapsta=1 to mask values
    ! if init0 is present and false, do not initialize values
    ! for mapsta<0. this prevents group 1 variables being set undef over
    ! ice. if init2 is present and true, apply a second initialization to
    ! a subset of variables for where mapsta==2

    character(len=*),           intent(in) :: vname
    real            ,           intent(in) :: var(:)
    real, optional  ,           intent(in) :: dir(:)
    character(len=*), optional, intent(in) :: usemask
    character(len=*), optional, intent(in) :: init0
    character(len=*), optional, intent(in) :: init2

    ! local variables
    real, dimension(nseal_cpl) :: varout
    logical                    :: lmask, linit0, linit2
    real                       :: varloc

    lmask = .false.
    if (present(usemask)) then
      lmask = (trim(usemask) == "true")
    end if
    linit0 = .true.
    if (present(init0)) then
      linit0 = (trim(init0) == "true")
    end if
    linit2 = .false.
    if (present(init2)) then
      linit2 = (trim(init2) == "true")
    end if

    ! DEBUG
    !write(*,'(a)')' writing variable ' //trim(vname)//' to history file '//trim(fname)

    varout = undef
    do jsea = 1,nseal_cpl
      call init_get_isea(isea, jsea)
      ! initialization
      !print '(a,5i6,f8.2)','WWvar2d ',jsea,isea,ix,iy,n,var(jsea)

      varloc = var(jsea)
      if (linit0) then
        if (mapsta(mapsf(isea,2),mapsf(isea,1)) < 0) varloc = undef
      end if
      if (linit2) then
        if (mapsta(mapsf(isea,2),mapsf(isea,1)) == 2) varloc = undef
      end if

      if (present(dir)) then
        if (varloc .ne. undef) then
          if (lmask) then
            if (mapsta(mapsf(isea,2),mapsf(isea,1)) == 1) then
              varout(jsea) = varloc*dir(jsea)
            end if
          else
            varout(jsea) = varloc*dir(jsea)
          end if
        end if
      else
        varout(jsea) = varloc
      end if
    end do

    ierr = pio_inq_varid(pioid,  trim(vname), vardesc)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, vardesc, int(1,kind=Pio_Offset_Kind))
    call pio_write_darray(pioid, vardesc, iodesc2d, varout, ierr)
    call handle_err(ierr, 'put variable '//trim(vname))

  end subroutine write_var2d

  !===============================================================================
  subroutine write_var3d(iodesc, vname, var, init2)
    ! write (nseal,:) array as (nx,ny,:)
    ! if init2 is present and true, apply a second initialization to
    ! a subset of variables for where mapsta==2

    type(io_desc_t),         intent(inout) :: iodesc
    character(len=*),           intent(in) :: vname
    real            ,           intent(in) :: var(:,:)
    character(len=*), optional, intent(in) :: init2

    ! local variables
    real, allocatable, dimension(:) :: varloc
    logical                         :: linit2
    integer                         :: lb, ub

    linit2 = .false.
    if (present(init2)) then
      linit2 = (trim(init2) == "true")
    end if

    lb = lbound(var,2)
    ub = ubound(var,2)
    allocate(varloc(lb:ub))

    ! DEBUG
    ! write(nds(1),'(a,2i6)')' writing variable ' //trim(vname)//' to history file ' &
    !    //trim(fname)//' with bounds ',lb,ub
    !print *,'XYXY ',trim(vname),minval(var3d,1),maxval(var3d,1),minval(var3d,2),maxval(var3d,2)

    var3d = undef
    do jsea = 1,nseal_cpl
      call init_get_isea(isea, jsea)
      ! initialization
      varloc(:) = var(jsea,:)
      print '(a,2i6,2f8.2)','WWvar2d ',jsea,isea,minval(varloc),maxval(varloc)
      if (mapsta(mapsf(isea,2),mapsf(isea,1)) < 0) varloc(:) = undef
      if (linit2) then
        if (mapsta(mapsf(isea,2),mapsf(isea,1)) == 2) varloc(:) = undef
      end if
      var3d(jsea,:) = varloc(:)
      !print *,trim(vname),jsea,lb,ub,varloc
    end do

    ierr = pio_inq_varid(pioid,  trim(vname), vardesc)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, vardesc, int(1,kind=PIO_OFFSET_KIND))
    call pio_write_darray(pioid, vardesc, iodesc, var3d, ierr)

    deallocate(varloc)
  end subroutine write_var3d

  !===============================================================================
  subroutine handle_err(ierr,string)
    use w3odatmd  , only : ndse
    use w3servmd  , only : extcde

    ! input/output variables
    integer         , intent(in) :: ierr
    character(len=*), intent(in) :: string

    integer :: strerror_status
    character(len=pio_max_name) :: err_msg

    if (ierr /= PIO_NOERR) then
      strerror_status = pio_strerror(ierr, err_msg)
      write(ndse,*) "*** WAVEWATCH III netcdf error: ",trim(string),':',trim(err_msg)
      call extcde ( 49 )
    end if
  end subroutine handle_err

  !===============================================================================
  subroutine wav_initdecomp_2d(iodesc)

    type(io_desc_t), intent(out) :: iodesc

    integer          :: n
    integer, pointer :: dof2d(:)

    allocate(dof2d(nseal_cpl))
    dof2d = 0
    n = 0
    do jsea = 1,nseal_cpl
      call init_get_isea(isea, jsea)
      ix = mapsf(isea,1)                 ! global ix
      iy = mapsf(isea,2)                 ! global iy
      n = n+1
      dof2d(n) = (iy-1)*nx + ix          ! local index : global index
      !print '(a,4i6,2f10.2)','PIO ',jsea,ix,iy,(iy-1)*nx + ix,xgrd(iy,ix),ygrd(iy,ix)
    end do
    call pio_initdecomp(wav_pio_subsystem, PIO_REAL, (/nx,ny/), dof2d, iodesc)

    deallocate(dof2d)
  end subroutine wav_initdecomp_2d

  !===============================================================================
  subroutine wav_initdecomp_3d(nz, iodesc)

    integer        , intent(in)  :: nz
    type(io_desc_t), intent(out) :: iodesc

    integer          :: k,n
    integer, pointer :: dof3d(:)

    allocate(dof3d(nz*nseal_cpl))
    dof3d = 0
    n = 0
    do k = 1,nz
      do jsea = 1,nseal_cpl
        call init_get_isea(isea, jsea)
        ix = mapsf(isea,1)     ! global ix
        iy = mapsf(isea,2)     ! global iy
        n = n+1
        dof3d(n) = ((iy-1)*nx + ix) + (k-1)*nx*ny ! local index : global index
      end do
    end do
    call pio_initdecomp(wav_pio_subsystem, PIO_REAL, (/nx,ny,nz/), dof3d, iodesc)

    deallocate(dof3d)
  end subroutine wav_initdecomp_3d
end module w3iogoncmd_pio
