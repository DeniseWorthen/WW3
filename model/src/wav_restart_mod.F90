!> @file wav_restart_mod
!!
!> @brief Handle WW3 restart files as netCDF using PIO
!!
!> @author Denise.Worthen@noaa.gov
!> @date 08-26-2024
module wav_restart_mod

  use w3parall      , only : init_get_isea
  use w3adatmd      , only : nsealm
  use w3gdatmd      , only : nth, nk, nx, ny, mapsf, nspec, nseal, nsea
  use w3odatmd      , only : ndso, iaproc, addrstflds, rstfldlist, rstfldcnt
  use w3wdatmd      , only : ice
  use wav_pio_mod   , only : pio_iotype, pio_ioformat, wav_pio_subsystem
  use wav_pio_mod   , only : handle_err, wav_pio_initdecomp
#ifdef W3_PDLIB
    use yowNodepool , only : ng
#endif
  use pio
  use netcdf

  implicit none

  private

  type(file_desc_t) :: pioid
  type(var_desc_t)  :: varid
  type(io_desc_t)   :: iodesc2dint
  type(io_desc_t)   :: iodesc2d
  type(io_desc_t)   :: iodesc3dk

  integer(kind=Pio_Offset_Kind) :: frame

  public :: write_restart
  public :: read_restart

  ! used/reused in module
  character(len=12) :: vname
  integer           :: ik, ith, ix, iy, kk, nseal_cpl, isea, jsea, ierr, i

  !===============================================================================
contains
  !===============================================================================
  !> Write a WW3 restart file
  !!
  !! @details Called by w3wavemd to write a restart file at a given frequency or
  !! time
  !!
  !! @param[in]     fname    the time-stamped file name
  !! @param[in]     va       the va array
  !! @param[in]     mapsta   the mapsta + 8*mapst2 array
  !!
  !> author DeniseWorthen@noaa.gov
  !> @date 08-26-2024
  subroutine write_restart (fname, va, mapsta)

    use w3odatmd , only : time_origin, calendar_name, elapsed_secs

    real            , intent(in) :: va(1:nspec,0:nsealm)
    integer         , intent(in) :: mapsta(ny,nx)
    character(len=*), intent(in) :: fname

    ! local variables
    integer              :: timid, xtid, ytid, ztid
    integer              :: nmode
    integer              :: dimid(4)
    real   , allocatable :: lva(:,:)
    integer, allocatable :: lmap(:)
    real   , allocatable :: lvar(:)
    !-------------------------------------------------------------------------------

#ifdef W3_PDLIB
    nseal_cpl = nseal - ng
#else
    nseal_cpl = nseal
#endif
    allocate(lmap(1:nseal_cpl))
    allocate(lva(1:nseal_cpl,1:nspec))
    allocate(lvar(1:nseal_cpl))

    ! create the netcdf file
    frame = 1
    pioid%fh = -1
    nmode = pio_clobber
    ! only applies to classic NETCDF files.
    if (pio_iotype == PIO_IOTYPE_NETCDF .or. pio_iotype == PIO_IOTYPE_PNETCDF) then
      nmode = ior(nmode,pio_ioformat)
    endif
    ierr = pio_createfile(wav_pio_subsystem, pioid, pio_iotype, trim(fname), nmode)
    call handle_err(ierr, 'pio_create')
    if (iaproc == 1) write(ndso,'(a)')' Writing restart file '//trim(fname)

    ierr = pio_def_dim(pioid,    'nx',    nx, xtid)
    ierr = pio_def_dim(pioid,    'ny',    ny, ytid)
    ierr = pio_def_dim(pioid, 'nspec', nspec, ztid)
    ierr = pio_def_dim(pioid,  'time', PIO_UNLIMITED, timid)

    !define the time variable
    ierr = pio_def_var(pioid, 'time', PIO_DOUBLE, (/timid/), varid)
    call handle_err(ierr,'def_timevar')
    ierr = pio_put_att(pioid, varid, 'units', trim(time_origin))
    call handle_err(ierr,'def_time_units')
    ierr = pio_put_att(pioid, varid, 'calendar', trim(calendar_name))
    call handle_err(ierr,'def_time_calendar')

    vname = 'va'
    dimid = (/xtid, ytid, ztid, timid/)
    ierr = pio_def_var(pioid, trim(vname), PIO_REAL, dimid, varid)
    call handle_err(ierr, 'define variable '//trim(vname))
    ierr = pio_put_att(pioid, varid, '_FillValue', nf90_fill_float)
    call handle_err(ierr, 'define _FillValue '//trim(vname))

    vname = 'mapsta'
    ierr = pio_def_var(pioid, trim(vname), PIO_INT, (/xtid, ytid, timid/), varid)
    call handle_err(ierr, 'define variable '//trim(vname))
    ierr = pio_put_att(pioid, varid, '_FillValue', nf90_fill_int)
    call handle_err(ierr, 'define _FillValue '//trim(vname))

    ! define any requested additional fields
    if (addrstflds) then
      do i = 1,rstfldcnt
        vname = trim(rstfldlist(i))
        ierr = pio_def_var(pioid, trim(vname), PIO_REAL, (/xtid, ytid, timid/), varid)
        call handle_err(ierr, 'define variable '//trim(vname))
        ierr = pio_put_att(pioid, varid, '_FillValue', nf90_fill_float)
        call handle_err(ierr, 'define _FillValue '//trim(vname))
      end do
    end if
    ! end variable definitions
    ierr = pio_enddef(pioid)
    call handle_err(ierr, 'end variable definition')

    ! initialize the decomp
    call wav_pio_initdecomp(iodesc2dint, use_int=.true.)
    if (addrstflds) call wav_pio_initdecomp(iodesc2d)
    call wav_pio_initdecomp(nspec, iodesc3dk)

    ! write the time
    ierr = pio_inq_varid(pioid,  'time', varid)
    call handle_err(ierr, 'inquire variable time ')
    ierr = pio_put_var(pioid, varid, (/1/), real(elapsed_secs,8))
    call handle_err(ierr, 'put time')

    ! mapsta is global
    lmap(:) = 0
    do jsea = 1,nseal_cpl
      call init_get_isea(isea, jsea)
      ix = mapsf(isea,1)
      iy = mapsf(isea,2)
      lmap(jsea) = mapsta(iy,ix)
    end do

    ! write PE local map
    vname = 'mapsta'
    ierr = pio_inq_varid(pioid,  trim(vname), varid)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, varid, int(1,kind=Pio_Offset_Kind))
    call pio_write_darray(pioid, varid, iodesc2dint, lmap, ierr)
    call handle_err(ierr, 'put variable '//trim(vname))

    ! write va
    lva(:,:) = 0.0
    do jsea = 1,nseal_cpl
      kk = 0
      do ik = 1,nk
        do ith = 1,nth
          kk = kk + 1
          lva(jsea,kk) = va(kk,jsea)
        end do
      end do
    end do

    vname = 'va'
    ierr = pio_inq_varid(pioid,  trim(vname), varid)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, varid, int(1,kind=PIO_OFFSET_KIND))
    call pio_write_darray(pioid, varid, iodesc3dk, lva, ierr)
    call handle_err(ierr, 'put variable '//trim(vname))

    ! write requested additional fields
    if (addrstflds) then
      do i = 1,rstfldcnt
        vname = trim(rstfldlist(i))
        ! TODO: make generic routine (in=ice, out=lvar)
        if (vname == 'ice') then
          lvar(:) = 0.0
          do jsea = 1,nseal_cpl
            call init_get_isea(isea, jsea)
            lvar(jsea) = ice(isea)
          end do
        end if

        ! write PE local field
        ierr = pio_inq_varid(pioid,  trim(vname), varid)
        call handle_err(ierr, 'inquire variable '//trim(vname))
        call pio_setframe(pioid, varid, int(1,kind=Pio_Offset_Kind))
        call pio_write_darray(pioid, varid, iodesc2d, lvar, ierr)
        call handle_err(ierr, 'put variable '//trim(vname))
      end do
    end if

    call pio_syncfile(pioid)
    if (addrstflds) call pio_freedecomp(pioid, iodesc2d)
    call pio_freedecomp(pioid, iodesc2dint)
    call pio_freedecomp(pioid, iodesc3dk)
    call pio_closefile(pioid)

  end subroutine write_restart

  !===============================================================================
  !> Read a WW3 restart file
  !!
  !> @details Called by w3init to read a restart file which is known to exist or to
  !! initialize a set of variables when the filename is "none".
  !!
  !! @param[in]     fname     the time-stamped file name
  !! @param[out]    va        the va array, optional
  !! @param[out]    mapsta    the mapsta array, optional
  !! @param[inout]  mapst2    the mapst2 array, optional
  !!
  !> author DeniseWorthen@noaa.gov
  !> @date 08-26-2024
  subroutine read_restart (fname, va, mapsta, mapst2)

    use mpi_f08
    use w3adatmd    , only : mpi_comm_wave
    use w3gdatmd    , only : sig
    use w3wdatmd    , only : time, tlev, tice, trho, tic1, tic5, wlv, asf, fpis

    character(len=*)  , intent(in)    :: fname
    real   , optional , intent(out)   :: va(1:nspec,0:nsealm)
    integer, optional , intent(out)   :: mapsta(ny,nx)
    integer, optional , intent(inout) :: mapst2(ny,nx)

    ! local variables
    type(MPI_Comm)       :: wave_communicator  ! needed for mpi_f08
    real                 :: global_input(nsea), global_output(nsea)
    integer              :: ifill
    real                 :: rfill
    real   , allocatable :: lva(:,:)
    real   , allocatable :: lvar(:)
    integer, allocatable :: lmap(:)
    integer, allocatable :: lmap2d(:,:)
    integer, allocatable :: st2init(:,:)
    !-------------------------------------------------------------------------------

    ! cold start, set initial values and return.
    if (trim(fname)  == 'none') then
      tlev(1) = -1
      tlev(2) =  0
      tice(1) = -1
      tice(2) =  0
      trho(1) = -1
      trho(2) =  0
      tic1(1) = -1
      tic1(2) =  0
      tic5(1) = -1
      tic5(2) =  0
      wlv     =  0.
      ice     =  0.
      asf     =  1.
      fpis    =  sig(nk)
      if (iaproc == 1) write(ndso,'(a)')' Initializing WW3 at rest '
      return
    end if

    ! read a netcdf restart
    wave_communicator%mpi_val = MPI_COMM_WAVE
#ifdef W3_PDLIB
    nseal_cpl = nseal - ng
#else
    nseal_cpl = nseal
#endif
    allocate(lva(1:nseal_cpl,1:nspec))
    allocate(lvar(1:nseal_cpl))
    allocate(lmap2d(1:ny,1:nx))
    allocate(st2init(1:ny,1:nx))
    allocate(lmap(1:nseal_cpl))
    lva(:,:) = 0.0
    lvar(:) = 0.0
    lmap2d(:,:) = 0
    lmap(:) = 0
    ! save a copy of initial mapst2 from mod_def
    st2init = mapst2

    ! all times are restart times
    tlev = time
    tice = time
    trho = time
    tic1 = time
    tic5 = time
    frame = 1
    ierr = pio_openfile(wav_pio_subsystem, pioid, pio_iotype, trim(fname), pio_nowrite)
    call handle_err(ierr, 'open file '//trim(fname))
    if (iaproc == 1) write(ndso,'(a)')' Reading restart file '//trim(fname)

    ! initialize the decomp
    call wav_pio_initdecomp(iodesc2dint, use_int=.true.)
    if (addrstflds) call wav_pio_initdecomp(iodesc2d)
    call wav_pio_initdecomp(nspec, iodesc3dk)

    vname = 'va'
    ierr = pio_inq_varid(pioid, trim(vname), varid)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, varid, frame)
    call pio_read_darray(pioid, varid, iodesc3dk, lva, ierr)
    call handle_err(ierr, 'get variable '//trim(vname))
    ierr = pio_get_att(pioid, varid, "_FillValue", rfill)
    call handle_err(ierr, 'get variable _FillValue'//trim(vname))

    va = 0.0
    do jsea = 1,nseal_cpl
      kk = 0
      do ik = 1,nk
        do ith = 1,nth
          kk = kk + 1
          if (lva(jsea,kk) .ne. rfill) then
            va(kk,jsea) = lva(jsea,kk)
          end if
        end do
      end do
    end do

    vname = 'mapsta'
    ierr = pio_inq_varid(pioid, trim(vname), varid)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, varid, frame)
    call pio_read_darray(pioid, varid, iodesc2dint, lmap, ierr)
    call handle_err(ierr, 'get variable '//trim(vname))
    ierr = pio_get_att(pioid, varid, "_FillValue", ifill)
    call handle_err(ierr, 'get variable _FillValue'//trim(vname))

    ! fill global array with PE local values
    global_input = 0.0
    global_output = 0.0
    do jsea = 1,nseal_cpl
      call init_get_isea(isea, jsea)
      if (lmap(jsea) .ne. ifill) then
        global_input(isea) = real(lmap(jsea))
      end if
    end do
    ! reduce across all PEs to create global array
    call MPI_AllReduce(global_input, global_output, nsea, MPI_REAL, MPI_SUM, wave_communicator, ierr)

    ! fill global array on each PE
    lmap2d = 0
    do isea = 1,nsea
      ix = mapsf(isea,1)
      iy = mapsf(isea,2)
      lmap2d(iy,ix) = int(global_output(isea))
    end do

    mapsta = mod(lmap2d+2,8) - 2
    mapst2 = st2init + (lmap2d-mapsta)/8

    ! read additional restart fields
    if (addrstflds) then
      do i = 1,rstfldcnt
        vname = trim(rstfldlist(i))
        ierr = pio_inq_varid(pioid, trim(vname), varid)
        call handle_err(ierr, 'inquire variable '//trim(vname))
        call pio_setframe(pioid, varid, frame)
        call pio_read_darray(pioid, varid, iodesc2d, lvar, ierr)
        call handle_err(ierr, 'get variable '//trim(vname))
        ierr = pio_get_att(pioid, varid, "_FillValue", rfill)
        call handle_err(ierr, 'get variable _FillValue'//trim(vname))

        ! fill global array with PE local values
        global_input = 0.0
        global_output = 0.0
        do jsea = 1,nseal_cpl
          call init_get_isea(isea, jsea)
          if (lvar(jsea) .ne. rfill) then
            global_input(isea) = lvar(jsea)
          end if
        end do
        ! reduce across all PEs to create global array
        call MPI_AllReduce(global_input, global_output, nsea, MPI_REAL, MPI_SUM, wave_communicator, ierr)

        if (vname == 'ice') then
          ! fill global array on each PE
          ! TODO : make generic routine (in=global_ouput, out=ice)
          ice = 0.0
          do isea = 1,nsea
            ice(isea) = global_output(isea)
          end do
        end if
      end do
    end if

    call pio_syncfile(pioid)
    if (addrstflds) call pio_freedecomp(pioid, iodesc2d)
    call pio_freedecomp(pioid, iodesc2dint)
    call pio_freedecomp(pioid, iodesc3dk)
    call pio_closefile(pioid)

  end subroutine read_restart
end module wav_restart_mod
