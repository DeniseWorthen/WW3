module wav_restart_mod

  use w3parall          , only : init_get_isea
  use w3adatmd          , only : mpi_comm_wave
  use w3gdatmd          , only : nth, nk, nx, ny, nspec, mapsta, nseal, nsea
  use wav_import_export , only : nseal_cpl
  use wav_pio_mod       , only : pio_iotype, wav_pio_subsystem
  use wav_pio_mod       , only : handle_err, wav_pio_initdecomp
  use pio
  use netcdf

  implicit none

  private

  type(file_desc_t) :: pioid
  type(var_desc_t)  :: varid
  type(io_desc_t)   :: iodesc2dint
  type(io_desc_t)   :: iodesc3dk

  public :: write_restart
  public :: read_restart

  !===============================================================================
contains
  !===============================================================================

  subroutine write_restart (fname, a)

    use w3gdatmd , only : mapsf
    use w3odatmd , only : time_origin, calendar_name, elapsed_secs

    real            , intent(in) :: a(nth,nk,0:nseal)
    character(len=*), intent(in) :: fname

    ! local variables
    character(len=12) :: vname
    integer           :: timid, xtid, ytid, ztid, ierr
    integer           :: ik, ith, ix, iy, kk
    integer           :: isea, jsea
    integer           :: dimid(4)
    integer           :: lmap(1:nseal_cpl)
    real              :: varout(1:nseal_cpl,1:nspec)

    ! create the netcdf file
    pioid%fh = -1
    ierr = pio_createfile(wav_pio_subsystem, pioid, pio_iotype, trim(fname), pio_clobber)
    call handle_err(ierr, 'pio_create')

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

    vname = 'mapsta'
    ierr = pio_def_var(pioid, trim(vname), PIO_INT, (/xtid, ytid, timid/), varid)
    call handle_err(ierr, 'define variable '//trim(vname))
    ierr = pio_put_att(pioid, varid, '_FillValue', nf90_fill_int)

    ! end variable definitions
    ierr = pio_enddef(pioid)
    call handle_err(ierr, 'end variable definition')

    ! initialize the decomp
    call wav_pio_initdecomp(iodesc2dint, use_int=.true.)
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

    ! write mapsta
    vname = 'mapsta'
    ierr = pio_inq_varid(pioid,  trim(vname), varid)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, varid, int(1,kind=Pio_Offset_Kind))
    call pio_write_darray(pioid, varid, iodesc2dint, lmap, ierr)
    call handle_err(ierr, 'put variable '//trim(vname))

    ! write va
    varout(:,:) = 0.0
    do jsea = 1,nseal_cpl
      kk = 0
      do ik = 1,nk
        do ith = 1,nth
          kk = kk + 1
          varout(jsea,kk) = a(ith,ik,jsea)
        end do
      end do
    end do

    vname = 'va'
    ierr = pio_inq_varid(pioid,  trim(vname), varid)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, varid, int(1,kind=PIO_OFFSET_KIND))
    call pio_write_darray(pioid, varid, iodesc3dk, varout, ierr)
    call handle_err(ierr, 'put variable '//trim(vname))

    !call pio_syncfile(pioid)
    call pio_freedecomp(pioid, iodesc2dint)
    call pio_freedecomp(pioid, iodesc3dk)
    call pio_closefile(pioid)

  end subroutine write_restart

  subroutine read_restart (fname, vaout)

    real            , intent(out) :: vaout(1:nspec,1:nsea)
    character(len=*), intent(in)  :: fname

    ! local variables
    character(len=12) :: vname
    integer           :: ierr
    logical           :: exists
    ! decompositions are real, need to make an integer one to write mapsta as int
    real              :: maptmp(ny,nx)
    real              :: atmp(nth,nk,0:nseal)

    ! open the netcdf file
    inquire(file = trim(fname), exist=exists)
    if (exists) then
      pioid%fh = -1
      ierr = pio_openfile(wav_pio_subsystem, pioid, pio_iotype, trim(fname), pio_write)
      call handle_err(ierr, 'open file '//trim(fname))
    else
      !error out
    end if

    ! initialize the decomp
    call wav_pio_initdecomp(iodesc2dint)
    call wav_pio_initdecomp(nspec, iodesc3dk)

    vname = 'va'
    ierr = pio_inq_varid(pioid, trim(vname), varid)
    call pio_read_darray(pioid, varid, iodesc3dk, atmp, ierr)
    call handle_err(ierr, 'get variable '//trim(vname))

    vname = 'mapsta'
    ierr = pio_inq_varid(pioid, trim(vname), varid)
    call pio_read_darray(pioid, varid, iodesc2dint, maptmp, ierr)
    call handle_err(ierr, 'get variable '//trim(vname))

    !mapsta = int(maptmp,4)

    call pio_freedecomp(pioid, iodesc2dint)
    call pio_freedecomp(pioid, iodesc3dk)
    call pio_closefile(pioid)

  end subroutine read_restart
end module wav_restart_mod
