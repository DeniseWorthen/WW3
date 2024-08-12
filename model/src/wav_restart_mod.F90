module wav_restart_mod

  use w3parall          , only : init_get_isea
  use w3adatmd          , only : mpi_comm_wave
  use w3gdatmd          , only : nth, nk, nx, ny, nspec, mapsta, nseal
  use wav_import_export , only : nseal_cpl
  use wav_pio_mod       , only : pio_iotype, wav_pio_subsystem
  use wav_pio_mod       , only : handle_err, wav_pio_initdecomp
  use pio
  use netcdf

  implicit none

  private

  type(file_desc_t) :: pioid
  type(var_desc_t)  :: varid
  type(io_desc_t)   :: iodesc2d
  type(io_desc_t)   :: iodesc3dk

  public :: write_restart
  !public :: read_restart

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
    real, allocatable :: varout(:,:)
    real, allocatable :: lmap(:)

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
    ierr = pio_put_att(pioid, varid, 'units', 'unitless')
    ierr = pio_put_att(pioid, varid, 'long_name', 'map status')

    ! end variable definitions
    ierr = pio_enddef(pioid)
    call handle_err(ierr, 'end variable definition')

    ! initialize the decomp
    call wav_pio_initdecomp(iodesc2d)
    call wav_pio_initdecomp(nspec, iodesc3dk)
    !print '(a,8i8)','XXX ',nsea,nseal,nsealm,nseal_cpl,lbound(va,1),ubound(va,1),lbound(va,2),ubound(va,2)
    !print '(a,5i8)','XXX ',nsea,nseal,nsealm,nseal_cpl,ubound(a,3)
    !WRITEBUFF(1:NSPEC) = VA(1:NSPEC,JSEA)

    ! write the time
    ierr = pio_inq_varid(pioid,  'time', varid)
    call handle_err(ierr, 'inquire variable time ')
    ierr = pio_put_var(pioid, varid, (/1/), real(elapsed_secs,8))
    call handle_err(ierr, 'put time')

    ! mapsta is global
    allocate(lmap(1:nseal_cpl))
    lmap = 0.0
    do jsea = 1,nseal_cpl
      call init_get_isea(isea, jsea)
      ix = mapsf(isea,1)
      iy = mapsf(isea,2)
      lmap(jsea) = real(mapsta(iy,ix),4)
    end do

    ! write mapsta
    vname = 'mapsta'
    ierr = pio_inq_varid(pioid,  trim(vname), varid)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, varid, int(1,kind=Pio_Offset_Kind))
    call pio_write_darray(pioid, varid, iodesc2d, lmap, ierr)
    call handle_err(ierr, 'put variable '//trim(vname))
    deallocate(lmap)

    ! write va
    allocate(varout(1:nseal_cpl,1:nspec))
    varout = 0.0
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
    deallocate(varout)

    !call pio_syncfile(pioid)
    call pio_freedecomp(pioid, iodsec2d)
    call pio_freedecomp(pioid, iodesc3dk)
    call pio_closefile(pioid)

  end subroutine write_restart

  ! subroutine read_restart (fname, vaout)

  !   real            , intent(in) :: vaout(nspec,0:nseal)
  !   character(len=*), intent(in) :: fname


  !   ! ! decompositions are real, need to make an integer one to write mapsta as int
  !   ! real, allocatable :: lmap(:)

  !   ! vname = 'va'
  !   ! dimid = (/xtid, ytid, ztid, timid/)
  !   ! ierr = pio_def_var(pioid, trim(vname), PIO_REAL, dimid, varid)
  !   ! call handle_err(ierr, 'define variable '//trim(vname))

  !   ! vname = 'mapsta'
  !   ! ierr = pio_def_var(pioid, trim(vname), PIO_INT, (/xtid, ytid, timid/), varid)
  !   ! call handle_err(ierr, 'def_mapsta')
  !   ! ierr = pio_put_att(pioid, varid, 'units', 'unitless')
  !   ! ierr = pio_put_att(pioid, varid, 'long_name', 'map status')

  !   ! ! end variable definitions
  !   ! ierr = pio_enddef(pioid)
  !   ! call handle_err(ierr, 'end variable definition')

  ! end subroutine read_restart

end module wav_restart_mod
