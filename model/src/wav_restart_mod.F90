module wav_restart_mod

  use w3parall          , only : init_get_isea
  use w3adatmd          , only : mpi_comm_wave, nsealm
  use w3gdatmd          , only : nth, nk, nx, ny, nspec, nseal, nsea
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

  subroutine write_restart (fname, va, mapsta)

    use w3gdatmd , only : mapsf
    use w3odatmd , only : time_origin, calendar_name, elapsed_secs
    ! eventually, and mapsta will come from module but for now, I want
    ! to test a read/write cycle w/o changing the internal values
    !use w3wdatmd, only : va

    real            , intent(in) :: va(1:nspec,0:nsealm)
    integer         , intent(in) :: mapsta(ny,nx)
    character(len=*), intent(in) :: fname

    ! local variables
    character(len=12) :: vname
    integer           :: timid, xtid, ytid, ztid, ierr
    integer           :: ik, ith, ix, iy, kk
    integer           :: isea, jsea
    integer           :: dimid(4)
    integer           :: lmap(1:nseal_cpl)
    real              :: va_out(1:nseal_cpl,1:nspec)

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
    ierr = pio_put_att(pioid, varid, '_FillValue', nf90_fill_float)

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

    ! write mapsta PE local mapsta
    vname = 'mapsta'
    ierr = pio_inq_varid(pioid,  trim(vname), varid)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, varid, int(1,kind=Pio_Offset_Kind))
    call pio_write_darray(pioid, varid, iodesc2dint, lmap, ierr)
    call handle_err(ierr, 'put variable '//trim(vname))

    ! write va
    va_out(:,:) = 0.0
    do jsea = 1,nseal_cpl
      kk = 0
      do ik = 1,nk
        do ith = 1,nth
          kk = kk + 1
          va_out(jsea,kk) = va(kk,jsea)
        end do
      end do
    end do

    vname = 'va'
    ierr = pio_inq_varid(pioid,  trim(vname), varid)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, varid, int(1,kind=PIO_OFFSET_KIND))
    call pio_write_darray(pioid, varid, iodesc3dk, va_out, ierr)
    call handle_err(ierr, 'put variable '//trim(vname))

    call pio_syncfile(pioid)
    call pio_freedecomp(pioid, iodesc2dint)
    call pio_freedecomp(pioid, iodesc3dk)
    call pio_closefile(pioid)

  end subroutine write_restart

  subroutine read_restart (fname, va_out, mapsta_out)

    use mpi
    !use mpi_f08  !? why doesn't this work
    use w3adatmd , only : mpi_comm_wave
    use w3gdatmd , only : mapsf, sig
    use w3wdatmd , only : time, tlev, tice, trho, tic1, tic5, wlv, asf, ice, fpis

    ! debug
    use w3odatmd , only : iaproc

    real            , intent(out) :: va_out(1:nspec,0:nsealm)
    integer         , intent(out) :: mapsta_out(ny,nx)
    character(len=*), intent(in)  :: fname

    ! local variables
    integer           :: ik, ith, kk
    integer           :: isea, jsea
    character(len=12) :: vname
    integer           :: ierr
    logical           :: exists
    real              :: vatmp(1:nseal_cpl,1:nspec)
    integer           :: global_input(nsea), global_output(nsea)
    ! debug
    integer :: ix, iy

    ! open the netcdf file; should this be different routine?
    if (trim(fname)  == 'none') then
      !fill needed fields and return
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
      return
    else
      ! all times are restart times
      tlev = time
      tice = time
      trho = time
      tic1 = time
      tic5 = time
      inquire(file=trim(fname), exist=exists)
      if (exists) then
        pioid%fh = -1
        ierr = pio_openfile(wav_pio_subsystem, pioid, pio_iotype, trim(fname), pio_nowrite)
        call handle_err(ierr, 'open file '//trim(fname))
      else
        !error out
      end if
    end if
    !debug
    va_out = -999.0
    mapsta_out = -99
    ! initialize the decomp
    call wav_pio_initdecomp(iodesc2dint, use_int=.true.)
    call wav_pio_initdecomp(nspec, iodesc3dk)

    vname = 'va'
    ierr = pio_inq_varid(pioid, trim(vname), varid)
    call pio_read_darray(pioid, varid, iodesc3dk, vatmp, ierr)
    call handle_err(ierr, 'get variable '//trim(vname))

    do jsea = 1,nseal_cpl
      kk = 0
      do ik = 1,nk
        do ith = 1,nth
          kk = kk + 1
          va_out(kk,jsea) = vatmp(jsea,kk)
        end do
      end do
    end do
    ! mapsta is global
    ! lmap(:) = 0
    ! do jsea = 1,nseal_cpl
    !   call init_get_isea(isea, jsea)
    !   ix = mapsf(isea,1)
    !   iy = mapsf(isea,2)
    !   lmap(jsea) = mapsta(iy,ix)
    ! end do

    vname = 'mapsta'
    ierr = pio_inq_varid(pioid, trim(vname), varid)
    call pio_read_darray(pioid, varid, iodesc2dint, mapsta_out, ierr)
    call handle_err(ierr, 'get variable '//trim(vname))
    ! do jsea = 1,nseal_cpl
    !   call init_get_isea(isea, jsea)
    !   ix = mapsf(isea,1)
    !   iy = mapsf(isea,2)
    !   write(100+iaproc,*)jsea,isea,iy,ix,mapsta_out(iy,ix)
    ! end do

    ! fill global array with PE local values
    global_input = 0
    global_output = 0
    do jsea = 1,nseal_cpl
      call init_get_isea(isea, jsea)
      ix = mapsf(isea,1)
      iy = mapsf(isea,2)
      if (mapsta_out(iy,ix) .ne. nf90_fill_int) then
        global_input(isea) = mapsta_out(iy,ix)
      end if
      !write(200+iaproc,*)jsea,isea,iy,ix,mapsta_out(iy,ix),global_input(isea)
    end do

    ! do isea = 1,nsea
    !   ix = mapsf(isea,1)
    !   iy = mapsf(isea,2)
    !   write(300+iaproc,*)isea,iy,ix,global_input(isea)
    ! end do

    ! reduce across all PEs
    call MPI_AllReduce(global_input, global_output, nsea, MPI_INTEGER, MPI_SUM, MPI_COMM_WAVE, ierr)

    ! fill global array on each PE
    mapsta_out = 0
    do isea = 1,nsea
      ix = mapsf(isea,1)
      iy = mapsf(isea,2)
      mapsta_out(iy,ix) = global_output(isea)
    end do

    ! do isea = 1,nsea
    !   ix = mapsf(isea,1)
    !   iy = mapsf(isea,2)
    !   write(400+iaproc,*)isea,iy,ix,mapsta_out(iy,ix)
    ! end do

    call pio_syncfile(pioid)
    call pio_freedecomp(pioid, iodesc2dint)
    call pio_freedecomp(pioid, iodesc3dk)
    call pio_closefile(pioid)

  end subroutine read_restart
end module wav_restart_mod
