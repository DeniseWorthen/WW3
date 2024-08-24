module wav_restart_mod

  use w3parall      , only : init_get_isea
  use w3adatmd      , only : mpi_comm_wave, nsealm
  use w3gdatmd      , only : nth, nk, nx, ny, nspec, nseal, nsea
  use wav_pio_mod   , only : pio_iotype, wav_pio_subsystem
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
  type(io_desc_t)   :: iodesc3dk

  integer(kind=Pio_Offset_Kind) :: lframe

  public :: write_restart
  public :: read_restart

  !===============================================================================
contains
  !===============================================================================

  subroutine write_restart (fname, va_in, map_in)

    use w3gdatmd , only : mapsf
    use w3odatmd , only : time_origin, calendar_name, elapsed_secs
    !debug
    use w3odatmd , only : iaproc

    real            , intent(in) :: va_in(1:nspec,0:nsealm)
    integer         , intent(in) :: map_in(ny,nx)
    character(len=*), intent(in) :: fname

    ! local variables
    character(len=12)    :: vname
    integer              :: timid, xtid, ytid, ztid, ierr
    integer              :: ik, ith, ix, iy, kk, nseal_cpl
    integer              :: isea, jsea
    integer              :: dimid(4)
    integer, allocatable :: locmap(:)
    real, allocatable    :: locva(:,:)
    !-------------------------------------------------------------------------------

#ifdef W3_PDLIB
    nseal_cpl = nseal - ng
#else
    nseal_cpl = nseal
#endif
    allocate(locmap(1:nseal_cpl))
    allocate(locva(1:nseal_cpl,1:nspec))

    ! create the netcdf file
    lframe = 1
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
    call handle_err(ierr, 'define _FillValue '//trim(vname))

    vname = 'mapsta'
    ierr = pio_def_var(pioid, trim(vname), PIO_INT, (/xtid, ytid, timid/), varid)
    call handle_err(ierr, 'define variable '//trim(vname))
    ierr = pio_put_att(pioid, varid, '_FillValue', nf90_fill_int)
    call handle_err(ierr, 'define _FillValue '//trim(vname))

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
    locmap(:) = 0
    do jsea = 1,nseal_cpl
      call init_get_isea(isea, jsea)
      ix = mapsf(isea,1)
      iy = mapsf(isea,2)
      locmap(jsea) = map_in(iy,ix)
    end do

    ! write PE local map
    vname = 'mapsta'
    ierr = pio_inq_varid(pioid,  trim(vname), varid)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, varid, int(1,kind=Pio_Offset_Kind))
    call pio_write_darray(pioid, varid, iodesc2dint, locmap, ierr)
    call handle_err(ierr, 'put variable '//trim(vname))

    ! write va
    locva(:,:) = 0.0
    do jsea = 1,nseal_cpl
      kk = 0
      do ik = 1,nk
        do ith = 1,nth
          kk = kk + 1
          locva(jsea,kk) = va_in(kk,jsea)
        end do
      end do
    end do

    vname = 'va'
    ierr = pio_inq_varid(pioid,  trim(vname), varid)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, varid, int(1,kind=PIO_OFFSET_KIND))
    call pio_write_darray(pioid, varid, iodesc3dk, locva, ierr)
    call handle_err(ierr, 'put variable '//trim(vname))

    call pio_syncfile(pioid)
    call pio_freedecomp(pioid, iodesc2dint)
    call pio_freedecomp(pioid, iodesc3dk)
    call pio_closefile(pioid)

  end subroutine write_restart

  !===============================================================================
  subroutine read_restart (fname, va_out, map_out)

    use mpi
    !use mpi_f08  !? why doesn't this work
    use w3adatmd    , only : mpi_comm_wave
    use w3gdatmd    , only : mapsf, mapst2, sig, nseal
    use w3wdatmd    , only : time, tlev, tice, trho, tic1, tic5, wlv, asf, ice, fpis
    ! debug
    !use w3odatmd , only : iaproc
    !use w3gdatmd , only : xgrd, ygrd

    real            , intent(out) :: va_out(1:nspec,0:nsealm)
    integer         , intent(out) :: map_out(ny,nx)
    character(len=*), intent(in)  :: fname

    ! local variables
    integer              :: ik, ith, ix, iy, kk, nseal_cpl
    integer              :: isea, jsea
    character(len=12)    :: vname
    integer              :: ierr
    logical              :: exists
    integer              :: global_input(nsea), global_output(nsea)
    integer              :: ifill
    real                 :: rfill
    real, allocatable    :: valoc(:,:)
    integer, allocatable :: maploc2d(:,:)
    integer, allocatable :: maploc(:)
    ! debug
    !integer :: ix, iy
    !-------------------------------------------------------------------------------

#ifdef W3_PDLIB
    nseal_cpl = nseal - ng
#else
    nseal_cpl = nseal
#endif
    allocate(valoc(1:nseal_cpl,1:nspec))
    allocate(maploc2d(1:ny,1:nx))
    allocate(maploc(1:nseal_cpl))
    valoc(:,:) = 0.0
    maploc2d(:,:) = 0
    maploc(:) = 0

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
        lframe = 1
        ierr = pio_openfile(wav_pio_subsystem, pioid, pio_iotype, trim(fname), pio_nowrite)
        call handle_err(ierr, 'open file '//trim(fname))
      else
        !error out
      end if
    end if

    ! initialize the decomp
    call wav_pio_initdecomp(iodesc2dint, use_int=.true.)
    call wav_pio_initdecomp(nspec, iodesc3dk)

    vname = 'va'
    ierr = pio_inq_varid(pioid, trim(vname), varid)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, varid, lframe)
    call pio_read_darray(pioid, varid, iodesc3dk, valoc, ierr)
    call handle_err(ierr, 'get variable '//trim(vname))
    ierr = pio_get_att(pioid, varid, "_FillValue", rfill)
    call handle_err(ierr, 'get variable _FillValue'//trim(vname))

    va_out = 0.0
    do jsea = 1,nseal_cpl
      kk = 0
      do ik = 1,nk
        do ith = 1,nth
          kk = kk + 1
          if (valoc(jsea,kk) .ne. rfill) then
            va_out(kk,jsea) = valoc(jsea,kk)
          end if
        end do
      end do
    end do
    !if (trim(fname) .eq. 'ufs.cpld.ww3.r.2021-03-22-64800.nc') then
    !  print *,'XXX ',minval(va_out),maxval(va_out)
    !end if

    vname = 'mapsta'
    ierr = pio_inq_varid(pioid, trim(vname), varid)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, varid, lframe)
    call pio_read_darray(pioid, varid, iodesc2dint, maploc, ierr)
    call handle_err(ierr, 'get variable '//trim(vname))
    ierr = pio_get_att(pioid, varid, "_FillValue", ifill)
    call handle_err(ierr, 'get variable _FillValue'//trim(vname))
    !if (trim(fname) .eq. 'ufs.cpld.ww3.r.2021-03-22-64800.nc') then
    !  print *,'YYY ',minval(maploc),maxval(maploc)
    !end if

    ! mapsta is local on this PE
    ! locmap(:) = 0
    ! do jsea = 1,nseal_cpl
    !   call init_get_isea(isea, jsea)
    !   ix = mapsf(isea,1)
    !   iy = mapsf(isea,2)
    !   locmap(jsea) = maploc(iy,ix)
    ! end do

    ! !if (trim(fname) .eq. 'ufs.cpld.ww3.r.2021-03-22-64800.nc') then
    !   do jsea = 1,nseal_cpl
    !     call init_get_isea(isea, jsea)
    !     ix = mapsf(isea,1)
    !     iy = mapsf(isea,2)
    !     write(100+iaproc,'(5i8,2f12.5)')jsea,isea,iy,ix,maploc(jsea),xgrd(iy,ix),ygrd(iy,ix)
    !   end do
    !   !do iy = 1,ny
    !   !  do ix = 1,nx
    !   !    write(100+iaproc,'(3i8,2f12.5)')iy,ix,maploc(iy,ix),xgrd(iy,ix),ygrd(iy,ix)
    !   !  end do
    !   !end do
    ! end if

    ! fill global array with PE local values
    global_input = 0
    global_output = 0
    do jsea = 1,nseal_cpl
      call init_get_isea(isea, jsea)
      ix = mapsf(isea,1)
      iy = mapsf(isea,2)
      if (maploc(jsea) .ne. ifill) then
        global_input(isea) = maploc(jsea)
      end if
      !if (trim(fname) .eq. 'ufs.cpld.ww3.r.2021-03-22-64800.nc') then
      !  write(200+iaproc,'(5i6,2f12.5)')jsea,isea,iy,ix,global_input(isea),xgrd(iy,ix),ygrd(iy,ix)
      !end if
    end do

    !do isea = 1,nsea
    !  ix = mapsf(isea,1)
    !  iy = mapsf(isea,2)
    !  global_input(isea) = maploc(iy,ix)
    !end do

    ! if (trim(fname) .eq. 'ufs.cpld.ww3.r.2021-03-22-64800.nc') then
    !   do isea = 1,nsea
    !     ix = mapsf(isea,1)
    !     iy = mapsf(isea,2)
    !     write(300+iaproc,'(4i8,2f12.5)')isea,iy,ix,global_input(isea),xgrd(iy,ix),ygrd(iy,ix)
    !   end do
    ! end if

    ! reduce across all PEs to create global array
    call MPI_AllReduce(global_input, global_output, nsea, MPI_INTEGER, MPI_SUM, MPI_COMM_WAVE, ierr)

    ! fill global array on each PE
    maploc2d = 0
    do isea = 1,nsea
      ix = mapsf(isea,1)
      iy = mapsf(isea,2)
      maploc2d(iy,ix) = global_output(isea)
    end do

    ! if (trim(fname) .eq. 'ufs.cpld.ww3.r.2021-03-22-64800.nc') then
    !   do isea = 1,nsea
    !     ix = mapsf(isea,1)
    !     iy = mapsf(isea,2)
    !     write(400+iaproc,*)isea,iy,ix,maploc2d(iy,ix)
    !   end do
    ! end if

    map_out = mod(maploc2d+2,8) - 2
    mapst2 = (maploc2d-map_out)/8

    call pio_syncfile(pioid)
    call pio_freedecomp(pioid, iodesc2dint)
    call pio_freedecomp(pioid, iodesc3dk)
    call pio_closefile(pioid)

  end subroutine read_restart
end module wav_restart_mod
