module wav_restart_mod

  use w3parall          , only : init_get_isea
  use w3gdatmd          , only : nth, nk, nx, ny, nspec, mapsf, mapsta, nsea, nseal
  use w3adatmd          , only : nsealm, mpi_comm_wave
  use w3odatmd          , only : iaproc, naproc
  use wav_import_export , only : nseal_cpl
  use w3iogoncmd_pio    , only : handle_err       !TODO: move the wav_pio
  use pio
  use netcdf

  implicit none

  private
  ! used/reused in module
  integer             :: isea, jsea, ix, iy, ierr
  character(len=1024) :: fname

  type(iosystem_desc_t) :: wav_pio_subsystem
  type(file_desc_t)     :: pioid
  type(var_desc_t)      :: varid
  type(io_desc_t)       :: iodesc

  public :: write_restart

  !===============================================================================
contains
  !===============================================================================

  subroutine write_restart (fname, a)

    use w3odatmd   , only : time_origin, calendar_name, elapsed_secs

    ! input/output variables
    real            , intent(in) :: a(nth,nk,0:nseal)
    character(len=*), intent(in) :: fname

    ! local variables
    character(len=12) :: vname
    integer           :: timid, xtid, ytid, ztid, ierr
    integer           :: ik, ith, ix, iy, kk
    integer           :: dimid(4)
    real,allocatable  :: varout(:,:)

    ! pio
    integer :: pio_iotype
    integer :: nprocs
    integer :: istride
    integer :: basetask
    integer :: numiotasks
    integer :: rearranger
    integer :: my_task
    integer :: master_task
    integer :: status

    ! TODO: for now, hardwire  the io system
    pioid%fh = -1
    pio_iotype = PIO_IOTYPE_PNETCDF
    nprocs = naproc
    my_task = iaproc - 1
    master_task = 0
    istride = 4
    basetask = 1
    numiotasks = max((nprocs-basetask)/istride,1)
    !numiotasks = 2
    rearranger = PIO_REARR_BOX
    !rearranger = PIO_REARR_SUBSET

    call pio_init(my_task, MPI_COMM_WAVE, numiotasks, master_task, istride, rearranger, &
         wav_pio_subsystem, base=basetask)
    call pio_seterrorhandling(wav_pio_subsystem, PIO_RETURN_ERROR)

    ! create the netcdf file
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

    ! end variable definitions
    ierr = pio_enddef(pioid)
    call handle_err(ierr, 'end variable definition')

    ! initialize the decomp
    call wav_initdecomp_3d(nz=nspec, iodesc=iodesc)
    !print '(a,8i8)','XXX ',nsea,nseal,nsealm,nseal_cpl,lbound(va,1),ubound(va,1),lbound(va,2),ubound(va,2)
    !print '(a,5i8)','XXX ',nsea,nseal,nsealm,nseal_cpl,ubound(a,3)
    !WRITEBUFF(1:NSPEC) = VA(1:NSPEC,JSEA)

    ! write the time
    ierr = pio_inq_varid(pioid,  'time', varid)
    call handle_err(ierr, 'inquire variable time ')
    ierr = pio_put_var(pioid, varid, (/1/), real(elapsed_secs,8))
    call handle_err(ierr, 'put time')

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
    call pio_write_darray(pioid, varid, iodesc, varout, ierr)
    deallocate(varout)

    !call pio_syncfile(pioid)
    call pio_freedecomp(pioid, iodesc)
    call pio_closefile(pioid)

  end subroutine write_restart

  !===============================================================================
  subroutine wav_initdecomp_3d(nz, iodesc)

    type(io_desc_t)      :: iodesc
    integer , intent(in) :: nz

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

end module wav_restart_mod
