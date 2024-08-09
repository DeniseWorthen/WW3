module wav_restart_mod

  use w3parall          , only : init_get_isea
  use w3gdatmd          , only : nk, nx, ny, nspec, mapsf, mapsta, nsea
  use w3adatmd          , only : mpi_comm_wave
  use w3wdatmd          , only : va
  use w3odatmd          , only : iaproc, naproc
  use wav_import_export , only : nseal_cpl
  use w3iogoncmd_pio    , only : wav_initdecomp, handle_err       !TODO: move the wav_pio
  use pio
  use netcdf

  implicit none

  private
  ! used/reused in module
  integer             :: isea, jsea, ix, iy, ierr
  integer             :: len_s, len_m, len_p, len_k
  character(len=1024) :: fname

  type(iosystem_desc_t) :: wav_pio_subsystem
  type(file_desc_t)     :: pioid
  type(var_desc_t)      :: varid
  type(io_desc_t)       :: iodesc3d

  integer :: pio_iotype

  public :: write_restart

  !===============================================================================
contains
  !===============================================================================

  subroutine write_restart (fname)
    use w3odatmd   , only : time_origin, calendar_name, elapsed_secs

    character(len=*), intent(in) :: fname

    ! local variables
    character(len=12)    :: vname

    integer :: timid, xtid, ytid, ztid, j, ierr
    integer :: dimid(3)
    real,allocatable  :: varout(:,:)

    ! pio
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

    call pio_init(my_task, MPI_COMM_WAVE, numiotasks, master_task, istride, rearranger, &
         wav_pio_subsystem, base=basetask)
    call pio_seterrorhandling(wav_pio_subsystem, PIO_RETURN_ERROR)

    ! create the netcdf file
    ierr = pio_createfile(wav_pio_subsystem, pioid, pio_iotype, trim(fname), pio_clobber)
    call handle_err(ierr, 'pio_create')

    ierr = pio_def_dim(pioid,    'nx',    nx, xtid)
    ierr = pio_def_dim(pioid,    'ny',    ny, ytid)
    ierr = pio_def_dim(pioid, 'nspec', nspec, ztid)
    ierr = pio_def_dim(pioid, 'time', PIO_UNLIMITED, timid)

    ! define the time variable
    ierr = pio_def_var(pioid, 'time', PIO_DOUBLE, (/timid/), varid)
    call handle_err(ierr,'def_timevar')
    ierr = pio_put_att(pioid, varid, 'units', trim(time_origin))
    call handle_err(ierr,'def_time_units')
    ierr = pio_put_att(pioid, varid, 'calendar', trim(calendar_name))
    call handle_err(ierr,'def_time_calendar')

    vname = 'va'
    dimid = (/xtid, ytid, ztid/)
    ierr = pio_def_var(pioid, trim(vname), PIO_REAL, dimid, varid)
    call handle_err(ierr, 'define variable '//trim(vname))

    ! end variable definitions
    ierr = pio_enddef(pioid)
    call handle_err(ierr, 'end variable definition')

    call wav_initdecomp(nz=nspec, iodesc=iodesc3d)

    !WRITEBUFF(1:NSPEC) = VA(1:NSPEC,JSEA)

    allocate(varout(nspec,nseal_cpl))
    varout = 0.0
    do j = 1,nseal_cpl
      call init_get_isea(isea, jsea)
      varout(:,jsea) = va(:,jsea)
    end do

    ierr = pio_inq_varid(pioid,  trim(vname), varid)
    call handle_err(ierr, 'inquire variable '//trim(vname))
    call pio_setframe(pioid, varid, int(1,kind=PIO_OFFSET_KIND))
    call pio_write_darray(pioid, varid, iodesc3d, varout, ierr)
    deallocate(varout)

    call pio_freedecomp(pioid,iodesc3d)
    call pio_closefile(pioid)
  end subroutine write_restart

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

end module wav_restart_mod
