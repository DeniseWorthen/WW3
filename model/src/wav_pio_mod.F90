!> @file wav_pio
!!
!> @brief Manage WAV PIO
!!
!> @author Denise.Worthen@noaa.gov
!> @date 08-02-2024

module wav_pio_mod

  use w3gdatmd          , only : nk, nx, ny, mapsf
  use w3parall          , only : init_get_isea
  use wav_import_export , only : nseal_cpl
  use pio
  use netcdf

  implicit none

  private

  interface wav_pio_initdecomp
    module procedure wav_pio_initdecomp_2d
    module procedure wav_pio_initdecomp_3d
  end interface wav_pio_initdecomp

  type(iosystem_desc_t) :: wav_pio_subsystem
  integer               :: pio_iotype

  public :: wav_pio_init
  public :: pio_iotype
  public :: wav_pio_subsystem
  public :: wav_pio_initdecomp
  public :: handle_err

  !===============================================================================
contains
  !===============================================================================

  subroutine wav_pio_init()

    use w3adatmd , only : mpi_comm_wave
    use w3odatmd , only : naproc, iaproc

    ! pio
    integer :: nprocs
    integer :: istride
    integer :: basetask
    integer :: numiotasks
    integer :: rearranger
    integer :: my_task
    integer :: master_task

    !-------------------------------------------------------------------------------

    ! TODO: for now, hardwire  the io system
    pio_iotype = PIO_IOTYPE_PNETCDF
    nprocs = naproc
    my_task = iaproc - 1
    master_task = 0
    istride = 4
    basetask = 1
    numiotasks = max((nprocs-basetask)/istride,1)
    !numiotasks = 2
    rearranger = PIO_REARR_BOX
    print '(a,4i8)','SETUP ',nprocs,iaproc,my_task,numiotasks, nseal_cpl

    call pio_init(my_task, MPI_COMM_WAVE, numiotasks, master_task, istride, rearranger, &
         wav_pio_subsystem, base=basetask)
    call pio_seterrorhandling(wav_pio_subsystem, PIO_RETURN_ERROR)

  end subroutine wav_pio_init

  !===============================================================================
  subroutine wav_pio_initdecomp_2d(iodesc, use_int)

    type(io_desc_t),           intent(out) :: iodesc
    logical        , optional, intent(in)  :: use_int

    ! local variables
    integer          :: n, isea, jsea, ix, iy
    integer, pointer :: dof2d(:)
    logical          :: luse_int

    luse_int = .false.
    if (present(use_int)) luse_int = use_int

    allocate(dof2d(nseal_cpl))
    dof2d = 0
    n = 0
    do jsea = 1,nseal_cpl
      call init_get_isea(isea, jsea)
      ix = mapsf(isea,1)                 ! global ix
      iy = mapsf(isea,2)                 ! global iy
      n = n+1
      dof2d(n) = (iy-1)*nx + ix          ! local index : global index
    end do
    if (luse_int) then
      call pio_initdecomp(wav_pio_subsystem, PIO_INT,  (/nx,ny/), dof2d, iodesc)
    else
      call pio_initdecomp(wav_pio_subsystem, PIO_REAL, (/nx,ny/), dof2d, iodesc)
    end if

  end subroutine wav_pio_initdecomp_2d

  !===============================================================================
  subroutine wav_pio_initdecomp_3d(nz, iodesc)

    integer ,         intent(in)  :: nz
    type(io_desc_t) , intent(out) :: iodesc

    ! local variables
    integer          :: n, k, isea, jsea, ix, iy
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
  end subroutine wav_pio_initdecomp_3d

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

end module wav_pio_mod
