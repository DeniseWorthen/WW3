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

  subroutine wav_pio_init(gcomp, rc)
#ifdef CESMCOUPLED
   use shr_pio_mod, only: shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat
#endif
    use ESMF         , only : ESMF_GridComp, ESMF_UtilStringUpperCase, ESMF_VM, ESMF_FAILURE
    use ESMF         , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_ERROR
    use NUOPC        , only : NUOPC_CompAttributeGet
    use wav_kind_mod , only : CL=>SHR_KIND_CL, CS=>SHR_KIND_CS
    use w3adatmd     , only : mpi_comm_wave
    use w3odatmd     , only : naproc, iaproc
    use wav_shr_mod  , only : chkerr

    ! input/output arguments
    type(ESMF_GridComp), intent(in)    :: gcomp
    integer            , intent(out)   :: rc

    integer           :: pio_ioformat
    integer           :: pio_numiotasks
    integer           :: pio_stride
    integer           :: pio_rearranger
    integer           :: pio_root
    character(len=CS) :: cvalue
    logical           :: isPresent, isSet
    integer           :: my_task
    integer           :: master_task
    ! TODO: get right logunit
    integer           :: logunit = 6
    character(len=CS) :: subname='wav_pio_init'
    character(*), parameter :: u_FILE_u = &                  !< a character string for an ESMF log message
       __FILE__
    !-------------------------------------------------------------------------------
    rc = ESMF_SUCCESS

#ifdef CESMCOUPLED
    ! TODO: needs testing
    wav_pio_subsystem => shr_pio_getiosys(inst_name)
    pio_iotype =  shr_pio_getiotype(inst_name)
    if ((pio_iotype==PIO_IOTYPE_NETCDF).or.(pio_iotype==PIO_IOTYPE_PNETCDF)) then
      nmode0 = shr_pio_getioformat(inst_name)
    else
      nmode0 = 0
    endif

    call pio_seterrorhandling(wav_pio_subsystem, PIO_RETURN_ERROR)
#else
    my_task = iaproc - 1
    master_task = 0

    ! code lifted from CMEPS med_io_mod.F90
    ! query component specific PIO attributes
    ! pio_netcdf_format
    call NUOPC_CompAttributeGet(gcomp, name='pio_netcdf_format', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       cvalue = ESMF_UtilStringUpperCase(cvalue)
       if (trim(cvalue) .eq. 'CLASSIC') then
         pio_ioformat = 0
       else if (trim(cvalue) .eq. '64BIT_OFFSET') then
         pio_ioformat = PIO_64BIT_OFFSET
       else if (trim(cvalue) .eq. '64BIT_DATA') then
         pio_ioformat = PIO_64BIT_DATA
       else
         call ESMF_LogWrite(trim(subname)//': need to provide valid option for pio_ioformat ' &
              //'(CLASSIC|64BIT_OFFSET|64BIT_DATA)', ESMF_LOGMSG_ERROR)
         rc = ESMF_FAILURE
         return
       end if
    else
       cvalue = '64BIT_OFFSET'
       pio_ioformat = PIO_64BIT_OFFSET
    end if
    if (my_task == 0) write(logunit,*) trim(subname), ' : pio_netcdf_format = ', trim(cvalue), pio_ioformat

    ! pio_typename
    call NUOPC_CompAttributeGet(gcomp, name='pio_typename', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       cvalue = ESMF_UtilStringUpperCase(cvalue)
       if (trim(cvalue) .eq. 'NETCDF') then
         pio_iotype = PIO_IOTYPE_NETCDF
       else if (trim(cvalue) .eq. 'PNETCDF') then
         pio_iotype = PIO_IOTYPE_PNETCDF
       else if (trim(cvalue) .eq. 'NETCDF4C') then
         pio_iotype = PIO_IOTYPE_NETCDF4C
       else if (trim(cvalue) .eq. 'NETCDF4P') then
         pio_iotype = PIO_IOTYPE_NETCDF4P
       else
         call ESMF_LogWrite(trim(subname)//': need to provide valid option for pio_typename ' &
              //'(NETCDF|PNETCDF|NETCDF4C|NETCDF4P)', ESMF_LOGMSG_ERROR)
         rc = ESMF_FAILURE
         return
       end if
    else
       cvalue = 'NETCDF'
       pio_iotype = PIO_IOTYPE_NETCDF
    end if
    if (my_task == 0) write(logunit,*) trim(subname), ' : pio_typename = ', trim(cvalue), pio_iotype

    ! pio_root
    call NUOPC_CompAttributeGet(gcomp, name='pio_root', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_root
       if (pio_root < 0) then
          pio_root = 1
       endif
       pio_root = min(pio_root, naproc-1)
    else
       pio_root = 1
    end if
    if (my_task == 0) write(logunit,*) trim(subname), ' : pio_root = ', pio_root

    ! pio_stride
    call NUOPC_CompAttributeGet(gcomp, name='pio_stride', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_stride
    else
       pio_stride = -99
    end if
    if (my_task == 0) write(logunit,*) trim(subname), ' : pio_stride = ', pio_stride

    ! pio_numiotasks
    call NUOPC_CompAttributeGet(gcomp, name='pio_numiotasks', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_numiotasks
    else
       pio_numiotasks = -99
    end if
    if (my_task == 0) write(logunit,*) trim(subname), ' : pio_numiotasks = ', pio_numiotasks

    ! check for parallel IO, it requires at least two io pes
    if (naproc > 1 .and. pio_numiotasks == 1 .and. &
       (pio_iotype .eq. PIO_IOTYPE_PNETCDF .or. pio_iotype .eq. PIO_IOTYPE_NETCDF4P)) then
       pio_numiotasks = 2
       pio_stride = min(pio_stride, naproc/2)
       if (my_task == 0) then
          write(logunit,*) ' parallel io requires at least two io pes - following parameters are updated:'
          write(logunit,*) trim(subname), ' : pio_stride = ', pio_stride
          write(logunit,*) trim(subname), ' : pio_numiotasks = ', pio_numiotasks
       end if
    endif

    ! check/set/correct io pio parameters
    if (pio_stride > 0 .and. pio_numiotasks < 0) then
       pio_numiotasks = max(1, naproc/pio_stride)
       if (my_task == 0) write(logunit,*) trim(subname), ' : update pio_numiotasks = ', pio_numiotasks
    else if(pio_numiotasks > 0 .and. pio_stride < 0) then
       pio_stride = max(1, naproc/pio_numiotasks)
       if (my_task == 0) write(logunit,*) trim(subname), ' : update pio_stride = ', pio_stride
    else if(pio_numiotasks < 0 .and. pio_stride < 0) then
       pio_stride = max(1,naproc/4)
       pio_numiotasks = max(1,naproc/pio_stride)
       if (my_task == 0) write(logunit,*) trim(subname), ' : update pio_numiotasks = ', pio_numiotasks
       if (my_task == 0) write(logunit,*) trim(subname), ' : update pio_stride = ', pio_stride
    end if
    if (pio_stride == 1) then
       pio_root = 0
    endif

    if (pio_root + (pio_stride)*(pio_numiotasks-1) >= naproc .or. &
       pio_stride <= 0 .or. pio_numiotasks <= 0 .or. pio_root < 0 .or. pio_root > naproc-1) then
       if (naproc < 100) then
          pio_stride = max(1, naproc/4)
       else if(naproc < 1000) then
          pio_stride = max(1, naproc/8)
       else
          pio_stride = max(1, naproc/16)
       end if
       if(pio_stride > 1) then
          pio_numiotasks = naproc/pio_stride
          pio_root = min(1, naproc-1)
       else
          pio_numiotasks = naproc
          pio_root = 0
       end if
       if (my_task == 0) then
          write(logunit,*) 'pio_stride, iotasks or root out of bounds - resetting to defaults:'
          write(logunit,*) trim(subname), ' : pio_root = ', pio_root
          write(logunit,*) trim(subname), ' : pio_stride = ', pio_stride
          write(logunit,*) trim(subname), ' : pio_numiotasks = ', pio_numiotasks
       end if
    end if

    ! pio_rearranger
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearranger', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       cvalue = ESMF_UtilStringUpperCase(cvalue)
       if (trim(cvalue) .eq. 'BOX') then
         pio_rearranger = PIO_REARR_BOX
       else if (trim(cvalue) .eq. 'SUBSET') then
         pio_rearranger = PIO_REARR_SUBSET
       else
         call ESMF_LogWrite(trim(subname)//': need to provide valid option for pio_rearranger (BOX|SUBSET)', ESMF_LOGMSG_ERROR)
         rc = ESMF_FAILURE
         return
       end if
    else
       cvalue = 'SUBSET'
       pio_rearranger = PIO_REARR_SUBSET
    end if
    if (my_task == 0) write(logunit,*) trim(subname), ' : pio_rearranger = ', trim(cvalue), pio_rearranger

    ! init PIO
    if (my_task == 0) then
      write(logunit,*) trim(subname),' calling pio init'
      write(logunit,*) trim(subname), ' : pio_root = ', pio_root
      write(logunit,*) trim(subname), ' : pio_stride = ', pio_stride
      write(logunit,*) trim(subname), ' : pio_numiotasks = ', pio_numiotasks
    end if

    call pio_init(my_task, MPI_COMM_WAVE, pio_numiotasks, master_task, pio_stride, pio_rearranger, &
         wav_pio_subsystem, base=pio_root)
    call pio_seterrorhandling(wav_pio_subsystem, PIO_RETURN_ERROR)
#endif
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
