!> @file wav_pio
!!
!> @brief Manage PIO for WW3
!!
!> @author Denise.Worthen@noaa.gov
!> @date 08-02-2024
module wav_pio_mod

  use w3gdatmd    , only : nk, nx, ny, mapsf
  use w3parall    , only : init_get_isea
  use w3gdatmd    , only : nseal
  use pio
  use netcdf
#ifdef W3_PDLIB
  use yowNodepool , only : ng
#endif
  implicit none

  private

  interface wav_pio_initdecomp
    module procedure wav_pio_initdecomp_2d
    module procedure wav_pio_initdecomp_3d
  end interface wav_pio_initdecomp

  integer                        :: pio_iotype
  integer                        :: pio_ioformat
  type(iosystem_desc_t), pointer :: wav_pio_subsystem

  ! Restart write options, set from WAV attributes in wav_pio_init.
  ! The defaults keep the original behavior.
  logical :: pio_restart_nofill     = .false. !< set PIO_NOFILL on restart files
  integer :: pio_restart_syncfreq   = 0       !< pio_syncfile every N va variables (0 = never)
  logical :: pio_restart_final_sync = .true.  !< pio_syncfile after all restart fields are written

  public :: wav_pio_init
  public :: pio_iotype
  public :: pio_ioformat
  public :: wav_pio_subsystem
  public :: pio_restart_nofill
  public :: pio_restart_syncfreq
  public :: pio_restart_final_sync
  public :: wav_pio_initdecomp
  public :: handle_err

  !===============================================================================
contains
  !===============================================================================
  !> Configure PIO for WW3
  !!
  !> @details Use either CESM shr code or configuration variables to configure PIO.
  !! This configuration code is lifted from CMEPS.
  !!
  !! @param       gcomp             an ESMF_GridComp object
  !! @param       mpi_comm          the MPI communicator
  !! @param[in]   stdout            the logfile unit on the root_task
  !! @param[in]   numprocs          naproc/nthrds
  !! @param[out]  rc                a return code
  !!
  !> @author Denise.Worthen@noaa.gov
  !> @date 08-02-2024
  subroutine wav_pio_init(gcomp, mpi_comm, stdout, numprocs, rc)

#ifdef CESMCOUPLED
    use shr_pio_mod, only : shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat
#endif
    use ESMF         , only : ESMF_GridComp, ESMF_UtilStringUpperCase, ESMF_VM, ESMF_FAILURE
    use ESMF         , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_ERROR
    use NUOPC        , only : NUOPC_CompAttributeGet
    use wav_kind_mod , only : CL=>SHR_KIND_CL, CS=>SHR_KIND_CS
    use w3odatmd     , only : iaproc
    use wav_shr_mod  , only : chkerr

    ! input/output arguments
    type(ESMF_GridComp), intent(in)    :: gcomp
    integer            , intent(in)    :: mpi_comm
    integer            , intent(in)    :: stdout
    integer            , intent(in)    :: numprocs
    integer            , intent(out)   :: rc

    integer           :: pio_numiotasks
    integer           :: pio_stride
    integer           :: pio_rearranger
    integer           :: pio_root
    integer           :: pio_debug_level
    integer           :: pio_buffer_limit                     ! MB
    integer           :: comm_type, fcd
    logical           :: enable_hs_c2i, enable_isend_c2i
    logical           :: enable_hs_i2c, enable_isend_i2c
    integer           :: max_pend_req_c2i, max_pend_req_i2c
    integer           :: ierr
    character(len=CS) :: cvalue
    logical           :: isPresent, isSet
    integer           :: my_task, master_task
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
    if (my_task == 0) write(stdout,*) trim(subname), ' : pio_netcdf_format = ', trim(cvalue), pio_ioformat

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
    if (my_task == 0) write(stdout,*) trim(subname), ' : pio_typename = ', trim(cvalue), pio_iotype

    ! pio_root
    call NUOPC_CompAttributeGet(gcomp, name='pio_root', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_root
       if (pio_root < 0) then
          pio_root = 1
       endif
       pio_root = min(pio_root, numprocs-1)
    else
       pio_root = 1
    end if
    if (my_task == 0) write(stdout,*) trim(subname), ' : pio_root = ', pio_root

    ! pio_stride
    call NUOPC_CompAttributeGet(gcomp, name='pio_stride', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_stride
    else
       pio_stride = -99
    end if
    if (my_task == 0) write(stdout,*) trim(subname), ' : pio_stride = ', pio_stride

    ! pio_numiotasks
    call NUOPC_CompAttributeGet(gcomp, name='pio_numiotasks', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
       read(cvalue,*) pio_numiotasks
    else
       pio_numiotasks = -99
    end if
    if (my_task == 0) write(stdout,*) trim(subname), ' : pio_numiotasks = ', pio_numiotasks

    ! check for parallel IO, it requires at least two io pes
    if (numprocs > 1 .and. pio_numiotasks == 1 .and. &
       (pio_iotype .eq. PIO_IOTYPE_PNETCDF .or. pio_iotype .eq. PIO_IOTYPE_NETCDF4P)) then
       pio_numiotasks = 2
       pio_stride = min(pio_stride, numprocs/2)
       if (my_task == 0) then
          write(stdout,*) ' parallel io requires at least two io pes - following parameters are updated:'
          write(stdout,*) trim(subname), ' : pio_stride = ', pio_stride
          write(stdout,*) trim(subname), ' : pio_numiotasks = ', pio_numiotasks
       end if
    endif

    ! check/set/correct io pio parameters
    if (pio_stride > 0 .and. pio_numiotasks < 0) then
       pio_numiotasks = max(1, numprocs/pio_stride)
       if (my_task == 0) write(stdout,*) trim(subname), ' : update pio_numiotasks = ', pio_numiotasks
    else if(pio_numiotasks > 0 .and. pio_stride < 0) then
       pio_stride = max(1, numprocs/pio_numiotasks)
       if (my_task == 0) write(stdout,*) trim(subname), ' : update pio_stride = ', pio_stride
    else if(pio_numiotasks < 0 .and. pio_stride < 0) then
       pio_stride = max(1,numprocs/4)
       pio_numiotasks = max(1,numprocs/pio_stride)
       if (my_task == 0) write(stdout,*) trim(subname), ' : update pio_numiotasks = ', pio_numiotasks
       if (my_task == 0) write(stdout,*) trim(subname), ' : update pio_stride = ', pio_stride
    end if
    if (pio_stride == 1) then
       pio_root = 0
    endif

    if (pio_root + (pio_stride)*(pio_numiotasks-1) >= numprocs .or. &
       pio_stride <= 0 .or. pio_numiotasks <= 0 .or. pio_root < 0 .or. pio_root > numprocs-1) then
       if (numprocs < 100) then
          pio_stride = max(1, numprocs/4)
       else if(numprocs < 1000) then
          pio_stride = max(1, numprocs/8)
       else
          pio_stride = max(1, numprocs/16)
       end if
       if(pio_stride > 1) then
          pio_numiotasks = numprocs/pio_stride
          pio_root = min(1, numprocs-1)
       else
          pio_numiotasks = numprocs
          pio_root = 0
       end if
       if (my_task == 0) then
          write(stdout,*) 'pio_stride, iotasks or root out of bounds - resetting to defaults:'
          write(stdout,*) trim(subname), ' : pio_root = ', pio_root
          write(stdout,*) trim(subname), ' : pio_stride = ', pio_stride
          write(stdout,*) trim(subname), ' : pio_numiotasks = ', pio_numiotasks
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
    if (my_task == 0) write(stdout,*) trim(subname), ' : pio_rearranger = ', trim(cvalue), pio_rearranger

    ! init PIO
    if (my_task == 0) then
      write(stdout,*) trim(subname),' calling pio init'
      write(stdout,*) trim(subname), ' : pio_root = ', pio_root
      write(stdout,*) trim(subname), ' : pio_stride = ', pio_stride
      write(stdout,*) trim(subname), ' : pio_numiotasks = ', pio_numiotasks
    end if

    allocate(wav_pio_subsystem)
    call pio_init(my_task, mpi_comm, pio_numiotasks, master_task, pio_stride, pio_rearranger, &
         wav_pio_subsystem, base=pio_root)

    ! pio_buffer_limit (MB): amount of data PIO collects before rearranging and writing.
    ! If not set, PIO's default is used.
    call NUOPC_CompAttributeGet(gcomp, name='pio_buffer_limit', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
      read(cvalue,*) pio_buffer_limit
      if (pio_buffer_limit <= 0) then
        call ESMF_LogWrite(trim(subname)//': pio_buffer_limit must be > 0 (MB)', ESMF_LOGMSG_ERROR)
        rc = ESMF_FAILURE
        return
      end if
      call pio_set_buffer_size_limit(int(pio_buffer_limit,PIO_OFFSET_KIND)*1024_PIO_OFFSET_KIND*1024_PIO_OFFSET_KIND)
      if (my_task == 0) write(stdout,*) trim(subname), ' : pio_buffer_limit (MB) = ', pio_buffer_limit
    else
      if (my_task == 0) write(stdout,*) trim(subname), ' : pio_buffer_limit not set, using the PIO default'
    end if

    ! pio_rearr_comm_type (P2P|COLL). If not set, PIO's default rearranger options
    ! are used and pio_set_rearr_opts is not called. If set, the other options below
    ! default to PIO's initial values.
    call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_type', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
      cvalue = ESMF_UtilStringUpperCase(cvalue)
      if (trim(cvalue) .eq. 'P2P') then
        comm_type = PIO_REARR_COMM_P2P
      else if (trim(cvalue) .eq. 'COLL') then
        comm_type = PIO_REARR_COMM_COLL
      else
        call ESMF_LogWrite(trim(subname)//': need to provide valid option for pio_rearr_comm_type (P2P|COLL)', ESMF_LOGMSG_ERROR)
        rc = ESMF_FAILURE
        return
      end if
      if (my_task == 0) write(stdout,*) trim(subname), ' : pio_rearr_comm_type = ', trim(cvalue)

      ! flow control direction
      call NUOPC_CompAttributeGet(gcomp, name='pio_rearr_comm_fcd', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (isPresent .and. isSet) then
        cvalue = ESMF_UtilStringUpperCase(cvalue)
      else
        cvalue = '2DDISABLE'
      end if
      if (trim(cvalue) .eq. '2DENABLE') then
        fcd = PIO_REARR_COMM_FC_2D_ENABLE
      else if (trim(cvalue) .eq. 'COMP2IO') then
        fcd = PIO_REARR_COMM_FC_1D_COMP2IO
      else if (trim(cvalue) .eq. 'IO2COMP') then
        fcd = PIO_REARR_COMM_FC_1D_IO2COMP
      else if (trim(cvalue) .eq. '2DDISABLE') then
        fcd = PIO_REARR_COMM_FC_2D_DISABLE
      else
        call ESMF_LogWrite(trim(subname)//': need to provide valid option for pio_rearr_comm_fcd '// &
             '(2DENABLE|COMP2IO|IO2COMP|2DDISABLE)', ESMF_LOGMSG_ERROR)
        rc = ESMF_FAILURE
        return
      end if
      if (my_task == 0) write(stdout,*) trim(subname), ' : pio_rearr_comm_fcd = ', trim(cvalue)

      ! flow control details, compute tasks to io tasks and io tasks to compute tasks
      call get_logical_attribute(gcomp, 'pio_rearr_comm_enable_hs_comp2io', .false., enable_hs_c2i, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call get_logical_attribute(gcomp, 'pio_rearr_comm_enable_isend_comp2io', .false., enable_isend_c2i, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call get_integer_attribute(gcomp, 'pio_rearr_comm_max_pend_req_comp2io', 0, max_pend_req_c2i, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call get_logical_attribute(gcomp, 'pio_rearr_comm_enable_hs_io2comp', .false., enable_hs_i2c, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call get_logical_attribute(gcomp, 'pio_rearr_comm_enable_isend_io2comp', .false., enable_isend_i2c, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call get_integer_attribute(gcomp, 'pio_rearr_comm_max_pend_req_io2comp', 0, max_pend_req_i2c, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (my_task == 0) then
        write(stdout,*) trim(subname), ' : pio_rearr_comm comp2io hs, isend, max_pend_req = ', &
             enable_hs_c2i, enable_isend_c2i, max_pend_req_c2i
        write(stdout,*) trim(subname), ' : pio_rearr_comm io2comp hs, isend, max_pend_req = ', &
             enable_hs_i2c, enable_isend_i2c, max_pend_req_i2c
      end if

      ierr = pio_set_rearr_opts(wav_pio_subsystem, comm_type, fcd, &
           enable_hs_c2i, enable_isend_c2i, max_pend_req_c2i, &
           enable_hs_i2c, enable_isend_i2c, max_pend_req_i2c)
      if (ierr /= PIO_NOERR) then
        call ESMF_LogWrite(trim(subname)//': pio_set_rearr_opts failed', ESMF_LOGMSG_ERROR)
        rc = ESMF_FAILURE
        return
      end if
    else
      if (my_task == 0) write(stdout,*) trim(subname), ' : pio_rearr_comm_type not set, using the PIO default rearranger options'
    end if

    ! restart write options
    call get_logical_attribute(gcomp, 'pio_restart_nofill', .false., pio_restart_nofill, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call get_integer_attribute(gcomp, 'pio_restart_syncfreq', 0, pio_restart_syncfreq, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (pio_restart_syncfreq < 0) then
      call ESMF_LogWrite(trim(subname)//': pio_restart_syncfreq must be >= 0', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
    end if
    call get_logical_attribute(gcomp, 'pio_restart_final_sync', .true., pio_restart_final_sync, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == 0) then
      write(stdout,*) trim(subname), ' : pio_restart_nofill = ', pio_restart_nofill
      write(stdout,*) trim(subname), ' : pio_restart_syncfreq = ', pio_restart_syncfreq
      write(stdout,*) trim(subname), ' : pio_restart_final_sync = ', pio_restart_final_sync
    end if

    ! PIO debug related options
    ! pio_debug_level
    call NUOPC_CompAttributeGet(gcomp, name='pio_debug_level', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
      read(cvalue,*) pio_debug_level
      if (pio_debug_level < 0 .or. pio_debug_level > 6) then
        call ESMF_LogWrite(trim(subname)//': need to provide valid option for pio_debug_level (0-6)', ESMF_LOGMSG_ERROR)
        rc = ESMF_FAILURE
        return
      end if
    else
      pio_debug_level = 0
    end if
    if (my_task == 0) write(stdout,*) trim(subname), ' : pio_debug_level = ', pio_debug_level

    ! set PIO debug level
    call pio_setdebuglevel(pio_debug_level)

    call pio_seterrorhandling(wav_pio_subsystem, PIO_RETURN_ERROR)
#endif
  end subroutine wav_pio_init

  !===============================================================================
  !> Get an optional logical WAV attribute (TRUE|FALSE)
  !!
  !! @param[in]   gcomp     an ESMF_GridComp object
  !! @param[in]   aname     the attribute name
  !! @param[in]   default   the value used if the attribute is not set
  !! @param[out]  avalue    the value
  !! @param[out]  rc        a return code
  subroutine get_logical_attribute(gcomp, aname, default, avalue, rc)

    use ESMF         , only : ESMF_GridComp, ESMF_UtilStringUpperCase
    use ESMF         , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LogWrite, ESMF_LOGMSG_ERROR
    use NUOPC        , only : NUOPC_CompAttributeGet
    use wav_kind_mod , only : CS=>SHR_KIND_CS
    use wav_shr_mod  , only : chkerr

    type(ESMF_GridComp), intent(in)  :: gcomp
    character(len=*)   , intent(in)  :: aname
    logical            , intent(in)  :: default
    logical            , intent(out) :: avalue
    integer            , intent(out) :: rc

    character(len=CS) :: cvalue
    logical           :: isPresent, isSet
    character(*), parameter :: u_FILE_u = &                  !< a character string for an ESMF log message
         __FILE__
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    avalue = default
    call NUOPC_CompAttributeGet(gcomp, name=trim(aname), value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
      cvalue = ESMF_UtilStringUpperCase(cvalue)
      if (trim(cvalue) .eq. 'TRUE') then
        avalue = .true.
      else if (trim(cvalue) .eq. 'FALSE') then
        avalue = .false.
      else
        call ESMF_LogWrite('wav_pio_init: need to provide valid option for '//trim(aname)//' (TRUE|FALSE)', &
             ESMF_LOGMSG_ERROR)
        rc = ESMF_FAILURE
        return
      end if
    end if

  end subroutine get_logical_attribute

  !===============================================================================
  !> Get an optional integer WAV attribute
  !!
  !! @param[in]   gcomp     an ESMF_GridComp object
  !! @param[in]   aname     the attribute name
  !! @param[in]   default   the value used if the attribute is not set
  !! @param[out]  avalue    the value
  !! @param[out]  rc        a return code
  subroutine get_integer_attribute(gcomp, aname, default, avalue, rc)

    use ESMF         , only : ESMF_GridComp, ESMF_SUCCESS
    use NUOPC        , only : NUOPC_CompAttributeGet
    use wav_kind_mod , only : CS=>SHR_KIND_CS
    use wav_shr_mod  , only : chkerr

    type(ESMF_GridComp), intent(in)  :: gcomp
    character(len=*)   , intent(in)  :: aname
    integer            , intent(in)  :: default
    integer            , intent(out) :: avalue
    integer            , intent(out) :: rc

    character(len=CS) :: cvalue
    logical           :: isPresent, isSet
    character(*), parameter :: u_FILE_u = &                  !< a character string for an ESMF log message
         __FILE__
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    avalue = default
    call NUOPC_CompAttributeGet(gcomp, name=trim(aname), value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) read(cvalue,*) avalue

  end subroutine get_integer_attribute

  !===============================================================================
  !> Define a decomposition for a 2d variable in WW3
  !!
  !! @param[out]  iodesc   the PIO decomposition handle
  !! @param[out]  use_int  define a decomposition for an integer array
  !!
  !> @author Denise.Worthen@noaa.gov
  !> @date 08-02-2024
  subroutine wav_pio_initdecomp_2d(iodesc, use_int)

    type(io_desc_t),           intent(out) :: iodesc
    logical        , optional, intent(in)  :: use_int

    ! local variables
    integer :: n, isea, jsea, ix, iy, nseal_cpl
    logical :: luse_int
    integer(kind=PIO_OFFSET_KIND) :: lnx,lny
    integer(kind=PIO_OFFSET_KIND), allocatable :: dof2d(:)
#ifdef W3_PDLIB
    nseal_cpl = nseal - ng
#else
    nseal_cpl = nseal
#endif
    luse_int = .false.
    if (present(use_int)) luse_int = use_int

    allocate(dof2d(nseal_cpl))
    dof2d = 0
    lnx = int(nx,PIO_OFFSET_KIND)
    lny = int(ny,PIO_OFFSET_KIND)

    n = 0
    do jsea = 1,nseal_cpl
      call init_get_isea(isea, jsea)
      ix = mapsf(isea,1)                 ! global ix
      iy = mapsf(isea,2)                 ! global iy
      n = n+1
      dof2d(n) = (iy-1)*lnx + ix         ! local index : global index
    end do

    ! optional: write the decomposition to a file for the standalone PIO benchmark
    call dump_decomp(dof2d)

    if (luse_int) then
      call pio_initdecomp(wav_pio_subsystem, PIO_INT,  (/nx,ny/), dof2d, iodesc)
    else
      call pio_initdecomp(wav_pio_subsystem, PIO_REAL, (/nx,ny/), dof2d, iodesc)
    end if
    deallocate(dof2d)

  end subroutine wav_pio_initdecomp_2d

  !===============================================================================
  !> Write the 2d restart decomposition to a file, once per run, if the environment
  !! variable WW3_PIO_DUMP_DECOMP=1
  !!
  !! @details Task 0 gathers each task's global indices (dof2d, in task order) and
  !! writes one unformatted stream file, ww3_restart_decomp_<ntasks>.bin, with
  !! ntasks, nx, ny, counts(ntasks), then the indices of each task in turn
  !! (all default integers). The file is read by the standalone PIO benchmark.
  !!
  !! @param[in]  dof2d   this task's global indices, as passed to pio_initdecomp
  subroutine dump_decomp(dof2d)

    use mpi_f08
    use w3adatmd , only : mpi_comm_wave
    use w3odatmd , only : ndso

    integer(kind=PIO_OFFSET_KIND), intent(in) :: dof2d(:)

    ! local variables
    logical, save        :: checked = .false.     ! the switch is checked once per run
    type(MPI_Comm)       :: wave_communicator      ! needed for mpi_f08
    integer              :: ntasks, my_task, ierr
    integer              :: idump, envstatus, nlocal, n, iunit
    character(len=16)    :: envvalue
    character(len=64)    :: fname
    integer, allocatable :: counts(:), displs(:)
    integer, allocatable :: ldof(:), gdof(:)
    !-------------------------------------------------------------------------------

    if (checked) return
    checked = .true.

    wave_communicator = MPI_COMM_WAVE
    call MPI_Comm_rank(wave_communicator, my_task, ierr)
    call MPI_Comm_size(wave_communicator, ntasks, ierr)

    ! task 0 reads the switch and passes it on
    idump = 0
    if (my_task == 0) then
      call get_environment_variable('WW3_PIO_DUMP_DECOMP', envvalue, status=envstatus)
      if (envstatus == 0 .and. trim(envvalue) == '1') idump = 1
    end if
    call MPI_Bcast(idump, 1, MPI_INTEGER, 0, wave_communicator, ierr)
    if (idump == 0) return

    ! the indices are written as default integers
    if (int(nx,PIO_OFFSET_KIND)*int(ny,PIO_OFFSET_KIND) > int(huge(1),PIO_OFFSET_KIND)) then
      if (my_task == 0) write(ndso,'(a)') ' dump_decomp: nx*ny too large for default integers, no file written'
      return
    end if

    nlocal = size(dof2d)
    allocate(ldof(max(1,nlocal)))
    if (nlocal > 0) ldof(1:nlocal) = int(dof2d(1:nlocal))

    ! counts per task, then the indices in task order
    allocate(counts(ntasks), displs(ntasks))
    call MPI_Gather(nlocal, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, wave_communicator, ierr)
    if (my_task == 0) then
      displs(1) = 0
      do n = 2,ntasks
        displs(n) = displs(n-1) + counts(n-1)
      end do
      allocate(gdof(max(1,sum(counts))))
    else
      allocate(gdof(1))
    end if
    call MPI_Gatherv(ldof, nlocal, MPI_INTEGER, gdof, counts, displs, MPI_INTEGER, 0, wave_communicator, ierr)

    if (my_task == 0) then
      write(fname,'(a,i0,a)') 'ww3_restart_decomp_', ntasks, '.bin'
      open(newunit=iunit, file=trim(fname), access='stream', form='unformatted', status='replace')
      write(iunit) ntasks, nx, ny
      write(iunit) counts
      write(iunit) gdof(1:sum(counts))
      close(iunit)
      write(ndso,'(a,i0,a,i0,a)') ' dump_decomp: wrote '//trim(fname)//' (', ntasks, ' tasks, ', &
           sum(counts), ' indices)'
    end if

    deallocate(counts, displs, ldof, gdof)

  end subroutine dump_decomp

  !===============================================================================
  !> Define a decomposition for a 3d variable in WW3
  !!
  !! @param[in]   nz       the non-spatial dimension
  !! @param[out]  iodesc   the PIO decomposition handle
  !!
  !> @author Denise.Worthen@noaa.gov
  !> @date 08-02-2024
  subroutine wav_pio_initdecomp_3d(nz, iodesc)

    integer ,         intent(in)  :: nz
    type(io_desc_t) , intent(out) :: iodesc

    ! local variables
    integer :: n, k, isea, jsea, ix, iy, nseal_cpl
    integer(kind=PIO_OFFSET_KIND) :: lnx,lny
    integer(kind=PIO_OFFSET_KIND), allocatable :: dof3d(:)
#ifdef W3_PDLIB
    nseal_cpl = nseal - ng
#else
    nseal_cpl = nseal
#endif
    allocate(dof3d(nz*nseal_cpl))

    dof3d = 0
    lnx = int(nx,PIO_OFFSET_KIND)
    lny = int(ny,PIO_OFFSET_KIND)

    n = 0
    do k = 1,nz
      do jsea = 1,nseal_cpl
        call init_get_isea(isea, jsea)
        ix = mapsf(isea,1)                           ! global ix
        iy = mapsf(isea,2)                           ! global iy
        n = n+1
        dof3d(n) = ((iy-1)*lnx + ix) + (k-1)*lnx*lny ! local index : global index
      end do
    end do

    call pio_initdecomp(wav_pio_subsystem, PIO_REAL, (/nx,ny,nz/), dof3d, iodesc)
    deallocate(dof3d)

  end subroutine wav_pio_initdecomp_3d

  !===============================================================================
  !> Handle errors
  !!
  !! @param[in]  ierr        the error code
  !! @param[in]  string      the error message
  !!
  !> @author Denise.Worthen@noaa.gov
  !> @date 08-02-2024
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
