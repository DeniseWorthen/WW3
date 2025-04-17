 !===============================================================================
  !> Read configuration attributes and advertise the import/export fields

  !> @details Called by NUOPC to read configuration attributes and to advertise the
  !! import and export fields. The configuration attributes are used to control run
  !! time settings, such as ESMF memory profiling, additional debug logging
  !! and character strings for specific use cases. A set of configuration attributes
  !! is also read to describe any scalar fields to be added to a state. For coupling
  !! with the wave model, only a scalar field for the dimensions of the wave model
  !! is required. The scalar field is added to the export state to communicate to the
  !! CMEPS mediator the domain dimensions of the wave model in order to write
  !! mediator history and restart files. The attribute ScalarFieldName sets the name
  !! of the scalar field in the export state, the ScalarFieldCount sets the
  !! dimensionality of the scalar field and the ScalarFieldIdxGridNX (NY) set the
  !! index of the NX or NY dimension in the scalar field.
  !!
  !! @param[in]    gcomp             an ESMF_GridComp object
  !! @param[in]    importState       an ESMF_State object for import fields
  !! @param[in]    exportState       an ESMF_State object for export fields
  !! @param[in]    clock             an ESMF_Clock object
  !! @param[out]   rc                return code
  !!
  !> @author mvertens@ucar.edu, Denise.Worthen@noaa.gov
  !> @date 01-05-2022
  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    use wav_shr_flags, only : w3_pdlib_flag

    ! input/output arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=CL) :: logmsg
    logical           :: isPresent, isSet
    character(len=CL) :: cvalue
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    call ufs_settimer(wtime)
    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//' called', ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! advertise fields
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
      flds_scalar_name = trim(cvalue)
      call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogWrite(trim(subname)//'Need to set attribute ScalarFieldName',&
           ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
      rc = ESMF_FAILURE
      return
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
      read(cvalue, *) flds_scalar_num
      write(logmsg,*) flds_scalar_num
      call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogWrite(trim(subname)//'Need to set attribute ScalarFieldCount',&
           ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
      rc = ESMF_FAILURE
      return
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
      read(cvalue,*) flds_scalar_index_nx
      write(logmsg,*) flds_scalar_index_nx
      call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogWrite(trim(subname)//'Need to set attribute ScalarFieldIdxGridNX',&
           ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
      rc = ESMF_FAILURE
      return
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
      read(cvalue,*) flds_scalar_index_ny
      write(logmsg,*) flds_scalar_index_ny
      call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogWrite(trim(subname)//'Need to set attribute ScalarFieldIdxGridNY',&
           ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
      rc = ESMF_FAILURE
      return
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ProfileMemory", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
      read(cvalue,*) profile_memory
      call ESMF_LogWrite(trim(subname)//': profile_memory = '//trim(cvalue), ESMF_LOGMSG_INFO)
    end if

    call NUOPC_CompAttributeGet(gcomp, name="merge_import", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
      if (trim(cvalue) == '.true.') then
        merge_import = .true.
      end if
    end if
    if (merge_import) then
      if (w3_pdlib_flag) then
        call ESMF_LogWrite('Merge_import is not valid with PDLIB', ESMF_LOGMSG_INFO)
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
    end if

    call NUOPC_CompAttributeGet(gcomp, name='dbug_flag', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
      read(cvalue,*) dbug_flag
    end if
    write(logmsg,'(A,i6)') trim(subname)//': Wave cap dbug_flag is ',dbug_flag
    call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_INFO)

    ! Get casename
    call NUOPC_CompAttributeGet(gcomp, name="case_name", value=casename, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(logmsg,'(A)') trim(subname)//': Wave casename setting : '//trim(casename)
    call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_INFO)

    ! Get component instance
    call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
      call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", value=inst_suffix, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      cvalue = inst_suffix(2:)
      read(cvalue, *) inst_index
    else
      inst_suffix = ""
      inst_index=1
    endif

    ! Determine wave-ice coupling
    wav_coupling_to_cice = .false.
    call NUOPC_CompAttributeGet(gcomp, name='wav_coupling_to_cice', value=cvalue, isPresent=isPresent, &
         isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) wav_coupling_to_cice
    end if
    write(logmsg,'(A,l)') trim(subname)//': Wave wav_coupling_to_cice setting is ',wav_coupling_to_cice
    call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_INFO)

    ! Determine Runtime logging
    call NUOPC_CompAttributeGet(gcomp, name="RunTimeLog", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) runtimelog=(trim(cvalue)=="true")
    write(logmsg,*) runtimelog
    call ESMF_LogWrite('WW3_cap:RunTimeLog = '//trim(logmsg), ESMF_LOGMSG_INFO)
    if (runtimelog) then
      call ufs_file_setLogUnit('./log.ww3.timer',nu_timer,runtimelog)
    end if

    ! Determine verbose native WW3 logging
    call NUOPC_CompAttributeGet(gcomp, name="verboselog", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) verboselog=(trim(cvalue)=="true")
    write(logmsg,*) verboselog
    call ESMF_LogWrite('WW3_cap: Verbose WW3 native logging is = '//trim(logmsg), ESMF_LOGMSG_INFO)

    call advertise_fields(importState, exportState, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeAdvertise
