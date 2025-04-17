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

    !use wav_shr_flags, only : w3_pdlib_flag

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
    ! retrieve configuration settings
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

    if (cesmcoupled) then
      if (len_trim(inst_suffix) > 0) then
        user_restfname = trim(casename)//'.ww3'//trim(inst_suffix)//'.r.'
        user_histfname = trim(casename)//'.ww3'//trim(inst_suffix)//'.hi.'
      else
        user_restfname = trim(casename)//'.ww3.r.'
        user_histfname = trim(casename)//'.ww3.hi.'
      endif

      ! netcdf is used for CESM history and restart
      use_historync = .true.
      use_restartnc = .true.
    else
      call NUOPC_CompAttributeGet(gcomp, name='use_restartnc', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (isPresent .and. isSet) then
        use_restartnc=(trim(cvalue)=="true")
      end if
      if (root_task) write(stdout,'(a,l4)') trim(subname)//': Wave use_restartnc setting is ',use_restartnc

      ! user filenaming is required with netcdf restarts or restart_from_binary. If netcdf restarts are not used,
      ! only native WW3 file naming is possible
      if (use_restartnc) then
        user_restfname = trim(casename)//'.ww3.r.'
        if (root_task) write(stdout,'(a)') trim(subname)//': Custom restart prefix is '//trim(user_restfname)
      end if

      call NUOPC_CompAttributeGet(gcomp, name='use_historync', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (isPresent .and. isSet) then
        use_historync=(trim(cvalue)=="true")
      end if
      if (root_task) write(stdout,'(a,l4)') trim(subname)//': Wave use_historync setting is ',use_historync

      ! user filenaming is optional with netcdf output. If netcdf history is not used, only native WW3
      ! naming is possible
      if (use_historync) then
        call NUOPC_CompAttributeGet(gcomp, name='user_histname', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        if (trim(cvalue)=="true") then
          user_histfname = trim(casename)//'.ww3.hi.'
          if (root_task) write(stdout,'(a)') trim(subname)//': Custom history prefix is '//trim(user_histfname)
        else
          user_histfname = ''
        end if
      end if
    end if ! if (cesmcoupled)

    ! allow startup from binary restarts as special case
    if (use_restartnc) then
      call NUOPC_CompAttributeGet(gcomp, name='restart_from_binary', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (isPresent .and. isSet) then
        restart_from_binary=(trim(cvalue)=="true")
      end if
      if (root_task) write(stdout,'(a,l4)') trim(subname)//': Wave restart_from_binary setting is ',restart_from_binary
    end if

    !XXXX
    !XXXX

    !--------------------------------------------------------------------
    ! Set up data structures
    !--------------------------------------------------------------------

    call w3nmod ( 1, 6, 6 )
    call w3ndat (    6, 6 )
    call w3naux (    6, 6 )
    call w3nout (    6, 6 )
    call w3ninp (    6, 6 )

    call w3setg ( 1, 6, 6 )
    call w3setw ( 1, 6, 6 )
    call w3seta ( 1, 6, 6 )
    call w3seto ( 1, 6, 6 )
    call w3seti ( 1, 6, 6 )

    !----------------------------------------------------------------------------
    ! Generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=mpi_comm, peCount=petcount, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_InfoGetFromHost(gcomp, info=info, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_InfoGet(info, key="/NUOPC/Hint/PePerPet/MaxCount", value=num_threads, default=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    naproc = petcount
    iaproc = iam + 1
    napout = 1
    naperr = 1
    if (iaproc == napout) root_task = .true.

    !--------------------------------------------------------------------
    ! IO set-up
    !--------------------------------------------------------------------

    if (cesmcoupled) then
      shrlogunit = 6
      if ( root_task ) then
        call NUOPC_CompAttributeGet(gcomp, name="diro", value=diro, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeGet(gcomp, name="logfile", value=logfile, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        open (newunit=stdout, file=trim(diro)//"/"//trim(logfile))
        logfile_is_assigned = .true.
      else
        stdout = 6
      endif
    else
      if ( root_task ) then
        open (newunit=stdout, file='log.ww3')
        logfile_is_assigned = .true.
      else
        stdout = 6
      end if
    end if

    call set_shel_io(stdout,mds,ntrace)

    if ( root_task ) then
      write(stdout,'(a)')'      *** WAVEWATCH III Program shell ***      '
      write(stdout,'(a)')'==============================================='
      write(stdout,'(/)')
      write(stdout,'(a,l)')' Wave wav_coupling_to_cice setting is ',wav_coupling_to_cice
    end if

    !--------------------------------------------------------------------
    ! Initialize run type
    !--------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=starttype, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if ( trim(starttype) == trim('startup')) then
      runtype = "initial"
    else if (trim(starttype) == trim('continue') ) then
      runtype = "continue"
    else if (trim(starttype) == trim('branch')) then
      runtype = "branch"
    end if
    if ( root_task ) then
      write(stdout,'(a)') ' WW3 runtype is '//trim(runtype)
    end if
    call ESMF_LogWrite('WW3 runtype is '//trim(runtype), ESMF_LOGMSG_INFO)

    !--------------------------------------------------------------------
    ! Time initialization
    !--------------------------------------------------------------------

    ! TIME0 = from ESMF clock
    ! NOTE - are not setting TIMEN here

    if ( root_task ) then
      write(stdout,'(a)')'  Time interval : '
      write(stdout,'(a)')'--------------------------------------------------'
    end if

    call ESMF_ClockPrint(clock, options="startTime", preString="Model Start Time: ", &
         unit=msgString, rc=rc)
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    call ESMF_ClockPrint(clock, options="currTime", preString="Model Current Time: ", &
         unit=msgString, rc=rc)
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    call ESMF_ClockGet( clock, startTime=startTime, currTime=currTime, rc=rc)
    TimeOffset = currTime - startTime
    call ESMF_TimeIntervalGet(TimeOffset, h_r8=toff, rc=rc)
    write(msgstring,'(a,g14.7)')'TimeOffset: CurrTime - StartTime = ',toff
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    ! Initial run or restart run
    if ( runtype == "initial") then
      call ESMF_ClockGet( clock, startTime=esmfTime, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
#ifndef W3_CESMCOUPLED
      esmfTime = esmfTime + TimeOffset
#endif
    else
      call ESMF_ClockGet( clock, currTime=esmfTime, rc=rc )
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
    ! Determine time attributes for history output
    call ESMF_TimeGet( startTime, timeString=time_origin, calendar=calendar, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    time_origin = 'seconds since '//time_origin(1:10)//' '//time_origin(12:19)
    !call ESMF_ClockGet(clock, calendar=calendar)
    if (calendar == ESMF_CALKIND_GREGORIAN) then
      calendar_name = 'standard'
    else if (calendar == ESMF_CALKIND_NOLEAP) then
      calendar_name = 'noleap'
    end if
    call ESMF_TimeGet( esmfTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ymd2date(yy, mm, dd, start_ymd)

    hh = start_tod/3600
    mm = (start_tod - (hh * 3600))/60
    ss = start_tod - (hh*3600) - (mm*60)

    time0(1) = start_ymd
    time0(2) = hh*10000 + mm*100 + ss

    call ESMF_ClockGet( clock, stopTime=stopTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( stopTime, yy=yy, mm=mm, dd=dd, s=stop_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ymd2date(yy, mm, dd, stop_ymd)

    hh = stop_tod/3600
    mm = (stop_tod - (hh * 3600))/60
    ss = stop_tod - (hh*3600) - (mm*60)

    timen(1) = stop_ymd
    timen(2) = hh*10000 + mm*100 + ss

    call stme21 ( time0 , dtme21 )
    if ( root_task ) then
      write (stdout,'(a)')' Starting time : '//trim(dtme21)
      write (stdout,'(a,i8,2x,i8)') ' start_ymd, stop_ymd = ',start_ymd, stop_ymd
    end if

    ! !--------------------------------------------------------------------
    ! ! Initialize PIO. This needs to be done prior to the call to w3init
    ! ! in order to read the restart file. The filename strings must also
    ! ! be available
    ! !--------------------------------------------------------------------

    ! if (cesmcoupled) then
    !   if (len_trim(inst_suffix) > 0) then
    !     user_restfname = trim(casename)//'.ww3'//trim(inst_suffix)//'.r.'
    !     user_histfname = trim(casename)//'.ww3'//trim(inst_suffix)//'.hi.'
    !   else
    !     user_restfname = trim(casename)//'.ww3.r.'
    !     user_histfname = trim(casename)//'.ww3.hi.'
    !   endif

    !   ! netcdf is used for CESM history and restart
    !   use_historync = .true.
    !   use_restartnc = .true.
    ! else
    !   call NUOPC_CompAttributeGet(gcomp, name='use_restartnc', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !   if (isPresent .and. isSet) then
    !     use_restartnc=(trim(cvalue)=="true")
    !   end if
    !   if (root_task) write(stdout,'(a,l4)') trim(subname)//': Wave use_restartnc setting is ',use_restartnc

    !   ! user filenaming is required with netcdf restarts or restart_from_binary. If netcdf restarts are not used,
    !   ! only native WW3 file naming is possible
    !   if (use_restartnc) then
    !     user_restfname = trim(casename)//'.ww3.r.'
    !     if (root_task) write(stdout,'(a)') trim(subname)//': Custom restart prefix is '//trim(user_restfname)
    !   end if

    !   call NUOPC_CompAttributeGet(gcomp, name='use_historync', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !   if (isPresent .and. isSet) then
    !     use_historync=(trim(cvalue)=="true")
    !   end if
    !   if (root_task) write(stdout,'(a,l4)') trim(subname)//': Wave use_historync setting is ',use_historync

    !   ! user filenaming is optional with netcdf output. If netcdf history is not used, only native WW3
    !   ! naming is possible
    !   if (use_historync) then
    !     call NUOPC_CompAttributeGet(gcomp, name='user_histname', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    !     if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !     if (trim(cvalue)=="true") then
    !       user_histfname = trim(casename)//'.ww3.hi.'
    !       if (root_task) write(stdout,'(a)') trim(subname)//': Custom history prefix is '//trim(user_histfname)
    !     else
    !       user_histfname = ''
    !     end if
    !   end if
    ! end if ! if (cesmcoupled)

    ! ! allow startup from binary restarts as special case
    ! if (use_restartnc) then
    !   call NUOPC_CompAttributeGet(gcomp, name='restart_from_binary', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !   if (isPresent .and. isSet) then
    !     restart_from_binary=(trim(cvalue)=="true")
    !   end if
    !   if (root_task) write(stdout,'(a,l4)') trim(subname)//': Wave restart_from_binary setting is ',restart_from_binary
    ! end if

     if (use_restartnc .or. use_historync) then
       call wav_pio_init(gcomp, mpi_comm, stdout, naproc/num_threads, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
     end if

    !--------------------------------------------------------------------
    ! Wave model initialization
    !--------------------------------------------------------------------

#ifndef W3_CESMCOUPLED
    call waveinit_ufs(gcomp, stdout, ntrace, mpi_comm, mds, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
#else
    !time = time0
    call ESMF_ClockGet( clock, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call waveinit_cesm(gcomp, ntrace, mpi_comm, mds, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
#endif
    !call mpi_barrier ( mpi_comm, ierr )
    if ( root_task ) then
      inquire(unit=stdout, name=logfile)
      write(*,'(a)')'WW3 log written to '//trim(logfile)
    end if

    if (wav_coupling_to_cice) then
      if (nwav_elev_spectrum .gt. nk) then
        call ESMF_LogWrite('nwav_elev_spectrum is greater than nk ', ESMF_LOGMSG_INFO)
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
    end if

    !--------------------------------------------------------------------
    ! Intialize the list of requested output variables for netCDF output.
    ! This needs to occur after mod_def has been read in w3init since
    ! some variables are available only if they are defined in the mod_def
    !--------------------------------------------------------------------

    if (use_historync) then
      call wav_history_init(stdout)
    end if

    call advertise_fields(importState, exportState, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeAdvertise
