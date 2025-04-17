  !========================================================================
  !> Realize the import and export fields.

  !> @details Called by NUOPC to realize the import and export fields
  !! for the wave model. After the wave model initializes, the global index
  !! for all sea points is retrieved using the WW3 mapsf array. A global index
  !! array is then constructed which contains both land and sea points, with
  !! the land points at the end of the array. An ESMF Distgrid object is created
  !! using this global index array. The distgrid is then transfered to the ESMF
  !! Mesh provided for the wave model domain. If the provided Mesh does not contain
  !! a grid mask, then the internal WW3 mask is transfered to the Mesh, otherwise
  !! the mask provided with the mesh file will be used. This mask is used by
  !! CMEPS to map to and from the wave model. Once the mesh has been created, the
  !! advertised fields are realized on the mesh.
  !!
  !! @param[in]    gcomp           an ESMF_GridComp object
  !! @param[in]    importState     an ESMF_State object for import fields
  !! @param[in]    exportState     an ESMF_State object for export fields
  !! @param[in]    clock           an ESMF_Clock object
  !! @param[out]   rc              return code
  !!
  !> @author mvertens@ucar.edu, Denise.Worthen@noaa.gov
  !> @date 01-05-2022
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

!    use w3odatmd        , only : w3nout, w3seto, naproc, naperr
!    use w3timemd        , only : stme21
!    use w3adatmd        , only : w3naux, w3seta
!    use w3idatmd        , only : w3seti, w3ninp
!    use w3gdatmd        , only : nk, nseal, nsea, nx, ny, mapsf, w3nmod, w3setg
!    use w3gdatmd        , only : rlgtype, ungtype, gtype
!    use w3wdatmd        , only : va, time, w3ndat, w3dimw, w3setw
!    use w3parall        , only : init_get_isea
!    use wav_shel_inp    , only : set_shel_io
!    use wav_history_mod , only : wav_history_init
!    use wav_pio_mod     , only : wav_pio_init
!    use wav_shr_mod     , only : diagnose_mesh, write_meshdecomp, wav_loginit
#ifdef W3_PDLIB
!    use yowNodepool     , only : ng
#endif

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! ! local variables
    ! type(ESMF_DistGrid)            :: distGrid
    ! type(ESMF_Mesh)                :: Emesh
    ! type(ESMF_Array)               :: elemMaskArray
    ! type(ESMF_VM)                  :: vm
    ! type(ESMF_Time)                :: esmfTime, startTime, currTime, stopTime
    ! type(ESMF_TimeInterval)        :: TimeOffset
    ! type(ESMF_TimeInterval)        :: TimeStep
    ! type(ESMF_Calendar)            :: calendar
    ! type(ESMF_Info)                :: info
    ! character(CL)                  :: cvalue
    ! integer                        :: shrlogunit
    ! integer                        :: yy,mm,dd,hh,ss
    ! integer                        :: start_ymd         ! start date (yyyymmdd)
    ! integer                        :: start_tod         ! start time of day (sec)
    ! integer                        :: stop_ymd          ! stop date (yyyymmdd)
    ! integer                        :: stop_tod          ! stop time of day (sec)
    ! integer                        :: ix, iy
    ! character(CL)                  :: starttype
    ! integer                        :: ntrace(2)
    ! integer                        :: n, jsea,isea, ncnt
    ! integer                        :: nlnd, nlnd_global, nlnd_local
    ! integer                        :: my_lnd_start, my_lnd_end
    ! integer, allocatable, target   :: mask_global(:)
    ! integer, allocatable, target   :: mask_local(:)
    ! integer, allocatable           :: gindex_lnd(:)
    ! integer, allocatable           :: gindex_sea(:)
    ! integer, allocatable           :: gindex(:)
    ! integer(i4)                    :: maskmin
    ! integer(i4), pointer           :: meshmask(:)
    ! character(23)                  :: dtme21
    ! integer                        :: iam, mpi_comm, num_threads
    ! character(ESMF_MAXSTR)         :: msgString
    ! character(ESMF_MAXSTR)         :: diro
    ! character(CL)                  :: logfile
    ! logical                        :: local
    ! integer                        :: imod, idsi, idso, idss, idst, idse
    ! integer                        :: mds(15) ! Note that nds is set to this in w3initmod
    ! integer                        :: stdout
    ! integer                        :: petcount
    ! real(r8)                       :: toff
    ! logical                        :: isPresent, isSet
    character(len=*), parameter    :: subname = '(wav_comp_nuopc:InitializeRealize)'
    ! -------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) call ESMF_LogWrite(trim(subname)//' called', ESMF_LOGMSG_INFO)

    call ufs_settimer(wtime)
!     !--------------------------------------------------------------------
!     ! Set up data structures
!     !--------------------------------------------------------------------

!     call w3nmod ( 1, 6, 6 )
!     call w3ndat (    6, 6 )
!     call w3naux (    6, 6 )
!     call w3nout (    6, 6 )
!     call w3ninp (    6, 6 )

!     call w3setg ( 1, 6, 6 )
!     call w3setw ( 1, 6, 6 )
!     call w3seta ( 1, 6, 6 )
!     call w3seto ( 1, 6, 6 )
!     call w3seti ( 1, 6, 6 )

!     !----------------------------------------------------------------------------
!     ! Generate local mpi comm
!     !----------------------------------------------------------------------------

!     call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
!     if (ChkErr(rc,__LINE__,u_FILE_u)) return

!     call ESMF_VMGet(vm, mpiCommunicator=mpi_comm, peCount=petcount, localPet=iam, rc=rc)
!     if (ChkErr(rc,__LINE__,u_FILE_u)) return

!     call ESMF_InfoGetFromHost(gcomp, info=info, rc=rc)
!     if (ChkErr(rc,__LINE__,u_FILE_u)) return
!     call ESMF_InfoGet(info, key="/NUOPC/Hint/PePerPet/MaxCount", value=num_threads, default=1, rc=rc)
!     if (ChkErr(rc,__LINE__,u_FILE_u)) return

!     naproc = petcount
!     iaproc = iam + 1
!     napout = 1
!     naperr = 1
!     if (iaproc == napout) root_task = .true.

!     !--------------------------------------------------------------------
!     ! IO set-up
!     !--------------------------------------------------------------------

!     if (cesmcoupled) then
!       shrlogunit = 6
!       if ( root_task ) then
!         call NUOPC_CompAttributeGet(gcomp, name="diro", value=diro, rc=rc)
!         if (chkerr(rc,__LINE__,u_FILE_u)) return
!         call NUOPC_CompAttributeGet(gcomp, name="logfile", value=logfile, rc=rc)
!         if (chkerr(rc,__LINE__,u_FILE_u)) return
!         open (newunit=stdout, file=trim(diro)//"/"//trim(logfile))
!         logfile_is_assigned = .true.
!       else
!         stdout = 6
!       endif
!     else
!       if ( root_task ) then
!         open (newunit=stdout, file='log.ww3')
!         logfile_is_assigned = .true.
!       else
!         stdout = 6
!       end if
!     end if

!     call set_shel_io(stdout,mds,ntrace)

!     if ( root_task ) then
!       write(stdout,'(a)')'      *** WAVEWATCH III Program shell ***      '
!       write(stdout,'(a)')'==============================================='
!       write(stdout,'(/)')
!       write(stdout,'(a,l)')' Wave wav_coupling_to_cice setting is ',wav_coupling_to_cice
!     end if

!     !--------------------------------------------------------------------
!     ! Initialize run type
!     !--------------------------------------------------------------------

!     call NUOPC_CompAttributeGet(gcomp, name='start_type', value=starttype, rc=rc)
!     if (ChkErr(rc,__LINE__,u_FILE_u)) return
!     if ( trim(starttype) == trim('startup')) then
!       runtype = "initial"
!     else if (trim(starttype) == trim('continue') ) then
!       runtype = "continue"
!     else if (trim(starttype) == trim('branch')) then
!       runtype = "branch"
!     end if
!     if ( root_task ) then
!       write(stdout,'(a)') ' WW3 runtype is '//trim(runtype)
!     end if
!     call ESMF_LogWrite('WW3 runtype is '//trim(runtype), ESMF_LOGMSG_INFO)

!     !--------------------------------------------------------------------
!     ! Time initialization
!     !--------------------------------------------------------------------

!     ! TIME0 = from ESMF clock
!     ! NOTE - are not setting TIMEN here

!     if ( root_task ) then
!       write(stdout,'(a)')'  Time interval : '
!       write(stdout,'(a)')'--------------------------------------------------'
!     end if

!     call ESMF_ClockPrint(clock, options="startTime", preString="Model Start Time: ", &
!          unit=msgString, rc=rc)
!     call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
!     call ESMF_ClockPrint(clock, options="currTime", preString="Model Current Time: ", &
!          unit=msgString, rc=rc)
!     call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
!     call ESMF_ClockGet( clock, startTime=startTime, currTime=currTime, rc=rc)
!     TimeOffset = currTime - startTime
!     call ESMF_TimeIntervalGet(TimeOffset, h_r8=toff, rc=rc)
!     write(msgstring,'(a,g14.7)')'TimeOffset: CurrTime - StartTime = ',toff
!     call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
!     ! Initial run or restart run
!     if ( runtype == "initial") then
!       call ESMF_ClockGet( clock, startTime=esmfTime, rc=rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
! #ifndef W3_CESMCOUPLED
!       esmfTime = esmfTime + TimeOffset
! #endif
!     else
!       call ESMF_ClockGet( clock, currTime=esmfTime, rc=rc )
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!     endif
!     ! Determine time attributes for history output
!     call ESMF_TimeGet( startTime, timeString=time_origin, calendar=calendar, rc=rc )
!     if (ChkErr(rc,__LINE__,u_FILE_u)) return
!     time_origin = 'seconds since '//time_origin(1:10)//' '//time_origin(12:19)
!     !call ESMF_ClockGet(clock, calendar=calendar)
!     if (calendar == ESMF_CALKIND_GREGORIAN) then
!       calendar_name = 'standard'
!     else if (calendar == ESMF_CALKIND_NOLEAP) then
!       calendar_name = 'noleap'
!     end if
!     call ESMF_TimeGet( esmfTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
!     if (ChkErr(rc,__LINE__,u_FILE_u)) return
!     call ymd2date(yy, mm, dd, start_ymd)

!     hh = start_tod/3600
!     mm = (start_tod - (hh * 3600))/60
!     ss = start_tod - (hh*3600) - (mm*60)

!     time0(1) = start_ymd
!     time0(2) = hh*10000 + mm*100 + ss

!     call ESMF_ClockGet( clock, stopTime=stopTime, rc=rc)
!     if (ChkErr(rc,__LINE__,u_FILE_u)) return
!     call ESMF_TimeGet( stopTime, yy=yy, mm=mm, dd=dd, s=stop_tod, rc=rc )
!     if (ChkErr(rc,__LINE__,u_FILE_u)) return
!     call ymd2date(yy, mm, dd, stop_ymd)

!     hh = stop_tod/3600
!     mm = (stop_tod - (hh * 3600))/60
!     ss = stop_tod - (hh*3600) - (mm*60)

!     timen(1) = stop_ymd
!     timen(2) = hh*10000 + mm*100 + ss

!     call stme21 ( time0 , dtme21 )
!     if ( root_task ) then
!       write (stdout,'(a)')' Starting time : '//trim(dtme21)
!       write (stdout,'(a,i8,2x,i8)') ' start_ymd, stop_ymd = ',start_ymd, stop_ymd
!     end if

!     !--------------------------------------------------------------------
!     ! Initialize PIO. This needs to be done prior to the call to w3init
!     ! in order to read the restart file. The filename strings must also
!     ! be available
!     !--------------------------------------------------------------------

!     if (cesmcoupled) then
!       if (len_trim(inst_suffix) > 0) then
!         user_restfname = trim(casename)//'.ww3'//trim(inst_suffix)//'.r.'
!         user_histfname = trim(casename)//'.ww3'//trim(inst_suffix)//'.hi.'
!       else
!         user_restfname = trim(casename)//'.ww3.r.'
!         user_histfname = trim(casename)//'.ww3.hi.'
!       endif

!       ! netcdf is used for CESM history and restart
!       use_historync = .true.
!       use_restartnc = .true.
!     else
!       call NUOPC_CompAttributeGet(gcomp, name='use_restartnc', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!       if (isPresent .and. isSet) then
!         use_restartnc=(trim(cvalue)=="true")
!       end if
!       if (root_task) write(stdout,'(a,l4)') trim(subname)//': Wave use_restartnc setting is ',use_restartnc

!       ! user filenaming is required with netcdf restarts or restart_from_binary. If netcdf restarts are not used,
!       ! only native WW3 file naming is possible
!       if (use_restartnc) then
!         user_restfname = trim(casename)//'.ww3.r.'
!         if (root_task) write(stdout,'(a)') trim(subname)//': Custom restart prefix is '//trim(user_restfname)
!       end if

!       call NUOPC_CompAttributeGet(gcomp, name='use_historync', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!       if (isPresent .and. isSet) then
!         use_historync=(trim(cvalue)=="true")
!       end if
!       if (root_task) write(stdout,'(a,l4)') trim(subname)//': Wave use_historync setting is ',use_historync

!       ! user filenaming is optional with netcdf output. If netcdf history is not used, only native WW3
!       ! naming is possible
!       if (use_historync) then
!         call NUOPC_CompAttributeGet(gcomp, name='user_histname', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
!         if (ChkErr(rc,__LINE__,u_FILE_u)) return
!         if (trim(cvalue)=="true") then
!           user_histfname = trim(casename)//'.ww3.hi.'
!           if (root_task) write(stdout,'(a)') trim(subname)//': Custom history prefix is '//trim(user_histfname)
!         else
!           user_histfname = ''
!         end if
!       end if
!     end if ! if (cesmcoupled)

!     ! allow startup from binary restarts as special case
!     if (use_restartnc) then
!       call NUOPC_CompAttributeGet(gcomp, name='restart_from_binary', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!       if (isPresent .and. isSet) then
!         restart_from_binary=(trim(cvalue)=="true")
!       end if
!       if (root_task) write(stdout,'(a,l4)') trim(subname)//': Wave restart_from_binary setting is ',restart_from_binary
!     end if

!     if (use_restartnc .or. use_historync) then
!       call wav_pio_init(gcomp, mpi_comm, stdout, naproc/num_threads, rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!     end if

!     !--------------------------------------------------------------------
!     ! Wave model initialization
!     !--------------------------------------------------------------------

! #ifndef W3_CESMCOUPLED
!     call waveinit_ufs(gcomp, stdout, ntrace, mpi_comm, mds, rc)
!     if (ChkErr(rc,__LINE__,u_FILE_u)) return
! #else
!     time = time0
!     call ESMF_ClockGet( clock, timeStep=timeStep, rc=rc)
!     if (ChkErr(rc,__LINE__,u_FILE_u)) return
!     call waveinit_cesm(gcomp, ntrace, mpi_comm, mds, rc)
!     if (ChkErr(rc,__LINE__,u_FILE_u)) return
! #endif
!     !call mpi_barrier ( mpi_comm, ierr )
!     if ( root_task ) then
!       inquire(unit=stdout, name=logfile)
!       write(*,'(a)')'WW3 log written to '//trim(logfile)
!     end if

!     if (wav_coupling_to_cice) then
!       if (nwav_elev_spectrum .gt. nk) then
!         call ESMF_LogWrite('nwav_elev_spectrum is greater than nk ', ESMF_LOGMSG_INFO)
!         call ESMF_Finalize(endflag=ESMF_END_ABORT)
!       end if
!     end if

    !--------------------------------------------------------------------
    ! Mesh initialization
    !--------------------------------------------------------------------

    if (gtype .eq. ungtype) then
      unstr_mesh = .true.
    else
      unstr_mesh = .false.
    end if

    ! Create a  global index array for sea points.
    !
    ! Note that nsea is the global number of sea points - and nseal is the local
    ! number of sea points. For the unstr mesh, the nsea points are on mesh nodes.
    ! We will use the gindex to set the element distgrid of a dual mesh. A dual mesh
    ! contains the mesh nodes at the center of each element. For the domain decomposition
    ! case (PDLIB), set a value of the local sea points on this processor minus the
    ! ghost points.
#ifdef W3_PDLIB
    nseal_cpl = nseal - ng
#else
    nseal_cpl = nseal
#endif
    allocate(gindex_sea(nseal_cpl))
    do jsea=1, nseal_cpl
      call init_get_isea(isea, jsea)
      ix = mapsf(isea,1)
      iy = mapsf(isea,2)
      gindex_sea(jsea) = ix + (iy-1)*nx
    end do

    if (unstr_mesh) then
      ! create distGrid from global index array of sea points with no ghost points
      DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex_sea, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      ! create a global index array for non-sea (i.e. land points)
      allocate(mask_global(nx*ny), mask_local(nx*ny))
      mask_local(:) = 0
      mask_global(:) = 0
      do jsea=1, nseal_cpl
        call init_get_isea(isea, jsea)
        ix = mapsf(isea,1)
        iy = mapsf(isea,2)
        mask_local(ix + (iy-1)*nx) = 1
      end do
      call ESMF_VMAllReduce(vm, sendData=mask_local, recvData=mask_global, count=nx*ny, &
           reduceflag=ESMF_REDUCE_MAX, rc=rc)

      nlnd_global = nx*ny - nsea
      nlnd_local = nlnd_global / naproc
      my_lnd_start = nlnd_local*iam + min(iam, mod(nlnd_global, naproc)) + 1
      if (iam < mod(nlnd_global, naproc)) then
        nlnd_local = nlnd_local + 1
      end if
      my_lnd_end = my_lnd_start + nlnd_local - 1

      allocate(gindex_lnd(my_lnd_end - my_lnd_start + 1))
      ncnt = 0
      do n = 1,nx*ny
        if (mask_global(n) == 0) then ! this is a land point
          ncnt = ncnt + 1
          if (ncnt >= my_lnd_start .and. ncnt <= my_lnd_end) then
            gindex_lnd(ncnt - my_lnd_start + 1) = n
          end if
        end if
      end do
      deallocate(mask_global)
      deallocate(mask_local)

      ! create a global index that includes both sea and land - but put land at the end
      nlnd = (my_lnd_end - my_lnd_start + 1)
      allocate(gindex(nlnd + nseal_cpl))
      do ncnt = 1,nlnd + nseal
        if (ncnt <= nseal_cpl) then
          gindex(ncnt) = gindex_sea(ncnt)
        else
          gindex(ncnt) = gindex_lnd(ncnt-nseal_cpl)
        end if
      end do

      ! create distGrid from global index array
      DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! get the mesh file name
    call NUOPC_CompAttributeGet(gcomp, name='mesh_wav', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! read in the mesh with the the DistGrid
    EMesh = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, &
         elementDistgrid=Distgrid,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      if (unstr_mesh) then
        call diagnose_mesh(EMesh, size(gindex_sea), 'EMesh', rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        deallocate(gindex_sea)
      else
        call diagnose_mesh(EMesh, size(gindex), 'EMesh', rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        deallocate(gindex)
        deallocate(gindex_sea)
        deallocate(gindex_lnd)
      end if
    end if

    if (.not. unstr_mesh) then
      ! obtain the mesh mask and find the minimum value across all PEs
      call ESMF_MeshGet(EMesh, elementDistgrid=Distgrid, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_DistGridGet(Distgrid, localDe=0, elementCount=ncnt, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      allocate(meshmask(ncnt))
      elemMaskArray = ESMF_ArrayCreate(Distgrid, farrayPtr=meshmask, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_MeshGet(Emesh, elemMaskArray=elemMaskArray, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_VMAllFullReduce(vm, sendData=meshmask, recvData=maskmin, count=ncnt, &
           reduceflag=ESMF_REDUCE_MIN, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      if (maskmin == 1) then
        ! replace mesh mask with internal mask
        meshmask(:) = 0
        meshmask(1:nseal_cpl) = 1
        call ESMF_MeshSet(mesh=EMesh, elementMask=meshmask, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      end if

      if (dbug_flag > 5) then
        call ESMF_ArrayWrite(elemMaskArray, 'meshmask.nc', variableName = 'mask', &
             overwrite=.true., rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
      end if
      deallocate(meshmask)
    end if

    if (dbug_flag > 5) then
      call write_meshdecomp(Emesh, 'emesh', rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !--------------------------------------------------------------------
    ! Realize the actively coupled fields
    !--------------------------------------------------------------------
    call realize_fields(gcomp, mesh=Emesh, flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! !--------------------------------------------------------------------
    ! ! Intialize the list of requested output variables for netCDF output.
    ! ! This needs to occur after mod_def has been read in w3init since
    ! ! some variables are available only if they are defined in the mod_def
    ! !--------------------------------------------------------------------

    ! if (use_historync) then
    !   call wav_history_init(stdout)
    ! end if

    !--------------------------------------------------------------------
    ! Write the header string for WW3 native logging
    !--------------------------------------------------------------------

    if (root_task) then
      if (verboselog) call wav_loginit(stdout)
    end if

    if (root_task) call ufs_logtimer(nu_timer,time,start_tod,'InitializeRealize time: ',runtimelog,wtime)

    if (dbug_flag > 5) call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeRealize
