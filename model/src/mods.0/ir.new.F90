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

    use w3odatmd        , only : w3seto, naproc
    use w3adatmd        , only : w3seta
    use w3idatmd        , only : w3seti
    use w3gdatmd        , only : w3setg
    use w3wdatmd        , only : va, time, w3setw
    use wav_history_mod , only : wav_history_init
    use wav_shr_mod     , only : wav_loginit

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    integer                        :: start_tod         ! start time of day (sec)
    logical                        :: local
    integer                        :: imod
    logical                        :: isPresent, isSet
    character(ESMF_MAXSTR)         :: preamb = './'
    character(ESMF_MAXSTR)         :: ifname = 'ww3_multi.inp'
    character(len=*), parameter    :: subname = '(wav_comp_nuopc:InitializeRealize)'
    ! -------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) call ESMF_LogWrite(trim(subname)//' called', ESMF_LOGMSG_INFO)

    call ufs_settimer(wtime)

    !--------------------------------------------------------------------
    ! Realize the actively coupled fields
    !--------------------------------------------------------------------
    call realize_fields(gcomp, mesh=Emesh, flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

#ifndef W3_CESMCOUPLED
    !TODO: when is this required?
    if (multigrid) then
      do imod = 1,nrgrd
        call w3setg ( imod, mdse, mdst )
        call w3setw ( imod, mdse, mdst )
        call w3seta ( imod, mdse, mdst )
        call w3seti ( imod, mdse, mdst )
        call w3seto ( imod, mdse, mdst )
        call wmsetm ( imod, mdse, mdst )
        local = iaproc .gt. 0 .and. iaproc .le. naproc
        if ( local .and. flcold .and. fliwnd ) call w3uini( va )
      enddo
    end if
#endif
    !--------------------------------------------------------------------
    ! Intialize the list of requested output variables for netCDF output.
    ! This needs to occur after mod_def has been read in w3init since
    ! some variables are available only if they are defined in the mod_def
    !--------------------------------------------------------------------

    if (use_historync) then
      call wav_history_init(stdout)
    end if

    !--------------------------------------------------------------------
    ! Write the header string for WW3 native logging
    !--------------------------------------------------------------------

    if (root_task) then
      if (verboselog) call wav_loginit(stdout)
    end if

    if (root_task) call ufs_logtimer(nu_timer,time,start_tod,'InitializeRealize time: ',runtimelog,wtime)

    if (dbug_flag > 5) call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeRealize
