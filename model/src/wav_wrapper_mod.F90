module wav_wrapper_mod

#ifdef CESMCOUPLED
  use perf_mod      , only : t_startf, t_stopf, t_barrierf
  use shr_file_mod  , only : shr_file_getlogunit, shr_file_setlogunit
  use wav_kind_mod  , only : r8 => shr_kind_r8, r4 => shr_kind_r4, i4 => shr_kind_i4
  use wav_kind_mod  , only : CL => shr_kind_cl, CS => shr_kind_cs

  implicit none

  real(r8) :: wtime = 0.0
contains
  ! Define stub routines that do nothing - they are just here to avoid
  ! having cppdefs in the main program
  subroutine ufs_settimer(timevalue)
    real(r8),    intent(inout) :: timevalue
  end subroutine ufs_settimer
  subroutine ufs_logtimer(nunit,elapsedsecs,tod,string,runtimelog,wtime0)
    integer,          intent(in) :: nunit
    integer(i4),      intent(in) :: elapsedsecs(2), tod
    character(len=*), intent(in) :: string
    logical,          intent(in) :: runtimelog
    real(r8),         intent(in) :: wtime0
  end subroutine ufs_logtimer
  subroutine ufs_file_setLogUnit(filename,nunit,runtimelog)
    character(len=*),  intent(in)  :: filename
    logical,           intent(in)  :: runtimelog
    integer,           intent(out) :: nunit
  end subroutine ufs_file_setLogUnit
  subroutine ufs_logfhour(msg,hour)
    character(len=*),  intent(in)  :: msg
    real(r8),          intent(in)  :: hour
  end subroutine ufs_logfhour
#else

  use wav_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4, i4 => shr_kind_i4
  use wav_kind_mod , only : CL => shr_kind_cl, CS => shr_kind_cs

  implicit none

  real(r8) :: wtime = 0.0
contains
  subroutine ufs_settimer(timevalue)
    real(r8),    intent(inout) :: timevalue
    real(r8)                   :: MPI_Wtime
    timevalue = MPI_Wtime()
  end subroutine ufs_settimer

  subroutine ufs_logtimer(nunit,elapsedsecs,tod,string,runtimelog,wtime0)
    integer,          intent(in)    :: nunit
    integer(i4),      intent(in)    :: elapsedsecs(2),tod
    character(len=*), intent(in)    :: string
    logical,          intent(in)    :: runtimelog
    real(r8),         intent(in)    :: wtime0
    real(r8)                        :: MPI_Wtime, timevalue
    if (.not. runtimelog) return
    if (wtime0 > 0.) then
      timevalue = MPI_Wtime()-wtime0
      write(nunit,'(3i8,a,g14.7)')elapsedsecs,tod,' WW3 '//trim(string),timevalue
    end if
  end subroutine ufs_logtimer

  subroutine ufs_file_setLogUnit(filename,nunit,runtimelog)
    character(len=*),  intent(in)    :: filename
    logical,           intent(in)    :: runtimelog
    integer,           intent(out)   :: nunit
    if (.not. runtimelog) return
    open (newunit=nunit, file=trim(filename))
  end subroutine ufs_file_setLogUnit

  subroutine ufs_logfhour(msg,hour)
    character(len=*), intent(in) :: msg
    real(r8),         intent(in) :: hour
    character(len=CS)            :: filename
    integer(r4)                  :: nunit
    write(filename,'(a,i3.3)')'log.ice.f',int(hour)
    open(newunit=nunit,file=trim(filename))
    write(nunit,'(a)')'completed: ww3'
    write(nunit,'(a,f10.3)')'forecast hour:',hour
    write(nunit,'(a)')'valid time: '//trim(msg)
    close(nunit)
  end subroutine ufs_logfhour

  ! Define stub routines that do nothing - they are just here to avoid
  ! having cppdefs in the main program
  subroutine shr_file_setLogUnit(nunit)
    integer, intent(in) :: nunit
  end subroutine shr_file_setLogUnit
  subroutine shr_file_getLogUnit(nunit)
    integer, intent(in) :: nunit
  end subroutine shr_file_getLogUnit
  subroutine t_startf(string)
    character(len=*) :: string
  end subroutine t_startf
  subroutine t_stopf(string)
    character(len=*) :: string
  end subroutine t_stopf
  subroutine t_barrierf(string, comm)
    character(len=*) :: string
    integer:: comm
  end subroutine t_barrierf
#endif

end module wav_wrapper_mod
