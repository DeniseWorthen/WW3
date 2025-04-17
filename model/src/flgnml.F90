    if (flgnml) then
      open(newunit=ndsi, file=trim(fnmpre)//"ww3_shel.nml", status='old', iostat=ierr)

      !--------------------
      ! Read namelist
      !--------------------

      call w3nmlshel (mpi_comm, ndsi, trim(fnmpre)//'ww3_shel.nml', nml_domain, nml_input, &
           nml_output_type, nml_output_date, nml_output_path, nml_homog_count, nml_homog_input, ierr)

      !--------------------
      ! 2.1 forcing flags
      !--------------------

      flh(-7:10)  = .false.
      flagtfc(-7) = trim(nml_input%forcing%ice_param1)
      flagtfc(-6) = trim(nml_input%forcing%ice_param2)
      flagtfc(-5) = trim(nml_input%forcing%ice_param3)
      flagtfc(-4) = trim(nml_input%forcing%ice_param4)
      flagtfc(-3) = trim(nml_input%forcing%ice_param5)
      flagtfc(-2) = trim(nml_input%forcing%mud_density)
      flagtfc(-1) = trim(nml_input%forcing%mud_thickness)
      flagtfc(0)  = trim(nml_input%forcing%mud_viscosity)
      flagtfc(1)  = trim(nml_input%forcing%water_levels)
      flagtfc(2)  = trim(nml_input%forcing%currents)
      flagtfc(3)  = trim(nml_input%forcing%winds)
      flagtfc(4)  = trim(nml_input%forcing%ice_conc)
      flagtfc(5)  = trim(nml_input%forcing%atm_momentum)
      flagtfc(6)  = trim(nml_input%forcing%air_density)
      flagtfc(7)  = trim(nml_input%assim%mean)
      flagtfc(8)  = trim(nml_input%assim%spec1d)
      flagtfc(9)  = trim(nml_input%assim%spec2d)

      if (trim(nml_input%forcing%ice_param1) .eq. 'H') then
        flagtfc(-7)='T'
        flh(-7)=.true.
      end if
      if (trim(nml_input%forcing%ice_param2) .eq. 'H') THEN
        flagtfc(-6)='T'
        flh(-6)=.true.
      end if
      if (trim(nml_input%forcing%ice_param3) .eq. 'H') THEN
        flagtfc(-5)='T'
        flh(-5)=.true.
      end if
      if (trim(nml_input%forcing%ice_param4) .eq. 'H') THEN
        flagtfc(-4)='T'
        flh(-4)=.true.
      end if
      if (trim(nml_input%forcing%ice_param5) .eq. 'H') THEN
        flagtfc(-3)='T'
        flh(-3)=.true.
      end if
      if (trim(nml_input%forcing%mud_density) .eq. 'H') THEN
        flagtfc(-2)='T'
        flh(-2)=.true.
      end if
      if (trim(nml_input%forcing%mud_thickness) .eq. 'H') THEN
        flagtfc(-1)='T'
        flh(-1)=.true.
      end if
      if (trim(nml_input%forcing%mud_viscosity) .eq. 'H') THEN
        flagtfc(0)='T'
        flh(0)=.true.
      end if
      if (trim(nml_input%forcing%water_levels) .eq. 'H') THEN
        flagtfc(1)='T'
        flh(1)=.true.
      end if
      if (trim(nml_input%forcing%currents) .eq. 'H') THEN
        flagtfc(2)='T'
        flh(2)=.true.
      end if
      if (trim(nml_input%forcing%winds) .eq. 'H') THEN
        flagtfc(3)='T'
        flh(3)=.true.
      end if
      if (trim(nml_input%forcing%ice_conc) .eq. 'H') THEN
        flagtfc(4)='T'
        flh(4)=.true.
      end if
      if (trim(nml_input%forcing%atm_momentum) .eq. 'H') THEN
        flagtfc(5)='T'
        flh(5)=.true.
      end if
      if (trim(nml_input%forcing%air_density) .eq. 'H') THEN
        flagtfc(6)='T'
        flh(6)=.true.
      end if

      if ( iaproc .eq. napout ) write (ndso, 920)
      DO J=JFIRST, 9
        if (flagtfc(j).eq.'T') THEN
          inflags1(j)=.true.
          flagsc(j)=.false.
        end if
        if (flagtfc(j).eq.'F') THEN
          inflags1(j)=.false.
          flagsc(j)=.false.
        end if
        if (flagtfc(j).eq.'C') THEN
          inflags1(j)=.true.
          flagsc(j)=.true.
        end if
        if ( j .le. 6 ) then
          flh(j) = flh(j) .and. inflags1(j)
        end if
        if ( inflags1(j) ) then
          yesxno = 'YES/--'
        else
          yesxno = '---/NO'
        end IF
        if ( flh(j) ) then
          strng  = '(homogeneous field) '
        else if ( flagsc(j) ) then
          strng  = '(coupling field) '
        else
          strng  = '                    '
        end if
        if ( iaproc .eq. napout ) write (ndso,921) idflds(j), yesxno, strng
      end do
      if (w3_cou_flag) then
        if (flagsc(1) .and. inflags1(2) .and. .not. flagsc(2)) goto 2102
        if (flagsc(2) .and. inflags1(1) .and. .not. flagsc(1)) goto 2102
      end if

      inflags1(10) = .false.
      if (w3_mgw_flag .or. w3_mgp_flag) then
        inflags1(10) = .true.
        flh(10)   = .true.
      end if
      if ( inflags1(10) .and. iaproc.eq.napout ) write (ndso,921) idflds(10), 'yes/--', ' '

      flflg  = inflags1(-7) .or. inflags1(-6) .or. inflags1(-5) .or. inflags1(-4) &
                            .or. inflags1(-3) .or. inflags1(-2) .or. inflags1(-1) &
                            .or. inflags1(0)  .or. inflags1(1)  .or. inflags1(2)  &
                            .or. inflags1(3)  .or. inflags1(4)  .or. inflags1(5)  &
                            .or. inflags1(6)  .or. inflags1(7)  .or. inflags1(8)  &
                            .or. inflags1(9)
      flhom  = flh(-7) .or. flh(-6) .or. flh(-5) .or. flh(-4)              &
                       .or. flh(-3) .or. flh(-2) .or. flh(-1) .or. flh(0)  &
                       .or. flh(1)  .or. flh(2)  .or. flh(3)  .or. flh(4)  &
                       .or. flh(5)  .or. flh(6)  .or. flh(10)

      if ( iaproc .eq. napout ) write (ndso,922)
      ! inflags2 is just "initial value of inflags1", i.e. does *not* get changed when
      ! model reads last record of ice.ww3
      inflags2=inflags1
      if (w3_t_flag) then
        write (ndst,9020) flflg, inflags1, flhom, flh
      end if

      !--------------------
      ! 2.2 Time setup
      !--------------------

      read (nml_domain%start,*) time0
      call t2d(time0,startdate,ierr)
      call d2j(startdate,startjulday,ierr)
      read(nml_domain%stop,*) timen
      call t2d(timen,stopdate,ierr)
      call d2j(stopdate,stopjulday,ierr)

      !--------------------
      ! 2.3 Domain setup
      !--------------------

      iostyp = nml_domain%iostyp
      if (w3_pdlib_flag) then
        if (iostyp .gt. 1) then
          write(*,*) 'iostyp not supported in domain decomposition mode'
          call extcde ( 6666 )
        endif
      end if

      call w3iogr ( 'GRID', ndsm )
      if ( flagll ) then
        factor = 1.
      else
        factor = 1.e-3
      end if

      !--------------------
      ! 2.4 Output dates
      !--------------------

      read(nml_output_date%field%start,    *)   odat(1), odat(2)
      read(nml_output_date%field%stride,   *)   odat(3)
      read(nml_output_date%field%stop,     *)   odat(4), odat(5)

      read(nml_output_date%field%outffile, *)  ofiles(1)

      read(nml_output_date%point%start,    *)   odat(6), odat(7)
      read(nml_output_date%point%stride,   *)   odat(8)
      read(nml_output_date%point%stop,     *)   odat(9), odat(10)

      read(nml_output_date%point%outffile, *)  ofiles(2)

      read(nml_output_date%track%start,      *) odat(11), odat(12)
      read(nml_output_date%track%stride,     *) odat(13)
      read(nml_output_date%track%stop,       *) odat(14), odat(15)

      read(nml_output_date%restart%start,    *) odat(16), odat(17)
      read(nml_output_date%restart%stride,   *) odat(18)
      read(nml_output_date%restart%stop,     *) odat(19), odat(20)

      read(nml_output_date%restart2%start,   *) odat(36), odat(37)
      read(nml_output_date%restart2%stride,  *) odat(38)
      read(nml_output_date%restart2%stop,    *) odat(39), odat(40)

      read(nml_output_date%boundary%start,   *) odat(21), odat(22)
      read(nml_output_date%boundary%stride,  *) odat(23)
      read(nml_output_date%boundary%stop,    *) odat(24), odat(25)

      read(nml_output_date%partition%start,  *) odat(26), odat(27)
      read(nml_output_date%partition%stride, *) odat(28)
      read(nml_output_date%partition%stop,   *) odat(29), odat(30)

      read(nml_output_date%coupling%start,   *) odat(31), odat(32)
      read(nml_output_date%coupling%stride,  *) odat(33)
      read(nml_output_date%coupling%stop,    *) odat(34), odat(35)

      ! set the time stride at 0 or more
      odat(3) = max ( 0 , odat(3) )
      odat(8) = max ( 0 , odat(8) )
      odat(13) = max ( 0 , odat(13) )
      odat(18) = max ( 0 , odat(18) )
      odat(23) = max ( 0 , odat(23) )
      odat(28) = max ( 0 , odat(28) )
      odat(33) = max ( 0 , odat(33) )
      odat(38) = max ( 0 , odat(38) )

      if (w3_cou_flag) then
        ! test the validity of the coupling time step
        if (odat(33) == 0) then
          if ( iaproc .eq. napout ) then
            write(ndso,1010) odat(33), int(dtmax)
          end if
          odat(33) = int(dtmax)
        else if (mod(odat(33),int(dtmax)) .ne. 0) then
          goto 2009
        end if
      end if

      !--------------------
      ! 2.5 Output types
      !--------------------

      npts   = 0
      notype = 6
      if (w3_cou_flag) then
        notype = 7
      end if
      do j = 1, notype

        ! outpts(i)%ofiles(j)=ofiles(j)
        if ( odat(5*(j-1)+3) .ne. 0 ) then

          if ( j .eq. 1 ) then

            ! type 1: fields of mean wave parameters
            fldout = nml_output_type%field%list
            call w3flgrdflag ( ndso, ndso, ndse, fldout, flgd, flgrd, iaproc, napout, ierr )
            if ( ierr .ne. 0 ) goto 2222

          else if ( j .eq. 2 ) then

            ! type 2: point output
            open (newunit=ndsl, file=trim(fnmpre)//trim(nml_output_type%point%file), &
                 form='formatted', status='old', err=2104, iostat=ierr)
            print *,iaproc,ndso,napout,' ZZZ0 '//trim(fnmpre)//trim(nml_output_type%point%file)

            ! first loop to count the number of points
            ! second loop to allocate the array and store the points
            ipts = 0
            do iloop=1,2
              rewind (ndsl)

              if ( iloop.eq.2) then
                npts = ipts
                if ( npts.gt.0 ) then
                  allocate ( x(npts), y(npts), pnames(npts) )
                  ipts = 0 ! reset counter to be reused for next do loop
                else
                  allocate ( x(1), y(1), pnames(1) )
                  goto 2054
                end if
              end if

              do
                read (ndsl,*,err=2004,iostat=ierr) tmpline
                ! if end of file or stopstring, then exit
                if ( ierr.ne.0 .or. index(tmpline,"STOPSTRING").ne.0 ) exit

                ! leading blanks removed and placed on the right
                test = adjustl ( tmpline )
                if ( test(1:1).eq.comstr .or. len_trim(test).eq.0 ) then
                  ! if comment or blank line, then skip
                  cycle
                else
                  ! otherwise, backup to beginning of line
                  backspace ( ndsl, err=2004, iostat=ierr)
                  read (ndsl,*,err=2004,iostat=ierr) xx, yy, pn
                end if
                ipts = ipts + 1
                if ( iloop .eq. 1 ) cycle
                if ( iloop .eq. 2 ) then
                  x(ipts)      = xx
                  y(ipts)      = yy
                  pnames(ipts) = pn
                  if ( iaproc .eq. napout ) then
                    if ( flagll ) then
                      if ( ipts .eq. 1 ) then
                        write (ndso,2945) factor*xx, factor*yy, pn
                      else
                        write (ndso,2946) ipts, factor*xx, factor*yy, pn
                      end if
                    else
                      if ( ipts .eq. 1 ) then
                        write (ndso,2955) factor*xx, factor*yy, pn
                      else
                        write (ndso,2956) ipts, factor*xx, factor*yy, pn
                      end if
                    end if
                  end if
                end if ! iloop.eq.2
              end do ! end of file
            end do ! iloop
            close(ndsl)

          else if ( j .eq. 3 ) then

            ! Type 3: track output
            tflagi = nml_output_type%track%format
            if ( .not. tflagi ) nds(11) = -nds(11)
            if ( iaproc .eq. napout ) then
              if ( .not. tflagi ) then
                write (ndso,3945) 'input', 'unformatted'
              else
                write (ndso,3945) 'input', 'formatted'
              end if
            end if

          else if ( j .eq. 6 ) then

            ! Type 6: partitioning
            iprt(1) = nml_output_type%partition%x0
            iprt(2) = nml_output_type%partition%xn
            iprt(3) = nml_output_type%partition%nx
            iprt(4) = nml_output_type%partition%y0
            iprt(5) = nml_output_type%partition%yn
            iprt(6) = nml_output_type%partition%ny
            prtfrm = nml_output_type%partition%format

            if ( iaproc .eq. napout ) then
              if ( prtfrm ) then
                yesxno = 'YES/--'
              else
                yesxno = '---/NO'
              end if
              write (ndso,6945) iprt, yesxno
            end if

          else if ( j .eq. 7 ) then
#ifdef W3_COU
            ! Type 7: coupling
            fldout = nml_output_type%coupling%sent
            call w3flgrdflag ( ndso, ndso, ndse, fldout, flg2, flgr2, iaproc, napout, ierr )
            if ( ierr .ne. 0 ) goto 2222
            fldin = nml_output_type%coupling%received
            cplt0 = nml_output_type%coupling%couplet0
#endif

          end if ! j
        end if ! odat
      end do ! j

      ! Extra fields to be written in the restart
      fldrst = nml_output_type%restart%extra
      call w3flgrdflag ( ndso, ndso, ndse, fldrst, flogr, flogrr, iaproc, napout, ierr )
      if (present(rstfldlist)) then
        if (trim(fldrst) .ne. 'unset')then
          rstfldlist = trim(fldrst)
        else
          rstfldlist = ' '
        end if
      end if
      if ( ierr .ne. 0 ) goto 2222

      ! force minimal allocation to avoid memory seg fault
      if ( .not.allocated(x) .and. npts.eq.0 ) allocate ( x(1), y(1), pnames(1) )

      !--------------------
      ! 2.6 Homogeneous field data
      !--------------------

      if ( flhom ) then
        if ( iaproc .eq. napout ) write (ndso,951)                   &
             'Homogeneous field data (and moving grid) ...'

        nh(-7) = nml_homog_count%n_ic1
        nh(-6) = nml_homog_count%n_ic2
        nh(-5) = nml_homog_count%n_ic3
        nh(-4) = nml_homog_count%n_ic4
        nh(-3) = nml_homog_count%n_ic5
        nh(-2) = nml_homog_count%n_mdn
        nh(-1) = nml_homog_count%n_mth
        nh(0)  = nml_homog_count%n_mvs
        nh(1)  = nml_homog_count%n_lev
        nh(2)  = nml_homog_count%n_cur
        nh(3)  = nml_homog_count%n_wnd
        nh(4)  = nml_homog_count%n_ice
        nh(5)  = nml_homog_count%n_tau
        nh(6)  = nml_homog_count%n_rho
        nh(10)  = nml_homog_count%n_mov

        n_tot = nml_homog_count%n_tot

        do j=jfirst,10
          if ( nh(j) .gt. nhmax ) goto 2006
        end do

        ! Store homogeneous fields
        if ( n_tot .gt. 0 ) then
          ihh(:)=0
          do ih=1,n_tot
            read(nml_homog_input(ih)%name,*) idtst
            select case (idtst)
            case ('IC1')
              j=-7
            case ('IC2')
              j=-6
            case ('IC3')
              j=-5
            case ('IC4')
              j=-4
            case ('IC5')
              j=-3
            case ('MDN')
              j=-2
            case ('MTH')
              j=-1
            case ('MVS')
              j=0
            case ('LEV')
              j=1
            case ('CUR')
              j=2
            case ('WND')
              j=3
            case ('ICE')
              j=4
            case ('TAU')
              j=5
            case ('RHO')
              j=6
            case ('MOV')
              j=10
            case DEFAULT
              goto 2062
            end SELECT
            ihh(j)=ihh(j)+1
            read(nml_homog_input(ih)%date,*) tho(:,j,ihh(j))
            ha(ihh(j),j) = nml_homog_input(ih)%value1
            hd(ihh(j),j) = nml_homog_input(ih)%value2
            hs(ihh(j),j) = nml_homog_input(ih)%value3
          end do
        end if

        if (w3_o7_flag) then
          do j=jfirst, 10
            if ( flh(j) .and. iaproc.eq.napout ) then
              write (ndso,952) nh(j), idflds(j)
              do i=1, nh(j)
                if ( ( j .le. 1 ) .or. ( j .eq. 4 ) .or. ( j .eq. 6 ) ) then
                  write (ndso,953) i, tho(1,j,i), tho(2,j,i), ha(i,j)
                else if ( ( j .eq. 2 ) .or. ( j .eq. 5 ) .or. ( j .eq. 10 ) ) then
                  write (ndso,953) i, tho(1,j,i), tho(2,j,i), ha(i,j), hd(i,j)
                else if ( j .eq. 3 ) then
                  write (ndso,953) i, tho(1,j,i), tho(2,j,i), ha(i,j), hd(i,j), hs(i,j)
                end if
              end do
            end if
          end do
        end if

        if ( ( flh(-7) .and. (nh(-7).eq.0) ) .or. &
             ( flh(-6) .and. (nh(-6).eq.0) ) .or. &
             ( flh(-5) .and. (nh(-5).eq.0) ) .or. &
             ( flh(-4) .and. (nh(-4).eq.0) ) .or. &
             ( flh(-3) .and. (nh(-3).eq.0) ) .or. &
             ( flh(-2) .and. (nh(-2).eq.0) ) .or. &
             ( flh(-1) .and. (nh(-1).eq.0) ) .or. &
             ( flh(0)  .and. (nh(0).eq.0)  ) .or. &
             ( flh(1)  .and. (nh(1).eq.0)  ) .or. &
             ( flh(2)  .and. (nh(2).eq.0)  ) .or. &
             ( flh(3)  .and. (nh(3).eq.0)  ) .or. &
             ( flh(4)  .and. (nh(4).eq.0)  ) .or. &
             ( flh(5)  .and. (nh(5).eq.0)  ) .or. &
             ( flh(6)  .and. (nh(6).eq.0)  ) .or. &
             ( flh(10) .and. (nh(10).eq.0) ) ) goto 2007

      end if ! flhom

      !--------------------
      ! 2.7 User-defined directory
      !--------------------

      FNMGRD = trim(nml_output_path%grd_out)
      if (FNMGRD(len_trim(FNMGRD):len_trim(FNMGRD)) /= '/') then
        FNMGRD = trim(FNMGRD) // '/'
      end if

      FNMPNT = trim(nml_output_path%pnt_out)
      if (FNMPNT(len_trim(FNMPNT):len_trim(FNMPNT)) /= '/') then
        FNMPNT = trim(FNMPNT) // '/'
      end if

      FNMRST = trim(nml_output_path%rst_out)
      if (FNMRST(len_trim(FNMRST):len_trim(FNMRST)) /= '/') then
        FNMRST = trim(FNMRST) // '/'
      end if

    end if ! flgnml