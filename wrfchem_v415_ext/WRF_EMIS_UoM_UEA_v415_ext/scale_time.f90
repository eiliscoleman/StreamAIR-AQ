
! The RETRO data are already seasonally scaled, only scale NAEI and INERIS
! This routine reads the monthly factors from files into an array

  SUBROUTINE scale_season()
    USE emission
    IMPLICIT NONE

    INTEGER, PARAMETER     :: ncomps_tno=12
    INTEGER, PARAMETER     :: ncomps_naei=7
    INTEGER            :: n, m, ilun, rstat, country, source, ncomps
    REAL               :: temp(nmonths)
    CHARACTER(len=5), DIMENSION(ncomps_tno):: comps_tno
    CHARACTER(len=5), DIMENSION(ncomps_naei):: comps_naei
    CHARACTER(len=160) :: filename
    CHARACTER(len=5)   :: component_name

!'CO','NOx','SOx','NH3','NMVOC','P10','P25','BC_1','EC_1_25','EC_25_10','OC_25','OC_25_10'
!'CO','NOx','SOx','NH3','NMVOC','P10','P25','P25','P25','P10','P25','P10'

    DATA comps_naei /'CO','NOx','SOx','NH3','NMVOC','P10','P25'/ ! Must be same as TNO, INERIS and NAEI!
    DATA comps_tno /'CO','NOx','SOx','NH3','NMVOC','P10','P25','P25','P25','P10','P25','P10'/ ! Must be same as TNO, INERIS and NAEI!
    ilun = 81
    season_fac(:,:,:,:) = 1.
    ! which inventory? if NAEI lets do stuff for the first 7 components, if TNO lets do it for all 12 
    if(naei(1))then
    ncomps=7
    else if (retrotno(1)) then
    ncomps=12
    end if

    DO n = 1, ncomps ! there are 12 components (CO, NOx, SOx,NH3, NMVOC, P10, P25)

       ! Construct filename
       component_name = comps_tno(n)
       filename = TRIM(path)//'scaling_factors/MonthlyFac_'//TRIM(component_name)//'.dat'
       
       WRITE(*,FMT='(A,A)') '   - Reading from file: ', TRIM(filename)

       ! Read from the text file
       OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',STATUS='OLD',IOSTAT=rstat)

       ! Read every line in the file
       DO WHILE (.true.)
          READ(ilun,'(I2,I3,12F7.3)', END=998) country,source,temp(1),temp(2),temp(3),temp(4),temp(5), &
             &temp(6),temp(7),temp(8),temp(9),temp(10),temp(11),temp(12)
!         convert source number to appropriate number as i have added other sub-sources such as tra1, tra2, tra3 etc
          if(source.eq.8) THEN
          source=12
          else if(source.eq.9) THEN
          source=13
          else if(source.eq.10) THEN 
          source=14
          else if(source.eq.11) THEN
          source=15
          else if(source.eq.12) THEN
          source=16
          else if(source.eq.13)THEN 
          source=17
          end if
          ! Debug:
!          WRITE(*,*) country,source,temp
          if(country.ne.0)then
          DO m = 1, nmonths
             season_fac(n,country,source,m) = temp(m)
          END DO
          end if
       END DO
998    CONTINUE
       CLOSE(ilun)
       
       ! Use country code 10 for the whole of Germany (code 60)
!       season_fac(n,60,:,:) = season_fac(n,10,:,:)

    END DO

    WRITE(*,*) '   - Done reading seasonal factors'

  END SUBROUTINE scale_season

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE scale_daily_init()
    USE emission
    IMPLICIT NONE
    INTEGER, PARAMETER     :: ncomps_tno=12
    INTEGER, PARAMETER     :: ncomps_naei=7
    INTEGER            :: ncomps
    INTEGER            :: n, m, ilun, rstat, country, source
    REAL               :: temp(ndays)
    CHARACTER(len=5), DIMENSION(ncomps_tno):: comps_tno
    CHARACTER(len=5), DIMENSION(ncomps_naei):: comps_naei
    CHARACTER(len=160) :: filename
    CHARACTER(len=5)   :: component_name

    DATA comps_naei /'CO','NOx','SOx','NH3','NMVOC','P10','P25'/ ! Must be same as TNO, INERIS and NAEI!
    DATA comps_tno /'CO','NOx','SOx','NH3','NMVOC','P10','P25','P25','P25','P10','P25','P10'/ ! Must be same as TNO, INERIS and NAEI!
    write(6,*)'ini scale_daily_init'
    ilun = 81
    daily_fac(:,:,:,:) = 1.
    ! which inventory? if NAEI lets do stuff for the first 7 components, if TNO lets do it for all 12 
    if(naei(1))then
    ncomps=7
    else if (retrotno(1)) then
    ncomps=12
    end if


    DO n = 1, ncomps ! there are 12 components (CO, NOx, SOx, NH3, NMVOC, P10, P25)

       ! Construct filename
       component_name = comps_tno(n)
       filename = TRIM(path)//'scaling_factors/DailyFac_'//TRIM(component_name)//'.dat'
       
       WRITE(*,FMT='(A,A)') '   - Reading from file: ', TRIM(filename)

       ! Read from the text file
       OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',STATUS='OLD',IOSTAT=rstat)
!       daily_fac(:,:,:,:) = 1.
       ! Read every line in the file
       DO WHILE (.true.)
          READ(ilun,'(I2,I3,F8.4,6F9.4)', END=997) country,source,temp(1),temp(2),temp(3),temp(4),temp(5),temp(6),temp(7)
!         convert source number to appropriate number as i have added other sub-sources such as tra1, tra2, tra3 etc
          if(source.eq.8) THEN
          source=12
          else if(source.eq.9) THEN
          source=13
          else if(source.eq.10) THEN 
          source=14
          else if(source.eq.11) THEN
          source=15
          else if(source.eq.12) THEN
          source=16
          else if(source.eq.13)THEN 
          source=17
          end if
          ! Debug:
!          WRITE(*,*) country,source,temp
          if(country.ne.0)then
          DO m = 1, ndays
             daily_fac(n,country,source,m) = temp(m)
!             write(6,*)'daily_fac',daily_fac(n,country,source,m),n,country,source,m
          END DO
          end if
       END DO
997    CONTINUE
       CLOSE(ilun)
       
    END DO

    WRITE(*,*) '   - Done reading daily factors'

  END SUBROUTINE scale_daily_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE scale_daily(old_value,day_of_week,s,n,i,j)
    USE emission
    IMPLICIT NONE

    ! I/O variables
    INTEGER, INTENT(IN) :: i,j,s,n,day_of_week
    REAL, INTENT(INOUT) :: old_value

    ! Local variables
    INTEGER      :: comp, country_index

    ! Find component
    IF (chem_scheme.eq.4) THEN !CRIMECH scheme

      IF(n .EQ. 4) THEN
         comp = 3 ! SOx
      ELSE IF(n .EQ. 2) THEN  ! NO
         comp = 2 ! NOx
      ELSE IF(n .EQ. 3) THEN  ! NO2
         comp = 2 ! NOx
      ELSE IF(n .EQ. 1) THEN
         comp = 1 ! CO
      ELSE IF(n .EQ. 5) THEN
         comp = 4 ! NH3
      ELSE IF(n .GE. 6 .AND. n .LE. 31) THEN  ! please note this is for listing in CRIMECH_UK speciation!!!!!!
         comp = 5 ! NMVOC
      ELSE IF(n .GE. 32.AND. n .LE. 36) THEN  ! please note this is for listing in CRIMECH_UK speciation!!!!!!
         comp = 7 ! PM25
      ELSE IF(n .GE. 37) THEN  ! please note this is for listing in CRIMECH_UK speciation!!!!!!
         comp = 6 ! PM10
      ELSE
         comp = 0
      END IF 

    ELSE IF (chem_scheme.eq.5) THEN! CBMZ_UK
 
      IF(n .EQ. 2) THEN
         comp = 3 ! SOx
      ELSE IF(n .EQ. 3) THEN  ! NO
         comp = 2 ! NOx
      ELSE IF(n .EQ. 20) THEN  ! NO2
         comp = 2 ! NOx
      ELSE IF(n .EQ. 4) THEN
         comp = 1 ! CO
      ELSE IF(n .EQ. 19) THEN
         comp = 4 ! NH3
      ELSE IF(n .GE. 5 .AND. n .LE. 18) THEN  ! please note this is for listing in CBMZ_UK speciation!!!!!!
         comp = 5 ! NMVOC
      ELSE IF(n .GE. 21 .AND. n .LE. 22) THEN  ! please note this is for listing in CBMZ_UK speciation!!!!!!
         comp = 5 ! NMVOC
      ELSE IF(n .GE. 23.AND. n .LE. 27) THEN  ! please note this is for listing in CBMZ_UK speciation!!!!!!
         comp = 7 ! PM25
      ELSE IF(n .GE. 28) THEN  ! please note this is for listing in CRIMECH_UK speciation!!!!!!
         comp = 6 ! PM10
      ELSE
         comp = 0
      END IF 
    
    ELSE IF (chem_scheme.eq.6) THEN! SAPRC99
      IF(n .EQ. 4) THEN
         comp = 3 ! SOx
      ELSE IF(n .EQ. 2) THEN  ! NO
         comp = 2 ! NOx
      ELSE IF(n .EQ. 3) THEN  ! NO2
         comp = 2 ! NOx
      ELSE IF(n .EQ. 1) THEN
         comp = 1 ! CO
      ELSE IF(n .EQ. 5) THEN
         comp = 4 ! NH3
      ELSE IF(n .GE. 6 .AND. n .LE. 40) THEN  ! please note this is for listing in SAPRC99_UK speciation!!!!!!
         comp = 5 ! NMVOC
      ELSE IF(n .GE. 41.AND. n .LE. 45) THEN  ! please note this is for listing in SAPRC99_UK speciation!!!!!!
         comp = 7 ! PM25
      ELSE IF(n .GE. 46) THEN  ! please note this is for listing in SAPRC99_UK speciation!!!!!!
         comp = 6 ! PM10
      ELSE
         comp = 0
      END IF     
    ELSE  ! NOT CRIMECH NOR CBMZ NOR SAPRC99 STOP!!!
!       write(*,*) 'HOLD ON A MINUTE! DAILY SCALING NOT IMPLEMENTED FOR CHEM_SCHEME ', chem_scheme
!       STOP
    END IF



!    write(6,*)'in scale_daily'
    ! Find country code
    if(retrotno(1))then 
    country_index = wrf_countries(i,j)
    else  ! naei
    country_index = 27
    end if
!    write(6,*)'country_index=',country_index
    ! Do the scaling

    IF(comp .NE. 0 .AND. country_index .LE. 63) THEN
       old_value = old_value*daily_fac(comp,country_index,s,day_of_week)
!     write(6,*)'daily_fac=',i,j,daily_fac(comp,country_index,s,day_of_week),comp,country_index,s,day_of_week
    END IF

  END SUBROUTINE scale_daily

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE scale_hour(old_value,new_value,hour,s,i,j)
    USE emission
    IMPLICIT NONE

    ! I/O variables
    INTEGER, INTENT(IN) :: i,j,s,hour
    REAL, INTENT(IN)    :: old_value
    REAL, INTENT(OUT)   :: new_value

    ! Local variables
    INTEGER      :: localtime, dst
    REAL         :: factor(nsources,24)

    ! Public power stations
    DATA factor(1,:) / 0.79,0.72,0.72,0.71,0.74,0.8,0.92,1.08,1.19,1.22,1.21,1.21, &
         1.17,1.15,1.14,1.13,1.1,1.07,1.04,1.02,1.02,1.01,0.96,0.88 /
    ! Comm./inst. combustion
    DATA factor(2,:) / 0.38,0.36,0.36,0.36,0.37,0.5,1.19,1.53,1.57,1.56,1.35,1.16, &
         1.07,1.06,1,0.98,0.99,1.12,1.41,1.52,1.39,1.35,1,0.42 /
    ! Industrial combustion
    DATA factor(3,:) / 0.75,0.75,0.78,0.82,0.88,0.95,1.02,1.09,1.16,1.22,1.28,1.3, &
         1.22,1.24,1.25,1.16,1.08,1.01,0.95,0.9,0.85,0.81,0.78,0.75 /
    ! Production processes
    DATA factor(4,:) / 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /
    ! Extraction fossil fuels
    DATA factor(5,:) / 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /
    ! Solvents
    DATA factor(6,:) / 0.5,0.35,0.2,0.1,0.1,0.2,0.75,1.25,1.4,1.5,1.5,1.5,1.5,1.5, &
         1.5,1.5,1.5,1.4,1.25,1.1,1,0.9,0.8,0.7 /
    ! Road traffic gasoline
    DATA factor(7,:) / 0.19,0.09,0.06,0.05,0.09,0.22,0.86,1.84,1.86,1.41,1.24,1.2, &
         1.32,1.44,1.45,1.59,2.03,2.08,1.51,1.06,0.74,0.62,0.61,0.44 /
    ! Road traffic diesel
    DATA factor(8,:) / 0.19,0.09,0.06,0.05,0.09,0.22,0.86,1.84,1.86,1.41,1.24,1.2, &
         1.32,1.44,1.45,1.59,2.03,2.08,1.51,1.06,0.74,0.62,0.61,0.44 /
    ! Road traffic LPG
    DATA factor(9,:) / 0.19,0.09,0.06,0.05,0.09,0.22,0.86,1.84,1.86,1.41,1.24,1.2, &
         1.32,1.44,1.45,1.59,2.03,2.08,1.51,1.06,0.74,0.62,0.61,0.44 /
    ! Road traffic non-exhaust (volatilisation)
    DATA factor(10,:) / 0.19,0.09,0.06,0.05,0.09,0.22,0.86,1.84,1.86,1.41,1.24,1.2, &
         1.32,1.44,1.45,1.59,2.03,2.08,1.51,1.06,0.74,0.62,0.61,0.44 /
    ! Road traffic non-exhaust (road, tire wear)
    DATA factor(11,:) / 0.19,0.09,0.06,0.05,0.09,0.22,0.86,1.84,1.86,1.41,1.24,1.2, &
         1.32,1.44,1.45,1.59,2.03,2.08,1.51,1.06,0.74,0.62,0.61,0.44 /
    ! Other mobile (trains, plains etc.)
    DATA factor(12,:) / 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /
    ! Waste
    DATA factor(13,:) / 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /
    ! Agriculture
    DATA factor(14,:) / 0.6,0.6,0.6,0.6,0.6,0.65,0.75,0.9,1.1,1.35,1.45,1.6,1.8, &
         1.75,1.7,1.55,1.35,1.1,0.9,0.75,0.65,0.6,0.6,0.6 /
    ! Nature
    DATA factor(15,:) / 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /


    ! Calculate local time (hour is UTC, timezone is difference from UTC)
    IF(daylight_savings) THEN
       dst = 1
    ELSE
       dst = 0
    END IF

    IF(retrotno(1).and. scaling_emep_hour) THEN
       localtime = hour + timediff_from_utc(i,j)
    ELSE IF(scaling_hour) THEN
       localtime = hour + timezone + dst
    END IF

    IF(localtime .GT. 24) THEN
       localtime = localtime - 24
    ELSE IF(localtime .LT. 1) THEN
       localtime = localtime + 24
    END IF
    new_value = old_value*factor(s,localtime)

  END SUBROUTINE scale_hour

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE get_wrfcountries_tno()
    USE emission
    USE module_tno
    IMPLICIT NONE

    INTEGER :: i, j, it, jt, country_index, dst, timezones(73)
    REAL    :: lat, lon

    ! Values taken from EMEP Unified model
    DATA timezones / &
         1,1,1,2,1,1,2,1,1,1,2,1,0,0,1,1,1,1,1,1,2,1,1,1,2, & ! countries 1-25
         3,0,1,1,1,1,1,1,1,1,3,3,3,2,2,2,4,2,2,2,1,1,1,1,1, & ! countries 26-50
         1,1,6,4,2,4,1,0,1,1,3,1,1,1,1,1,0,6,3,1,4,1,1 /      ! countries 51-73

    WRITE(*,*) ' - Find country numbers for wrf domain'

    ! Go through every WRF grid box
    DO j = 1, jx
       DO i = 1, ix

          lat = wrf_lats(i,j)
          lon = wrf_lons(i,j)

          ! Calculate indices of TNO inventory
          it = NINT((lon*8.) - t_istart + 0.5)
          jt = NINT((lat*16.) - t_jstart + 0.5)

          IF(it .GE. 1 .AND. it .LE. t_iemis .AND. jt .GE. 1 .AND. jt .LE. t_jemis) THEN
             wrf_countries(i,j) = tno_countries(it,jt)
          END IF
       END DO
    END DO

    ! Calculate time difference from UTC
    IF(scaling_emep_hour) THEN
       IF(daylight_savings) THEN
          dst = 1
       ELSE
          dst = 0
       END IF

       ! Go through every WRF grid box
       DO j = 1, jx
          DO i = 1, ix
             country_index = wrf_countries(i,j)
             IF(country_index .EQ. 67 .OR. country_index .GT. 73) THEN
                ! 67 is undefined, use timezone from namelist instead
!                write(6,*)'oops country index=',country_index
                timediff_from_utc(i,j) = timezone + dst
             ELSE
                ! Find timezone value and consider daylight saving time
!                write(6,*)'ok country index and timezone=',country_index,timezones(country_index),j,i
                timediff_from_utc(i,j) = timezones(country_index) + dst
             END IF
          END DO
       END DO

    END IF

  END SUBROUTINE get_wrfcountries_tno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE get_wrfcountries_ineris()
    USE emission
    USE module_ineris
    IMPLICIT NONE

    INTEGER :: i, j, ii, ji, country_index, dst, timezones(73)
    REAL    :: lat, lon

    ! Values taken from EMEP Unified model
    DATA timezones / &
         1,1,1,2,1,1,2,1,1,1,2,1,0,0,1,1,1,1,1,1,2,1,1,1,2, & ! countries 1-25
         3,0,1,1,1,1,1,1,1,1,3,3,3,2,2,2,4,2,2,2,1,1,1,1,1, & ! countries 26-50
         1,1,6,4,2,4,1,0,1,1,3,1,1,1,1,1,0,6,3,1,4,1,1 /      ! countries 51-73

    WRITE(*,*) ' - Find country numbers for wrf domain'

    ! Go through every WRF grid box
    DO j = 1, jx
       DO i = 1, ix

          lat = wrf_lats(i,j)
          lon = wrf_lons(i,j)

          ! Calculate indices of INERIS inventory
          ii = NINT((lon*10.) - i_istart + 1)
          ji = NINT((lat*10.) - i_jstart + 1)

          IF(ii .GE. 1 .AND. ii .LE. i_iemis .AND. ji .GE. 1 .AND. ji .LE. i_jemis) THEN
             wrf_countries(i,j) = ineris_countries(ii,ji)
          END IF

       END DO
    END DO

    ! Calculate time difference from UTC
    IF(scaling_emep_hour) THEN
       IF(daylight_savings) THEN
          dst = 1
       ELSE
          dst = 0
       END IF

       ! Go through every WRF grid box
       DO j = 1, jx
          DO i = 1, ix
             country_index = wrf_countries(i,j)
             IF(country_index .EQ. 67 .OR. country_index .GT. 73) THEN
                ! 67 is undefined, use timezone from namelist instead
                timediff_from_utc(i,j) = timezone + dst
             ELSE
                ! Find timezone value and consider daylight saving time
                timediff_from_utc(i,j) = timezones(country_index) + dst
             END IF
          END DO
       END DO

    END IF

  END SUBROUTINE get_wrfcountries_ineris

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
