MODULE module_ineris
  USE emission

!  INTEGER, PARAMETER :: i_iemis = 601, i_jemis = 451, i_comps = 4, &
!  INTEGER, PARAMETER :: i_iemis = 601, i_jemis = 451, i_comps = 7, &
!       i_istart = -200, i_iend = 400, i_jstart = 300, i_jend = 750
  INTEGER, PARAMETER :: i_iemis = 601, i_jemis = 471, i_comps = 4, &
       i_istart = -200, i_iend = 400, i_jstart = 280, i_jend = 750
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ineris_plot
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ineris_readdata
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ineris_data
  REAL, ALLOCATABLE, DIMENSION(:,:) :: i_longitude
  REAL, ALLOCATABLE, DIMENSION(:,:) :: i_latitude
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: ineris_covered
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: ineris_lanuv
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: ineris_naei
  INTEGER            :: ineris_countries(i_iemis,i_jemis)
  CHARACTER(len=4)   :: syear
  CHARACTER(len=5), DIMENSION(i_comps):: components_ineris

  DATA components_ineris /'CO','NOx','SOx','NMVOC'/ ! NMVOC must always be last
!  DATA components_ineris /'CO','NH3','NOx','PM25','PMco','SOx','NMVOC'/

CONTAINS

  SUBROUTINE allocate_ineris()
  ! Allocate arrays for INERIS emissions
  ! Using allocatable arrays to avoid compiler limits on array size
  ! Note that if array is too big for available RAM on the system an out-of-memory error will be generated.
  ! T. Pugh, 20/09/10
    USE emission
    IMPLICIT NONE

    ALLOCATE(ineris_plot(i_iemis,i_jemis,i_comps))
    ALLOCATE(ineris_readdata(i_iemis,i_jemis,i_comps,nsources))
    ALLOCATE(ineris_data(i_iemis,i_jemis,nchem,nsources))
    ALLOCATE(i_longitude(i_iemis,i_jemis))
    ALLOCATE(i_latitude(i_iemis,i_jemis))
    ALLOCATE(ineris_covered(i_iemis,i_jemis))
    ALLOCATE(ineris_lanuv(i_iemis,i_jemis))
    ALLOCATE(ineris_naei(i_iemis,i_jemis))

  END SUBROUTINE allocate_ineris

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_ineris()
    USE emission
    IMPLICIT NONE

    INTEGER            :: i,j,n,c,s,ilun,rstat,country
    REAL               :: lat,lon,source(nsources),sum_sources
    CHARACTER(len=160) :: filename
    CHARACTER(len=5)   :: component_name

    ! Now we will read the variables from the INERIS inventory
    PRINT *, ' - Reading INERIS emissions from txt files'

    ilun = 21
    ineris_covered  = .false.
    ineris_plot     = 0.D0
    ineris_readdata = 0.D0
    ineris_countries= 67   ! initialize everything to 67 (undefined)

    DO n = 1, i_comps      ! there are 4 components (CO, NMVOC, NOX, SOx)

       ! Construct filename
       WRITE(syear,'(I4)') year
       component_name = components_ineris(n)
       filename = TRIM(path)//'ineris/'//syear//'/grid'//TRIM(component_name)//'_emis'//syear//'_INERIS_10km_latlon'

       WRITE(*,FMT='(A,A)') '   - Reading from file: ', TRIM(filename)

       ! Read from the text file. The unit is Mg/(0.1deg x 0.1deg)/year
       OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',STATUS='OLD',IOSTAT=rstat)

       ! Read every point in the inventory:
       DO WHILE (.true.)
!       DO c = 1,nlines(n)
          READ(ilun,FMT=*,END=999) country,lat,lon,source(1),source(2), &
source(3),source(4),source(5),source(6),source(7),source(8), &
source(9),source(10),source(11)

          ! Calculate index
          i = NINT((lon*10.) - i_istart + 1)
          j = NINT((lat*10.) - i_jstart + 1)

          ineris_covered(i,j) = .true.

          ! Season scaling (pr. month)
          IF(scaling_season .AND. country .LE. ncountries) THEN
             DO s=1,nsources
                ! Debug:
!                WRITE(*,*) 'Before seasonal factor applied: ',source(s)
                source(s) = source(s)*season_fac(n,country,s,month)
!                WRITE(*,*) 'After seasonal factor applied: ',source(s)
             END DO
          END IF

          sum_sources = 0.
          DO s=1,nsources
             sum_sources = sum_sources + source(s)
          END DO

          ! BUG FIX: The codes for maritime zones have both 2 and 3 digits for
          ! the years 2006 and 2007. We skip the 3 digit zones to avoid doubled
          ! ship emissions
          IF(country .GT. 300 .AND. country .LT. 350) THEN
             source(:) = 0.D0
             sum_sources = 0.
          END IF

          ! Store the sum in an array to be used by plot_ineris()
          ineris_plot(i,j,n) = ineris_plot(i,j,n) + sum_sources

          ! Store the sources in an array
          DO s=1,nsources
             ineris_readdata(i,j,n,s) = ineris_readdata(i,j,n,s) + source(s)
          END DO

          ! Store the country number in an array
          ineris_countries(i,j) = country
       END DO
999    CONTINUE
       CLOSE(ilun)

    END DO ! n = 1, i_comps


    ! Short summary
    DO n = 1, i_comps
       sum_sources = 0.D0
       DO j = 1, i_jemis
          DO i = 1, i_iemis
             sum_sources = sum_sources + ineris_plot(i,j,n)
          END DO
       END DO
       WRITE(*,FMT='(3A,F15.3,A)') 'INERIS: Sum of ',components_ineris(n),' = ',sum_sources,' t/yr'
    END DO

    WRITE(*,*) ' - Done reading from txt files'

  END SUBROUTINE read_ineris

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE latlon_ineris()
    IMPLICIT NONE

    ! We want to write lats and longs to file for the whole INERIS grid
    INTEGER            :: i,j,ilun,rstat
    REAL               :: lon,lat
    CHARACTER(len=160) :: filename

    PRINT *, ' - Write latitudes and longitudes for INERIS grid'

    ! Construct filename
    filename = './plotfiles/ineris_latlon.dat'
    WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

    ! Open file for writing
    ilun = 41
    OPEN(ilun,FILE=filename,FORM='FORMATTED',STATUS='UNKNOWN',IOSTAT=rstat)

    ! Loop through INERIS grid
    DO j = i_jstart, i_jend
       DO i = i_istart, i_iend
          lon = i/10.
          lat = j/10.

          WRITE(ilun,FMT='(2F12.1)') j/10.,i/10.  ! Write to file for plotting
          i_longitude(i-i_istart+1,j-i_jstart+1) = lon
          i_latitude(i-i_istart+1,j-i_jstart+1) = lat
       END DO
    END DO

    CLOSE(ilun)

    WRITE(*,*) ' - Done writing latitudes and longitudes'

  END SUBROUTINE latlon_ineris

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE plot_ineris()
    IMPLICIT NONE
    CHARACTER (len=80)          :: component_name,filename
    CHARACTER (len=2)           :: smonth
    INTEGER                     :: i,j,n,ilun,rstat

    ilun = 51

    IF(month .LT. 10) THEN
       WRITE(smonth,'(I1)') month
    ELSE
       WRITE(smonth,'(I2)') month
    END IF

    ! Units are in t/year if we run this routine before convert_units()
    DO n = 1, i_comps
       component_name = components_ineris(n)

       IF(scaling_season) THEN
          filename = './plotfiles/ineris_'//syear//'_'//TRIM(smonth)//'_'//TRIM(component_name)//'.dat'
       ELSE
          filename = './plotfiles/ineris_'//syear//'_'//TRIM(component_name)//'.dat'
       END IF

       WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

       OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',IOSTAT=rstat)
       DO j=1,i_jemis
          DO i=1,i_iemis
             WRITE(UNIT=ilun,FMT='(F16.4)',IOSTAT=rstat) ineris_plot(i,j,n)
          END DO
       END DO
       CLOSE(UNIT=ilun)
    END DO

    ! Also write country number to plot file
    filename = './plotfiles/ineris_countries_'//syear//'.dat'
    WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

    OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',IOSTAT=rstat)
    DO j=1,i_jemis
       DO i=1,i_iemis
          WRITE(UNIT=ilun,FMT='(I10)',IOSTAT=rstat) ineris_countries(i,j)
       END DO
    END DO
    CLOSE(UNIT=ilun)

  END SUBROUTINE plot_ineris

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE convert_ineris()
    IMPLICIT NONE

    ! Convert the gas emissions from tons/(0.1dx0.1d)/yr to mole/(h*km^2):
    ! Have to multiply by 1e6 to convert from tons to g.
    ! Have to divide by 365*24 to convert from pr. year to pr. hour.
    ! Have to divide by gridbox_area to get pr. square km.
    ! Have to divide by molar mass to convert from g to mole.
    REAL        :: constant,earth_radius,pi,lat,gridbox_area,mw(3)
    INTEGER ineris_to_radm(3)
    INTEGER     :: i,j,n, s
    DATA mw / 28.0101, 46.0055, 64.054 / ! CO, NOx (NO2), SOx (SO2)
    DATA ineris_to_radm / 11, 2, 1 /     ! CO, NOx, SO2

    ineris_data = 0.d0
    earth_radius = 6371.     ! Earth's radius in km
    pi = 3.14159265         ! Pi
    constant = 1000000./(365.*24.)

    DO n = 1, i_comps-1  ! Skip the VOCs for now
       DO j = 1, i_jemis
          DO i = 1, i_iemis

             ! Calculate area of 0.1 deg x 0.1 deg grid box
             lat = i_latitude(i,j)
             gridbox_area = (earth_radius**2)*((0.1*pi/180.)**2)*COS((lat*pi)/180.)
             
             ! Debug:
!             WRITE(*,*) 'lat: ',lat,', gridbox_area: ',gridbox_area
!             WRITE(*,*) 'n=',n,', i=',i,', j=',j,', before conversion: ',SUM(ineris_readdata(i,j,n,:)),' tons/yr'

             ! Convert from tons/(0.1degx0.1deg)/yr to mole/(km^2*h)
             ineris_data(i,j,ineris_to_radm(n),:) = &
                  (constant*ineris_readdata(i,j,n,:))/(gridbox_area*mw(n))

             ! Debug:
!             WRITE(*,*) 'n=',n,', i=',i,', j=',j,', after conversion: ',SUM(ineris_data(i,j,ineris_to_radm(n),:)),' mole/km2/yr'

          END DO
       END DO
    END DO

  END SUBROUTINE convert_ineris

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE vocs_ineris()
    USE voc_factors
    IMPLICIT NONE

    ! Divide VOCs into RADM-components and convert the units

    REAL              :: mw(nvocs),sum_vocs
    REAL              :: constant,earth_radius,pi,lat,gridbox_area,value
    INTEGER           :: i,j,n,s,radm_nr(nvocs)

    ! Molecular weight of 47 of the 50 most significant NMVOC species:
    DATA mw / &
         46.06,58.12,30.07,44.09,92.14,32.04,28.05,72.15,72.15,58.09, &
         86.18,106.16,58.12,30.03,78.11,42.08,72.11,116.16,142.29,106.17, &
         120.19,60.10,100.21,88.11,106.16,106.16,114.23,100.16,26.04,128.20, &
         56.11,156.31,74.08,74.12,44.05,86.20,118.18,120.19,60.09,54.09, &
         136.24,56.10,90.12,120.20,120.20,156.31,70.14 /

    DATA radm_nr / &
         7,7,10,7,15,7,12,8,8,17,9,16,7,4,15,13,17,5,9,15,15,7,9,5,16,16,9,17, &
         7,9,13,9,5,7,3,8,7,15,7,14,14,14,7,15,15,9,14 /

    sum_vocs = 0.
    earth_radius = 6371.     ! Earth's radius in km
    pi = 3.14159265         ! Pi
    constant = 1000000./(365.*24.)

    WRITE(*,*) ' - Dividing VOC species based on sources'
    ! HVA MED SOURCE 10 og 11??
    ! For now we treat SNAP 10 and 11 as SNAP sector 1: Combustion in energy and transformation industries

    ! SJEKK AT DETTE BLIR RIKTIG!!!!!!!!

    DO j=1,i_jemis
       DO i=1,i_iemis
          ! Calculate area of 0.1 deg x 0.1 deg grid box
          lat = i_latitude(i,j)
          gridbox_area = (earth_radius**2)*((0.1*pi/180.)**2)*COS((lat*pi)/180.)

          DO n=1,nvocs
             DO s=1,nvocsrc
                ! Divide VOCs
                value = ineris_readdata(i,j,i_comps,s)*factor(n,s)
                sum_vocs = sum_vocs + value

                ! Convert units
                value = (constant*value)/(gridbox_area*mw(n))
                ineris_data(i,j,radm_nr(n),s) = &
                     ineris_data(i,j,radm_nr(n),s) + value
             END DO

             ! Then treat agriculture+other (SNAP 10+11). For now treated as SNAP 1
             DO s=10,11
                ! Divide VOCs
                value = ineris_readdata(i,j,i_comps,s)*factor(n,1) ! OBS: factor(n,1)
                sum_vocs = sum_vocs + value

                ! Convert units
                value = (constant*value)/(gridbox_area*mw(n))
                ineris_data(i,j,radm_nr(n),s) = &
                     ineris_data(i,j,radm_nr(n),s) + value
             END DO

          END DO
       END DO
    END DO

    WRITE(*,*) '  - Sum VOC before speciation: ', SUM(ineris_readdata(:,:,i_comps,:)), ' tons/yr'
    WRITE(*,*) '  - Sum VOC after speciation:  ', sum_vocs, ' tons/yr'

    WRITE(*,*) ' - Done dividing VOC species'

  END SUBROUTINE vocs_ineris

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE overlap_ineris_lanuv()
    USE module_lanuv
    IMPLICIT NONE

    INTEGER                 :: i,j,il,jl,ilun
    REAL                    :: low_lon,high_lon,low_lat,high_lat,lon,lat
    LOGICAL                 :: covered, covered_completely
    CHARACTER(len=160)      :: filename

    ilun = 51
    ineris_lanuv    = .false.    ! .true. if the whole INERIS box is covered by LANUV

    filename = './plotfiles/ineris_lanuv.dat'  ! Construct filename
    WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

    OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED')

    DO j=1,i_jemis
       DO i=1,i_iemis
          
          ! Find the boundaries of the INERIS grid boxes
          low_lon = i_longitude(i,j)-0.05
          high_lon = i_longitude(i,j)+0.05
          low_lat = i_latitude(i,j)-0.05
          high_lat = i_latitude(i,j)+0.05

          ! Then loop through all LANUV grid boxes
          covered = .false.
          covered_completely = .true.

          DO jl=1,l_jemis
!             DO il=1,l_iemis
             DO il=1,l_iemis_short
                lon = l_longitude(il,jl)
                lat = l_latitude(il,jl)

                ! If within INERIS grid box
                IF((lon .GE. low_lon) .AND. (lon .LT. high_lon) .AND. &
                     (lat .GE. low_lat) .AND. (lat .LT. high_lat)) THEN
                   covered = .true.

                   IF(.NOT. lanuv_covered(il,jl)) THEN
                      covered_completely = .false.
                   END IF
                END IF

             END DO
          END DO

          IF(covered .AND. covered_completely) THEN
             ineris_lanuv(i,j) = .true.
             WRITE(ilun,'(I1)') 1
          ELSE
             ineris_lanuv(i,j) = .false.
             WRITE(ilun,'(I1)') 0
          END IF
       END DO
    END DO

    CLOSE(ilun)

  END SUBROUTINE overlap_ineris_lanuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_overlap_ineris_lanuv()
    USE module_lanuv
    IMPLICIT NONE

    INTEGER                 :: i,j,icovered,ilun
    CHARACTER(len=160)      :: filename

    ilun = 51
    ineris_lanuv    = .false.    ! .true. if the whole INERIS box is covered by LANUV

    filename = './plotfiles/ineris_lanuv.dat'  ! Construct filename
    WRITE(*,FMT='(A,A)') '   - Reading from file: ', TRIM(filename)

    OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED')

    DO j=1,i_jemis
       DO i=1,i_iemis
          READ(ilun,*) icovered
          IF(icovered .EQ. 1) THEN
             ineris_lanuv(i,j) = .true.
          ELSE
             ineris_lanuv(i,j) = .false.
          END IF
       END DO
    END DO

    CLOSE(ilun)

  END SUBROUTINE read_overlap_ineris_lanuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE overlap_ineris_naei()
    USE module_naei
    IMPLICIT NONE

    INTEGER                 :: i,j,in,jn,ilun
    REAL                    :: low_lon,high_lon,low_lat,high_lat,lon,lat
    LOGICAL                 :: covered, covered_partly
    CHARACTER(len=160)      :: filename

    ilun = 51
    ineris_naei    = .false.    ! .true. if some of the INERIS box is covered by NAEI

    filename = './plotfiles/ineris_naei.dat'  ! Construct filename
    WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

    OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED')

    DO j=1,i_jemis
       DO i=1,i_iemis
          
          ! Find the boundaries of the INERIS grid boxes
          low_lon = i_longitude(i,j)-0.05
          high_lon = i_longitude(i,j)+0.05
          low_lat = i_latitude(i,j)-0.05
          high_lat = i_latitude(i,j)+0.05

          ! Then loop through all NAEI grid boxes
          covered = .false.
          covered_partly = .false.

          DO jn=1,n_jemis
             DO in=1,n_iemis
                lon = n_longitude(in,jn)
                lat = n_latitude(in,jn)

                ! If within INERIS grid box
                IF((lon .GE. low_lon) .AND. (lon .LT. high_lon) .AND. &
                     (lat .GE. low_lat) .AND. (lat .LT. high_lat)) THEN
                   covered = .true.

                   IF(naei_covered(in,jn)) THEN
                      covered_partly = .true.
                   END IF
                END IF

             END DO
          END DO

          IF(covered .AND. covered_partly) THEN
             ineris_naei(i,j) = .true.
             WRITE(ilun,'(I1)') 1
          ELSE
             ineris_naei(i,j) = .false.
             WRITE(ilun,'(I1)') 0
          END IF
       END DO
    END DO

    CLOSE(ilun)

  END SUBROUTINE overlap_ineris_naei

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_overlap_ineris_naei()
    USE module_naei
    IMPLICIT NONE

    INTEGER                 :: i,j,icovered,ilun
    CHARACTER(len=160)      :: filename

    ilun = 51
    ineris_naei    = .false.    ! .true. if some of the INERIS box is covered by NAEI

    filename = './plotfiles/ineris_naei.dat'  ! Construct filename
    WRITE(*,FMT='(A,A)') '   - Reading from file: ', TRIM(filename)

    OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED')

    DO j=1,i_jemis
       DO i=1,i_iemis
          READ(ilun,*) icovered
          IF(icovered .EQ. 1) THEN
             ineris_naei(i,j) = .true.
          ELSE
             ineris_naei(i,j) = .false.
          END IF
       END DO
    END DO

    CLOSE(ilun)

  END SUBROUTINE read_overlap_ineris_naei

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE merge_lanuv_ineris()
    USE module_lanuv
    IMPLICIT NONE

    INTEGER, PARAMETER:: nspecies = 13
    REAL              :: lon, lat, factor_ineris(nradm,nsources),earth_radius, &
         pi, gridbox_area,factor_lanuv(nradm,nsources),factor,total_voc(3) ! TEMPORARY!!!!
    INTEGER           :: i,j,k,n,s,ii,ji,species(nspecies)

    ! The following SNAP categories are missing from the LANUV data: 5,6,9,10,11
    ! VOC species:
!    DATA species / 3,4,5,7,8,9,10,12,13,14,15,16,17 /  ! TEMPORARY!!!!

    factor_ineris = 0.d0
    factor_lanuv = 0.d0
    earth_radius = 6371.     ! Earth's radius in km
    pi = 3.14159265          ! Pi

    WRITE(*,*) ' - Merging LANUV data with INERIS data'


    ! Calculate INERIS values for each species and sector in area covered by LANUV data
    DO s=1,nsources
       DO n=1,nradm
          DO j=1,i_jemis
             DO i=1,i_iemis
                IF(ineris_lanuv(i,j)) THEN

                   ! Calculate area of 0.1 deg x 0.1 deg grid box
                   lat = i_latitude(i,j)
                   gridbox_area = (earth_radius**2)*((0.1*pi/180.)**2)*COS((lat*pi)/180.)

                   ! Add to total in unit mole/hr
                   factor_ineris(n,s) = factor_ineris(n,s) + ineris_data(i,j,n,s)*gridbox_area

                END IF
             END DO
          END DO
       END DO
    END DO

    ! Write total INERIS data to screen
!    WRITE(*,*) 'Total INERIS data (mole/hr) in area where we have LANUV coverage:'
!    DO n=1,nradm
!       WRITE(*,*) ename(n)
!       WRITE(*,*) factor_ineris(n,:)
!    END DO


! TEMPORARY:
!    total_voc = 0.d0

    ! Write total for sectors 1,3, and 4
!    DO n=1,nspecies
!       total_voc(1) = total_voc(1) + factor_ineris(species(n),1)
!       total_voc(2) = total_voc(2) + factor_ineris(species(n),3)
!       total_voc(3) = total_voc(3) + factor_ineris(species(n),4)
!    END DO

!    WRITE(*,*) 'total_voc:      ',total_voc
!    WRITE(*,*) 'SUM(total_voc): ',SUM(total_voc)

! END OF TEMPORARY



    ! Add INERIS data to missing sectors in LANUV - go through all LANUV grid boxes
    DO s = 1, nsources
       DO n = 1, nradm
          DO j=1,l_jemis
             DO i=1,l_iemis_short

                ! Find the INERIS grid box that covers the LANUV grid box
                ! First find lat/lon
                lon = l_longitude(i,j)
                lat = l_latitude(i,j)

                ! Then calculate indices of INERIS inventory
                ii = NINT((lon*10.) - i_istart + 1)
                ji = NINT((lat*10.) - i_jstart + 1)

                ! Some sectors and SO2 is not included in LANUV - use INERIS data
                IF((n .EQ. 1) .OR. (s .EQ. 5) .OR. (s .EQ. 6) .OR. &
                     (s .EQ. 9) .OR. (s .EQ. 10) .OR. (s .EQ. 11)) THEN
                   ! Add INERIS value to LANUV arrays
                   lanuv_data(i,j,n,s) = ineris_data(ii,ji,n,s)
                END IF
                
             END DO
          END DO
       END DO
    END DO


    ! Calculate factors for LANUV data

    gridbox_area = 1. ! LANUV data is given in per sq. km.

    ! Calculate LANUV values for each species and sector
    DO s=1,nsources
       DO n=1,nradm
          DO j=1,l_jemis
             DO i=1,l_iemis

                ! Find the INERIS grid box that covers the LANUV grid box
                ! First find lat/lon
                lon = l_longitude(i,j)
                lat = l_latitude(i,j)

                ! Then calculate indices of INERIS inventory
                ii = NINT((lon*10.) - i_istart + 1)
                ji = NINT((lat*10.) - i_jstart + 1)
                
                ! Then check if INERIS and LANUV overlap completely
                IF(ineris_lanuv(ii,ji)) THEN

                   ! Add to total in unit mole/hr
                   factor_lanuv(n,s) = factor_lanuv(n,s) + lanuv_data(i,j,n,s)*gridbox_area

                   ! Then add industry
                   IF(s .LE. 4) THEN
                      DO k=1,kx
                         factor_lanuv(n,s) = factor_lanuv(n,s) + lanuv_industry(i,k,j,n,s)*gridbox_area
                      END DO
                   END IF

                END IF
             END DO
          END DO
       END DO
    END DO

    ! Write total LANUV data to screen
!    WRITE(*,*)
!    WRITE(*,*) 'Total LANUV data (mole/hr) in area where we have complete LANUV coverage in INERIS:'
!    DO n=1,nradm
!       WRITE(*,*) ename(n)
!       WRITE(*,*) factor_lanuv(n,:)
!    END DO


    ! Scale LANUV data based on INERIS data

    ! Calculate LANUV values for each species and sector
    DO s=1,nsources
       DO n=1,nradm

          IF(factor_lanuv(n,s) .GT. 0.) THEN
             factor = factor_ineris(n,s)/factor_lanuv(n,s)
          ELSE
             factor = 1.
          END IF

          DO j=1,l_jemis
             DO i=1,l_iemis

                ! Find the INERIS grid box that covers the LANUV grid box
                ! First find lat/lon
                lon = l_longitude(i,j)
                lat = l_latitude(i,j)

                ! Then calculate indices of INERIS inventory
                ii = NINT((lon*10.) - i_istart + 1)
                ji = NINT((lat*10.) - i_jstart + 1)
                
                ! Then check if INERIS and LANUV overlap completely
                IF(ineris_lanuv(ii,ji)) THEN

                   ! Then scale LANUV data based on INERIS values
                   lanuv_data(i,j,n,s) = lanuv_data(i,j,n,s)*factor

                   ! Then scale LANUV industry
                   IF(s .LE. 4) THEN
                      DO k=1,kx
                         lanuv_industry(i,k,j,n,s) = lanuv_industry(i,k,j,n,s)*factor
                      END DO
                   END IF

                END IF
             END DO
          END DO
       END DO
    END DO



    ! DEBUGGING:

    ! Calculate factors for LANUV data
    factor_lanuv = 0.d0
    gridbox_area = 1. ! LANUV data is given in per sq. km.

    ! Calculate LANUV values for each species and sector
    DO s=1,nsources
       DO n=1,nradm
          DO j=1,l_jemis
             DO i=1,l_iemis

                ! Find the INERIS grid box that covers the LANUV grid box
                ! First find lat/lon
                lon = l_longitude(i,j)
                lat = l_latitude(i,j)

                ! Then calculate indices of INERIS inventory
                ii = NINT((lon*10.) - i_istart + 1)
                ji = NINT((lat*10.) - i_jstart + 1)
                
                ! Then chekc if INERIS and LANUV overlap completely
                IF(ineris_lanuv(ii,ji)) THEN

                   ! Add to total in unit mole/hr
                   factor_lanuv(n,s) = factor_lanuv(n,s) + lanuv_data(i,j,n,s)*gridbox_area

                   ! Then add industry
                   IF(s .LE. 4) THEN
                      DO k=1,kx
                         factor_lanuv(n,s) = factor_lanuv(n,s) + lanuv_industry(i,k,j,n,s)*gridbox_area
                      END DO
                   END IF

                END IF
             END DO
          END DO
       END DO
    END DO

    ! Write new total LANUV data to screen
!    WRITE(*,*)
!    WRITE(*,*) 'NEW total LANUV data (mole/hr) in area where we have complete LANUV coverage in INERIS:'
!    DO n=1,nradm
!       WRITE(*,*) ename(n)
!       WRITE(*,*) factor_lanuv(n,:)
!    END DO

    ! END OF DEBUGGING


    ! Remove Ruhr emissions
    IF(ruhr .EQV. .false.) THEN

       lanuv_data = 0.d0
       lanuv_industry = 0.d0

    END IF


  END SUBROUTINE merge_lanuv_ineris


!!! REMOVE THIS SUBROUTINE:

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE overlap_lanuv_ineris()
    USE module_lanuv
    IMPLICIT NONE

    INTEGER                 :: i,j,il,jl,ilun
    REAL                    :: low_lon,high_lon,low_lat,high_lat,lon,lat
    LOGICAL                 :: lanuv_ineris(l_iemis,l_jemis)
    LOGICAL                 :: covered, covered_completely
    CHARACTER(len=160)      :: filename

    ilun = 51
!    ineris_lanuv    = .false.    ! .true. if the whole INERIS box is covered by LANUV
    lanuv_ineris    = .false.

    filename = './plotfiles/lanuv_ineris.dat'  ! Construct filename
    WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

    OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED')

    DO j=1,i_jemis
       DO i=1,i_iemis
          
          ! Find the boundaries of the INERIS grid boxes
          low_lon = i_longitude(i,j)-0.05
          high_lon = i_longitude(i,j)+0.05
          low_lat = i_latitude(i,j)-0.05
          high_lat = i_latitude(i,j)+0.05

          ! Then loop through all LANUV grid boxes
          covered = .false.
          covered_completely = .true.

          DO jl=1,l_jemis
!             DO il=1,l_iemis
             DO il=1,l_iemis_short
                lon = l_longitude(il,jl)
                lat = l_latitude(il,jl)

                ! If within INERIS grid box
                IF((lon .GE. low_lon) .AND. (lon .LT. high_lon) .AND. &
                     (lat .GE. low_lat) .AND. (lat .LT. high_lat)) THEN
                   covered = .true.

                   IF(.NOT. lanuv_covered(il,jl)) THEN
                      covered_completely = .false.
                   END IF
                END IF

             END DO
          END DO

          IF(covered .AND. covered_completely) THEN

             ! G� gjennom p� nytt for � fylle inn verdier i array
             DO jl=1,l_jemis
!             DO il=1,l_iemis
                DO il=1,l_iemis_short
                   lon = l_longitude(il,jl)
                   lat = l_latitude(il,jl)

                   ! If within INERIS grid box
                   IF((lon .GE. low_lon) .AND. (lon .LT. high_lon) .AND. &
                        (lat .GE. low_lat) .AND. (lat .LT. high_lat)) THEN

                      lanuv_ineris(il,jl) = .true.

                   END IF

                END DO
             END DO
          END IF
       END DO
    END DO


    ! G� gjennom p� nytt for � skrive til fil
    DO jl=1,l_jemis
       DO il=1,l_iemis_short
          IF(lanuv_ineris(il,jl)) THEN
             WRITE(ilun,'(I1)') 1
          ELSE
             WRITE(ilun,'(I1)') 0
          END IF
       END DO
    END DO

    CLOSE(ilun)

  END SUBROUTINE overlap_lanuv_ineris

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE deallocate_ineris()
  ! Deallocate arrays for INERIS emissions which are not required later on in program
  ! T. Pugh, 18/09/10
    USE emission
    IMPLICIT NONE

    DEALLOCATE(ineris_plot)
    DEALLOCATE(ineris_readdata)

  END SUBROUTINE deallocate_ineris

END MODULE module_ineris

