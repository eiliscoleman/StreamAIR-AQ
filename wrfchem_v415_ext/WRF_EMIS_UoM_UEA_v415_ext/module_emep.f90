MODULE module_emep
  USE emission

!  INTEGER, PARAMETER :: e_iemis = 601, e_jemis = 451, e_comps = 4, &
!       e_istart = -200, e_iend = 400, e_jstart = 300, e_jend = 750
  INTEGER, PARAMETER :: e_iemis = 601, e_jemis = 471, e_comps = 4, &
       e_istart = -200, e_iend = 400, e_jstart = 280, e_jend = 750
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: emep_plot
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: emep_readdata
! REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: emep_data
  REAL, ALLOCATABLE, DIMENSION(:,:)     :: e_longitude
  REAL, ALLOCATABLE, DIMENSION(:,:)     :: e_latitude
  LOGICAL, ALLOCATABLE, DIMENSION(:,:)  :: emep_covered
!  LOGICAL, ALLOCATABLE, DIMENSION(:,:)  :: emep_lanuv
!  LOGICAL, ALLOCATABLE, DIMENSION(:,:)  :: emep_naei
  CHARACTER(len=4)   :: syear
  CHARACTER(len=5), DIMENSION(e_comps):: components_emep

  DATA components_emep /'CO','NOx','SOx','NMVOC'/ ! NMVOC must always be last
!  DATA components_emep /'CO','NH3','NOx','PM25','PMco','SOx','NMVOC'/

CONTAINS

  SUBROUTINE allocate_emep()
  ! Allocate arrays for EMEP emissions
  ! Using allocatable arrays to avoid compiler limits on array size
  ! Note that if array is too big for available RAM on the system an out-of-memory error will be generated.
  ! T. Pugh, 20/09/10
    USE emission
    IMPLICIT NONE

    ALLOCATE(emep_plot(e_iemis,e_jemis,e_comps))
    ALLOCATE(emep_readdata(e_iemis,e_jemis,e_comps,nsources))
!    ALLOCATE(emep_data(e_iemis,e_jemis,nchem,nsources))
    ALLOCATE(e_longitude(e_iemis,e_jemis))
    ALLOCATE(e_latitude(e_iemis,e_jemis))
    ALLOCATE(emep_covered(e_iemis,e_jemis))
!    ALLOCATE(emep_lanuv(e_iemis,e_jemis))
!    ALLOCATE(emep_naei(e_iemis,e_jemis))

  END SUBROUTINE allocate_emep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_emep()
    USE emission
    IMPLICIT NONE

    INTEGER            :: i,j,n,c,s,ilun,rstat,country
    REAL               :: lat,lon,source(nsources),sum_sources
    CHARACTER(len=160) :: filename, buffer
    CHARACTER(len=5)   :: component_name

!!! TEMP:
    INTEGER :: min_i,max_i,min_j,max_j
    min_i = 1000
    max_i = 0
    min_j = 1000
    max_j = 0
!!!

    ! Now we will read the variables from the EMEP inventory
    PRINT *, ' - Reading EMEP emissions from txt files'

    ilun = 21
    emep_covered  = .false.
    emep_plot     = 0.D0
    emep_readdata = 0.D0

    DO n = 1, e_comps      ! there are 4 components (CO, NMVOC, NOX, SOx)

       ! Construct filename
       WRITE(syear,'(I4)') year
       component_name = components_emep(n)
       filename = TRIM(path)//'EMEP/'//syear//'/grid'//TRIM(component_name)//'_emis'//syear//'_10km_latlon'

       WRITE(*,FMT='(A,A)') '   - Reading from file: ', TRIM(filename)

       ! Read from the text file. The unit is Mg/(0.1deg x 0.1deg)/year
       OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',STATUS='OLD',IOSTAT=rstat)

       ! Skip the first five lines
       READ(ilun,FMT='(A)') buffer
!       WRITE(*,*) buffer
       READ(ilun,FMT='(A)') buffer
!       WRITE(*,*) buffer
       READ(ilun,FMT='(A)') buffer
!       WRITE(*,*) buffer
       READ(ilun,FMT='(A)') buffer
!       WRITE(*,*) buffer
       READ(ilun,FMT='(A)') buffer
!       WRITE(*,*) buffer

       ! Read every point in the inventory:
       DO WHILE (.true.)
!       DO c = 1,nlines(n)
          READ(ilun,FMT=*,END=999) country,lat,lon,source(1), & 
source(2),source(3),source(4),source(5),source(6),source(7), & 
source(8),source(9),source(10),source(11)

          ! Debug:
!          WRITE(*,*) country,lat,lon,source
          !

          ! Calculate index
          i = NINT((lon*10.) - e_istart + 1)
          j = NINT((lat*10.) - e_jstart + 1)

!!! TEMP:
          IF(i .LT. min_i) THEN
             min_i = i
          END IF
          IF(i .GT. max_i) THEN
             max_i = i
          END IF
          IF(j .LT. min_j) THEN
             min_j = j
          END IF
          IF(j .GT. max_j) THEN
             max_j = j
          END IF
!!!

          emep_covered(i,j) = .true.

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


! SJEKK OM SKIPSUTSLIPPENE ER DOBLET OGSÅ I EMEP, ELLER BARE I INERIS!!!!!!!!!!!!!!!!!!!!

          ! BUG FIX: Ship emissions are doubled in 2006. Remove half:
!          IF(year .EQ. 2006 .AND. country .GT. 300 .AND. country .LT. 350) THEN
!             source(:) = 0.D0
!             sum_sources = 0.
!          END IF

          ! Store the sum in an array to be used by plot_emep()
          emep_plot(i,j,n) = emep_plot(i,j,n) + sum_sources

          ! Store the sources in an array
          DO s=1,nsources
             emep_readdata(i,j,n,s) = emep_readdata(i,j,n,s) + source(s)
          END DO

       END DO
999    CONTINUE
       CLOSE(ilun)

    END DO ! n = 1, e_comps

!!! TEMP:
    WRITE(*,*) 'EMEP boundary information:'
    WRITE(*,*) 'i: ', min_i, max_i
    WRITE(*,*) 'j: ', min_j, max_j
!!!


    ! Short summary
    DO n = 1, e_comps
       sum_sources = 0.D0
       DO j = 1, e_jemis
          DO i = 1, e_iemis
             sum_sources = sum_sources + emep_plot(i,j,n)
          END DO
       END DO
       WRITE(*,FMT='(3A,F15.3,A)') 'EMEP: Sum of ',components_emep(n),' = ',sum_sources,' t/yr'
    END DO

    WRITE(*,*) ' - Done reading from txt files'

  END SUBROUTINE read_emep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE latlon_emep()
    IMPLICIT NONE

    ! We want to write lats and longs to file for the whole EMEP grid
    INTEGER            :: i,j,ilun,rstat
    REAL               :: lon,lat
    CHARACTER(len=160) :: filename

    PRINT *, ' - Write latitudes and longitudes for EMEP grid'

    ! Construct filename
    filename = './plotfiles/emep_latlon.dat'
    WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

    ! Open file for writing
    ilun = 41
    OPEN(ilun,FILE=filename,FORM='FORMATTED',STATUS='UNKNOWN',IOSTAT=rstat)

    ! Loop through EMEP grid
    DO j = e_jstart, e_jend
       DO i = e_istart, e_iend
          lon = i/10.
          lat = j/10.

          WRITE(ilun,FMT='(2F12.1)') j/10.,i/10.  ! Write to file for plotting
          e_longitude(i-e_istart+1,j-e_jstart+1) = lon
          e_latitude(i-e_istart+1,j-e_jstart+1) = lat
       END DO
    END DO

    CLOSE(ilun)

    WRITE(*,*) ' - Done writing latitudes and longitudes'

  END SUBROUTINE latlon_emep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE plot_emep()
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
    DO n = 1, e_comps
       component_name = components_emep(n)

       IF(scaling_season) THEN
          filename = './plotfiles/emep_'//syear//'_'//TRIM(smonth)//'_'//TRIM(component_name)//'.dat'
       ELSE
          filename = './plotfiles/emep_'//syear//'_'//TRIM(component_name)//'.dat'
       END IF

       WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

       OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',IOSTAT=rstat)
       DO j=1,e_jemis
          DO i=1,e_iemis
             WRITE(UNIT=ilun,FMT='(F16.4)',IOSTAT=rstat) emep_plot(i,j,n)
          END DO
       END DO
       CLOSE(UNIT=ilun)
    END DO

  END SUBROUTINE plot_emep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE deallocate_emep()
  ! Deallocate arrays for EMEP emissions which are not required later on in program
  ! T. Pugh, 18/09/10
    USE emission
    IMPLICIT NONE

    DEALLOCATE(emep_plot)
    DEALLOCATE(emep_readdata)

  END SUBROUTINE deallocate_emep

END MODULE module_emep
