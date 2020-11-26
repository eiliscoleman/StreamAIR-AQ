MODULE module_retro
  USE emission

  INTEGER, PARAMETER :: r_sources = 7, r_iemis = 720, r_jemis = 360, &
       r_files = 27 ! Changed from 44->27 - removed fire emissions
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: retro_plot
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: retro_data
  REAL, ALLOCATABLE, DIMENSION(:,:)     :: r_longitude
  REAL, ALLOCATABLE, DIMENSION(:,:)     :: r_latitude
  LOGICAL, ALLOCATABLE, DIMENSION(:,:)  :: retro_ineris
  CHARACTER(len=160) :: retropath
  INTEGER            :: sector(r_sources)
  DATA sector / 1,3,2,6,7,10,11 / ! From RETRO to SNAP

  ! Make our own component type:
  TYPE component
     CHARACTER(len=80) :: name
     INTEGER           :: number
     INTEGER           :: nfields
     CHARACTER(len=10), POINTER :: fields(:)
     CHARACTER(len=80) :: file
     INTEGER           :: low_res    ! 1= we have 360x180 instead of 720x360
     REAL              :: mw         ! Molecular Weight (g/mol)
  END TYPE component

  TYPE(component), DIMENSION(r_files) :: components

CONTAINS

  SUBROUTINE allocate_retro()
  ! Allocate arrays for RETRO emissions
  ! Using allocatable arrays to avoid compiler limits on array size
  ! Note that if array is too big for available RAM on the system an out-of-memory error will be generated.
  ! T. Pugh, 20/09/10
    USE emission
    IMPLICIT NONE

    ALLOCATE(retro_plot(r_iemis,r_jemis,nchem))
    ALLOCATE(retro_data(r_iemis,r_jemis,nchem,nsources))
    ALLOCATE(r_longitude(r_iemis,r_jemis))
    ALLOCATE(r_latitude(r_iemis,r_jemis))
    ALLOCATE(retro_ineris(r_iemis,r_jemis))

  END SUBROUTINE allocate_retro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE read_retro()
  IMPLICIT NONE

  ! First we declare the variables we need:
  INTEGER            :: i,j,k,c,n,s,source
  REAL*8             :: retro_read(r_iemis,r_jemis),tempdata(360,180),value
  CHARACTER(len=160) :: filename
  CHARACTER(len=4)   :: long_nc,lat_nc,time_nc

  ! Initialize variables:
  long_nc = 'lon'
  lat_nc = 'lat'
  time_nc = 'time'
  retro_read = 0.D0
!  retropath = '../retro/RETROEMIS/'
!  retropath = '/mn/sunspot/d1/oivinho/EMIS/RETROEMIS/'
  retropath= TRIM(path)//'RETROEMIS/'  !- TP 23/08/10

  CALL readcomponents()

  ! Start with reading the variables from the RETRO inventory:
  PRINT *, ' - Reading RETRO emissions from nc files'
  PRINT *, ''

  DO n = 1, r_files      ! Read every netCDF file
     filename = TRIM(retropath)//components(n)%file  ! Construct filename

     ! Print some info first:
     WRITE(*,FMT='(A,A)') 'Name of component: ', components(n)%name
     WRITE(*,FMT='(A,I3)') 'RADM2 number:      ', components(n)%number
     WRITE(*,FMT='(A,F7.3)') 'Mole mass:         ', components(n)%mw
     WRITE(*,FMT='(A,A)') 'Filename:          ', filename

     ! We want to get data from all the sources we have chosen:
     DO c = 1, components(n)%nfields

        ! Emissions from ships have lower resolution:
        IF (components(n)%low_res .EQ. 1) THEN
           tempdata = 0.D0
           CALL readnc_2d_2(tempdata,components(n)%fields(c), &
                long_nc,lat_nc,time_nc,filename,360,180,month)
           ! Convert the 360x180 grid in to a 720x360 grid
           DO j = 1, 180
              DO i = 1, 360
                 DO k = 1, 2
                    DO s = 1, 2
                       retro_read(2*(i-1)+k,2*(j-1)+s) = tempdata(i,j)
                    END DO
                 END DO
              END DO
           END DO
        ELSE
           CALL readnc_2d_2(retro_read,components(n)%fields(c), &
                long_nc,lat_nc,time_nc,filename,r_iemis,r_jemis,month)
        END IF

        IF(components(n)%fields(c) .EQ. 'pow') THEN
           source = 1
        ELSE IF(components(n)%fields(c) .EQ. 'inc') THEN
           source = 2
        ELSE IF(components(n)%fields(c) .EQ. 'res') THEN
           source = 3
        ELSE IF(components(n)%fields(c) .EQ. 'sol') THEN
           source = 4
        ELSE IF(components(n)%fields(c) .EQ. 'tra') THEN
           source = 5
        ELSE IF(components(n)%fields(c) .EQ. 'agr') THEN
           source = 6
        ELSE
           source = 7
        END IF

        ! Convert the gas emissions from kg/(s*m^2) to mole/(km^2*h).
        retro_read(:,:) = (retro_read(:,:)*3.6*1.E12)/components(n)%mw

        ! Store the sum in an array to be used by plot_retro()
        retro_plot(:,:,components(n)%number) = retro_plot(:,:,components(n)%number) + retro_read(:,:)

        ! Add the emission data from different sources into a big array:
        retro_data(:,:,components(n)%number,sector(source)) = &
             retro_data(:,:,components(n)%number,sector(source)) + retro_read(:,:)
        retro_read = 0.D0
     END DO

     WRITE(*,*) '-------------------------------------------------------------'

  END DO ! n = 1, r_files

  WRITE(*,*) ' - Done reading from nc files'

END SUBROUTINE read_retro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE readcomponents()
  IMPLICIT NONE

  ! Local variables
  INTEGER            :: i,j,n,ilun
  CHARACTER(len=160) :: filename
  CHARACTER(len=5)   :: buffer

  filename = TRIM(path)//'components.txt'

  WRITE(*,FMT='(A,A)') '  - Reading ', filename
  WRITE(*,*)

  ilun = 12
  OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',STATUS='OLD')
  
  DO n = 1, r_files
     READ(UNIT=ilun,FMT='(A)') buffer
     READ(UNIT=ilun,FMT='(A)') components(n)%name
     READ(UNIT=ilun,FMT='(I4)') components(n)%number  ! Width specifier added to format statement to avert compiler error - TP 25/08/10
     READ(UNIT=ilun,FMT='(I4)') components(n)%nfields ! Width specifier added to format statement to avert compiler error - TP 25/08/10
     j = components(n)%nfields
     ALLOCATE(components(n)%fields(j))
     DO i = 1, j
        READ(UNIT=ilun,FMT='(A)') components(n)%fields(i)
     END DO
     READ(UNIT=ilun,FMT='(A)') components(n)%file
     READ(UNIT=ilun,FMT='(I4)') components(n)%low_res ! Width specifier added to format statement to avert compiler error - TP 25/08/10
     READ(UNIT=ilun,FMT='(F10.5)') components(n)%mw
  END DO
  CLOSE(UNIT=ilun)

END SUBROUTINE readcomponents

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE latlon_retro()
    IMPLICIT NONE

    ! We want to write lats and longs to file for the whole RETRO grid
    INTEGER            :: i,j,ilun,rstat
    REAL               :: lon,lat
    CHARACTER(len=160) :: filename

    PRINT *, ' - Write latitudes and longitudes for RETRO grid'

    ! Construct filename
    filename = './plotfiles/retro_latlon.dat'
    WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

    ! Open file for writing
    ilun = 41
    OPEN(ilun,FILE=filename,FORM='FORMATTED',STATUS='UNKNOWN',IOSTAT=rstat)

    ! Calculate longitudes
!    lon = -0.25
!    DO i = 1, r_iemis
!       lon = lon + 0.5
!       r_longitude(i,:) = lon
!    END DO

    ! Calculate longitudes
    lon = -0.25
    DO i = 1, r_iemis/2
       lon = lon + 0.5
       r_longitude(i,:) = lon
    END DO
    lon = -180.25
    DO i = r_iemis/2 + 1, r_iemis
       lon = lon + 0.5
       r_longitude(i,:) = lon
    END DO

    ! Calculate latitudes
    lat = -90.25
    DO j = 1, r_jemis
       lat = lat + 0.5
       r_latitude(:,j) = lat
    END DO

    ! Write lats and longs to file for plotting
    DO j = 1, r_jemis
       DO i = 1, r_iemis
          WRITE(ilun,FMT='(2F12.2)') r_latitude(i,j),r_longitude(i,j)
       END DO
    END DO


!    filename_nc = TRIM(retropath)//'RETRO_TNO_CO_2000.nc'  ! Construct filename
!    CALL readnc_1d_2(lons,'lon','lon','lat',filename_nc,r_iemis)
!    WRITE(*,*) lons
!    CALL readnc_1d_2(lats,'lat','lon','lat',filename_nc,r_jemis)
!    WRITE(*,*) lats

    CLOSE(ilun)

    WRITE(*,*) ' - Done writing latitudes and longitudes'

  END SUBROUTINE latlon_retro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE plot_retro()
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

    ! Unit is mole/km2/hour
    DO n = 1, nchem
       component_name = ename(n)
       filename = './plotfiles/retro_'//TRIM(smonth)//'_'//TRIM(component_name(3:))//'.dat'

       WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

       OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',IOSTAT=rstat)
       DO j=1,r_jemis
          DO i=1,r_iemis
             WRITE(UNIT=ilun,FMT='(F16.4)',IOSTAT=rstat) retro_plot(i,j,n)
          END DO
       END DO
       CLOSE(UNIT=ilun)
    END DO

  END SUBROUTINE plot_retro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE overlap_retro_ineris()
    USE module_ineris
    IMPLICIT NONE

    INTEGER                 :: i,j,ii,ji,ilun
    REAL                    :: low_lon,high_lon,low_lat,high_lat,lon,lat
    LOGICAL                 :: covered, covered_partly
    CHARACTER(len=160)      :: filename

    ilun = 51
    retro_ineris = .false.    ! .true. if some of the RETRO box is covered by INERIS

    filename = './plotfiles/retro_ineris.dat'  ! Construct filename
    WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

    OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED')

    DO j=1,r_jemis
       DO i=1,r_iemis
          
          ! Find the boundaries of the RETRO grid boxes
          low_lon = r_longitude(i,j)-0.25
          high_lon = r_longitude(i,j)+0.25
          low_lat = r_latitude(i,j)-0.25
          high_lat = r_latitude(i,j)+0.25

          ! Then loop through all INERIS grid boxes
          covered = .false.
          covered_partly = .false.

          DO ji=1,i_jemis
             DO ii=1,i_iemis
                lon = i_longitude(ii,ji)
                lat = i_latitude(ii,ji)

                ! If within RETRO grid box
                IF((lon .GE. low_lon) .AND. (lon .LT. high_lon) .AND. &
                     (lat .GE. low_lat) .AND. (lat .LT. high_lat)) THEN
                   covered = .true.

                   IF(ineris_covered(ii,ji)) THEN
                      covered_partly = .true.
                   END IF
                END IF

             END DO
          END DO

          IF(covered .AND. covered_partly) THEN
             retro_ineris(i,j) = .true.
             WRITE(ilun,'(I1)') 1 ! Width specifier added to format statement to avert compiler error - TP 25/08/10
          ELSE
             retro_ineris(i,j) = .false.
             WRITE(ilun,'(I1)') 0 ! Width specifier added to format statement to avert compiler error - TP 25/08/10
          END IF
       END DO
    END DO

    CLOSE(ilun)

  END SUBROUTINE overlap_retro_ineris

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_overlap_retro_ineris()
    USE module_ineris
    IMPLICIT NONE

    INTEGER                 :: i,j,icovered,ilun
    CHARACTER(len=160)      :: filename

    ilun = 51
    retro_ineris  = .false.  ! .true. if some of the RETRO box is covered by INERIS

    filename = './plotfiles/retro_ineris.dat'  ! Construct filename
    WRITE(*,FMT='(A,A)') '   - Reading from file: ', TRIM(filename)

    OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED')

    DO j=1,r_jemis
       DO i=1,r_iemis
          READ(ilun,*) icovered ! Format statement changed from I to * to avert compiler error - TP 25/08/10
          IF(icovered .EQ. 1) THEN
             retro_ineris(i,j) = .true.
          ELSE
             retro_ineris(i,j) = .false.
          END IF
       END DO
    END DO

    CLOSE(ilun)

  END SUBROUTINE read_overlap_retro_ineris

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE deallocate_retro()
  ! Deallocate arrays for RETRO emissions which are not required later on in program
  ! T. Pugh, 18/09/10
    USE emission
    IMPLICIT NONE

    DEALLOCATE(retro_plot)

  END SUBROUTINE deallocate_retro

END MODULE module_retro
