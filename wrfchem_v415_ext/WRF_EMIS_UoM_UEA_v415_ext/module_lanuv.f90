MODULE module_lanuv
  USE emission

  INTEGER, PARAMETER :: l_iemis = 385, l_jemis = 250, l_comps = 3, &
       l_istart = 2490, l_iend = 2738, l_jstart = 5576, l_jend = 5825,&
       l_iemis_diff = 654, l_iemis_short = 249, l_istart2 = 3393, &
       l_iend2 = 3528
  ! Have reduced the array sizes because of memory constraints

!  INTEGER, PARAMETER :: l_iemis = 249, l_jemis = 250, l_comps = 3, &
!       l_istart = 2490, l_iend = 2738, l_jstart = 5576, l_jend = 5825,&
!       l_iemis2 = 136, l_istart2 = 3393, l_iend2 = 3528
!  INTEGER, PARAMETER :: l_iemis = 1039, l_jemis = 250, l_comps = 3, &
!       l_istart = 2490, l_iend = 3528, l_jstart = 5576, l_jend = 5825
  REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: lanuv_plot
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)   :: lanuv_readdata
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)   :: lanuv_data
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)   :: lanuv_readindustry
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: lanuv_industry ! Only 4 SNAP sources for industry
  REAL, ALLOCATABLE, DIMENSION(:,:)       :: l_longitude
  REAL, ALLOCATABLE, DIMENSION(:,:)       :: l_latitude
  REAL                                    :: lanuv_ind_fac(4,l_comps)
  LOGICAL, ALLOCATABLE, DIMENSION(:,:)    :: lanuv_covered
  CHARACTER (len=3), DIMENSION(l_comps)   :: components_lanuv

  DATA components_lanuv /'CO','NOX','VOC'/
!  DATA lanuv_ind_fac / &  ! Temporary factors
!       0.333, 0., 0.333, 0.333, &
!       0.333, 0., 0.333, 0.333, &
!       0.333, 0., 0.333, 0.333 &
!       /
  DATA lanuv_ind_fac / &  ! These factors are derived from INERIS data for 2003
       0.0550274908, 0., 0.4406814255, 0.5042910837, &  ! CO
       0.6425936281, 0., 0.3296879367, 0.0277184351, &  ! NOX
       0.0970776612, 0., 0.0800661209, 0.8228562179 &   ! VOC
       /

CONTAINS

  SUBROUTINE allocate_lanuv()
  ! Allocate arrays for LANUV emissions
  ! Using allocatable arrays to avoid compiler limits on array size
  ! Note that if array is too big for available RAM on the system an out-of-memory error will be generated.
  ! T. Pugh, 20/09/10
    USE emission
    IMPLICIT NONE

    ALLOCATE(lanuv_plot(l_iemis,l_jemis,l_comps))
    ALLOCATE(lanuv_readdata(l_iemis,l_jemis,l_comps,nsources))
    ALLOCATE(lanuv_data(l_iemis,l_jemis,nchem,nsources))
    ALLOCATE(lanuv_readindustry(l_iemis,max_k,l_jemis,l_comps))
    ALLOCATE(lanuv_industry(l_iemis,max_k,l_jemis,nchem,4))
    ALLOCATE(l_longitude(l_iemis,l_jemis))
    ALLOCATE(l_latitude(l_iemis,l_jemis))
    ALLOCATE(lanuv_covered(l_iemis,l_jemis))

  END SUBROUTINE allocate_lanuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_lanuv()
    IMPLICIT NONE

    INTEGER            :: i,j,k,n,x,y,c,s,ilun,rstat,ibuffer,height, &
         vert_layer
    REAL               :: lat,long,source,sum_sources,x_temp,y_temp
    REAL               :: nox_road,nox_air,nox_offroad,nox_rail,nox_ship, &
         co_road,co_air,co_offroad,co_rail,co_ship,voc_road,voc_air, &
         voc_offroad,voc_rail,voc_ship,nox_total,co_total,voc_total
    CHARACTER(len=160) :: filename
    CHARACTER(len=300) :: buffer
    CHARACTER(len=3)   :: component_name
    INTEGER, DIMENSION(l_comps)          :: nlines_industry
    DATA nlines_industry / 5381,5482,11186 /

    ! Now we will read the variables from the LANUV inventory
    PRINT *, ' - Reading LANUV Ruhr emissions from txt files'

    ilun = 21
    lanuv_covered  = .false.
    lanuv_plot     = 0.D0
    lanuv_readdata     = 0.D0
    lanuv_readindustry = 0.D0


    ! *** Heating ***
    DO n = 1, l_comps      ! there are 3 components (CO, NOX, VOC)

       ! Construct filename
       component_name = components_lanuv(n)
       filename = TRIM(path)//'lanuv/heating_'//TRIM(component_name)//'_2004.txt'

       WRITE(*,FMT='(A,A)') '   - Reading from file: ', TRIM(filename)

       ! Read from the '.txt' file. The unit is kg/km2/year
       OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',STATUS='OLD',IOSTAT=rstat)
       READ(ilun,FMT='(A)') buffer

       ! Read every point in the inventory:
       DO c = 1,11548
          READ(ilun,FMT=*) x,y,source

          i = (x/1000) - l_istart + 1
          j = (y/1000) - l_jstart + 1

          ! Skip empty array fields to reduce memory requirements
          IF(i .GT. 900) THEN
             i = i - l_iemis_diff
          END IF

          lanuv_covered(i,j) = .true.

          ! Convert from kg/year to t/year
          source = source/1000.

          ! Store the sum in an array to be used by plot_lanuv()
          lanuv_plot(i,j,n) = lanuv_plot(i,j,n) + source

          ! Store the sources in an array
          lanuv_readdata(i,j,n,2) = lanuv_readdata(i,j,n,2) + source   ! For now, store all in SNAP 2

       END DO
       CLOSE(ilun)

    END DO ! n = 1, n_lanuv


    ! *** Industry ***
    DO n = 1, l_comps      ! there are 3 components (CO, NOX, VOC)

       ! Construct filename
       component_name = components_lanuv(n)
       filename = TRIM(path)//'lanuv/industry_'//TRIM(component_name)//'_2004.txt'
       WRITE(*,FMT='(A,A)') '   - Reading from file: ', TRIM(filename)

       ! Read from the '.txt' file. The unit is kg/year
       OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',STATUS='OLD',IOSTAT=rstat)
       READ(ilun,FMT='(A)') buffer

       ! Read every point in the inventory:
       DO c = 1,nlines_industry(n)
          READ(ilun,FMT=*) ibuffer,x_temp,y_temp,height,source     ! SJEKK SNAP-KODER!!!!

          i = NINT((x_temp/1000.) - l_istart + 1)
          j = NINT((y_temp/1000.) - l_jstart + 1)

!          IF(i .GE. 220 .AND. i .LT. 904) THEN
!             WRITE(*,*) 'i = ', i
!          END IF

          ! Skip empty array fields to reduce memory requirements
          IF(i .GT. 900) THEN
!             WRITE(*,*) 'i before = ', i
             i = i - l_iemis_diff
!             WRITE(*,*) 'i after = ', i
          END IF

          lanuv_covered(i,j) = .true.

          ! Convert from kg/year to t/year
          source = source/1000.

          ! Store the sum in an array to be used by plot_lanuv()
          lanuv_plot(i,j,n) = lanuv_plot(i,j,n) + source

          ! For now, we just put the emission in the WRF vertical 
          ! layer that is closest to the stack height
          IF(height .GE. heights(kx+1)) THEN
             vert_layer = kx
          ELSE
             DO k = 1, kx
                IF(height .GE. heights(k) .AND. height .LT. heights(k+1)) THEN
                   vert_layer = k
                END IF
             END DO
          END IF

          ! Store the industry emissions in a separate array because of height information
          lanuv_readindustry(i,vert_layer,j,n) = lanuv_readindustry(i,vert_layer,j,n) + source   ! SJEKK SNAP-KODER!!!

       END DO
       CLOSE(ilun)

    END DO ! n = 1, n_lanuv


    ! *** Transport ***

    filename = TRIM(path)//'lanuv/transport_2000_2007.txt'   ! Construct filename
    WRITE(*,FMT='(A,A)') '   - Reading from file: ', TRIM(filename)

    ! Read from the '.txt' file. The unit is kg/km2/year
    OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',STATUS='OLD',IOSTAT=rstat)
    READ(ilun,FMT='(A)') buffer

    ! Read every point in the inventory:
    DO c = 1,35072
       READ(ilun,FMT=*) ibuffer,x,y,nox_road,nox_air,nox_offroad,nox_rail, &
            nox_ship,co_road,co_air,co_offroad,co_rail,co_ship,voc_road,voc_air, &
            voc_offroad,voc_rail,voc_ship,nox_total,co_total,voc_total

       i = x - l_istart + 1
       j = y - l_jstart + 1

       ! Skip empty array fields to reduce memory requirements
       IF(i .GT. 900) THEN
          i = i - l_iemis_diff
       END IF

       lanuv_covered(i,j) = .true.

       ! Convert from kg/year to t/year
       co_total = co_total/1000.
       nox_total = nox_total/1000.
       voc_total = voc_total/1000.
       co_road = co_road/1000.
       nox_road = nox_road/1000.
       voc_road = voc_road/1000.
       co_air = co_air/1000.
       nox_air = nox_air/1000.
       voc_air = voc_air/1000.
       co_offroad = co_offroad/1000.
       nox_offroad = nox_offroad/1000.
       voc_offroad = voc_offroad/1000.
       co_rail = co_rail/1000.
       nox_rail = nox_rail/1000.
       voc_rail = voc_rail/1000.
       co_ship = co_ship/1000.
       nox_ship = nox_ship/1000.
       voc_ship = voc_ship/1000.

       ! Store the sums in an array to be used by plot_lanuv()
       lanuv_plot(i,j,1) = lanuv_plot(i,j,1) + co_total
       lanuv_plot(i,j,2) = lanuv_plot(i,j,2) + nox_total
       lanuv_plot(i,j,3) = lanuv_plot(i,j,3) + voc_total

       ! Store the sources in an array
       lanuv_readdata(i,j,1,7) = lanuv_readdata(i,j,1,7) + co_road
       lanuv_readdata(i,j,1,8) = lanuv_readdata(i,j,1,8) + co_air+co_offroad+co_rail+co_ship
       lanuv_readdata(i,j,2,7) = lanuv_readdata(i,j,2,7) + nox_road
       lanuv_readdata(i,j,2,8) = lanuv_readdata(i,j,2,8) + nox_air+nox_offroad+nox_rail+nox_ship
       lanuv_readdata(i,j,3,7) = lanuv_readdata(i,j,3,7) + voc_road
       lanuv_readdata(i,j,3,8) = lanuv_readdata(i,j,3,8) + voc_air+voc_offroad+voc_rail+voc_ship

    END DO
    CLOSE(ilun)


    ! Short summary
    DO n = 1, l_comps      ! there are 3 components (CO, NOX, VOC)
       sum_sources = 0.D0
       DO j = 1, l_jemis
          DO i = 1, l_iemis
             sum_sources = sum_sources + lanuv_plot(i,j,n)
          END DO
       END DO
       WRITE(*,FMT='(3A,F15.3,A)') 'LANUV: Sum of ',components_lanuv(n),' = ',sum_sources,' t/yr'
    END DO

    WRITE(*,*) ' - Done reading from txt files'

  END SUBROUTINE read_lanuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE latlon_lanuv()
    IMPLICIT NONE

    ! We want to calculate lats and longs for the whole LANUV-grid
    INTEGER            :: i,j,ilun,rstat
    REAL               :: phi,lam,gx,gy
    CHARACTER(len=160) :: filename

    PRINT *, ' - Calculate latitudes and longitudes for LANUV grid'

    ! Construct filename
    filename = './plotfiles/lanuv_latlon.dat'
    WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

    ! Open file for writing
    ilun = 41
    OPEN(ilun,FILE=filename,FORM='FORMATTED',STATUS='UNKNOWN',IOSTAT=rstat)

    phi = 0.
    lam = 0.

    ! Loop through LANUV grid
    DO j = l_jstart, l_jend
       DO i = l_istart, l_iend

          gx = (i*1000.)+500.
          gy = (j*1000.)+500.

          CALL G2PHILA(gx,gy,phi,lam)   ! Convert from Gauss-Krueger to lat and long
          
          phi=phi*180./3.14159265       ! Convert from radians to degrees
          lam=lam*180./3.14159265

          WRITE(ilun,FMT='(2F12.6)') phi,lam     ! Write to file for plotting

          l_longitude(i-l_istart+1,j-l_jstart+1) = lam
          l_latitude(i-l_istart+1,j-l_jstart+1) = phi

       END DO

       ! Second part of grid
       DO i = l_istart2, l_iend2

          gx = (i*1000.)+500.
          gy = (j*1000.)+500.

          CALL G2PHILA(gx,gy,phi,lam)   ! Convert from Gauss-Krueger to lat and long
          
          phi=phi*180./3.14159265       ! Convert from radians to degrees
          lam=lam*180./3.14159265

          WRITE(ilun,FMT='(2F12.6)') phi,lam     ! Write to file for plotting

          l_longitude(i-l_istart+1-l_iemis_diff,j-l_jstart+1) = lam
          l_latitude(i-l_istart+1-l_iemis_diff,j-l_jstart+1) = phi

       END DO

    END DO

    CLOSE(ilun)

    WRITE(*,*) ' - Done calculating latitudes and longitudes'

  END SUBROUTINE latlon_lanuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE plot_lanuv()
    IMPLICIT NONE
    CHARACTER (len=80)          :: component_name,filename
    INTEGER                     :: i,j,n,ilun,rstat

    ilun = 51

    ! Units are in kg/year if we run this routine before convert_units()
    DO n = 1, l_comps
       component_name = components_lanuv(n)
       filename = './plotfiles/lanuv_'//TRIM(component_name)//'.dat'

       WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

       OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',IOSTAT=rstat)
       DO j=1,l_jemis
          DO i=1,l_iemis
             WRITE(UNIT=ilun,FMT='(F16.4)',IOSTAT=rstat) lanuv_plot(i,j,n)
          END DO
       END DO
       CLOSE(UNIT=ilun)
    END DO

  END SUBROUTINE plot_lanuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE convert_lanuv()
    IMPLICIT NONE

    ! Convert the gas emissions from tons/yr to mole/(h*km^2):
    ! Have to multiply by 1e6 to convert from tons to g.
    ! Have to divide by 365*24 to convert from pr. year to pr. hour.
    ! Have to divide by gridbox_area to get pr. square km.
    ! Have to divide by molar mass to convert from g to mole.
    REAL        :: constant,gridbox_area,mw(2)
    INTEGER     :: lanuv_to_radm(2)
    INTEGER     :: i, j, n, s
    DATA mw / 28.0101, 46.0055 / ! CO, NOX (NO2)
    DATA lanuv_to_radm / 11, 2 / ! CO, NOx

    lanuv_data = 0.d0
    lanuv_industry = 0.d0
    gridbox_area = 1.  ! LANUV data is per sq. km
    constant = 1000000./(365.*24.)

    DO n = 1,l_comps-1 ! Skip the VOCs for now
       DO j=1,l_jemis
          DO i=1,l_iemis
             lanuv_data(i,j,lanuv_to_radm(n),:) = &
                  (constant*lanuv_readdata(i,j,n,:))/(gridbox_area*mw(n))

             ! Industry is in another array because of height information
             DO s = 1, 4
                lanuv_industry(i,:,j,lanuv_to_radm(n),s) = &
                     (constant*lanuv_readindustry(i,:,j,n)*lanuv_ind_fac(s,n))/(gridbox_area*mw(n))
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE convert_lanuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE vocs_lanuv()
    USE voc_factors
    IMPLICIT NONE

    ! Divide VOCs into RADM-components and convert the units

    REAL              :: mw(nvocs)
    REAL              :: constant,gridbox_area,value,sum_vocs
    INTEGER           :: i,j,k,n,s,radm_nr(nvocs)

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
    gridbox_area = 1.  ! LANUV data is per sq. km
    constant = 1000000./(365.*24.)

    WRITE(*,*) ' - Dividing VOC species based on sources'

    ! All the (47*9 = 423) factors are in this file
    !INCLUDE 'input/voc_factors.f90'

    DO s=1,nvocsrc
       DO n=1,nvocs
          DO j=1,l_jemis
             DO i=1,l_iemis
                ! Divide VOCs
                value = lanuv_readdata(i,j,l_comps,s)*factor(n,s)
                sum_vocs = sum_vocs + value

                ! Convert units
                value = (constant*value)/(gridbox_area*mw(n))
                lanuv_data(i,j,radm_nr(n),s) = &
                     lanuv_data(i,j,radm_nr(n),s) + value


                ! Industry must be done seperately:

                DO k=1,kx
                   IF(s .LE. 4) THEN
                      ! Divide VOC category into SNAP sources and then into VOC species
                      value = lanuv_readindustry(i,k,j,l_comps)*lanuv_ind_fac(s,l_comps)*factor(n,s)
                      sum_vocs = sum_vocs + value

                      ! Convert units
                      value = (constant*value)/(gridbox_area*mw(n))
                      lanuv_industry(i,k,j,radm_nr(n),s) = &
                           lanuv_industry(i,k,j,radm_nr(n),s) + value
                   END IF
                END DO

             END DO
          END DO
       END DO
    END DO

    WRITE(*,*) '  - Sum VOC before speciation: ', &
 SUM(lanuv_readdata(:,:,l_comps,:))+SUM(lanuv_readindustry(:,:,:,l_comps)), ' tons/yr'
    WRITE(*,*) '  - Sum VOC after speciation:  ', sum_vocs, ' tons/yr'

    WRITE(*,*) ' - Done dividing VOC species'

  END SUBROUTINE vocs_lanuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Find nearest LANUV-point from given latitude and longitude
  SUBROUTINE find_lanuv_index(lon,lat,istart,iend,jstart,jend,i_save,j_save)
    IMPLICIT NONE

    INTEGER :: ifn, jfn, i_save, j_save
    INTEGER :: istart, iend, jstart, jend
    REAL    :: lon, lat, smallest_diff, diff_x, diff_y, diff_sum

    smallest_diff = 10.
    DO jfn = jstart, jend
       DO ifn = istart, iend
          diff_x = lon-l_longitude(ifn,jfn)
          diff_y = lat-l_latitude(ifn,jfn)
          diff_sum = SQRT(diff_x**2+diff_y**2)
          IF(diff_sum .LT. smallest_diff) THEN
             smallest_diff = diff_sum
             ! Save the coordinates of the closest LANUV point
             i_save = ifn
             j_save = jfn
          END IF
       END DO
    END DO
  END SUBROUTINE find_lanuv_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Find nearest LANUV-point from given latitude and longitude
  SUBROUTINE find_lanuv_index_dist(lon,lat,istart,iend,jstart,jend,&
       i_save,j_save,delta_x,delta_y)
    IMPLICIT NONE

    INTEGER :: ifn, jfn, i_save, j_save
    INTEGER :: istart, iend, jstart, jend
    REAL    :: lon, lat, lat_save, smallest_diff, diff_x, diff_y, &
         save_diff_x, save_diff_y, diff_sum, earth_radius, pi, &
         delta_x, delta_y

    earth_radius = 6371.     ! Earth's radius in km
    pi = 3.14159265         ! Pi

    smallest_diff = 10.
    DO jfn = jstart, jend
       DO ifn = istart, iend
          diff_x = lon-l_longitude(ifn,jfn)
          diff_y = lat-l_latitude(ifn,jfn)
          diff_sum = SQRT(diff_x**2+diff_y**2)
          IF(diff_sum .LT. smallest_diff) THEN
             smallest_diff = diff_sum
             ! Save the coordinates of the closest LANUV point
             i_save = ifn
             j_save = jfn
             save_diff_x = (diff_x*pi)/180 ! x-distance in radians
             save_diff_y = (diff_y*pi)/180 ! y-distance in radians
             lat_save = (l_latitude(ifn,jfn)*pi)/180 ! latitude in radians
          END IF
       END DO
    END DO

    ! Calculate distance in km
    delta_x = earth_radius*save_diff_x*COS(lat_save)
    delta_y = earth_radius*save_diff_y

!    distance = SQRT(delta_x**2+delta_y**2)

  END SUBROUTINE find_lanuv_index_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE deallocate_lanuv()
  ! Deallocate arrays for LANUV emissions which are not required later on in program
  ! T. Pugh, 18/09/10
    USE emission
    IMPLICIT NONE

    DEALLOCATE(lanuv_plot)
    DEALLOCATE(lanuv_readdata)
    DEALLOCATE(lanuv_readindustry)

  END SUBROUTINE deallocate_lanuv

END MODULE module_lanuv
