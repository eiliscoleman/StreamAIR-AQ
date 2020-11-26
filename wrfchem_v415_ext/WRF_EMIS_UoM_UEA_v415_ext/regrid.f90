  ! Read latitudes and longitudes from geo_em* or wrfinput* file
  SUBROUTINE get_wrflatlon()
    USE emission
    IMPLICIT NONE
    
    INTEGER            :: i,j
    REAL*8             :: temp_read(ix,jx)
    REAL               :: wrf_i,wrf_j,lat,lon
    CHARACTER(len=160) :: geo_em_file, wrfinput_file
    
    ! Read from geo_em file
    IF(read_geofile) THEN
       geo_em_file = 'geo_em.d'//id_nr(d)//'.nc' ! Create geo_em filename

       WRITE(*,*) ' - Will read latitudes and longitudes from geo_em file: '
       WRITE(*,FMT='(A,A)') '   -  ',geo_em_file

       CALL readnc_2d_2(temp_read,'XLAT_M','west_east','south_north','Time', &
            geo_em_file,ix,jx,1)
       wrf_lats = temp_read
       temp_read = 0.D0
       CALL readnc_2d_2(temp_read,'XLONG_M','west_east','south_north','Time', &
            geo_em_file,ix,jx,1)
       wrf_lons = temp_read
       temp_read = 0.D0

    ! Read from wrfinput file
    ELSEIF(read_wrfinput) THEN
       wrfinput_file = 'wrfinput_d'//id_nr(d) ! Create wrfinput filename

       WRITE(*,*) ' - Will read latitudes and longitudes from wrfinput file: '
       WRITE(*,FMT='(A,A)') '   -  ',wrfinput_file

       CALL readnc_2d_2(temp_read,'XLAT','west_east','south_north','Time', &
            wrfinput_file,ix,jx,1)
       wrf_lats = temp_read
       temp_read = 0.D0
       CALL readnc_2d_2(temp_read,'XLONG','west_east','south_north','Time', &
            wrfinput_file,ix,jx,1)
       wrf_lons = temp_read
       temp_read = 0.D0

write(6,*)'Back from reading wrfinput'

    ! Calculate lats and longs (only approximate values)
    ELSE
       WRITE(*,*) ' - Will calculate latitudes and longitudes for wrf domain'
       DO j = 1, jx
          DO i = 1, ix
             wrf_i = FLOAT(i)+0.5
             wrf_j = FLOAT(j)+0.5
             CALL MAPCF(wrf_i,wrf_j,lat,lon)
             wrf_lats(i,j) = lat
             wrf_lons(i,j) = lon
          END DO
       END DO

    END IF

  END SUBROUTINE get_wrflatlon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE wrf_grid()
    USE emission
    USE module_retro
    USE module_tno
    USE module_ineris
    USE module_naei
    USE module_lanuv
    IMPLICIT NONE

    INTEGER :: i,j,ir,jr,ii,ji,in_big,jn_big,in,jn,i_s,j_s,n,k,nboxes, &
         il_big,jl_big,il,jl,it_big,jt_big,it,jt
    REAL    :: wrf_i,wrf_j,small_wrf_i,small_wrf_j,lat,lon, &
         lat_temp,lon_temp,lat_mid,lon_mid,lat_small,lon_small, &
         small_interval,start_value
    REAL    :: ineris_lat_start,ineris_lat_end, &
         ineris_lon_start,ineris_lon_end,naei_lat_start,&
         naei_lat_end,naei_lon_start,naei_lon_end, &
         lanuv_lat_start,lanuv_lat_end,lanuv_lon_start,lanuv_lon_end, &
         tno_lat_start,tno_lat_end,tno_lon_start,tno_lon_end
    REAL    :: mean_value(max_k,nchem,nsources), value(max_k,nchem,nsources)

write(6,*)'In wrf_grid'

    ! Define approximate lat/lon intervals of emission inventories
    IF(ineris(d)) THEN
    ineris_lat_start = MINVAL(i_latitude) - 1
    ineris_lat_end   = MAXVAL(i_latitude) + 1
    ineris_lon_start = MINVAL(i_longitude) - 1
    ineris_lon_end   = MAXVAL(i_longitude) + 1
    ENDIF
    IF(naei(d)) THEN
    naei_lat_start   = MINVAL(n_latitude) - 1
    naei_lat_end     = MAXVAL(n_latitude) + 1
    naei_lon_start   = MINVAL(n_longitude) - 1
    naei_lon_end     = MAXVAL(n_longitude) + 1
    ENDIF
    IF(lanuv(d)) THEN
    lanuv_lat_start  = MINVAL(l_latitude) - 1
    lanuv_lat_end    = MAXVAL(l_latitude) + 1
    lanuv_lon_start  = MINVAL(l_longitude) - 1
    lanuv_lon_end    = MAXVAL(l_longitude) + 1
    ENDIF
    IF(retrotno(d)) THEN
    tno_lat_start  = MINVAL(t_latitude) - 1
    tno_lat_end    = MAXVAL(t_latitude) + 1
    tno_lon_start  = MINVAL(t_longitude) - 1
    tno_lon_end    = MAXVAL(t_longitude) + 1
    ENDIF

write(6,*)'A bit futher in wrf_grid'

    WRITE(*,*) 'INERIS: ',ineris_lat_start,ineris_lat_end, &
         ineris_lon_start,ineris_lon_end
    WRITE(*,*) 'NAEI:   ',naei_lat_start,naei_lat_end,naei_lon_start, &
         naei_lon_end
    WRITE(*,*) 'LANUV:  ',lanuv_lat_start,lanuv_lat_end,lanuv_lon_start, &
         lanuv_lon_end
    WRITE(*,*) 'TNO:  ',tno_lat_start,tno_lat_end,tno_lon_start, &
         tno_lon_end

    WRITE(*,*)
    WRITE(*,*) 'retro_interval  = ',retro_interval  
    WRITE(*,*) 'ineris_interval = ',ineris_interval
    WRITE(*,*) 'naei_interval   = ',naei_interval
    WRITE(*,*) 'lanuv_interval  = ',lanuv_interval
    WRITE(*,*) 'retro_tno       = ',tno_interval


    ! Regridding from emission inventories to WRF grid
    WRITE(*,*) ' - Start regridding of emissions to WRF grid'
    DO j = 1, jx
       WRITE(*,*) 'j = ',j
       DO i = 1, ix

          wrf_i = FLOAT(i)+0.5
          wrf_j = FLOAT(j)+0.5

          lat = wrf_lats(i,j)
          lon = wrf_lons(i,j)
          
          ! Calculate lat and lon of WRF grid box
          ! (it will differ from wrf_lats and wrf_lon when read_wrfinput=.true.)
          CALL MAPCF(wrf_i,wrf_j,lat_mid,lon_mid)

          ! Divide into smaller gridboxes based on emission inventory resolution
          ! First find out how many small grid boxes we need
          IF(lanuv(d) .AND. (lat .GT. lanuv_lat_start) .AND. (lat .LT. &
               lanuv_lat_end) .AND. (lon .GT. lanuv_lon_start) .AND. &
               (lon .LT. lanuv_lon_end)) THEN
             nboxes = MAX(NINT(spacing/lanuv_interval),1) ! WRF-res(km)/interval(km)
          ELSE IF(naei(d) .AND. (lat .GT. naei_lat_start) .AND. (lat .LT. &
               naei_lat_end) .AND. (lon .GT. naei_lon_start) .AND. &
               (lon .LT. naei_lon_end)) THEN
             nboxes = MAX(NINT(spacing/naei_interval),1) ! WRF-res(km)/interval(km)
          ELSE IF(ineris(d) .AND. (lat .GT. ineris_lat_start) .AND. (lat .LT. &
               ineris_lat_end) .AND. (lon .GT. ineris_lon_start) .AND. &
               (lon .LT. ineris_lon_end)) THEN
             nboxes = MAX(NINT(spacing/ineris_interval),1)! WRF-res(km)/interval(km)
          ELSE IF(retro(d)) THEN
             nboxes = MAX(NINT(spacing/retro_interval),1) ! WRF-res(km)/interval(km)
          ELSE IF(retrotno(d) .AND. (lat .GT. tno_lat_start) .AND. (lat .LT. &
               tno_lat_end) .AND. (lon .GT. tno_lon_start) .AND. &
               (lon .LT. tno_lon_end)) THEN
             nboxes = MAX(NINT(spacing/tno_interval),1) ! WRF-res(km)/interval(km)
          ELSE
!             WRITE(*,FMT='(A,I6,A,I6)') 'WARNING: No emission inventory covers the chosen WRF grid point! i = ',i,', j = ',j
          END IF

          small_interval = 1./nboxes ! This is the real interval that is used
          start_value = (-1)*((nboxes-1)/2)*small_interval

          ! DEBUG:
          !WRITE(*,*) 'nboxes: ',nboxes,', small_interval: ',small_interval, ', start_value: ',start_value


          ! An example to describe the method
          ! *********************************
          !  -----
          ! | | | |
          ! |- - -
          ! | |x| |
          ! |- - -
          ! |o| | |
          !  -----
          !
          ! x = wrf_i = 1.5
          !
          ! spacing        = 5 km, and within ineris emission inventory interval
          ! nboxes         = 5 km / 2 km = 2.5 ~ 3
          ! small_interval = 1/3 km = 0.333
          ! start_value    = (-1)*((3-1)/2)*0.333 = -0.333
          !
          ! o = small_wrf_i = 1.5 - 0.333 + 0 = 1.167
          !

          mean_value(:,:,:) = 0.d0

          ! Then loop through all the small grid boxes
          DO j_s = 1, nboxes
             DO i_s = 1, nboxes

                value(:,:,:) = 0.d0

                small_wrf_i = wrf_i + start_value + small_interval*(i_s-1)
                small_wrf_j = wrf_j + start_value + small_interval*(j_s-1)

                ! Calculate new lat and lon, of small grid box
                CALL MAPCF(small_wrf_i,small_wrf_j,lat_temp,lon_temp)

                ! Use the difference from mid point to find the most
                ! correct lat/long
                lat_small = wrf_lats(i,j) - (lat_mid-lat_temp)
                lon_small = wrf_lons(i,j) - (lon_mid-lon_temp)

                ! Now find the emission inventory that is valid for that point:

                ! *** RETRO *** (most probable choice)
                IF(retro(d)) THEN

                   ! Calculate indices of RETRO inventory
                   lat_temp = lat_small + 90.25
                   lon_temp = lon_small + 0.25
                   IF(lon_temp .LT. 0.25) THEN
                      lon_temp = lon_temp + 360.
                   END IF
                   ! Multiply by 2 and round off to the nearest integer:
                   ir = NINT(lon_temp*2) ! Calculate longitude
                   jr = NINT(lat_temp*2) ! Calculate latitude

                   ! Check for array out of bounds
                   ir = MIN(ir,720)
                   ir = MAX(ir,1)
                   jr = MIN(jr,360)
                   jr = MAX(jr,1)

                   IF(ineris(d) .AND. retro_ineris(ir,jr)) THEN
                      ! Calculate indices of INERIS inventory
                      ii = NINT((lon*10.) - i_istart + 1)
                      ji = NINT((lat*10.) - i_jstart + 1)

                      ! Check boundaries
                      ii = MAX(ii,1)
                      ji = MAX(ji,1)
                      ii = MIN(ii,i_iemis)
                      ji = MIN(ji,i_jemis)

                      IF(lanuv(d) .AND. ineris_lanuv(ii,ji)) THEN

                         ! Find approximate LANUV indices to speed things up later:
                         IF(i_s .EQ. 1 .AND. j_s .EQ. 1) THEN
                            CALL find_lanuv_index(lon,lat,1,l_iemis_short,1,l_jemis,il_big,jl_big)
                         END IF

                         ! Calculate indices of LANUV inventory
                         CALL find_lanuv_index(lon_small,lat_small, &
                              MAX(il_big-NINT(spacing),1), &
                              MIN(il_big+NINT(spacing),l_iemis_short), &
                              MAX(jl_big-NINT(spacing),1), &
                              MIN(jl_big+NINT(spacing),l_jemis), &
                              il,jl)
                          write(6,*)'here1'
                         CALL get_lanuv(lon_small,lat_small,il,jl,value) ! Use LANUV

                      ELSE IF(naei(d) .AND. ineris_naei(ii,ji)) THEN

                         ! Find approximate NAEI indices to speed things up later:
                         IF(i_s .EQ. 1 .AND. j_s .EQ. 1) THEN
                            CALL find_naei_index(lon,lat,1,n_iemis,1,n_jemis,in_big,jn_big)
                         END IF

                         ! Calculate indices of NAEI inventory
                         CALL find_naei_index(lon_small,lat_small, &
                              MAX(in_big-NINT(spacing),1), &
                              MIN(in_big+NINT(spacing),n_iemis), &
                              MAX(jn_big-NINT(spacing),1), &
                              MIN(jn_big+NINT(spacing),n_jemis), &
                              in,jn)
                         write(6,*)'here2'
                         CALL get_naei(in,jn,value) ! Use NAEI

                      ELSE
                        write(6,*)'here3'
                         CALL get_ineris(ii,ji,value) ! Use INERIS
                      END IF

!                   ELSEIF(naei(d)) THEN 
!                   ! Added to allow use of NAEI with RETRO. Use NAEI preferably. THIS DOES NOT WORK PROPERLY YET. PROBLEMS WITH UK LANDMASK.
!                   ! If outside NAEI range use RETRO. TP 06/09/10
!                      IF(lon_small.GE.naei_lon_start .AND. lon_small.LE.naei_lon_end .AND. &
!                         lat_small.GE.naei_lat_start .AND. lat_small.LE.naei_lat_end) THEN !lat and lon within extent of NAEI grid
!
!                         ! Find approximate NAEI indices to speed things up later:
!                         IF(used_naei.EQV..false.) THEN
!                            CALL find_naei_index(lon,lat,1,n_iemis,1,n_jemis,in_big,jn_big)
!                         END IF
!                         used_naei=.true. ! NAEI has now been used on this loop through the small grid boxes
!
!                         ! Calculate indices of NAEI inventory
!                         CALL find_naei_index(lon_small,lat_small, &
!                              MAX(in_big-NINT(spacing),1), &
!                              MIN(in_big+NINT(spacing),n_iemis), &
!                              MAX(jn_big-NINT(spacing),1), &
!                              MIN(jn_big+NINT(spacing),n_jemis), &
!                              in,jn)
!                         CALL get_naei(in,jn,value) ! Use NAEI
!
!                      ELSE
!                         ! Do nothing and get value from RETRO instead
!                         CALL get_retro(ir,jr,value) ! Use RETRO
!                      ENDIF

                   ELSE
                     write(6,*)'here4'
                      CALL get_retro(ir,jr,value) ! Use RETRO
                   END IF

                ! *** INERIS ***
                ELSE IF(ineris(d)) THEN

                   ! Calculate indices of INERIS inventory
                   ii = NINT((lon*10.) - i_istart + 1)
                   ji = NINT((lat*10.) - i_jstart + 1)
                   
                   ! Check if WRF grid point is covered by INERIS emissions
                   IF((ii .GE. 1) .AND. (ii .LE. i_iemis) .AND. (ji .GE. 1) &
                        .AND. (ji .LE. i_jemis)) THEN

                      ! Check boundaries
                      ii = MAX(ii,1)
                      ji = MAX(ji,1)
                      ii = MIN(ii,i_iemis)
                      ji = MIN(ji,i_jemis)

                      IF(lanuv(d) .AND. ineris_lanuv(ii,ji)) THEN

                         ! Find approximate LANUV indices to speed things up later:
                         IF(i_s .EQ. 1 .AND. j_s .EQ. 1) THEN
                            CALL find_lanuv_index(lon,lat,1,l_iemis_short,1,l_jemis,il_big,jl_big)
                         END IF

                         ! Calculate indices of LANUV inventory
                         CALL find_lanuv_index(lon_small,lat_small, &
                              MAX(il_big-NINT(spacing),1), &
                              MIN(il_big+NINT(spacing),l_iemis_short), &
                              MAX(jl_big-NINT(spacing),1), &
                              MIN(jl_big+NINT(spacing),l_jemis), &
                              il,jl)
                          write(6,*)'here5'
                         CALL get_lanuv(lon_small,lat_small,il,jl,value) ! Use LANUV

                      ELSE IF(naei(d) .AND. ineris_naei(ii,ji)) THEN

                         ! Find approximate NAEI indices to speed things up later:
                         IF(i_s .EQ. 1 .AND. j_s .EQ. 1) THEN
                            CALL find_naei_index(lon,lat,1,n_iemis,1,n_jemis,in_big,jn_big)
                         END IF

                         ! Calculate indices of NAEI inventory
                         CALL find_naei_index(lon_small,lat_small, &
                              MAX(in_big-NINT(spacing),1), &
                              MIN(in_big+NINT(spacing),n_iemis), &
                              MAX(jn_big-NINT(spacing),1), &
                              MIN(jn_big+NINT(spacing),n_jemis), &
                              in,jn)
!                         write(6,*)'here6'
                         CALL get_naei(in,jn,value) ! Use NAEI

                      ELSE
                         CALL get_ineris(ii,ji,value) ! Use INERIS
                      END IF
                   ELSE
                      WRITE(*,*) 'WARNING: INERIS emission'  
!     &inventory does not cover this WRF grid point! Set retro = .true.  
!     &in namelist to fix this problem.'
                   END IF

                ! *** NAEI ***
                ELSE IF(naei(d)) THEN
                   ! Find approximate NAEI indices to speed things up later:
                   IF(i_s .EQ. 1 .AND. j_s .EQ. 1) THEN
                      CALL find_naei_index(lon,lat,1,n_iemis,1,n_jemis,in_big,jn_big)
                   END IF

                   ! Calculate indices of NAEI inventory
                   CALL find_naei_index(lon_small,lat_small, &
                        MAX(in_big-NINT(spacing),1), &
                        MIN(in_big+NINT(spacing),n_iemis), &
                        MAX(jn_big-NINT(spacing),1), &
                        MIN(jn_big+NINT(spacing),n_jemis), &
                        in,jn)
!                       write(6,*)'here7',in,jn
                   CALL get_naei(in,jn,value) ! Use NAEI
                   
!                ! *** TNO ***


                ELSE IF(retrotno(d)) THEN
                   ! Calculate indices of TNO inventory
                   it = NINT((lon*8.) - t_istart + 0.5)
                   jt = NINT((lat*16.) - t_jstart + 0.5)

                   ! Check if WRF grid point is covered by TNO emissions
                   IF((it .GE. 1) .AND. (it .LE. t_iemis) .AND. (jt .GE. 1) &
                        .AND. (jt .LE. t_jemis)) THEN

                      ! Check boundaries
                      it = MAX(it,1)
                      jt = MAX(jt,1)
                      it = MIN(it,t_iemis)
                      jt = MIN(jt,t_jemis)

                       CALL get_tno(it,jt,value) ! Use TNO
                   ELSE
                      WRITE(*,*) 'WARNING: TNO emission inventory does not cover this WRF grid point!'
                      WRITE(*,*) ' Set retro = .true. in namelist to fix this problem.'
                   END IF




!                ! *** LANUV ***
                ELSE IF(lanuv(d)) THEN
                   WRITE(*,*) 'WARNING: Option not supported! ...ineris'
! In order  &
! to use LANUV emission inventory you must set ineris = .true.!'
                ELSE
                   WRITE(*,*) 'WARNING: No emission inventory chosen'
! or  &
! none of the chosen emission inventories covers this WRF grid point!'
                END IF

                ! Adding value from small gridbox to the mean value of the large gridbox
                mean_value(:,:,:) = mean_value(:,:,:) + (value(:,:,:)/(nboxes*nboxes))

             END DO ! i_s
          END DO    ! j_s

          ! Then we add the right value to the large matrix:
          DO n = 1, nchem
             DO k = 1, kx
                emi5d(i,k,j,n,:) = mean_value(k,n,:)
!             if(n.eq.6)write(6,*)'nh3=',mean_value(k,n,:),k,n
             END DO
          END DO

       END DO ! i
    END DO    ! j

    WRITE(*,*) ' - Done regridding'

  END SUBROUTINE wrf_grid
