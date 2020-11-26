  SUBROUTINE init()
    USE emission
    IMPLICIT NONE

    ! Moved namelist specification to top in line with standard to avert compiler error - TP 25/08/10
    namelist/general/path,year,month,io_style_emissions,run_days,chem_scheme
    namelist/mapinfo/max_dom,id_dom,ix2,jx2,kx,DX2,XLATC2,XLONC2, &
         CLAT1,CLAT2,IPROJ,heights,read_geofile,read_wrfinput
    namelist/reading/lanuv,naei,ineris,retro,retrotno
    namelist/writing/wrfchemi,binary,ascii,plot
    namelist/scaling/scaling_hour,timezone,scaling_emep_hour, &
         daylight_savings,scaling_emep_daily,day_of_week_init, &
         scaling_season
    namelist/intervals/lanuv_interval,naei_interval,ineris_interval, &
         retro_interval,tno_interval
    namelist/misc/naei_pointsrc,london,london_pointsrc,ruhr,coarser, &
coarser2_ratio,coarser_city,coarser2_city_ratio,removing_area,city_region

    INTEGER :: i,k,k2,low_k, high_k
    REAL    :: overlap

    ! We initialize the variables:
    namelist = 'namelist'
    DXBIGDO = 25000.
    ! INEST1 = Nesting or no nesting flag  (no nesting => inest1 = 0)
    INEST1 = 0
    XNESSTR = 1.
    YNESSTR = 1.
    REKM = 6370.  ! REKM = Earth radius (km)

    ! Then we read variables from the namelist:
    OPEN(13, FILE=namelist, FORM='FORMATTED')
    READ(13, NML=general)
    READ(13, NML=mapinfo)
    READ(13, NML=reading)
    READ(13, NML=writing)
    READ(13, NML=scaling)
    READ(13, NML=intervals)
    READ(13, NML=misc)
    CLOSE(13)

    ! Print info from namelist:
    WRITE(*,*) 'INPUT OPTIONS:'
    WRITE(*,*) '**************'
    WRITE(*,*) 'Month:             ', month
    WRITE(*,*) 'Year:              ', year
    WRITE(*,*) 'Time zone:         ', timezone
    IF(chem_scheme.EQ.1) THEN
       WRITE(*,*) 'Using RADM2 chemical scheme'
    ELSEIF(chem_scheme.EQ.2) THEN
       WRITE(*,*) 'Using CBMZ chemical scheme. Note CMBZ currently only supported for NAEI.'
    ELSEIF(chem_scheme.EQ.3) THEN
       WRITE(*,*) 'Using CRIMECH chemical scheme'
    ELSEIF(chem_scheme.EQ.4) THEN
       WRITE(*,*) 'Using CRIMECH chemical scheme UK' 
    ELSEIF(chem_scheme.EQ.5) THEN
       WRITE(*,*) 'Using CBMZ chemical scheme UK' 
    ENDIF
    WRITE(*,*) 'Read LANUV data:                ', lanuv
    WRITE(*,*) 'Read NAEI data:                 ', naei
    WRITE(*,*) 'Read INERIS data:               ', ineris
    WRITE(*,*) 'Read RETRO data:                ', retro
    WRITE(*,*) 'Read TNO data:                  ', retrotno
    WRITE(*,*) 'Scale hour:                     ', scaling_hour
    WRITE(*,*) 'Scale season:                   ', scaling_season
    WRITE(*,*) 'Read NAEI point sources:        ', naei_pointsrc
    WRITE(*,*) 'Read London NAEI data:          ', london
    WRITE(*,*) 'Read London NAEI point sources: ', london_pointsrc
    WRITE(*,*) 'Read Ruhr LANUV data:           ', ruhr
    WRITE(*,*) 'Coarser background emissions:   ', coarser
    WRITE(*,*) 'Coarser city emissions:         ', coarser_city
    WRITE(*,*) 'Remove emissions in city:       ', removing_area
    WRITE(*,*) 'Number of domains: ', max_dom
    WRITE(*,*)
    DO i = 1, max_dom
       WRITE(*,*) 'DOMAIN ', id_dom(i)
       WRITE(*,*) 'x-dim: ',ix2(i)
       WRITE(*,*) 'y-dim: ',jx2(i)
       WRITE(*,*) 'z-dim: ',kx
       WRITE(*,*) 'Resolution  : ', DX2(i)
       WRITE(*,*) 'Centre point: ', XLATC2(i), XLONC2(i)
       WRITE(*,*)
    END DO

    ! Calculate vertical interpolation factors between WRF and NAEI
    !hfac_wrf_emi(k,k2) 10,4
    ! Skal gå fra emi-grid til WRF-grid
    ! Trenger heights og emi_heights

    ! Må legge inn en test som fikser hvis høyeste wrf-lag er for lavt

    WRITE(*,*) '- Calculating height approximations'
    
    DO k = 1, kx ! kx=5??
       
       DO k2 = 1, emi_k ! emi_k = 4

          low_k = MAX(heights(k),emi_heights(k2))
          high_k = MIN(heights(k+1),emi_heights(k2+1))
          overlap = high_k - low_k

          IF(overlap .GT. 0) THEN
             hfac_wrf_emi(k,k2) = overlap/(emi_heights(k2+1)-emi_heights(k2))
             WRITE(*,*) k,k2,overlap,hfac_wrf_emi(k,k2)
          END IF
       END DO

    END DO

!    WRITE(*,*) hfac_wrf_emi

  END SUBROUTINE init
