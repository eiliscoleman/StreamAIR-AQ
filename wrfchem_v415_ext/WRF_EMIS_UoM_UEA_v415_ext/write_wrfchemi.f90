
! Write to netCDF files: wrfchemi_00z_d0X and wrfchemi_12z_d0X
  SUBROUTINE write_wrfchemi()
    USE netcdf
    USE emission
    IMPLICIT NONE

    INTEGER      :: i,k,j,n,s,start_hour,end_hour,hour,month_wrfinput,start_date,date, &
         comp,time_loop,ntime_loops,day_of_week
    INTEGER      :: nc2ename_radm(nradm), nc2ename_cbmz(ncbmz)
    INTEGER      :: nc2ename_crimech(ncrimech)
    INTEGER      :: nc2ename_cbmz_uk(ncbmz_uk)
    INTEGER      :: nc2ename_crimech_uk(ncrimech_uk)
    INTEGER      :: nc2ename_saprc99_uk(saprc99_uk)
    INTEGER,ALLOCATABLE :: nc2ename(:)
    INTEGER      :: days_in_month(12)
    REAL         :: old_value,new_value,sum_value
    REAL         :: emi3d_out(ix, jx, kx)
    CHARACTER(len=2) :: smonth, sdate
    CHARACTER(len=4) :: syear
    CHARACTER(len=12) :: e_varnames_radm(nradm)
    CHARACTER(len=12) :: e_varnames_cbmz(ncbmz)
    CHARACTER(len=12) :: e_varnames_cbmz_uk(ncbmz_uk)
    CHARACTER(len=12) :: e_varnames_crimech(ncrimech)
    CHARACTER(len=12) :: e_varnames_crimech_uk(ncrimech_uk)
    CHARACTER(len=12) :: e_varnames_saprc99_uk(saprc99_uk)
    CHARACTER(len=12), ALLOCATABLE :: e_varnames(:)
    CHARACTER(len=9) :: day_of_week_char(7)
    DATA days_in_month / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
    ! Change ename order to correct order in nc-file:
    DATA nc2ename_radm / 19, 1, 2, 11, 10, 7, 8, 9, 16, 12, 13, 14, &
         15, 18, 4, 3, 17, 5, 6, 30, 20, 21, 28, 29, 26, 27, 22, 23, 24, 25 /

    DATA nc2ename_cbmz / 19, 1, 2, 11, 10, 7, 8, 9, 16, 12, 13, 14, &
         15, 18, 4, 3, 17, 5, 6, 30, 20, 21, 28,29,26,27, &
         22,23,24,25,31,32,33,34,35,36,37 /

    DATA nc2ename_crimech /1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, &
         14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, &
         29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42,43,44 /

    DATA nc2ename_crimech_uk /1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, &
         14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, &
         29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,40,41 /

    DATA nc2ename_saprc99_uk /1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, &
         14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, &
         29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,40,41, 42, 43, 44, &
         45, 46, 47, 48, 49, 50/


    DATA nc2ename_cbmz_uk / 19, 1, 2, 11, 10, 7, 8, 9, 16, 12, 13, 14, &
         15, 18, 4, 3, 17, 5, 6,20,21,22,23,24,25,26,27,28,29,30 /


    DATA e_varnames_radm /    &
       'E_ISO ','E_SO2 ','E_NO  ','E_CO  ','E_ETH ','E_HC3 ','E_HC5 ','E_HC8 ', &
       'E_XYL ','E_OL2 ','E_OLT ','E_OLI ','E_TOL ','E_CSL ','E_HCHO','E_ALD ', &
       'E_KET ','E_ORA2','E_NH3 ','E_PM_10','E_PM25I','E_PM25J','E_ECI', &
       'E_ECJ','E_ORGI','E_ORGJ','E_SO4I','E_SO4J','E_NO3I','E_NO3J'/ 
    DATA e_varnames_cbmz /    &
       'E_ISO ','E_SO2 ','E_NO  ','E_CO  ','E_ETH ','E_HC3 ','E_HC5 ','E_HC8 ', &
       'E_XYL ','E_OL2 ','E_OLT ','E_OLI ','E_TOL ','E_CSL ','E_HCHO','E_ALD ', &
       'E_KET ','E_ORA2','E_NH3 ','E_PM_10','E_PM25I','E_PM25J','E_ECI', &
       'E_ECJ','E_ORGI','E_ORGJ','E_SO4I','E_SO4J','E_NO3I','E_NO3J','E_NO2', &
       'E_CH3OH','E_C2H5OH','E_SO4C','E_NO3C','E_ORGC','E_ECC'/
    DATA e_varnames_cbmz_uk /    &
       'E_ISO ','E_SO2 ','E_NO  ','E_CO  ','E_ETH ','E_HC3 ','E_HC5 ','E_HC8 ', &
       'E_XYL ','E_OL2 ','E_OLT ','E_OLI ','E_TOL ','E_CSL ','E_HCHO','E_ALD ', &
       'E_KET ','E_ORA2','E_NH3 ','E_NO2', &
       'E_CH3OH','E_C2H5OH','E_PM25    ', &
       'E_OC_DOM', 'E_OC_TRA','E_BC_1', 'E_EC_1_25', 'E_PM10', &
       'E_OC_25_10', 'E_EC_25_10'/
 

    DATA e_varnames_crimech /  &
       'E_CO','E_NO','E_NO2','E_SO2','E_NH3','E_C2H6','E_C2H4', &
       'E_C5H8','E_TM123B','E_TM124B','E_TM135B','E_OETHTOL','E_METHTOL', &
       'E_PETHTOL','E_DIME35EB','E_HCHO','E_CH3CHO','E_C2H5CHO','E_KET', &
       'E_MEK','E_CH3OH','E_C2H5OH','E_C3H6','E_C2H2','E_BENZENE', &
       'E_NC4H10','E_TOLUENE', &
       'E_OXYL','E_C3H8','E_TBUT2ENE','E_CH3CO2H','E_PM_10', &
        'E_PM25I','E_PM25J','E_ECI', &
       'E_ECJ','E_ORGI','E_ORGJ','E_SO4I','E_SO4J','E_NO3I','E_NO3J', &
       'E_ORGC','E_ECC'/
    DATA e_varnames_crimech_uk /  &
       'E_CO','E_NO','E_NO2','E_SO2','E_NH3','E_C2H6','E_C2H4','E_C5H8', &
       'E_TM123B','E_TM124B','E_TM135B','E_OETHTOL','E_METHTOL', &
       'E_PETHTOL','E_DIME35EB','E_HCHO','E_CH3CHO','E_C2H5CHO','E_KET', &
       'E_MEK','E_CH3OH','E_C2H5OH','E_C3H6','E_C2H2','E_BENZENE', &
       'E_NC4H10','E_TOLUENE', &
       'E_OXYL','E_C3H8','E_TBUT2ENE','E_CH3CO2H','E_PM25','E_OC_DOM', &
       'E_OC_TRA','E_BC_1','E_EC_1_25','E_PM_10', &
       'E_OC_25_10','E_EC_25_10','E_OIN_25','E_OIN_10'/
       
    DATA e_varnames_saprc99_uk /    &
       'E_CO','E_NO      ', 'E_NO2     ', 'E_SO2     ', &
       'E_NH3     ', 'E_C2H6    ','E_ETHENE','E_OLE1','E_OLE2   ',  &  
       'E_C3H8','E_C2H2','E_ALK3','E_ALK4','E_ALK5','E_C3H6',  & 
       'E_ARO1','E_ARO2','E_HCHO','E_CCHO','E_RCHO', &
       'E_ACET','E_MEK','E_ISOPRENE','E_TERP','E_SESQ','E_PHEN', &
       'E_CRES','E_MEOH','E_GLY','E_MGLY','E_BACL','E_ISOPROD',  &
       'E_METHACRO','E_MVK','E_PROD2','E_CH4','E_BALD','E_HCOOH', &
       'E_CCO_OH','E_RCO_OH','E_PM25','E_OC_DOM', &
       'E_OC_TRA','E_BC_1','E_EC_1_25','E_PM_10', &
       'E_OC_25_10','E_EC_25_10','E_OIN_25','E_OIN_10'/


    DATA day_of_week_char /    &
         'Monday','Tuesday','Wednesday','Thursday','Friday', &
         'Saturday','Sunday' /

    ! NetCDF declarations for wrfinput and wrfchemi files
    CHARACTER(len=80) :: file_nc, times_var_value
    INTEGER           :: ncid, x_dimid, y_dimid, x_dim, y_dim
    INTEGER           :: rec_dimid, datestrlen_dimid, west_east_dimid, south_north_dimid, &
         emissions_zdimid, times_varid, e_varids(nchem)

    INTEGER           :: dimids(4), dimids_time(2)
    INTEGER           :: start(4), count(4), start_times(2), count_times(2)
    INTEGER, PARAMETER:: datestrlen_dimlen = 19, fieldtype = 104
    CHARACTER(len=*), PARAMETER :: rec_dim = "Time"
    CHARACTER(len=*), PARAMETER :: datestrlen_dim = "DateStrLen"
    CHARACTER(len=*), PARAMETER :: west_east_dim = "west_east"
    CHARACTER(len=*), PARAMETER :: south_north_dim = "south_north"
    CHARACTER(len=*), PARAMETER :: emissions_zdim = "emissions_zdim"
    CHARACTER(len=*), PARAMETER :: times_var = "Times"

    ! Global attributes names (gan)
    CHARACTER(len=*), PARAMETER :: gan_title = "TITLE"
    CHARACTER(len=*), PARAMETER :: gan_start_date = "START_DATE"
    CHARACTER(len=*), PARAMETER :: gan_we_grid = "WEST-EAST_GRID_DIMENSION"
    CHARACTER(len=*), PARAMETER :: gan_sn_grid = "SOUTH-NORTH_GRID_DIMENSION"
    CHARACTER(len=*), PARAMETER :: gan_bt_grid = "BOTTOM-TOP_GRID_DIMENSION"
    CHARACTER(len=*), PARAMETER :: gan_dx = "DX"
    CHARACTER(len=*), PARAMETER :: gan_dy = "DY"
    CHARACTER(len=*), PARAMETER :: gan_gridtype = "GRIDTYPE"
    CHARACTER(len=*), PARAMETER :: gan_diff_opt = "DIFF_OPT"
    CHARACTER(len=*), PARAMETER :: gan_km_opt = "KM_OPT"
    CHARACTER(len=*), PARAMETER :: gan_damp_opt = "DAMP_OPT"
    CHARACTER(len=*), PARAMETER :: gan_dampcoef = "DAMPCOEF"
    CHARACTER(len=*), PARAMETER :: gan_khdif = "KHDIF"
    CHARACTER(len=*), PARAMETER :: gan_kvdif = "KVDIF"
    CHARACTER(len=*), PARAMETER :: gan_mp_physics = "MP_PHYSICS"
    CHARACTER(len=*), PARAMETER :: gan_ra_lw_physics = "RA_LW_PHYSICS"
    CHARACTER(len=*), PARAMETER :: gan_ra_sw_physics = "RA_SW_PHYSICS"
    CHARACTER(len=*), PARAMETER :: gan_sf_sfclay_physics = "SF_SFCLAY_PHYSICS"
    CHARACTER(len=*), PARAMETER :: gan_sf_surface_physics = "SF_SURFACE_PHYSICS"
    CHARACTER(len=*), PARAMETER :: gan_bl_pbl_physics = "BL_PBL_PHYSICS"
    CHARACTER(len=*), PARAMETER :: gan_cu_physics = "CU_PHYSICS"
    CHARACTER(len=*), PARAMETER :: gan_surface_input = "SURFACE_INPUT_SOURCE"
    CHARACTER(len=*), PARAMETER :: gan_sst_update = "SST_UPDATE"
    CHARACTER(len=*), PARAMETER :: gan_grid_fdda = "GRID_FDDA"
    CHARACTER(len=*), PARAMETER :: gan_gfdda_interval_m = "GFDDA_INTERVAL_M"
    CHARACTER(len=*), PARAMETER :: gan_gfdda_end_h = "GFDDA_END_H"
    CHARACTER(len=*), PARAMETER :: gan_grid_sfdda = "GRID_SFDDA"
    CHARACTER(len=*), PARAMETER :: gan_sgfdda_interval_m = "SGFDDA_INTERVAL_M"
    CHARACTER(len=*), PARAMETER :: gan_sgfdda_end_h = "SGFDDA_END_H"
    CHARACTER(len=*), PARAMETER :: gan_we_start_unstag = "WEST-EAST_PATCH_START_UNSTAG"
    CHARACTER(len=*), PARAMETER :: gan_we_end_unstag = "WEST-EAST_PATCH_END_UNSTAG"
    CHARACTER(len=*), PARAMETER :: gan_we_start_stag = "WEST-EAST_PATCH_START_STAG"
    CHARACTER(len=*), PARAMETER :: gan_we_end_stag = "WEST-EAST_PATCH_END_STAG"
    CHARACTER(len=*), PARAMETER :: gan_sn_start_unstag = "SOUTH-NORTH_PATCH_START_UNSTAG"
    CHARACTER(len=*), PARAMETER :: gan_sn_end_unstag = "SOUTH-NORTH_PATCH_END_UNSTAG"
    CHARACTER(len=*), PARAMETER :: gan_sn_start_stag = "SOUTH-NORTH_PATCH_START_STAG"
    CHARACTER(len=*), PARAMETER :: gan_sn_end_stag = "SOUTH-NORTH_PATCH_END_STAG"
    CHARACTER(len=*), PARAMETER :: gan_bt_start_unstag = "BOTTOM-TOP_PATCH_START_UNSTAG"
    CHARACTER(len=*), PARAMETER :: gan_bt_end_unstag = "BOTTOM-TOP_PATCH_END_UNSTAG"
    CHARACTER(len=*), PARAMETER :: gan_bt_start_stag = "BOTTOM-TOP_PATCH_START_STAG"
    CHARACTER(len=*), PARAMETER :: gan_bt_end_stag = "BOTTOM-TOP_PATCH_END_STAG"
    CHARACTER(len=*), PARAMETER :: gan_grid_id = "GRID_ID"
    CHARACTER(len=*), PARAMETER :: gan_parent_id = "PARENT_ID"
    CHARACTER(len=*), PARAMETER :: gan_i_parent_start = "I_PARENT_START"
    CHARACTER(len=*), PARAMETER :: gan_j_parent_start = "J_PARENT_START"
    CHARACTER(len=*), PARAMETER :: gan_parent_grid_ratio = "PARENT_GRID_RATIO"
    CHARACTER(len=*), PARAMETER :: gan_dt = "DT"
    CHARACTER(len=*), PARAMETER :: gan_cen_lat = "CEN_LAT"
    CHARACTER(len=*), PARAMETER :: gan_cen_lon = "CEN_LON"
    CHARACTER(len=*), PARAMETER :: gan_truelat1 = "TRUELAT1"
    CHARACTER(len=*), PARAMETER :: gan_truelat2 = "TRUELAT2"
    CHARACTER(len=*), PARAMETER :: gan_moad_cen_lat = "MOAD_CEN_LAT"
    CHARACTER(len=*), PARAMETER :: gan_stand_lon = "STAND_LON"
!    CHARACTER(len=*), PARAMETER :: gan_pole_lat = "POLE_LAT"
!    CHARACTER(len=*), PARAMETER :: gan_pole_lon = "POLE_LON"
    CHARACTER(len=*), PARAMETER :: gan_gmt = "GMT"
    CHARACTER(len=*), PARAMETER :: gan_julyr = "JULYR"
    CHARACTER(len=*), PARAMETER :: gan_julday = "JULDAY"
    CHARACTER(len=*), PARAMETER :: gan_map_proj = "MAP_PROJ"
    CHARACTER(len=*), PARAMETER :: gan_mminlu = "MMINLU"
    CHARACTER(len=*), PARAMETER :: gan_num_land_cat = "NUM_LAND_CAT"
    CHARACTER(len=*), PARAMETER :: gan_iswater = "ISWATER"
    CHARACTER(len=*), PARAMETER :: gan_islake = "ISLAKE"
    CHARACTER(len=*), PARAMETER :: gan_isice = "ISICE"
    CHARACTER(len=*), PARAMETER :: gan_isurban = "ISURBAN"
    CHARACTER(len=*), PARAMETER :: gan_isoilwater = "ISOILWATER"

    ! Global attributes values (gav)
!    CHARACTER(len=*), PARAMETER :: gav_title = "OUTPUT FROM WRF-CHEM V3.2.1 EMISSIONS PREPROCESSOR"
    CHARACTER(len=80) :: gav_start_date, gav_gridtype, gav_mminlu
    INTEGER          :: gav_we_grid, gav_sn_grid, gav_bt_grid, gav_diff_opt
    INTEGER          ::  gav_km_opt, gav_damp_opt, gav_mp_physics, gav_ra_lw_physics
    INTEGER          ::  gav_ra_sw_physics, gav_sf_sfclay_physics, gav_sf_surface_physics
    INTEGER          ::  gav_bl_pbl_physics, gav_cu_physics, gav_surface_input
    INTEGER          ::  gav_sst_update, gav_grid_fdda, gav_gfdda_interval_m
    INTEGER          ::  gav_gfdda_end_h, gav_grid_sfdda, gav_sgfdda_interval_m
    INTEGER          ::  gav_sgfdda_end_h, gav_we_start_unstag, gav_we_end_unstag
    INTEGER          ::  gav_we_start_stag, gav_we_end_stag, gav_sn_start_unstag
    INTEGER          ::  gav_sn_end_unstag, gav_sn_start_stag, gav_sn_end_stag
    INTEGER          ::  gav_bt_start_unstag, gav_bt_end_unstag, gav_bt_start_stag
    INTEGER          ::  gav_bt_end_stag, gav_grid_id, gav_parent_id, gav_i_parent_start
    INTEGER          ::  gav_j_parent_start, gav_parent_grid_ratio, gav_julyr, gav_julday
    INTEGER          ::  gav_map_proj, gav_num_land_cat, gav_iswater
    INTEGER          ::  gav_islake, gav_isice, gav_isurban, gav_isoilwater
    REAL             :: gav_dx, gav_dy, gav_dampcoef, gav_khdif, gav_kvdif
    REAL             ::  gav_dt, gav_cen_lat, gav_cen_lon, gav_truelat1, gav_truelat2
    REAL             ::  gav_moad_cen_lat, gav_stand_lon, gav_gmt !, gav_pole_lat, gav_pole_lon

  !What chemical scheme are we using? Set appropriate parameters for array size and output species. - TP 01/12/10
  IF(chem_scheme.EQ.1) THEN
     ALLOCATE(nc2ename(nradm))
     nc2ename=nc2ename_radm
     ALLOCATE(e_varnames(nradm))
     e_varnames=e_varnames_radm
  ELSEIF(chem_scheme.EQ.2) THEN
     ALLOCATE(nc2ename(ncbmz))
     nc2ename=nc2ename_cbmz
     ALLOCATE(e_varnames(ncbmz))
     e_varnames=e_varnames_cbmz
  ELSEIF(chem_scheme.EQ.3) THEN
     ALLOCATE(nc2ename(ncrimech))
     nc2ename=nc2ename_crimech
     ALLOCATE(e_varnames(ncrimech))
     e_varnames=e_varnames_crimech
  ELSEIF(chem_scheme.EQ.4) THEN
     ALLOCATE(nc2ename(ncrimech_uk))
     nc2ename=nc2ename_crimech_uk
     ALLOCATE(e_varnames(ncrimech_uk))
     e_varnames=e_varnames_crimech_uk
  ELSEIF(chem_scheme.EQ.5) THEN
     ALLOCATE(nc2ename(ncbmz_uk))
     nc2ename=nc2ename_cbmz_uk
     ALLOCATE(e_varnames(ncbmz_uk))
     e_varnames=e_varnames_cbmz_uk
  ELSEIF(chem_scheme.EQ.6) THEN
     ALLOCATE(nc2ename(saprc99_uk))
     nc2ename=nc2ename_saprc99_uk
     ALLOCATE(e_varnames(saprc99_uk))
     e_varnames=e_varnames_saprc99_uk
  ELSE
     WRITE(*,*)'Supported chemical scheme not allocated'
     STOP
  ENDIF

    ! Open wrfinput to read dimensions and global attributes
    file_nc = 'wrfinput_d'//id_nr(d)
    WRITE(*,*) ' - Reading wrfinput_d'//id_nr(d)//' - checking dimensions'
    call check( nf90_open(file_nc, nf90_nowrite, ncid) )

    ! Get the dimids of the west_east and south_north coordinate variables.
    call check( nf90_inq_dimid(ncid, 'west_east', x_dimid) )
    call check( nf90_inq_dimid(ncid, 'south_north', y_dimid) )

    ! Read the dimension data.
    call check( nf90_Inquire_Dimension(ncid, x_dimid, len = x_dim) )
    call check( nf90_Inquire_Dimension(ncid, y_dimid, len = y_dim) )

    ! Check that dimensions are the same in wrfinput_d0X
    WRITE(*,*) 'WARNING: kx not being checked! Verify that kemit in namelist.input equals kx in namelist. kx = ',kx
    IF(ix .NE. x_dim) THEN
       WRITE(*,*) 'ERROR: x dimensions do not agree!'
       WRITE(*,*) 'From namelist: ix        = ',ix
       WRITE(*,*) 'From wrfinput: west_east = ',x_dim
       STOP
    END IF
    IF(jx .NE. y_dim) THEN
       WRITE(*,*) 'ERROR: y dimensions do not agree!'
       WRITE(*,*) 'From namelist: jx          = ',jx
       WRITE(*,*) 'From wrfinput: south_north = ',y_dim
       STOP
    END IF

    ! Read global attribute values from wrfinput file
    call check( nf90_get_att(ncid, nf90_global, gan_start_date, gav_start_date) )
    WRITE(*,*) 'Start time in wrfinput file: ', gav_start_date
    call check( nf90_get_att(ncid, nf90_global, gan_we_grid, gav_we_grid) )
    call check( nf90_get_att(ncid, nf90_global, gan_sn_grid, gav_sn_grid) )
    call check( nf90_get_att(ncid, nf90_global, gan_bt_grid, gav_bt_grid) )
    call check( nf90_get_att(ncid, nf90_global, gan_dx, gav_dx) )
    call check( nf90_get_att(ncid, nf90_global, gan_dy, gav_dy) )
    call check( nf90_get_att(ncid, nf90_global, gan_gridtype, gav_gridtype) )
    call check( nf90_get_att(ncid, nf90_global, gan_diff_opt, gav_diff_opt) )
    call check( nf90_get_att(ncid, nf90_global, gan_km_opt, gav_km_opt) )
    call check( nf90_get_att(ncid, nf90_global, gan_damp_opt, gav_damp_opt) )
    call check( nf90_get_att(ncid, nf90_global, gan_dampcoef, gav_dampcoef) )
    call check( nf90_get_att(ncid, nf90_global, gan_khdif, gav_khdif) )
    call check( nf90_get_att(ncid, nf90_global, gan_kvdif, gav_kvdif) )
    call check( nf90_get_att(ncid, nf90_global, gan_mp_physics, gav_mp_physics) )
    call check( nf90_get_att(ncid, nf90_global, gan_ra_lw_physics, gav_ra_lw_physics) )
    call check( nf90_get_att(ncid, nf90_global, gan_ra_sw_physics, gav_ra_sw_physics) )
    call check( nf90_get_att(ncid, nf90_global, gan_sf_sfclay_physics, gav_sf_sfclay_physics) )
    call check( nf90_get_att(ncid, nf90_global, gan_sf_surface_physics, gav_sf_surface_physics) )
    call check( nf90_get_att(ncid, nf90_global, gan_bl_pbl_physics, gav_bl_pbl_physics) )
    call check( nf90_get_att(ncid, nf90_global, gan_cu_physics, gav_cu_physics) )
    call check( nf90_get_att(ncid, nf90_global, gan_surface_input, gav_surface_input) )
    call check( nf90_get_att(ncid, nf90_global, gan_sst_update, gav_sst_update) )
    call check( nf90_get_att(ncid, nf90_global, gan_grid_fdda, gav_grid_fdda) )
    call check( nf90_get_att(ncid, nf90_global, gan_gfdda_interval_m, gav_gfdda_interval_m) )
    call check( nf90_get_att(ncid, nf90_global, gan_gfdda_end_h, gav_gfdda_end_h) )
    call check( nf90_get_att(ncid, nf90_global, gan_grid_sfdda, gav_grid_sfdda) )
    call check( nf90_get_att(ncid, nf90_global, gan_sgfdda_interval_m, gav_sgfdda_interval_m) )
    call check( nf90_get_att(ncid, nf90_global, gan_sgfdda_end_h, gav_sgfdda_end_h) )
    call check( nf90_get_att(ncid, nf90_global, gan_we_start_unstag, gav_we_start_unstag) )
    call check( nf90_get_att(ncid, nf90_global, gan_we_end_unstag, gav_we_end_unstag) )
    call check( nf90_get_att(ncid, nf90_global, gan_we_start_stag, gav_we_start_stag) )
    call check( nf90_get_att(ncid, nf90_global, gan_we_end_stag, gav_we_end_stag) )
    call check( nf90_get_att(ncid, nf90_global, gan_sn_start_unstag, gav_sn_start_unstag) )
    call check( nf90_get_att(ncid, nf90_global, gan_sn_end_unstag, gav_sn_end_unstag) )
    call check( nf90_get_att(ncid, nf90_global, gan_sn_start_stag, gav_sn_start_stag) )
    call check( nf90_get_att(ncid, nf90_global, gan_sn_end_stag, gav_sn_end_stag) )
    call check( nf90_get_att(ncid, nf90_global, gan_bt_start_unstag, gav_bt_start_unstag) )
    call check( nf90_get_att(ncid, nf90_global, gan_bt_end_unstag, gav_bt_end_unstag) )
    call check( nf90_get_att(ncid, nf90_global, gan_bt_start_stag, gav_bt_start_stag) )
    call check( nf90_get_att(ncid, nf90_global, gan_bt_end_stag, gav_bt_end_stag) )
    call check( nf90_get_att(ncid, nf90_global, gan_grid_id, gav_grid_id) )
    call check( nf90_get_att(ncid, nf90_global, gan_parent_id, gav_parent_id) )
    call check( nf90_get_att(ncid, nf90_global, gan_i_parent_start, gav_i_parent_start) )
    call check( nf90_get_att(ncid, nf90_global, gan_j_parent_start, gav_j_parent_start) )
    call check( nf90_get_att(ncid, nf90_global, gan_parent_grid_ratio, gav_parent_grid_ratio) )
    call check( nf90_get_att(ncid, nf90_global, gan_dt, gav_dt) )
    call check( nf90_get_att(ncid, nf90_global, gan_cen_lat, gav_cen_lat) )
    call check( nf90_get_att(ncid, nf90_global, gan_cen_lon, gav_cen_lon) )
    call check( nf90_get_att(ncid, nf90_global, gan_truelat1, gav_truelat1) )
    call check( nf90_get_att(ncid, nf90_global, gan_truelat2, gav_truelat2) )
    call check( nf90_get_att(ncid, nf90_global, gan_moad_cen_lat, gav_moad_cen_lat) )
    call check( nf90_get_att(ncid, nf90_global, gan_stand_lon, gav_stand_lon) )
!    call check( nf90_get_att(ncid, nf90_global, gan_pole_lat, gav_pole_lat) )
!    call check( nf90_get_att(ncid, nf90_global, gan_pole_lon, gav_pole_lon) )
    call check( nf90_get_att(ncid, nf90_global, gan_gmt, gav_gmt) )
    call check( nf90_get_att(ncid, nf90_global, gan_julyr, gav_julyr) )
    call check( nf90_get_att(ncid, nf90_global, gan_julday, gav_julday) )
    call check( nf90_get_att(ncid, nf90_global, gan_map_proj, gav_map_proj) )
    call check( nf90_get_att(ncid, nf90_global, gan_mminlu, gav_mminlu) )
    call check( nf90_get_att(ncid, nf90_global, gan_num_land_cat, gav_num_land_cat) )
    call check( nf90_get_att(ncid, nf90_global, gan_iswater, gav_iswater) )
    call check( nf90_get_att(ncid, nf90_global, gan_islake, gav_islake) )
    call check( nf90_get_att(ncid, nf90_global, gan_isice, gav_isice) )
    call check( nf90_get_att(ncid, nf90_global, gan_isurban, gav_isurban) )
    call check( nf90_get_att(ncid, nf90_global, gan_isoilwater, gav_isoilwater) )

    ! Close the wrfinput file.
    call check( nf90_close(ncid) )

    
    ! Use two 12-h emission data files
    IF(io_style_emissions .EQ. 1) THEN
       ntime_loops = 2
    ELSEIF(io_style_emissions .EQ. 2) THEN
       READ(gav_start_date(6:7),'(I2)') month_wrfinput
       READ(gav_start_date(9:10),'(I2)') start_date
       IF(month_wrfinput .NE. month) THEN
          start_date = 1
       END IF
       ntime_loops = MIN0(run_days,days_in_month(month)-start_date+1)
       WRITE(syear,'(I4)') year
       smonth = '00'
       IF(month .LT. 10) THEN
          WRITE(smonth(2:2),'(I1)') month
       ELSE
          WRITE(smonth(1:2),'(I2)') month
       END IF
       day_of_week = day_of_week_init - 1
    ELSE
       WRITE(*,*) 'ERROR: Wrong value for io_style_emissions!'
       STOP
    END IF

    DO time_loop = 1,ntime_loops

       ! Open wrfchemi to write dimensions, variables and attributes
       IF(io_style_emissions .EQ. 1) THEN
          times_var_value = gav_start_date
          IF(time_loop .EQ. 1) THEN
             start_hour = 1
             end_hour = 12
             file_nc = 'wrfchemi_00z_d'//id_nr(d)
             WRITE(*,*) ' - Writing to netcdf file wrfchemi_00z_d'//id_nr(d)
          ELSEIF(time_loop .EQ. 2) THEN
             start_hour = 13
             end_hour = 24
             file_nc = 'wrfchemi_12z_d'//id_nr(d)
             WRITE(*,*) ' - Writing to netcdf file wrfchemi_12z_d'//id_nr(d)
             WRITE(gav_start_date(12:13),'(A2)') '12'
          END IF
       ELSEIF(io_style_emissions .EQ. 2) THEN
          start_hour = 1
          end_hour = 24
          date = start_date + time_loop - 1
          sdate = '00'
          IF(date .LT. 10) THEN
             WRITE(sdate(2:2),'(I1)') date
          ELSE
             WRITE(sdate(1:2),'(I2)') date
          END IF
          day_of_week = day_of_week + 1
          IF(day_of_week .GT. 7) THEN
             day_of_week = day_of_week - 7
          END IF
          times_var_value = syear//'-'//smonth//'-'//sdate//'_00:00:00'
          file_nc = 'wrfchemi_d'//id_nr(d)//'_'//syear//'-'//smonth//'-'//sdate//'_00:00:00'
          IF((naei(d).OR.retrotno(d)) .AND. scaling_emep_daily) THEN
             WRITE(*,*) ' - Writing to netcdf file ',TRIM(file_nc), &
                  ' (',TRIM(day_of_week_char(day_of_week)),')'
          ELSE
             WRITE(*,*) ' - Writing to netcdf file ',TRIM(file_nc)
          END IF
       ELSE
          WRITE(*,*) 'ERROR: Wrong value for io_style_emissions!'
          STOP
       END IF

       call check( nf90_create(file_nc, nf90_clobber, ncid) )

       ! Define the dimensions.
       call check( nf90_def_dim(ncid, rec_dim, nf90_unlimited, rec_dimid) )
       call check( nf90_def_dim(ncid, datestrlen_dim, datestrlen_dimlen, datestrlen_dimid) )
       call check( nf90_def_dim(ncid, west_east_dim, ix, west_east_dimid) )
       call check( nf90_def_dim(ncid, south_north_dim, jx, south_north_dimid) )
       call check( nf90_def_dim(ncid, emissions_zdim, kx, emissions_zdimid) )

       dimids_time = (/ datestrlen_dimid, rec_dimid /)
       dimids = (/ west_east_dimid, south_north_dimid, emissions_zdimid, rec_dimid /)

       ! Define the netCDF variables
       call check( nf90_def_var(ncid, times_var, NF90_CHAR, dimids_time, times_varid) )
       DO comp=1,nchem
          call check( nf90_def_var(ncid, TRIM(e_varnames(comp)), NF90_REAL, dimids, e_varids(comp)) )
       END DO

       ! Assign units attributes to the netCDF variables.
       DO comp=1,nchem
          call check( nf90_put_att(ncid, e_varids(comp), 'FieldType', fieldtype) )
          call check( nf90_put_att(ncid, e_varids(comp), 'MemoryOrder', 'XYZ') )
          IF(comp .EQ. 1) THEN
             call check( nf90_put_att(ncid, e_varids(comp), 'description', &
                  'Isoprene EMISSIONS (Anth. for RADM/RACM, Anth+Bio for CBMZ)') )
          ELSE
             call check( nf90_put_att(ncid, e_varids(comp), 'description', 'EMISSIONS') )
          END IF
          call check( nf90_put_att(ncid, e_varids(comp), 'units', 'mol km^-2 hr^-1') )
          call check( nf90_put_att(ncid, e_varids(comp), 'stagger', '') )
          call check( nf90_put_att(ncid, e_varids(comp), 'coordinates', 'XLONG XLAT') )
       END DO

       ! Write global attributes (we probably don't need all the attributes from wrfinput)
       call check( nf90_put_att(ncid, nf90_global, 'TITLE', 'Created by Steven Utembe') )
       call check( nf90_put_att(ncid, nf90_global, gan_start_date, gav_start_date) )
       call check( nf90_put_att(ncid, nf90_global, gan_we_grid, gav_we_grid) )
       call check( nf90_put_att(ncid, nf90_global, gan_sn_grid, gav_sn_grid) )
       call check( nf90_put_att(ncid, nf90_global, gan_bt_grid, gav_bt_grid) )
       call check( nf90_put_att(ncid, nf90_global, gan_dx, gav_dx) )
       call check( nf90_put_att(ncid, nf90_global, gan_dy, gav_dy) )
       call check( nf90_put_att(ncid, nf90_global, gan_cen_lat, gav_cen_lat) )
       call check( nf90_put_att(ncid, nf90_global, gan_cen_lon, gav_cen_lon) )
       call check( nf90_put_att(ncid, nf90_global, gan_truelat1, gav_truelat1) )
       call check( nf90_put_att(ncid, nf90_global, gan_truelat2, gav_truelat2) )
       call check( nf90_put_att(ncid, nf90_global, gan_moad_cen_lat, gav_moad_cen_lat) )
       call check( nf90_put_att(ncid, nf90_global, gan_stand_lon, gav_stand_lon) )
       call check( nf90_put_att(ncid, nf90_global, gan_map_proj, gav_map_proj) )
       call check( nf90_put_att(ncid, nf90_global, gan_mminlu, gav_mminlu) )

       ! End define mode.
       call check( nf90_enddef(ncid) )

       ! These settings tell netcdf to write one timestep of data. (The
       ! setting of start(4) inside the loop below tells netCDF which
       ! timestep to write.)
       count_times = (/ datestrlen_dimlen, 1 /)
       start_times = (/ 1, 1 /)

       count = (/ ix, jx, kx, 1 /)
       start = (/ 1, 1, 1, 1 /)

       ! Write data for each time step
       do hour = start_hour, end_hour
       
          start_times(2) = hour - start_hour + 1
          start(4) = hour - start_hour + 1

          ! First write to Times variable
          IF(hour .LE. 10) THEN
             WRITE(times_var_value(13:13),'(I1.1)') hour-1
          ELSE
             WRITE(times_var_value(12:13),'(I2.2)') hour-1
          END IF
          call check( nf90_put_var(ncid, times_varid, times_var_value, &
               start = start_times, count = count_times) )
!          write(6,*)'whats nchem=',nchem
          ! Then write the emission data
!       DO n=1,nradm_all
          DO comp=1,nchem
!       DO comp=1,1
             n = nc2ename(comp)
!             write(6,*)'n=',n
             emi3d_out(:,:,:) = 0.d0
!             IF(n.GT.0) THEN !Added to account for if a species is not passed to a WRF input - TP 02/12/10
             ! Write out 3-D emission arrays to netCDF wrfchemi* file
             DO i=1,ix
                DO k=1,kx
                   DO j=1,jx
                      sum_value = 0.
                      DO s=1,nsources
                         old_value = emi5d(i,k,j,n,s)
!                         if(n.eq.6)write(6,*)'nh3=',emi5d(i,k,j,n,s)
                         IF(old_value .GT. 0.0) THEN

                            ! Scale daily
!                            IF((io_style_emissions .EQ. 2) .AND. ineris(d) &
!                                 .AND. scaling_ineris_daily) THEN
                            IF((io_style_emissions .EQ. 2) .AND. (naei(d).OR.retrotno(d)) &
                                 .AND. scaling_emep_daily) THEN
                               
                               ! DEBUG:
!                               IF(n .LT. 6 .AND. i .EQ. 31 .AND. j .EQ. 54 .AND. k .EQ. 1) THEN
!                                  WRITE(*,'(4I3,A,ES14.7)') day_of_week,hour,n,s,'. Before scale_daily: ',old_value
!                               END IF
!                               if(i.eq.1.and.k.eq.1.and.j.eq.1)write(6,*)'calling scale daily'
                               CALL scale_daily(old_value,day_of_week,s,n,i,j)

                               ! DEBUG:
!                               IF(n .LT. 6 .AND. i .EQ. 31 .AND. j .EQ. 54 .AND. k .EQ. 1) THEN
!                                  WRITE(*,'(4I3,A,ES14.7)') day_of_week,hour,n,s,'. After scale_daily:  ',old_value
!                               END IF
!                               write(6,*)'back from scale_daily',s
                            END IF
                         
                            ! Scale hourly
!                            IF((ineris(d) .AND. scaling_ineris_hour) .OR. &
!                                 (scaling_hour)) THEN
                            IF(((naei(d).OR.retrotno(d)) .AND. scaling_emep_hour) .OR. &
                                 (scaling_hour)) THEN

                               CALL scale_hour(old_value,new_value,hour,s,i,j)
                               sum_value = sum_value + new_value
!                               write(6,*)'back from scale_hour'
                            ELSE
                               sum_value = sum_value + old_value
                            END IF

                         END IF

                      END DO
                      emi3d_out(i,j,k) = sum_value
                   END DO
                END DO
             END DO
!             ENDIF
             ! Put variable in netCDF file
             call check( nf90_put_var(ncid, e_varids(comp), emi3d_out, &
                  start = start, count = count) )
          END DO
       end do

       ! Close the wrfchemi file.
       call check( nf90_close(ncid) )

    END DO ! time_loop


  END SUBROUTINE write_wrfchemi

  subroutine check(status)
    use netcdf
    implicit none
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  
