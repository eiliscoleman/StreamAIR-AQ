MODULE emission
  IMPLICIT NONE

  INTEGER, PARAMETER :: nradm = 30, ncbmz=37, ncrimech=44,ndom = 5, &
       max_k = 10, emi_k = 4, nsources = 15, ncountries = 63, &
       nmonths = 12, ndays = 7, ncrimech_uk=41,ncbmz_uk=30,saprc99_uk=50
  INTEGER            :: nchem
  REAL               :: height_fac(nsources,emi_k), spacing, &
       hfac_wrf_emi(max_k,emi_k),season_fac(12,ncountries,nsources,12), &
       daily_fac(12,ncountries,nsources,7)
  REAL, ALLOCATABLE  :: emi3d(:,:,:), emi5d(:,:,:,:,:), &
       emi5d_orig(:,:,:,:,:),wrf_lats(:,:),wrf_lons(:,:)
  INTEGER, ALLOCATABLE::wrf_countries(:,:), timediff_from_utc(:,:)
  CHARACTER (len=80) :: path,namelist,city_region
  INTEGER            :: d,month,year,timezone,day_of_week_init, &
       io_style_emissions,run_days,emi_heights(emi_k+1),max_dom,INEST1, IPROJ
  INTEGER            :: ix,jx,kx,ILX,JLX, coarser_ratio, &
       coarser_city_ratio, chem_scheme
  REAL               :: lanuv_interval,naei_interval,ineris_interval, &
       retro_interval,tno_interval
  LOGICAL            :: read_wrfinput,read_geofile, &
       wrfchemi,binary,ascii,plot,scaling_hour,scaling_season, &
       scaling_emep_hour,daylight_savings,scaling_emep_daily, &
       removing_area,naei_pointsrc,london,london_pointsrc,ruhr,coarser, &
       coarser_city
  LOGICAL :: lanuv(ndom),naei(ndom),ineris(ndom),retro(ndom),retrotno(ndom)
  LOGICAL, ALLOCATABLE :: city_regions(:,:)
  CHARACTER (len=2), DIMENSION(ndom)     :: id_dom,id_nr
  CHARACTER (len=12), DIMENSION(nradm):: ename_radm
  CHARACTER (len=12), DIMENSION(ncbmz):: ename_cbmz
  CHARACTER (len=12), DIMENSION(ncbmz_uk):: ename_cbmz_uk
  CHARACTER (len=12), DIMENSION(ncrimech):: ename_crimech
  CHARACTER (len=12), DIMENSION(ncrimech_uk):: ename_crimech_uk
  CHARACTER (len=12), DIMENSION(saprc99_uk):: ename_saprc99_uk
  CHARACTER (len=12), ALLOCATABLE, DIMENSION(:) :: ename

  DATA id_nr /'01','02','03','04','05'/
  DATA ename_radm /    &
       'e_so2 ','e_no  ','e_ald ','e_hcho','e_ora2','e_nh3 ','e_hc3 ','e_hc5 ', &
       'e_hc8 ','e_eth ','e_co  ','e_ol2 ','e_olt ','e_oli ','e_tol ','e_xyl ', &
       'e_ket ','e_csl ','e_iso ','e_pm25i','e_pm25j','e_so4i','e_so4j',        &
       'e_no3i','e_no3j','e_orgi','e_orgj','e_eci','e_ecj','e_pm10'/
  DATA ename_cbmz /    &
       'e_so2 ','e_no  ','e_ald ','e_hcho','e_ora2','e_nh3 ','e_hc3 ','e_hc5 ', &
       'e_hc8 ','e_eth ','e_co  ','e_ol2 ','e_olt ','e_oli ','e_tol ','e_xyl ', &
       'e_ket ','e_csl ','e_iso ','e_pm25i','e_pm25j','e_so4i','e_so4j',        &
       'e_no3i','e_no3j','e_orgi','e_orgj','e_eci','e_ecj','e_pm10','e_no2',    &
       'e_ch3oh','e_c2h5oh','e_so4c','e_no3c','e_orgc','e_ecc'/
  DATA ename_crimech /    &
       'e_co      ', 'e_no      ', 'e_no2     ', 'e_so2     ', &
       'e_nh3     ', 'e_c2h6    ', 'e_c2h4    ', &
       'e_c5h8    ', 'e_tm123b  ', 'e_tm124b  ', &
       'e_tm135b  ', 'e_oethtol ', 'e_methtol ', 'e_pethtol ', &
       'e_dime35eb', 'e_hcho    ', 'e_ch3cho  ', 'e_c2h5cho ', &
       'e_ket', 'e_mek     ', 'e_ch3oh   ', 'e_c2h5oh  ', &
       'e_c3h6    ', 'e_c2h2    ', 'e_benzene ', 'e_nc4h10  ', &
       'e_toluene ', 'e_oxyl    ', 'e_c3h8    ', 'e_tbut2ene', &
       'e_ch3co2h ', 'e_pm25i   ', 'e_pm25j   ', 'e_so4i    ', &
       'e_so4j    ', 'e_no3i    ', 'e_no3j    ', 'e_orgi    ', &
       'e_orgj    ', 'e_eci     ', 'e_ecj     ', 'e_pm10    ', &
       'e_orgc    ', 'e_ecc     '/     
  DATA ename_cbmz_uk /    &   
       'e_so2 ','e_no  ','e_ald ','e_hcho','e_ora2','e_nh3 ','e_hc3 ','e_hc5 ', &
       'e_hc8 ','e_eth ','e_co  ','e_ol2 ','e_olt ','e_oli ','e_tol ','e_xyl ', &
       'e_ket ','e_csl ','e_iso ','e_no2','e_ch3oh','e_c2h5oh','e_pm25    ', &
       'e_oc_dom  ', 'e_oc_tra  ','e_bc_1    ', 'e_ec_1_25 ', 'e_pm10    ', &
       'e_oc_25_10', 'e_ec_25_10'/


  DATA ename_crimech_uk /    &
       'e_co      ', 'e_no      ', 'e_no2     ', 'e_so2     ', &
       'e_nh3     ', 'e_c2h6    ', 'e_c2h4    ', &
       'e_c5h8    ', 'e_tm123b  ', 'e_tm124b  ', &
       'e_tm135b  ', 'e_oethtol ', 'e_methtol ', 'e_pethtol ', &
       'e_dime35eb', 'e_hcho    ', 'e_ch3cho  ', 'e_c2h5cho ', &
       'e_ket', 'e_mek     ', 'e_ch3oh   ', 'e_c2h5oh  ', &
       'e_c3h6    ', 'e_c2h2    ', 'e_benzene ', 'e_nc4h10  ', &
       'e_toluene ', 'e_oxyl    ', 'e_c3h8    ', 'e_tbut2ene', &
       'e_ch3co2h ', 'e_pm25    ', 'e_oc_dom  ', 'e_oc_tra  ', &
       'e_bc_1    ', 'e_ec_1_25 ', 'e_pm10    ', &
       'e_oc_25_10', 'e_ec_25_10', 'e_oin_25  ', 'e_oin_10  '/  


  DATA ename_saprc99_uk /    &
       'e_co','e_no      ', 'e_no2     ', 'e_so2     ', &
       'e_nh3     ', 'e_c2h6    ','e_ethene','e_ole1','e_ole2   ',  &  
       'e_c3h8','e_c2h2','e_alk3','e_alk4','e_alk5','e_c3h6',  & 
       'e_aro1','e_aro2','e_hcho','e_ccho','e_rcho', &
       'e_acet','e_mek','e_isoprene','e_terp','e_sesq','e_phen', &
       'e_cres','e_meoh','e_gly','e_mgly','e_bacl','e_isoprod',  &
       'e_methacro','e_mvk','e_prod2','e_ch4','e_bald','e_hcooh', &
       'e_cco_oh','e_rco_oh','e_pm25    ', 'e_oc_dom  ', 'e_oc_tra  ', &
       'e_bc_1    ', 'e_ec_1_25 ', 'e_pm10    ', &
       'e_oc_25_10', 'e_ec_25_10', 'e_oin_25  ', 'e_oin_10  '/ 
       
  DATA emi_heights / 0, 50, 150, 250, 400 / ! Emission upper height limits
  !Array height_fac details the fraction of each NAEI source that should be attributed
  !to each NAEI emission level, as defined in emi_heights. - Deduced TP 27/09/10
!  DATA height_fac / &
!       0.1, 0.5, 0.5, 0.9, 0.9, 1.0, 1.0, 1.0, 0.8, 1.0, 1.0, &
!       0.2, 0.5, 0.5, 0.1, 0.1, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, &
!       0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
!       0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  &
!       /

  DATA height_fac / &
       0.1, 0.5, 0.5, 0.9, 0.9, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.8, 1.0, 1.0, &
       0.2, 0.5, 0.5, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, &
       0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
       0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  &
       /

  REAL       :: XNESSTR, YNESSTR, REKM, CLAT1, CLAT2, &
       DX, DXBIGDO, XLATC, XLONC 
!, wrf_i, wrf_j, smallest_diff_naei
  INTEGER, DIMENSION(ndom):: ix2,jx2,coarser2_ratio, &
       coarser2_city_ratio
  INTEGER, DIMENSION(max_k+1) :: heights
  REAL, DIMENSION(ndom)   :: DX2,XLATC2,XLONC2

END MODULE emission
