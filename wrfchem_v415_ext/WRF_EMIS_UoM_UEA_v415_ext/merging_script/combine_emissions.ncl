;
;  Script for merging TNO and NAEI emission files
;
;  This is based on the matlab script to do the same task
;    - treatment of the emission variables will be hardcoded, 
;      so that we can check for existence of NAEI emissions.
;
;   Rules for merging:
;     1) TNO emissions will only be taken where there are no
;         NAEI emissions (this will be based on a total summation of 
;         NAEI emissions, rather than done on a species by species basis).
;     2) TNO emissions are used for all OC_DOM, OC_TRA, BC_1, EC_1_25,
;               OC_25_10 and EC_25_10 emission variables
;     3) A multiplying factor of 1.6 is applied to OC_DOM, OC_TRA, and OC_25_10
;               emissions (to roughly account for the Carbon to "everything" mass ratio for OM)
;     4) OIN 2.5 PM species will be either the difference between
;               E_PM25 and the sum of BC_1 + EC_1_25 + OC_DOM + OC_TRA
;           or 10% of the E_PM25 mass (whichever is smaller)
;          (we'll make sure that all emissions are >0)
;     5) OIN 10 PM species will be the difference between
;               E_PM_10 and the sum of OC_25_10 + EC_25_10
;          (we'll make sure that all emissions are >0)
;   
;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

; define input and output files
;;;;;;;;;; line 33 ;;;;;;; change the calling script if this moves!!!!! ;;;;;;;;;;;
file_name = "wrfchemi_d01_2014-01-01_00:00:00"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
tno_dir  = "TNO_2016_files/"
naei_dir = "NAEI_2016_files/"
out_dir  = "combined_2016_files/"


; full list of variables (for calculating the total NAEI emissions)
var_full = (/"E_CO","E_NO","E_NO2","E_SO2","E_NH3",\
			"E_C2H6","E_C2H4","E_C5H8","E_TM123B","E_TM124B", \
			"E_TM135B","E_OETHTOL","E_METHTOL","E_PETHTOL","E_DIME35EB", \
			"E_HCHO","E_CH3CHO","E_C2H5CHO","E_KET","E_MEK", \
			"E_CH3OH","E_C2H5OH","E_C3H6","E_C2H2","E_BENZENE", \
			"E_NC4H10","E_TOLUENE","E_OXYL","E_C3H8","E_TBUT2ENE",\
			"E_CH3CO2H", \
			"E_PM25","E_OC_DOM","E_OC_TRA","E_BC_1","E_EC_1_25",\
			"E_PM_10","E_OC_25_10","E_EC_25_10","E_OIN_25","E_OIN_10"/)

; list of variables for which we filter TNO emissions to avoid clashs with
; NAEI emissions
var_filter = (/"E_CO","E_NO","E_NO2","E_SO2","E_NH3",\
			"E_C2H6","E_C2H4","E_C5H8","E_TM123B","E_TM124B", \
			"E_TM135B","E_OETHTOL","E_METHTOL","E_PETHTOL","E_DIME35EB", \
			"E_HCHO","E_CH3CHO","E_C2H5CHO","E_KET","E_MEK", \
			"E_CH3OH","E_C2H5OH","E_C3H6","E_C2H2","E_BENZENE", \
			"E_NC4H10","E_TOLUENE","E_OXYL","E_C3H8","E_TBUT2ENE",\
			"E_CH3CO2H", \
			"E_PM25", \
			"E_PM_10" /)

; OC->OM scaling factor
scale_oc = 1.6


;;;;;;;;;;;;;; operational section of script ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; open the files for processing
F_TNO  = addfile(tno_dir+file_name+".nc","r")
F_NAEI = addfile(naei_dir+file_name+".nc","r")
F_OUT  = addfile(out_dir+file_name+".nc","w")

;;; loop through variables to pull out for making full NAEI summation
;;;   (units don't matter here - it's a purely binary check)
var_tot = F_NAEI->$var_full(0)$
do ivar = 1,dimsizes(var_full)-1
	var_temp = F_NAEI->$var_full(ivar)$
	var_tot = var_tot + var_temp
	delete(var_temp)
end do

;; create data mask to apply to TNO data
naei_empty = var_tot * 0.0 ; initialise mask to correct dimensions
naei_empty = where(var_tot.gt.0.0, 0.0, 1.0) ; set mask to zero where there's NAEI data - and 1.0 where there isn't

;; loop through the variables to be combined (in a straightforward manner)
do ivar = 0,dimsizes(var_filter)-1
	; load data
	var_tno = F_TNO->$var_filter(ivar)$
	var_naei= F_NAEI->$var_filter(ivar)$
	; merge data files - applying the filter
	var_naei = var_naei + naei_empty*var_tno
	; save data
	F_OUT->$var_filter(ivar)$ = var_naei
	; delete variables
	delete(var_tno)
	delete(var_naei)
end do

;; process the aerosol variables (this is done explicitly, rather than with loops)
; load aerosol variables
tno_oc_dom = F_TNO->E_OC_DOM     ; these are the TNO data - as NAEI doesn't carry this info
tno_oc_tra = F_TNO->E_OC_TRA
tno_oc_25_10 = F_TNO->E_OC_25_10
tno_bc_1 = F_TNO->E_BC_1
tno_ec_1_25 = F_TNO->E_EC_1_25
tno_ec_25_10 = F_TNO->E_EC_25_10
comb_pm25 = F_OUT->E_PM25       ; these are the combined NAEI & TNO data
comb_pm10 = F_OUT->E_PM_10 	 	;   (so need to take these from the output file, not TNO or NAEI)

; apply the OC -> OM scaling factor
tno_oc_dom = tno_oc_dom * scale_oc
tno_oc_tra = tno_oc_tra * scale_oc
tno_oc_25_10 = tno_oc_25_10 * scale_oc

; determine OIN_25 & OIN_10
comb_oin25 = comb_pm25 - tno_bc_1 - tno_ec_1_25 - tno_oc_dom - tno_oc_tra
comb_oin10 = comb_pm10 - tno_oc_25_10 - tno_ec_25_10 

; use either the difference between PM2.5 emissions and emitted BC/OC mass,
;         or 10% of the PM2.5 emission mass (whichever is smaller)
comb_oin25 = where(comb_oin25.lt.(comb_pm25*0.1),comb_oin25,comb_pm25*0.1)

; make sure we've not got any negative emissions of OIN
comb_oin25 = where(comb_oin25.gt.0.0, comb_oin25, 0.0)
comb_oin10 = where(comb_oin10.gt.0.0, comb_oin10, 0.0)

; write the aerosol emission data into the combined files
F_OUT->E_OC_DOM    = tno_oc_dom
F_OUT->E_OC_TRA    = tno_oc_tra
F_OUT->E_OC_25_10  = tno_oc_25_10
F_OUT->E_BC_1      = tno_bc_1
F_OUT->E_EC_1_25   = tno_ec_1_25
F_OUT->E_EC_25_10  = tno_ec_25_10
F_OUT->E_PM25      = comb_pm25
F_OUT->E_PM_10     = comb_pm10
F_OUT->E_OIN_25    = comb_oin25
F_OUT->E_OIN_10    = comb_oin10

end


