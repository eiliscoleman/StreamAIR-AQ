!**************************************
!* Author:  Øivind Hodnebrog          *
!* E-mail:  oivinho@geo.uio.no        *
!* Updated: November, 2010            *
!**************************************
! Updated T. Pugh, Dec 2010           *
!**************************************
!
! The function of this program is to read emissions data from
! one or several inventories, and produce emission input files for 
! use with the WRF-Chem model. See README for more info.
!
!*****************************************************************

PROGRAM main
  IMPLICIT NONE

  PRINT *, '************************************************'
  PRINT *, '*                                              *'
  PRINT *, '*   WRF/Chem emission preprocessor program     *'
  PRINT *, '*   --------------------------------------     *'
  PRINT *, '*                                              *'
  PRINT *, '* This program makes anthropogenic emissions   *'
  PRINT *, '* from several inventories available as input  *'
  PRINT *, '* to the WRF/Chem model.                       *'
  PRINT *, '*                                              *'
  PRINT *, '*        - See README for more info -          *'
  PRINT *, '*                                              *'
  PRINT *, '************************************************'
  PRINT *, ''

  PRINT *, ' - Reading namelist'
  PRINT *, ''
  CALL init()
  CALL main_prog()

END PROGRAM main

! This is the main program
SUBROUTINE main_prog()
  USE emission
  USE module_lanuv
  USE module_naei
  USE module_ineris
  USE module_emep
  USE module_retro
  USE module_tno
  IMPLICIT NONE

!Local variables
  INTEGER :: i
  LOGICAL :: data_init

  !What chemical scheme are we using? Set appropriate parameters for array size and output species. - TP 23/09/10
  IF(chem_scheme.EQ.1) THEN
     nchem=nradm
     ALLOCATE(ename(nradm))
     ename=ename_radm
  ELSEIF(chem_scheme.EQ.2) THEN
     nchem=ncbmz
     ALLOCATE(ename(ncbmz))
     ename=ename_cbmz
  ELSEIF(chem_scheme.EQ.3) THEN
     nchem=ncrimech
     ALLOCATE(ename(ncrimech))
     ename=ename_crimech
  ELSEIF(chem_scheme.EQ.4) THEN
     nchem=ncrimech_uk
     ALLOCATE(ename(ncrimech_uk))
     ename=ename_crimech_uk
  ELSEIF(chem_scheme.EQ.5) THEN
     nchem=ncbmz_uk
     ALLOCATE(ename(ncbmz_uk))
     ename=ename_cbmz_uk
  ELSEIF(chem_scheme.EQ.6) THEN
     nchem=saprc99_uk
     ALLOCATE(ename(saprc99_uk))
     ename=ename_saprc99_uk
  ELSE
     WRITE(*,*)'Supported chemical scheme not allocated'
     STOP
  ENDIF

  IF(scaling_season) THEN
     CALL scale_season()
  END IF
  IF(scaling_emep_daily) THEN
     CALL scale_daily_init()
  END IF

  ! *** LANUV ***
  data_init=.false.
  DO i=1,ndom !If any of the nests required LANUV then initialise the dataset - TP 29/11/10
  IF((lanuv(i).EQV..true.).AND.(data_init.EQV..false.)) THEN
     CALL allocate_lanuv() ! Allocate array size for naei variables - TP 18/09/10
     CALL read_lanuv()    ! Read variables from the LANUV inventory
     CALL latlon_lanuv()  ! Calculate lats and longs for the LANUV grid
     CALL plot_lanuv()    ! Write to plot files (plot with e.g. MATLAB)
     CALL convert_lanuv() ! Convert units to mole/(km^2*h)
     CALL vocs_lanuv()    ! Split VOCs based on sources and convert units
     CALL deallocate_lanuv()
     data_init=.true.
  END IF
  ENDDO

  ! *** NAEI ***
  data_init=.false.
  DO i=1,ndom
  IF((naei(i).EQV..true.).AND.(data_init.EQV..false.)) THEN
     CALL allocate_naei() ! Allocate array size for naei variables - TP 18/09/10
     CALL read_naei()     ! Read variables from the NAEI inventory
     CALL latlon_naei()   ! Read lats and longs for the NAEI grid
     CALL plot_naei()     ! Write to plot files (plot with e.g. MATLAB)
     CALL convert_naei()  ! Convert units to mole/(km^2*h)
     CALL vocs_naei()     ! Split VOCs based on sources and convert units
     CALL pm_naei()       ! Convert units for PM emissions to ug m-3 m s-1 - TP 10/09/10
     CALL deallocate_naei()
     data_init=.true.
  END IF
  ENDDO

  ! *** INERIS ***
  data_init=.false.
  DO i=1,ndom
  IF((ineris(i).EQV..true.).AND.(data_init.EQV..false.)) THEN
     CALL allocate_ineris() ! Allocate array size for naei variables - TP 18/09/10
     CALL read_ineris()    ! Read variables from the INERIS inventory
     CALL latlon_ineris()  ! Write lats and longs for the INERIS grid
     CALL plot_ineris()    ! Write to plot files (plot with e.g. MATLAB)
     CALL convert_ineris() ! Convert units to mole/(km^2*h)
     CALL vocs_ineris()    ! Split VOCs based on sources and convert units
     IF(lanuv(i)) THEN
!        CALL overlap_ineris_lanuv()  ! Find grid boxes that overlap
        CALL read_overlap_ineris_lanuv()  ! Find grid boxes that overlap
!!! REMOVE:
!        CALL overlap_lanuv_ineris()  ! Find grid boxes that overlap
!!!
        CALL merge_lanuv_ineris() ! Replace missing sectors with INERIS data
     END IF
     IF(naei(i)) THEN
!        CALL overlap_ineris_naei()   ! Find grid boxes that overlap
        CALL read_overlap_ineris_naei()  ! Find grid boxes that overlap
     END IF
     CALL deallocate_ineris()
     data_init=.true.
  END IF
  ENDDO

!!! NOT FINISHED YET:
  ! *** EMEP ***
!  data_init=.false.
!  DO i=1,ndom
!  IF(emep(i).EQV..true..AND.data_init.EQV..false.) THEN
!     CALL allocate_emep() ! Allocate array size for naei variables - TP 18/09/10
!     CALL read_emep()    ! Read variables from the INERIS inventory
!     CALL latlon_emep()  ! Write lats and longs for the INERIS grid
!     CALL plot_emep()    ! Write to plot files (plot with e.g. MATLAB)
!     IF(lanuv(i)) THEN
!        CALL overlap_emep_lanuv()  ! Find grid boxes that overlap
!        CALL read_overlap_emep_lanuv()  ! Find grid boxes that overlap
!!! REMOVE:
!        CALL overlap_lanuv_emep()  ! Find grid boxes that overlap
!!!
!        CALL merge_lanuv_ineris() ! Replace missing sectors with INERIS data
!     END IF
!     IF(naei(i)) THEN
!        CALL overlap_ineris_naei()   ! Find grid boxes that overlap
!        CALL read_overlap_ineris_naei()  ! Find grid boxes that overlap
!     END IF
!     CALL deallocate_emep()
!     data_init=.true.
!  END IF
!  ENDDO
!!!

  ! *** RETRO ***
  data_init=.false.
  DO i=1,ndom
  IF((retro(i).EQV..true.).AND.(data_init.EQV..false.)) THEN
     CALL allocate_retro() ! Allocate array size for naei variables - TP 18/09/10
     CALL read_retro()     ! Read variables from the RETRO inventory
     CALL latlon_retro()   ! Write lats and longs for the RETRO grid
     CALL plot_retro()     ! Write to plot files (plot with e.g. MATLAB)
     IF(ineris(i)) THEN
!        CALL overlap_retro_ineris()  ! Find grid boxes that overlap
        CALL read_overlap_retro_ineris()  ! Find grid boxes that overlap
     END IF
     CALL deallocate_retro()
     data_init=.true.
  END IF
  ENDDO

  ! *** TNO ***
  data_init=.false.
  DO i=1,ndom
  IF((retrotno(i).EQV..true.).AND.(data_init.EQV..false.)) THEN
  
     CALL allocate_tno() ! Allocate array size for tno variables - SU 16/02/11
     CALL read_tno()     ! Read variables from the tno inventory
     CALL latlon_tno()   ! Read lats and longs for the tno grid
     CALL plot_tno()     ! Write to plot files (plot with e.g. MATLAB)
     CALL convert_tno()  ! Convert units to mole/(km^2*h)
     CALL vocs_tno()     ! Split VOCs based on sources and convert units
     CALL pm_tno()       ! Convert units for PM emissions to ug m-3 m s-1 - TP 10/09/10
     CALL deallocate_tno()
     data_init=.true.
  END IF
  ENDDO


  ! This loop goes through the rest of the program once for each domain
  DO d = 1, max_dom
     ix = ix2(d)
     jx = jx2(d)
     DX = DX2(d)
     DXBIGDO = DX
     XLATC = XLATC2(d)
     XLONC = XLONC2(d)
     ILX = ix+1
     JLX = jx+1
     coarser_ratio = coarser2_ratio(d)
     coarser_city_ratio = coarser2_city_ratio(d)

     spacing = DX/1000.

     ! Have to allocate matrices now that we know the dimensions:
     ALLOCATE(emi5d(ix,kx,jx,nchem,nsources))
     emi5d = 0.D0
     ALLOCATE(emi3d(ix,kx,jx))
     emi3d = 0.D0

     ALLOCATE(city_regions(ix,jx))
     CALL define_regions()

     ALLOCATE(wrf_lats(ix,jx))
     wrf_lats = 0.D0
     ALLOCATE(wrf_lons(ix,jx))
     wrf_lons = 0.D0
     CALL get_wrflatlon()

     ! Find country number of each WRF grid box
     ALLOCATE(wrf_countries(ix,jx))
     wrf_countries = 67 ! initialize everything to 67 (undefined)
     ALLOCATE(timediff_from_utc(ix,jx))
     IF(retrotno(d).or.naei(d)) THEN
        write(6,*)'calling get_wrfcountries'
        CALL get_wrfcountries_tno()
     ELSEIF(ineris(d)) THEN
        CALL get_wrfcountries_ineris()
     END IF

write(6,*)'Before wrf_grid if'
     IF(.NOT. coarser .OR. d .EQ. 1) THEN
write(6,*)'Just about to call wrf_grid'
        CALL wrf_grid()
write(6,*)'Back from wrf_grid'
     END IF


     !!! REMOVE:
     ! Save London emissions to file
!     IF(d .EQ. 1) THEN
!        CALL save_london()
!     END IF
     !!!

     !!! REMOVE:
     ! Read London emissions from file
!     CALL read_london()
     !!!


     ! Make emissions coarser by specific ratios
     IF(coarser) THEN
        IF(d .EQ. 1) THEN
           ALLOCATE(emi5d_orig(ix,kx,jx,nchem,nsources))
           emi5d_orig = 0.D0
           emi5d_orig = emi5d
!        ELSEIF(d .GT. 1) THEN
        END IF

        CALL make_coarser()

        IF(d .EQ. max_dom) THEN
           DEALLOCATE(emi5d_orig)
        END IF
     END IF

     ! Remove emissions in a selected area
     IF (removing_area) THEN
        CALL remove_area()
     END IF

     IF (wrfchemi) THEN
        ! Now we write all the emissions data to two binary files
        ! that can be read by convert_emiss.exe in WRFV3/chem.
        CALL write_wrfchemi()
     END IF

     IF (binary) THEN
        ! Now we write all the emissions data to two binary files
        ! that can be read by convert_emiss.exe in WRFV3/chem.
        WRITE(*,*) ' - Writing to binary file wrfem_00to12z_d'//id_nr(d)
        CALL write_binary()
     END IF

     IF (ascii) THEN
        ! We also write the data to two ASCII-files so we can check
        ! the output.
        WRITE(*,*) ' - Writing to ASCII file wrfem_00to12z_d'//id_nr//'.txt'
        CALL write_ascii()
     END IF

     IF (plot) THEN
        WRITE(*,*) ' - Writing to ASCII-files for plotting'
        CALL write_plot()
     END IF

     CALL summary()

     ! Deallocate matrices:
     DEALLOCATE(wrf_lats)
     DEALLOCATE(wrf_lons)
     DEALLOCATE(emi5d)
     DEALLOCATE(emi3d)
     DEALLOCATE(city_regions)
     DEALLOCATE(wrf_countries)
     DEALLOCATE(timediff_from_utc)
  END DO ! d = 1, max_dom

END SUBROUTINE main_prog
