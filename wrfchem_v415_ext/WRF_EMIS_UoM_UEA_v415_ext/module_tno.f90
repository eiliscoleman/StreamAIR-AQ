
  MODULE module_tno
  USE emission

  INTEGER, PARAMETER :: t_sources = 15, t_iemis = 720, t_jemis = 672,&  
       t_istart = -240, t_iend = 479, t_jstart = 480, t_jend = 1151, &
    t_comps = 16 ! Changed from 44->27 - removed fire emissions
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: tno_plot
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: tno_data
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: tno_readdata
  REAL, ALLOCATABLE, DIMENSION(:,:) :: tno_readdata_ocdom
  REAL, ALLOCATABLE, DIMENSION(:,:)     :: t_longitude
  REAL, ALLOCATABLE, DIMENSION(:,:)     :: t_latitude
  CHARACTER(len=160) :: tnopath
  CHARACTER (len=8), DIMENSION(t_comps):: components_tno
  INTEGER            :: sector(t_sources)
  INTEGER            :: tno_countries(t_iemis,t_jemis)
  DATA sector / 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 / ! From tno to SNAP
  DATA components_tno /'CO','CO','NOX','NOX','SO2','SO2','NH3', &
  'VOC','VOC','PM10','PM25','BC_1','EC_1_25','EC_25_10','OC_25','OC_25_10'/     ! the entries for some are duplicated to account for fact that in the input file there are also shipping entries

  ! Make our own component type:
  TYPE component
     CHARACTER(len=80) :: name
     INTEGER           :: number
     INTEGER           :: nfields
     CHARACTER(len=10), POINTER :: fields(:)
     CHARACTER(len=80) :: file
     INTEGER           :: low_res    ! 1= we have 1442x436 instead of 1442x872
     REAL              :: mw         ! Molecular Weight (g/mol)
  END TYPE component

  TYPE(component), DIMENSION(t_comps) :: components

CONTAINS

  SUBROUTINE allocate_tno()
  ! Allocate arrays for tno emissions
  ! Using allocatable arrays to avoid compiler limits on array size
  ! Note that if array is too big for available RAM on the system an out-of-memory error will be generated.
  ! T. Pugh, 20/09/10
    USE emission
    IMPLICIT NONE

    ALLOCATE(tno_data(t_iemis,t_jemis,nchem,t_sources))
    ALLOCATE(tno_readdata_ocdom(t_iemis,t_jemis))
    ALLOCATE(tno_readdata(t_iemis,t_jemis,t_comps,t_sources))
    ALLOCATE(t_longitude(t_iemis,t_jemis))
    ALLOCATE(t_latitude(t_iemis,t_jemis))

  END SUBROUTINE allocate_tno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE read_tno()
  IMPLICIT NONE

  ! First we declare the variables we need:
  INTEGER            :: i,j,k,c,n,s,source,lons,lats,country
  REAL*8             :: tno_read(t_iemis,t_jemis),tempdata(1442,436),value
  CHARACTER(len=160) :: filename
  CHARACTER(len=4)   :: long_nc,lat_nc,time_nc
  CHARACTER(len=3)   :: country_tno
  ! Initialize variables:
  long_nc = 'lon'
  lat_nc = 'lat'
  time_nc = 'time'
  tno_read = 0.D0
  tno_readdata_ocdom = 0.D0
  tno_readdata = 0.D0
  tno_countries= 67   ! initialize everything to 67 (undefined)
  tnopath= TRIM(path)//'TNO_EMISS/'  !- TP 23/08/10

    open(12,file='countrys.txt')
    i=1
    do while(.true.)
    read(12,13,end=14)lats,lons,country_tno  ! Read lat lon and country name
! Find EMEP country code from TNO country code
       CALL country_tno2emep(country_tno,country)
! Store the country number in an array
       tno_countries(lons,lats) = country

!    country_index(lons,lats)=country_tno
    i=i+1
    end do
14  continue
    close(12)
13  format(i4,1x,i4,1x,a3)

  CALL readcomponents()

  ! Start with reading the variables from the tno inventory:
  PRINT *, ' - Reading tno emissions from nc files'
  PRINT *, ''

  DO n = 1, t_comps      ! Read every netCDF file
     filename = TRIM(tnopath)//components(n)%file  ! Construct filename

     ! Print some info first:
     WRITE(*,FMT='(A,A)') 'Name of component: ', components(n)%name
     WRITE(*,FMT='(A,I3)') ' number:      ', components(n)%number
     WRITE(*,FMT='(A,F7.3)') 'Mole mass:         ', components(n)%mw
     WRITE(*,FMT='(A,A)') 'Filename:          ', filename

     ! We want to get data from all the sources we have chosen:
     DO c = 1, components(n)%nfields

        ! Emissions from ships have lower resolution:
        IF (components(n)%low_res .EQ. 1) THEN
           tempdata = 0.D0
           CALL readnc_2d_2(tempdata,components(n)%fields(c), &
                long_nc,lat_nc,time_nc,filename,1442,436,month)
           ! Convert the 1442x436 grid in to a 1442x872 grid
           DO j = 1, 436
              DO i = 1, 1442
                    DO s = 1, 2
                       tno_read(i,2*(j-1)+s) = tempdata(i,j)
                    END DO
              END DO
           END DO
        ELSE
           CALL readnc_2d_2(tno_read,components(n)%fields(c), &
                long_nc,lat_nc,time_nc,filename,t_iemis,t_jemis,month)
        write(6,*)'done reading nc'
        END IF

        IF(components(n)%fields(c) .EQ. 'pow') THEN
           source = 1
        ELSE IF(components(n)%fields(c) .EQ. 'res') THEN
           source = 2
        ELSE IF(components(n)%fields(c) .EQ. 'inc') THEN
           source = 3
        ELSE IF(components(n)%fields(c) .EQ. 'pei') THEN
           source = 4
        ELSE IF(components(n)%fields(c) .EQ. 'exf') THEN
           source = 5
        ELSE IF(components(n)%fields(c) .EQ. 'sol') THEN
           source = 6
        ELSE IF(components(n)%fields(c) .EQ. 'tra1') THEN
           source = 7
        ELSE IF(components(n)%fields(c) .EQ. 'tra2') THEN
           source = 8
        ELSE IF(components(n)%fields(c) .EQ. 'tra3') THEN
           source = 9
        ELSE IF(components(n)%fields(c) .EQ. 'tra4') THEN
           source = 10
        ELSE IF(components(n)%fields(c) .EQ. 'tra5') THEN
           source = 11
        ELSE IF(components(n)%fields(c) .EQ. 'nrt') THEN
           source = 12
        ELSE IF(components(n)%fields(c) .EQ. 'was') THEN
           source = 13
        ELSE IF(components(n)%fields(c) .EQ. 'agr') THEN
           source = 14
        ELSE IF(components(n)%fields(c) .EQ. 'nat') THEN
           source = 15
 
        END IF
!        if(components(n)%number.eq.4)write(6,*)'tno_read:',tno_read(:,:)
        ! Convert the gas emissions from kg/(s*m^2) to mole/(km^2*h).
!        tno_read(:,:) = (tno_read(:,:)*3.6*1.E12)/components(n)%mw
        ! Store the sum in an array to be used by plot_tno()
!        tno_plot(:,:,components(n)%number) = tno_plot(:,:,components(n)%number) + tno_read(:,:)
        ! Add the emission data from different sources into a big array:
!         write(6,*)'stuff',sector(source),components(n)%number
!         write(6,*)'hey',season_fac(:,:,:,month)       

        ! Season scaling (pr. month)
       ! Go through TNO grid and apply season_factor 
       DO j = 1, 672
          DO i = 1, 720
        IF(scaling_season .AND. tno_countries(i,j) .LE. ncountries) THEN
           tno_read(i,j)=tno_read(i,j)*season_fac(components(n)%number,tno_countries(i,j),sector(source),month)
!          if(tno_countries(i,j).eq.27)write(6,*)'seasonfac',season_fac(components(n)%number,tno_countries(i,j),sector(source),month),month, &
!                           i,j,tno_countries(i,j),sector(source),components(n)%number
        END IF
          END DO
       END DO

        if(components(n)%number.eq.2) then
!        if(components(n)%number.eq.2.or.components(n)%number.eq.5) then  !if components are NOX or VOC lets reduce them by 30% DEFRA
        if(sector(source).eq.7.or.sector(source).eq.8.or.sector(source).eq.9) then
        tno_readdata(:,:,components(n)%number,sector(source)) = &
             tno_readdata(:,:,components(n)%number,sector(source)) + tno_read(:,:)*1.0  !0.7   ! ADDED THIS LINE FOR TRAFFIC NOX FOR SENSITIVITY STUDIES!!!!!
        else   ! OTHER SOURCES OF NOX eg POW not doubled or halved
        tno_readdata(:,:,components(n)%number,sector(source)) = &
             tno_readdata(:,:,components(n)%number,sector(source)) + tno_read(:,:) 
       end if

        else ! not NOx or VOC
        tno_readdata(:,:,components(n)%number,sector(source)) = &
             tno_readdata(:,:,components(n)%number,sector(source)) + tno_read(:,:)
        end if
        tno_read = 0.D0

!        if(filename.eq.'./input/TNO_EMISS/tno_ecoc_ec_25_10.nc'.and.c.eq.2)then   ! get the domestic combustion emissions from oc
!         tno_readdata_ocdom(:,:)=tno_read(:,:)
!        write(6,*)'got bc1 domestic combustion emissions'
!        end if  

     END DO
     WRITE(*,*) '-------------------------------------------------------------'

  END DO ! n = 1, t_comps

  WRITE(*,*) ' - Done reading from nc files'

END SUBROUTINE read_tno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE convert_tno()
    IMPLICIT NONE

    ! Convert the gas emissions from kg/(0.125dx0.0625d)/yr to mole/(h*km^2):
    ! Have to multiply by 1e3 to convert from kg to g.
    ! Have to divide by 365*24 to convert from pr. year to pr. hour.
    ! Have to divide by gridbox_area to get pr. square km.
    ! Have to divide by molar mass to convert from g to mole.
    REAL        :: constant,earth_radius,pi,lat,gridbox_area,mw(4)
    INTEGER tno_to_radm(4),tno_to_cri(4), tno_to_saprc(4)
    INTEGER     :: i,j,n, s
    DATA mw / 28.0101, 46.0055,64.054,17.0 / ! CO, NOx (NO2), SOx (SO2),NH3
    DATA tno_to_radm / 11,2,1,6 /     ! CO, NOx, SO2, NH3
    DATA tno_to_cri / 1,2,4,5 /   ! CO, NOx, SO2, NH3
    DATA tno_to_cri / 1,2,4,5 /   ! CO, NOx, SO2, NH3

    tno_data = 0.d0
    earth_radius = 6371.     ! Earth's radius in km
    pi = 3.14159265         ! Pi
    constant = 1000./(365.*24.)

! Skip the VOCs for now 
       DO j = 1, t_jemis 
          DO i = 1, t_iemis 

          ! Calculate area of 0.125 deg x 0.0625 deg grid box
          lat = t_latitude(i,j)
          gridbox_area = (earth_radius**2)*((0.125*pi/180.)*(0.0625*pi/180.))*COS((lat*pi)/180.)

             IF(chem_scheme.EQ.1.OR.chem_scheme.EQ.2.OR.chem_scheme.EQ.5) THEN
             
! CO add area/point and shipping

                tno_data(i,j,11,:) = &
                    (constant*tno_readdata(i,j,1,1:t_sources))/(gridbox_area*mw(1)) !+ &
                    
! NOX add area/point and shipping

                tno_data(i,j,2,:) = &
                    (constant*tno_readdata(i,j,2,1:t_sources)*1.0)/(gridbox_area*mw(2)) !+ & 80%NO 

! NOX add area/point and shipping

!                tno_data(i,j,3,:) = &
!                    (constant*tno_readdata(i,j,2,1:t_sources)*0.2)/(gridbox_area*mw(2)) !+ &  and 20% NO2

                    
! SO2 add area/point and shipping

                tno_data(i,j,1,:) = &
                    (constant*tno_readdata(i,j,3,1:t_sources))/(gridbox_area*mw(3)) !+ &

! NH3 add area/point emissions        

                tno_data(i,j,6,:) = & 
                    (constant*tno_readdata(i,j,4,1:t_sources))/(gridbox_area*mw(4)) 

             ELSE   
             
! CO add area/point and shipping   

                tno_data(i,j,1,:) = & 
                    (constant*tno_readdata(i,j,1,1:t_sources))/(gridbox_area*mw(1))  !+ &
!                write(6,*)'CO data:',tno_data(i,j,1,:)                    
! NOx add area/point and shipping

                tno_data(i,j,2,:) = & 
                    (constant*tno_readdata(i,j,2,1:t_sources)*1.0)/(gridbox_area*mw(2))  !+ &
! NOx add area/point and shipping

!                tno_data(i,j,3,:) = &
!                    (constant*tno_readdata(i,j,2,1:t_sources)*0.20)/(gridbox_area*mw(2))  !+ &
                    
! SO2 add area/point and shipping                    

                tno_data(i,j,4,:) = & 
                    (constant*tno_readdata(i,j,3,1:t_sources))/(gridbox_area*mw(3))  !+ &

! NH3 add area/point emissions        

                tno_data(i,j,5,:) = & 
                    (constant*tno_readdata(i,j,4,1:t_sources))/(gridbox_area*mw(4))  
                 
            ENDIF   
          END DO  
       END DO  


!    DO n = 1, t_comps-2  ! Skip the VOCs for now
!       DO j = 1, t_jemis
!          DO i = 1, t_iemis
!
!             ! Calculate area of 0.125 deg x 0.0625 deg grid box
!             lat = t_latitude(i,j)
!             gridbox_area = (earth_radius**2)*((0.125*pi/180.)*(0.0625*pi/180.))*COS((lat*pi)/180.)
!             write(6,*)'mw and chemscheme',mw(n),chem_scheme
!             IF(chem_scheme.EQ.1.OR.chem_scheme.EQ.2) THEN
!                tno_data(i,j,tno_to_radm(n),:) = &
!                    (constant*tno_readdata(i,j,n,1:nsources))/(gridbox_area*mw(n))
!             ELSE
!                tno_data(i,j,tno_to_cri(n),:) = & 
!                    (constant*tno_readdata(i,j,n,1:nsources))/(gridbox_area*mw(n))
!            ENDIF
!          END DO
!       END DO
!    END DO
!    write(6,*)'tno_data=',tno_data(:,:,11,:)
  END SUBROUTINE convert_tno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  
  DO n = 1, t_comps
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
  SUBROUTINE vocs_tno()
    ! All the (47*9 = 423) factors are in this file
    USE voc_factors
    IMPLICIT NONE

    ! Divide VOCs into RADM or CBMZ components and convert the units

    REAL              :: mw(nvocs)
    REAL              :: constant,gridbox_area,value,sum_vocs
    REAL              :: earth_radius,pi,lat
    INTEGER           :: i,j,n,s,radm_nr(nvocs),cbmz_nr(nvocs)
    INTEGER           :: crimech_nr(nvocs),crimech_nr_uk(25)
    INTEGER           :: saprc_nr(18)
    INTEGER           :: cbmz_nr_uk(25)
    REAL              :: HFRAC(25),MOLWTHC(25), hfrac_saprc(18), molwt_saprc(18)

       DATA  (HFRAC(I), I=1,25) / 2.470E-02,4.290E-02,3.985E-01, &
        3.150E-02,1.480E-02,3.371E-02,1.200E-02,1.330E-02,1.300E-03, &
        2.092E-03,1.150E-02,1.659E-02,1.155E-02,1.340E-01, &
        1.375E-03,2.020E-02,5.820E-02,8.540E-02,3.550E-03,1.180E-02, &
        4.080E-03,3.240E-03,5.140E-03,4.520E-03, &
        9.260E-03 /
! 
      DATA  (MOLWTHC(I), I=1,25) / 30.0,44.0,58.0,28.0,42.0,56.0, &
        26.0,30.0,44.0,58.0,58.0,72.0,32.0,46.0,60.0,78.0, &
        92.0,106.0,120,120,120,120,120,120,134 /

! saprc vocs C2H6,C3H8,ACET ,ALK3 ,ALK4,ALK5 ,ARO1 ,ARO2 ,C2H2 ,C3H6 ,CCHO ,CRES ,ETHENE ,HCHO ,MEK 
!OLE1,OLE2,RCHO 

      DATA  (molwt_saprc(I), I=1,18) / 30.07,44.10,58.08,55.46,75.71,106.64,102.62, &
      114.02,26.04,42.08,44.05,102.87,28.05,30.03, &
      72.11,62.24,79.10,67.32 /
      
     DATA (crimech_nr_uk (i), I=1,25) /6,29,26,7,23,30,24,16,17,18,19,20,21,22, &
        31,25,27,28,9,10,11,12,13,14,15/ 
!     DATA (cbmz_nr_uk(i), I=1,25) / 10,7,8,12,7,14,7,4,3,3,17,17,21,22, &
!        5,13,15,16,15,15,15,9,9,9,15 /
     DATA (cbmz_nr_uk(i), I=1,25) /10,7,8,12,13,14,7,4,3,3,17,17,21,22, &
    5,15,15,16,16,16,16,16,16,16,16/

     DATA (saprc_nr (i), I=1,18) /6,10,21,12,13,14,16,17,11,15,19,27,7,18,22,8,9,20 / 


      DATA  (hfrac_saprc(I), I=1,18) /0.1498,0.0921,0.0117,0.0767 , &
      0.0982 ,0.0945,0.0718,0.0516,0.0454,0.0342,0.0048,0.0006, &
      0.1455,0.0344,0.0085,0.0293,0.0476,0.0033 /
      
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

    DATA cbmz_nr / &
         33,7,10,7,15,32,12,8,8,17,9,16,7,4,15,13,17,5,9,15,15,7,9,5,16,16,9,17, &
         7,9,13,9,5,7,3,8,7,15,7,14,14,14,7,15,15,9,14 / ! Note all NOx emissions are assumed to be NO

    DATA crimech_nr / &
        22,26,6,29,27,21,7,26,26,19,26,28,26,16,19,25,20,19,26,28,10,22,26,19,28, &
        28,26,19,24,26,23,26,19,22,17,26,22,11,22,30,30,30,22,9,12,26,30 /
    DATA saprc_nr / &
         6,10,21,12,13,14,16,17,11,15,19,27,7,18,22,8,9,20 /

    sum_vocs = 0.
    earth_radius = 6371.     ! Earth's radius in km
    pi = 3.14159265         ! Pi
    constant = 1000./(365.*24.)

    WRITE(*,*) ' - Dividing VOC species based on sources'


    DO j=1,t_jemis
       DO i=1,t_iemis
          ! Calculate area of 0.125 deg x 0.0625 deg grid box
          lat = t_latitude(i,j)
          gridbox_area = (earth_radius**2)*((0.125*pi/180.)*(0.0625*pi/180.))*COS((lat*pi)/180.)
          
          IF(chem_scheme.EQ.1.or.chem_scheme.EQ.2.or.chem_scheme.EQ.3.or.chem_scheme.EQ.4.or. &
               chem_scheme.EQ.5) then
          
          DO n=1,25
             DO s=1,13  ! nvocsrc
                ! Divide VOCs
!                write(6,*)'tno_readdata',tno_readdata(i,j,t_comps,s)
                value = tno_readdata(i,j,5,s)*hfrac(n)  
                sum_vocs = sum_vocs + value
!                if(n.eq.1)write(6,*)'voc number1:',j,i,s,sum_vocs
                ! Convert units
                value = (constant*value)/(gridbox_area*molwthc(n))
                
                IF(chem_scheme.EQ.1) THEN !RADM2
!                  WRITE(6,*)'STOP!!!!!, uk speciation not implemented for this'
!                  STOP
                   tno_data(i,j,radm_nr(n),s) = &
                        tno_data(i,j,radm_nr(n),s) + value
                ELSEIF(chem_scheme.EQ.2) THEN !CBMZ WITH OLD SPECICIATION INHERITED FROM OIVIND
                   tno_data(i,j,cbmz_nr_uk(n),s) = &
                        tno_data(i,j,cbmz_nr_uk(n),s) + value
                ELSEIF(chem_scheme.EQ.3.or.chem_scheme.EQ.4) THEN !CRIMECH ...4 IS UK SPECIATION
                   tno_data(i,j,crimech_nr_uk(n),s) = & 
                        tno_data(i,j,crimech_nr_uk(n),s) + value
                ELSEIF(chem_scheme.EQ.5) THEN !CBMZ WITH UK EMISSION SPECIATION
                   tno_data(i,j,cbmz_nr_uk(n),s) = & 
                        tno_data(i,j,cbmz_nr_uk(n),s) + value 
                ENDIF
                     
            END DO
          END DO
             
             else   !saprc99 scheme
             
          DO n=1,18
             DO s=1,13  ! nvocsrc
                ! Divide VOCs
!                write(6,*)'tno_readdata',tno_readdata(i,j,t_comps,s)
                value = tno_readdata(i,j,5,s)*hfrac_saprc(n)  
                sum_vocs = sum_vocs + value

                ! Convert units
                value = (constant*value)/(gridbox_area*molwt_saprc(n))
                

                   tno_data(i,j,saprc_nr(n),s) = & 
                        tno_data(i,j,saprc_nr(n),s) + value

                     
             END DO
         END DO
             end if

          
       END DO
    END DO

    WRITE(*,*) '  - Sum VOC before speciation: ', SUM(tno_readdata(:,:,5,1:13)), ' tons/yr'
    WRITE(*,*) '  - Sum VOC after speciation:  ', sum_vocs, ' tons/yr'

    WRITE(*,*) ' - Done dividing VOC species'

  END SUBROUTINE vocs_tno


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE latlon_tno()
    IMPLICIT NONE

    ! We want to write lats and longs to file for the whole TNO grid
    INTEGER            :: i,j,ilun,rstat
    REAL               :: lon,lat
    CHARACTER(len=160) :: filename

    PRINT *, ' - Write latitudes and longitudes for TNO grid'

    ! Construct filename
    filename = './plotfiles/tno_latlon.dat'
    WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

    ! Open file for writing
    ilun = 41
    OPEN(ilun,FILE=filename,FORM='FORMATTED',STATUS='UNKNOWN',IOSTAT=rstat)

    ! Loop through TNO grid
    DO j = t_jstart, t_jend
       DO i = t_istart, t_iend
          lon = (i/8.) + (1./16.)
          lat = (j/16.) + (1./32.)

          WRITE(ilun,FMT='(2F12.5)') lat,lon  ! Write to file for plotting
          t_longitude(i-t_istart+1,j-t_jstart+1) = lon
          t_latitude(i-t_istart+1,j-t_jstart+1) = lat
       END DO
    END DO

    CLOSE(ilun)

    WRITE(*,*) ' - Done writing latitudes and longitudes'

  END SUBROUTINE latlon_tno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE plot_tno()
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
!    DO n = 1, nchem
!       component_name = ename(n)
!       filename = './plotfiles/tno_'//TRIM(smonth)//'_'//TRIM(component_name(3:))//'.dat'
!
!       WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)
!
!       OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',IOSTAT=rstat)
!       DO j=1,t_jemis
!          DO i=1,t_iemis
!             WRITE(UNIT=ilun,FMT='(F16.4)',IOSTAT=rstat) tno_plot(i,j,n)
!          END DO
!       END DO
!       CLOSE(UNIT=ilun)
!    END DO

  END SUBROUTINE plot_tno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE pm_tno()
    IMPLICIT NONE
    ! Convert the particle emissions from kilogram km-1 yr-1 to ug m-3 m s-1
    !
    ! Same routine for both RADM2 and CBMZ as position of each species
    ! in ename array is the same in both cases - TP 23/09/10

    REAL :: s_per_year, ug_per_ton, m2_per_km2
    REAL :: temp,temp2
    REAL              :: constant,gridbox_area
    REAL              :: earth_radius,pi,lat
    INTEGER :: tno_to_radm(7),tno_to_cri(9),tno_to_cbmz(9),tno_to_cbmz_uk(8)
    INTEGER :: tno_to_cri_uk(10),tno_to_saprc(9),tno_to_saprc_uk(10)
    INTEGER :: n,nn,j,i,s
    INTEGER :: npm10,npm25i,npm25j,norgi,norgj,neci,necj,norgc,necc
    INTEGER :: npm10_cri, npm25i_cri,npm25j_cri,norgi_cri,norgj_cri,neci_cri,necj_cri
    INTEGER :: norgc_cri, necc_cri
    
    INTEGER :: npm10_saprc, npm25i_saprc,npm25j_saprc,norgi_saprc,norgj_saprc,neci_saprc,necj_saprc
    INTEGER :: norgc_saprc, necc_saprc
    
    INTEGER :: npm10_cri_uk, npm25_cri_uk
    INTEGER :: nocdom_cri_uk,noctra_cri_uk,nbc1_cri_uk,nec125_cri_uk
    INTEGER :: norgc_cri_uk, necc_cri_uk,oin_25_cri_uk,oin_10_cri_uk

    INTEGER :: npm10_saprc_uk, npm25_saprc_uk
    INTEGER :: nocdom_saprc_uk,noctra_saprc_uk,nbc1_saprc_uk,nec125_saprc_uk
    INTEGER :: norgc_saprc_uk, necc_saprc_uk,oin_25_saprc_uk,oin_10_saprc_uk

    INTEGER :: npm10_cbmz_uk, npm25_cbmz_uk
    INTEGER :: nocdom_cbmz_uk,noctra_cbmz_uk,nbc1_cbmz_uk,nec125_cbmz_uk
    INTEGER :: norgc_cbmz_uk, necc_cbmz_uk
    REAL    ::  tpm25,tbc1,tec125,tocdom,toctra,temp5
    INTEGER :: lats,lons
    CHARACTER (LEN=3) :: country_index(720,672),country
    INTEGER :: nnbc1,nnec125,nnec2510,nnoc25,nnoc2510,nnpm25,nnpm10



   DATA tno_to_radm / 30, 20, 21, 26, 27, 28, 29 /     ! PM10, PM2.5 (<0.01 um), PM2.5 (>0.01 um, accumulation mode)
! and e_orgi, eorgj,e_eci,e_ecj
    DATA tno_to_cbmz / 30, 20, 21, 26, 27, 28, 29, 36, 37 /  ! PM10, pm2.5, pm2.5
! last two are e_orgc,e_ecc
    DATA tno_to_cri / 42, 32, 33,38,39,40,41,43,44 /  ! PM10, pm2.5, pm2.5
    DATA tno_to_cri_uk / 32, 33, 34,35,36,37,38,39,40,41 /
    DATA tno_to_saprc / 51, 41, 42,47,48,49,50,52,53 / 
    DATA tno_to_saprc_uk / 41, 42, 43,44,45,46,47,48,49,50 /
    DATA tno_to_cbmz_uk / 23, 24, 25,26,27,28,29,30 /
    earth_radius = 6371.     ! Earth's radius in km
    pi = 3.14159265         ! Pi
    constant = 1000./(365.*24.)
  
    npm10 = 1
    npm25i= 2
    npm25j= 3
    norgi= 4
    norgj= 5
    neci= 6
    necj= 7
    norgc= 8
    necc= 9
    npm10_cri=1
    npm25i_cri=2
    npm25j_cri=3
    norgi_cri=4
    norgj_cri=5
    neci_cri=6
    necj_cri=7
    norgc_cri=8
    necc_cri=9

    npm25_saprc_uk=1
    nocdom_saprc_uk=2
    noctra_saprc_uk=3
    nbc1_saprc_uk=4
    nec125_saprc_uk=5
    npm10_saprc_uk=6
    norgc_saprc_uk=7
    necc_saprc_uk=8
    oin_25_saprc_uk=9
    oin_10_saprc_uk=10
    
    npm25_cri_uk=1
    nocdom_cri_uk=2
    noctra_cri_uk=3
    nbc1_cri_uk=4
    nec125_cri_uk=5
    npm10_cri_uk=6
    norgc_cri_uk=7
    necc_cri_uk=8
    oin_25_cri_uk=9
    oin_10_cri_uk=10

    npm25_cbmz_uk=1
    nocdom_cbmz_uk=2
    noctra_cbmz_uk=3
    nbc1_cbmz_uk=4
    nec125_cbmz_uk=5
    npm10_cbmz_uk=6
    norgc_cbmz_uk=7
    necc_cbmz_uk=8

    npm10_saprc=1
    npm25i_saprc=2
    npm25j_saprc=3
    norgi_saprc=4
    norgj_saprc=5
    neci_saprc=6
    necj_saprc=7
    norgc_saprc=8
    necc_saprc=9

    s_per_year = 365.0*86400.0 !Seconds per year
!+hrm tonnes per year  
!   ug_per_ton = 1.0D+12
!   for TNO kg per year 
    ug_per_ton = 1.0D+9
    m2_per_km2 = 1.0D+6
    

    WRITE(*,*)'  - Allocating PM emissions'

    ! For PM2.5. Need to split between Aitken and Accumulation modes.
    ! Split currently based on ratio of total PM2.5 to total PM0.01 in tno UK summary statistics
    DO n=1,t_comps
       IF(components_tno(n).EQ.'BC_1') THEN
          nnbc1=n-4   ! subtracted the indices referred to shipping files
          write(6,*)'nnbc1=',nnbc1
       ENDIF
       IF(components_tno(n).EQ.'EC_1_25') THEN
          nnec125=n-4  
          write(6,*)'nnec125=',nnec125
       ENDIF   
       IF(components_tno(n).EQ.'EC_25_10') THEN
          nnec2510=n-4  
          write(6,*)'nnec2510=',nnec2510
       ENDIF   
       IF(components_tno(n).EQ.'OC_25') THEN
          nnoc25=n-4  
          write(6,*)'nnoc25=',nnoc25
       ENDIF   
       IF(components_tno(n).EQ.'OC_25_10') THEN
          nnoc2510=n-4  
          write(6,*)'nnoc2510=',nnoc2510
       ENDIF   
       IF(components_tno(n).EQ.'PM25') THEN
          nnpm25=n-4
          write(6,*)'nnpm25=',nnpm25
       ENDIF
       IF(components_tno(n).EQ.'PM10') THEN
          nnpm10=n-4
          write(6,*)'nnpm10=',nnpm10
       ENDIF
    ENDDO
  
    
   IF(chem_scheme.EQ.1) THEN 
    DO j=1,t_jemis
       DO i=1,t_iemis
              ! Calculate area of 0.125 deg x 0.0625 deg grid box

          lat = t_latitude(i,j)
          gridbox_area = (earth_radius**2)*((0.125*pi/180.)*(0.0625*pi/180.))*COS((lat*pi)/180.)
     
          DO s=1,t_sources
! split BC_1 between eci (25%) and ecj accumulation (75%) and add all EC_1_25 to ecj

          temp = ((tno_readdata(i,j,nnbc1,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)           
          tno_data(i,j,tno_to_radm(neci),s)=temp*0.25 ! Size range 0.01-2.5 um
          temp2 = temp*0.75
          temp=((tno_readdata(i,j,nnec125,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_radm(necj),s)=temp+temp2
          
! now lets split oc25 between orgi and orgj

          temp = ((tno_readdata(i,j,nnoc25,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)                  
          tno_data(i,j,tno_to_radm(norgi),s)=temp*0.2 ! Size range 0.01-2.5 um
          tno_data(i,j,tno_to_radm(norgj),s)=temp*0.8 ! Size range <0.01 um
          
! Now lets assign pm25 to pm25i and pm25j

          temp = ((tno_readdata(i,j,nnpm25,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_radm(npm25i),s)=temp*0.2 ! Size range <0.01 um
          tno_data(i,j,tno_to_radm(npm25j),s)=temp*0.8 ! Size range 0.01-2.5 um

! Now lets assign pm10          
          temp = ((tno_readdata(i,j,nnpm10,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_radm(npm10),s)=temp ! all of pm10 emissions


          temp=0.0
          ENDDO
       END DO
    END DO
    nn=0


   ELSE IF (chem_scheme.EQ.2) THEN !CBMZ
!     write(6,*)'whats this',tno_readdata(:,:,nnpm25,:)
    DO j=1,t_jemis
       DO i=1,t_iemis
              ! Calculate area of 0.125 deg x 0.0625 deg grid box

          lat = t_latitude(i,j)
          gridbox_area = (earth_radius**2)*((0.125*pi/180.)* &
               (0.0625*pi/180.))*COS((lat*pi)/180.)

          DO s=1,t_sources
! split BC_1 between eci (25%) and ecj accumulation (75%) and add all EC_1_25 to ecj

          temp = ((tno_readdata(i,j,nnbc1,s)/s_per_year)* &
                 ug_per_ton)/(m2_per_km2*gridbox_area)
!          write(6,*)'eci data',temp,s,nnbc1,neci,gridbox_area
          tno_data(i,j,tno_to_cbmz(neci),s)=temp*0.25 ! Size range 0.01-2.5 um
          temp2 = temp*0.75
          temp=((tno_readdata(i,j,nnec125,s)/s_per_year)* &
               ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cbmz(necj),s)=temp+temp2

! now lets split oc25 between orgi and orgj

          temp = ((tno_readdata(i,j,nnoc25,s)/s_per_year)* &
                 ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cbmz(norgi),s)=temp*0.2 ! Size range 0.01-2.5 um
          tno_data(i,j,tno_to_cbmz(norgj),s)=temp*0.8 ! Size range <0.01 um

! now lets assign oc_25_10 to orgc

          temp = ((tno_readdata(i,j,nnoc2510,s)/s_per_year)* &
                 ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cbmz(norgc),s)=temp

! now lets assign ec_25_10 to ecc

          temp = ((tno_readdata(i,j,nnec2510,s)/s_per_year)* &
                 ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cbmz(necc),s)=temp

! Now lets assign pm25 to pm25i and pm25j

       temp = ((tno_readdata(i,j,nnpm25,s)/s_per_year)* &
              ug_per_ton)/(m2_per_km2*gridbox_area)
       tno_data(i,j,tno_to_cbmz(npm25i),s)=temp*0.2 ! Size range <0.01 um
          tno_data(i,j,tno_to_cbmz(npm25j),s)=temp*0.8 ! Size range 0.01-2.5 um
!              write(6,*)'yesss',tno_data(i,j,tno_to_cbmz(npm25i),s), &
!                    tno_data(i,j,tno_to_cbmz(npm25j),s), tno_to_cbmz(npm25i), &
!                    tno_to_cbmz(npm25j),s,nnpm25,i,j
! Now lets assign pm10
          temp = ((tno_readdata(i,j,nnpm10,s)/s_per_year)* &
                 ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cbmz(npm10),s)=temp ! all of pm10 emissions


          temp=0.0
          temp2=0.0
          ENDDO
       END DO
    END DO
    nn=0


    ELSE IF (chem_scheme.EQ.3) THEN !CRIMECH SCHEME WITH OLD SPECIATION (containing eci,ecj,orgi,orgj,orgc,pm25i,pm25j etc)

    DO j=1,t_jemis
       DO i=1,t_iemis
              ! Calculate area of 0.125 deg x 0.0625 deg grid box

          lat = t_latitude(i,j)
          gridbox_area = (earth_radius**2)*((0.125*pi/180.)*(0.0625*pi/180.))*COS((lat*pi)/180.)

          DO s=1,t_sources
! split BC_1 between eci (25%) and ecj accumulation (75%) and add all EC_1_25 to ecj

          temp = ((tno_readdata(i,j,nnbc1,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cri(neci_cri),s)=temp*0.25 ! Size range 0.01-2.5 um
          temp2 = temp*0.75
          temp=((tno_readdata(i,j,nnec125,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cri(necj_cri),s)=temp+temp2

! now lets split oc25 between orgi and orgj

          temp = ((tno_readdata(i,j,nnoc25,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cri(norgi_cri),s)=temp*0.2 ! Size range 0.01-2.5 um
          tno_data(i,j,tno_to_cri(norgj_cri),s)=temp*0.8 ! Size range <0.01 um

! now lets assign oc_25_10 to orgc

          temp = ((tno_readdata(i,j,nnoc2510,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cri(norgc_cri),s)=temp

! now lets assign ec_25_10 to ecc

          temp = ((tno_readdata(i,j,nnec2510,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cri(necc_cri),s)=temp

! Now lets assign pm25 to pm25i and pm25j

          temp = ((tno_readdata(i,j,nnpm25,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cri(npm25i_cri),s)=temp*0.2 ! Size range <0.01 um
          tno_data(i,j,tno_to_cri(npm25j_cri),s)=temp*0.8 ! Size range 0.01-2.5 um

! Now lets assign pm10
          temp = ((tno_readdata(i,j,nnpm10,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cri(npm10_cri),s)=temp ! all of pm10 emissions

          temp=0.0
          temp2=0.0   
          ENDDO
       END DO
    END DO
    nn=0


    ELSE IF (chem_scheme.EQ.4) THEN !CRIMECH WITH UK SPECIATION FOR RONOCO WORK
    
    DO j=1,t_jemis
       DO i=1,t_iemis
              ! Calculate area of 0.125 deg x 0.0625 deg grid box
              
          lat = t_latitude(i,j)
          gridbox_area = (earth_radius**2)*((0.125*pi/180.)*(0.0625*pi/180.))*COS((lat*pi)/180.)
          
          DO s=1,t_sources
!  get all BC_1
          tbc1=tno_readdata(i,j,nnbc1,s)
          temp = ((tno_readdata(i,j,nnbc1,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cri_uk(nbc1_cri_uk),s)=temp

! now lets get OC_25 emissions from domestic combustion. we are getting domestic comb only which is from
! source s=2. Because its ocdom only lets set all the rest to zero
           
          if(s.eq.2)then
          tocdom=tno_readdata(i,j,nnoc25,s)
          temp = ((tno_readdata(i,j,nnoc25,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cri_uk(nocdom_cri_uk),s)=temp
          else
          tocdom=0.0
          tno_data(i,j,tno_to_cri_uk(nocdom_cri_uk),s)=0.0
          end if 

! now lets get OC_25 emissions from traffic and the rest of remaining categories
! we have alread removed domestic comb so lets not double count. we will then set
! value =0 when s=2
         if(s.ne.2)then  
         toctra=tno_readdata(i,j,nnoc25,s)        
         temp=((tno_readdata(i,j,nnoc25,s)/s_per_year)*ug_per_ton) &
           /(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cri_uk(noctra_cri_uk),s)=temp
!          write(6,*)'tnodata oc traffic:',temp,s
          else
          toctra=0.0
          tno_data(i,j,tno_to_cri_uk(noctra_cri_uk),s)=0.0
          end if 

! now lets assign ec_1_25
         tec125=tno_readdata(i,j,nnec125,s)
         temp=((tno_readdata(i,j,nnec125,s)/s_per_year)*ug_per_ton) &
           /(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cri_uk(nec125_cri_uk),s)=temp

! now lets assign oc_25_10 to orgc

          temp = ((tno_readdata(i,j,nnoc2510,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cri_uk(norgc_cri_uk),s)=temp

! now lets assign ec_25_10 to ecc

          temp = ((tno_readdata(i,j,nnec2510,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cri_uk(necc_cri_uk),s)=temp

! Now lets assign pm25
          temp5=tno_readdata(i,j,nnpm25,s)
          temp = ((tno_readdata(i,j,nnpm25,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cri_uk(npm25_cri_uk),s)=temp 

! Now lets assign pm10 making sure to subtract pm25 so we dont double count
          temp = ((tno_readdata(i,j,nnpm10,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
       tno_data(i,j,tno_to_cri_uk(npm10_cri_uk),s)=temp-tno_data(i,j,tno_to_cri_uk(npm25_cri_uk),s)
         
! diagnostic on other inorganics
          
       !   tno_data(i,j,tno_to_cri_uk(oin_25_cri_uk),s) &
        !      =   tno_data(i,j,tno_to_cri_uk(npm25_cri_uk),s) &
        !         - tno_data(i,j,tno_to_cri_uk(nbc1_cri_uk),s) &
        !         - tno_data(i,j,tno_to_cri_uk(nec125_cri_uk),s) &
         !        - tno_data(i,j,tno_to_cri_uk(nocdom_cri_uk),s) &
          !       - tno_data(i,j,tno_to_cri_uk(noctra_cri_uk),s)
! now oin_10
        !  tno_data(i,j,tno_to_cri_uk(oin_10_cri_uk),s) &
         !    =    tno_data(i,j,tno_to_cri_uk(npm10_cri_uk),s) &
         !        - tno_data(i,j,tno_to_cri_uk(norgc_cri_uk),s) &
          !       - tno_data(i,j,tno_to_cri_uk(necc_cri_uk),s)

!          write(6,*)'whats tno_data for pm10',tno_data(i,j,tno_to_cri_uk(oin_10_cri_uk),s), &
!                     oin_10_cri_uk,tno_to_cri_uk(oin_10_cri_uk),s
!          tpm25=tno_data(i,j,tno_to_cri_uk(npm25_cri_uk),s)
!          tbc1=tno_data(i,j,tno_to_cri_uk(nbc1_cri_uk),s)
!          tec125=tno_data(i,j,tno_to_cri_uk(nec125_cri_uk),s)
!          tocdom=tno_data(i,j,tno_to_cri_uk(nocdom_cri_uk),s)
!          toctra=tno_data(i,j,tno_to_cri_uk(noctra_cri_uk),s)

!          calc_PM25=temp5-tbc1-tec125-tocdom-toctra
          

!           write(6,*)'source:',s,j,i,t_longitude(i,1),t_latitude(1,j),'temp=',temp, &
!                temp5, tbc1+tec125+tocdom+toctra,country_index(i,j)
          temp=0.0
          temp2=0.0
          temp5=0.0
          ENDDO
       END DO
    END DO

		

    nn=0

    ELSE IF (chem_scheme.EQ.5) THEN !CBMZ WITH UK SPECIATION
    
    DO j=1,t_jemis
       DO i=1,t_iemis
              ! Calculate area of 0.125 deg x 0.0625 deg grid box
    
          lat = t_latitude(i,j)
          gridbox_area = (earth_radius**2)*((0.125*pi/180.)*(0.0625*pi/180.))*COS((lat*pi)/180.)
    
          DO s=1,t_sources
!  get all BC_1

          temp = ((tno_readdata(i,j,nnbc1,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cbmz_uk(nbc1_cbmz_uk),s)=temp

! now lets get OC_25 emissions from domestic combustion. we are getting domestic comb only which is from
! source s=2. Because its ocdom only lets set all the rest to zero
    
          if(s.eq.2)then
          temp = ((tno_readdata(i,j,nnoc25,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cbmz_uk(nocdom_cbmz_uk),s)=temp
          else
          tno_data(i,j,tno_to_cbmz_uk(nocdom_cbmz_uk),s)=0.0
          end if

! now lets get OC_25 emissions from traffic and the rest of remaining categories
! we have alread removed domestic comb so lets not double count. we will then set
! value =0 when s=2
         if(s.ne.2)then
         temp=((tno_readdata(i,j,nnoc25,s)/s_per_year)*ug_per_ton) &
           /(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cbmz_uk(noctra_cbmz_uk),s)=temp
!          write(6,*)'tnodata oc traffic:',temp,s
          else
          tno_data(i,j,tno_to_cbmz_uk(noctra_cbmz_uk),s)=0.0
          end if
! now lets assign ec_1_25

         temp=((tno_readdata(i,j,nnec125,s)/s_per_year)*ug_per_ton) &
           /(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cbmz_uk(nec125_cbmz_uk),s)=temp

! now lets assign oc_25_10 to orgc

          temp = ((tno_readdata(i,j,nnoc2510,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cbmz_uk(norgc_cbmz_uk),s)=temp

! now lets assign ec_25_10 to ecc

          temp = ((tno_readdata(i,j,nnec2510,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cbmz_uk(necc_cbmz_uk),s)=temp

! Now lets assign pm25

          temp = ((tno_readdata(i,j,nnpm25,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cbmz_uk(npm25_cbmz_uk),s)=temp

! Now lets assign pm10
          temp = ((tno_readdata(i,j,nnpm10,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_cbmz_uk(npm10_cbmz_uk),s)=temp

          temp=0.0
          temp2=0.0
          ENDDO
       END DO
    END DO
    nn=0


    ELSE IF (chem_scheme.EQ.6) THEN !SAPRC99 UK spec

        DO j=1,t_jemis
       DO i=1,t_iemis
              ! Calculate area of 0.125 deg x 0.0625 deg grid box
              
          lat = t_latitude(i,j)
          gridbox_area = (earth_radius**2)*((0.125*pi/180.)*(0.0625*pi/180.))*COS((lat*pi)/180.)
          
          DO s=1,t_sources
!  get all BC_1
          tbc1=tno_readdata(i,j,nnbc1,s)
          temp = ((tno_readdata(i,j,nnbc1,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_saprc_uk(nbc1_saprc_uk),s)=temp

! now lets get OC_25 emissions from domestic combustion. we are getting domestic comb only which is from
! source s=2. Because its ocdom only lets set all the rest to zero
           
          if(s.eq.2)then
          tocdom=tno_readdata(i,j,nnoc25,s)
          temp = ((tno_readdata(i,j,nnoc25,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_saprc_uk(nocdom_saprc_uk),s)=temp
          else
          tocdom=0.0
          tno_data(i,j,tno_to_saprc_uk(nocdom_saprc_uk),s)=0.0
          end if 

! now lets get OC_25 emissions from traffic and the rest of remaining categories
! we have alread removed domestic comb so lets not double count. we will then set
! value =0 when s=2
         if(s.ne.2)then  
         toctra=tno_readdata(i,j,nnoc25,s)        
         temp=((tno_readdata(i,j,nnoc25,s)/s_per_year)*ug_per_ton) &
           /(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_saprc_uk(noctra_saprc_uk),s)=temp
!          write(6,*)'tnodata oc traffic:',temp,s
          else
          toctra=0.0
          tno_data(i,j,tno_to_saprc_uk(noctra_saprc_uk),s)=0.0
          end if 

! now lets assign ec_1_25
         tec125=tno_readdata(i,j,nnec125,s)
         temp=((tno_readdata(i,j,nnec125,s)/s_per_year)*ug_per_ton) &
           /(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_saprc_uk(nec125_saprc_uk),s)=temp

! now lets assign oc_25_10 to orgc

          temp = ((tno_readdata(i,j,nnoc2510,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_saprc_uk(norgc_saprc_uk),s)=temp

! now lets assign ec_25_10 to ecc

          temp = ((tno_readdata(i,j,nnec2510,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_saprc_uk(necc_saprc_uk),s)=temp

! Now lets assign pm25
          temp5=tno_readdata(i,j,nnpm25,s)
          temp = ((tno_readdata(i,j,nnpm25,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
          tno_data(i,j,tno_to_saprc_uk(npm25_saprc_uk),s)=temp 

! Now lets assign pm10 making sure to subtract pm25 so we dont double count
          temp = ((tno_readdata(i,j,nnpm10,s)/s_per_year)*ug_per_ton)/(m2_per_km2*gridbox_area)
       tno_data(i,j,tno_to_saprc_uk(npm10_saprc_uk),s)=temp-tno_data(i,j,tno_to_saprc_uk(npm25_saprc_uk),s)
         
! diagnostic on other inorganics
          
          tno_data(i,j,tno_to_saprc_uk(oin_25_saprc_uk),s) &
              =   tno_data(i,j,tno_to_saprc_uk(npm25_saprc_uk),s) &
                 - tno_data(i,j,tno_to_saprc_uk(nbc1_saprc_uk),s) &
                 - tno_data(i,j,tno_to_saprc_uk(nec125_saprc_uk),s) &
                 - tno_data(i,j,tno_to_saprc_uk(nocdom_saprc_uk),s) &
                 - tno_data(i,j,tno_to_saprc_uk(noctra_saprc_uk),s)
! now oin_10
          tno_data(i,j,tno_to_saprc_uk(oin_10_saprc_uk),s) &
             =    tno_data(i,j,tno_to_saprc_uk(npm10_saprc_uk),s) &
                 - tno_data(i,j,tno_to_saprc_uk(norgc_saprc_uk),s) &
                 - tno_data(i,j,tno_to_saprc_uk(necc_saprc_uk),s)

!          write(6,*)'whats tno_data for pm10',tno_data(i,j,tno_to_saprc_uk(oin_10_saprc_uk),s), &
!                     oin_10_saprc_uk,tno_to_saprc_uk(oin_10_saprc_uk),s
!          tpm25=tno_data(i,j,tno_to_saprc_uk(npm25_saprc_uk),s)
!          tbc1=tno_data(i,j,tno_to_saprc_uk(nbc1_saprc_uk),s)
!          tec125=tno_data(i,j,tno_to_saprc_uk(nec125_saprc_uk),s)
!          tocdom=tno_data(i,j,tno_to_saprc_uk(nocdom_saprc_uk),s)
!          toctra=tno_data(i,j,tno_to_saprc_uk(noctra_saprc_uk),s)

!          calc_PM25=temp5-tbc1-tec125-tocdom-toctra
          

!           write(6,*)'source:',s,j,i,t_longitude(i,1),t_latitude(1,j),'temp=',temp, &
!                temp5, tbc1+tec125+tocdom+toctra,country_index(i,j)
          temp=0.0
          temp2=0.0
          temp5=0.0
          ENDDO
       END DO
    END DO

    nn=0
    END IF
!  write(6,*)'whats this',tno_data(:,:,6,:)
  END SUBROUTINE pm_tno


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Find nearest TNO-point from given latitude and longitude
  SUBROUTINE find_tno_index(lon,lat,istart,iend,jstart,jend,i_save,j_save)
    IMPLICIT NONE

    INTEGER :: ifn, jfn, i_save, j_save
    INTEGER :: istart, iend, jstart, jend
    REAL    :: lon, lat, smallest_diff, diff_x, diff_y, diff_sum

    smallest_diff = 10.
    DO jfn = jstart, jend
       DO ifn = istart, iend
          diff_x = lon-t_longitude(ifn,jfn)
          diff_y = lat-t_latitude(ifn,jfn)
          diff_sum = SQRT(diff_x**2+diff_y**2)
          IF(diff_sum .LT. smallest_diff) THEN
             smallest_diff = diff_sum
             ! Save the coordinates of the closest TNO point
             i_save = ifn
             j_save = jfn
          END IF
       END DO
    END DO
  END SUBROUTINE find_tno_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE country_tno2emep(country_tno,country)
    IMPLICIT NONE

    CHARACTER(len=3), INTENT(IN) :: country_tno
    INTEGER, INTENT(INOUT)         :: country

    ! Find EMEP country code from TNO country code
    IF (country_tno .EQ. 'ALB') THEN
       country =        1
    ELSE IF (country_tno .EQ. 'AUT') THEN
       country =    2
    ELSE IF (country_tno .EQ. 'BEL') THEN
       country =        3
    ELSE IF (country_tno .EQ. 'BGR') THEN
       country =        4
    ELSE IF (country_tno .EQ. 'DNK') THEN
       country =        6
    ELSE IF(country_tno .EQ. 'FIN') THEN
       country =        7
    ELSE IF(country_tno .EQ. 'FRA') THEN
       country =        8
    ELSE IF(country_tno .EQ. 'GRC') THEN
       country =        11
    ELSE IF(country_tno .EQ. 'HUN') THEN
       country =        12
    ELSE IF(country_tno .EQ. 'ISL') THEN
       country =        13
    ELSE IF(country_tno .EQ. 'IRL') THEN
       country =        14
    ELSE IF(country_tno .EQ. 'ITA') THEN
       country =        15
    ELSE IF(country_tno .EQ. 'LUX') THEN
       country =        16
    ELSE IF(country_tno .EQ. 'NLD') THEN
       country =        17
    ELSE IF(country_tno .EQ. 'NOR') THEN
       country =        18
    ELSE IF(country_tno .EQ. 'POL') THEN
       country =        19
    ELSE IF(country_tno .EQ. 'PRT') THEN
       country =        20
!    ELSE IF(country_tno .EQ. 'ROM') THEN 
    ELSE IF(country_tno .EQ. 'ROU') THEN
       country =        21
    ELSE IF(country_tno .EQ. 'ESP') THEN
       country =        22
    ELSE IF(country_tno .EQ. 'SWE') THEN
       country =        23
    ELSE IF(country_tno .EQ. 'CHE') THEN
       country =        24
    ELSE IF(country_tno .EQ. 'TUR') THEN
       country =        25
    ELSE IF(country_tno .EQ. 'RUS') THEN
       country =        26
    ELSE IF(country_tno .EQ. 'GBR') THEN
       country =        27
    ELSE IF(country_tno .EQ. 'BLR') THEN
       country =        39
    ELSE IF(country_tno .EQ. 'UKR') THEN
       country =        40
    ELSE IF(country_tno .EQ. 'MDA') THEN
       country =        41
    ELSE IF(country_tno .EQ. 'EST') THEN
       country =        43
    ELSE IF(country_tno .EQ. 'LVA') THEN
       country =        44
    ELSE IF(country_tno .EQ. 'LTU') THEN
       country =        45
    ELSE IF(country_tno .EQ. 'CZE') THEN
       country =        46
    ELSE IF(country_tno .EQ. 'SVK') THEN
       country =        47
    ELSE IF(country_tno .EQ. 'SVN') THEN
       country =        48
    ELSE IF(country_tno .EQ. 'HRV') THEN
       country =        49
    ELSE IF(country_tno .EQ. 'BIH') THEN
       country =        50
    ELSE IF(country_tno .EQ. 'GEO') THEN
       country =        54
    ELSE IF(country_tno .EQ. 'CYP') THEN
       country =        55
    ELSE IF(country_tno .EQ. 'ARM') THEN
       country =        56
    ELSE IF(country_tno .EQ. 'MLT') THEN
       country =        57
    ELSE IF(country_tno .EQ. 'LIE') THEN
       country =        59
    ELSE IF(country_tno .EQ. 'DEU') THEN
!       country =     60
       country =     10
    ELSE IF(country_tno .EQ. 'MCO') THEN
       country =        62
    ELSE IF(country_tno .EQ. 'KGZ') THEN
       country =        68
    ELSE IF(country_tno .EQ. 'AZE') THEN
       country =        69
    ELSE IF(country_tno .EQ. 'YUG') THEN
       country =        72
    ELSE IF(country_tno .EQ. 'KAZ') THEN
       country =        92
    ELSE IF(country_tno .EQ. 'TKM') THEN
       country =        95
    ELSE ! 67 = undefined
       country = 67
    ENDIF

  END SUBROUTINE country_tno2emep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE deallocate_tno()
  ! Deallocate arrays for tno emissions which are not required later on in program
  ! T. Pugh, 18/09/10
    USE emission
    IMPLICIT NONE

!    DEALLOCATE(tno_plot)

  END SUBROUTINE deallocate_tno

END MODULE module_tno
