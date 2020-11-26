MODULE module_naei
! Modified from the original routine by Øivind Hodnebrog to read direct
! from the data on the NAEI website by by T.Pugh
! 24/08/10
  USE emission

  INTEGER, PARAMETER :: n_iemis = 812, n_jemis = 1378, &
       n_comps = 7, n_sources = 17
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: naei_plot
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: naei_readdata
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: naei_data
  REAL, ALLOCATABLE, DIMENSION(:,:) :: n_longitude
  REAL, ALLOCATABLE, DIMENSION(:,:) :: n_latitude
  LOGICAL            :: naei_covered(n_iemis,n_jemis)
  CHARACTER (len=5), DIMENSION(n_comps):: components_naei, components_naei_l
  INTEGER            :: tno_countries(720,672)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: n_x
  INTEGER, ALLOCATABLE, DIMENSION(:) :: n_y !Now declared here for new lat/lon conversion routine - TP 26/08/10

  DATA components_naei /'CO','NOX','SO2','NH3','VOC','PM10','PM2_5'/
  DATA components_naei_l /'co','nox','so2','nh3','voc','pm10','pm2_5'/ !Lowercase for creating filenames

CONTAINS

  SUBROUTINE allocate_naei()
  ! Allocate arrays for NAEI emissions
  ! Using allocatable arrays to avoid compiler limits on array size
  ! Note that if array is too big for available RAM on the system an out-of-memory error will be generated.
  ! T. Pugh, 18/09/10
    USE emission
    IMPLICIT NONE

    ALLOCATE(naei_plot(n_iemis,n_jemis,n_comps))
    ALLOCATE(naei_readdata(n_iemis,n_jemis,n_comps,nsources))
    ALLOCATE(naei_data(n_iemis,n_jemis,nchem,nsources))
    ALLOCATE(n_longitude(n_iemis,n_jemis))
    ALLOCATE(n_latitude(n_iemis,n_jemis))
    ALLOCATE(n_x(n_iemis*n_jemis))
    ALLOCATE(n_y(n_iemis*n_jemis))

  END SUBROUTINE allocate_naei

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_naei()
    USE emission
    IMPLICIT NONE

    INTEGER            :: i,j,n,c,s,ilun,rstat,ci,cj,ii
    REAL               :: source2d(n_sources,n_jemis,n_iemis)
    REAL               :: source1d(n_sources,(n_jemis*n_iemis))
    REAL               :: sum_source,tot_sum_source
    CHARACTER(len=100) :: filename
    CHARACTER(len=5)   :: component_name
    CHARACTER(len=4)   :: syear
    CHARACTER(len=20), DIMENSION(n_sources) :: source_name
    INTEGER   :: lons,lats,country
    CHARACTER(len=3)   :: country_tno
    ! NAEI datafile prefixes
    DATA source_name /'energyprod','domcom','indcom','indproc','offshore', &
         'solvents','roadtrans','tra2','tra3','tra4','tra5', &
         'othertrans','waste','agric','nature', &
         'totarea','total' /

    ! Now we will read the variables from the NAEI inventory
    PRINT *, ' - Reading NAEI UK emissions from dat files'

    naei_covered  = .false.
    naei_plot     = 0.D0
    naei_readdata     = 0.D0

    ! Copy integer variable year to character variable
    WRITE(syear,'(I4)')year
! first lets read all country indices from whole tno datatbase (even though we are dealing with naei here)

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




    ! We want to read every component
    DO n = 1, n_comps

       DO s = 1,n_sources

          ! Construct filename
          component_name = components_naei(n)
          filename = TRIM(path)//'NAEI/'//TRIM(syear)//'/'//TRIM(component_name)//'/'//&
              TRIM(source_name(s))//TRIM(components_naei_l(n))//'15.asc'
     
          ! Print some info first:
          WRITE(*,FMT='(A,A)') '   - Reading from file: ', TRIM(filename)

          ! Read in from the '.txt'-file. The unit is t/yr !!!!! CHECK UNITS!
          ilun = 21
          OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',STATUS='OLD',IOSTAT=rstat)
          ! If able to open the file then read the data inside
          IF(rstat.EQ.0) THEN
             ! Dump the first 6 lines which contain the header
             DO cj=1,6
                READ(UNIT=ilun,FMT=*)
             ENDDO
             ! Read every point in the inventory:
             DO cj = 1,n_jemis
                IF(component_name.EQ.'CO'.OR.component_name.EQ.'NOX'.OR.component_name.EQ.'SO2' &
                   .OR.component_name.EQ.'PM10'.OR.component_name.EQ.'PM2_5' &
                   .OR.component_name.EQ.'NH3') THEN
!                   if(component_name.EQ.'NOX'.and.s.EQ.7)then
                    if(component_name.EQ.'NOX')then
                   READ(UNIT=ilun,FMT=*) source2d(s,n_jemis+1-cj,1:n_iemis)
                   source2d(s,n_jemis+1-cj,1:n_iemis)=source2d(s,n_jemis+1-cj,1:n_iemis)*1.0 ! this line added NOX for sensitivity runs!!!!!!
                   else
                   READ(UNIT=ilun,FMT=*) source2d(s,n_jemis+1-cj,1:n_iemis) ! Note: flipping north-south for correct orientation - TP 09/09/10
                   end if

                ELSEIF(component_name.EQ.'VOC') THEN
                   READ(UNIT=ilun,FMT=*) source2d(s,n_jemis+1-cj,1:n_iemis)
                   source2d(s,n_jemis+1-cj,1:n_iemis)=source2d(s,n_jemis+1-cj,1:n_iemis)*1.0 ! reduced VOC for defra case study
                   IF(s.EQ.15) THEN ! Set natural source of VOC to zero
                      source2d(15,n_jemis+1-cj,:) = 0.0D0
                   ENDIF
                ENDIF
             ENDDO
          !Not all species have emission from all sources. Therefore set source to zero if unable to open file.
          ELSE
             WRITE(6,*)'    File does not exist for this source. Setting source to zero'
             source2d(s,:,:)=0.0D0
          ENDIF

          !Rearrange into 1D array
          DO ci=1,n_iemis
             DO cj=1,n_jemis
                c=((ci*n_jemis)-n_jemis)+cj
                source1d(s,c)=source2d(s,cj,ci)
                IF(source1d(s,c).lt.0.0) THEN ! Set no-data elements to zero
                   source1d(s,c)=0.0D0
                ENDIF
                n_x(c)=((ci*1000)-500)-50000 ! x-coordinate in m relative to the origin of the British National Grid.
                n_y(c)=((cj*1000)-500)-50000 ! x-coordinate in m relative to the origin of the British National Grid.
             ENDDO
          ENDDO

          CLOSE(UNIT=ilun)
       ENDDO ! n_sources

       ! Calculate point source emissions by subtracting 'total area sources' (16) from 'total emissions' (17)
       ! Put point sources into source 16, and set 17 to zero.
       DO c=1,(n_iemis*n_jemis)
          source1d(16,c)=source1d(17,c)-source1d(16,c)
          source1d(17,c)=0.0D0
       ENDDO

       ! Carry out other processing
       DO c=1,(n_iemis*n_jemis)

          ! Check if inside 81 km x 81 km London square
          IF((n_x(c) .GE. 488500+50000) .AND. (n_x(c) .LE. 568500+50000) .AND. (n_y(c) .GE. 139500+50000) .AND. &
               (n_y(c) .LE. 219500+50000)) THEN ! +50000 because origin of x and y is not same as British National Grid. TP 07/09/10
             IF(london_pointsrc .EQV. .false.) THEN
                ! Set point source to zero
                source1d(16,c) = 0.0D0
!                WRITE(*,FMT='(A,I7,A,I7)') &
!                  'Removing London NAEI point source from: x=',x,', y=',y
             END IF
             IF(london .EQV. .false.) THEN
                ! Set all emission sources to zero
                source1d(:,c) = 0.0D0
!                WRITE(*,FMT='(A,I7,A,I7)') &
!                  'Removing London NAEI emissions from: x=',x,', y=',y
             END IF

             ! If only NOx-emissions in London:
!             IF(component_name .NE. 'NOX') THEN
!                source1d(:,c) = 0.0D0
!             END IF
          END IF

          IF(naei_pointsrc) THEN
             ! Tilbury is a power plant
             IF((n_x(c) .EQ. 566500+50000) .AND. (n_y(c) .EQ. 175500+50000)) THEN
                source1d(1,c) = source1d(16,c)
             ELSE
                ! We split the point sources equally between industry and power
                source1d(1,c) = source1d(1,c) + 0.5*source1d(16,c)
                source1d(3,c) = source1d(3,c) + 0.5*source1d(16,c)
             END IF
          END IF
          source1d(16,c) = 0.0D0
          
          i = ((n_x(c)+500)/1000) + 50 ! Must add 50 to avoid negative indices because emissions grid starts 50km west and south of British National Grid origin
          j = ((n_y(c)+500)/1000) + 50
          naei_covered(i,j) = .true.

          ! Season scaling (pr. month)
          IF(scaling_season) THEN
             DO s=1,nsources ! Note: nsources is specified in module emission rather than locally
                source1d(s,c) = source1d(s,c)*season_fac(n,27,s,month) ! UK is no. 27
             END DO
          END IF

          ! Add the sources from all the categories
          sum_source = 0.0D0
          DO s = 1,16
             sum_source = sum_source + source1d(s,c)
          END DO

          ! Store the sum_source in an array to be used by plot_naei()
          naei_plot(i,j,n) = sum_source

          ! Store the sources in an array
         
          naei_readdata(i,j,n,:) = source1d(1:15,c)


       END DO ! c lines in file
        
    END DO ! n = 1, n_comps


    ! Short summary
    DO n = 1, n_comps
       tot_sum_source = 0.0D0
       DO j = 1, n_jemis
          DO i = 1, n_iemis
             tot_sum_source = tot_sum_source + naei_plot(i,j,n)
          END DO
       END DO
       WRITE(*,FMT='(3A,F15.3,A)') 'NAEI: Sum of ',components_naei(n),' = ',tot_sum_source,' t/yr'
    END DO

    WRITE(*,*) ' - Done reading from input .txt files'

  END SUBROUTINE read_naei

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE latlon_naei()
    USE emission
    IMPLICIT NONE
    ! We want to get latitudes and longitudes for the whole NAEI-grid
    ! Re-written to calculate latitude and longitude directly from British
    ! National Grid coordinates.
    ! Currently finds lat and lon of centre of grid square
    ! T. Pugh
    ! 26/08/10

    INTEGER :: ci,cj,c
    REAL(KIND=KIND(0.0)) :: lat,lon

    WRITE(6,*)' - Calculating latitudes and longitudes for NAEI grid'

    DO ci=1,n_iemis
       DO cj=1,n_jemis
          c=((ci*n_jemis)-n_jemis)+cj
          CALL BNGtolatlonconv(REAL(n_x(c)),REAL(n_y(c)),lat,lon)
          n_longitude(ci,cj)=lon
          n_latitude(ci,cj)=lat
       ENDDO
    ENDDO

    WRITE(*,*) 'Min and max lon',MINVAL(n_longitude), MAXVAL(n_longitude)
    WRITE(*,*) 'Min and max lat',MINVAL(n_latitude), MAXVAL(n_latitude)

    WRITE(*,*) ' - Done calculating latitudes and longitudes'

  END SUBROUTINE latlon_naei

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE plot_naei()
    IMPLICIT NONE
    CHARACTER (len=80)          :: component_name,filename
    CHARACTER(len=4)            :: syear
    CHARACTER (len=2)           :: smonth
    INTEGER                     :: i,j,n,ilun,rstat

    ilun = 51

    IF(year .LE. 2003) THEN
       syear = '2003'
    ELSE
       syear = '2013'
    END IF

    IF(month .LT. 10) THEN
       WRITE(smonth,'(I1)') month
    ELSE
       WRITE(smonth,'(I2)') month
    END IF

    ! Units are in tons/year/km² if we run this routine before convert_units()
    DO n = 1, n_comps
       component_name = components_naei(n)
       IF(scaling_season) THEN
          filename = './plotfiles/naei_'//syear//'_'//TRIM(smonth)//'_'//TRIM(component_name)//'.dat'
       ELSE
          filename = './plotfiles/naei_'//syear//'_'//TRIM(component_name)//'.dat'
       END IF

       WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

       OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',IOSTAT=rstat)
       DO j=1,n_jemis
          DO i=1,n_iemis
             WRITE(UNIT=ilun,FMT='(F12.4)',IOSTAT=rstat) naei_plot(i,j,n)
          END DO
!          WRITE(UNIT=ilun,FMT='(A)',IOSTAT=rstat)
       END DO
       CLOSE(UNIT=ilun)
    END DO

  END SUBROUTINE plot_naei

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE convert_naei()
    IMPLICIT NONE
    ! Convert the gas emissions from tons/yr to mole/(h*km^2):
    ! Have to multiply by 1e6 to convert from tons to g.
    ! Have to divide by 365*24 to convert from pr. year to pr. hour.
    ! Have to divide by gridbox_area to get pr. square km.
    ! Have to divide by molar mass to convert from g to mole.
    !
    ! Same routine for both RADM2 and CBMZ as position of each species
    ! in ename array is the same in both cases - TP 23/09/10

    REAL        :: constant,gridbox_area,mw(4)
    INTEGER     :: naei_to_radm(4),naei_to_crimech(4),naei_to_cbmz(4),naei_to_saprc(4)
    INTEGER     :: i, j, n, s
    DATA mw / 28.0101, 46.0055, 64.054,17.03 / ! CO, NOX (NO2), SO2,NH3
    DATA naei_to_radm / 11, 2, 1,6 /     ! CO, NOx, SO2,NH3
    DATA naei_to_crimech / 1, 2, 4,5 /   ! CO, NOx, SO2,NH3
    DATA naei_to_saprc / 1,2,4,5 /   ! CO, NOx, SO2, NH3
    DATA naei_to_cbmz / 11, 2, 1,6 /   ! CO, NOx, SO2,NH3
    naei_data = 0.d0
    gridbox_area = 1.  ! NAEI data is per sq. km
    constant = 1000000./(365.*24.)
    
    IF(chem_scheme.EQ.1.OR.chem_scheme.EQ.2) THEN
    DO n = 1,n_comps-3 ! Skip the VOCs and PM for now
       DO j=1,n_jemis
          DO i=1,n_iemis
             naei_data(i,j,naei_to_radm(n),:) = &
                  (constant*naei_readdata(i,j,n,1:nsources))/(gridbox_area*mw(n))
          END DO
       END DO
    END DO


   ELSEIF(chem_scheme.EQ.4) THEN    !crimech
     DO n = 1,n_comps-3 ! Skip the VOCs and PM for now 
       DO j=1,n_jemis
          DO i=1,n_iemis
             naei_data(i,j,naei_to_crimech(n),:) = & 
                  (constant*naei_readdata(i,j,n,1:nsources))/(gridbox_area*mw(n))
!             if(n.eq.2)then
!             naei_data(i,j,naei_to_crimech(n),:) = &
!                  (constant*naei_readdata(i,j,n,1:nsources)*0.80)/(gridbox_area*mw(n))   ! 80% of NOx is NO
!             naei_data(i,j,3,:) = &
!                  (constant*naei_readdata(i,j,n,1:nsources)*0.20)/(gridbox_area*mw(n))   ! 20% of NOx is NO2 for now
!             end if 

          END DO
       END DO
    END DO
    
    ELSEIF(chem_scheme.EQ.5) THEN
    DO n = 1,n_comps-3 ! Skip the VOCs and PM for now
       DO j=1,n_jemis
          DO i=1,n_iemis
             naei_data(i,j,naei_to_cbmz(n),:) = &
                  (constant*naei_readdata(i,j,n,1:nsources))/(gridbox_area*mw(n))
          END DO
       END DO
    END DO

   ELSEIF(chem_scheme.EQ.6) THEN    !saprc99
     DO n = 1,n_comps-3 ! Skip the VOCs and PM for now 
       DO j=1,n_jemis
          DO i=1,n_iemis
             naei_data(i,j,naei_to_saprc(n),:) = & 
                  (constant*naei_readdata(i,j,n,1:nsources))/(gridbox_area*mw(n))
!             if(n.eq.2)then
!             naei_data(i,j,naei_to_saprc(n),:) = &
!                  (constant*naei_readdata(i,j,n,1:nsources)*0.80)/(gridbox_area*mw(n))   ! 80% of NOx is NO
!             naei_data(i,j,3,:) = &
!                  (constant*naei_readdata(i,j,n,1:nsources)*0.20)/(gridbox_area*mw(n))   ! 20% of NOx is NO2 for now
!             end if 

          END DO
       END DO
    END DO
   END IF 

  END SUBROUTINE convert_naei

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE vocs_naei()
    ! All the (47*9 = 423) factors are in this file
    USE voc_factors
    IMPLICIT NONE

    ! Divide VOCs into RADM or CBMZ components and convert the units

    REAL              :: mw(nvocs)
    REAL              :: constant,gridbox_area,value,sum_vocs
    INTEGER           :: i,j,n,s,radm_nr(nvocs),cbmz_nr(nvocs)
    INTEGER           :: crimech_nr(nvocs)
    INTEGER           :: saprc_nr(18)
    INTEGER           :: nn
    INTEGER           :: crimech_nr_uk(25)
    INTEGER           :: cbmz_nr_uk(25)
    REAL              :: HFRAC(25),MOLWTHC(25), hfrac_saprc(18), molwt_saprc(18)

       DATA  (HFRAC(I), I=1,25) / 2.470E-02,4.290E-02,3.985E-01, &
        3.150E-02,1.480E-02,3.371E-02,1.200E-02,1.330E-02,1.300E-03, &
        2.092E-03,1.150E-02,1.659E-02,1.155E-02,1.340E-01, &
        1.375E-03,2.020E-02,5.820E-02,8.540E-02,3.550E-03,1.180E-02, &
        4.080E-03,3.240E-03,5.140E-03,4.520E-03, &
        9.260E-03 /  ! THIS IS THE HFRACS FOR THE 25 VOCs USED IN UK SPECIATION
! 
! saprc vocs C2H6,C3H8,ACET ,ALK3 ,ALK4,ALK5 ,ARO1 ,ARO2 ,C2H2 ,C3H6 ,CCHO ,CRES ,ETHENE ,HCHO ,MEK 
!OLE1,OLE2,RCHO 

       DATA  (molwt_saprc(I), I=1,18) / 30.07,44.10,58.08,55.46,75.71,106.64,102.62, &
      114.02,26.04,42.08,44.05,102.87,28.05,30.03, &
      72.11,62.24,79.10,67.32 /
      DATA  (MOLWTHC(I), I=1,25) / 30.0,44.0,58.0,28.0,42.0,56.0, &
        26.0,30.0,44.0,58.0,58.0,72.0,32.0,46.0,60.0,78.0, &
        92.0,106.0,120,120,120,120,120,120,134 /
     DATA (crimech_nr_uk (i), I=1,25) /6,29,26,7,23,30,24,16,17,18,19,20,21,22, &
        31,25,27,28,9,10,11,12,13,14,15/ 
!     DATA (cbmz_nr_uk(i), I=1,25) / 10,7,8,12,7,14,7,4,3,3,17,17,32,33, &
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
    gridbox_area = 1.  ! NAEI data is per sq. km
    constant = 1000000./(365.*24.)


    !Find which index of the components_naei array contains VOCs
    DO i=1,n_comps
        IF(components_naei(i).eq.'VOC') THEN
            nn=i
        ENDIF
    ENDDO


    WRITE(*,*) ' - Dividing VOC species based on sources'

    IF(chem_scheme.EQ.1.or.chem_scheme.EQ.2.or.chem_scheme.EQ.3.or.chem_scheme.EQ.4.or. &
               chem_scheme.EQ.5) then

    DO s=1,nvocsrc
       DO n=1,25
          DO j=1,n_jemis
             DO i=1,n_iemis
                ! Divide VOCs
                value = naei_readdata(i,j,nn,s)*hfrac(n) 
                sum_vocs = sum_vocs + value

                ! Convert units
                value = (constant*value)/(gridbox_area*molwthc(n))

                IF(chem_scheme.EQ.1) THEN !RADM2
!                  write(6,*)'STOP!!!  NOT MAINTENED FOR THIS YET!'
!                  STOP
                   naei_data(i,j,radm_nr(n),s) = &
                        naei_data(i,j,radm_nr(n),s) + value
                ELSEIF(chem_scheme.EQ.2.or.chem_scheme.eq.5) THEN !CBMZ
                   naei_data(i,j,cbmz_nr_uk(n),s) = &
                        naei_data(i,j,cbmz_nr_uk(n),s) + value
                ELSEIF(chem_scheme.EQ.3.or.chem_scheme.eq.4) THEN !CRIMECH
                   naei_data(i,j,crimech_nr_uk(n),s) = & 
                        naei_data(i,j,crimech_nr_uk(n),s) + value
!                   write(6,*)'check naei:',naei_data(i,j,crimech_nr_uk(n),s),s
                ENDIF
             END DO
          END DO
       END DO
    END DO
    
    
                 else   !saprc99 scheme
      DO s=1,13
       DO n=1,18
          DO j=1,n_jemis
             DO i=1,n_iemis
                ! Divide VOCs
                value = naei_readdata(i,j,nn,s)*hfrac_saprc(n) 
                sum_vocs = sum_vocs + value  
                 ! Convert units
                value = (constant*value)/(gridbox_area*molwt_saprc(n))
                naei_data(i,j,saprc_nr(n),s) = & 
                        naei_data(i,j,saprc_nr(n),s) + value
             END DO
          END DO
       END DO
      END DO                
            
    end if    
                        
!          DO n=1,18
!             DO s=1,13  ! nvocsrc
!                ! Divide VOCs
!                write(6,*)'naei_readdata',naei_readdata(i,j,t_comps,s)
!                value = naei_readdata(i,j,5,s)*hfrac_saprc(n)  
!                sum_vocs = sum_vocs + value
!
!                ! Convert units
!                value = (constant*value)/(gridbox_area*molwt_saprc(n))
!                
!
!                   naei_data(i,j,saprc_nr(n),s) = & 
!                        naei_data(i,j,saprc_nr(n),s) + value

          


    WRITE(*,*) ' - Sum VOC before speciation: ',SUM(naei_readdata(:,:,nn,:)), ' tons/yr'
    WRITE(*,*) ' - Sum VOC after speciation:  ',sum_vocs, ' tons/yr'

    WRITE(*,*) ' - Done dividing VOC species'

  END SUBROUTINE vocs_naei

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE pm_naei()
    IMPLICIT NONE
    ! Convert the particle emissions from tons km-1 yr-1 to ug m-3 m s-1
    !
    ! Same routine for both RADM2 and CBMZ as position of each species
    ! in ename array is the same in both cases - TP 23/09/10

    REAL :: s_per_year, ug_per_ton, m2_per_km2
    REAL :: temp
    INTEGER :: naei_to_radm(3),naei_to_crimech(3)
    INTEGER :: naei_to_cbmz_uk(2),naei_to_crimech_uk(2),naei_to_saprc_uk(2)
    INTEGER :: n,nn,j,i,s
    DATA naei_to_radm / 30, 20, 21 /     ! PM10, PM2.5 (<0.01 um), PM2.5 (>0.01 um, accumulation mode)
    DATA naei_to_crimech / 42, 32, 33 /  ! PM10, pm2.5, pm2.5
    DATA naei_to_crimech_uk / 37,32 /  ! PM10, pm2.5
    DATA naei_to_saprc_uk / 46,41 /  ! PM10, pm2.5
    DATA naei_to_cbmz_uk / 28,23 /  ! PM10, pm2.5

    s_per_year = 365.0*86400.0 !Seconds per year
    ug_per_ton = 1.0D+12
    m2_per_km2 = 1.0D+6

    WRITE(*,*)'  - Allocating PM emissions'

    ! For PM2.5. Need to split between Aitken and Accumulation modes.
    ! Split currently based on ratio of total PM2.5 to total PM0.01 in NAEI UK summary statistics
    DO n=1,n_comps
       IF(components_naei(n).EQ.'PM2_5') THEN
          nn=n
       ENDIF
    ENDDO
  
    
   IF(chem_scheme.EQ.1.OR.chem_scheme.EQ.2) THEN 
    DO j=1,n_jemis
       DO i=1,n_iemis
          DO s=1,nsources
          temp = ((naei_readdata(i,j,nn,s)/s_per_year)*ug_per_ton)/m2_per_km2
          naei_data(i,j,naei_to_radm(2),s)=temp*0.2 ! Size range <0.01 um  !!!THESE FACTORS MAY NEED TO BE REVISITED - TP 09/09/10
          naei_data(i,j,naei_to_radm(3),s)=temp*0.8 ! Size range 0.01-2.5 um
          temp=0.0
          ENDDO
       END DO
    END DO
    nn=0

    ! For PM10
    DO n=1,n_comps
       IF(components_naei(n).EQ.'PM10') THEN
          nn=n
          !write(6,*)'component is PM10',nn
       ENDIF
    ENDDO
   
    DO j=1,n_jemis
       DO i=1,n_iemis
          DO s=1,nsources
          temp = ((naei_readdata(i,j,nn,s)/s_per_year)*ug_per_ton)/m2_per_km2
          !Need to subtract PM2.5 total from PM10 total to get PM in the range 2.5 to 10 um
          temp = temp - naei_data(i,j,naei_to_radm(2),s) - naei_data(i,j,naei_to_radm(3),s)
          naei_data(i,j,naei_to_radm(1),s) = temp
          temp=0.0
          ENDDO
       END DO
    END DO
    nn=0

    ELSE IF (chem_scheme.EQ.3) THEN !CRIMECH

    DO j=1,n_jemis
       DO i=1,n_iemis
          DO s=1,nsources
          temp = ((naei_readdata(i,j,nn,s)/s_per_year)*ug_per_ton)/m2_per_km2          
          naei_data(i,j,naei_to_crimech(2),s)=temp*0.2 ! Size range <0.01 um  !!!THESE FACTORS MAY NEED TO BE REVISITED - TP 09/09/10
          naei_data(i,j,naei_to_crimech(3),s)=temp*0.8 ! Size range 0.01-2.5 um
          temp=0.0
          ENDDO
       END DO
    END DO
    nn=0

    ! For PM10
    DO n=1,n_comps
       IF(components_naei(n).EQ.'PM10') THEN
          nn=n
          !write(6,*)'component is PM10',nn
       ENDIF
    ENDDO
   
    DO j=1,n_jemis
       DO i=1,n_iemis
          DO s=1,nsources
          temp = ((naei_readdata(i,j,nn,s)/s_per_year)*ug_per_ton)/m2_per_km2
          !Need to subtract PM2.5 total from PM10 total to get PM in the range 2.5 to 10 um
      temp=temp-naei_data(i,j,naei_to_crimech(2),s)-naei_data(i,j,naei_to_crimech(3),s)
          naei_data(i,j,naei_to_crimech(1),s) = temp
          temp=0.0
          ENDDO
       END DO
    END DO
    nn=0

    ELSE IF (chem_scheme.EQ.4) THEN !CRIMECH


! OK lets do assign NaEI PM2.5 first
    DO n=1,n_comps
       IF(components_naei(n).EQ.'PM_25') THEN
          nn=n
          !write(6,*)'component is PM_25',nn
       ENDIF
    ENDDO


    DO j=1,n_jemis
       DO i=1,n_iemis
          DO s=1,nsources
          temp = ((naei_readdata(i,j,nn,s)/s_per_year)*ug_per_ton)/m2_per_km2        
          naei_data(i,j,naei_to_crimech_uk(2),s)=temp 
          temp=0.0
          ENDDO   
       END DO     
    END DO

!  assign bc_1,OC_DOM,OC_TRA,BC_1,OC_25_10,EC_25_10,OIN to zero in naei (as we dont have them, we will add then in later from TNO emision file merging)
    DO j=1,n_jemis
       DO i=1,n_iemis
          DO s=1,nsources
          naei_data(i,j,33,s)=0.0   ! OC_DOM
          naei_data(i,j,34,s)=0.0   ! OC_TRA
          naei_data(i,j,35,s)=0.0   ! BC_1
          naei_data(i,j,36,s)=0.0   ! EC_1_25
          naei_data(i,j,38,s)=0.0   ! OC_25_10
          naei_data(i,j,39,s)=0.0   ! EC_25_10
          naei_data(i,j,40,s)=0.0   ! OIN_25
          naei_data(i,j,41,s)=0.0   ! OIN_10
          ENDDO
       END DO
    END DO



    ! For PM10    
    DO n=1,n_comps
       IF(components_naei(n).EQ.'PM10') THEN
          nn=n 
          !write(6,*)'component is PM10',nn
       ENDIF      
    ENDDO      
   
    DO j=1,n_jemis
       DO i=1,n_iemis
          DO s=1,nsources
          temp = ((naei_readdata(i,j,nn,s)/s_per_year)*ug_per_ton)/m2_per_km2
!            if(temp.gt.0.0)write(6,*)'pm10=',temp,s,j,i
          !Need to subtract PM2.5 total from PM10 total to get PM in the range 2.5 to 10 um
          temp=temp-naei_data(i,j,naei_to_crimech_uk(2),s)
          naei_data(i,j,naei_to_crimech_uk(1),s) = temp
          temp=0.0
          ENDDO   
       END DO     
    END DO     
    nn=0    
    write(6,*)'finished assigning pm10 data'

    ELSE IF (chem_scheme.EQ.5) THEN !cbmz with UK speciation

! first lets assign PM25

    DO n=1,n_comps
       IF(components_naei(n).EQ.'PM_25') THEN
          nn=n
          !write(6,*)'component is PM_25',nn
       ENDIF
    ENDDO


    DO j=1,n_jemis
       DO i=1,n_iemis
          DO s=1,nsources
          temp = ((naei_readdata(i,j,nn,s)/s_per_year)*ug_per_ton)/m2_per_km2
          naei_data(i,j,naei_to_cbmz_uk(2),s)=temp
          temp=0.0
          ENDDO
       END DO
    END DO
     nn=0

!  assign bc_1,OC_DOM,OC_TRA,BC_1,OC_25_10,EC_25_10, to zero in naei (as we dont have them, we will add then in later from TNO emision file merging)
    DO j=1,n_jemis
       DO i=1,n_iemis
          DO s=1,nsources
          naei_data(i,j,24,s)=0.0   ! OC_DOM
          naei_data(i,j,25,s)=0.0   ! OC_TRA
          naei_data(i,j,26,s)=0.0   ! BC_1
          naei_data(i,j,27,s)=0.0   ! EC_1_25
          naei_data(i,j,29,s)=0.0   ! OC_25_10
          naei_data(i,j,30,s)=0.0   ! EC_25_10
          ENDDO
       END DO
    END DO



    ! For PM10    
    DO n=1,n_comps
       IF(components_naei(n).EQ.'PM10') THEN
          nn=n
          !write(6,*)'component is PM10',nn
       ENDIF
    ENDDO

    DO j=1,n_jemis
       DO i=1,n_iemis
          DO s=1,nsources
          temp = ((naei_readdata(i,j,nn,s)/s_per_year)*ug_per_ton)/m2_per_km2
          !Need to subtract PM2.5 total from PM10 total to get PM in the range 2.5 to 10 um
          temp=temp-naei_data(i,j,naei_to_cbmz_uk(2),s)
          naei_data(i,j,naei_to_cbmz_uk(1),s) = temp
          temp=0.0
          ENDDO
       END DO
    END DO
    nn=0

    ELSE IF (chem_scheme.EQ.6) THEN !SAPRC99 UK

! OK lets do assign NaEI PM2.5 first
    DO n=1,n_comps
       IF(components_naei(n).EQ.'PM_25') THEN
          nn=n
          !write(6,*)'component is PM_25',nn
       ENDIF
    ENDDO


    DO j=1,n_jemis
       DO i=1,n_iemis
          DO s=1,nsources
          temp = ((naei_readdata(i,j,nn,s)/s_per_year)*ug_per_ton)/m2_per_km2        
          naei_data(i,j,naei_to_saprc_uk(2),s)=temp 
          temp=0.0
          ENDDO   
       END DO     
    END DO

!  assign bc_1,OC_DOM,OC_TRA,BC_1,OC_25_10,EC_25_10,OIN to zero in naei (as we dont have them, we will add then in later from TNO emision file merging)
    DO j=1,n_jemis
       DO i=1,n_iemis
          DO s=1,nsources
          naei_data(i,j,42,s)=0.0   ! OC_DOM
          naei_data(i,j,43,s)=0.0   ! OC_TRA
          naei_data(i,j,44,s)=0.0   ! BC_1
          naei_data(i,j,45,s)=0.0   ! EC_1_25
          naei_data(i,j,47,s)=0.0   ! OC_25_10
          naei_data(i,j,48,s)=0.0   ! EC_25_10
          naei_data(i,j,49,s)=0.0   ! OIN_25
          naei_data(i,j,50,s)=0.0   ! OIN_10
          ENDDO
       END DO
    END DO



    ! For PM10    
    DO n=1,n_comps
       IF(components_naei(n).EQ.'PM10') THEN
          nn=n 
          !write(6,*)'component is PM10',nn
       ENDIF      
    ENDDO      
   
    DO j=1,n_jemis
       DO i=1,n_iemis
          DO s=1,nsources
          temp = ((naei_readdata(i,j,nn,s)/s_per_year)*ug_per_ton)/m2_per_km2
!            if(temp.gt.0.0)write(6,*)'pm10=',temp,s,j,i
          !Need to subtract PM2.5 total from PM10 total to get PM in the range 2.5 to 10 um
          temp=temp-naei_data(i,j,naei_to_saprc_uk(2),s)
          naei_data(i,j,naei_to_saprc_uk(1),s) = temp
          temp=0.0
          ENDDO   
       END DO     
    END DO     
    nn=0    
    write(6,*)'finished assigning pm10 data'

    END IF

  END SUBROUTINE pm_naei

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Find nearest NAEI-point from given latitude and longitude
  SUBROUTINE find_naei_index(lon,lat,istart,iend,jstart,jend,i_save,j_save)
    IMPLICIT NONE

    INTEGER :: ifn, jfn, i_save, j_save
    INTEGER :: istart, iend, jstart, jend
    REAL    :: lon, lat, smallest_diff, diff_x, diff_y, diff_sum

    smallest_diff = 10.
    DO jfn = jstart, jend
       DO ifn = istart, iend
          diff_x = lon-n_longitude(ifn,jfn)
          diff_y = lat-n_latitude(ifn,jfn)
          diff_sum = SQRT(diff_x**2+diff_y**2)
          IF(diff_sum .LT. smallest_diff) THEN
             smallest_diff = diff_sum
             ! Save the coordinates of the closest NAEI point
             i_save = ifn
             j_save = jfn
          END IF
       END DO
    END DO
  END SUBROUTINE find_naei_index

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

  SUBROUTINE deallocate_naei()
  ! Deallocate arrays for NAEI emissions which are not required later on in program
  ! T. Pugh, 18/09/10
    USE emission
    IMPLICIT NONE

    DEALLOCATE(naei_plot)
    DEALLOCATE(naei_readdata)
    DEALLOCATE(n_x)
    DEALLOCATE(n_y)

  END SUBROUTINE deallocate_naei

END MODULE module_naei
