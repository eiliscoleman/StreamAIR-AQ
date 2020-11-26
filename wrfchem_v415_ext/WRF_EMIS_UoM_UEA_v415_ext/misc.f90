
  SUBROUTINE define_regions()
    USE emission
    IMPLICIT NONE

    ! Define megacity / hot spot regions
    WRITE(*,*) '- Defining megacity / hot spot region'

    city_regions = .false.

    IF(city_region .EQ. 'london_t1') THEN
       WRITE(*,*) ' - Defining region: LONDON - T1'
       city_regions(82:90,37:45) = .true.

    ELSEIF(city_region .EQ. 'london_t2') THEN
       WRITE(*,*) ' - Defining region: LONDON - T2'
       city_regions(127:135,82:90) = .true.

!    IF(city_region .EQ. 'london') THEN
!       WRITE(*,*) ' - Defining region: LONDON'
!       city_regions(82:90,46:54) = .true.

    ELSE IF(city_region .EQ. 'ruhr_t1') THEN
       WRITE(*,*) ' - Defining region: RUHR - T1'
       city_regions(55:72,91:108) = .true.

    ELSE IF(city_region .EQ. 'ruhr_t2') THEN
       WRITE(*,*) ' - Defining region: RUHR - T2'
       city_regions(163:180,82:99) = .true.

! Cairo is located at approximately 30.05 deg N, 31.25 deg E
    ELSE IF(city_region .EQ. 'cairo_t1') THEN
       WRITE(*,*) ' - Defining region: CAIRO - T1'
       city_regions(73:81,181:189) = .true.

    ELSE IF(city_region .EQ. 'cairo_t2') THEN
       WRITE(*,*) ' - Defining region: CAIRO - T2'
       city_regions(118:126,181:189) = .true.

!    ELSE IF(city_region .EQ. 'cairo') THEN
!       WRITE(*,*) ' - Defining region: CAIRO'
!       city_regions(127:135,127:135) = .true.

    ELSE
       WRITE(*,*) ' - No region defined!'

    END IF

  END SUBROUTINE define_regions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE make_coarser()
    USE emission
    IMPLICIT NONE

    ! Make nested emissions as coarse as specified by coarser_ratio

    INTEGER     :: i,k,j,n,s,coarser_ix,coarser_jx,j_coarse,i_coarse,i_fine_start, &
         j_fine_start,i_fine_end,j_fine_end,n_squares
    REAL        :: value

    WRITE(*,*) ' - Making the background emission resolution coarser by ratio ',coarser_ratio

!    emi5d = 0.D0        ! Erase old values, if any

    coarser_ix = ix/coarser_ratio
    coarser_jx = jx/coarser_ratio

    DO j_coarse = 1, coarser_jx
       DO i_coarse = 1, coarser_ix
          i_fine_start = (i_coarse-1)*coarser_ratio + 1
          j_fine_start = (j_coarse-1)*coarser_ratio + 1

          i_fine_end = i_fine_start + coarser_ratio - 1
          j_fine_end = j_fine_start + coarser_ratio - 1        

          ! If not within megacity / hot spot region
          IF(city_regions(i_fine_start,j_fine_start) .EQV. .false.) THEN

             DO k = 1, kx
                DO n=1,nchem
                   DO s=1,nsources
                      value = 0.
                      n_squares = 0

                      ! Find average value in coarser_ratio*coarser_ratio area
                      DO j = j_fine_start, j_fine_end
                         DO i = i_fine_start, i_fine_end
                            value = value + emi5d_orig(i,k,j,n,s)
                            n_squares = n_squares + 1
                         END DO
                      END DO
                   
                      IF(n_squares .EQ. 0) THEN
                         WRITE(*,*) 'WARNING: Zero-division during interpolation!'
                      END IF

                      ! Calculate average:
                      value = value/n_squares

                      ! Replace old value by average value
                      DO j = j_fine_start, j_fine_end
                         DO i = i_fine_start, i_fine_end
                            emi5d(i,k,j,n,s) = value
                         END DO
                      END DO

                   END DO
                END DO
             END DO

          END IF

       END DO
    END DO


    IF(coarser_city .AND. coarser_city_ratio .NE. 1) THEN
       WRITE(*,*) ' - Making the city region emission resolution coarser by ratio ',coarser_city_ratio
    
       coarser_ix = ix/coarser_city_ratio
       coarser_jx = jx/coarser_city_ratio

       DO j_coarse = 1, coarser_jx
          DO i_coarse = 1, coarser_ix
             i_fine_start = (i_coarse-1)*coarser_city_ratio + 1
             j_fine_start = (j_coarse-1)*coarser_city_ratio + 1

             i_fine_end = i_fine_start + coarser_city_ratio - 1
             j_fine_end = j_fine_start + coarser_city_ratio - 1        

             ! If within megacity / hot spot region
             IF(city_regions(i_fine_start,j_fine_start)) THEN

                DO k = 1, kx
                   DO n=1,nchem
                      DO s=1,nsources
                         value = 0.
                         n_squares = 0

                         ! Find average value in coarser_city_ratio*coarser_city_ratio area
                         DO j = j_fine_start, j_fine_end
                            DO i = i_fine_start, i_fine_end
                               value = value + emi5d_orig(i,k,j,n,s)
                               n_squares = n_squares + 1
                            END DO
                         END DO
                   
                         IF(n_squares .EQ. 0) THEN
                            WRITE(*,*) 'WARNING: Zero-division during interpolation!'
                         END IF

                         ! Calculate average:
                         value = value/n_squares

                         ! Replace old value by average value
                         DO j = j_fine_start, j_fine_end
                            DO i = i_fine_start, i_fine_end
                               emi5d(i,k,j,n,s) = value
                            END DO
                         END DO

                      END DO
                   END DO
                END DO

             END IF

             END DO
          END DO

       END IF

  END SUBROUTINE make_coarser

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE remove_area()
    USE emission
    IMPLICIT NONE

    REAL :: sum_co, sum_nox, sum_voc
    INTEGER :: i,j

    ! Remove emissions in a defined megacity / hot spot region

    WRITE(*,*) ' - Removing emissions from region: '//city_region

    DO j = 1, jx
       DO i = 1, ix
          IF(city_regions(i,j)) THEN
             emi5d(i,:,j,:,:) = 0.d0
          END IF
       END DO
    END DO

    WRITE(*,*) ' - Done removing emissions from megacity / hot spot region'

  END SUBROUTINE remove_area

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE summary()
    USE emission
    IMPLICIT NONE

    INTEGER :: n
    CHARACTER(len=9) :: comp
    REAL             :: sum_emis,average_emis

    WRITE(*,*)
    WRITE(*,*) 'Summary for domain ',id_nr(d)
    WRITE(*,*) '*********************'
    WRITE(*,*) 'Note: Values are without diurnal scaling.'
    DO n = 1, nchem
       sum_emis = SUM(emi5d(:,:,:,n,:))
       comp = ename(n)
       average_emis = sum_emis/(ix*jx)
       WRITE(*,*) 'Average emissions of ',comp(3:),' = ',average_emis, &
            ' mole/(km^2 hr)'
    END DO
    WRITE(*,*)

  END SUBROUTINE summary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! REMOVE:
  SUBROUTINE save_london
    USE emission
    IMPLICIT NONE

    INTEGER                 :: ilun
    CHARACTER(len=160)      :: filename
    CHARACTER (len=2)       :: smonth
    CHARACTER (len=4)       :: syear
    ilun = 60

    WRITE(syear,'(I4)') year
    IF(month .LT. 10) THEN
       WRITE(smonth,'(I1)') month
    ELSE
       WRITE(smonth,'(I2)') month
    END IF

    filename = TRIM(path)//'london_'//syear//TRIM(smonth)//'.dat'  ! Construct filename
    WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)

    ! Choose to save an area covering 9 times 81x81 km2 squares
    OPEN(UNIT=ilun,FILE=filename,FORM='UNFORMATTED')
! Time period 1:
!    WRITE(UNIT=ilun) emi5d(73:99,:,28:54,:,:)
! Time period 2:
    WRITE(UNIT=ilun) emi5d(118:144,:,73:99,:,:)
!    WRITE(UNIT=ilun) emi5d(73:99,:,37:63,:,:)
    CLOSE(UNIT=ilun)

  END SUBROUTINE save_london
!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! REMOVE:
  SUBROUTINE read_london
    USE emission
    IMPLICIT NONE

    INTEGER                 :: ilun
    CHARACTER(len=160)      :: filename
    CHARACTER (len=2)       :: smonth
    CHARACTER (len=4)       :: syear
    ilun = 60

    WRITE(syear,'(I4)') year
    IF(month .LT. 10) THEN
       WRITE(smonth,'(I1)') month
    ELSE
       WRITE(smonth,'(I2)') month
    END IF

    filename = TRIM(path)//'london_'//syear//TRIM(smonth)//'.dat'  ! Construct filename
    WRITE(*,FMT='(A,A)') '   - Reading from file: ', TRIM(filename)

    ! Choose to save an area covering 9 times 81x81 km2 squares
    OPEN(UNIT=ilun,FILE=filename,FORM='UNFORMATTED')
! Time period 1:
!    READ(UNIT=ilun) emi5d(64:90,:,172:198,:,:)    ! CHOOSE WHERE TO PUT LONDON EMISSIONS!!
! Time period 2:
    READ(UNIT=ilun) emi5d(109:135,:,172:198,:,:)    ! CHOOSE WHERE TO PUT LONDON EMISSIONS!!
!    READ(UNIT=ilun) emi5d(118:144,:,118:144,:,:)    ! CHOOSE WHERE TO PUT LONDON EMISSIONS!!
    CLOSE(UNIT=ilun)

  END SUBROUTINE read_london
!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
