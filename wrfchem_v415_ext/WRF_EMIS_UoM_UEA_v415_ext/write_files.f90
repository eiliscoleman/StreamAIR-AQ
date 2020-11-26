  SUBROUTINE write_binary()
    USE emission
    IMPLICIT NONE

    ! Now we write all the emissions data to two binary files
    ! that can be read by convert_emiss.exe in WRFV3/chem.
    INTEGER      :: i,k,j,n,s,hour
    REAL         :: old_value,new_value,sum_value

    OPEN(11,FILE='wrfem_00to12z_d'//id_nr(d),FORM='UNFORMATTED')
    WRITE(11) nchem
    WRITE(11) ename
    
    DO hour=1, 24
       IF(hour .EQ. 13) THEN
          CLOSE(11)
          WRITE(*,*) ' - Writing to binary file wrfem_12to24z_d'//id_nr(d)
          OPEN(11,FILE='wrfem_12to24z_d'//id_nr(d),FORM='UNFORMATTED')
          WRITE(11) nchem
          WRITE(11) ename
       END IF
       WRITE(11) hour
       
       DO n=1,nchem
          ! Write out 3-D emission arrays to unformatted file
          DO i=1,ix
             DO k=1,kx
                DO j=1,jx
!                   WRITE(*,FMT='(A,I,A,I,A,I,A,I)') 'n=',n,' i=',i,' j=',j,' k=',k
                   sum_value = 0.
                   DO s=1,nsources
                      old_value = emi5d(i,k,j,n,s)
!                      WRITE(*,FMT='(A,I,A,F20.10)') 's=',s,' old_value=',old_value

                      ! Scale hourly
                      IF(((ineris(d) .OR. retrotno(d)) .AND. scaling_emep_hour) .OR. &
                           (scaling_hour)) THEN
                         CALL scale_hour(old_value,new_value,hour,s,i,j)
                         sum_value = sum_value + new_value
                      ELSE
                         sum_value = sum_value + old_value
                      END IF

                   END DO
                   emi3d(i,k,j) = sum_value
                   
!                   IF((i.EQ.28).AND.(j.EQ.37).AND.(k.EQ.1).AND.(n.EQ.11)) THEN
!                      WRITE(*,FMT='(A,I,A,I,A,I,A,I)') 'n=',n,' i=',i,' j=',j,' k=',k
!                      WRITE(*,FMT='(A,I,A,F20.10)') 'hour=',hour,' emi3d(i,k,j)',emi3d(i,k,j)
!                   END IF
                END DO
             END DO
          END DO
          WRITE(11) emi3d
       END DO

    END DO
    CLOSE(11)
  END SUBROUTINE write_binary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE write_ascii()
    USE emission
    IMPLICIT NONE

    INTEGER      :: i,k,j,n,s,hour  
    REAL         :: old_value,new_value,sum_value

    ! TESTUTSKRIFT TIL ASCII:
    OPEN(15,FILE='wrfem_00to12z_d'//id_nr(d)//'.txt',FORM='FORMATTED')
    WRITE(15,*) nchem
    WRITE(15,*) ename
    
    DO hour=1, 24
       IF(hour .EQ. 13) THEN
          CLOSE(15)
          WRITE(*,*) ' - Writing to ASCII file wrfem_12to24z_d01.txt'
          OPEN(15,FILE='wrfem_12to24z_d'//id_nr(d)//'.txt',FORM='FORMATTED')
          WRITE(15,*) nchem
          WRITE(15,*) ename
       END IF
       WRITE(15,*) hour
       
       DO n=1,nchem
          ! Write out 3-D emission arrays to unformatted file
          DO i=1,ix
             DO k=1,kx
                DO j=1,jx
                   sum_value = 0.
                   DO s=1,nsources
                      old_value = emi5d(i,k,j,n,s)

                      ! Scale hourly
                      IF(((ineris(d) .OR. retrotno(d)) .AND. scaling_emep_hour) .OR. &
                           (scaling_hour)) THEN
                         CALL scale_hour(old_value,new_value,hour,s,i,j)
                         sum_value = sum_value + new_value
                      ELSE
                         sum_value = sum_value + old_value
                      END IF

                   END DO
                   emi3d(i,k,j) = sum_value
                END DO
             END DO
          END DO
          WRITE(15,*) emi3d
       END DO

    END DO
    CLOSE(15)

  END SUBROUTINE write_ascii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Kan bruke gnuplot eller matlab for å plotte filene,
  ! men matlab egner seg nok best.
  SUBROUTINE write_plot()
    USE emission
    IMPLICIT NONE

    CHARACTER (len=80)          :: component_name,filename,file_lat,file_lon
    CHARACTER (len=2)           :: hour_string
    REAL                        :: old_value,new_value,sum_value
    INTEGER                     :: i,k,j,n,s,width,ilun,ilun_lat,ilun_lon, &
         rstat,comp,hour,hour_end
!    INTEGER, PARAMETER          :: n_plots = 17
    INTEGER, PARAMETER          :: n_plots = 3
    INTEGER, DIMENSION(n_plots) :: radm_numbers

!    DATA radm_numbers / 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17 /
    DATA radm_numbers / 2,7,11 /
    ilun = 31
    ilun_lat = 32
    ilun_lon = 33

    IF(scaling_hour) THEN
       hour_end = 24
    ELSE
       hour_end = 1
    END IF

    DO n = 1, n_plots
       component_name = ename(radm_numbers(n))
       width = LEN_TRIM(ename(radm_numbers(n)))

       DO hour = 1, hour_end
          CALL i2a(hour,hour_string)
          IF(hour .LT. 10) THEN
             hour_string = '0'//hour_string
          END IF

          IF(n .EQ. 1 .AND. hour .EQ. 1) THEN
             file_lat = './plotfiles/wrf/wrf_lat_d'//id_nr(d)//'.dat'
             file_lon = './plotfiles/wrf/wrf_lon_d'//id_nr(d)//'.dat'
             WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(file_lat)
             WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(file_lon)
             OPEN(UNIT=ilun_lat,FILE=file_lat,FORM='FORMATTED',IOSTAT=rstat)
             OPEN(UNIT=ilun_lon,FILE=file_lon,FORM='FORMATTED',IOSTAT=rstat)
          END IF

          filename = './plotfiles/wrf/wrf_'//component_name(3:width)//'_'//hour_string//'_d'//id_nr(d)//'.dat'
          WRITE(*,FMT='(A,A)') '   - Writing to file: ', TRIM(filename)
          OPEN(UNIT=ilun,FILE=filename,FORM='FORMATTED',IOSTAT=rstat)
          DO i=1,ix
             DO j=1,jx
                DO k=1,kx
                   comp = radm_numbers(n)
                   sum_value = 0.
                   DO s=1,nsources
                      old_value = emi5d(i,k,j,comp,s)

                      ! Scale hourly
                      IF(((ineris(d) .OR. retrotno(d)) .AND. scaling_emep_hour) .OR. &
                           (scaling_hour)) THEN
                         CALL scale_hour(old_value,new_value,hour,s,i,j)
                         sum_value = sum_value + new_value
                      ELSE
                         sum_value = sum_value + old_value
                      END IF

                   END DO
                   WRITE(UNIT=ilun,FMT='(F20.10)',IOSTAT=rstat) sum_value
                END DO

                IF(n .EQ. 1 .AND. hour .EQ. 1) THEN
!                   CALL MAPCF(FLOAT(i)+0.5,FLOAT(j)+0.5,XLAT,XLON)
                   WRITE(UNIT=ilun_lat,FMT='(F20.10)',IOSTAT=rstat) wrf_lats(i,j)
                   WRITE(UNIT=ilun_lon,FMT='(F20.10)',IOSTAT=rstat) wrf_lons(i,j)
                END IF

             END DO
             WRITE(UNIT=ilun,FMT='(A)',IOSTAT=rstat)
          END DO
          CLOSE(UNIT=ilun)
          IF(n .EQ. 1 .AND. hour .EQ. 1) THEN
             CLOSE(UNIT=ilun_lat)
             CLOSE(UNIT=ilun_lon)
          END IF
       END DO
    END DO

    WRITE(*,*) ' - Done writing to files for plotting'

  END SUBROUTINE write_plot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
