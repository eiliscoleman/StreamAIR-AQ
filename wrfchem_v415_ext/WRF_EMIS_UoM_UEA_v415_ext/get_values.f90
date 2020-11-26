
  SUBROUTINE get_retro(ir,jr,value)
    USE emission
    USE module_retro
    IMPLICIT NONE

    INTEGER                 :: ir, jr, s, n, k, k2
    REAL                    :: value(max_k,nchem,nsources)

    ! Distribute emissions vertically depending on emission source
    value(:,:,:) = 0.d0
    DO s = 1, nsources
       DO n = 1, nchem
          DO k = 1, kx
             DO k2 = 1, emi_k
                value(k,n,s) = value(k,n,s) + retro_data(ir,jr,n,s)* &
                     hfac_wrf_emi(k,k2)*height_fac(s,k2)
             END DO
          END DO
       END DO
    END DO


  END SUBROUTINE get_retro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE get_tno(ir,jr,value)
    USE emission
    USE module_tno
    IMPLICIT NONE

    INTEGER                 :: ir, jr, s, n, k, k2
    REAL                    :: value(max_k,nchem,nsources)
!    write(6,*)'hey nchem=',nchem,nsources,max_k
    ! Distribute emissions vertically depending on emission source

      value(:,:,:) = 0.d0
    DO s = 1, nsources
       DO n = 1, nchem
          DO k = 1, kx
             DO k2 = 1, emi_k
                value(k,n,s) = value(k,n,s) + & 
                     tno_data(ir,jr,n,s)*hfac_wrf_emi(k,k2)*height_fac(s,k2)
!        if(n.eq.41)write(6,*)'whats',hfac_wrf_emi(k,k2),height_fac(s,k2), &
!             value(k,n,s),tno_data(ir,jr,n,s),ir,jr,s,k,k2
             END DO
          END DO
       END DO
    END DO


  END SUBROUTINE get_tno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE get_ineris(ii,ji,value)
    USE emission
    USE module_ineris
    IMPLICIT NONE

    INTEGER                 :: ii, ji, s, n, k, k2
    REAL                    :: value(max_k,nchem,nsources)

    ! Distribute emissions vertically depending on emission source
    value(:,:,:) = 0.d0
    DO s = 1, nsources
       DO n = 1, nchem
          DO k = 1, kx
             DO k2 = 1, emi_k
                value(k,n,s) = value(k,n,s) + ineris_data(ii,ji,n,s)* &
                     hfac_wrf_emi(k,k2)*height_fac(s,k2)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE get_ineris

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE get_naei(in,jn,value)
    USE emission
    USE module_naei
    IMPLICIT NONE

    INTEGER                 :: in, jn, s, n, k, k2
    REAL                    :: value(max_k,nchem,nsources)

    ! Distribute emissions vertically depending on emission source
    value(:,:,:) = 0.d0
    DO s = 1, nsources
       DO n = 1, nchem
          DO k = 1, kx
             DO k2 = 1, emi_k
                value(k,n,s) = value(k,n,s) + naei_data(in,jn,n,s)* &
                     hfac_wrf_emi(k,k2)*height_fac(s,k2)
!                 if(value(k,n,s).gt.0.)write(6,*)'get naei',value(k,n,s),k,n,s
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE get_naei

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE get_lanuv(lon,lat,il,jl,value)
    USE emission
    USE module_lanuv
    IMPLICIT NONE

    INTEGER                 :: il,jl, il2, jl2, s, n, k, k2
    REAL                    :: value(max_k,nchem,nsources), &
         lon, lat, delta_x, delta_y

    ! Distribute emissions vertically depending on emission source
    value(:,:,:) = 0.d0
    DO s = 1, nsources
       DO n = 1, nchem
          DO k = 1, kx
             DO k2 = 1, emi_k
                value(k,n,s) = value(k,n,s) + lanuv_data(il,jl,n,s)* &
                     hfac_wrf_emi(k,k2)*height_fac(s,k2)
             END DO
          END DO
       END DO
    END DO

    ! Then check for additional emissions (il=250:385) - calculate new indices
    CALL find_lanuv_index_dist(lon,lat,l_iemis_short+1,l_iemis,1,l_jemis,il2,jl2,&
         delta_x,delta_y)
 
    ! If less than 0.5 km away in each direction; then add the emission
    IF(ABS(delta_x) .LT. 0.5 .AND. ABS(delta_y) .LT. 0.5) THEN
       DO s = 1, nsources
          DO n = 1, nchem
             DO k = 1, kx
                DO k2 = 1, emi_k
                   value(k,n,s) = value(k,n,s) + lanuv_data(il2,jl2,n,s)* &
                        hfac_wrf_emi(k,k2)*height_fac(s,k2)
                END DO
             END DO
          END DO
       END DO
    END IF


    ! Then add emissions from industry
    DO s = 1, 4 ! Only the first 4 SNAP sectors
       DO n = 1, nchem
          DO k = 1, kx
             value(k,n,s) = value(k,n,s) + lanuv_industry(il,k,jl,n,s)
          END DO
       END DO
    END DO

    ! Again check for additional emissions (il=250:385)
    ! If less than 0.5 km away in each direction; then add the emission
    IF(ABS(delta_x) .LT. 0.5 .AND. ABS(delta_y) .LT. 0.5) THEN
       DO s = 1, 4 ! Only the first 4 SNAP sectors
          DO n = 1, nchem
             DO k = 1, kx
                value(k,n,s) = value(k,n,s) + lanuv_industry(il2,k,jl2,n,s)
             END DO
          END DO
       END DO
    END IF

  END SUBROUTINE get_lanuv

