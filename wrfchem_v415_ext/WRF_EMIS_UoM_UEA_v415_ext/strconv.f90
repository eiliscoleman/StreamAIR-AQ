FUNCTION front_trim(buf) RESULT(res)
  IMPLICIT NONE
  !// Arguments
  CHARACTER(LEN=*), INTENT(IN)   :: buf
  INTEGER                        :: res
  !// Local variables
  INTEGER                        :: i, lng
  lng = LEN_TRIM(buf)
  DO i = 1, lng
    IF (buf(i:i) .NE. ' ') THEN
      res = i
      RETURN
    END IF
  END DO
END FUNCTION front_trim

FUNCTION a2i(buf) RESULT(res)
  IMPLICIT NONE
  !// Arguments
  CHARACTER(LEN=*), INTENT(IN)   :: buf
  INTEGER                        :: res
  !// Local variables
  INTEGER                        :: i, foffs, lng
  INTEGER, EXTERNAL              :: front_trim
  LOGICAL                        :: neg
  neg = .FALSE.
  lng = LEN_TRIM(buf)
  foffs = front_trim(buf)
  res = 0
  DO i = foffs, lng
    IF(buf(i:i) .EQ. '-') THEN
      neg = .TRUE.
      CONTINUE
    END IF
    IF(buf(i:i) .GE. '0' .AND. buf(i:i) .LE. '9' .AND. &
      buf(i:i) .NE. ' ') THEN
      res = res * 10
      res = res + (IACHAR(buf(i:i)) - 48)
    END IF
  END DO
  IF(neg) THEN
    res = res * (-1)
  END IF
END FUNCTION a2i

FUNCTION a2f(buf) RESULT(res)
  IMPLICIT NONE
  !// Arguments
  CHARACTER(LEN=*)               :: buf
  REAL                           :: res
  !// Local variables
  INTEGER                        :: i,j,k,itmp, foffs, lng
  INTEGER                        :: p, q
  REAL                           :: rtmp
  INTEGER, EXTERNAL              :: front_trim
  INTEGER, EXTERNAL              :: a2i
  LOGICAL                        :: neg
  LOGICAL                        :: exponential
  neg = .FALSE.
  exponential = .FALSE.
  lng = LEN_TRIM(buf)
  foffs = front_trim(buf)
  k = INDEX(buf,'.')
  p = INDEX(buf,'e')
  q = INDEX(buf,'E')
  IF(p /= -1 .AND. p /= 0) THEN
    exponential = .TRUE.
  END IF
  IF(q /= -1 .AND. q /= 0) THEN
    exponential = .TRUE.
    p = q
  END IF
  IF (k /= -1 .AND. k /= 0) THEN
  !// We have a floating point number
    itmp = 0
    !// Get the integer part of the number
    DO i = foffs, k - 1
      IF(buf(i:i) .EQ. '-') THEN
        neg = .TRUE.
        CONTINUE
      END IF
      IF(buf(i:i) .GE. '0' .AND. buf(i:i) .LE. '9' .AND. &
        buf(i:i) .NE. ' ') THEN
        itmp = itmp * 10
        itmp = itmp + (IACHAR(buf(i:i)) - 48)
      END IF
    END DO
    res = DBLE(itmp)
    rtmp = 0.
    q = 0
    IF(.NOT. exponential) THEN
      !// We do not have an exponential number
      DO i = LEN_TRIM(buf), k+1, -1
        itmp = (IACHAR(buf(i:i)) - 48)
        rtmp = rtmp + FLOAT(itmp)
        rtmp = rtmp / 10.
      END DO
    ELSE
      !// We have an exponential number
      DO i = p-1, k+1, -1
        itmp = (IACHAR(buf(i:i)) - 48)
        rtmp = rtmp + FLOAT(itmp)
        rtmp = rtmp / 10.
      END DO
      q = a2i(buf(p+1:LEN_TRIM(buf)))
    END IF
    res = res + rtmp
    IF(exponential) THEN
      res = res * 10.**q
    END IF
    IF(neg) THEN
      res = res * (-1)
    END IF
  ELSE
  !// We have an integer
    res = REAL(a2i(buf))
  END IF
END FUNCTION a2f

FUNCTION a2d(buf) RESULT(res)
  IMPLICIT NONE
  !// Arguments
  CHARACTER(LEN=*)               :: buf
  DOUBLE PRECISION               :: res
  !// Local variables
  INTEGER                        :: i,j,k,itmp, foffs, lng
  INTEGER                        :: p, q
  DOUBLE PRECISION               :: dtmp
  INTEGER, EXTERNAL              :: front_trim
  INTEGER, EXTERNAL              :: a2i
  LOGICAL                        :: neg
  LOGICAL                        :: exponential
  neg = .FALSE.
  exponential = .FALSE.
  lng = LEN_TRIM(buf)
  foffs = front_trim(buf)
  k = INDEX(buf,'.')
  p = INDEX(buf,'e')
  q = INDEX(buf,'E')
  IF(p /= -1 .AND. p /= 0) THEN
    exponential = .TRUE.
  END IF
  IF(q /= -1 .AND. q /= 0) THEN
    exponential = .TRUE.
    p = q
  END IF
  IF (k /= -1 .AND. k /= 0) THEN
  !// We have a floating point number
    itmp = 0
    !// Get the integer part of the number
    DO i = foffs, k - 1
      IF(buf(i:i) .EQ. '-') THEN
        neg = .TRUE.
        CONTINUE
      END IF
      IF(buf(i:i) .GE. '0' .AND. buf(i:i) .LE. '9' .AND. &
        buf(i:i) .NE. ' ') THEN
        itmp = itmp * 10
        itmp = itmp + (IACHAR(buf(i:i)) - 48)
      END IF
    END DO
    res = DBLE(itmp)
    dtmp = 0.
    q = 0
    IF(.NOT. exponential) THEN
      !// We do not have an exponential number
      DO i = LEN_TRIM(buf), k+1, -1
        itmp = (IACHAR(buf(i:i)) - 48)
        dtmp = dtmp + FLOAT(itmp)
        dtmp = dtmp / 10.
      END DO
    ELSE
      !// We have an exponential number
      DO i = p-1, k+1, -1
        itmp = (IACHAR(buf(i:i)) - 48)
        dtmp = dtmp + FLOAT(itmp)
        dtmp = dtmp / 10.
      END DO
      q = a2i(buf(p+1:LEN_TRIM(buf)))
    END IF
    res = res + dtmp
    IF(exponential) THEN
      res = res * 10.**q
    END IF
    IF(neg) THEN
      res = res * (-1)
    END IF
  ELSE
  !// We have an integer
    res = DBLE(a2i(buf))
  END IF
END FUNCTION a2d

SUBROUTINE i2a(int, res)

   IMPLICIT NONE
   !// Arguments
   INTEGER, INTENT(IN)            :: int
   CHARACTER(LEN=*), INTENT(OUT)  :: res
   !// Local variables
   INTEGER                        :: i,j,k
   CHARACTER(LEN=80)              :: sbuf
   sbuf = ' '
   res = ' '
   i = int
   k = 1
   DO
      j = MOD(i,10)
      sbuf(k:k) = ACHAR(j+48)
      i = i / 10
      IF (i <= 0) THEN
         EXIT
      END IF
      k = k + 1
   END DO
   k = LEN_TRIM(sbuf)
   IF (k>1) THEN
      j = k
      i = 1
      DO i = 1, k
         res(i:i) = sbuf(j:j)
         j = j - 1
      END DO
   ELSE
      res(1:1) = sbuf(1:1)
   END IF

END SUBROUTINE i2a
