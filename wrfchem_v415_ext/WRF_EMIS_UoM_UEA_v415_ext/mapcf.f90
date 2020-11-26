  ! This subroutine is copied and slightly modified from emiss_v03.F:

     SUBROUTINE MAPCF (XI,YJ,XLATI,XLONG)
       USE emission
       IMPLICIT NONE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                     C
!                                                                     C
!     THIS SUBROUTINE COMPUTES THE LATITUDE AND LONGITUDE FROM MODEL  C
!     INDEXES OF A POINT.                                             C
!                                                                     C
!     INPUT :                                                         C
!                                                                     C
!        XI : X COORDINATE OF THE POINT IN MODEL INDEX.               C
!                                                                     C
!        YJ : Y COORDINATE OF THE POINT IN MODEL INDEX.               C
!                                                                     C
!     OUTPUT :                                                        C
!                                                                     C
!        XLAT : LATITUDE OF THE POINT.                                C
!                                                                     C
!        XLON : LONGITUDE OF THE POINT.                               C
!     NOTE ***                                                        C
!                                                                     C
!     A WEST LONGITUDE IS GIVEN BY A NEGATIVE NUMBER; POSITIVE        C
!     ANGLES DENOTE EAST LONGITUDE.                                   C
!                                                                     C
!     A NORTH LATITUDE IS GIVEN BY A POSITIVE NUMBER, AND A NEGATIVE  C
!     NUMBER FOR A SOUTH LATITUDE.                                    C
!                                                                     C
!                                                                     C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     REAL        ::  XN, PSI1, FRTY5D, CONV, RLAT1, RLAT2, POLE
     REAL        ::  PSX, CELL, CELL2, R, C2, XCNTR, PHICTR, YCNTR
     REAL        ::  CNTRJ, CNTRI, XIBD, YJBD, X, Y, FLP, FLPP, RXN
     REAL        ::  CEL1, CEL2
     REAL        ::  XI, YJ, XLATI, XLONG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        XN = 0.
!  9/5/03 Modified for more general Lambert Conformal with 2 standard parallels,CLAT1,CLAT2
        PSI1=CLAT1
        FRTY5D=ATAN(1.)
        CONV=45./FRTY5D
        IF (IPROJ.EQ.1)THEN
!23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
           RLAT1=CLAT1/CONV
           RLAT2=CLAT2/CONV
           IF(CLAT1.NE.CLAT2)THEN
              XN=ALOG(COS(RLAT1)/COS(RLAT2))/ALOG(TAN(FRTY5D+.5*RLAT2)/ &
                   TAN(FRTY5D+.5*RLAT1))
           ELSE
              XN=SIN(RLAT1)
           ENDIF
!           WRITE(*,*)'IPROJ,XN,PSI1,CONV,REKM,DX= ',IPROJ,XN,PSI1,CONV,REKM,DX
        ENDIF
        IF (IPROJ.EQ.2) XN = 1.0
!
!-----PSI1 IS LATITUDE WHERE CONE OR PLANE INTERSECTS EARTH
!
        POLE = 90.
        PSI1 = PSI1/CONV
        IF (XLATC.LT.0.)POLE=-90.
!
!-----REKM IS RADIUS OF EARTH IN KM
!     SUBTRACT XROT DEGREES (modify XCNTR,YCNTR) TO ROTATE LAMBERT CONFORMAL PROJECTION
!     CALCULATE R
!
        IF (IPROJ.NE.3) THEN
           PSX = (POLE-XLATC)/CONV
           IF (IPROJ.EQ.1) THEN
              IF(CLAT1.EQ.CLAT2)THEN
                 CELL = REKM/tan(psi1)
                 CELL2=1.
              ELSE
                 CELL =REKM*SIN(PSI1)/XN
                 CELL2 = (TAN(PSX/2.))/(TAN(PSI1/2.))
              ENDIF
           ENDIF
           IF (IPROJ.EQ.2) THEN
              CELL = REKM*SIN(PSX)/XN
              CELL2 = (1.+COS(PSI1))/(1.+COS(PSX))
           ENDIF
           R = CELL*(CELL2)**XN
           XCNTR = 0.0
           YCNTR = -R
        ENDIF
!
!-----FOR MERCATOR PROJECTION, THE PROJECTION IS TRUE AT LATITUDE PSI1
!
        IF (IPROJ.EQ.3) THEN
           C2 = REKM*COS(PSI1)
           XCNTR = 0.0
           PHICTR = XLATC/CONV
           CELL = (1.0+SIN(PHICTR))/COS(PHICTR)
           YCNTR = C2*ALOG(CELL)
        ENDIF
!     WRITE (6,1010) XCNTR,YCNTR
!1010 FORMAT (1X,'X COORD GRID CNTR = ',F8.1,' Y COORD = ',F8.1)
!
!-----CALCULATE X AND Y POSITIONS OF GRID
!
        CNTRJ = (JLX+1)/2.
        CNTRI = (ILX+1)/2.
        IF(INEST1.EQ.1)THEN
           XIBD=XNESSTR+(XI-1)*DX/DXBIGDO
           YJBD=YNESSTR+(YJ-1)*DX/DXBIGDO
        ELSE
           XIBD=XI
           YJBD=YJ
        ENDIF
        X = XCNTR+(XIBD-CNTRI)*DXBIGDO/1000.
        Y = YCNTR+(YJBD-CNTRJ)*DXBIGDO/1000.
!     WRITE (6,1020) X,Y
!1020 FORMAT (1X,'X INDEX = ',F8.1,' Y INDEX = ',F8.1)
!
!-----NOW CALC LAT AND LONG OF THIS POINT
!
        IF (IPROJ.NE.3) THEN
           IF (XLATC.LT.0.0) THEN
              FLP = ATAN2(X,Y)
           ELSE
              FLP = ATAN2(X,-Y)
           ENDIF
           FLPP = (FLP/XN)*CONV+XLONC
           IF (FLPP.LT.-180.) FLPP = FLPP+360.
           IF (FLPP.GT.180.) FLPP = FLPP-360.
           XLONG = FLPP
!
!-----NOW SOLVE FOR LATITUDE
!
           R = SQRT(X*X+Y*Y)
           IF (XLATC.LT.0.0) R = -R
           IF (IPROJ.EQ.1) THEN
              IF(CLAT1.EQ.CLAT2)THEN
                 XLATI = (1./tan(psi1) + psi1 - r/rekm)*conv
              ELSE
                 CELL = (R*XN)/(REKM*SIN(PSI1))
                 RXN = 1.0/XN
                 CEL1 = TAN(PSI1/2.)*(CELL)**RXN
                 CEL2 = ATAN(CEL1)
                 PSX = 2.*CEL2*CONV
                 XLATI = POLE-PSX
              ENDIF
           ENDIF
           IF (IPROJ.EQ.2) THEN
              CELL = R/REKM
              CEL1 = CELL/(1.0+COS(PSI1))
              CEL2 = ATAN(CEL1)
              PSX = 2.*CEL2*CONV
              XLATI = POLE-PSX
           ENDIF
        ENDIF
!
!     CALCULATIONS FOR MERCATOR LAT, LONG AND MAP SCALES...
!
        IF (IPROJ.EQ.3) THEN
           XLONG = XLONC+((X-XCNTR)/C2)*CONV
           CELL = EXP(Y/C2)
           XLATI = 2.0*(CONV*ATAN(CELL))-90.0
        ENDIF
        RETURN
!
     END SUBROUTINE  MAPCF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
