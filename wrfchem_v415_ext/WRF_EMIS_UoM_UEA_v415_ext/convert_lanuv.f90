! Comment from Oivind Hodnebrog:
! Received this program from Michael Memmesheimer and Christoph Kessler
! on May 26, 2009.

SUBROUTINE G2PHILA(GX,GY,PHI,LAM)
!     *******************************************************************
!     Rechnet Gauss-Krueger-Koordinaten in Laenge und Breite um.
!     change Gauss-Krueger to latitude longitude 
!     programmed probably by Wenzel Bruecher, formerly at
!     Institut Meteorologie und Geophysik University Cologne
!     some translation be christoph kessler May 2009 ck
!     Mail: ch.kessler@t-online.de 
!     more at
!     http://en.wikipedia.org/wiki/Gauss-Kr%C3%BCger_coordinate_system
!     and the much longer german page
!     possibly http://earth-info.nima.mil/GandG/geotrans/
!     is a public domain site for conversion (I did not test it, ck)
!     testdata May 2009
!     gx   4369976.     gy  5540166.    
!     phi-lat  0.8723996     lam-lon 0.1777755     
!     lat  49.98481  lon  10.18578    
!     gx can also have other first digits
!     *******************************************************************
!     ->  GX     = Rechtswert in m 
!     ->  GX     = West/East value in meter e.g. 4369976;  4 is degree 
!         (4x3 = 12 degree from Greenwich) 369976 is distance to latitude 12
!     ->  GY     = Sout/North in meter e.g. 5540166 5 000 000 is offset 
!         to something equatorlike 540km 166 meter is distance to offset
!     <-  PHI    = Breite (Bogenmass) 
!     <-  PHI    = Latitude (rad) 
!     <-  LAM    = Laenge (Bogenmass)
!     <-  LAM    = Longitude (rad)
!     *******************************************************************

  IMPLICIT NONE

  INTEGER KZ
  REAL BF,E0,ETA,F2,F4,F6,LH,KOEFF1,KOEFF2,KOEFF3,KOEFF4, KOEFF5,KOEFF6, &
       NF,RHO,SIGMA,TF,X,Y,phi,lam,gx,gy
     
!     6 Grad Null bei gx=2500000
!     6 Degree Zero at x=2500000
!     kz = 2 als Angabe des GK-Systemes KZ is first digit of west-east
!     number and is multiplied by 3 degrees in LH

  KZ = INT(gx/1000000.)
  X  = gy
  Y  = gx-(KZ*1000000.+ 500000.)
!     also gibt LH den Mittelmeridian (6., 9. Grad etc. an)
!     thus LH indicates middle meridian (6., 9. degrees etc.)
  LH = KZ * 3
!
!     Berechnung der geographischen Breite des Lotfusspunktes
!     Calculation of geo. latitude as computation of latitude of base of
!     perpendicular
!
  RHO = 180/3.141592654
  E0  = 111120.619680
  F2  = 0.143885358/RHO
  F4  = 0.000210790/RHO
  F6  = 0.000000423/RHO
  SIGMA = X/E0/RHO
  BF  = SIGMA+F2*SIN(2*SIGMA)+F4*SIN(4*SIGMA)+F6*SIN(6*SIGMA)
  LH  = LH/RHO
!     Berechnung der KOeffizienten computation of coefficients
!
  TF  = TAN(BF)
  ETA = 0.006719219*((COS(BF))**2)
  NF  = 6398786.849/SQRT(1+ETA)
  KOEFF2 = -1/(2 * NF**2)*TF*(1+ETA)
  KOEFF4 = 1/(24 * NF**4)*TF*(5+3*TF**2+6*ETA*(1-TF**2))
  KOEFF6 = -1/(720*NF**6)*TF*(61+90*TF**2+45*TF**4)
  KOEFF1 = 1/(NF*COS(BF))
  KOEFF3 = 1/(6*NF**3*COS(BF))*(1+2*TF**2+ETA)
  KOEFF5 = 1/(120*NF**5*COS(BF))*(5+28*TF**2+24*TF**4)
!
!     Berechnung von Breite und Laenge computation of lat long
!
  PHI = BF+KOEFF2*Y**2+KOEFF4*Y**4+KOEFF6*Y**6
  LAM = LH+KOEFF1*Y+KOEFF3*Y**3+KOEFF5*Y**5
  
  RETURN
END SUBROUTINE G2PHILA
