SUBROUTINE BNGtolatlonconv(E,N,lat,lon)
! Subroutine to convert British National Grid coordinates into latitude and longitude
! Calculations from Annexe C of 'A guide to coordinate systems in Great Britain' by
! Ordnance Survey. Available online at http://www.ordnancesurvey.co.uk/oswebsite/gps/
! information/coordinatesystemsinfo/guidecontents/
! Note that the datum used is WGS 1934. A transformation will need to be applied if
! working with another datum
! This subroutine must be compiled in double precision
! T. Pugh
! 26/08/10
IMPLICIT NONE

!Set reals to double precision
INTEGER, PARAMETER :: DP=KIND(0.0D0)
INTEGER, PARAMETER :: SP=KIND(0.0)
REAL, PARAMETER :: pi=3.14159265358979323846

! Eastings and Northings to be converted. In single precision for compatibility 
! with rest of program. 
REAL(KIND=SP) :: E,N
! Intermediate parameters
REAL(KIND=DP) :: a,b,e2
REAL(KIND=DP) :: N0,E0,F0,phi0,lambda0
REAL(KIND=DP) :: ns, phi_prime
REAL(KIND=DP) :: M1,M2,M3,M4,M
REAL(KIND=DP) :: v,rho,eta2
REAL(KIND=DP) :: VII,VIII,IX,X,XI,XII,XIIA
REAL(KIND=DP) :: phi,lambda
! Output variables
REAL(KIND=SP) :: lat,lon

!WRITE(6,*)'E and N are',E,N

! Set elipsoid constants for Airy 1830 ellipsoid (used by British National Grid)
a=6377563.396
b=6356256.910
e2=((a**2.0)-(b**2.0))/(a**2.0)

! Set projection constants for British National Grid
N0=-100000.0
E0=400000.0
F0=0.9996012717
phi0=49.0*(pi/180.0) ! 49 degrees N
lambda0=-2.0*(pi/180.0) ! 2 degrees W

! Do the calculations

phi_prime=((N-N0)/(a*F0))+phi0

ns=(a-b)/(a+b)

!Calculate M using a interative procedure until the condition
!(N-N0-M)<0.01mm is satisfied
DO

   M1=(1.0+ns+((5.0/4.0)*(ns**2.0))+((5.0/4.0)*(ns**3.0)))*(phi_prime-phi0)

   M2=((3.0*ns)+(3.0*(ns**2.0))+((21.0/8.0)*(ns**3.0)))*(sin(phi_prime-phi0))*(cos(phi_prime+phi0))

   M3=(((15.0/8.0)*(ns**2.0))+((15.0/8.0)*(ns**3.0)))*(sin(2.0*(phi_prime-phi0)))*(cos(2.0*(phi_prime+phi0)))

   M4=((35.0/24.0)*(ns**3))*(sin(3.0*(phi_prime-phi0)))*(cos(3.0*(phi_prime+phi0)))

   M=b*F0*(M1-M2+M3-M4)

   IF((N-N0-M)<0.01E-3) EXIT !Remember 0.01 is in mm

   phi_prime=((N-N0-M)/(a*F0))+phi_prime

ENDDO

!Compute v, rho and eta**2
v=a*F0*((1.0-(e2*(sin(phi_prime)**2.0)))**(-0.5))

rho=a*F0*(1.0-e2)*((1.0-(e2*(sin(phi_prime)**2.0)))**(-1.5))

eta2=(v/rho)-1.0

!Compute phi and lambda
VII=(tan(phi_prime))/(2.0*rho*v)

VIII=((tan(phi_prime))/(24.0*rho*(v**3.0))) * (5.0+(3.0*(tan(phi_prime)**2.0))+eta2-(9.0*(tan(phi_prime)**2.0)*eta2))

IX=((tan(phi_prime))/(720.0*rho*(v**5.0))) * (61.0+(90.0*(tan(phi_prime)**2.0))+(45.0*(tan(phi_prime)**4.0)))

X=(1.0/(cos(phi_prime)))/v

XI=((1.0/(cos(phi_prime)))/(6.0*(v**3.0))) * ((v/rho)+(2.0*(tan(phi_prime)**2.0)))

XII=((1.0/(cos(phi_prime)))/(120.0*(v**5.0))) * (5.0+(28.0*(tan(phi_prime)**2.0))+(24.0*(tan(phi_prime)**4.0)))

XIIA=((1.0/(cos(phi_prime)))/(5040.0*(v**7.0))) * (61.0+(662.0*(tan(phi_prime)**2.0)) &
+(1320.0*(tan(phi_prime)**4.0))+(720.0*(tan(phi_prime)**6.0)))

phi=phi_prime-(VII*((E-E0)**2.0))+(VIII*((E-E0)**4.0))-(IX*((E-E0)**6.0))
lambda=lambda0+(X*(E-E0))-(XI*((E-E0)**3.0))+(XII*((E-E0)**5.0))-(XIIA*((E-E0)**7.0))

! Convert from radians to degrees
lat=REAL(phi*(180.0/pi))
lon=REAL(lambda*(180.0/pi))

! Write out output
!WRITE(6,*)'Lat and lon are',lat,lon

END SUBROUTINE
