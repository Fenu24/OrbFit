MODULE planet_masses
IMPLICIT NONE

PUBLIC
! Max number of bodies in binary ephemeris file
INTEGER,PARAMETER :: nbepx=70
! Binary (random access) ephemeris
!
! filbep   -  Name of the binary ephemeris file
! nbep     -  Number of bodies
! masbep   -  Masses
!
CHARACTER*80 filbep
INTEGER nbep
DOUBLE PRECISION, DIMENSION(nbepx) :: masbep
! file for asteroid constants and names
character*80 filbec

! selection of asteroids; mostly to avoid self-perturbations
! note that asteroid are identified by number, since only
! numbered asteroids can have known mass
INTEGER astid(nbepx),iatrue
! strings with asteroid names                                          
character*30 astnam(nbepx) 

! maximum number of planets (from JPL epehemrides)
INTEGER,PARAMETER ::  nplax=10
! maximum number of perturbing bodies (including massive asteroids)
INTEGER,PARAMETER :: nmassx=nplax+nbepx
! mass of sun, ?, of planets as used in integration, sun+earth+moon, earth
DOUBLE PRECISION gm0,gmu,gm(nmassx),gmse,gmearth,gmoon
!  masses of the planets, flags to find them in order
INTEGER npla,nmass,itarg(nplax),listpl(12)
CHARACTER*30,DIMENSION(nmassx) :: ordnam !planet names
!  distance below which close approach is stored
DOUBLE PRECISION,DIMENSION(nmassx) :: dmin

! Close approach control
DOUBLE PRECISION dmea,dmoon,dmjup,dmast,dter

! Dynamical Model
INTEGER ilun,imerc,iplut,irel,icrel,iast,icast,iaber,istat, iclun
PUBLIC iast ! for selast
PUBLIC ilun,imerc,iplut,irel,iaber,istat,iclun ! for rmodel


END MODULE planet_masses
