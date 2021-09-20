MODULE massmod
  IMPLICIT NONE
  PUBLIC
! ==================================
! {\bf pardim.h} ORBIT8V
! max dimensions of state vectors
  INTEGER,PARAMETER :: nbox=10
  INTEGER,PARAMETER :: nastx=1100
  INTEGER,PARAMETER :: nvzx=nastx
  INTEGER,PARAMETER :: norbx=nbox-1+2*nastx
  INTEGER,PARAMETER :: nyx=nbox-1+nastx+nvzx
  INTEGER,PARAMETER :: nvarx=nyx*6
  INTEGER,PARAMETER :: nvar2x=nyx*3
  INTEGER,PARAMETER :: nangx=3
! ==============================================
! {\bf commas.h} : pardim.h ORBIT8V
! list of  masses: the active masses are the gm,
! the ones used to compute the position of the sun in
! barycentric coordinates are the rm
! also J2 and relativistic correction constants
  DOUBLE PRECISION, DIMENSION(nbox) :: gm,rm,sm,pmu,rmt,smp
  DOUBLE PRECISION :: cj2,scw
! planet names as read from input
  CHARACTER*9 nompla(nbox)  
! added 9/9/09: Yarkovsky parameters da/dt in au/My
  DOUBLE PRECISION dadt9(nastx)
END MODULE massmod
