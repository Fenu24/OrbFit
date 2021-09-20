      SUBROUTINE ekcc1(elem,type,xv,gm,dt)
      IMPLICIT NONE

      DOUBLE PRECISION elem(6),xv(6),gm,dt
      CHARACTER*(*) type

      DOUBLE PRECISION eq(6)
      double precision dxde(6,6),ddxde(3,6,6)

      CALL  kepequ(elem,eq)
      CALL  prop2b(0.d0,eq,dt,xv,gm,0,dxde,ddxde)
      RETURN
      END SUBROUTINE ekcc1
