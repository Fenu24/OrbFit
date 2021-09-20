MODULE virtual_impactor
  USE orbit_elements
  USE tp_trace
  USE dyn_param
  IMPLICIT NONE
  PRIVATE

!  INTEGER, PARAMETER, PUBLIC :: nprevx=5

! object oriented definition
  TYPE virt_imp
     TYPE (orbit_elem) :: ele ! orbital elements of one represntative of the VI
     TYPE(orb_uncert)  :: unc ! covariance and normal matrices
     TYPE (tp_point) :: tp     ! TP point and auxiliary variables
     TYPE(dyn_par) :: dyn    ! nongravitational parameters
!     TYPE (tp_point) :: prev_tp(nprevx) ! TP point on previous encounters
  END TYPE virt_imp

  PUBLIC virt_imp

! storage area
  INTEGER, PARAMETER :: nvix=100
  TYPE(virt_imp), DIMENSION(nvix) :: vis ! list of VIs found so far
  INTEGER num_vi
! individual VI used by riskchecktp as temporary storage
  TYPE(virt_imp) :: curr_vi 
  LOGICAL :: stored_vi

  PUBLIC nvix, vis, num_vi,curr_vi, stored_vi

  PUBLIC store_vi, vi_draw

CONTAINS

  SUBROUTINE store_vi()
    IF(.not.stored_vi) RETURN
    IF(num_vi.lt.nvix)THEN
       num_vi=num_vi+1
       vis(num_vi)=curr_vi
    ELSE
       WRITE(*,*)' store_vi: too many VI, increase max from ', nvix 
    ENDIF
  END SUBROUTINE store_vi

  SUBROUTINE vi_draw(del)
    USE fund_const
    USE planet_masses, ONLY: gmearth
    USE propag_state , ONLY: pro_ele
    INTEGER, PARAMETER :: nj=20
    DOUBLE PRECISION, INTENT(OUT) :: del(2,nj)
!==========correspondent ellipse in orbital elements space=============== 
    DOUBLE PRECISION, DIMENSION(2,2) :: dtpdt,dtpdi
    DOUBLE PRECISION, DIMENSION(2) :: dtpc,dd,ddd
    TYPE(tp_point)   :: tp_curr
    DOUBLE PRECISION dtpde0(2,6), dee(6,6)
    DOUBLE PRECISION :: v_infi,v_infty        ! vel at infinity      
    DOUBLE PRECISION tafter,tbefore           ! time of the final elements
    TYPE(orbit_elem) elkep, elcor, el1
    TYPE(orb_uncert) unckep, unc1
    INTEGER fail_flag,j, nit
    DOUBLE PRECISION U,mu,b_e,dthe, det, rindex
    INTEGER, PARAMETER :: nitmax=5
    DOUBLE PRECISION, PARAMETER :: epsnew=1.d-6
! compute impact cross section
    U=curr_vi%tp%opik%coord(1)
    mu=gmearth/reau**3
    b_e=sqrt(1.d0+(2.d0*mu)/U**2)
! convert to keplerian
    CALL coo_cha(curr_vi%ele,'KEP',elkep,fail_flag,dee)
    CALL convertunc(curr_vi%unc,dee,unckep)
! derivative of TP w.r. to keplerian elements
! assuming jc=1, will kill us if multiple minima
    dtpde0=MATMUL(dtpde(1:2,1:6,1),dee)
! select omega, mean_an
    dtpdt=dtpde0(1:2,5:6)
    det=dtpdt(1,1)*dtpdt(2,2)-dtpdt(1,2)*dtpdt(2,1)
    dtpdi(1,1)=dtpdt(2,2)/det
    dtpdi(2,2)=dtpdt(1,1)/det
    dtpdi(1,2)=-dtpdt(1,2)/det
    dtpdi(2,1)=-dtpdt(2,1)/det    
! loop on polygon
    elcor=elkep
    dthe=dpig/nj
    DO j=1,nj
       tp_curr=curr_vi%tp
       dd(1)=b_e*cos(j*dthe)
       dd(2)=b_e*sin(j*dthe)
       DO nit=1,nitmax
          dtpc=dd-tp_curr%opik%coord(4:5)
          del(1:2,j) =MATMUL(dtpdi,dd)*reau+elkep%coord(5:6)
          elcor%coord(5:6)=del(1:2,j)
! time up to which to search for close approaches                       
          v_infi=v_infty(elcor) 
          CALL aftclov(tp_curr%iplam,elkep%t,tp_curr%tcla,v_infi,tbefore,tafter) 
! make new covariance matrix available for target plane analysis        
          CALL cov_avai(unckep,elkep%coo,elkep%coord) 
! reset storage of close encounters                                     
          njc=0 
          iplam=0 
! propagation to search for new target plane point                      
          CALL pro_ele(elcor,tafter,el1,unckep,unc1)
         ! check presence of target plane                                        
          IF(njc.eq.0.or.iplam.ne.3)THEN 
             WRITE(*,*)' fold found, no TP'
             del=1.d50 
          ELSEIF(tcla(njc).lt.tbefore.or.tcla(1).gt.tafter)THEN 
             WRITE(*,*)' fold found, no TP in the required time span'
             del=1.d50
          ENDIF
! get data on TP analysys of the new close approach
          rindex=tp_curr%rindex                     
          CALL arrloadtp(tp_curr,rindex) 
          ddd=tp_curr%opik%coord(4:5)-dd ! by how much we missed
          IF(sqrt(ddd(1)**2+ddd(2)**2).lt.epsnew) EXIT
       ENDDO
    ENDDO

  END SUBROUTINE vi_draw
END MODULE virtual_impactor
