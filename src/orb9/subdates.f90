!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                         D A T J M 1                         *      
!  *                                                             *      
!  *        Calcolo della data (giorno, mese, anno, ora)         *      
!  *                dal Tempo Giuliano Modificato                *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
!  INPUT:     TJM     (R*8)  -  Tempo Giuliano Modificato               
!                               (= Tempo Giuliano - 2 400 000.5)        
!                                                                       
!  OUTPUT:    GIORNO  (I*4)  -  Giorno del mese ( 1 <= GIORNO <= 31 )   
!             MESE    (I*4)  -  Mese dell'anno  ( 1 <=  MESE  <= 12 )   
!             ANNO    (I*4)  -  Anno (es: 1987)                         
!             ORA     (R*8)  -  Ora del giorno  ( 0. <= ORA < 24. )     
!             CMESE   (CH*3) -  Mese in caratteri (es: 'DIC')           
!                                                                       
SUBROUTINE datjm1(tjm,giorno,mese,anno,ora,cmese) 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: tjm
  INTEGER, INTENT(OUT) :: giorno,anno, mese 
  CHARACTER*3, INTENT(OUT) :: cmese
  DOUBLE PRECISION, INTENT(OUT) :: ora
! END INTERFACE
  CHARACTER*3 cmv(12) 
  DATA cmv/'GEN','FEB','MAR','APR','MAG','GIU','LUG','AGO','SET',   &
     &         'OTT','NOV','DIC'/     
  DOUBLE PRECISION a
  INTEGER ia, ic, ib, id, ie, ig  
  a=tjm+2400001.d0 
  ia=a 
  ora=(a-ia)*24.d0 
  if(ia.lt.2299161)then 
     ic=ia+1524 
  else 
     ib=FLOOR((ia-1867216.25d0)/36524.25d0) 
     ic=ia+ib-FLOOR(ib/4.d0)+1525 
  end if
  id=FLOOR((ic-122.1d0)/365.25d0) 
  ie=FLOOR(365.25d0*id) 
  ig=FLOOR((ic-ie)/30.6001d0) 
  giorno=ic-ie-FLOOR(30.6001d0*ig) 
  mese=ig-1-12*FLOOR(ig/14.d0) 
  anno=id-4715-FLOOR((7+mese)/10.d0) 
  cmese=cmv(mese) 
END subroutine datjm1
