! conversion of observations from new data type to old multi vector
    SUBROUTINE obs_convert(obs,obsw,m,aln,den,tau,tut,idsta,iobs,objid,smag, &
    & rmsa,rmsd,rmsmag,sel)
    USE astrometric_observations
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: m ! number of observations
! new data types
    TYPE(ast_obs),INTENT(IN),DIMENSION(m) :: obs
    TYPE(ast_wbsr),INTENT(IN),DIMENSION(m) :: obsw
! observations both arcs: alpha, delta, time (ET and UT), station code, 
    DOUBLE PRECISION, INTENT(OUT) ::  aln(m),den(m),tau(m),tut(m) 
    INTEGER, INTENT(OUT) ::  idsta(m),iobs(m) 
! asteroid identifier, apparent magnitude                               
    CHARACTER*9, INTENT(OUT) :: objid(m) 
    CHARACTER*6, INTENT(OUT) :: smag(m) 
! observation rms, magnitude rms,weight                                 
    DOUBLE PRECISION, INTENT(OUT) :: rmsd(m),rmsa(m),rmsmag(m) 
! selection flags
    INTEGER, INTENT(OUT) :: sel(m) 
! =====================================
    INTEGER i

    aln(1:m)=obs(1:m)%coord(1)
    den(1:m)=obs(1:m)%coord(2)
    tau(1:m)=obs(1:m)%time_tdt
    tut(1:m)=obs(1:m)%time_utc
    idsta(1:m)=obs(1:m)%obscod_i
    objid(1:m)=obs(1:m)%objdes
    smag(1:m)=obs(1:m)%mag_str   
    rmsa(1:m)=obsw(1:m)%rms_coord(1)
    rmsd(1:m)=obsw(1:m)%rms_coord(2)
    rmsmag(1:m)=obsw(1:m)%rms_mag
    sel(1:m)=obsw(1:m)%sel_coord
! iobs as a function of obs%type, obs%tech

! this is wrong just to try
    DO i=1,m
      IF(obs(i)%type.eq.'O')THEN
         iobs(i)=1000
      ELSEIF(obs(i)%type.eq.'R')THEN
         iobs(i)=2002
      ELSEIF(obs(i)%type.eq.'V')THEN
         iobs(i)=2003
      ELSEIF(obs(i)%type.eq.'S')THEN
         iobs(i)=1000
      ENDIF
    ENDDO

    RETURN
    END subroutine obs_convert
