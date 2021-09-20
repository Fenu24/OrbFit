MODULE controlmod
  IMPLICIT NONE
  PUBLIC
! ===============================================
! pardiv.h LONG7
!  maximum number of divisors
  INTEGER, PARAMETER :: ndex=100
!==================================
!  control.h
  DOUBLE PRECISION dd(ndex)
!  control for perturbing planets included and order of theory
  INTEGER ipla(8),nitmax,isec,nplin,nplou,nplin0,nplou0
  INTEGER npli4,nplo4,iqcr
  DOUBLE PRECISION fconte,fconti,divrat,divm
! version number
  CHARACTER *20 vers
!  maximum number of nonlinear forced terms                             
  INTEGER, PARAMETER :: nnlx=20 
!  maximum number of terms in gener.lng
  INTEGER, PARAMETER :: igx=1500
!  integer lists for generic function
  INTEGER ::  ig(20,igx),nig
! =======================parameters=====================================
!  maximum number of iterations                                         
  INTEGER, PARAMETER :: maxit=50
!  accuracy of the computation of Laplace coefficients                  
  DOUBLE PRECISION, PARAMETER :: excrit=1.d-6 
! max number of terms in generic Hamiltonian (degree 4)                 
  INTEGER, PARAMETER :: ntx=200 


END MODULE controlmod
