subroutine do_calc(npt,temper,press,nu0,radius,gamma,fac,omega, &
     vortlambda, ov, inialpha, inibeta, alpha, beta)

  USE general
  USE energies
  USE text
  USE spectra
  USE glob
  USE velocities

  IMPLICIT NONE

  INTEGER :: npt
  REAL (KIND=dp), DIMENSION(0:npt) :: inialpha,inibeta,alpha,beta
  REAL (KIND=dp) :: temper,press,nu0,radius,gamma,fac,omega

  nmax = npt
  t = temper
  p = press
  r = radius

  
