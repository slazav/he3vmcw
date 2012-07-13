MODULE glob

  USE general

  REAL (KIND=dp), SAVE :: t,p,h,r,lo

  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: evrp,evfp,evzp
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: evrm,evfm,evzm
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: elrp,elfp,elzp
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: elrm,elfm,elzm
  REAL (KIND=dp), DIMENSION(0:maxnpt-1), SAVE :: ewm,ewp

END MODULE glob


