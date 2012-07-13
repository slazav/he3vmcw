MODULE velocities

USE general
USE glob
USE energies

IMPLICIT NONE

CONTAINS

  SUBROUTINE clusterprofile(r,omega,ov)
    ! Vortex-free velocity profile
    IMPLICIT NONE
    INTEGER :: i
    REAL (KIND=dp) :: r,omega,ov,rr
    DO i=0,nmax-1
       rr=(i+(3-SQRT(3.))/6)*dx*r
       evrm(i)=0._dp
       evfm(i)=0._dp
       evzm(i)=0._dp
       elrm(i)=0._dp
       elfm(i)=0._dp
       elzm(i)=1._dp
       ewm(i)=2*omega
       IF (rr > r*SQRT(ov/omega)) THEN
          evfm(i)=omega*rr-ov*r*r/rr
          ewm(i)=0._dp
       END IF
       rr=(i+(3+SQRT(3.))/6)*dx*r
       evrp(i)=0._dp
       evfp(i)=0._dp
       evzp(i)=0._dp
       elrp(i)=0._dp
       elfp(i)=0._dp
       elzp(i)=1._dp
       ewp(i)=2*omega
       IF (rr > r*SQRT(ov/omega)) THEN
          evfp(i)=omega*rr-ov*r*r/rr
          ewp(i)=0._dp
       END IF
    END DO
  END SUBROUTINE clusterprofile

  SUBROUTINE uniformvortcluster(r,omega,ov)
    ! uniform vortex cluster
    IMPLICIT NONE
    INTEGER :: i
    REAL (KIND=dp) :: r,omega,ov,rr
    DO i=0,nmax-1
       rr=(i+(3-SQRT(3.))/6)*dx*r
       evrm(i)=0._dp
       evfm(i)=(omega-ov)*rr
       evzm(i)=0._dp
       elrm(i)=0._dp
       elfm(i)=0._dp
       elzm(i)=1._dp
       ewm(i)=2*ov
       rr=(i+(3+SQRT(3.))/6)*dx*r
       evrp(i)=0._dp
       evfp(i)=(omega-ov)*rr
       evzp(i)=0._dp
       elrp(i)=0._dp
       elfp(i)=0._dp
       elzp(i)=1._dp
       ewp(i)=2*ov
    END DO
  END SUBROUTINE uniformvortcluster


  SUBROUTINE twistedstate(r,omega,kr)
    ! Twisted-state velocity profile
    IMPLICIT NONE
    INTEGER :: i
    REAL (KIND=dp) :: r,omega,kr,rr
    DO i=0,nmax-1
       rr=(i+(3-SQRT(3.))/6)*dx*r
       evrm(i)=0._dp
       evfm(i)=omega*rr-(kr**2/LOG(1+kr**2))*omega*rr/(1+(kr*rr/r)**2)
       evzm(i)=-omega*r*(kr/(LOG(1+kr**2)*(1+(kr*rr/r)**2))-1/kr)
       elrm(i)=0._dp
       elfm(i)=(kr*rr/r)/SQRT(1+(kr*rr/r)**2)
       elzm(i)=1._dp/SQRT(1+(kr*rr/r)**2)
       ewm(i)=2*omega*(kr**2/LOG(1+kr**2))*(1+(kr*rr/r)**2)**(-1.5)
       rr=(i+(3+SQRT(3.))/6)*dx*r
       evrp(i)=0._dp
       evfp(i)=omega*rr-(kr**2/LOG(1+kr**2))*omega*rr/(1+(kr*rr/r)**2)
       evzp(i)=-omega*r*(kr/(LOG(1+kr**2)*(1+(kr*rr/r)**2))-1/kr)
       elrp(i)=0._dp
       elfp(i)=(kr*rr/r)/SQRT(1+(kr*rr/r)**2)
       elzp(i)=1._dp/SQRT(1+(kr*rr/r)**2)
       ewp(i)=2*omega*(kr**2/LOG(1+kr**2))*(1+(kr*rr/r)**2)**(-1.5)
    END DO
  END SUBROUTINE twistedstate




END MODULE velocities
