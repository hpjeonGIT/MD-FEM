SUBROUTINE FORCE(NS, q, param, sys, bulk, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(PTCL):: q(NS%Npt)
TYPE(PRMT):: param
TYPE(SSTM):: sys
TYPE(MATL):: bulk(NS%Nmat)
REAL*8:: dt
!
INTEGER:: i, j, k, id_a, id_c
REAL*8 :: r2, xr(Ndim), eps, ecut, sigma, r2i, rc2, r02, r6i, lj, ee, chi, d,&
     rspl2, a2, a3, Ecoh, rmax2, r0, Ncoord, rho0
!
!
Ecoh = param%Ecoh
r0 = param%r0
chi = param%chi
d = REAL(Ndim)
Ncoord = d*(d+1.)
eps = Ecoh*2./Ncoord
sigma = r0/2.**(1./6.)
r02 = r0**2
rspl2 = (sigma*1.244455)**2
rmax2 = (sigma*1.711238)**2
ee = EXP(1.0)
a2 = 0.5424494
a3 = 0.09350527
rho0 = 1./ee !0.3678794411717238
!
DO i=1, NS%Npt
   q(i)%rho = 0.0
   DO j=1, Ndim
      q(i)%ff(j) = 0.0
      q(i)%drho(j) = 0.0
   END DO
END DO
!
! [[[[[[[[[[[[[[[[[[[[[[[[[[[ LJ 6-12 potential ]]]]]]]]]]]]]]]]]]]]]]]]]]]
DO i=1, NS%Npt-1
   id_c = q(i)%id
   DO j = i+1, NS%Npt
      id_a = q(j)%id
      r2 = 0.0
      DO k = 1, Ndim
         xr(k) = q(i)%xx(k) - q(j)%xx(k)
         r2 = r2 + xr(k)**2
      END DO
      !rc2 = param%rc2(id_c, id_a)
      IF (r2 < rmax2) THEN
         !eps = param%eps(id_c, id_a)
         !ecut = param%ecut(id_c, id_a)
         !sigma = param%sig(id_c, id_a)     
         !r02 = sigma**2*2.**(1./3.)
         !rspl2 = (sigma*1.244455)**2
         IF ( r2 < rspl2) THEN
            r2i= 1./r2
            r6i = sigma**6*r2i**3
            lj = 48.*r2i*r6i*(r6i - .5)
            DO k = 1, Ndim
               q(i)%ff(k) = q(i)%ff(k) + lj*xr(k)*eps*chi
               q(j)%ff(k) = q(j)%ff(k) - lj*xr(k)*eps*chi               
            END DO
            sys%Epot = sys%Epot + 4.*r6i*(r6i - 1.)*eps*chi
         ELSE
            lj = (rmax2 - r2)*(-4.*a2 + 6.*a3*(rmax2 - r2)/sigma**2)/sigma**4
            DO k = 1, Ndim
               q(i)%ff(k) = q(i)%ff(k) + lj*xr(k)*eps*chi
               q(j)%ff(k) = q(j)%ff(k) - lj*xr(k)*eps*chi
            END DO
            sys%Epot = sys%Epot + (-a2*(rmax2-r2)**2/sigma**4 &
                 +a3*(rmax2-r2)**3/sigma**6)*eps*chi
         END IF            
         !
         ! Density update for EAM
         q(i)%rho = q(i)%rho + ((rmax2 - r2)/(rmax2 - r02))**2 
         q(j)%rho = q(j)%rho + ((rmax2 - r2)/(rmax2 - r02))**2 
         DO k = 1, Ndim
            q(i)%drho(k) = q(i)%drho(k) + &
                 (rmax2 - r2)*4.*xr(k)/(rmax2 - r02)**2
            q(j)%drho(k) = q(j)%drho(k) - &
                 (rmax2 - r2)*4.*xr(k)/(rmax2 - r02)**2
         END DO                
      END IF
   END DO
END DO
!
! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[ EAM potential ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
DO i=1, NS%Npt
   q(i)%rho = q(i)%rho*rho0/Ncoord
   sys%Epot = sys%Epot + (1.-chi)*Ecoh*q(i)%rho*LOG(q(i)%rho)/rho0
   DO k=1, Ndim
      q(i)%drho(k) = q(i)%drho(k)*rho0/Ncoord
      q(i)%ff(k) = q(i)%ff(k) + (1.-chi)*Ecoh*(LOG(q(i)%rho)+1.)* &
           q(i)%drho(k)/rho0
   END DO
END DO
!
!
!DO i=1,17
!q(i)%ff(2) = q(i)%ff(2)  + 0.1
!END DO
!q(41)%ff(2) = q(41)%ff(1) + 0.1
!q(81)%ff(2) = q(81)%ff(1) + 0.1
!q(122)%ff(1) = q(122)%ff(1) + 0.1
!q(162)%ff(1) = q(162)%ff(1) + 0.1
!q(203)%ff(1) = q(203)%ff(1) + 0.1
!q(243)%ff(1) = q(243)%ff(1) + 0.1
!q(284)%ff(1) = q(284)%ff(1) + 0.1
!q(324)%ff(1) = q(324)%ff(1) + 0.1
!q(365)%ff(1) = q(365)%ff(1) + 0.1
!
100 RETURN
!
END SUBROUTINE FORCE
