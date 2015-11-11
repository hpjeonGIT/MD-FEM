!
! ######################## Initialization routine ############################
!
! Initialization before loop.
!
SUBROUTINE INIT(NS, n, e, q, bulk, IC, param, dt)
USE DATASET
IMPLICIT NONE
INTERFACE
   FUNCTION fluct(x)
     IMPLICIT NONE
     REAL*8:: fluct, x, l0, l, r, v1, v2
     REAL*8:: rand1, rand2, ran2
   END FUNCTION fluct
END INTERFACE
!
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(PTCL):: q(NS%Npt)
TYPE(MATL):: bulk(NS%Nmat)
TYPE(COND):: IC(NS%Nic)
TYPE(PRMT):: param
REAL*8    :: dt
!
INTEGER:: i, j, k, l1, l2, Nr
REAL*8 :: x, y, Le, Cd, xv(NS%Npt, Ndim), xm, mv2, kb, T0, lambda
!
!
kb = 8.617343E-5
IF (param%thermo == 'STOCH ') THEN
   Nr = 0
ELSE
   Nr = Nref
END IF
!
! Node initialization
DO i=1, NS%Nnode
   DO j=1, Ndim
      n(i)%xx(j) = n(i)%ix(j)
      n(i)%xv(j) = 0.0
      n(i)%xa(j) = 0.0
      n(i)%extF(j) = 0.0
   END DO
END DO
!
! Particle initialization
IF (param%temp > 1.E-6) THEN
   mv2 = 0.0
   !
   ! Distribute velocity using normal distribution with
   ! standard deviation of the given temperature
   DO i=1, NS%Npt
      xm = param%xm(q(i)%id)
      q(i)%rho = 0.0
      DO j=1, Ndim
         xv(i,j) = fluct(SQRT(param%temp))
         mv2 = mv2 + xm*xv(i,j)**2
         q(i)%drho(j) = 0.0
      END DO
   END DO
   !
   ! Rescale all velocity using isokinetic scaling
   T0 = mv2/kb/REAL(3*NS%Npt - Nr)
   !
   lambda = SQRT(param%temp/T0)
   DO i = 1, NS%Npt
      xm = param%xm(q(i)%id)
      DO j = 1,Ndim
         q(i)%xv(j) = xv(i,j)*lambda
      END DO
   END DO
ELSE
   DO i = 1, NS%Npt
      DO j = 1,Ndim
         q(i)%xv(j) = 0.0
      END DO   
   END DO
END IF
!
DO i=1, NS%Nic
   k = IC(i)%node
   DO j=1, Ndim
      n(k)%xv(j) = IC(i)%pp(j)
   END DO
END DO
!
RETURN
END SUBROUTINE INIT
!
! ############# Initialization routine for FEM only simulation ################
!
! Initialization before loop.
!
SUBROUTINE INIT_FEM(NS, n, e, bulk, IC, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(MATL):: bulk(NS%Nmat)
TYPE(COND):: IC(NS%Nic)
REAL*8:: dt
!
INTEGER:: i, j, k, l1, l2
!
! Node initialization
DO i=1, NS%Nnode
   DO j=1, Ndim
      n(i)%xx(j) = n(i)%ix(j)
      n(i)%xv(j) = 0.0
      n(i)%xa(j) = 0.0
      n(i)%extF(j) = 0.0
   END DO
END DO
!
! Element initialization
DO i = 1, NS%Nelem
   DO j=1, Nint
      e(i)%sigma(j,1) = 0.
      e(i)%sigma(j,2) = 0.
      e(i)%sigma(j,3) = 0.
   END DO
END DO
!
DO i=1, NS%Nic
   k = IC(i)%node
   DO j=1, Ndim
      n(k)%xv(j) = IC(i)%pp(j)
   END DO
END DO
!
RETURN
END SUBROUTINE INIT_FEM
!
! ############# Initialization routine for MD only simulation  ###############
!
! Initialization before loop.
!
SUBROUTINE INIT_MD(NS, q, param, dt)
USE DATASET
IMPLICIT NONE
INTERFACE
   FUNCTION fluct(x)
     IMPLICIT NONE
     REAL*8:: fluct, x, l0, l, r, v1, v2
     REAL*8:: rand1, rand2, ran2
   END FUNCTION fluct
END INTERFACE
!
TYPE(NMBR):: NS
TYPE(PTCL):: q(NS%Npt)
TYPE(PRMT):: param
REAL*8    :: dt
!
INTEGER:: i, j, k, l1, l2, Nr
REAL*8 :: x, y, Le, Cd, xv(NS%Npt, Ndim), xm, mv2, kb, T0, lambda
!
!
kb = 8.617343E-5
IF (param%thermo == 'STOCH ') THEN
   Nr = 0
ELSE
   Nr = Nref
END IF
!
! Particle initialization
IF ( param%temp > 1.E-6 ) THEN
   mv2 = 0.0
   !
   ! Distribute velocity using normal distribution with
   ! standard deviation of the given temperature
   DO i=1, NS%Npt
      xm = param%xm(q(i)%id)
      q(i)%rho = 0.0
      DO j=1, Ndim
         xv(i,j) = fluct(SQRT(param%temp))
         mv2 = mv2 + xm*xv(i,j)**2
         q(i)%drho(j) = 0.0
      END DO
   END DO
   !
   ! Rescale all velocity using isokinetic scaling
   T0 = mv2/kb/REAL(3*NS%Npt - Nr)
   !
   lambda = SQRT(param%temp/T0)
   DO i = 1, NS%Npt
      xm = param%xm(q(i)%id)
      DO j = 1,Ndim
         q(i)%xv(j) = xv(i,j)*lambda
      END DO
   END DO
ELSE
   DO i = 1, NS%Npt
      DO j = 1,Ndim
         q(i)%xv(j) = 0.0
      END DO
   END DO
END IF
!
RETURN
END SUBROUTINE INIT_MD
!
! ###################### Velocity prediction routine #########################
!
! Central difference time integration routine 1 - equivalent to velocity Verlet
! scheme in MD.
!
SUBROUTINE VERLET1(NS, n, e, q, NPpair, param, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(PTCL):: q(NS%Npt)
INTEGER   :: NPpair(NS%Npair,5)
TYPE(PRMT):: param
REAL*8:: dt
!
INTEGER:: i, j, el, p1, p2, p3, p4, n1, n2, n3, n4
REAL*8 :: xm
!
! Node update
DO i=1, NS%Nnode
   DO j=1, Ndim
      n(i)%xx(j) = n(i)%xx(j) + dt*n(i)%xv(j) + 0.5*n(i)%xa(j)*dt*dt
      n(i)%xv(j) = n(i)%xv(j) + 0.5*dt*n(i)%xa(j)
      n(i)%extF(j) = 0.0
   END DO   
END DO
!
! Particle update
DO i=1, NS%Npt
   xm = param%xm(q(i)%id)
   DO j=1, Ndim
      q(i)%xv(j) = q(i)%xv(j) + 0.5*dt*q(i)%ff(j)/xm
      q(i)%xx(j) = q(i)%xx(j) + dt*q(i)%xv(j)
   END DO
END DO
!
! Pair particle update
DO i=1, NS%Npair
   el = NPpair(i,1)
   p1 = NPpair(i,2)
   n1 = e(el)%l(1)
   n2 = e(el)%l(2)
   n3 = e(el)%l(3)
   n4 = e(el)%l(4)
   DO j=1, Ndim
      q(p1)%xx(j) = (n(n1)%xx(j) + n(n2)%xx(j) + n(n4)%xx(j) + n(n3)%xx(j))/4.
   END DO
END DO
!
RETURN
END SUBROUTINE VERLET1
!
SUBROUTINE VERLET_FEM1(NS, n, e, param, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(PRMT):: param
REAL*8:: dt
!
INTEGER:: i, j
!
! Node update
DO i=1, NS%Nnode
   DO j=1, Ndim
      n(i)%xx(j) = n(i)%xx(j) + dt*n(i)%xv(j) + 0.5*n(i)%xa(j)*dt*dt
      n(i)%xv(j) = n(i)%xv(j) + 0.5*dt*n(i)%xa(j)
      n(i)%extF(j) = 0.0
   END DO   
END DO
!
RETURN
END SUBROUTINE VERLET_FEM1
!
SUBROUTINE VERLET_MD1(NS, q, param, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(PTCL):: q(NS%Npt)
TYPE(PRMT):: param
REAL*8:: dt
!
INTEGER:: i, j
REAL*8 :: xm
!
! Particle update
DO i=1, NS%Npt
   xm = param%xm(q(i)%id)
   DO j=1, Ndim
      q(i)%xv(j) = q(i)%xv(j) + 0.5*dt*q(i)%ff(j)/xm
      q(i)%xx(j) = q(i)%xx(j) + dt*q(i)%xv(j)
   END DO
END DO
!
RETURN
END SUBROUTINE VERLET_MD1
!
! ###################### Velocity correction routine #########################
!
! Central difference time integration routine 2 - 
SUBROUTINE VERLET2(NS, n, q, NPpair, sys, param, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(PTCL):: q(NS%Npt)
INTEGER   :: NPpair(NS%Npair,5)
TYPE(SSTM):: sys
TYPE(PRMT):: param
REAL*8:: dt
!
INTEGER:: i, j
REAL*8 :: xm
!
! Node update
DO i=1, NS%Nnode
   DO j=1, Ndim
      n(i)%xv(j) = n(i)%xv(j) + 0.5*dt*n(i)%xa(j)
   END DO   
END DO
!
! Particle update
sys%mv2 = 0.0
DO i=1, NS%Npt
   xm = param%xm(q(i)%id)
   DO j=1, Ndim
      q(i)%xv(j) = q(i)%xv(j) + 0.5*dt*q(i)%ff(j)/xm
      sys%mv2 = sys%mv2 + xm*q(i)%xv(j)**2
   END DO
END DO
!
RETURN
END SUBROUTINE VERLET2
!
SUBROUTINE VERLET_FEM2(NS, n, sys, param, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(SSTM):: sys
TYPE(PRMT):: param
REAL*8:: dt
!
INTEGER:: i, j
REAL*8 :: xm
!
! Node update
DO i=1, NS%Nnode
   DO j=1, Ndim
      n(i)%xv(j) = n(i)%xv(j) + 0.5*dt*n(i)%xa(j)
   END DO   
END DO
!
RETURN
END SUBROUTINE VERLET_FEM2
!
SUBROUTINE VERLET_MD2(NS, q, sys, param, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(PTCL):: q(NS%Npt)
TYPE(SSTM):: sys
TYPE(PRMT):: param
REAL*8:: dt
!
INTEGER:: i, j
REAL*8 :: xm
!
! Particle update
sys%mv2 = 0.0
DO i=1, NS%Npt
   xm = param%xm(q(i)%id)
   DO j=1, Ndim
      q(i)%xv(j) = q(i)%xv(j) + 0.5*dt*q(i)%ff(j)/xm
      sys%mv2 = sys%mv2 + xm*q(i)%xv(j)**2
   END DO
END DO
!
RETURN
END SUBROUTINE VERLET_MD2
!
! ##################### Acceleration estimation routine #######################
!
! 
SUBROUTINE ACCEL(NS, n, e, q, NPpair, bulk, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(PTCL):: q(NS%Npt)
INTEGER   :: NPpair(NS%Npair,5)
TYPE(MATL):: bulk(NS%Nmat)
REAL*8:: dt
!
INTEGER:: i, j, el, p1, p2, p3, p4, n1, n2, n3, n4
!
CALL INTFORCE(NS, n, e, bulk, dt)
!goto 10
!
! Pair node update
DO i=1, NS%Npair
   el = NPpair(i,1)
   p1 = NPpair(i,2)
   n1 = e(el)%l(1)
   n2 = e(el)%l(2)
   n3 = e(el)%l(3)
   n4 = e(el)%l(4)
   DO j=1, Ndim
      n(n1)%extF(j) = n(n1)%extF(j) + q(p1)%ff(j)/4.
      n(n2)%extF(j) = n(n2)%extF(j) + q(p1)%ff(j)/4.
      n(n3)%extF(j) = n(n3)%extF(j) + q(p1)%ff(j)/4.
      n(n4)%extF(j) = n(n4)%extF(j) + q(p1)%ff(j)/4.
   END DO
END DO
!
! Acceleration of nodes
10 DO i=1, NS%Nnode
   DO j=1, Ndim
      n(i)%xa(j) = (n(i)%extF(j) - n(i)%intF(j))/n(i)%mass
   END DO
END DO
!
RETURN
END SUBROUTINE ACCEL
! 
SUBROUTINE ACCEL_FEM(NS, n, e, bulk, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(MATL):: bulk(NS%Nmat)
REAL*8:: dt
!
INTEGER:: i, j
!
CALL INTFORCE(NS, n, e, bulk, dt)
!
! Acceleration of nodes
10 DO i=1, NS%Nnode
   DO j=1, Ndim
      n(i)%xa(j) = (n(i)%extF(j) - n(i)%intF(j))/n(i)%mass
   END DO
END DO
!
RETURN
END SUBROUTINE ACCEL_FEM
!
! ##################### Force Boundary condition routine ######################
!
SUBROUTINE BCFORCE(NS, n, e, FBC)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(COND):: FBC(NS%Nfbc)
!
INTEGER:: i, j, k
!
!
DO i=1, NS%Nfbc
   k = FBC(i)%node
   DO j=1, Ndim
      n(k)%extF(j) = FBC(i)%pp(j)
   END DO
END DO
!
RETURN
END SUBROUTINE BCFORCE
!
! ##################### Fixed Boundary condition routine ######################
!
SUBROUTINE BCFIXED(NS, n, e, XBC)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(COND):: XBC(NS%Nxbc)
!
INTEGER:: i, j, k
!
!
DO i=1, NS%Nxbc
   k = XBC(i)%node
   DO j=1, Ndim
      IF (XBC(i)%pp(j) < 99998.) THEN         
         n(k)%xx(j) = XBC(i)%pp(j)
      END IF
   END DO
END DO
!
RETURN
END SUBROUTINE BCFIXED
! ############## Random Gaussian(Normal) Distribution Function ################
!
! For stochastic thermostat, fluctuation dissipation theorem is implemented.
! Basically, random number generator which follows Gaussian distribution is
! needed to implement this thermal noise force.
! Random number is generated using FORTRAN intrinsic fucntion - RANDOM_SEED
! and RANDOM_NUMBER. But during the implementation, it is found that those
! intrinsic functions may not work well under false seeding - to use this
! routine on new machine or new compiler, please make sure that the
! distribution follows zero-mean with suitable deviation.
!
! This function provides a random number with zero-mean and deviation x
! along Gaussian distribution.
! <p> = 0.0
! <p**2> = x**2
!
FUNCTION fluct(x)
  IMPLICIT NONE
  REAL*8:: fluct, x, r, v1, v2
  REAL*8:: rand1, rand2, ran2
  REAL:: ranf
  !
  ! Initialization
  r=1.
  DO WHILE (r.ge.1.)
     CALL RANDOM_NUMBER(rand1)
     CALL RANDOM_NUMBER(rand2)
     v1 = 2.*rand1 - 1.
     v2 = 2.*rand2 - 1.
     r = v1*v1+v2*v2
  END DO
  fluct = v1*SQRT(-2.*log(r)/r)*x
  RETURN
END FUNCTION fluct
