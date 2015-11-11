!
! ################# FATE ######################################################
!
! FEM and MD hybrid simulation program.
!
! For macro and meso scale analysis, this program is developed using plane
! strain with explicit time integration. Integrating stress wave along given
! element sets, stress and strain will be calculated using conventional
! solid mechanics constitutive equations. 
! For MD part, conventional scheme with two dimensional DOF wil be implemented.
! For coupling MD and FEM, nodal point and particle will be joined.
! Displacement of node will be imposed on particle and force of particle will
! be applied to node as external force term. Then whole process will be
! straightforward.
!
! Byoungseon Jeon
! University of California, Department of Applied Science, graduate student
! Los Alamos National Laboratory, T-12, GRA
!
! 12.28.2005 
! As postporcessor, HAYATE was prepared.
!
! 12.30.2005 
! As preprocessor, NANOHA wass prepared.
!
! 01.05.2006 
! Programming began. 
! 
! 01.12.2006
! Modify FEM-MD coupling routine and NANOHA
!
PROGRAM FATE
USE DATASET
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NODE), POINTER:: n(:)
TYPE(ELEM), POINTER:: e(:)
TYPE(PTCL), POINTER:: q(:)
TYPE(MATL), POINTER:: bulk(:)
TYPE(COND), POINTER:: FBC(:), XBC(:), IC(:)
INTEGER,    POINTER:: NPpair(:,:)
TYPE(PRMT):: param
TYPE(SSTM):: sys
TYPE(NMBR):: NS
TYPE(DATE):: time
REAL*8    :: dt
!
! INTERNAL VARIABLES
REAL:: time0, time1, time2, secnds
time1 = secnds(time0)
!
CALL PARSOR(NS, n, e, q, bulk, FBC, XBC, IC, NPpair, param, time)
!
time%now = 0.0
time%Nloop = 0
time%Ndump = 0
!
! Hybrid analysis
IF ( NS%Npt /= 0 .AND. NS%Nnode /= 0) THEN
   !
   ! Initialize velocity and mass
   CALL INIT(NS, n, e, q, bulk, IC, param, dt)
   !
   ! Open file for time data print
   OPEN(UNIT=15, file="timeplot.dat")
   DO WHILE (time%now < time%max+dt*0.5)
      !
      ! Velocity prediction
      CALL VERLET1(NS, n, e, q, NPpair, param, dt)
      !
      ! Mass of each node and shape function derivative
      CALL MASS(NS, n, e, bulk)
      !
      ! Force for MD part
      CALL FORCE(NS, q, param, sys, dt)
      !
      ! Acceleration for FEM part
      IF (NS%fbc /= 'OFF') CALL BCFORCE(NS, n, e, FBC)   
      CALL ACCEL(NS, n, e, q, NPpair, bulk, dt)
      ! 
      ! Boundary condition
      IF (NS%xbc /= 'OFF') CALL BCFIXED(NS, n, e, XBC)
      ! 
      ! Velocity correction(update0
      CALL VERLET2(NS, n, q, NPpair, sys, param, dt)
      !
      ! Postprocessing intermediate results
      IF ( INT((time%now+.5*dt)/time%dump) > time%Ndump-1) THEN      
         CALL POST(NS, n, e, q, NPpair, time, param, dt)
         PRINT 100,  dt, time%now, time%max
      END IF
      CALL TPOST(NS, n, e, q, time)
      time%now = time%now + dt
      time%Nloop = time%Nloop + 1
   END DO
   !
   DEALLOCATE(n, e, q, bulk, IC, NPpair)
   CLOSE(15)
   IF (NS%fbc /= 'OFF') DEALLOCATE(FBC)
   IF (NS%xbc /= 'OFF') DEALLOCATE(XBC)
!
! Only FEM analysis
ELSE IF ( NS%Npt == 0 .AND. NS%Nnode /= 0) THEN
   CALL INIT_FEM(NS, n, e, bulk, IC, dt)
   OPEN(UNIT=15, file="timeplot.dat")
   DO WHILE (time%now < time%max+dt*0.5)
      CALL VERLET_FEM1(NS, n, e, param, dt)
      CALL MASS(NS, n, e, bulk)
      IF (NS%fbc /= 'OFF') CALL BCFORCE(NS, n, e, FBC)   
      CALL ACCEL_FEM(NS, n, e, bulk, dt)
      IF (NS%xbc /= 'OFF') CALL BCFIXED(NS, n, e, XBC)
      CALL VERLET_FEM2(NS, n, sys, param, dt)
      IF ( INT((time%now+.5*dt)/time%dump) > time%Ndump-1) THEN      
         CALL POST_FEM(NS, n, e, time, param, dt)
         PRINT 100,  dt, time%now, time%max
      END IF
      CALL TPOST_FEM(NS, n, e, time)
      time%now = time%now + dt
      time%Nloop = time%Nloop + 1
   END DO
   DEALLOCATE(n, e, bulk, IC)
   CLOSE(15)
   IF (NS%fbc /= 'OFF') DEALLOCATE(FBC)
   IF (NS%xbc /= 'OFF') DEALLOCATE(XBC)
!
! Only MD analysis
ELSE IF ( NS%Npt /= 0 .AND. NS%Nnode == 0) THEN
   CALL INIT_MD(NS, q, dt)
   OPEN(UNIT=15, file="timeplot.dat")
   DO WHILE (time%now < time%max+dt*0.5)
      CALL VERLET_MD1(NS, q, param, dt)
      CALL FORCE(NS, q, param, sys, dt)
      CALL VERLET_MD2(NS, q, sys, param, dt)
      IF ( INT((time%now+.5*dt)/time%dump) > time%Ndump-1) THEN      
         CALL POST_MD(NS, q, time, param, dt)
         PRINT 100,  dt, time%now, time%max
      END IF
      CALL TPOST_MD(NS, q, time)
      time%now = time%now + dt
      time%Nloop = time%Nloop + 1
   END DO
   DEALLOCATE(q, bulk, IC)
   CLOSE(15)
ELSE
   PRINT *, "No problem set given"
END IF
!
time2 = secnds(time1)
PRINT *, "Wall time is", time2
100 FORMAT("dt = ", ES10.4, " at ", ES8.2, " out of ", ES8.2)
!
CONTAINS
!
! ###################### INPUT DATA parsing routine ###########################
!
! This subroutine reads positions of all nodes and component node sets for each
! element. Besides, initial/boundary condition and material properties will
! be available through this routine.
! 
! Unit normalization
! Length: 1. means 1Angstrom = 10E-10 m
! Mass: 1. means 1.6605E-27 kg
! Energy: 1. means 1 eV = 1.602E-19 J
! Time: 1. means 10.2192 fs = 1.02192E-14 sec
! 1 Pa = 1kg/m/s/s = 6.289193E-12
! 200GPa = 1.2578386
!
SUBROUTINE PARSOR(NS, n, e, q, bulk, FBC, XBC, IC, NPpair, param, time)
USE DATASET
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NODE), POINTER:: n(:)
TYPE(ELEM), POINTER:: e(:)
TYPE(PTCL), POINTER:: q(:)
TYPE(MATL), POINTER:: bulk(:)
TYPE(COND), POINTER:: FBC(:), XBC(:), IC(:)
INTEGER,    POINTER:: NPpair(:,:)
TYPE(PRMT):: param
TYPE(SSTM):: sys
TYPE(NMBR):: NS
TYPE(DATE):: time
!
! INTERNAL VARIABLES
INTEGER:: AllocateStatus, i, j
CHARACTER*8:: dummy
!
![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[ Geometric data ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
!
OPEN(UNIT=11, file="model.inp")
!
! HEADER
READ(11,*) dummy
READ(11,*) dummy ! 
READ(11,*) time%max, time%dump, dt, param%temp
!
! Node data
READ(11,*) dummy
READ(11,*) NS%Nnode
IF (NS%Nnode /= 0) THEN
   ALLOCATE(n(NS%Nnode), stat=AllocateStatus)
   IF (AllocateStatus /= 0) STOP " == Not Enough Memory for Node Allocation =="
   DO i=1, NS%Nnode
      READ(11,*) dummy, (n(i)%ix(j), j=1, Ndim)
   END DO
END IF
!
! Element data
READ(11,*) dummy 
READ(11,*) NS%Nelem
IF (NS%Nelem /= 0) THEN
   ALLOCATE(e(NS%Nelem), stat=AllocateStatus)
   IF (AllocateStatus /= 0) STOP " == Not Enough Memory for El. Allocation =="
   DO i=1, NS%Nelem
      READ(11,*) dummy, (e(i)%l(j), j=1, 4)
   END DO
END IF
!
! Particle data
READ(11,*) dummy
READ(11,*) NS%Npt
IF (NS%Npt /=0) THEN
   ALLOCATE(q(NS%Npt), stat=AllocateStatus)
   IF (AllocateStatus /=0) STOP " == Not Enough Memory for Pt. Allocation =="
   DO i=1, NS%Npt
      !
      ! Read particle position
      READ(11,*) dummy, (q(i)%xx(j), j=1, Ndim), q(i)%id
   END DO
END IF
!
! Node-particle pair
READ(11,*) dummy
READ(11,*) NS%Npair
IF (NS%Npair /= 0) THEN
   ALLOCATE(NPpair(NS%Npair, 5), stat=AllocateStatus)
   !
   ! NPpair(*,1) = element
   ! NPpair(*,2+) = particle
   IF (AllocateStatus /= 0) STOP " == Not Enough Memory for Node Allocation =="
   DO i=1, NS%Npair
      !
      ! 
      READ(11,*) NPpair(i,1), NPpair(i,2)
   END DO
END IF
!
! Force boundary condition
READ(11,*) dummy 
READ(11,*) NS%Nfbc
IF (NS%Nfbc /=0) THEN
   NS%fbc = 'ON '
   ALLOCATE(FBC(NS%Nfbc))
   DO i=1, NS%Nfbc
      !
      ! READ node and force
      READ(11,*) FBC(i)%node, (FBC(i)%pp(j), j=1, Ndim)
   END DO
ELSE
   NS%fbc = 'OFF'
END IF
!
! Fixed boundary condition
READ(11,*) dummy 
READ(11,*) NS%Nxbc
IF (NS%Nxbc /= 0) THEN
   NS%xbc = 'ON '
   ALLOCATE(XBC(NS%Nxbc))
   DO i=1, NS%Nxbc
      !
      ! READ node and displacement
      READ(11,*) XBC(i)%node, (XBC(i)%pp(j), j=1, Ndim)
   END DO
ELSE
   NS%xbc = 'OFF'
END IF
!
! Initial condition
READ(11,*) dummy
READ(11,*) NS%Nic
ALLOCATE(IC(NS%Nic))
DO i=1, NS%Nic
   !
   ! READ node and velocity
   READ(11,*) IC(i)%node, (IC(i)%pp(j), j=1, Ndim)
END DO
!
CLOSE(11)
!
![[[[[[[[[[[[[[[[[[[  Material property and potential data ]]]]]]]]]]]]]]]]]]]]
!
OPEN(UNIT=12, file="potential.prm")
!
READ(12,*) dummy ! Material data
READ(12,*) NS%Nmat
ALLOCATE(bulk(NS%Nmat))
DO i=1, NS%Nmat
   !
   !          density      thickness     elastic     poisson ratio
   READ(12,*) bulk(i)%rho, bulk(i)%t, bulk(i)%C(1), bulk(i)%C(2)
END DO
!
! Cohesive energy
READ(12,*) dummy
READ(12,*) param%Ecoh
!
! NN distance for lattice
READ(12,*) dummy
READ(12,*) param%r0
!
! LJ/EAM ratio
READ(12,*) dummy
READ(12,*) param%chi
!
! Mass for particle kind
READ(12,*) dummy
READ(12,*) param%xm(Nparam)
!
! Thermostat
READ(12,*) dummy
READ(12,*) param%thermo
READ(12,*) param%tau, param%alpha
!
CLOSE(12)
!
END SUBROUTINE PARSOR
!
END PROGRAM FATE
