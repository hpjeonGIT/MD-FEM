!
! ################ Pure elastic constitutive routine  ####d###################
!
SUBROUTINE ELAST(NS, n, e, bulk)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(MATL):: bulk(NS%Nmat)
!
INTEGER:: i, j, k, l
REAL*8 :: weight, Eyoung, nu, const
!
weight = 1.0
Eyoung = bulk(1)%C(1)
nu = bulk(1)%C(2)
!
DO i=1, NS%Nelem
   DO j=1, Nint
      e(i)%eps(j,1) = 0.0
      e(i)%eps(j,2) = 0.0
      e(i)%eps(j,3) = 0.0
      DO k=1, Nmesh
         l = e(i)%l(k)
         e(i)%eps(j,1) = e(i)%eps(j,1)+(n(l)%xx(1)-n(l)%ix(1))* &
              e(i)%BjI(1,k,j)/e(i)%J(j)
         e(i)%eps(j,2) = e(i)%eps(j,2)+(n(l)%xx(2)-n(l)%ix(2))* &
              e(i)%BjI(2,k,j)/e(i)%J(j)
         e(i)%eps(j,3) = e(i)%eps(j,3)+ &
              (n(l)%xx(2)-n(l)%ix(2))*e(i)%BjI(1,k,j)*0.5/e(i)%J(j) + &
              (n(l)%xx(1)-n(l)%ix(1))*e(i)%BjI(2,k,j)*0.5/e(i)%J(j)
!         print *, j, k, e(i)%BjI(1,k,j)/e(i)%J(j),e(i)%BjI(2,k,j)/e(i)%J(j)
      END DO
   END DO
END DO         
!
const = Eyoung/(1.+ nu)/(1.-2.*nu)
DO i=1, NS%Nelem
   DO j=1, Nint
      e(i)%sigma(j,1) = const*((1.-nu)*e(i)%eps(j,1) + nu*e(i)%eps(j,2))
      e(i)%sigma(j,2) = const*(nu*e(i)%eps(j,1) + (1.-nu)*e(i)%eps(j,2))
      e(i)%sigma(j,3) = Eyoung*e(i)%eps(j,3)/(1.+nu)
   END DO
END DO
!
RETURN
END SUBROUTINE ELAST
!
! ############ Green strain elastic constitutive routine  ####d###############
! 
! Strain - Green strain
! Stress - 2nd PK stress
! Total lagrangian
! non rate dependent constitutive rule
! For internal force routine, 1st PK is needed.
! 
SUBROUTINE ELAST_GREEN(NS, n, e, bulk)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(MATL):: bulk(NS%Nmat)
!
INTEGER:: i, j, k, l
REAL*8 :: weight, Eyoung, nu, const, F(Ndim,Ndim)
!
weight = 1.0
Eyoung = bulk(1)%C(1)
nu = bulk(1)%C(2)
!
DO i=1, NS%Nelem
   DO j=1, Nint
      !
      ! Initialization
      DO k=1,2
         DO l=1,2
            e(i)%F(j,k,l) = 0.0
         END DO
      END DO
      !
      ! Deformation gradient
      DO k=1, Nmesh
         l = e(i)%l(k)
         e(i)%F(j,1,1) = e(i)%F(j,1,1) + &
              (n(l)%xx(1)-n(l)%ix(1))*e(i)%BjI(1,k,j)/e(i)%J(j)
         e(i)%F(j,1,2) = e(i)%F(j,1,2) + &
              (n(l)%xx(1)-n(l)%ix(1))*e(i)%BjI(2,k,j)/e(i)%J(j)
         e(i)%F(j,2,1) = e(i)%F(j,2,1) + &
              (n(l)%xx(2)-n(l)%ix(2))*e(i)%BjI(1,k,j)/e(i)%J(j)
         e(i)%F(j,2,2) = e(i)%F(j,2,2) + &
              (n(l)%xx(2)-n(l)%ix(2))*e(i)%BjI(2,k,j)/e(i)%J(j)
      END DO
   END DO
   !
   ! Green strain
   DO j=1, Nint      
      e%E(j,1) = 0.5*(e(i)%F(1,1)*e(i)%F(1,1) + e(i)%F(1,2)*e(i)%F(2,1) - 1.0)
      e%E(j,2) = 0.5*(e(i)%F(2,1)*e(i)%F(1,2) + e(i)%F(2,2)*e(i)%F(2,2) - 1.0)
      e%E(j,3) = 0.5*(e(i)%F(1,1)*e(i)%F(1,2) + e(i)%F(2,1)*e(i)%F(1,1))
   END DO
END DO         
!
const = Eyoung/(1.+ nu)/(1.-2.*nu)
DO i=1, NS%Nelem
   DO j=1, Nint
      !
      ! 2nd PK stress
      e(i)%S(j,1) = const*((1.-nu)*e(i)%E(j,1) + nu*e(i)%E(j,2))
      e(i)%S(j,2) = const*(nu*e(i)%E(j,1) + (1.-nu)*e(i)%E(j,2))
      e(i)%S(j,3) = Eyoung*e(i)%E(j,3)/(1.+nu)
   END DO
END DO
!
RETURN
END SUBROUTINE ELAST_GREEN
!
! ###### Rate of deformation for  elastic constitutive routine  ####d##########
! 
! Strain rate - Rate of deformation
! Stress - Cauchy stress
! Updated lagrangian
! 
SUBROUTINE ELAST_JAUMANN(NS, n, e, bulk, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(MATL):: bulk(NS%Nmat)
REAL*8    :: dt
!
INTEGER:: i, j, k, l
REAL*8 :: Eyoung, nu, const, Lvg(2,2), D(3), W(3), cauchy(3), dcauchy(3)
!
Eyoung = bulk(1)%C(1)
nu = bulk(1)%C(2)
const = Eyoung/(1.+ nu)/(1.-2.*nu)
!
DO i=1, NS%Nelem
   DO j=1, Nint
      !
      ! Initialization
      DO k=1,2
         DO l=1,2
            Lvg(k,l) = 0.0
         END DO
      END DO
      !
      ! Deformation gradient
      DO k=1, Nmesh
         l = e(i)%l(k)
         !
         ! L:velocity gradient
         Lvg(1,1) = Lvg(1,1) + n(l)%xv(1)*e(i)%BjI(1,k,j)/e(i)%J(j)
         Lvg(1,2) = Lvg(1,2) + n(l)%xv(1)*e(i)%BjI(2,k,j)/e(i)%J(j) 
         Lvg(2,1) = Lvg(2,1) + n(l)%xv(2)*e(i)%BjI(1,k,j)/e(i)%J(j)
         Lvg(2,2) = Lvg(2,2) + n(l)%xv(2)*e(i)%BjI(2,k,j)/e(i)%J(j)
      END DO
      !
      ! rate of deformation(symmetric)
      D(1) = Lvg(1,1)
      D(2) = Lvg(2,2)
      D(3) = 0.5*(Lvg(1,2) + Lvg(2,1))
      !
      ! Spin tensor(anti-symmetric)
      W(1) = 0.0
      W(2) = 0.0
      W(3) = 0.5*(Lvg(1,2) - Lvg(2,1))
      !
      ! cauchy stress rate
      dcauchy(1) = const*((1.-nu)*D(1) + nu*D(2))
      dcauchy(2) = const*(nu*D(1) + (1.-nu)*D(2))
      dcauchy(3) = Eyoung*D(3)/(1.+nu)
      !
      ! old cauchy stress
      cauchy(1) = e(i)%sigma(j,1)
      cauchy(2) = e(i)%sigma(j,2)
      cauchy(3) = e(i)%sigma(j,3)
      !
      ! Jaumann rate
      e(i)%sigma(j,1) = dcauchy(1)*dt + cauchy(1)
      e(i)%sigma(j,2) = dcauchy(2)*dt + cauchy(2)
      e(i)%sigma(j,3) = (dcauchy(3) + W(3)*cauchy(1) - cauchy(2)*W(3))*dt + &
           cauchy(3)
   END DO
END DO         
!
RETURN
END SUBROUTINE ELAST_JAUMANN
