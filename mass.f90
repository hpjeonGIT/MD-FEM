!
! ###################### Mass of node calculation ###########################
!
! This subroutine calculates mass of each node and store Jacobian of each
! element for Gauss quadrature.
! 
SUBROUTINE MASS(NS, n, e, bulk)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(MATL):: bulk(NS%Nmat)
!
INTEGER:: i, j, k
REAL*8 :: xi(Nint), eta(Nint), x(Nmesh), y(Nmesh), edge(Nmesh, Ndim), weight, &
     dxdxi(4,Nint)
!
!
weight = 1./SQRT(3.)
xi = RESHAPE((/-weight, weight, weight, -weight/), (/4/))
eta = RESHAPE((/-weight, -weight, weight, weight/), (/4/))
edge = RESHAPE((/-1., 1., 1., -1.,  -1.,-1.,1.,1. /), (/4,2/))
weight = 1.0
!
! Initialize Jacobian of all elements
DO i=1, NS%Nnode
   n(i)%mass = 0.0
END DO
!
!
DO i=1, NS%Nelem
   DO j=1, Nmesh
      x(j) = n(e(i)%l(j))%xx(1)
      y(j) = n(e(i)%l(j))%xx(2)
   END DO
   DO j=1, Nint
      !
      ! dxdxi(1,#) = dx/dxi
      ! dxdxi(2,#) = dy/dxi
      ! dxdxi(3,#) = dx/deta
      ! dxdxi(4,#) = dy/deta
      dxdxi(1,j) = ((x(2)-x(1))*(1.-eta(j)) + (x(3)-x(4))*(1.+eta(j)))/4.
      dxdxi(2,j) = ((y(2)-y(1))*(1.-eta(j)) + (y(3)-y(4))*(1.+eta(j)))/4.
      dxdxi(3,j) = ((x(4)-x(1))*(1.-xi(j)) + (x(3)-x(2))*(1.+xi(j)))/4.
      dxdxi(4,j) = ((y(4)-y(1))*(1.-xi(j)) + (y(3)-y(2))*(1.+xi(j)))/4.
      e(i)%J(j) = dxdxi(1,j)*dxdxi(4,j) - dxdxi(2,j)*dxdxi(3,j)
      IF (e(i)%J(j) < 0.0) STOP "======== Negative Jacobian ======"
   END DO
   !
   DO k=1, Nmesh
      !
      ! Unit thickness(t=1.0)
      DO j=1, Nint
         n(e(i)%l(j))%mass = n(e(i)%l(j))%mass + bulk(1)%rho*0.25*bulk(1)%t* &
              weight*(1.+edge(k,1)*xi(j))*(1.+edge(k,2)*eta(j))*e(i)%J(j)
         ! 
         ! e(i)%BjI(x, SHAPEFTN, GaussPNT)
         e(i)%BjI(1,k,j) = edge(k,1)*(1.+edge(k,2)*eta(j))*dxdxi(4,j) &
              - edge(k,2)*(1.+edge(k,1)*xi(j))*dxdxi(2,j)
         e(i)%BjI(2,k,j) = - edge(k,1)*(1.+edge(k,2)*eta(j))*dxdxi(3,j) &
              + edge(k,2)*(1.+edge(k,1)*xi(j))*dxdxi(1,j)
!         print *, i,j,k,j, e(i)%BjI(1,k,j), e(i)%BjI(2,k,j)
      END DO
   END DO
END DO
!
!
RETURN
END SUBROUTINE MASS
!
! ###################### Internal force calculation ###########################
!
! Calculates stress at each Gauss point of 112.72627339521146all elements.
! 
SUBROUTINE INTFORCE(NS, n, e, bulk, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(MATL):: bulk(NS%Nmat)
REAL*8    :: dt
!
INTEGER:: i, j, k, l
REAL*8 :: weight
!
weight = 1.0
!
!CALL ELAST(NS, n, e, bulk)
CALL ELAST_JAUMANN(NS, n, e, bulk, dt) 
!
DO i=1, NS%Nnode
   DO j=1, Ndim
      n(i)%intF(j) = 0.0
   END DO
END DO
!
DO i=1, NS%Nelem
   ! 
   ! In plane strain mode, 
   ! sigma(#,1) = sigma_xx
   ! sigma(#,2) = sigma_yy
   ! sigma(#,3) = sigma_xy
   DO k=1, Nmesh
      l = e(i)%l(k)
      DO j=1, Nint
         n(l)%intF(1) = n(l)%intF(1) + bulk(1)%t*weight*( &
              e(i)%sigma(j,1)*e(i)%BjI(1,k,j) + &
              e(i)%sigma(j,3)*e(i)%BjI(2,k,j))
         n(l)%intF(2) = n(l)%intF(2) + bulk(1)%t*weight*( &
              e(i)%sigma(j,2)*e(i)%BjI(2,k,j) + &
              e(i)%sigma(j,3)*e(i)%BjI(1,k,j))
      END DO
   END DO
END DO
!
!
RETURN
END SUBROUTINE INTFORCE
