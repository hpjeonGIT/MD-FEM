MODULE DATASET
  IMPLICIT NONE
  !
  ! Ndim = dimension of problem
  ! Nint = order of integration
  ! Nmesh = order of mesh - trigonal/quad
  INTEGER, PARAMETER:: Ndim = 2, Nint = 4, Nmesh = 4, Nparam = 1, Nref = 3
  !
  ! NUMBER type
  ! Nnode = Number of all nodes
  ! Nelem = Number of all elements
  ! Nmat  = Number of all material kinds
  ! Nfbc  = Number of force boundary condition sets
  ! Nxbc  = Number of fixed boundary condition sets
  ! Nic   = Number of initial condition sets
  ! Npt   = Number of all particles
  TYPE NMBR
     INTEGER:: Nnode, Nelem, Nmat, Nfbc, Nxbc, Nic
     INTEGER:: Npt, Npair
     CHARACTER*3:: fbc, xbc
  END TYPE NMBR
  !
  ! Node data type
  ! ix = initial position
  ! xx = current position
  ! xv = current velocity
  ! xa = current acceleration
  ! mass = mass per node
  ! intF = internal force on node
  ! extF = external force on node
  TYPE NODE
     REAL*8 :: ix(Ndim), xx(Ndim), xv(Ndim), xa(Ndim), mass, &
          intF(Ndim), extF(Ndim)
  END TYPE NODE
  !
  ! Particle data type
  ! xx = current position
  ! xv = current velocity
  ! ff = current force
  ! rho = density for EAM
  ! drho = gradient of rho
  ! id = id of each particle
  TYPE PTCL
     REAL*8 :: xx(3), xv(3), ff(3), rho, drho(2)
     INTEGER:: id
  END TYPE PTCL
  !
  ! ELEMENT type
  ! l = Node number for component of element set
  ! J = Jacobian for Gauss quadrature mapping
  ! BjI = Shape function spatial derivative(without 1/|J|)
  ! eps = strain tensor 
  ! sigma = Cauchy stress tensor
  ! E = green strain
  ! S = 2nd PK stress 
  ! F = deformation gradient
  TYPE ELEM
     INTEGER:: l(Nint)
     REAL*8 :: J(Nint), BjI(Ndim,Nmesh,Nint), eps(Nint,3), sigma(Nint,3), &
          E(Nint,3), S(Nint,3), F(Nint, Ndim, Ndim)
  END TYPE ELEM
  !
  ! MATERIAL type
  ! rho = density
  ! t = thickness
  ! C(#) = material property depending on constitutive equations
  TYPE MATL
     REAL*8 :: rho, t, C(10)
  END TYPE MATL
  !
  ! PARAMETER type
  ! Ecoh = cohesive energy
  ! r0 = NN distance for lattice
  ! chi = LJ/EAM ratio
  ! xm = mass
  ! temp = given temperature
  TYPE PRMT
     REAL*8 :: Ecoh, r0, chi, xm(Nparam), temp, tau, alpha
     CHARACTER*6 :: thermo
  END TYPE PRMT
  !
  ! Boundary condition type
  TYPE COND
     INTEGER:: node
     REAL*8 :: pp(Ndim)
  END TYPE COND
  !
  ! System parameter
  TYPE SSTM
     REAL*8:: mv2, Epot
  END TYPE SSTM
  !
  ! DATE type
  ! Nloop = Number of iteration loop
  ! Ndump = Number of dump frequency
  ! now = current time
  ! dump = intermediate result dump frequency
  ! max = termination time
  TYPE DATE
     INTEGER:: Nloop, Ndump
     REAL*8 :: now, dump, max
  END TYPE DATE
  !
END MODULE DATASET
