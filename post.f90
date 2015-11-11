!
! ######################## POSTprocessing routine #############################
!
SUBROUTINE POST(NS, n, e, q, NPpair, time, param, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(PTCL):: q(NS%Npt)
TYPE(DATE):: time
TYPE(PRMT):: param
INTEGER   :: NPpair(NS%Npair,5)
REAL*8    :: dt
!
INTEGER:: i, j, k
REAL*8 :: kb
CHARACTER(LEN=14):: NUMBER, FILENAME, RSLTNAME
DATA NUMBER/'0123456789'/
FILENAME = "hybrid00.geo"
RSLTNAME = "hybrid00.res"
kb = 8.617343E-5
!
IF(INT(time%Ndump/10) < 1 ) THEN
   FILENAME(7:7) = NUMBER(1:1)
   FILENAME(8:8) = NUMBER(time%Ndump+1:time%Ndump+1)
   RSLTNAME(7:7) = NUMBER(1:1)
   RSLTNAME(8:8) = NUMBER(time%Ndump+1:time%Ndump+1)
ELSE
   FILENAME(7:7) = NUMBER(INT(time%Ndump/10)+1:INT(time%Ndump/10)+1)
   FILENAME(8:8) = NUMBER(time%Ndump-INT(time%Ndump/10)*10+1: &
        time%Ndump-INT(time%Ndump/10)*10+1)
   RSLTNAME(7:7) = NUMBER(INT(time%Ndump/10)+1:INT(time%Ndump/10)+1)
   RSLTNAME(8:8) = NUMBER(time%Ndump-INT(time%Ndump/10)*10+1: &
        time%Ndump-INT(time%Ndump/10)*10+1)
END IF
!
! ####################### MESH GEOMETRY PRINT #################################
!
OPEN (UNIT=20, FILE=FILENAME)
!
! Node output
WRITE(20, 10)
WRITE(20, 20) time%now
WRITE(20, 30) time%max, time%dump, dt
WRITE(20, 40)
WRITE(20, 50) NS%Nnode
DO i=1, NS%Nnode
   WRITE(20, 60) i, n(i)%xx(1), n(i)%xx(2)
END DO
!
! Element output 
WRITE(20, 70)
WRITE(20, 50) NS%Nelem
DO i=1, NS%Nelem   
   WRITE(20, 80)  i, (e(i)%l(j), j=1, 4)
END DO
!
! Particle output
WRITE(20, 90)
WRITE(20, 50) NS%Npt
DO i=1, NS%Npt
   WRITE(20, 60) i, q(i)%xx(1), q(i)%xx(2)
END DO 
!
! Pair
WRITE(20, 90)
WRITE(20, 50) NS%Npair
DO i=1, NS%Npair
   WRITE(20, 80) (NPpair(i,j), j=1,2)
END DO
!
!
10 FORMAT("# hybrid model geometric data")
20 FORMAT("# time max, dum, dt", ES14.4)
30 FORMAT(3(F10.3, 1X))
40 FORMAT("# Node data")
50 FORMAT(I5)
60 FORMAT(I5, 2(ES14.4, 1X))
70 FORMAT("# Element data")
80 FORMAT(5(I5, 1X))
90 FORMAT("# Particle data")
!
CLOSE(20)
!
! ####################### FEA RESULT PRINT #################################
!
OPEN (UNIT=30, FILE=RSLTNAME)
WRITE(30, 110)
WRITE(30, 120)
WRITE(30, 130)
WRITE(30, 140)
WRITE(30, 150)
!
! Stress at Gauss points
WRITE(30, 160)
WRITE(30, 170)
WRITE(30, 180)
DO i=1, NS%Nelem
   WRITE(30,190) i, (e(i)%sigma(1,j), j=1,3)
   WRITE(30,200)    (e(i)%sigma(2,j), j=1,3)
   WRITE(30,200)    (e(i)%sigma(3,j), j=1,3)   
   WRITE(30,200)    (e(i)%sigma(4,j), j=1,3)   
END DO
WRITE(30, 250)
!
! Nodal displacement
!WRITE(30, 210)
!WRITE(30, 220)
!WRITE(30, 230)
!DO i=1, NS%Nnode
!   WRITE(30,240) i,  n(i)%xx(1)-n(i)%ix(1), n(i)%xx(2)-n(i)%ix(2)
!END DO
!WRITE(30, 250)
!
! particle temperature
WRITE(30, 260)
WRITE(30, 270)
WRITE(30, 280)
DO i=1, NS%Npt
   WRITE(30, 290) i, param%xm(q(i)%id)*(q(i)%xv(1)**2+q(i)%xv(2)**2)/2./kb
END DO
WRITE(30, 300)
!
!
CLOSE(30)

110 FORMAT("GiD Post Results File 1.0")
120 FORMAT("GaussPoints \'2D_test\' ElemType Quadrilateral \'Multi\'")
130 FORMAT("  Number Of Gauss Points: 4")
140 FORMAT("Natural Coordinates: Internal")
150 FORMAT("End gausspoints")
!
160 FORMAT("Result \'Gauss Points Stresses\' \'LOAD ANALYSIS\' 1 Vector OnGaussPoints \'2D_test\'")
170 FORMAT("ComponentNames \'Sxx\', \'Syy\', \'Sxy\'")
180 FORMAT("Values")
190 FORMAT(I5, 3(1X, ES14.4))
200 FORMAT(5X, 3(1X, ES14.4))
!
210 FORMAT("Result \'DISPLACEMENT\' \'LOAD ANALYSIS\' 1 Vector OnNodes")
220 FORMAT("ComponentNames X-DISPLACEMENT, Y-DISPLACEMENT")
230 FORMAT("Values")
240 FORMAT(I5, 1X, ES14.4, 1X, ES14.4)
250 FORMAT("End values")
!
260 FORMAT("Result \'Particle\' \'temperature ANALYSIS\' 1 Vector OnNodes")
270 FORMAT("ComponentNames Temperature") 
280 FORMAT("Values")
290 FORMAT(I5, 1X, ES14.4)
300 FORMAT("End values")
!
time%Ndump = time%Ndump + 1
!
RETURN
END SUBROUTINE POST
!
!
! ######################## POSTprocessing routine #############################
!
SUBROUTINE POST_FEM(NS, n, e, time, param, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(DATE):: time
TYPE(PRMT):: param
REAL*8    :: dt
!
INTEGER:: i, j, k
REAL*8 :: kb
CHARACTER(LEN=14):: NUMBER, FILENAME, RSLTNAME
DATA NUMBER/'0123456789'/
FILENAME = "hybrid00.geo"
RSLTNAME = "hybrid00.res"
kb = 8.617343E-5
!
IF(INT(time%Ndump/10) < 1 ) THEN
   FILENAME(7:7) = NUMBER(1:1)
   FILENAME(8:8) = NUMBER(time%Ndump+1:time%Ndump+1)
   RSLTNAME(7:7) = NUMBER(1:1)
   RSLTNAME(8:8) = NUMBER(time%Ndump+1:time%Ndump+1)
ELSE
   FILENAME(7:7) = NUMBER(INT(time%Ndump/10)+1:INT(time%Ndump/10)+1)
   FILENAME(8:8) = NUMBER(time%Ndump-INT(time%Ndump/10)*10+1: &
        time%Ndump-INT(time%Ndump/10)*10+1)
   RSLTNAME(7:7) = NUMBER(INT(time%Ndump/10)+1:INT(time%Ndump/10)+1)
   RSLTNAME(8:8) = NUMBER(time%Ndump-INT(time%Ndump/10)*10+1: &
        time%Ndump-INT(time%Ndump/10)*10+1)
END IF
!
! ####################### MESH GEOMETRY PRINT #################################
!
OPEN (UNIT=20, FILE=FILENAME)
!
! Node output
WRITE(20, 10)
WRITE(20, 20) time%now
WRITE(20, 30) time%max, time%dump, dt
WRITE(20, 40)
WRITE(20, 50) NS%Nnode
DO i=1, NS%Nnode
   WRITE(20, 60) i, n(i)%xx(1), n(i)%xx(2)
END DO
!
! Element output 
WRITE(20, 70)
WRITE(20, 50) NS%Nelem
DO i=1, NS%Nelem   
   WRITE(20, 80)  i, (e(i)%l(j), j=1, 4)
END DO
!
! Particle output
WRITE(20, 90)
WRITE(20, 50) NS%Npt
!
! Pair
WRITE(20, 90)
WRITE(20, 50) NS%Npair
!
!
10 FORMAT("# hybrid model geometric data")
20 FORMAT("# time max, dum, dt", ES14.4)
30 FORMAT(3(F10.3, 1X))
40 FORMAT("# Node data")
50 FORMAT(I5)
60 FORMAT(I5, 2(ES14.4, 1X))
70 FORMAT("# Element data")
80 FORMAT(5(I5, 1X))
90 FORMAT("# Particle data")
!
CLOSE(20)
!
! ####################### FEA RESULT PRINT #################################
!
OPEN (UNIT=30, FILE=RSLTNAME)
WRITE(30, 110)
WRITE(30, 120)
WRITE(30, 130)
WRITE(30, 140)
WRITE(30, 150)
!
! Stress at Gauss points
WRITE(30, 160)
WRITE(30, 170)
WRITE(30, 180)
DO i=1, NS%Nelem
   WRITE(30,190) i, (e(i)%sigma(1,j), j=1,3)
   WRITE(30,200)    (e(i)%sigma(2,j), j=1,3)
   WRITE(30,200)    (e(i)%sigma(3,j), j=1,3)   
   WRITE(30,200)    (e(i)%sigma(4,j), j=1,3)   
END DO
WRITE(30, 250)
!
! particle temperature
WRITE(30, 260)
WRITE(30, 270)
WRITE(30, 280)
WRITE(30, 300)
!
!
CLOSE(30)

110 FORMAT("GiD Post Results File 1.0")
120 FORMAT("GaussPoints \'2D_test\' ElemType Quadrilateral \'Multi\'")
130 FORMAT("  Number Of Gauss Points: 4")
140 FORMAT("Natural Coordinates: Internal")
150 FORMAT("End gausspoints")
!
160 FORMAT("Result \'Gauss Points Stresses\' \'LOAD ANALYSIS\' 1 Vector OnGaussPoints \'2D_test\'")
170 FORMAT("ComponentNames \'Sxx\', \'Syy\', \'Sxy\'")
180 FORMAT("Values")
190 FORMAT(I5, 3(1X, ES14.4))
200 FORMAT(5X, 3(1X, ES14.4))
!
210 FORMAT("Result \'DISPLACEMENT\' \'LOAD ANALYSIS\' 1 Vector OnNodes")
220 FORMAT("ComponentNames X-DISPLACEMENT, Y-DISPLACEMENT")
230 FORMAT("Values")
240 FORMAT(I5, 1X, ES14.4, 1X, ES14.4)
250 FORMAT("End values")
!
260 FORMAT("Result \'Particle\' \'temperature ANALYSIS\' 1 Vector OnNodes")
270 FORMAT("ComponentNames Temperature") 
280 FORMAT("Values")
290 FORMAT(I5, 1X, ES14.4)
300 FORMAT("End values")
!
time%Ndump = time%Ndump + 1
!
RETURN
END SUBROUTINE POST_FEM
!
!
! ######################## POSTprocessing routine #############################
!
SUBROUTINE POST_MD(NS, q, time, param, dt)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(PTCL):: q(NS%Npt)
TYPE(DATE):: time
TYPE(PRMT):: param
REAL*8    :: dt
!
INTEGER:: i, j, k
REAL*8 :: kb
CHARACTER(LEN=14):: NUMBER, FILENAME, RSLTNAME
DATA NUMBER/'0123456789'/
FILENAME = "hybrid00.geo"
RSLTNAME = "hybrid00.res"
kb = 8.617343E-5
!
IF(INT(time%Ndump/10) < 1 ) THEN
   FILENAME(7:7) = NUMBER(1:1)
   FILENAME(8:8) = NUMBER(time%Ndump+1:time%Ndump+1)
   RSLTNAME(7:7) = NUMBER(1:1)
   RSLTNAME(8:8) = NUMBER(time%Ndump+1:time%Ndump+1)
ELSE
   FILENAME(7:7) = NUMBER(INT(time%Ndump/10)+1:INT(time%Ndump/10)+1)
   FILENAME(8:8) = NUMBER(time%Ndump-INT(time%Ndump/10)*10+1: &
        time%Ndump-INT(time%Ndump/10)*10+1)
   RSLTNAME(7:7) = NUMBER(INT(time%Ndump/10)+1:INT(time%Ndump/10)+1)
   RSLTNAME(8:8) = NUMBER(time%Ndump-INT(time%Ndump/10)*10+1: &
        time%Ndump-INT(time%Ndump/10)*10+1)
END IF
!
! ####################### MESH GEOMETRY PRINT #################################
!
OPEN (UNIT=20, FILE=FILENAME)
!
! Node output
WRITE(20, 10)
WRITE(20, 20) time%now
WRITE(20, 30) time%max, time%dump, dt
WRITE(20, 40)
WRITE(20, 50) NS%Nnode
!
! Element output 
WRITE(20, 70)
WRITE(20, 50) NS%Nelem
!
! Particle output
WRITE(20, 90)
WRITE(20, 50) NS%Npt
DO i=1, NS%Npt
   WRITE(20, 60) i, q(i)%xx(1), q(i)%xx(2)
END DO 
!
! Pair
WRITE(20, 90)
WRITE(20, 50) NS%Npair
!
!
10 FORMAT("# hybrid model geometric data")
20 FORMAT("# time max, dum, dt", ES14.4)
30 FORMAT(3(F10.3, 1X))
40 FORMAT("# Node data")
50 FORMAT(I5)
60 FORMAT(I5, 2(ES14.4, 1X))
70 FORMAT("# Element data")
80 FORMAT(5(I5, 1X))
90 FORMAT("# Particle data")
!
CLOSE(20)
!
! ####################### FEA RESULT PRINT #################################
!
OPEN (UNIT=30, FILE=RSLTNAME)
WRITE(30, 110)
WRITE(30, 120)
WRITE(30, 130)
WRITE(30, 140)
WRITE(30, 150)
!
! Stress at Gauss points
WRITE(30, 160)
WRITE(30, 170)
WRITE(30, 180)
WRITE(30, 250)
!
! particle temperature
WRITE(30, 260)
WRITE(30, 270)
WRITE(30, 280)
DO i=1, NS%Npt
   WRITE(30, 290) i, param%xm(q(i)%id)*(q(i)%xv(1)**2+q(i)%xv(2)**2)/2./kb
END DO
WRITE(30, 300)
!
!
CLOSE(30)

110 FORMAT("GiD Post Results File 1.0")
120 FORMAT("GaussPoints \'2D_test\' ElemType Quadrilateral \'Multi\'")
130 FORMAT("  Number Of Gauss Points: 4")
140 FORMAT("Natural Coordinates: Internal")
150 FORMAT("End gausspoints")
!
160 FORMAT("Result \'Gauss Points Stresses\' \'LOAD ANALYSIS\' 1 Vector OnGaussPoints \'2D_test\'")
170 FORMAT("ComponentNames \'Sxx\', \'Syy\', \'Sxy\'")
180 FORMAT("Values")
190 FORMAT(I5, 3(1X, ES14.4))
200 FORMAT(5X, 3(1X, ES14.4))
!
210 FORMAT("Result \'DISPLACEMENT\' \'LOAD ANALYSIS\' 1 Vector OnNodes")
220 FORMAT("ComponentNames X-DISPLACEMENT, Y-DISPLACEMENT")
230 FORMAT("Values")
240 FORMAT(I5, 1X, ES14.4, 1X, ES14.4)
250 FORMAT("End values")
!
260 FORMAT("Result \'Particle\' \'temperature ANALYSIS\' 1 Vector OnNodes")
270 FORMAT("ComponentNames Temperature") 
280 FORMAT("Values")
290 FORMAT(I5, 1X, ES14.4)
300 FORMAT("End values")
!
time%Ndump = time%Ndump + 1
!
RETURN
END SUBROUTINE POST_MD
!
! ################# POSTprocessing for time plot routine #####################
!
SUBROUTINE TPOST(NS, n, e, q, time)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(PTCL):: q(NS%Npt)
TYPE(DATE):: time
!
INTEGER:: i, j, k
!
!WRITE(15,3000) time%now, (n(j+364)%xv(2),j=1,5)
!WRITE(15,3000) time%now, (n(j+94)%xv(2),j=1,5)
WRITE(15,3000) time%now, (n(j+352)%xv(2),j=1,5)
3000 FORMAT(6(ES14.4, 1X))
!
RETURN
END SUBROUTINE TPOST
!
SUBROUTINE TPOST_FEM(NS, n, e, time)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(NODE):: n(NS%Nnode)
TYPE(ELEM):: e(NS%Nelem)
TYPE(DATE):: time
!
INTEGER:: i, j, k
!
WRITE(15,3000) time%now, n(353)%xv(2), n(354)%xv(2), n(355)%xv(2), n(356)%xv(2), n(357)%xv(2)
3000 FORMAT(6(ES14.4, 1X))
!
RETURN
END SUBROUTINE TPOST_FEM
!
SUBROUTINE TPOST_MD(NS, q, time)
USE DATASET
IMPLICIT NONE
TYPE(NMBR):: NS
TYPE(PTCL):: q(NS%Npt)
TYPE(DATE):: time
!
INTEGER:: i, j, k
!
WRITE(15,3000) time%now, q(958)%xv(2), q(960)%xv(2), q(962)%xv(2), q(964)%xv(2), q(966)%xv(2)
3000 FORMAT(6(ES14.4, 1X))
!
RETURN
END SUBROUTINE TPOST_MD
