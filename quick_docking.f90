MODULE gen_poses
    INTEGER, PARAMETER :: MAX_ATOMS = 10000

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! SUBROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!--------------------------------------------------    
! Rotate coord to X based on Golden ratio method
! Watson & Curtis, Appl. Crystal. (2013) 46, 1171-1177
!--------------------------------------------------

  SUBROUTINE rotate_coord(coords,X,Rot)
    IMPLICIT NONE
    REAL (kind = 8), DIMENSION(1:,1:), INTENT(INOUT) :: coords
    REAL (kind = 8), DIMENSION(1:,1:,1:), INTENT(INOUT) :: X
    REAL (kind = 8), DIMENSION(1:3) :: vec1, vec2, vec3
    REAL (kind = 8), DIMENSION(1:3,1:3) :: Umat
    REAL (kind = 8) :: pi, GR, theta, RotAngle
    INTEGER :: i, j, k, Matshape(1:2), Matshape2(1:3), Npose, m
    INTEGER :: start, end, counter, Natoms
    INTEGER, INTENT(IN) :: Rot
    CHARACTER(len=100) :: out_file_name,fs

    Matshape = Shape(coords)
    Matshape2 = Shape(X)
    X = 0

    ! vec1 is x-axis
    vec1(1) = 1
    vec1(2) = 0
    vec1(3) = 0
    
    Natoms = Matshape(1)
    Npose = Matshape2(1)

    ! If the flag Rot > 1 then rotate the poses about x-axis
    IF (Rot>1) Npose = Npose/Rot
    start = -(Npose - 1)/2
    end =    (Npose - 1)/2
    counter = 1
    pi = 3.1415
    GR = (1 + SQRT(5.0))/2

    DO i = start, end

       vec2(1) = COS(ASIN(2.0*i/Npose)) * COS(2.0*pi*i/GR)
       vec2(2) = COS(ASIN(2.0*i/Npose)) * SIN(2.0*pi*i/GR)
       vec2(3) = 2.0*i/Npose

       ! Umat to rotate vec1 to vec2
       CALL get_rot_matrix(vec1,vec2,Umat)

       DO k = 1, Natoms
          X(counter,k,:) = MATMUL(Umat, coords(k,:))
       END DO
       m = counter
       counter = counter + 1

       ! Rotate the current coordinates(m) about vec2 axis when flag is 1
       IF (Rot>1) THEN
          DO j = 1, Rot-1
             RotAngle = 360.0/Rot
             theta = (j*RotAngle)*3.1415/180.0
             ! Umat to rotate about vec1 by angle theta
             CALL get_self_rot_matrix(vec1,theta,Umat)
             DO k = 1, Natoms
                X(counter,k,:) = MATMUL(Umat, X(m,k,:))
             END DO
             counter = counter + 1
          END DO
       END IF

    END DO

  END SUBROUTINE rotate_coord

!--------------------------------------------------    
! Unit vec perpendicular to vec1, vec2
!--------------------------------------------------

  SUBROUTINE get_rot_axis(u,v,a)
     IMPLICIT NONE
     REAL (kind = 8), DIMENSION(1:), INTENT(IN) :: u, v
     REAL (kind = 8), DIMENSION(1:3), INTENT(OUT) :: a
     REAL (kind = 8), DIMENSION(1:3) :: temp_vec
     ! temp_vec = cross product of vec1, vec2
     temp_vec(1) =  u(2)*v(3) - v(2)*u(3)
     temp_vec(2) = -u(1)*v(3) + v(1)*u(3)
     temp_vec(3) =  u(1)*v(2) - v(1)*u(2)
     CALL get_unit_vector(temp_vec,a)

   END SUBROUTINE get_rot_axis

     
!--------------------------------------------------    
! Angle between vectors
!--------------------------------------------------

  SUBROUTINE get_rot_angle(u,v,theta)
     IMPLICIT NONE
     REAL (kind = 8), DIMENSION(1:), INTENT(IN) :: u, v
     REAL (kind = 8), DIMENSION(1:3) :: a, b
     REAL (kind = 8), INTENT(INOUT):: theta
     REAL :: cos_theta

     CALL get_unit_vector(u,a)
     CALL get_unit_vector(v,b)

     ! angle between u,v
     cos_theta = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
     theta = ACOS(cos_theta)

   END SUBROUTINE get_rot_angle

!--------------------------------------------------    
! get rotation matrix
! to rotate vec1 towards vec2
!--------------------------------------------------

  SUBROUTINE get_rot_matrix(vec1,vec2,Umat)
     IMPLICIT NONE
     REAL (kind = 8), DIMENSION(1:), INTENT(IN) :: vec1, vec2
     REAL (kind = 8), DIMENSION(1:3) :: a
     REAL (kind = 8), DIMENSION(1:,1:), INTENT(OUT) :: Umat
     REAL (kind = 8) :: theta, sin_theta, cos_theta
     REAL (kind = 8) :: u,v,w
     CALL get_rot_axis(vec1,vec2,a)
     CALL get_rot_angle(vec1,vec2,theta)
     CALL get_self_rot_matrix(a,theta,Umat)     

   END SUBROUTINE get_rot_matrix

     
!--------------------------------------------------    
! get self rotation matrix
! i.e. rotation about vector 'a' by theta degrees
!--------------------------------------------------

  SUBROUTINE get_self_rot_matrix(a,theta,Umat)
     IMPLICIT NONE
     REAL (kind = 8), DIMENSION(1:), INTENT(IN) :: a
     REAL (kind = 8), DIMENSION(1:,1:), INTENT(OUT) :: Umat
     REAL (kind = 8), INTENT(IN) :: theta
     REAL (kind = 8) :: sin_theta, cos_theta
     REAL (kind = 8) :: u,v,w

     u = a(1)
     v = a(2)
     w = a(3)

     cos_theta = COS(theta)
     sin_theta = SIN(theta)
     
     Umat(1,1) = u*u + (1 - u*u)*cos_theta
     Umat(1,2) = u*v*(1 - cos_theta) - w*sin_theta
     Umat(1,3) = u*w*(1 - cos_theta) + v*sin_theta

     Umat(2,1) = u*v*(1 - cos_theta) + w*sin_theta
     Umat(2,2) = v*v + (1 - v*v)*cos_theta
     Umat(2,3) = v*w*(1 - cos_theta) - u*sin_theta

     Umat(3,1) = u*w*(1 - cos_theta) - v*sin_theta
     Umat(3,2) = v*w*(1 - cos_theta) + u*sin_theta
     Umat(3,3) = w*w + (1 - w*w)*cos_theta

   END SUBROUTINE get_self_rot_matrix

     
!--------------------------------------------------    
!  Find minimum distance between atoms
!--------------------------------------------------

  SUBROUTINE get_min_dist(X1,X2,mdist)
    IMPLICIT NONE
    REAL (kind = 8), DIMENSION(1:,1:), INTENT(IN) :: X1, X2
    REAL (kind = 8), INTENT(OUT) :: mdist
    REAL (kind = 8) :: dist
    INTEGER :: i,j,Matshape1(1:2), Matshape2(1:2)

    Matshape1 = SHAPE(X1)
    Matshape2 = SHAPE(X2)

    mdist = 1000000.0

    DO i = 1, Matshape1(1)
       DO j = 1, Matshape2(1)
          dist = distance(X1(i,:),X2(j,:))
          IF (dist < mdist) mdist = dist
       END DO
    END DO

  END SUBROUTINE get_min_dist

!--------------------------------------------------    
! Align Coordinates X to coords
!--------------------------------------------------

  SUBROUTINE align(coords,X)
    IMPLICIT NONE
    REAL (kind = 8), DIMENSION(1:,1:), INTENT(INOUT) :: coords
    REAL (kind = 8), DIMENSION(1:,1:), INTENT(INOUT) :: X
    REAL (kind = 8), DIMENSION(1:3,1:3) :: Umat

    CALL rot_matrix(coords,X,Umat)
    X  = TRANSPOSE(MATMUL(Umat,  TRANSPOSE(X)))

  END SUBROUTINE align

!--------------------------------------------------    
! Get unit vector
!--------------------------------------------------    

  SUBROUTINE get_unit_vector(X,unit)
    IMPLICIT NONE
    REAL (kind = 8), DIMENSION(1:3), INTENT(IN) :: X
    REAL (kind = 8), DIMENSION(1:3), INTENT(OUT) :: unit
    REAL (kind = 8) :: DR
    
    DR = sqrt(X(1)*X(1) + X(2)*X(2) + X(3)*X(3))
    
    unit(1) = X(1)/DR
    unit(2) = X(2)/DR
    unit(3) = X(3)/DR
    IF (DR==0) unit = 0

  END SUBROUTINE get_unit_vector
!--------------------------------------------------    
! Get center of mass of coordinates
!--------------------------------------------------    

  SUBROUTINE get_COM(X,COM)
    IMPLICIT NONE
    REAL (kind = 8), DIMENSION(1:,1:), INTENT(IN) :: X
    REAL (kind = 8), DIMENSION(1:3), INTENT(OUT) :: COM
    INTEGER :: Matshape(1:2), i,j
    Matshape = SHAPE(X)
    COM(1) = 0
    COM(2) = 0
    COM(3) = 0
    
    DO i = 1, Matshape(1)
       DO j = 1,3
          COM(j) = COM(j) + X(i,j)
       END DO
    END DO
    COM(1) = COM(1)/Matshape(1)
    COM(2) = COM(2)/Matshape(1)
    COM(3) = COM(3)/Matshape(1)

  END SUBROUTINE get_COM
   
    
!--------------------------------------------------    
! Get aligned coordinates - Kabsch method
!--------------------------------------------------    

  SUBROUTINE rot_matrix(coords1,X,Umat)
    IMPLICIT NONE
    REAL (kind = 8), DIMENSION(1:,1:), INTENT(INOUT) :: coords1
    REAL (kind = 8), DIMENSION(1:,1:), INTENT(INOUT) :: X
    REAL (kind = 8), DIMENSION(1:3,1:3), INTENT(INOUT) :: Umat
    REAL (kind = 8), DIMENSION(1:3,1:3) :: u, vt
    REAL (kind = 8), DIMENSION(1:3,1:3) :: Amat
    REAL (kind = 8), DIMENSION(1:3) ::  s
    REAL (kind = 8) :: work(15)
    REAL (kind = 8) :: detu, detvt 
    INTEGER :: i,j, info

    ! Translate coordinates to origin
    CALL center_at_o(coords1)
    CALL center_at_o(X)

    Amat = MATMUL(TRANSPOSE(coords1), X)

    ! Singular value decomposition of Amat = u * [s] * vt
    ! This function is from standard BLAS/LAPACK(f90) module
    CALL dgesvd ('A','A', 3, 3, Amat, 3, s, u, 3, vt, 3, work, 15, info )

    detu = DETERMINANT(u)
    detvt = DETERMINANT(vt)
    detu = detu*detvt
    IF (detu < 0) THEN
       DO j = 1,3
          u(j,3) = -u(j,3)
       END DO
    END IF
    Umat = MATMUL(u,vt)

  END SUBROUTINE rot_matrix

!--------------------------------------------------    
! Calculate DETERMINANT of 3X3 matrix
!--------------------------------------------------    

  FUNCTION determinant(a)
    IMPLICIT NONE
    REAL (kind = 8), INTENT(IN), DIMENSION(1:3,1:3) :: a
    INTEGER :: i,j
    REAL (kind = 8) :: determinant

    determinant = 0.0
    determinant = determinant + a(1,1) *( a(2,2)*a(3,3) - a(3,2)*a(2,3) )
    determinant = determinant - a(1,2) *( a(2,1)*a(3,3) - a(3,1)*a(2,3) )
    determinant = determinant + a(1,3) *( a(2,1)*a(3,2) - a(3,1)*a(2,2) )

  END FUNCTION determinant

!--------------------------------------------------    
! Read .pdb file to get atom coordinates
!--------------------------------------------------    

  SUBROUTINE read_pdb(file_name,Y,Natoms)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: file_name
    REAL (kind = 8), DIMENSION(:,:), ALLOCATABLE :: X
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: Y
    INTEGER, INTENT(INOUT) :: Natoms
    CHARACTER(len=150) :: dummy
    INTEGER :: irslt, iatom


    ALLOCATE(X(1:MAX_ATOMS,1:3))

    ! atom_index
    iatom = 1 
    OPEN (unit=96, file=trim(file_name), status='old', action='read')
    !write(0,*) " "
    write(0,*) "Reading PDB file: ", file_name
    DO 

       READ (96, '(A100)', iostat = irslt) dummy
       IF (irslt .ne. 0) exit

       IF (dummy(:4) == "ATOM" .OR. dummy(:6) == "HETATM") THEN
          READ (dummy(32:55), '(F8.3,F8.3,F8.3)') X(iatom,1), X(iatom,2), X(iatom,3)
          iatom = iatom + 1
       END IF
    END DO


    iatom = iatom - 1
    CLOSE(96)
    ALLOCATE(Y(1:iatom,1:3))
    Y(1:iatom,1:3) = X(1:iatom,1:3)
    
    Natoms = iatom

    DEALLOCATE(X)
    write(0,'(A30, I5)') "Total atoms in PDB: ", Natoms
  END SUBROUTINE read_pdb

!--------------------------------------------------    
! WRITE .pdb file with X coordinates using input pdb
!--------------------------------------------------    

  SUBROUTINE write_pdb(file_name1,file_name2,out_file_name,coords1,coords2)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: file_name1,file_name2, out_file_name
    REAL (kind = 8), DIMENSION(1:,1:), INTENT(IN) :: coords1, coords2
    CHARACTER(len=100) :: dummy
    INTEGER :: irslt, iatom
    
    ! read and write PDB1
    iatom = 1 ! atom_index
    OPEN (unit=97, file=trim(out_file_name))
    OPEN (unit=96, file=trim(file_name1), status='old', action='read')
    DO 
       READ (96, '(A100)', iostat = irslt) dummy
       IF (irslt .ne. 0) exit

       IF (dummy(:4) == "ATOM" .OR. dummy(:6)=="HETATM") THEN
          WRITE (dummy(31:54), '(F8.3,F8.3,F8.3)') coords1(iatom,1), coords1(iatom,2), coords1(iatom,3)
          WRITE(97,'(A100)') dummy
          iatom = iatom + 1
       END IF
    END DO
    CLOSE(96)
    
    ! read and write PDB2
    iatom = 1 ! atom index
    OPEN (unit=96, file=trim(file_name2), status='old', action='read')
    DO 
       READ (96, '(A100)', iostat = irslt) dummy
       IF (irslt .ne. 0) exit

       IF (dummy(:4) == "ATOM" .OR. dummy(:6)=="HETATM") THEN
          WRITE (dummy(31:54), '(F8.3,F8.3,F8.3)') coords2(iatom,1), coords2(iatom,2), coords2(iatom,3)
          WRITE(97,'(A100)') dummy
          iatom = iatom + 1
       END IF
    END DO

    WRITE(97,'(A)') "TER"
    CLOSE(96)
    CLOSE(97)
    
  END SUBROUTINE write_pdb

!--------------------------------------------------    
! Center coordinates about the origin
!--------------------------------------------------    

  SUBROUTINE center_at_o(a)
    IMPLICIT NONE
    REAL (kind = 8), INTENT(INOUT), DIMENSION(1:,1:) :: a
    INTEGER :: i,j,Matshape(1:2)
    REAL (kind = 8) :: cen(1:3)
    Matshape = SHAPE(a)

    cen(1) = 0
    cen(2) = 0
    cen(3) = 0

    DO i = 1,Matshape(1)
       DO j = 1,Matshape(2)
          cen(j) = cen(j) + a(i,j)
       END DO
    END DO

    cen(1) = cen(1)/Matshape(1)
    cen(2) = cen(2)/Matshape(1)
    cen(3) = cen(3)/Matshape(1)
    !write(*,*) cen
    DO i = 1,Matshape(1)
       DO j = 1,Matshape(2)
          a(i,j) = a(i,j) - cen(j)
       END DO
    END DO
    
  END SUBROUTINE center_at_o

!--------------------------------------------------    
! Find translation distance using maximum dimensions
!--------------------------------------------------    

  SUBROUTINE Max_diagonal(X,Trans_dist)
    IMPLICIT NONE
    REAL (kind = 8), INTENT(IN), DIMENSION(1:,1:) :: X
    REAL (kind = 8), INTENT(OUT) :: Trans_dist
    INTEGER :: Natoms,iatom, Matshape(2)
    ! Maximum X, Y, Z in the PDB
    REAL (kind = 8) :: Xmax,Ymax,Zmax
    Xmax = 0
    Ymax = 0
    Zmax = 0
    Matshape = SHAPE(X)
    Natoms = Matshape(1)

    DO iatom = 1, Natoms
       IF (Xmax < abs(X(iatom,1))) Xmax = abs(X(iatom,1))
       IF (Ymax < abs(X(iatom,2))) Ymax = abs(X(iatom,2))
       IF (Zmax < abs(X(iatom,3))) Zmax = abs(X(iatom,3))
    END DO

    ! max diagonal distance
    Trans_dist = sqrt(Xmax*Xmax + Ymax*Ymax + Zmax*Zmax)

  END SUBROUTINE Max_diagonal

!--------------------------------------------------    
! Distance between two points
!--------------------------------------------------    

  FUNCTION distance(X,Y)
    IMPLICIT NONE
    REAL (kind = 8), DIMENSION(1:3), INTENT(IN) :: X,Y
    REAl (kind = 8) :: distance,r
    
    r =     (X(1)-Y(1))*(X(1)-Y(1))
    r = r + (X(2)-Y(2))*(X(2)-Y(2))
    r = r + (X(3)-Y(3))*(X(3)-Y(3))
    r = sqrt(r)
    distance = r
    END FUNCTION distance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 END MODULE gen_poses


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM MAIN
  USE gen_poses
  IMPLICIT NONE
  REAL (kind = 8), ALLOCATABLE, DIMENSION(:,:) :: coords1, coords2, XX2
  REAL (kind = 8),  ALLOCATABLE, DIMENSION(:,:,:) :: X1, X2
  REAL (kind =8) :: Trans, Trans1, Trans2
  INTEGER :: i, j, k, m, Natoms, irslt
  CHARACTER(len=100) :: inp_file_name1, inp_file_name2, out_file_name
  CHARACTER(len=100) :: cNpose, fs, cRot, cmdist
  REAL (kind = 8) :: cutoff_dist, min_dist
  INTEGER :: counter, Npose1, Rot, Npose2

  ! General output/help
  write(0,*) " "
  write(0,*) "********************************************************** "
  write(0,*) "QUICK DOCKING - Version 1.0 (10/1/2013)"
  write(0,*) "********************************************************** "
  write(0,*) " --- Krishnakumar M Ravikumar (yanglab.case.edu) "
  write(0,*) " "
  write(0,*) "USAGE: ./quick_docking pdb1 pdb2 [n] [rotations] [min_distance]"
  write(0,*) "pdb1 - PDB file of molecule A"
  write(0,*) "pdb2 - PDB file of molecule B"
  write(0,*) "n - odd integer (Default 25) independent rotations for A,B"
  write(0,*) "rotations - integer (Default: 5) number of relative rotations between A,B"
  write(0,*) "min_distance - float (Default: 4.0 Ang) minimum distance between A,B"
  write(0,*) " "
  write(0,*) "NOTES: "
  write(0,*) "A total of n*n*rotations starting conformations will be generated."
  write(0,*) "Output wil be stored as from 1.pdb to XXX.pdb in the current folder."
  write(0,*) "For all-atom docking, 0.0 < min_distance < 2.0 is recommended."
  write(0,*) "********************************************************** "
  write(0,*) " "

  ! If no inputs provided
  IF (iargc() < 2) THEN 
     write(0,*) "ERROR: 2 inputs PDB files needed!!"
     STOP
  END IF

  ! Read input file
  CALL getarg(1,inp_file_name1)
  CALL getarg(2,inp_file_name2)
  IF (iargc() > 2) THEN 
     CALL getarg(3,cNpose)
  ELSE 
     cNpose = "25"
  END IF
  IF (iargc() > 3) THEN 
     CALL getarg(4,cRot)
  ELSE 
     cRot = "5"
  END IF
  IF (iargc() > 4) THEN 
     CALL getarg(5,cmdist)
  ELSE 
     cmdist = "4.0"
  END IF

  READ(cNpose,*) Npose1
  READ(cRot,*) Rot
  READ(cmdist,*) cutoff_dist

  IF (MOD(Npose1,2) == 0) Npose1 = Npose1 + 1

  IF (Rot >= 1) Npose2 = Npose1*Rot
  write(0,'(A,I8,A)') "Generating ", Npose1*Npose2, " pdbs "

  ! read pdbs
  CALL read_pdb(trim(inp_file_name1),coords1,Natoms)
  CALL CENTER_AT_O(coords1)
  CALL Max_diagonal(coords1,Trans1)
  ALLOCATE (X1(1:Npose1,1:Natoms,1:3))

  CALL read_pdb(trim(inp_file_name2),coords2,Natoms)
  CALL CENTER_AT_O(coords2)
  CALL Max_diagonal(coords2,Trans2)
  ALLOCATE (X2(1:Npose2,1:Natoms,1:3))
  ALLOCATE (XX2(1:Natoms,1:3))

  !Rotate coordinates to generate conformations
  CALL rotate_coord(coords1,X1,1)
  CALL rotate_coord(coords2,X2,Rot)

  ! Translate protein B by Trans1 + Trans2 + delta
  Trans = Trans1 + Trans2 + 5.0
  !write(0,'(A,F8.3,A)') " Translating protein B by ", Trans, " Angstroms"
  DO i = 1, Npose2
     DO j = 1, Natoms
        X2(i,j,1) = X2(i,j,1) + Trans
     END DO
  END DO

  counter = 1
  DO i=1,Npose1
     DO j=1,Npose2
        IF (counter <10) fs = '(I1,A)'
        IF (counter >=10) fs = '(I2,A)'
        IF (counter >=100) fs = '(I3,A)'
        IF (counter >=1000) fs = '(I4,A)'
        IF (counter >=10000) fs = '(I5,A)'
        IF (counter >=100000) fs = '(I6,A)'
        XX2(:,:) = X2(j,:,:)
        ! Quick Docking
        DO k = 1, 5
           CALL get_min_dist(X1(i,:,:),XX2(:,:),min_dist)
           IF (min_dist < cutoff_dist) EXIT
           DO m = 1, Natoms
              XX2(m,1) = XX2(m,1) - (min_dist - cutoff_dist + 1.0)/2.0
           END DO
        END DO
        write(out_file_name,fs) counter,".pdb" 
        CALL write_pdb(inp_file_name1, inp_file_name2, out_file_name, X1(i,:,:), XX2(:,:))
        counter = counter + 1
     END DO
  END DO

END PROGRAM MAIN


