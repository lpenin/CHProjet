Program Projet
  use fonction
  Implicit NONE


! PARAMETRE DU SYSTEME
Integer:: Nx, Ny, Maxiter
Real*8:: Lx, Ly, dx, dy, dt, tfinal, D
Real*8,dimension(:),allocatable:: U_0, U, Mat_f
Real*8,dimension(:),allocatable::A,B1,B2,C1,C2
Real*8,Parameter:: coeff_a,coeff_b,coeff_c


! INITIALISATION
OPEN(unit,file='data', form='formatted', status='old')
READ(unit,*)
READ(unit,*) Lx, Ly
READ(unit,*)
READ(unit,*) D
READ(unit,*)
READ(unit,*) Nx, Ny
CLOSE(unit)

Allocate(U_0(1:Nx*Ny), U(1:Nx*Ny), Mat_f(1:Nx*Ny))
Allocate(A(1:Nx*Ny), B1(1:Nx*Ny), B2(1:Nx*Ny), C1(1:Nx*Ny), C2(1:Nx*Ny)) !expliquer d'ou vient la taille des vecteurs

coeff_a=1.0d0+ 2.0d0*D*dt/(dy*dy) + 2.0d0*D*dt/(dx*dx)
coeff_b= -1.0d0*D*dt/(dx*dx)
coeff_c= -1.0d0*D*dt/(dy*dy)

Maxiter= int(tfinal/dt)

!On donne la condition initiale à U
U(:)= 0.0d0

  DO kt = 1, Maxiter
    ! constructio n de Mat_f
    k=1
    DO j=1,Ny
      DO i=1,Nx
        Mat_f(k)=Mat_f(k)+f(i*dx,j*dy,kt))

      IF (i==1) THEN
        Mat_f(k)=Mat_f(k)-coeff_c*h((i-1)*dx,0,kt)

      ELSE if (i==Nx)THEN
        Mat_f(k)=Mat_f(k)-coeff_b*g((i-1)*dx,Ny*dy,kt)

      ELSE IF (j==0) THEN
        Mat_f(k)=Mat_f(k)-coeff_c*h(i*dx,0,kt)

      ELSE IF (j==0) THEN
        Mat_f(k)=Mat_f(k)-coeff_c*h(i*dx,(Ny-1)*dy,kt)

      k=k+1

    ENDIF

    END DO
  END DO


  END DO

Deallocate(U_0,U,Mat_f)


end program
