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
    DO i=1,Nx
      DO j=1,Ny
        Mat_f(k)=Mat_f()+U()+f(i*dx,j*dy,kt))
        k=K=&
      IF (i==1) DO

    END DO
  END DO


  END DO

Deallocate(U_0,U,Mat_f)


end program
