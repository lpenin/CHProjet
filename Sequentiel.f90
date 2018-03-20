Program Projet
  use fonction
  Implicit NONE


! PARAMETRE DU SYSTEME
Integer:: Nx, Ny
Real*8:: Lx, Ly, dx, dy, dt, tfinal, D
Real*8,dimension(:,:),allocatable:: U_0, U, Mat_f
Real*8,Parameter:: A,B,C


! INITIALISATION
OPEN(unit,file='data', form='formatted', status='old')
READ(unit,*)
READ(unit,*) Lx, Ly
READ(unit,*)
READ(unit,*) D
READ(unit,*)
READ(unit,*) Nx, Ny
CLOSE(unit)

Allocate(U_0(0:Nx+1,0:Ny+1), U(0:Nx+1,0:Ny+1), Mat_f(0:Nx+1,0:Ny+1))

A=1.0d0+ 2.0d0*D*dt/(dy*dy) + 2.0d0*D*dt/(dx*dx)
B= -1.0d0*D*dt/(dx*dx)
C= -1.0d0*D*dt/(dy*dy)




Deallocate(U_0,U,Mat_f)


end program
