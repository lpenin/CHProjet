Program Projet
  use fonction
  use Gradient
  Implicit NONE


! PARAMETRE DU SYSTEME
Integer:: Nx, Ny, i
Real*8:: Lx, Ly, dx, dy, dt, tfinal, D, t
Real*8,dimension(:,:),allocatable:: U_0, U, Mat_f
Real*8,dimension(:),allocatable::A,B1,B2,C1,C2
Real*8:: coeff_a,coeff_b,coeff_c


! INITIALISATION

OPEN(10,file='data', form='formatted', status='old')
READ(10,*)
READ(10,*) Lx, Ly
READ(10,*)
READ(10,*) D
READ(10,*)
READ(10,*) Nx, Ny
READ(10,*)
READ(10,*) tfinal
CLOSE(10)


Allocate(U_0(0:Nx+1,0:Ny+1), U(0:Nx+1,0:Ny+1), Mat_f(0:Nx+1,0:Ny+1))
Allocate(A(1:Nx*Ny), B1(1:Nx*Ny), B2(1:Nx*Ny), C1(1:Nx*Ny), C2(1:Nx*Ny)) !expliquer d'ou vient la taille des vecteurs

coeff_a=1.0d0+ 2.0d0*D*dt/(dy*dy) + 2.0d0*D*dt/(dx*dx)
coeff_b= -1.0d0*D*dt/(dx*dx)
coeff_c= -1.0d0*D*dt/(dy*dy)


do while (t<tfinal)
  !BIIIITTTEEEE
  !Inserez le GC


  t=t+dt
end do




!! matrice
A=coeff_a
do i=1,Nx*Ny
  if (i<Nx*Ny-Nx) then
    B1(i)=coeff_b
  end if
  if (i>Nx)then
    B2(i)=coeff_b
  end if
end do

C1=coeff_c
C2=coeff_c
do i=1,Nx
  C1(i*Nx)=0.

  C2((i-1)*Nx+1)=0.

end do
 !! fin INITIALISATION matrice


Deallocate(U_0,U,Mat_f)


end program
