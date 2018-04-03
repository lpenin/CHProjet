Program Projet
  use fonction
  use Gradient
  Implicit NONE


! PARAMETRE DU SYSTEME
Integer:: Nx, Ny, i,Maxiter,j,kt,k,n
Real*8:: Lx, Ly, dx, dy, dt, tfinal, D, t
Real*8,dimension(:),allocatable:: U_0, U, Mat_f,Ur
Real*8,dimension(:),allocatable::A,B1,B2,C1,C2
Real*8:: coeff_a,coeff_b,coeff_c
real:: t1,t2

call CPU_TIME( t1 )


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

Allocate(U_0(1:Nx*Ny), U(1:Nx*Ny), Mat_f(1:Nx*Ny))
Allocate(A(1:Nx*Ny), B1(1:Nx*Ny), B2(1:Nx*Ny), C1(1:Nx*Ny), C2(1:Nx*Ny),Ur(1:Nx*Ny)) !expliquer d'ou vient la taille des vecteurs

dx=Lx/(Nx+1)
dy=Ly/(Ny+1)
dt=0.9d0*dx*dx/(2.d0*D)
coeff_a=1.0d0+ 2.0d0*D*dt/(dy*dy) + 2.0d0*D*dt/(dx*dx)
coeff_b= -1.0d0*D*dt/(dx*dx)
coeff_c= -1.0d0*D*dt/(dy*dy)

n=Nx*Ny
!! matrice
B1=0.d0
B2=0.d0
A=coeff_a
do i=1,n
  if (i<n-Nx+1) then
    B1(i)=coeff_b
  end if
  if (i>Nx)then
    B2(i)=coeff_b
  end if
end do



C1=coeff_c
C2=coeff_c
do i=1,Nx
 C1(i*Nx)=0.d0
 C2((i-1)*Nx+1)=0.d0
end do


Maxiter= int(tfinal/dt)
print*, Maxiter
!On donne la condition initiale à U
U(:)= 0.0d0
mat_f=0.d0
  DO kt = 1, Maxiter
    ! construction de Mat_f

    k=1
    DO j=1,Ny
      DO i=1,Nx

        Mat_f(k)=f(i*dx,j*dy,kt*dt)

      IF (i==1) THEN
        Mat_f(k)=Mat_f(k)-coeff_b*h((i-1)*dx,0.d0,dt*kt)

      ELSE if (i==Nx)THEN
        Mat_f(k)=Mat_f(k)-coeff_b*g((i-1)*dx,Ny*dy,dt*kt)

      ELSE IF (j==0) THEN
        Mat_f(k)=Mat_f(k)-coeff_c*h(i*dx,0.d0,dt*kt)

      ELSE IF (j==0) THEN
        Mat_f(k)=Mat_f(k)-coeff_c*h(i*dx,(Ny-1)*dy,dt*kt)

      ENDIF
      k=k+1

      END DO
    END DO

  call GC(A,B1,B2,C1,C2,Mat_f,U,n)

  END DO
!  print*,u

  !!vecteur solution
  k=1
  ur=0.d0
  DO j=1,Ny
    DO i=1,Nx

      ur(k)=i*dx*(1.d0-i*dx)*j*dy*(1.d0-j*dy)
      k=k+1

    end do
  end do


  k=1
 DO j=1,Ny
   do i=1,Nx
   PRINT*, i*dx, j*dy, u(k), ur(k)
   !print*, C1(j),k

   k=k+1
 END DO
 end do





Deallocate(U_0,U,Mat_f)
call CPU_TIME( t2 )
! print *,t2 - t1

end program
