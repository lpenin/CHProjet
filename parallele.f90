PROGRAM para
  use Fonction_para
  Implicit NONE
  Include 'mpif.h'
  Integer::me,Np,statinfo,i1,iN
  Integer,dimension(MPI_STATUS_SIZE)::status

  ! PARAMETRE DU SYSTEME
  Integer:: Nx, Ny, i,Maxiter,j,kt,k,n,l, Nl
  Real*8:: Lx, Ly, dx, dy, dt, tfinal, D, t
  Real*8,dimension(:),allocatable:: U, Mat_f,Ur
  Real*8,dimension(:),allocatable::A,B1,B2,C1,C2, X1
  Real*8:: coeff_a,coeff_b,coeff_c
  Real*8:: t1,t2

  !!!   INITIALISATION    !!!

  ! Lis le fichier data
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

  ! Definition des matrices
  Allocate(U(1:Nx*Ny), Mat_f(1:Nx*Ny), Ur(1:Nx*Ny))
  Allocate(A(1:Nx*Ny), B1(1:Nx*Ny), B2(1:Nx*Ny), C1(1:Nx*Ny), C2(1:Nx*Ny))

  ! Creation des variables utiles
  dx=Lx/(Nx+1)
  dy=Ly/(Ny+1)
  dt=0.9d0*dx*dx/(4.d0*D)
  coeff_a=1.0d0 + 2.0d0*D*dt/(dy*dy) + 2.0d0*D*dt/(dx*dx)
  coeff_b= -1.0d0*D*dt/(dx*dx)
  coeff_c= -1.0d0*D*dt/(dy*dy)
  n=Nx*Ny

  ! Remplissage des 5 vecteurs representant la matrice
  A=coeff_a

  B1=0.d0
  B2=0.d0
  Do i=1,n
    If (i<n-Nx+1) then
      B1(i)=coeff_b
    End if
    If (i>Nx)then
      B2(i)=coeff_b
    End if
  End do

  C1=coeff_c
  C2=coeff_c
  Do i=1,Nx
    C1(i*Nx)=0.d0
    C2((i-1)*Nx+1)=0.d0
  End do

  ! Initialisation du vecteur U
  U(:)= 0.0d0
  !!!   FIN INITIALISATION    !!!


  ! Debut de parallelisation
  Call MPI_init(statinfo)
  call MPI_comm_rank(MPI_comm_world,me,Statinfo)
  call MPI_comm_size(MPI_comm_world,Np,Statinfo)

  ! Demmarage du chronometre
  If (me==0) then
    call CPU_TIME( t1 )
  End if


  !!!   BOUCLE EN TEMPS   !!!
  Maxiter= int(tfinal/dt)

  Do kt = 1, Maxiter
    ! Construction de Mat_f, le membre de droite
    mat_f=0.d0
    call charge(me,n,Np,i1,iN)

    Do k=i1,iN
      i=mod(k-1,Nx)+1
      j=k/Nx+1
      Mat_f(k)=U(k)+f(i*dx,j*dy,kt*dt)*dt
      If (i==1) then
        Mat_f(k)=Mat_f(k)-coeff_b*h((i-1)*dx,0.d0,dt*kt)
      Endif
      If (i==Nx)then
        Mat_f(k)=Mat_f(k)-coeff_b*g((i-1)*dx,Ny*dy,dt*kt)
      Endif
      If (j==1) then
        Mat_f(k)=Mat_f(k)-coeff_c*h(i*dx,0.d0,dt*kt)
      Endif
      If (j==Ny) then
        Mat_f(k)=Mat_f(k)-coeff_c*h(i*dx,(Ny-1)*dy,dt*kt)
      Endif
    End do

    call GC_para(A,B1,B2,C1,C2,Mat_f,U,n,Nx)

    t = t + dt
  End do
  !!!   FIN BOUCLE EN TEMPS   !!!


  ! ! Vecteur solution si f=y-y²+x-x² et g=h=0
  ! k=1
  ! ur=0.d0
  ! DO j=1,Ny
  !   DO i=1,Nx
  !     ur(k)=i*dx*(1.d0-i*dx)*j*dy*(1.d0-j*dy)
  !     k=k+1
  !   end do
  ! end do

  ! Fin du chronometre
  If (me==0) then
    call CPU_TIME( t2 )
    print *,t2 - t1
  End if

  Deallocate(U,Mat_f,Ur,A,B1,B2,C1,C2)
  call MPI_finalize(statinfo)
End program
