PROGRAM para
  use Fonction_para
  IMPLICIT NONE
  include 'mpif.h'
  integer::me,Np,statinfo,i1,iN
  integer,dimension(MPI_STATUS_SIZE)::status

  ! PARAMETRE DU SYSTEME
  Integer:: Nx, Ny, i,Maxiter,j,kt,k,n,l, Nl
  Real*8:: Lx, Ly, dx, dy, dt, tfinal, D, t
  Real*8,dimension(:),allocatable:: U, Mat_f,Ur
  Real*8,dimension(:),allocatable::A,B1,B2,C1,C2, X1
  Real*8:: coeff_a,coeff_b,coeff_c
  Real*8:: t1,t2

  !!!   INITIALISATION    !!!

  ! Lis fichier data
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

  ! Definition des matices
  Allocate(U(1:Nx*Ny), Mat_f(1:Nx*Ny), Ur(1:Nx*Ny))
  Allocate(A(1:Nx*Ny), B1(1:Nx*Ny), B2(1:Nx*Ny), C1(1:Nx*Ny), C2(1:Nx*Ny)) 

  ! Création des variables utiles
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

  ! Initialisation du vecteur U
  U(:)= 0.0d0
  !!!   FIN INITIALISATION    !!!



  !!!   BOUCLE EN TEMPS   !!!

  ! Debut de parallelisation
  Call MPI_init(statinfo)
  call MPI_comm_rank(MPI_comm_world,me,Statinfo)
  call MPI_comm_size(MPI_comm_world,Np,Statinfo)

  ! Chrono
  if (me==0) then
    call CPU_TIME( t1 )
  end if

  Maxiter= int(tfinal/dt)
  DO kt = 1, Maxiter
    ! Construction de Mat_f, le membre de droite
    mat_f=0.d0
    call charge(me,n,Np,i1,iN)

    do k=i1,iN
      i=mod(k-1,Nx)+1
      j=k/Nx
      j=int(j)+1
      Mat_f(k)=U(k)+f(i*dx,j*dy,kt*dt)*dt
      IF (i==1) THEN
        Mat_f(k)=Mat_f(k)-coeff_b*h((i-1)*dx,0.d0,dt*kt)
      end if
      if (i==Nx)THEN
        Mat_f(k)=Mat_f(k)-coeff_b*g((i-1)*dx,Ny*dy,dt*kt)
      endif
      IF (j==1) THEN
        Mat_f(k)=Mat_f(k)-coeff_c*h(i*dx,0.d0,dt*kt)
      endif
      IF (j==Ny) THEN
        Mat_f(k)=Mat_f(k)-coeff_c*h(i*dx,(Ny-1)*dy,dt*kt)
      end if
    END DO

    call GC_para(A,B1,B2,C1,C2,Mat_f,U,n,Nx)

    t = t + dt
  ENDDO

  !!vecteur solution si f=y-y²+x-x² et g=h=0
  ! k=1
  ! ur=0.d0
  ! DO j=1,Ny
  !   DO i=1,Nx
  !     ur(k)=i*dx*(1.d0-i*dx)*j*dy*(1.d0-j*dy)
  !     k=k+1
  !   end do
  ! end do

  ! if (me==0) then
  !   do i=1,n
  !     print*, u(i)
  !   end do
  ! end if

  Deallocate(U,Mat_f,Ur,A,B1,B2,C1,C2)
  if (me==0) then
    call CPU_TIME( t2 )
    print *,t2 - t1
  end if
  call MPI_finalize(statinfo)
End program
