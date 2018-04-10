PROGRAM para
  IMPLICIT NONE
  include 'mpif.h'
  integer::me,Np,statinfo,i1,iN
  integer,dimension(MPI_STATUS_SIZE)::status

  ! PARAMETRE DU SYSTEME
  Integer:: Nx, Ny, i,Maxiter,j,kt,k,n,l, Nl
  Real*8:: Lx, Ly, dx, dy, dt, tfinal, D, t
  Real*8,dimension(:),allocatable:: U_0, U, Mat_f,Ur, kappa, r, d1,W
  Real*8,dimension(:),allocatable::A,B1,B2,C1,C2, X1
  Real*8:: coeff_a,coeff_b,coeff_c
  real:: t1,t2
  REAL*8:: epsilon,alpha,beta,residu, drl, dwl, betal, residul, dr, dw, res_sum


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

  Allocate(U_0(1:Nx*Ny), U(1:Nx*Ny), Mat_f(1:Nx*Ny), W(1:Nx*Ny), kappa(1:Nx*Ny), r(1:Nx*Ny), d1(1:Nx*Ny))
  Allocate(A(1:Nx*Ny), B1(1:Nx*Ny), B2(1:Nx*Ny), C1(1:Nx*Ny), C2(1:Nx*Ny),Ur(1:Nx*Ny)) !expliquer d'ou vient la taille des vecteurs

  dx=Lx/(Nx+1)
  dy=Ly/(Ny+1)
  !dt=0.9d0*dx*dx/(4.d0*D)
  dt=1.d0
  coeff_a=1.0d0+ 2.0d0*D*dt/(dy*dy) + 2.0d0*D*dt/(dx*dx)
  coeff_b= -1.0d0*D*dt/(dx*dx)
  coeff_c= -1.0d0*D*dt/(dy*dy)
  n=Nx*Ny

  !! Matrice
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



  Call MPI_init(statinfo)
  call MPI_comm_rank(MPI_comm_world,me,Statinfo)
  call MPI_comm_size(MPI_comm_world,Np,Statinfo)

  if (me ==0) then
    call CPU_TIME( t1 )
  end if

  Maxiter= int(tfinal/dt)

  !On donne la condition initiale à U
  U(:)= 0.0d0
  mat_f=0.d0

  ! construction de Mat_f
  DO kt = 1, Maxiter
    k=1
    DO j=1,Ny
      DO i=1,Nx
        Mat_f(k)=U(k)+f(i*dx,j*dy,kt*dt)*dt

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


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!Gradient Conjugué
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !initialisation Gradient conjugue

    call charge(me,n,Np,i1,iN)
    call Mat_mul_para(A,B1,B2,C1,C2,W,U,Nx,i1,iN) !W=AU

    DO i=i1,iN
      kappa(i) = U(i)
      r(i)     = W(i) - Mat_f(i)
      d1(i)     = W(i) - Mat_f(i)
    END DO

    residu = 0.d0
    res_sum=0.d0
    DO i=i1,iN
      res_sum = res_sum + r(i)*r(i) !SUM(r*r)
    ENDDO
    call MPI_allreduce(res_sum,residu,1,MPI_integer,MPI_sum,MPI_comm_world,statinfo)


    ! boucle du Gradient conjugue

    l=1
    DO WHILE ((l<10000).AND.( SQRT(residu) .ge. 0.0001 ))
      call charge(me,n,Np,i1,iN)
      DO j=1,N
        call Mat_mul_para(A,B1,B2,C1,C2,W,d1,Nx,i1,iN)
      ENDDO

      drl = 0.0d0
      dwl = 0.0d0
      DO i = i1, iN
        drl = drl + d1(i)*r(i)
        dwl = dwl + d1(i)*w(i)
      ENDDO

      call MPI_allreduce(drl,dr,1,MPI_integer,MPI_sum,MPI_comm_world,statinfo)
      call MPI_allreduce(dwl,dw,1,MPI_integer,MPI_sum,MPI_comm_world,statinfo)

      alpha = dr/dw

      call charge(me,n,Np,i1,iN)
      DO i=i1,iN
        kappa(i) = kappa(i) - alpha*d1(i)
        r(i) = r(i) - alpha*W(i)
      END DO

      betal=0.d0
      DO i=i1,iN
        betal=betal+ r(i)*r(i)
      ENDDO
      call MPI_allreduce(betal/residu,beta,1,MPI_integer,MPI_sum,MPI_comm_world,statinfo)

      DO i=i1,iN
        d1(i) = r(i) + beta*d1(i)
      ENDDO

      residu = 0.d0
      res_sum=0.d0
      DO i=i1,iN
        res_sum    = res_sum+r(i)*r(i)
      ENDDO
      call MPI_allreduce(res_sum,residu,1,MPI_integer,MPI_sum,MPI_comm_world,statinfo)
      l=l+1
    ENDDO

    if (me==0) then
    print*,'l',l,'residu',residu
  end if
    DO i=1,N
      U(i)=kappa(i)
    ENDDO
    t = t + dt;

    ! Fin boucle en temps
  ENDDO

  !call GC(A,B1,B2,C1,C2,Mat_f,U,n,Nx)

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

  !
  k=1
  DO j=1,Ny
    do i=1,Nx
      if (me==0) then
      print*, u(k),ur(k)
      !print*, C1(j),k
    end if
      k=k+1
    END DO
  end do



  Deallocate(U_0,U,Mat_f)
  call CPU_TIME( t2 )
  ! print *,t2 - t1

  call MPI_finalize(statinfo)


Contains
  subroutine charge(me,N,Np,i1,iN)
    integer, intent(in):: me, N, Np
    integer,intent(out)::i1, iN
    i1=me*N/Np+1
    iN=(me+1)*N/Np
  end subroutine


  subroutine Mat_mul_para(A,B1,B2,C1,C2,X,U,Nx,i1,iN) !calcul AU=X
    real*8,dimension(:),intent(in)::A,B1,B2,C1,C2,U
    integer,intent(in)::Nx,i1,iN
    integer::n,i
    real*8,dimension(:), allocatable,intent(out):: X
    n=size(A)
    Allocate(X(n))
    X=0.d0
    do i=i1,iN

      if (i==1) then
        X(1)=A(1)*U(1)+C1(1)*U(2)+B1(1)*U(1+Nx)
      end if
      if (i>1 .and. i<Nx+1) then
        X(i)=A(i)*U(i)+C2(i)*U(i-1)+C1(i)*U(i+1)+B1(i)*U(i+Nx)
      end if
      if (i>Nx .and. i<N-Nx) then
        X(i)=A(i)*U(i)+C2(i)*U(i-1)+C1(i)*U(i+1)+B1(i)*U(i+Nx)+B2(i)*U(i-Nx)
      end if
      if (i>N-Nx-1 .and. i<N) then
        X(i)=A(i)*U(i)+C1(i)*U(i+1)+C2(i)*U(i-1)+B2(i)*U(i-Nx)
      end if
      if (i==n) then
        X(N)=A(N)*U(N)+C2(N-1)*U(N-1)+B2(N-Nx-1)*U(N-Nx-1)
      end if
    end do
  end subroutine

  function f(x,y,t)
    real*8,intent(in):: x, y, t
    real*8:: f
    f=2.0d0*(y - y*y + x - x*x)
    !f=sin(x)+cos(y)
    !f=exp( -1.0d0*pow((x-Lx/2.0d0),2) )*exp( -1.0d0*pow((y-Ly/2.0d0),2) )*cos( 2*atan(1.0D0)*t )
  end function

  function g(x,y,t)
    real*8,intent(in):: x, y, t
    real*8:: g
    g=0.0d0
    !g=sin(x)+cos(y)
  end function

  function h(x,y,t)
    real*8,intent(in):: x, y, t
    real*8:: h
    h=0.0d0
    !h=sin(x)+cos(y)
  end function

End program
