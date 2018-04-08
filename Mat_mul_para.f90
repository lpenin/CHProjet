PROGRAM matmul
  IMPLICIT NONE
  include 'mpif.h'
  integer::me,Np,statinfo,i1,iN
  integer,dimension(MPI_STATUS_SIZE)::status

  Integer:: Nx, Ny, i,Maxiter,j,kt,k,n
  Real*8:: Lx, Ly, dx, dy, dt, tfinal, D, t
  Real*8,dimension(:),allocatable:: U, x,x1
  Real*8,dimension(:),allocatable::A,B1,B2,C1,C2
  Real*8:: coeff_a,coeff_b,coeff_c

  Call MPI_init(statinfo)
  call MPI_comm_rank(MPI_comm_world,me,Statinfo)
  call MPI_comm_size(MPI_comm_world,Np,Statinfo)


  Nx=5
  Ny=5
  n=Nx*Ny
  Allocate(A(1:Nx*Ny), B1(1:Nx*Ny), B2(1:Nx*Ny), C1(1:Nx*Ny), C2(1:Nx*Ny), X(n), U(n))

  coeff_a=5.
  coeff_b=2.
  coeff_c=-1.
  do i=1,n
    U(i)=i
  end do

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


  call charge(me,n,Np,i1,iN)
  call Mat_mul_para(A,B1,B2,C1,C2,X,U,Nx,i1,iN)


  if (me==0) then
    do i=1,Np-1
      call charge(i,n,Np,i1,iN)
      Allocate(x1(1:iN-i1+1))
      call MPI_recv(x1,iN-i1+1,MPI_REAL8,i,100,MPI_COMM_WORLD,status,statinfo)
      X(i1:iN+1)=x1
      DEALLOCATE(x1)
    end do
  Else
    call MPI_send(X(i1:iN+1),iN-i1+1,MPI_REAL8,0,100,MPI_COMM_WORLD,statinfo)
    print*, 'je suis ', me, ' j''envoie', X
    print*,''
  end if


  call MPI_finalize(statinfo)

  if (me==0) then
    do i=1,n
      print*,x(i)
    end do
  end if

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

End program
