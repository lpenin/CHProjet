Module Fonction_para
  Implicit NONE

Contains
  !!!   Calcul la charge de tache, nbr de ligne de l'équation, de chaqu'un des proc   !!!
  subroutine charge(me,N,Np,i1,iN)
    Integer, intent(in):: me, N, Np
    Integer,intent(out)::i1, iN
    i1=me*N/Np+1
    iN=(me+1)*N/Np
  end subroutine


  !!!   Calcul les lignes i1 a iN du produit A*U a partir des 5 vecteurs diagonaux   !!!
  subroutine Mat_mul_para(A,B1,B2,C1,C2,X,U,Nx,i1,iN) !calcul AU=X
    Real*8,dimension(:),intent(in)::A,B1,B2,C1,C2,U
    Integer,intent(in)::Nx,i1,iN
    Integer::n,i
    Real*8,dimension(:),allocatable,intent(out):: X
    n=size(A)
    Allocate(X(n))
    X=0.d0
    Do i=i1,iN

      If (i==1) then
        X(1)=A(1)*U(1)+C1(1)*U(2)+B1(1)*U(1+Nx)
      End if
      If (i>1 .and. i<Nx+1) then
        X(i)=A(i)*U(i)+C2(i)*U(i-1)+C1(i)*U(i+1)+B1(i)*U(i+Nx)
      End if
      If (i>Nx .and. i<N-Nx) then
        X(i)=A(i)*U(i)+C2(i)*U(i-1)+C1(i)*U(i+1)+B1(i)*U(i+Nx)+B2(i)*U(i-Nx)
      End if
      If (i>N-Nx-1 .and. i<N) then
        X(i)=A(i)*U(i)+C1(i)*U(i+1)+C2(i)*U(i-1)+B2(i)*U(i-Nx)
      End if
      If (i==n) then
        X(N)=A(N)*U(N)+C2(N-1)*U(N-1)+B2(N-Nx-1)*U(N-Nx-1)
      End if
    End do
  End subroutine


  !!!   Gradient conjugue parallelise a partir de votre algorithme   !!!
  subroutine GC_para(A,B1,B2,C1,C2,Mat_f,U,n,Nx)
    include 'mpif.h'
    Integer,intent(in)::n, Nx
    Real*8,dimension(:),intent(in)::A,B1,B2,C1,C2, Mat_f
    Real*8,dimension(:),intent(inout)::U
    Real*8,dimension(:),allocatable:: U_0, r, d1, W, W1, Up
    Real*8:: epsilon, alpha, beta, residu, drl, dwl, betal, residul, dr, dw, res_sum
    Integer::i, j, l

    Integer::me,Np,statinfo,i1,iN
    Integer,dimension(MPI_STATUS_SIZE)::status

    call MPI_comm_rank(MPI_comm_world,me,Statinfo)
    call MPI_comm_size(MPI_comm_world,Np,Statinfo)

    Allocate(W(1:n), W1(1:n), U_0(1:n), r(1:n), d1(1:n), Up(1:n))

    call charge(me,n,Np,i1,iN)
    call Mat_mul_para(A,B1,B2,C1,C2,W,U,Nx,i1,iN) ! W=A*U

    Do i=i1,iN
      U_0(i) = U(i)
      r(i)     = W(i) - Mat_f(i)
      d1(i)     = W(i) - Mat_f(i)
    End do

    residu = 0.d0
    res_sum=0.d0
    Do i=i1,iN
      residu= residu + r(i)*r(i)
    End do

    If (me==0) then
      Do i=1,Np-1
        call MPI_recv(res_sum,1,MPI_real8,i,10,MPI_COMM_WORLD,status,statinfo)
        residu=residu+res_sum
      End do
    Else
      call MPI_send(residu,1,MPI_real8,0,10,MPI_COMM_WORLD,statinfo)
    End if
    call MPI_bcast(residu,1,MPI_real8,0,MPI_comm_world,statinfo)

    l=1
    Do while((l<10000).AND.( SQRT(residu) .ge. 0.0001 ))
      call MPI_allreduce(d1,w1,n,MPI_real8,MPI_sum,MPI_comm_world,statinfo)

      call charge(me,n,Np,i1,iN)
      call Mat_mul_para(A,B1,B2,C1,C2,w,w1,Nx,i1,iN)

      drl = 0.0d0
      dwl = 0.0d0
      dr=0.d0
      dw=0.d0
      Do i = i1, iN
        drl = drl + d1(i)*r(i)
        dwl = dwl + d1(i)*w(i)
      End do

      If (me==0) then
        Do i=1,Np-1
          call MPI_recv(dr,1,MPI_real8,i,100,MPI_COMM_WORLD,status,statinfo)
          call MPI_recv(dw,1,MPI_real8,i,200,MPI_COMM_WORLD,status,statinfo)
          drl=drl+dr
          dwl=dwl+dw
        End do
        alpha = drl/dwl
      Else
        call MPI_send(drl,1,MPI_real8,0,100,MPI_COMM_WORLD,statinfo)
        call MPI_send(dwl,1,MPI_real8,0,200,MPI_COMM_WORLD,statinfo)
      End if
      call MPI_bcast(alpha,1,MPI_real8,0,MPI_comm_world,statinfo)

      call charge(me,n,Np,i1,iN)
      Do i=i1,iN
        U_0(i) = U_0(i) - alpha*d1(i)
        r(i) = r(i) - alpha*W(i)
      End do
      betal=0.d0
      Do i=i1,iN
        beta=beta+ r(i)*r(i)
      End do

      If (me==0) then
        Do i=1,Np-1
          call MPI_recv(betal,1,MPI_real8,i,300,MPI_COMM_WORLD,status,statinfo)
          beta=betal+beta
        End do
      Else
        call MPI_send(beta/residu,1,MPI_real8,0,300,MPI_COMM_WORLD,statinfo)
      End if
      call MPI_bcast(beta,1,MPI_real8,0,MPI_comm_world,statinfo)

      Do i=i1,iN
        d1(i) = r(i) + beta*d1(i)
      End do
      residu = 0.d0
      res_sum=0.d0
      Do i=i1,iN
        res_sum    = res_sum+r(i)*r(i)
      End do

      If (me==0) then
        Do i=1,Np-1
          call MPI_recv(res_sum,1,MPI_real8,i,20,MPI_COMM_WORLD,status,statinfo)
          residu=residu+res_sum
        End do
      Else
        call MPI_send(residu,1,MPI_real8,0,20,MPI_COMM_WORLD,statinfo)
      End if
      call MPI_bcast(residu,1,MPI_real8,0,MPI_comm_world,statinfo)

      l=l+1
    End do

    u=0.d0
    Do i=i1,iN
      U(i)=U_0(i)
    End do

    call MPI_allreduce(u,up,n,MPI_real8,MPI_sum,MPI_comm_world,statinfo) !!! REGLER LE PB DE MATMUL PARA ET ON SUPRIME
    u=up
  end subroutine


  !!!   Fonction aux limites    !!!
  function f(x,y,t)
    Real*8,intent(in):: x, y, t
    Real*8:: f
    f=2.0d0*(y - y*y + x - x*x)
    !f=sin(x)+cos(y)
    !f=exp( -1.0d0*pow((x-Lx/2.0d0),2) )*exp( -1.0d0*pow((y-Ly/2.0d0),2) )*cos( 2*atan(1.0D0)*t )
  end function

  function g(x,y,t)
    Real*8,intent(in):: x, y, t
    Real*8:: g
    g=0.0d0
    !g=sin(x)+cos(y)
  end function

  function h(x,y,t)
    Real*8,intent(in):: x, y, t
    Real*8:: h
    h=0.0d0
    !h=sin(x)+cos(y)
  end function
end module
