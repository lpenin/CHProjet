Module fonction
  Implicit NONE

contains
  !!!   Calcul X=AU a partir des 5 vecteurs diagonaux   !!!
  Subroutine Mat_mul_creux(A,B1,B2,C1,C2,X,U,Nx)
    Real*8,dimension(:),intent(in)::A,B1,B2,C1,C2,U
    Integer,intent(in)::Nx
    Integer::n,i
    Real*8,dimension(:), allocatable,intent(out):: X

    n=size(A)
    Allocate(X(n))
    X(1)=A(1)*U(1)+C1(1)*U(2)+B1(1)*U(1+Nx)
    Do i=2,Nx
      X(i)=A(i)*U(i)+C2(i)*U(i-1)+C1(i)*U(i+1)+B1(i)*U(i+Nx)
    End do
    Do i=Nx+1,N-Nx-1
      X(i)=A(i)*U(i)+C2(i)*U(i-1)+C1(i)*U(i+1)+B1(i)*U(i+Nx)+B2(i)*U(i-Nx)
    End do
    Do i=N-Nx, N-1
      X(i)=A(i)*U(i)+C1(i)*U(i+1)+C2(i)*U(i-1)+B2(i)*U(i-Nx)
    End do
    X(N)=A(N)*U(N)+C2(N-1)*U(N-1)+B2(N-Nx-1)*U(N-Nx-1)
  End subroutine


  !!!   Notre gradient conjugue   !!!
  Subroutine GC_N(A,B1,B2,C1,C2,b,X,size,Nx)
    Integer,intent(in)::size,Nx
    Real*8,dimension(:),intent(in)::b,A,B1,B2,C1,C2
    Real*8,dimension(:),intent(inout)::X
    Real*8,dimension(:),allocatable::Ap,Ax,rk,rk1,p
    Real*8::alpha,beta
    Integer::k

    Allocate( Ax(size), Ap(size), rk(size), rk1(size), p(size))
    call Mat_mul_creux(A,B1,B2,C1,C2,Ax,X,Nx) ! Ax=A*x
    rk=b-Ax
    p=rk
    k=0
    call Mat_mul_creux(A,B1,B2,C1,C2,Ap,p,Nx) ! Ap=A*p

    Do while (sqrt(Dot_product(rk,rk))>0.001d0)
      call Mat_mul_creux(A,B1,B2,C1,C2,Ap,p,Nx) ! Ap=A*p
      alpha=Dot_product(p,rk)/(Dot_product(p,Ap))
      X=X+alpha*p
      rk1=rk-alpha*Ap
      beta=Dot_product(rk1,rk1)/Dot_product(rk,rk)
      p=rk1+beta*p
      k=k+1
      rk=rk1
    End do
    deallocate(Ax,Ap,rk,rk1,p)
  End subroutine


  !!!   Gradient conjugue repris de votre programme   !!!
  Subroutine GC(A,B1,B2,C1,C2,Mat_f,U,n,Nx)
    Integer,intent(in)::n, Nx
    Real*8,dimension(:),intent(in)::A,B1,B2,C1,C2, Mat_f
    Real*8,dimension(:),intent(inout)::U
    Real*8,dimension(:),allocatable:: U_0, r, d1, W
    Real*8:: epsilon, alpha, beta, residu, drl, dwl, betal, residul
    Integer::i, j, l

    Allocate(U_0(1:n), r(1:n), d1(1:n))
    call Mat_mul_creux(A,B1,B2,C1,C2,W,U,Nx) ! W=A*U

    Do i=1,N
      U_0(i) = U(i)
      r(i)     = W(i) - Mat_f(i)
      d1(i)     = W(i) - Mat_f(i)
    End do
    residu = 0.d0
    Do i=1,N
      residu = residu + r(i)*r(i)
    End do

    l=1
    Do WHILE ((l<10000).AND.( SQRT(residu) .ge. 0.0001 ))
      call Mat_mul_creux(A,B1,B2,C1,C2,W,d1,Nx) ! W=A*d1
      drl = 0.0d0
      dwl = 0.0d0
      Do i = 1, N
        drl = drl + d1(i)*r(i)
        dwl = dwl + d1(i)*w(i)
      End do
      alpha = drl/dwl
      Do i=1,N
        U_0(i) = U_0(i) - alpha*d1(i)
        r(i) = r(i) - alpha*W(i)
      End do
      betal=0.d0

      Do i=1,N
        betal=betal+ r(i)*r(i)
      End do
      beta = betal/residu
      Do i=1,N
        d1(i) = r(i) + beta*d1(i)
      End do
      residu = 0.d0
      Do i=1,N
        residu    = residu+r(i)*r(i)
      End do
      l=l+1
    End do
    Do i=1,N
      U(i)=U_0(i)
    End do
  End subroutine


  !!!   Fonction aux limites    !!!
  Function f(x,y,t)
    Real*8,intent(in):: x, y, t
    Real*8:: f
    f=2.0d0*(y - y*y + x - x*x)
    !f=sin(x)+cos(y)
    !f=exp( -1.0d0*pow((x-Lx/2.0d0),2) )*exp( -1.0d0*pow((y-Ly/2.0d0),2) )*cos( 2*atan(1.0D0)*t )
  End function

  Function g(x,y,t)
    Real*8,intent(in):: x, y, t
    Real*8:: g
    g=0.0d0
    !g=sin(x)+cos(y)
  End function

  Function h(x,y,t)
    Real*8,intent(in):: x, y, t
    Real*8:: h
    h=0.0d0
    !h=sin(x)+cos(y)
  End function
End module
