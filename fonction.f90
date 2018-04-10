Module fonction
  Implicit NONE

contains
  subroutine Mat_mul_creux(A,B1,B2,C1,C2,X,U,Nx) !calcul AU=X
    real*8,dimension(:),intent(in)::A,B1,B2,C1,C2,U
    integer,intent(in)::Nx
    integer::n,i
    real*8,dimension(:), allocatable,intent(out):: X

    n=size(A)
    Allocate(X(n))
    X(1)=A(1)*U(1)+C1(1)*U(2)+B1(1)*U(1+Nx)
    do i=2,Nx
      X(i)=A(i)*U(i)+C2(i)*U(i-1)+C1(i)*U(i+1)+B1(i)*U(i+Nx)
    end do
    do i=Nx+1,N-Nx-1
      X(i)=A(i)*U(i)+C2(i)*U(i-1)+C1(i)*U(i+1)+B1(i)*U(i+Nx)+B2(i)*U(i-Nx)
    end do
    do i=N-Nx, N-1
      X(i)=A(i)*U(i)+C1(i)*U(i+1)+C2(i)*U(i-1)+B2(i)*U(i-Nx)
    end do
    X(N)=A(N)*U(N)+C2(N-1)*U(N-1)+B2(N-Nx-1)*U(N-Nx-1)
  end subroutine

  subroutine GC_L(A,B1,B2,C1,C2,b,X,size,Nx)
    integer,intent(in)::size,Nx
    real*8,dimension(:),intent(in)::b,A,B1,B2,C1,C2
    real*8,dimension(:),intent(inout)::X
    real*8,dimension(:),allocatable::Ap,Ax,rk,rk1,p
    real*8::alpha,beta
    integer::k

    Allocate( Ax(size), Ap(size), rk(size), rk1(size), p(size))
    call Mat_mul_creux(A,B1,B2,C1,C2,Ax,X,Nx)
    rk=b-Ax
    p=rk
    k=0
    call Mat_mul_creux(A,B1,B2,C1,C2,Ap,p,Nx)

    do while (sqrt(dot_product(rk,rk))>0.001d0)
      call Mat_mul_creux(A,B1,B2,C1,C2,Ap,p,Nx)
      alpha=dot_product(p,rk)/(dot_product(p,Ap))   !A*p
      X=X+alpha*p
      rk1=rk-alpha*Ap
      beta=dot_product(rk1,rk1)/dot_product(rk,rk)
      p=rk1+beta*p
      k=k+1
      rk=rk1
    end do
    deallocate(Ax,Ap,rk,rk1,p)
  end subroutine

  subroutine GC(A,B1,B2,C1,C2,Mat_f,U,n,Nx)
    Real*8,dimension(:),intent(in)::A,B1,B2,C1,C2, Mat_f
    Real*8,dimension(:),intent(inout)::U
    Integer::n, Nx, i, j, l
    REAL*8:: epsilon, alpha, beta, residu, drl, dwl, betal, residul
    Real*8,dimension(:),allocatable:: U_0, r, d1, W

    Allocate(U_0(1:n), r(1:n), d1(1:n))
    call Mat_mul_creux(A,B1,B2,C1,C2,W,U,Nx) !W=AU

    DO i=1,N
      U_0(i) = U(i)
      r(i)     = W(i) - Mat_f(i)
      d1(i)     = W(i) - Mat_f(i)
    END DO
    residu = 0.d0
    DO i=1,N
      residu = residu + r(i)*r(i)
    ENDDO

    ! boucle du Gradient conjugue
    l=1
    DO WHILE ((l<10000).AND.( SQRT(residu) .ge. 0.0001 ))
      call Mat_mul_creux(A,B1,B2,C1,C2,W,d1,Nx)
      drl = 0.0d0
      dwl = 0.0d0
      DO i = 1, N
        drl = drl + d1(i)*r(i)
        dwl = dwl + d1(i)*w(i)
      ENDDO
      alpha = drl/dwl
      DO i=1,N
        U_0(i) = U_0(i) - alpha*d1(i)
        r(i) = r(i) - alpha*W(i)
      END DO
      betal=0.d0

      DO i=1,N
        betal=betal+ r(i)*r(i)
      ENDDO
      beta = betal/residu
      DO i=1,N
        d1(i) = r(i) + beta*d1(i)
      ENDDO
      residu = 0.d0
      DO i=1,N
        residu    = residu+r(i)*r(i)
      ENDDO
      l=l+1
      !Fin Gradient conjugue
    ENDDO
    DO i=1,N
      U(i)=U_0(i)
    ENDDO
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


end module
