Module fonction

contains
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

  subroutine Mat_mul_creux(A,B1,B2,C1,C2,X,U,Nx) !calcul AU=X
    real*8,dimension(:),intent(in)::A,B1,B2,C1,C2,U
    integer,intent(in)::Nx
    integer::n
    real*8,dimension(:), allocatable,intent(out):: X
    n=size(A)
    Allocate(X(n))

    do i=1,Nx
      X(i)=A(i)*U(i)+C1(i)*U(i+1)+C2(i)*U(i-1)+B1(i)*U(i+Nx)
    end do

    do i=Nx,n-Nx
      X(i)=A(i)*U(i)+C1(i)*U(i+1)+C2(i)*U(i-1)+B1(i)*U(i+Nx)+B2(i)*U(i-Nx)
    end do

    do i=n-Nx,n
      X(i)=A(i)*U(i)+C1(i)*U(i+1)+C2(i)*U(i-1)+B2(i)*U(i-Nx)
    end do

  end subroutine
end module
