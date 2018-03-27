Module fonction

contains
  function f(x,y,t)
    real*8,intent(in):: x, y
    integer, intent(in):: t
    real*8:: f
    f=2.0d0*(y - y*y + x - x*x)

    !f=sin(x)+cos(y)

    !f=exp( -1.0d0*pow((x-Lx/2.0d0),2) )*exp( -1.0d0*pow((y-Ly/2.0d0),2) )*cos( 2*atan(1.0D0)*t )
  end function

  function g(x,y,t)
    real*8,intent(in):: x, y
    integer, intent(in):: t
    real*8:: g
    g=0.0d0
    !g=sin(x)+cos(y)
  end function

  function h(x,y,t)
    real*8,intent(in):: x, y
    integer, intent(in):: t
    real*8:: h
    h=0.0d0
    !h=sin(x)+cos(y)
  end function

  subroutine Mat_mul_creux(A,B1,B2,C1,C2,X,U) !calcul AU=X
    real*8,dimension(:),intent(in)::A,B1,B2,C1,C2,U
    integer::n
    real*8,dimension(:), allocatable,intent(out):: X
    n=size(A)
    Allocate(X(n))

    do i=4,n-3
      X(i)=A(i)*U(i)+C1(i)*U(i+1)+C2(i)*U(i-1)+B1(i)*U(i+3)+B2(i)*U(i-3)
    end do
    X(1)=A(1)*U(1)+C1(1)*U(2)+B1(1)*U(4)
    X(2)=A(2)*U(2)+C1(2)*U(3)+C2(2)*U(1)+B1(2)*U(5)
    X(3)=A(3)*U(3)+C1(3)*U(4)+C2(3)*U(2)+B1(3)*U(6)

    X(n-2)=A(n-2)*U(n-2)+C1(n-2)*U(n-1)+C2(n-2)*U(n-3)+B2(n-2)*U(n-5)
    X(n-1)=A(n-1)*U(n-1)+C1(n-1)*U(n)+C2(n-1)*U(n-2)+B2(n-1)*U(n-4)
    X(n)=A(n)*U(n)+C2(n)*U(n-1)+B2(n)*U(n-3)
  end subroutine
end module
