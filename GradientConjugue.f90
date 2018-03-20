Module GC
contains
  function GC(A,B1,B2,C1,C2,x,b,size)
    real*8,intent(in)::size
    real*8,dimension(:),allocatable,intent(in)::x,b,A,B1,B2,C1,C2
    real*8,dimension(:),allocatable::Ap,Ax
    real*8::rk,rk1,p,alpha,beta
    integer::k
    Allocate
    Allocate(x(size), b(size), A(size), B1(size), B2(size), C1(size), C2(size), Ax(size), Ap(size))
    call Mat_mul_creux(A,B1,B2,C1,C2,Ax,x)
    r=b-Ax
    p=rk
    k=0
    call Mat_mul_creux(A,B1,B2,C1,C2,Ap,p)
    do while (sqrt(dotproduct(rk,rk))<0.000001)
      alpha=dotproduct(r,r)/(dotproduct(p,Ap))!A*p
      x=x+alpha*p
      rk1=rk-alpha*Ap
      beta=dotproduct(rk1,rk1)/dotproduct(rk,rk)
      p=rk1+beta*p
      k=k+1
      rk=rk1
    end do
  end function
end Module
