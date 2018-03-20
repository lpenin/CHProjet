Module GC
contains
  function GC(A,B1,B2,C1,C2,x,b,size)
    real*8,intent(in)::size
    real*8,dimension(:),allocatable,intent(in)::x,b,A,B1,B2,C1,C2
    real*8::rk,rk1,p,alpha,beta
    integer::k

    Allocate(x(size), b(size), A(size), B1(size), B2(size), C1(size), C2(size))
    r=b-!A*x0
    p=rk
    k=0
    do while (sqrt(dotproduct(rk,rk))<0.000001)
      alpha=dotproduct(r,r)/(dotproduct(p,))!A*p
      x=x+alpha*p
      rk1=rk-alpha!*A*p
      beta=dotproduct(rk1,rk1)/dotproduct(rk,rk)
      p=rk1+beta*p
      k=k+1
      rk=rk1
    end do
  end function
end Module
