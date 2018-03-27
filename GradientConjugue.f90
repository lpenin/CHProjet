Module Gradient
use fonction

contains
  subroutine GC(A,B1,B2,C1,C2,b,X,size)
    integer,intent(in)::size
    real*8,dimension(:),intent(in)::b,A,B1,B2,C1,C2
    real*8,dimension(:),intent(inout)::X
    real*8,dimension(:),allocatable::Ap,Ax,rk,rk1,p
    real*8::alpha,beta
    integer::k

    Allocate( Ax(size), Ap(size), rk(size), rk1(size), p(size))
    call Mat_mul_creux(A,B1,B2,C1,C2,Ax,X)
    rk=b-Ax
    p=rk
    k=0
    call Mat_mul_creux(A,B1,B2,C1,C2,Ap,p)

    do while (sqrt(dot_product(rk,rk))>0.001d0)
      call Mat_mul_creux(A,B1,B2,C1,C2,Ap,p)
      alpha=dot_product(rk,rk)/(dot_product(p,Ap))   !A*p
      X=X+alpha*p
      !print*,"alpha =", alpha
      rk1=rk-alpha*Ap

      beta=dot_product(rk1,rk1)/dot_product(rk,rk)
      p=rk1+beta*p
      k=k+1
      rk=rk1
    end do
    deallocate(Ax,Ap,rk,rk1,p)
  end subroutine



end Module
