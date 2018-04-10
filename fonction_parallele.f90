Module Fonction_para
  Implicit NONE

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

  subroutine GC_para(A,B1,B2,C1,C2,Mat_f,U,n,Nx)
    include 'mpif.h'
    Real*8,dimension(:),intent(in)::A,B1,B2,C1,C2, Mat_f
    Real*8,dimension(:),intent(inout)::U
    Integer::n, Nx, i, j, l
    REAL*8:: epsilon, alpha, beta, residu, drl, dwl, betal, residul, dr, dw, res_sum
    Real*8,dimension(:),allocatable:: U_0, r, d1, W, W1, Up
    integer::me,Np,statinfo,i1,iN
    integer,dimension(MPI_STATUS_SIZE)::status

    call MPI_comm_rank(MPI_comm_world,me,Statinfo)
    call MPI_comm_size(MPI_comm_world,Np,Statinfo)

    Allocate(W(1:n), W1(1:n), U_0(1:n), r(1:n), d1(1:n), Up(1:n))


    ! Initialisation Gradient conjugue
    call charge(me,n,Np,i1,iN)
    call Mat_mul_para(A,B1,B2,C1,C2,W,U,Nx,i1,iN) !W=AU

    DO i=i1,iN
      U_0(i) = U(i)
      r(i)     = W(i) - Mat_f(i)
      d1(i)     = W(i) - Mat_f(i)
    END DO

    residu = 0.d0
    res_sum=0.d0
    DO i=i1,iN
      res_sum = res_sum + r(i)*r(i) !SUM(r*r)
    ENDDO
    call MPI_allreduce(res_sum,residu,1,MPI_real8,MPI_sum,MPI_comm_world,statinfo)

    ! Boucle du Gradient conjugue
    l=1
    DO WHILE ((l<10000).AND.( SQRT(residu) .ge. 0.0001 ))
      call MPI_allreduce(d1,w1,n,MPI_real8,MPI_sum,MPI_comm_world,statinfo)

      call charge(me,n,Np,i1,iN)
      call Mat_mul_para(A,B1,B2,C1,C2,w,w1,Nx,i1,iN)

      drl = 0.0d0
      dwl = 0.0d0
      dr=0.d0
      dw=0.d0
      DO i = i1, iN
        drl = drl + d1(i)*r(i)
        dwl = dwl + d1(i)*w(i)
      ENDDO

      call MPI_allreduce(drl,dr,1,MPI_real8,MPI_sum,MPI_comm_world,statinfo)
      call MPI_allreduce(dwl,dw,1,MPI_real8,MPI_sum,MPI_comm_world,statinfo)

      alpha = dr/dw
      call charge(me,n,Np,i1,iN)
      DO i=i1,iN
        U_0(i) = U_0(i) - alpha*d1(i)
        r(i) = r(i) - alpha*W(i)
      END DO
      betal=0.d0
      DO i=i1,iN
        betal=betal+ r(i)*r(i)
      ENDDO
      call MPI_allreduce(betal/residu,beta,1,MPI_real8,MPI_sum,MPI_comm_world,statinfo)

      DO i=i1,iN
        d1(i) = r(i) + beta*d1(i)
      ENDDO
      residu = 0.d0
      res_sum=0.d0
      DO i=i1,iN
        res_sum    = res_sum+r(i)*r(i)
      ENDDO
      call MPI_allreduce(res_sum,residu,1,MPI_real8,MPI_sum,MPI_comm_world,statinfo)
      l=l+1
    ENDDO

    u=0.d0
    DO i=i1,iN
      U(i)=U_0(i)
    ENDDO

    ! Fin boucle en temps
    call MPI_allreduce(u,up,n,MPI_real8,MPI_sum,MPI_comm_world,statinfo)
    u=up
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
