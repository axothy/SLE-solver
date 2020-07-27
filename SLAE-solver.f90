module reflection_QR
contains

subroutine otr(A,B,X,n)
implicit none
    integer                      :: n
    integer                      :: i, j, k
    real                         :: Alpha, p, ts, t
    real(4), dimension (:,:)     :: A
    real(4), dimension (:)       :: B, X
    real(4), dimension (size(X)) :: s, w
    
    do i=1,n-1
        do k=i,n
            s(k-i+1)=A(k,i)
        end do
        
        Alpha=0
        
        do k=1,(n-i+1)
            Alpha=Alpha+s(k)*s(k)
        end do
        
        Alpha=sqrt(Alpha)
        
        do k=1,(n-i+1)
            w(k)=s(k)
        end do
        
        w(1)=w(1)-Alpha
        p=0
        
        do k=1,(n-i+1)
            p=p+w(k)*w(k)
        end do
        
        p=sqrt(p)
        
        do k=1,(n-i+1)
            w(k)=w(k)/p
        end do

        A(i,i)=Alpha

        do j=i+1,n
            ts=0
            do k=1,(n-i+1)
                ts=ts+w(k)*A(i+k-1,j)
            end do
            do k=1,(n-i+1)
                A(i+k-1,j)=A(i+k-1,j)-2*ts*w(k)
            end do
        end do

        ts=0
        do k=1,(n-i+1)
            ts=ts+w(k)*B(i+k-1)
        end do
        do k=1,(n-i+1)
            B(i+k-1)=B(i+k-1)-2*ts*w(k)
        end do
    end do

    X(n)=B(n)/A(n,n)
    do i=(n-1),1,-1
        t=0
        do k=(i+1),n
            t=t+A(i,k)*x(k)
        end do
        X(i)=(B(i)-t)/A(i,i)
    end do
end subroutine otr

end module reflection_QR

program SLE_Solver
use reflection_QR
implicit none
    integer              :: n
    integer              :: i,j
    real(4), allocatable :: A(:,:), B(:)
    real(4), allocatable :: X(:)
    character(len=10)    :: answer
    
    do while (.true.)
    print*, 'Vvedite razmer matritsy n x n:'
    print*, 'Input matrixe size n x n:'
    read*, n

    allocate(A(n,n),B(n),X(n))
    print*, 'Vvedite matritsu koefficientov A po strokam:'
    print*, 'Input matrix of koefficients A by lines:'
    do i=1,n
      read(*,*) (A(i,j),j=1,n)
    end do

    print*, 'Vvedite stolbec svobodnih chlenov B:'
    print*, 'Input the absoulte terms of column B (right side of equations):'
    read(*,*) (B(i),i=1,n)

    call otr(A,B,X,n)
    print*, '______________________________________________'
    print*, 'Reshenie SLAU metodom otrazheniy (Housholdera)'
    print*, 'Housholder reflection SLE solution'
    
    print*, 'Korni:'
    print*, 'Roots:'
    do i=1,n
        print*, X(i)
    end do
    
    print*, 'Hotite prodolzhit rabotu s programmoy? (y/n)'
    print*, 'Do you want to continue working with program? (y/n)'
    read*, answer 
    if (answer=='y') then
        cycle
    else
        stop
    end if
    deallocate(A,B,X)
    end do
end program SLE_Solver
