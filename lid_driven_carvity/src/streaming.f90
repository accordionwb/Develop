subroutine streaming(f,n,m)
real f(0:8,0:n,0:m)
! streaming
DO j=0,m
DO i=n,1,-1 !RIGHT TO LEFT
f(1,i,j)=f(1,i-1,j)
END DO
DO i=0,n-1 !LEFT TO RIGHT
f(3,i,j)=f(3,i+1,j)
END DO
END DO
DO j=m,1,-1 !TOP TO BOTTOM
DO i=0,n
f(2,i,j)=f(2,i,j-1)
END DO
DO i=n,1,-1
f(5,i,j)=f(5,i-1,j-1)
END DO
DO i=0,n-1
f(6,i,j)=f(6,i+1,j-1)
END DO
END DO
DO j=0,m-1 !BOTTOM TO TOP
DO i=0,n
f(4,i,j)=f(4,i,j+1)
END DO
DO i=0,n-1
f(7,i,j)=f(7,i+1,j+1)
END DO
DO i=n,1,-1
f(8,i,j)=f(8,i-1,j+1)
END DO
END DO
return
end

