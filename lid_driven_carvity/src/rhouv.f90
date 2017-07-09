subroutine rhouv(f,rho,u,v,cx,cy,n,m)
real f(0:8,0:n,0:m),rho(0:n,0:m),u(0:n,0:m),v(0:n,0:m),cx(0:8),cy(0:8)
do j=0,m
do i=0,n
ssum=0.0
do k=0,8
ssum=ssum+f(k,i,j)
end do
rho(i,j)=ssum
end do
end do
do i=1,n
rho(i,m)=f(0,i,m)+f(1,i,m)+f(3,i,m)+2.*(f(2,i,m)+f(6,i,m)+f(5,i,m))
end do
DO i=1,n
DO j=1,m-1
usum=0.0
vsum=0.0
DO k=0,8
usum=usum+f(k,i,j)*cx(k)
vsum=vsum+f(k,i,j)*cy(k)
END DO
u(i,j)=usum/rho(i,j)
v(i,j)=vsum/rho(i,j)
END DO
END DO
return
end

