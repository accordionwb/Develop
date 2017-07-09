subroutine result(u,v,rho,uo,n,m)
real u(0:n,0:m),v(0:n,0:m)
real rho(0:n,0:m),strf(0:n,0:m)
open(50, file='streamf')
! streamfunction calculations
strf(0,0)=0.
do i=0,n
rhoav=0.5*(rho(i-1,0)+rho(i,0))
if(i.ne.0) strf(i,0)=strf(i-1,0)-rhoav*0.5*(v(i-1,0)+v(i,0))
do j=1,m
rhom=0.5*(rho(i,j)+rho(i,j-1))
strf(i,j)=strf(i,j-1)+rhom*0.5*(u(i,j-1)+u(i,j))
end do
end do
! ———————————–
write(20,*) "VARIABLES =X, Y, U, V, S"
write(20,*) "ZONE ","I=",n+1,"J=",m+1,",","F=BLOCK"
do j=0,m
write(20,*)(i,i=0,n)
end do
do j=0,m
write(20,*)(j,i=0,n)
end do
do j=0,m
write(20,*)(u(i,j),i=0,n)
end do
do j=0,m
write(20,*)(v(i,j),i=0,n)
end do
do j=0,m
write(20,*)(strf(i,j),i=0,n)
end do
do j=0,m
write(30,*)j/float(m),u(n/2,j)/uo,u(n/4,j)/uo,u(3*n/4,j)/uo
end do
do i=0,n
write(40,*) i/float(n),v(i,m/2)/uo
end do
return
end

