subroutine sfbound(f,n,m,uo)
real f(0:8,0:n,0:m)
do j=0,m
! bounce back on west boundary
f(1,0,j)=f(3,0,j)
f(5,0,j)=f(7,0,j)
f(8,0,j)=f(6,0,j)
! bounce back on east boundary
f(3,n,j)=f(1,n,j)
f(7,n,j)=f(5,n,j)
f(6,n,j)=f(8,n,j)
end do
! bounce back on south boundary
do i=0,n
f(2,i,0)=f(4,i,0)
f(5,i,0)=f(7,i,0)
f(6,i,0)=f(8,i,0)
end do
! moving lid, north boundary
do i=1,n-1
rhon=f(0,i,m)+f(1,i,m)+f(3,i,m)+2.*(f(2,i,m)+f(6,i,m)+f(5,i,m))
f(4,i,m)=f(2,i,m)
f(8,i,m)=f(6,i,m)+rhon*uo/6.0
f(7,i,m)=f(5,i,m)-rhon*uo/6.0
end do
return
end

