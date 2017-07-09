! Computer code for lid-driven carvity flow
Program lid_cf
parameter (n=100,m=100)
real f(0:8,0:n,0:m)
real feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8), cx(0:8),cy(0:8)
real u(0:n,0:m), v(0:n,0:m)
integer i, n_p_sec, ia, ie
real :: t
character(2) :: percent = '%' //char(13)
open(unit=20,file='uvfield')
open(unit=30,file='uvely')
open(unit=40,file='vvelx')
open(unit=80,file='timeu')
!
write(*,*) 'Start Program:'
!
call system_clock(count_rate=n_p_sec)
call system_clock(count=ia)
!

uo=0.10
sumvelo=0.0
rhoo=5.00
dx=1.0
dy=dx
dt=1.0
alpha=0.01
Re=uo*m/alpha
write(*, *) "Re=", Re
omega=1.0/(3.*alpha+0.5)
mstep=40000
w(0)=4./9.
do i=1,4
w(i)=1./9.
end do
do i=5,8
w(i)=1./36.
end do
cx(0)=0
cx(1)=1
cx(2)=0
cx(3)=-1
cx(4)=0
cx(5)=1
cx(6)=-1
cx(7)=-1
cx(8)=1
cy(0)=0
cy(1)=0
cy(2)=1
cy(3)=0
cy(4)=-1
cy(5)=1
cy(6)=1
cy(7)=-1
cy(8)=-1
do j=0,m
do i=0,n
rho(i,j)=rhoo
u(i,j)=0.0
v(i,j)=0.0
end do
end do
do i=1,n-1
u(i,m)=uo
v(i,m)=0.0
end do

! ####### main loop #################
1 do kk=1,mstep
call collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m)
call streaming(f,n,m)
! ——————————–
call sfbound(f,n,m,uo)
call rhouv(f,rho,u,v,cx,cy,n,m)
!write(*,*) u(0,m/2),v(0,m/2),rho(0,m/2),u(n,m/2),v(n,m/2),rho(n,m/2)
write(80,*) kk,u(n/2,m/2),v(n/2,m/2)
!! time count
call sleepqq(10)
write(*,'(i2,a\)') kk,percent

END DO
! end of the main loop

call result(u,v,rho,uo,n,m)
call system_clock(count=ie)
t = (ie-ia)/real(n_p_sec)
write(*,*) "Time in second: ",t
stop
end
! end of the main program


