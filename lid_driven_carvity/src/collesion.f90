subroutine collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m)
real f(0:8,0:n,0:m)
real feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8), cx(0:8),cy(0:8)
real u(0:n,0:m), v(0:n,0:m)
DO i=0,n
DO j=0,m
t1=u(i,j)*u(i,j)+v(i,j)*v(i,j)
DO k=0,8
t2=u(i,j)*cx(k)+v(i,j)*cy(k)
feq(k,i,j)=rho(i,j)*w(k)*(1.0+3.0*t2+4.50*t2*t2-1.50*t1)
f(k,i,j)=omega*feq(k,i,j)+(1.-omega)*f(k,i,j)
END DO
END DO
END DO
return
end

