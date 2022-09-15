module global
implicit none

integer,parameter :: N=91
integer,parameter :: NX=45,NY=45
integer,parameter :: P=3 ! number of DOFs within a 1D element

real,parameter :: aaaa=6371229.,a=6371229.
real,parameter :: pi=3.14159265358979
!Output
real,parameter :: plevel=70.e3

! Grid

integer ids,ide,jds,jde
integer ims,ime,jms,jme

real xl,xr,yl,yr
real x(-4+1:2*NX+1+4),y(-4+1:2*NY+1+4),dx,dy


!Boundaries

integer bci(10,6*(8*(2*NX+1)+8*(2*NY+1)))
real bcr(6*(8*(2*NX+1)+8*(2*NY+1)))


end module global


program cubedsphere
use global
implicit none

integer i,j

ids = 1;     ide = 2*NX+1; jds = 1;     jde = 2*NY+1
ims = ids-4; ime = ide+4;  jms = jds-4; jme = jde+4

xl=-pi/4.*a
yl=-pi/4.*a

xr=pi/4.*a
yr=pi/4.*a


dx=(xr-xl)/float(NX)
dy=(yr-yl)/float(NY)


do i=-4+1,2*NX+1+4
x(i)=xl+float(i-1)*dx/2.
end do

do j=-4+1,2*NY+1+4
y(j)=yl+float(j-1)*dy/2.
end do

call ghostlocation()

stop
end program



subroutine pprop2sp(lambda,theta,x,y,k)
use global, only : aaaa
implicit none
real lambda,theta,x,y,x1,y1,a,b,pi,r
integer k

pi=2.*asin(1.)
r=aaaa

x1=x/r; y1=y/r
	
select case(k)	      
case(1:4)
lambda=x1+float(k-1)*pi/2.; theta=atan2(tan(y1)*cos(x1),1.)
case(5)
a=tan(x1); b=tan(y1)
lambda=atan2(a,-b); theta=atan2(1.,sqrt(a*a+b*b))
case(6)
a=tan(x1); b=tan(y1)
lambda=atan2(a,b); theta=-atan2(1.,sqrt(a*a+b*b))  
end select
	
return
end

subroutine pprosp2p(x,y,lambda,theta,k)
use global, only : aaaa
implicit none
real lambda,theta,x,y,x1,y1,a,b,pi,r
integer k

pi=2.*asin(1.)
r=aaaa

select case(k)
case(1:4)
x=atan(tan(lambda-float(k-1)*pi/2.))
y=atan(tan(theta)/cos(lambda-float(k-1)*pi/2.))
case(5:6)
x=atan((-1.)**(k+1)*sin(lambda)/tan(theta))
y=atan(-cos(lambda)/tan(theta))
end select

x=x*r; y=y*r
	
return
end

subroutine contravprosp2p(contrav1,contrav2,sv1,sv2,k,lambda,theta)
implicit none
real*8 contrav1,contrav2,sv1,sv2,lambda,theta,ia(2,2)
integer k

call matrixIA(ia,k,lambda,theta)
      
contrav1=ia(1,1)*sv1+ia(1,2)*sv2
contrav2=ia(2,1)*sv1+ia(2,2)*sv2

return
end

subroutine contravprop2sp(sv1,sv2,contrav1,contrav2,k,lambda,theta)
implicit none
real*8 contrav1,contrav2,sv1,sv2,lambda,theta,a(2,2)
integer k
      
call matrixA(a,k,lambda,theta)      
sv1=a(1,1)*contrav1+a(1,2)*contrav2
sv2=a(2,1)*contrav1+a(2,2)*contrav2	

return
end

subroutine matrixIA(ima,k,lambda,theta)
use global, only : aaaa
implicit none
real ima(2,2),lambda,theta,alambda,atheta,pi,a,b,c,d,temp,r
integer k     
      
pi=2.*asin(1.)
r=aaaa

if (k <= 4) then
alambda=lambda-float(k-1)*pi/2.; atheta=theta
a=sin(alambda); b=cos(alambda); c=sin(atheta); d=cos(atheta); temp=d*d*b*b+c*c
ima(1,1)=1./d; ima(1,2)=0.
ima(2,1)=a*c/temp; ima(2,2)=b/temp
else if (k==5) then
alambda=lambda; atheta=theta
a=sin(alambda); b=cos(alambda); c=sin(atheta); d=cos(atheta); temp=c+a*a*d*d/c
ima(1,1)=b/temp; ima(1,2)=-a/c/temp
temp=c+b*b*d*d/c; ima(2,1)=a/temp; ima(2,2)=b/c/temp
else
alambda=lambda; atheta=theta
a=sin(alambda); b=cos(alambda); c=sin(atheta); d=cos(atheta); temp=c+a*a*d*d/c
ima(1,1)=-b/temp; ima(1,2)=a/c/temp
temp=c+b*b*d*d/c; ima(2,1)=a/temp; ima(2,2)=b/c/temp
endif	

!ima=ima/r/r
	      
return
end

subroutine matrixA(ma,k,lambda,theta)
use global, only : aaaa
implicit none
real ma(2,2),lambda,theta,alambda,atheta,pi,a,b,c,d,temp,r
integer k
      
pi=2.*asin(1.)
r=aaaa

if (k <= 4) then
alambda=lambda-float(k-1)*pi/2.; atheta=theta
a=sin(alambda); b=cos(alambda); c=sin(atheta); d=cos(atheta)
ma(1,1)=d; ma(1,2)=0.d0; ma(2,1)=-c*d*a/b; ma(2,2)=b*d*d+c*c/b
else if (k==5) then
alambda=lambda; atheta=theta
a=sin(alambda); b=cos(alambda); c=sin(atheta); d=cos(atheta); temp=1.+a*a*d*d/c/c
ma(1,1)=b*c*temp; ma(2,1)=-c*c*a*temp; temp=1.+b*b*d*d/c/c
ma(1,2)=a*c*temp; ma(2,2)=b*c*c*temp
else
alambda=lambda; atheta=theta
a=sin(alambda); b=cos(alambda); c=sin(atheta); d=cos(atheta); temp=1.+a*a*d*d/c/c
ma(1,1)=-b*c*temp; ma(2,1)=c*c*a*temp; temp=1.+b*b*d*d/c/c
ma(1,2)=a*c*temp; ma(2,2)=b*c*c*temp
endif

!ma=ma*r*r

return
end

subroutine matrixG(mg,xx,yy)
use global, only : a
implicit none
real mg(2,2),x,y,rho,xx,yy
 
x=xx/a; y=yy/a
      
rho=dsqrt(1.+tan(x)**2.+tan(y)**2.)

mg(1,1)=(1.+tan(x)**2.)/(rho**4.*cos(x)**2.*cos(y)**2.)
mg(1,2)=-tan(x)*tan(y)/(rho**4.*cos(x)**2.*cos(y)**2.)
mg(2,1)=mg(1,2)
mg(2,2)=(1.+tan(y)**2.)/(rho**4.*cos(x)**2.*cos(y)**2.)

!mg=a*a*mg

return
end

subroutine matrixIG(img,xx,yy)
use global, only : a
implicit none
real img(2,2),x,y,rho,xx,yy

x=xx/a; y=yy/a

rho=sqrt(1.+tan(x)**2.+tan(y)**2.)

img(1,1)=(1.+tan(y)**2.)*(rho**2.*cos(x)**2.*cos(y)**2.)
img(1,2)=tan(x)*tan(y)*(rho**2.*cos(x)**2.*cos(y)**2.)
img(2,1)=img(1,2)
img(2,2)=(1.+tan(x)**2.)*(rho**2.*cos(x)**2.*cos(y)**2.)

!img=img/a/a

return
end

subroutine computejab(jab,xx,yy)
use global, only : a
implicit none
real jab,x,y,rho,xx,yy

x=xx/a; y=yy/a

rho=sqrt(1+tan(x)*tan(x)+tan(y)*tan(y))
jab=(1+tan(x)*tan(x))*(1+tan(y)*tan(y))/rho/rho/rho

!jab=a*a*(1+tan(x)*tan(x))*(1+tan(y)*tan(y))/rho/rho/rho

return
end subroutine computejab

subroutine ghostlocation()
use global, only: Nx,Ny,P,bci,bcr,xl,yl,x,y,dx,dy,a,ids,ide,jds,jde,ims,ime,jms,jme
implicit none
integer i,j,k,no,nos,noe
real tmp,dd

! ghost cell location

no=1

do k=1,6

nos=(k-1)*(8*(2*NX+1)+8*(2*NY+1))+1; noe=k*(8*(2*NX+1)+8*(2*NY+1))

bci(1,nos:noe)=k

! left boundary

do j=jds,jde
do i=ims,ids-1
  bci(2,no)=i; bci(3,no)=j
  select case(k)
    case(1)
      bci(6,no)=4; bci(7,no)=2*NX+i; bci(9,no)=2; tmp= 1.
    case(2)
      bci(6,no)=1; bci(7,no)=2*NX+i; bci(9,no)=2; tmp= 1.
    case(3)
      bci(6,no)=2; bci(7,no)=2*NX+i; bci(9,no)=2; tmp= 1.
    case(4)
      bci(6,no)=3; bci(7,no)=2*NX+i; bci(9,no)=2; tmp= 1.
    case(5)
      bci(6,no)=4; bci(8,no)=2*NY+i; bci(9,no)=1; tmp=-1.
    case(6)
      bci(6,no)=4; bci(8,no)=2-i;    bci(9,no)=1; tmp= 1.
  end select
  dd=a*tmp*atan(tan(y(j)/a)/abs(tan(x(i)/a)))
  if (bci(9,no)==1) then
    bci(7,no)=int((dd-xl)/dx)+1; bci(7,no)=max(min(bci(7,no),Nx),1)
    bcr(no)=dd !-x(bci(7,no),1)
  else
    bci(8,no)=int((dd-yl)/dy)+1; bci(8,no)=max(min(bci(8,no),Ny),1)
    bcr(no)=dd !-y(bci(8,no),1)
  end if
  no=no+1
end do
end do

! bottom boundary

do i=ids,ide
do j=jms,jds-1
  bci(2,no)=i; bci(3,no)=j
  select case(k)
    case(1)
      bci(6,no)=6; bci(8,no)=2*NY+j; bci(9,no)=1; tmp= 1.
    case(2)
      bci(6,no)=6; bci(7,no)=2*NX+j; bci(9,no)=2; tmp=-1.
    case(3)
      bci(6,no)=6; bci(8,no)=2-j;    bci(9,no)=1; tmp=-1.
    case(4)
      bci(6,no)=6; bci(7,no)=2-j;    bci(9,no)=2; tmp= 1.
    case(5)
      bci(6,no)=1; bci(8,no)=2*NY+j; bci(9,no)=1; tmp= 1.
    case(6)
      bci(6,no)=3; bci(8,no)=2-j;    bci(9,no)=1; tmp=-1.
  end select
  dd=a*tmp*atan(tan(x(i)/a)/abs(tan(y(j)/a)))
  if (bci(9,no)==1) then
    bci(7,no)=int((dd-xl)/dx)+1; bci(7,no)=max(min(bci(7,no),Nx),1)
    bcr(no)=dd !-x(bci(7,no),1)
  else
    bci(8,no)=int((dd-yl)/dy)+1; bci(8,no)=max(min(bci(8,no),Ny),1)
    bcr(no)=dd !-y(bci(8,no),1)
  end if
  no=no+1
end do
end do

! right boundary

do j=jds,jde
do i=ide+1,ime

  bci(2,no)=i; bci(3,no)=j
  select case(k)
    case(1)
      bci(6,no)=2; bci(7,no)=i-2*NX;              bci(9,no)=2; tmp= 1.
    case(2)
      bci(6,no)=3; bci(7,no)=i-2*NX;              bci(9,no)=2; tmp= 1.
    case(3)
      bci(6,no)=4; bci(7,no)=i-2*NX;              bci(9,no)=2; tmp= 1.
    case(4)
      bci(6,no)=1; bci(7,no)=i-2*NX;              bci(9,no)=2; tmp= 1.
    case(5)
      bci(6,no)=2; bci(8,no)=2*NY+1-(i-(2*NX+1)); bci(9,no)=1; tmp= 1.
    case(6)
      bci(6,no)=2; bci(8,no)=1+(i-(2*NX+1));      bci(9,no)=1; tmp=-1.
  end select
  dd=a*tmp*atan(tan(y(j)/a)/abs(tan(x(i)/a)))
  if (bci(9,no)==1) then
    bci(7,no)=int((dd-xl)/dx)+1; bci(7,no)=max(min(bci(7,no),Nx),1)
    bcr(no)=dd !-x(bci(7,no),1)
  else
    bci(8,no)=int((dd-yl)/dy)+1; bci(8,no)=max(min(bci(8,no),Ny),1)
    bcr(no)=dd !-y(bci(8,no),1)
  end if
  no=no+1
end do
end do

! top boundary

do i=ids,ide
do j=jde+1,jme
  bci(2,no)=i; bci(3,no)=j
  select case(k)
    case(1)
      bci(6,no)=5; bci(8,no)=1+(j-(2*NY+1));      bci(9,no)=1; tmp= 1.
    case(2)
      bci(6,no)=5; bci(7,no)=2*NX+1-(j-(2*NY+1)); bci(9,no)=2; tmp= 1.
    case(3)
      bci(6,no)=5; bci(8,no)=2*NY+1-(j-(2*NY+1)); bci(9,no)=1; tmp=-1.
    case(4)
      bci(6,no)=5; bci(7,no)=1+(j-(2*NY+1));      bci(9,no)=2; tmp=-1.
    case(5)
      bci(6,no)=3; bci(8,no)=2*NY+1-(j-(2*NY+1)); bci(9,no)=1; tmp=-1.
    case(6)
      bci(6,no)=1; bci(8,no)=1+(j-(2*NY+1));      bci(9,no)=1; tmp= 1.
  end select
  dd=a*tmp*atan(tan(x(i)/a)/abs(tan(y(j)/a)))
  if (bci(9,no)==1) then
    bci(7,no)=int((dd-xl)/dx)+1; bci(7,no)=max(min(bci(7,no),Nx),1)
    bcr(no)=dd !-x(bci(7,no),1)
  else
    bci(8,no)=int((dd-yl)/dy)+1; bci(8,no)=max(min(bci(8,no),Ny),1)
    bcr(no)=dd !-y(bci(8,no),1)
  end if
  no=no+1
end do
end do


end do


end subroutine ghostlocation
