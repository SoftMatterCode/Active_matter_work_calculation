!****************************************************************************************
! 
!          active particles in cirular box with time dependent steepess
!          
!             last modifed :: 20th november 2021
!****************************************************************************************
program test
implicit none
real*8, dimension (1:5000) ::x,y,vx,vy,dumx,dumy,direcx,direcy,vxdum,vydum,dumvx,dumvy
real*8, dimension (1:1000) ::const1,const2,r2,dw,w,w_config
real*8:: eta,L,l0,v0,t,dt,tau,r,pi,x2,y2,a,bound,alpha,fric,pot1,pot2,pot3
real*8:: eff,a0,f3x,f3y,direcabs,stef,hig,KBT,gasdev,ran1
real*8:: dstif1,dstif2,stif,coshb,w_sum,w_cycle,r_b,KBT1,KBT2
integer:: i,j,N,idum,time,k,interval,halftime,nconfig,z

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stef=2.50d0  ! stiffness
v0=20.0d0     ! activity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

N=1000                ! number of particles
tau=16000!8000        ! time
idum=967369025
alpha=1.0d0
fric=10.0d0           ! gamma (friction)
KBT=0.00010d0         ! Temperature
KBT1=0.0010d0
KBT2=0.00010d0 
a=1.0d0               ! Interparticle interaction potential parameter 
a0=a*(2.0**(1.0/6))   ! cut off for interparticle interaction (WCA)  
eff=1.0d0             ! Interparticle interaction potential parameter
hig=100.0d0           ! Trap height          
bound=80.0d0          ! Trap radius  
eta=1.0d0             ! Angular noise span 
L=bound-2.0d0         ! Initial condition parameter, so that particles are inside
l0=2.0d0               ! Vicsek radius
dt=0.001d0            ! Time increment
t=0.0d0               ! Initial Time
pi=4*atan(1.0d0)       ! pi value
interval=2000!200!2000
halftime=(int(tau/dt))/2
time = int(tau/dt)

nconfig=1

const1=0.0d0
const2=0.0d0
r2=0.0d0
w_config=0.0d0

dstif2 = (2.40d0/(halftime + 1))
dstif1 = (-2.40d0/halftime)
!w_cycle=0.0d0
!w_sum=0.0d0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 open(unit=1,file='input.dat',form='formatted')
 !open(unit=2,file='output.dat',form='formatted')
 !open(unit=3,file='work.dat',form='formatted')
 open(unit=4,file='total_work.dat',form='formatted')
   do i=1,N
       read(unit=1,fmt=*)x(i),y(i),vx(i),vy(i)
   end do

  ! do i=1,N
   !    write(unit=2,fmt=*)x(i),y(i),vx(i),vy(i)
   !end do

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do z = 1 , nconfig

  w_cycle=0.0d0
  w_sum=0.0d0

  do k = 1,time 
     !write(*,*)k 
     w=0.0d0
     dw=0.0d0
      write(*,*)k
     !$OMP PARALLEL
     !$OMP DO PRIVATE(i)
     do i=1,N
       direcx(i)=vx(i)
       direcy(i)=vy(i)
       f3x=0.0d0
       f3y=0.0d0
       r2=0.0d0
       
      
       do j=1,N
          if(j.ne.i)then
            x2=x(i)-x(j)
            y2=y(i)-y(j)

            r=sqrt(x2**2  + y2**2)
            pot1=0.0d0
            pot2=0.0d0
            pot3=0.0d0
            stif=0.0d0
            coshb=0.0d0
            r_b=0.0d0
! Direction choose
                    
            if(r.lt.l0)then
               direcx(i)=direcx(i)+vx(j)
               direcy(i)=direcy(i)+vy(j)             
            endif 
            
            
! WCA interaction

            
            if(r.lt.a0)then
               pot3=4*eff*(12*(a/r)**12-6*(a/r)**6)/r
               f3x = f3x + pot3*(x2/r)  
               f3y = f3y + pot3*(y2/r)
            endif

          endif
        end do
        
        direcabs =sqrt(direcx(i)**2 + direcy(i)**2)
        direcx(i)=direcx(i)/direcabs
        direcy(i)=direcy(i)/direcabs

        r2(i)=ran1(idum) 
        r2(i)=r2(i)-0.5d0
        r2(i)=r2(i)*eta

        vxdum(i)=(direcx(i)*cos(r2(i))-direcy(i)*sin(r2(i)))
        vydum(i)=(direcx(i)*sin(r2(i))+direcy(i)*cos(r2(i)))    

        r=sqrt(x(i)**2 +y(i)**2)

        const1(i)=gasdev(idum)
        const2(i)=gasdev(idum)

       if(k.le. halftime) then

        r_b = r-bound
        stif=((-2.40d0*k)/halftime + 2.50d0)
        coshb=((cosh((stif)*(r_b)))**2)
        pot1=-hig*(stif/coshb)        
        
        dumvx(i)=-fric*vx(i)*dt+v0*vxdum(i)*dt+alpha*(pot1*x(i)/r+f3x)*dt+sqrt(2.0d0*fric*KBT1*dt)*const1(i)
        dumvy(i)=-fric*vy(i)*dt+v0*vydum(i)*dt+alpha*(pot1*y(i)/r+f3y)*dt+sqrt(2.0d0*fric*KBT1*dt)*const2(i)
        
        dw(i) = (hig*dstif1*dt*(r_b))/coshb  
        
        else

        r_b = r-bound
        stif=(2.40d0*(k-halftime+1)/(halftime+1) + 0.10d0)
        coshb=((cosh((stif)*(r_b)))**2)
        pot2=-hig*(stif/coshb)
        
        dumvx(i)=-fric*vx(i)*dt+v0*vxdum(i)*dt+alpha*(pot2*x(i)/r+f3x)*dt+sqrt(2.0d0*fric*KBT2*dt)*const1(i)
        dumvy(i)=-fric*vy(i)*dt+v0*vydum(i)*dt+alpha*(pot2*y(i)/r+f3y)*dt+sqrt(2.0d0*fric*KBT2*dt)*const2(i)
        
        dw(i) = (hig*dstif2*dt*(r_b))/coshb
        
        end if

        dumx(i)=vx(i)*dt
        dumy(i)=vy(i)*dt
!             
      end do

      !$OMP END DO
      !$OMP END PARALLEL
      
     
     do i= 1,N
        
        w(i) = w(i) + dw(i)
        
        vx(i)=vx(i)+dumvx(i)
        vy(i)=vy(i)+dumvy(i)

        x(i)=x(i)+dumx(i)
        y(i)=y(i)+dumy(i)         

      end do

      w_sum = sum(w)/N
      w_cycle = w_cycle +  w_sum

      !if (mod(k,interval).eq. 0.0d0)then
      !   write(*,*)k 
      !   write(unit=3,fmt=*)k,w_sum      
      !   do i=1,N
      !        write(unit=2,fmt=*)x(i),y(i),vx(i),vy(i)
      !   end do
      !end if   

   
  end do !~~~~~~~~~~~~~ end of time (k) loop 

!   w_config(z) = w_cycle
   write(unit=4,fmt=*)z,w_cycle
   
end do   !~~~~~~~~~~~~~~~~~~~~~~~` end of config (z) loop  

end program test

	!********************************************
	!Gaussian Random Number generation
	!********************************************

	FUNCTION gasdev(idum)
	INTEGER idum
	REAL*8 gasdev
!      	USES ran1
	INTEGER iset

	real*8 fac,gset,rsq,v1,v2,ran1
	SAVE iset,gset
	DATA iset/0/
	if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
	v2=2.*ran1(idum)-1.
	rsq=v1**2+v2**2
	if(rsq.ge.1..or.rsq.eq.0.)goto 1
	fac=sqrt(-2.*log(rsq)/rsq)
	gset=v1*fac
	gasdev=v2*fac
	iset=1
	else
	gasdev=gset
	iset=0
	endif
	return
	END

	FUNCTION ran1(idum)
	INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
	real*8 ran1,AM,EPS,RNMX
	PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836)
	PARAMETER (NTAB=32)
	PARAMETER (NDIV=1+(IM-1)/NTAB)
	PARAMETER (EPS=1.2e-7)
	PARAMETER (RNMX=1.-EPS)

	INTEGER j,k,iv(NTAB),iy
	SAVE iv,iy
	DATA iv /NTAB*0/, iy /0/
	if (idum.le.0.or.iy.eq.0) then
	idum=max(-idum,1)
	do 11 j=NTAB+8,1,-1
	k=idum/IQ
	idum=IA*(idum-k*IQ)-IR*k
	if (idum.lt.0) idum=idum+IM
	if (j.le.NTAB) iv(j)=idum
11      continue
	iy=iv(1)
	endif
	k=idum/IQ
	idum=IA*(idum-k*IQ)-IR*k
	if (idum.lt.0) idum=idum+IM
	j=1+iy/NDIV
	iy=iv(j)
	iv(j)=idum
	ran1=min(AM*iy,RNMX)
	return
	END
