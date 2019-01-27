# Minghao Liu 

import numpy as geek 
import math as math
import matplotlib.pyplot as plt

# grid resolution
N=32
dx = 1./N

xa=0
xb=1

# discrete domain
x = geek.linspace(xa+0.5*dx,xb-0.5*dx,N)

#IC
u = [0 for i in range(N+2)]
counter = xa


for i in range(1,N+1):
	if counter <= 0.3:
		u[i] = 2
	elif (counter > 0.3) & (counter <= 0.6):
		u[i] = -1
	elif (counter > 0.6) & (counter <=1) :
		u[i] = 3
	counter = counter + dx


'''
for i in range(1,N+1):
	u[i] = 0.5*math.sin(2*math.pi*x[i-1])+100000 

print u
'''

# periodic BC
#u[0]=u[N];
#u[N+1]=u[1];

#outflow BC
u[0] = u[1]
u[N+1] = u[N]

# constant advection velocitypy
c=1.

# non-constant advection velocity
Cc = 1.


# CFL
Ca = 0.8

sList = []
for i in range(0,N+1):
	sList.append(abs(0.5*(u[i+1]+u[i])))


dt=Ca*dx/max(sList)
t=0;


# Ncycle revolutions for tmax
Ncycle = 1;
tmax=0.3; #Ncycle*2.*pi/abs(c);

#hold on;
plt.plot(x,u[1:N+1],'r*-')
plt.xlim([0, 1])
plt.grid(color='grey', linestyle='-', linewidth=0.5)
#plt.show()

 
u0 = []
ul = []
ur = []
for i in range(0,N+2):
	u0.append(u[i])
	ul.append(u[i])
	ur.append(u[i])


   
uNew = [0 for i in range(N+2)]

while t < tmax:
	for i in range(1,N+1):
		'''
		if u[i] > 0:
			uNew[i] = u[i] - dt/dx*((u[i]*u[i])-(u[i-1]*u[i-1]))/2
		else:
			uNew[i] = u[i] - dt/dx*((u[i+1]*u[i+1])-(u[i]*u[i]))/2   
		
		#u = u0 + au
		#uNew[i] = u[i] - 0.5*c*dt/dx*(u[i+1]-u[i-1]) + 0.5* (c*dt/dx)*(c*dt/dx)*(u[i+1]-2*u[i]+u[i-1])
		#uNew[i] = u[i] - 0.5*c*dt/dx*(3*u[i]-4*u[i-1]+u[i-2]) + 0.5* (c*dt/dx)* (c*dt/dx)*(u[i]-2*u[i-1]+u[i-2])
		#uNew[i] = u[i] - c*dt/dx* (u[i]-u[i-1]) -0.25*c*dt/dx*(1-c*dt/dx)*(u[i+1]-u[i])+0.25*c*dt/dx*(1-c*dt/dx)*(u[i-1]-u[i-2])
      
		#uNew[i] =  u[i] - dt/dx*u[i]*(u[i] - u[i-1])
		#uNew[i] = u[i] - dt/dx*(0.5*u[i]*u[i] - 0.5*u[i-1]*u[i-1])
		'''
		'''
		sr = 0.5 * (u[i+1] + u[i])
		sl = 0.5 * (u[i] + u[i-1])
		if u[i]>=u[i+1]:	
			if sr>0:
				FR = 0.5*u[i]*u[i]
			else:
				FR = 0.5*u[i+1]*u[i+1]

			if sl>0:
				FL = 0.5*u[i-1]*u[i-1]
			else:
				FL = 0.5*u[i]*u[i]

		else:
			if u[i] >= 0:
				FR = 0.5 * u[i]*u[i]
			elif u[i]<0 & 0< u[i+1]:
				FR = 0
			elif u[i+1] <=0:
				FR = 0.5 * u[i+1]*u[i+1]

			if u[i-1] >= 0:
				FL = 0.5 * u[i-1]*u[i-1]
			elif u[i-1]<0 & 0< u[i]:
				FL = 0
			elif u[i] <=0:
				FL = 0.5 * u[i]*u[i]	
		
		uNew[i] = u[i] - dt/dx*(FR-FL)
		'''
		
		#Delta for the slop, now we are using centered slop
		delta = 0.5*(u[i+1] - u[i-1])/dx   #delta i
		delta1 = 0.5*(u[(i+2)%(N+2)] - u[i])/dx    #delta i+1
		delta01 = 0.5*(u[i] - u[i-2])/dx   #delta i-1

		#ul[i] = u[i] - dx/2 * delta * (1 + Ca)
		#ur[i] = u[i] + dx/2 * delta * (1 - Ca)

		#F PLM:
		F_PLM_L = 0  #F plm of i - 0.5
		F_PLM_R = 0	 #F plm of i + 0.5	
		#This is for F_PLM_R
		#flux for shock solution when u[i] > =u [i+1].  
		if u[i] >= u[i+1]:
			s_half = (u[i]+u[i+1])/2
			if s_half >= 0:
				F_PLM_R = 0.5 * (u[i] + 0.5*dx*delta*(1-Ca)) * (u[i] + 0.5*dx*delta*(1-Ca))
			elif s_half < 0:
				F_PLM_R = 0.5 * (u[i+1] - 0.5*dx * delta1*(1+Ca)) * (u[i+1] - 0.5*dx * delta1*(1+Ca))
		#flux for rarefaction solution when u[i] < u [i+1].  
		elif u[i] < u[i+1]:
			if u[i] >= 0:
				F_PLM_R = 0.5 * (u[i] + 0.5*dx*delta*(1-Ca)) * (u[i] + 0.5*dx*delta*(1-Ca))
			elif u[i] < 0 and u[i+1] >0:
				F_PLM_R = 0
			elif u[i+1] <= 0:
				F_PLM_R = 0.5 * (u[i+1] - 0.5*dx * delta1*(1+Ca)) * (u[i+1] - 0.5*dx * delta1*(1+Ca))
		

		#This is for F_PLM_L
		#flux for shock solution when u[i-1] > =u [i].  
		if u[i-1] >= u[i-1]:
			s_half = (u[i-1]+u[i])/2
			if s_half >= 0:
				F_PLM_L = 0.5 * (u[i-1] + 0.5*dx*delta01*(1-Ca)) * (u[i-1] + 0.5*dx*delta01*(1-Ca))
			elif s_half < 0:
				F_PLM_L = 0.5 * (u[i] - 0.5*dx*delta*(1+Ca)) * (u[i] - 0.5*dx*delta*(1+Ca))
		#flux for rarefaction solution when u[i-1] < u [i].  
		elif u[i-1] < u[i]:
			if u[i-1] >= 0:
				F_PLM_L = 0.5 * (u[i-1] + 0.5*dx*delta01*(1-Ca)) * (u[i-1] + 0.5*dx*delta01*(1-Ca))
			elif u[i-1] < 0 and u[i] >0:
				F_PLM_L = 0
			elif u[i] <= 0:
				F_PLM_L = 0.5 * (u[i] - 0.5*dx * delta*(1+Ca)) * (u[i] - 0.5*dx * delta*(1+Ca))
		


		uNew[i] = u[i] - dt/dx * (F_PLM_R - F_PLM_L)

				


      
   #update t
	t = t+ dt
   
   
   #update BC
	uNew[0]=uNew[N]
	uNew[N+1]=uNew[1]
	#uNew[N+1] = uNew[N]

   #store updated solution
	u = []
	for i in range(0,N+2):
		u.append(uNew[i])

	sList = []
	for i in range(0,N+1):
		sList.append(abs(0.5*(u[i+1]+u[i])))

	dt=Ca*dx/abs(max(sList))

   #plot
	plt.clf()
	plt.plot(x,u[1:N+1],'r*-')
	plt.xlim([0, 1])
	plt.grid(True)
	plt.pause(0.1)


plt.plot(x,u0[1:N+1],'b*-')
plt.show()

   
   
   
   
#