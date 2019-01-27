# Minghao Liu 

import numpy as geek 
import math as math
import matplotlib.pyplot as plt


  
def minmod (x, y):
	return y
	
	if   abs(x) < abs(b) and x*y > 0:
		return x
	elif abs(y) < abs(x) and x*y > 0:
		return y
	elif x * y <= 0 :
		return 0
	


# grid resolution
N  = 32
dx = 1./N
xa = 0
xb = 1

# discrete domain
x  = geek.linspace(xa + 0.5 * dx, xb - 0.5 * dx, N)

#IC
u       = [0 for i in range(N+2)]
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
u[0]   = u[1]
u[N+1] = u[N]

# constant advection velocitypy
c  = 1.

# non-constant advection velocity
Cc = 1.


# CFL
Ca = 0.7


sList = []
for i in range(0,N+1):
	sList.append(abs(0.5 * (u[i+1] + u[i])))


dt = Ca * dx/max(sList)
t  = 0


# Ncycle revolutions for tmax
Ncycle = 1
tmax   = 0.3 #Ncycle*2.*pi/abs(c);

#hold on;
plt.plot(x,u[1:N+1],'r*-')
plt.xlim([0, 1])
plt.grid(color = 'grey', linestyle = '-', linewidth = 0.5)
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
		
		#method 0: centered slop.  1: upwind slop.  2:downwind slop
		#method 3: minmod          4: MC limiter    5:van Leer s limiter.

		# TODO: 2. Not sure why it would not stop
		# TODO: minmod. Def minmod doesn't work
		# TODO: 4. Not finished coding. Not sure about minmod 
		method = 3

		if method == 0:
			#Delta for the slop, now we are using centered slop
			delta    = 0.5 * (u[i+1]         - u[i-1] ) /dx   #delta i
			delta1   = 0.5 * (u[(i+2)%(N+2)] - u[i]   ) /dx    #delta i+1
			delta01  = 0.5 * (u[i]           - u[i-2] ) /dx   #delta i-1
		elif method == 1:
			delta    = (u[i]   - u[i-1] ) /dx   #delta i
			delta1   = (u[i+1] - u[i]   ) /dx   #delta i+1
			delta01  = (u[i-1] - u[i-2] ) /dx   #delta i-1
		elif method == 2:  #bug
			delta   = (u[i+1]         - u[i]  ) /dx   #delta i
			delta1  = (u[(i+2)%(N+2)] - u[i+1]) /dx   #delta i+1
			delta01 = (u[i]           - u[i-1]) /dx   #delta i-1
		elif method == 3:
			a   = (u[i]   - u[i-1]) / dx 
			a1  = (u[i+1] - u[i]  ) / dx 
			a01 = (u[i]   - u[i-1]) / dx 

			b    = (u[i+1]         - u[i]  	) / dx			
			b1   = (u[(i+2)%(N+2)] - u[i+1] ) / dx
			b01  = (u[i]           - u[i-1] ) / dx	

			
			if   abs(a) < abs(b) and a*b > 0:
				delta = a
			elif abs(b) < abs(a) and a*b > 0:
				delta = b
			elif a * b <= 0 :
				delta = 0

			if   abs(a1) < abs(b1) and a1*b1 > 0:
				delta1 = a1
			elif abs(b1) < abs(a1) and a1*b1 > 0:
				delta1 = b1
			elif a1*b1 <= 0 :
				delta1 = 0

			if   abs(a01) < abs(b01) and a01*b01 > 0:
				delta01 = a01
			elif abs(b01) < abs(a01) and a01*b01 > 0:
				delta01 = b01
			elif a01*b01 <= 0 :
				delta01 = 0
			
			#delta    = minmod (a  , b  )
			#delta1   = minmod (a1 , b1 )
			#delta01  = minmod (a01, b01)

		elif method == 4: 
			a   = 0.5 * (u[i+1]           - u[i-1]) / dx 
			a1  = 0.5 * (u[(i+2)%(N+2)]   - u[i]  ) / dx 
			a01 = 0.5 * (u[i]             - u[i-2]) / dx 

			b    = 2 * (u[i+1]         - u[i]  	) / dx			
			b1   = 2 * (u[(i+2)%(N+2)] - u[i+1] ) / dx
			b01  = 2 * (u[i]           - u[i-1] ) / dx	

			c   = 2 * (u[i]   - u[i-1]) / dx 
			c1  = 2 * (u[i+1] - u[i]  ) / dx 
			c01 = 2 * (u[i]   - u[i-1]) / dx 

			delta   = minmod(minmod (a, b), c)
			delta1  = minmod(minmod (a1, b1), c1)
			delta01 = minmod(minmod (a01, b01), c01)

		elif method == 5:
			a   = (u[i]   - u[i-1]) / dx 
			a1  = (u[i+1] - u[i]  ) / dx 
			a01 = (u[i]   - u[i-1]) / dx 

			b    = (u[i+1]         - u[i]  	) / dx			
			b1   = (u[(i+2)%(N+2)] - u[i+1] ) / dx
			b01  = (u[i]           - u[i-1] ) / dx	

			if a*b <= 0:
				delta = 0
			else:
				delta = (2 * a * b) / (a + b)
			
			if a1*b1 <= 0:
				delta1 = 0
			else:
				delta1 = (2 * a1 * b1) / (a1 + b1)

			if a01*b01 <= 0:
				delta01 = 0
			else:
				delta01 = (2 * a01 * b01) / (a01 + b01)

			#delta   = 1
			#delta1  = 1
			#delta01 = 1

		#ul[i] = u[i] - dx/2 * delta * (1 + Ca)
		#ur[i] = u[i] + dx/2 * delta * (1 - Ca)

		#F PLM:
		F_PLM_L = 0  #F plm of i - 0.5
		F_PLM_R = 0	 #F plm of i + 0.5	
		#This is for F_PLM_R
		#flux for shock solution when u[i] > =u [i+1].  
		if u[i] >= u[i+1]:
			s_half = (u[i] + u[i+1])/2
			if s_half >= 0:
				F_PLM_R = 0.5 * (u[i]   + 0.5 * dx * delta  * (1-Ca)) * (u[i]   + 0.5 * dx * delta  * (1 - Ca))
			elif s_half < 0:
				F_PLM_R = 0.5 * (u[i+1] - 0.5 * dx * delta1 * (1+Ca)) * (u[i+1] - 0.5 * dx * delta1 * (1 + Ca))
		#flux for rarefaction solution when u[i] < u [i+1].  
		elif u[i] < u[i+1]:
			if u[i] >= 0:
				F_PLM_R = 0.5 * (u[i]   + 0.5 * dx * delta  * (1 - Ca)) * (u[i]   + 0.5 * dx * delta  * (1 - Ca))
			elif u[i] < 0 and u[i+1] >0:
				F_PLM_R = 0
			elif u[i+1] <= 0:
				F_PLM_R = 0.5 * (u[i+1] - 0.5 * dx * delta1 * (1 + Ca)) * (u[i+1] - 0.5 * dx * delta1 * (1 + Ca))
		

		#This is for F_PLM_L
		#flux for shock solution when u[i-1] > =u [i].  
		if u[i-1] >= u[i-1]:
			s_half = (u[i-1] + u[i])/2
			if s_half  >= 0:
				F_PLM_L = 0.5 * (u[i-1] + 0.5 * dx * delta01 * (1 - Ca)) * (u[i-1] + 0.5 * dx * delta01 * (1 - Ca))
			elif s_half < 0:
				F_PLM_L = 0.5 * (u[i]   - 0.5 * dx * delta   * (1 + Ca)) * (u[i]   - 0.5 * dx * delta   * (1 + Ca))
		#flux for rarefaction solution when u[i-1] < u [i].  
		elif u[i-1] < u[i]:
			if u[i-1] >= 0:
				F_PLM_L = 0.5 * (u[i-1] + 0.5 * dx * delta01 * (1 - Ca)) * (u[i-1] + 0.5 * dx * delta01 * (1 - Ca))
			elif u[i-1] < 0 and u[i] >0:
				F_PLM_L = 0
			elif u[i] <= 0:
				F_PLM_L = 0.5 * (u[i]   - 0.5 * dx * delta   * (1 + Ca)) * (u[i]   - 0.5 * dx * delta   * (1 + Ca))
		

		# Calculation of u : time at n+1 space at i
		uNew[i] = u[i] - dt/dx * (F_PLM_R - F_PLM_L)

				


      
   #update t
	t = t + dt
   
   
   #update BC
	uNew[0]   = uNew[N]
	#uNew[N+1] = uNew[1]
	uNew[N+1] = uNew[N]

   #store updated solution
	u = []
	for i in range(0, N+2):
		u.append(uNew[i])

	sList = []
	for i in range(0, N+1):
		sList.append(abs(0.5 * (u[i+1] + u[i])))

	dt=Ca * dx/abs(max(sList))

   #plot
	plt.clf   ()
	plt.plot  (x,u[1:N+1], 'r*-')
	plt.xlim  ([0, 1])
	plt.grid  (True)
	plt.pause (0.1)


plt.plot(x, u0[1:N+1], 'b*-')
plt.show()

   
   
   
 