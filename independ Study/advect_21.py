# Minghao Liu 

import numpy as geek 
import math as math
import matplotlib.pyplot as plt



# grid resolution
N  = 32
dx = 1./N
xa = 0
xb = 1

# discrete domain
x  = geek.linspace(xa + 0.5 * dx, xb - 0.5 * dx, N)

#IC
u       = [0 for i in range(N+4)]
counter = xa

'''
for i in range(1,N+1):
	if counter <= 0.3:
		u[i] = -1
	elif (counter > 0.3) & (counter <= 0.6):
		u[i] = 2
	elif (counter > 0.6) & (counter <=1) :
		u[i] = -2
	counter = counter + dx
'''
for i in range(2,N+2):
	u[i] = geek.sin(2*math.pi*x[i-2])

#print u


# periodic BC
u[0]   = u[N]
u[1]   = u[N+1]
u[N+2] = u[2]
u[N+3] = u[3]

#outflow BC
#u[0]   = u[2]
#u[1]   = u[3]
#u[N+2] = u[N]
#u[N+3] = u[N+1]

# constant advection velocitypy
c  = 1.

# non-constant advection velocity
Cc = 1.


# CFL
Ca = 0.7


sList = []
for i in range(0,N+3):
	sList.append(abs(0.5 * (u[i+1] + u[i])))


dt = Ca * dx/max(sList)
t  = 0


# Ncycle revolutions for tmax
Ncycle = 1
tmax   = 0.7 #Ncycle*2.*pi/abs(c);

#hold on;
plt.plot(x,u[1:N+1],'r*-')
plt.xlim([0, 1])
plt.grid(color = 'grey', linestyle = '-', linewidth = 0.5)
#plt.show()

 
u0      = []
ul      = []
ur      = []
F_PLM = []

for i in range(0,N+4):
	u0.append(u[i])
	ul.append(0)
	ur.append(0)
	F_PLM.append(0)

uNew = [0 for i in range(N+4)]


# --------- CONTROL CENTER -------------   
#method  0: centered slop.    1: upwind slop.     2:downwind slop
#method  3: minmod            4: MC limiter       5:van Leer s limiter.
#boundry 0: periodic boundry  1: outflow boundry 

method  = 3
boundry = 0
# --------- CONTROL CENTER -------------

# TODO: 2. Not sure why it would not stop				 |  Still have trouble
# TODO: minmod. Def minmod doesn't work.                 |  Solved by hard coding minmod(), in method 3
# TODO: 4. Not finished coding. Not sure about minmod .  |  Solved by hard coding minmod(), in method 4
# TODO: The behaviour of the graph doesn t working ideal.

while t < tmax:
	#Loop for getting ur ul
	for i in range(1,N+3):
		delta = 0

		if method == 0:
			#Delta for the slop, now we are using centered slop
			delta    = 0.5 * (u[i+1]         - u[i-1] ) /dx   #delta i

		elif method == 1:
			delta    = (u[i]   - u[i-1] ) /dx   #delta i

		elif method == 2:  
			delta   = (u[i+1]         - u[i]  ) /dx   #delta i

		elif method == 3:
			a   = (u[i]   - u[i-1]) / dx 
			b    = (u[i+1]         - u[i]  	) / dx			
			delta = 0.5*(geek.sign(a) + geek.sign(b))  * min (abs(a) , abs(b))

		elif method == 4: 
			a   = 0.5 * (u[i+1]           - u[i-1]) / dx 
			b    = 2 * (u[i+1]         - u[i]  	) / dx		
			c   = 2 * (u[i]   - u[i-1]) / dx 

			delta = 0.5*(geek.sign(a) + geek.sign(b))  * min (abs(a) , abs(b))
			delta = 0.5*(geek.sign(a) + geek.sign(b))  * min (abs(delta) , abs(c))
			
		elif method == 5:
			a   = (u[i]   - u[i-1]) / dx 
			b    = (u[i+1]         - u[i]  	) / dx		

			if a*b <= 0:
				delta = 0
			else:
				delta = (2 * a * b) / (a + b)

		ul[i] = u[i] - dx/2 * delta * (1 + Ca)
		ur[i] = u[i] + dx/2 * delta * (1 - Ca)

	#Loop for getting F_PLM
	for i in range(2,N+3):
		
		if u[i-1] >= u[i]:
			s_half = (u[i-1]+u[i])/2
			if s_half >= 0:
				F_PLM[i] = 0.5 * ur[i-1] * ur[i-1]
			elif s_half < 0:
				F_PLM[i] = 0.5 * ul[i] * ul[i]
		else:
			if u[i-1] >= 0:
				F_PLM[i] = 0.5 * ur[i-1] * ur[i-1]
			elif u[i-1] < 0 and u[i] >0:
				F_PLM[i] = 0
			elif u[i] <= 0:
				F_PLM[i] = 0.5 * ul[i] * ul[i]	
			else: 
				F_PLM[i] = 0

	for i in range(2,N+2):
		uNew[i] = u[i] - dt/dx * (F_PLM[i+1] - F_PLM[i])

	#update t
	t = t + dt
   
   
   #update BC
	
	if   boundry == 0:
		uNew[0]   = uNew[N]
		uNew[1]   = uNew[N+1]
		uNew[N+2] = uNew[2]
		uNew[N+3] = uNew[3]  # periodic boundry 
	elif boundry == 1:
		uNew[0]    = uNew[2]
		uuNew[1]   = uNew[3]
		uNew[N+2]  = uNew[N]
		uNew[N+3]  = uNew[N+1] # outflow boundry

   #store updated solution
	u = []
	for i in range(0, N+4):
		u.append(uNew[i])

	sList = []
	for i in range(0, N+3):
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

