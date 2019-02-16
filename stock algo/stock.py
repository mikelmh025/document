# Minghao Liu 
# Buy and sell without algorithm
# Only buy first
# Only one transaction at a time


import numpy as geek 
import math as math
import matplotlib.pyplot as plt
import random
import sys

#Buy at given price
#No condition check
def buy (price):
	global flag_buy 
	global onbalance_price
	global onbalance_quant
	global onbalance_value
	global money
	flag_buy = 1
	onbalance_price = price
	onbalance_quant = fix_buy
	onbalance_value = onbalance_price * onbalance_quant
	money = money - onbalance_value
	return

#Sell at given price
#No Condition check
def sell (price):
	global flag_buy 
	global money
	global onbalance_price
	global onbalance_quant
	global onbalance_value
	global profit
	flag_buy = 0
	money = money + onbalance_quant * price
	profit = onbalance_quant * (price - onbalance_price)
	onbalance_price = 0
	onbalance_quant = 0
	onbalance_value = 0	
	return


# --------- CONTROL CENTER ------------- 
# Stock plot condition
N     =32      #grid
Npast = 2       #History Data, start with day 0
dx    = 1./N     
xa    = 0		  #Left bound
xb    = 1		  #Right bound
limit = 0.5      #rise and drop limitor

#Trading Account
money           = 100
flag_buy        = 0
fix_buy         = 10
onbalance_price = 0
onbalance_quant = 0
onbalance_value = onbalance_price * onbalance_quant

#Printer flags
flag_ini_plot     = 0	   #Plot initial graph
flag_day_plot     = 1		#Plot moving graph
flag_overall_plot = 1		#Plot overall graph
ez_overall_plot   = 0   	#Plot Ez overall graph 

flag_day_sum      = 1		#Print Each day data
flag_overall_sum  = 1		#Print overall data

time_day_plot     = 0.05   #Pause time for each day
time_ez_plot      = 1.5		#Pause time for ez graph


# --------- CONTROL CENTER ------------- 

# discrete domain
x  = geek.linspace(xa + 0.5 * dx, xb - 0.5 * dx, N)

#IC
u       = [0 for i in range(N)]

#Generate stock history
u[0] = 1 #Start price
for i in range(1,Npast):
	u[i] = u[i-1] + (random.random() * 2 - 1) * limit * u[i-1]
	if u[i] < 0:
		print "The stock just brocked up"
		sys.exit()



#initial graph
u0 =[]
for i in range(0,N):
	u0.append(u[i])
#plot initial graph
if flag_ini_plot == 1:
	plt.plot(x,u[0:N],'r*-')
	plt.xlim([0, 1])
	plt.grid(color = 'grey', linestyle = '-', linewidth = 0.5)
	plt.show()


#Generate new stock movements
t = Npast
while t < N:
	#update for next day
	u[t] = u[t-1] + (random.random() * 2 - 1) * limit * u[t-1]
	profit = 0
	#sell what we have
	if flag_buy == 1:
		sell (u[t])
	#buy new
	if flag_buy == 0:
		buy(u[t])


	#print day summary data
	if flag_day_sum == 1:
		print "-----------------day",t," -----------------"
		print "Today profit :", profit
		print "money :" , money
		print "onbalance value", onbalance_value
		print "sum :", money + onbalance_value	
	#plot day graph
	if flag_day_plot == 1:
		plt.clf   ()
		plt.plot  (x,u[0:N], 'r*-')
		plt.xlim  ([0, 1])
		plt.grid  (True)
		plt.pause (time_day_plot)
	t = t + 1


#plot overall graph
if flag_overall_plot == 1 or ez_overall_plot:
	plt.clf   ()
	plt.plot  (x,u[0:N], 'r*-')
	plt.xlim  ([0, 1])
	plt.grid  (True)
	plt.plot(x, u0[0:N+1], 'b*-')
	if ez_overall_plot == 1:
		plt.pause (time_ez_plot)
	else:
		plt.show()
#print overall summary data
if flag_overall_sum == 1:
	print "-----------------summary -----------------"
	print "money :" , money
	print "onbalance value", onbalance_value
	print "sum :", money + onbalance_value



















