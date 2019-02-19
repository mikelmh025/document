# Minghao Liu 
# Buy and sell without algorithm
# Only buy first
# Only one transaction at a time

from __future__ import division

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

#input two flux value 
#Counter for 4 types: 1. U U,  2. U D,  3.D U,  4. D D
#update type of this pattern to counter
def analyze (price1, price2):
	global counter_1 
	global counter_2 
	global counter_3 
	global counter_4 
	global counter_sum
	if price1 >= 0:
		if price2 >= 0:
			counter_1 = counter_1 + 1
		else:
			counter_2 = counter_2 + 1
	else:
		if price2 >= 0:
			counter_3 = counter_3 + 1 
		else:
			counter_4 = counter_4 + 1
	counter_sum = counter_sum + 1
	return

def erase(price1, price2):
	global counter_1 
	global counter_2 
	global counter_3 
	global counter_4 
	global counter_sum
	if price1 >= 0:
		if price2 >= 0:
			counter_1 = counter_1 - 1
		else:
			counter_2 = counter_2 - 1
	else:
		if price2 >= 0:
			counter_3 = counter_3 - 1 
		else:
			counter_4 = counter_4 - 1
	counter_sum = counter_sum - 1
	return


# --------- CONTROL CENTER ------------- 
# Stock plot condition
N     = 64      #grid
Npast = 32      #History Data, start with day 0
dx    = 1./N     
xa    = 0		  #Left bound
xb    = 1		  #Right bound
limit = 0.1      #rise and drop limitor

#analyze
lookback_range = Npast - 10

#Trading Account
money           = 100
flag_buy        = 0
fix_buy         = 10
onbalance_price = 0
onbalance_quant = 0
onbalance_value = onbalance_price * onbalance_quant

#Printer flags
flag_ini_plot     = 1	   #Plot initial graph
flag_day_plot     = 1		#Plot moving graph
flag_overall_plot = 1		#Plot overall graph
ez_overall_plot   = 0   	#Plot Ez overall graph 

flag_day_sum      = 0		#Print Each day data
flag_overall_sum  = 1		#Print overall data

time_day_plot     = 0.05   #Pause time for each day
time_ez_plot      = 1.5		#Pause time for ez graph


# --------- CONTROL CENTER ------------- 

# discrete domain
x  = geek.linspace(xa + 0.5 * dx, xb - 0.5 * dx, N)

#IC
u       = [0 for i in range(N)]
uC      = [0 for i in range(N)]
#Generate stock history
u[0] = 1 #Start price
for i in range(1,Npast):
	u[i] = u[i-1] + (random.random() * 2 - 1) * limit * u[i-1]
	if u[i] < 0:
		print "The stock just brocked up"
		sys.exit()


#Impliment 2 step look back analysis to history data
#Counter for 4 types: 1. U U,  2. U D,  3.D U,  4. D D
counter_1 = 0
counter_2 = 0
counter_3 = 0
counter_4 = 0
counter_sum = 0
for i in range (1,Npast):
	uC[i] = u[i] - u[i-1]

for i in range (1,Npast):
	analyze(uC[i-1],uC[i])
	if i >= lookback_range:
		#print "erase :", i-lookback_range, i+1-lookback_range
		erase(uC[i-lookback_range],uC[i+1-lookback_range])
#print uC


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
	u[t]  = u[t-1] + (random.random() * 2 - 1) * limit * u[t-1]
	uC[t] = u[t] - u[t-1]
	analyze(uC[t-1],uC[t])
	erase(uC[t-lookback_range],uC[t+1-lookback_range])
	profit = 0

	#Decision making
	
	p = 0
	
	'''
	# Print counter
	print "counter_sum : ", counter_sum
	print "counter_1: ", counter_1
	print "counter_2: ", counter_2
	print "counter_3: ", counter_3
	print "counter_4: ", counter_4	
	print "uC[t]: ",uC[t]
	'''
	flag_action = 0
	if uC[t] >= 0:
		if counter_1 == 0 and counter_2 == 0:
			flag_action = 0
		else:
			p = (counter_1 / counter_sum) / ((counter_1 / counter_sum)+(counter_2 / counter_sum))
	else:
		if counter_3 == 0 and counter_4 == 0:
			flag_action = 0
		else:
			p = (counter_4 / counter_sum) / ((counter_4 / counter_sum)+(counter_3 / counter_sum))
	#print "p: ", p
	
	if p >= 0.51:
		flag_action = 1

	#sell what we have
		if flag_buy == 1:
			sell (u[t])
	if flag_action == 1:
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




'''
#Show Counter
print counter_1
print counter_2
print counter_3
print counter_4

'''









