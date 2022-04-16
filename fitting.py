# -*- coding: utf-8 -*-
import scipy.optimize as so
import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import curve_fit,least_squares



#global k_b
#global j

H_hf0 = 415.46 # KOe
k_b =  1#(1.380649*10**(-23))/ (9.27*10**(-24))  #J⋅K−1
S_fe = 1
S_R = 2.5
g_fe = 2
j_FeR = 16.3
Z_RFei = 6# variable
gamma = 2.5
mu_R0 = 1.6#*9.27*10**(-24) #mu_B
beta = 1/150 #mu_B/KOe
mu_FeT= 2.2 #mu_B

H1 = []
H2 = []
mu_RT = []
y = []
z = []

# Experimental x and y data points
x = np.array([100,150,200,250,300,350,400,450,500,600,650,700])   #TData
y1 = np.array([415.46, 412.5, 407.6, 401.7, 388.8, 371.9, 358.1, 338.2, 324.3, 307.5, 258.7, 211.03])  #_6cData
y2 = np.array([379.6, 376.7, 377.7, 373.8, 362.9, 351.1, 334.2, 318.3, 306.4, 283.6,219.9,173.2])  #_9dData
y3 = np.array([356.7,354.7,350.8,345.9,334.04,319.2,310.3,296.4,273.6,246.7,218.9,159.24])  #_18fData
y4 = np.array([312.9, 308.9, 310, 307.1, 297.2, 281.34, 273.43,257.6,241.7,220.8,187.05,159.24])  #_18hData

def func(X, a, b, c, d):
    x,y1,y2,y3,y4 = X
    #for i in range(len(x)):
        #y1 = (1 + 0.5*S_fe) * (S_fe * np.tanh((2 * a *beta * y1 + 6 * b *beta * y2 + 12 * c * beta * y3 + 6 * d *beta * y4) * (S_fe + 0.5) / S_fe)) - 0.5 * (S_fe * np.tanh(0.5 * (2 * a *beta * y1 + 6 * b *beta * y2 + 12 * c *beta * y3 + 6 * d *beta * y4) / S_fe))
    H1 = -(S_R*j_FeR)/(k_b*x*g_fe*150)*(7*y1+6*y2+5*y3+5*y4) #Z_FeRi (7,6,5,5) # nbre des proches voisins(R) pour un Fe donné
    mu_RT = mu_R0*((1 + 0.5/S_R)*(1/np.tanh((1 + 0.5/S_R)*H1)) -(1/(2*S_R))*(1/np.tanh(H1/(2*S_R))))
    H2 = -((S_fe*beta)/(k_b*x*g_fe))*(2*a*y1+6*b*y2+12*c*y3+6*d*y4) + (S_fe/x)*(j_FeR*Z_RFei*gamma*mu_RT)
    y = H_hf0*(1 + 0.5/S_fe)*(1/np.tanh((1 + 0.5/S_fe)*H2)) -(H_hf0/(2*S_fe))*(1/np.tanh(H2/(2*S_fe)))
    return y




# initial guesses for a,b,c, d:
p0 = -200., 35., 15., 70.

popt,pcov = curve_fit(func,(x,y1,y2,y3,y4),y1, p0)
print(popt)

    
    
# Plot experimental data points
plt.plot(x, y1, 'ro', label='experimental-data')

# Plot the fitted function
z1 = func((x,y1,y2,y3,y4),  *popt)
plt.plot(x, z1 , 'r')



#print(res)

plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()



