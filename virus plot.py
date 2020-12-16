# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 10:17:09 2020

@author: matt hanson
"""
import numpy as np
import matplotlib.pyplot as plt

#function to generate points to graph
#finds the integral of the virus_dirivate
def rk4(func, x0, y0, xN, N, args=()):
    h = (xN-x0)/N
    t = np.linspace(x0, xN, N+1)
    Y = np.zeros((len(t),len(y0)))
    Y[0] =y0
    for n in range(N): 
        k1 = func(Y[n], *args)
        k2 = func(Y[n]+0.5*h*k1, *args)
        k3 = func(Y[n]+0.5*h*k2, *args)
        k4 = func(Y[n]+k3*h, *args)
        Y[n+1]=Y[n] + h*(k1+2*k2+2*k3+k4)/6
    return t,Y

#finds how the infection and susceptible populations are chaging at a certain point
def virus_dirivitive(y,b,a):
    i=y[0]
    s=y[1]
    ds=-b*s*i
    di=b*s*i-a*i
    dr=a*i
    return np.array([di,ds,dr])
    
#i0=initial infected population
i0=(8.2*10)**(-4)

#s0=initial susceptible population
s0=1-i0
r0=0
y0=np.array([i0,s0,r0])
t0=0
tN=1000
N=100000

#b_a=(β=.214,γ=1/14)
b_a=(.214,.0714)

#plot1
t,y1=rk4(virus_dirivitive, t0, y0, tN, N,(b_a))
s=y1[:,1]
i=y1[:,0]
r=y1[:,2]
plt.subplot(2,2,1)
plt.plot(t,i, 'g', label = 'infected ')
plt.plot(t,s, 'b', label = 'Susceptibles ')
plt.plot(t,r, 'r', label = 'recovered ')
plt.xlabel('time in days')
plt.ylabel('fraction of population')
plt.title('β=.214,γ=1/14, i0=8*10**-4')
plt.legend(loc='best')

#plot 2
b_a=(.114,.0714)
t,y1=rk4(virus_dirivitive, t0, y0, tN, N,(b_a))
s=y1[:,1]
i=y1[:,0]
r=y1[:,2]
plt.subplot(2,2,2)
plt.plot(t,i, 'g',)
plt.plot(t,s, 'b',)
plt.plot(t,r, 'r',)
plt.xlabel('time in days')
plt.ylabel('fraction of population')
plt.title('β=.114,γ=1/14, i0=8*10**-4')

#plot 3
b_a=(.18,.0714)
t,y1=rk4(virus_dirivitive, t0, y0, tN, N,(b_a))
s=y1[:,1]
i=y1[:,0]
r=y1[:,2]
plt.subplot(2,2,3)
plt.plot(t,i, 'g',)
plt.plot(t,s, 'b',)
plt.plot(t,r, 'r',)
plt.xlabel('time in days')
plt.ylabel('fraction of population')
plt.title('β=.18,γ=1/14, i0=8*10**-4')

#plot 4
i0=(8.2*10)**(-4)
s0=1-i0-.9
r0=0
y0=np.array([i0,s0,r0])
b_a=(.114,.0714)
t,y1=rk4(virus_dirivitive, t0, y0, tN, N,(b_a))
s=y1[:,1]
i=y1[:,0]
r=y1[:,2]
plt.subplot(2,2,4)
plt.plot(t,i, 'g',)
plt.plot(t,s, 'b',)
plt.plot(t,r, 'r',)
plt.xlabel('time in days')
plt.ylabel('fraction of population')
plt.title('β=.214,γ=1/14, i0=8*10**-4,s0=.1')
plt.tight_layout()
