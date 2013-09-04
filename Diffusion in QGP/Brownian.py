# -*- coding: utf-8 -*-
"""
Created on Mon Sep 02 16:06:46 2013

@author: Wiliam
"""

from pylab import *
import tests


def pupdate(p,dW,mu,m,dt):
    return p - dt*mu*p/m + dW
    
def xupdate(x,xdt):
    return x + xdt
    
def kin(p,m):
    return (p/m)**2
    
def dist(x,t,n):
    return ((x**2)/n)/(2*t)
    
    

  
vecp = vectorize(pupdate)
vecx = vectorize(xupdate)
veck = vectorize(kin)
vecdist = vectorize(dist)

    

dim = 3
Kb =1.38e-23   #Bolztman const
T = 3000.      #  "Temperature
R = 3.5e-9     #"Radius
m = 1.54e-11 # "Mass
visc = 0.001  # "Viscosity
mu = 6*pi*visc*R #"Stokes law
D = mu*Kb*T  #"Einstein relation

#"Initial conditions and parameters"
t =100.  
j= 3000
N = int(10**3) 
dt = t/N
dv = 1e-21
vf = 0.6e-18

p = zeros([j,N],dtype=(float,dim))
#x = zeros([j,N],dtype=(float,dim))
sqrtD = sqrt(2*D*dt)
dW = ones([j,N],dtype=(float,dim))*sqrtD*randn(j,N,dim)
#kinetic =  zeros([1,N],dtype=(float,1))


for i in range(1,int(N)):
    p[:,i] = vecp(p[:,i-1],dW[:,i-1],mu,m,dt)
    #x[:,i] = vecx(x[:,-1],(dt/m)*p[:,i])
#    kinetic[0,i] = sum(veck(p[:,i],m))/j

    
    
#tests.MSD(p,dv,vf,T,m,j,N)
#tests.MMD(p,dv,vf,T,m,j,N)
tests.MED(p,dv,vf,T,m,j,N)   

#    
    








