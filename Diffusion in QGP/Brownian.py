# -*- coding: utf-8 -*-
"""
Created on Mon Sep 02 16:06:46 2013

@author: Wiliam
"""

from pylab import *


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

    

dim = 2
Kb =1.38e-23   #Bolztman const
T =300000.      #  "Temperature
R = 3.5e-9     #"Radius
m = 1.54e-11 # "Mass
visc = 0.001  # "Viscosity
mu = 6*pi*visc*R #"Stokes law
D = mu*Kb*T  #"Einstein relation

#"Initial conditions and parameters"
t = 1.  
j= 100
N = (10**4)  
dt = t/N

p = zeros([j,N],dtype=(float,dim))
x = zeros([j,N],dtype=(float,dim))
dW = ones([j,N],dtype=(float,dim))*sqrt(2*D*dt)*randn(j,N,dim)
kinetic =  zeros([1,N],dtype=(float,1))
Dist = zeros([1,N],dtype=(float,1))
kinetic[0,0] = sum(veck(p[:,0],m))/j
variance = zeros([1,N],dtype=(float,1))



for i in range(1,int(N)):
    p[:,i] = vecp(p[:,i-1],dW[:,i-1],mu,m,dt)
    x[:,i] = vecx(x[:,-1],(dt/m)*p[:,i])
    kinetic[0,i] = sum(veck(p[:,i],m))/j
    Dist[0,i] = sum(vecdist(x[:,i],i*dt,j))
    variance[0,i] = var(x[:,i])
    

    
    

    
    








