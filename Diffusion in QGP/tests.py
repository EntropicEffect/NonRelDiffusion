# -*- coding: utf-8 -*-
"""
Created on Wed Sep 04 11:33:42 2013

@author: Wiliam
"""
from pylab import *

k =1.38e-23
def MSD(p,dv,vf,T,m,n,N):
    c = 4*pi*sqrt(((m)/(2*pi*k*T))**3)
    speeds = zeros([1,n])
    spdrange = arange(0,vf,dv)
    Maxwell = zeros([1,len(spdrange)])
    for i in range(n):
        speeds[0,i] = sqrt((1/m)*p[i][N-1].dot((1/m)*p[i][N-1]))
    g = 0
    for s in spdrange:
        Maxwell[0,g] = c*(s**2)*(e**((-m*s**2)/(2*k*T)))*dv
        g += 1
    weights = ones([1,n])/(n)
    hist(speeds[0],range=(0,vf),bins=vf/dv,weights=weights[0],hold=True)
    plot(spdrange,Maxwell[0])
    
def MMD(p,dm,mf,T,m,n,N):
    print("Check that drag term is non-zero")
    c= sqrt((1/(2*pi*m*k*T))**3)
    momenta = zeros([1,n])
    momrange = arange(0,mf,dm)
    Maxwell = zeros([1,len(momrange)])
    for i in range(n):
        momenta[0,i] = sqrt(p[i][N-1].dot(p[i][N-1]))
    g= 0
    for mo in momrange:
        Maxwell[0,g] = c*e**((-mo**2)/(2*m*k*T))*dm
        g += 1
    weights = ones([1,n])/(n)
    hist(momenta[0],range=(0,mf),bins=mf/dm,weights=weights[0],hold=True)
    plot(momrange,Maxwell[0]/1e30)  
    
    
def vartest(p,D,n,N,dt):
    print("Check that drag term is zero")
    variance = zeros([1,N])
    for i in range(N):
        variance[0,i] = var(p[:,i])
    xaxis = array(range(N))
    varofvar = 2.*(variance**2)/(n-1)
    variance = variance/(2*dt*xaxis)
    varofvar = varofvar*(1/(2*dt*xaxis)**2)
    plot(xaxis,D*ones([1,N])[0],hold=True)
    plot(xaxis,variance[0],hold=True)
    plot(xaxis,variance[0]+(sqrt(varofvar[0])),'--r',hold=True)
    plot(xaxis,variance[0]-(sqrt(varofvar[0])),'--b')
    
def MED(p,de,ef,T,m,n,N):  
    print("Check that drag term is non-zero")
    c = 2*sqrt((1/k*T)**3)
    energy = zeros([1,n])
    enrange = arange(0,ef,de)
    Maxwell = zeros([1,len(enrange)])
    for i in range(n):
        energy[0,i] = (p[i][N-1].dot(p[i][N-1]))/(2*m)
    g= 0 
    for en in enrange:
        Maxwell[0,g] = c*sqrt(en/pi)*e**(-en/(k*T))*de
        g+= 1
    weights = ones([1,n])*de
    hist(energy[0],range=(0,ef),bins=ef/de,normed=1,hold=True)    
    plot(enrange,Maxwell[0]/(2*m))   
    
    
     