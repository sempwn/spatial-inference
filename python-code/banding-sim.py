#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 16:59:56 2017

@author: m.irvine
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import convolve

def shortenKernel(l,N):
    k = np.floor(np.sqrt(2*np.log(1/0.025))*l**2)+1; #this selects the size of the kernel such that the mass of a unit outside of it is less than 0.025
    k = k + (np.mod(k,2)); #increses by one if k1 odd in order to do next operation
    if k > N-2: k = N-2

    return int(k)

def genKernel(l,xx,yy,k,N):
    #create kernel shortened to k.
    ker = np.exp(-1*(.5/l**2)* (xx**2+yy**2) );
    ker = ker/np.sum(ker[:]);
    mink,maxk = int(np.floor((N-k)/2)),int(np.ceil((N+k)/2))+1
    ker = ker[mink:maxk,mink:maxk];
    return ker
    
def calcBD(A,N,l1,l2,r,theta,k,xx,yy,eps=0.01):
    k1 = shortenKernel(l1,N)
    ker1 = genKernel(l1,xx,yy,k1,N)

    k2 = shortenKernel(l2,N)
    ker2 = genKernel(l2,xx,yy,k2,N)

    P1 = convolve(A,ker1,mode='wrap')

    P2 = np.maximum(0,convolve(A,ker2,mode='wrap')) 
    oxs = int(r*np.cos(theta))
    oys = int(r*np.sin(theta))
    P2 = np.roll(np.roll(P2,oxs,axis=0),oys,axis=1)
    B = 1-np.exp(-eps*P1);#prob of birth event
    D = 1-np.exp(-eps*k*P2);#prob of death event
    return B,D

def generateModelData(l1,l2,r,theta,k,N=100,T=1):
    A = (np.random.rand(N,N)<0.1).astype(float)
    spacek = np.linspace(-N,N,N)
    yy,xx = np.meshgrid(spacek,spacek)

    for t in np.arange(T):
        B,D = calcBD(A,N,l1,l2,r,theta,k,xx,yy,eps=0.01)
        births = np.random.rand(N,N)<B
        deaths = np.random.rand(N,N)<D
        A = A + births*(1-A) - deaths*A
        #print('births {}, deaths {}'.format(np.sum(births),np.sum(deaths)))
        #A = A - A*(A<0) + (1-A)*(A>1) 

        
    return A,B,D
    
if __name__ == '__main__':
    l1 = 1.0
    l2 = 1.0
    k = 10.0
    r = 10.
    theta = np.pi/4
    A,B,D = generateModelData(l1,l2,r,theta,k,N=80,T=500)
    plt.imshow(A,cmap='Greys')
    B,D = calcBD(A,N,l1,l2,r,theta,k,xx,yy,eps=0.01)
    plt.figure()
    plt.imshow(B)
    plt.title('birth rate')
    
    plt.figure()
    plt.imshow(D)
    plt.title('death rate')
    