#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Spatial inference on a single snapshot using a pair-wise correlation derived 
pseudo-likelihood.
Created on Sat Apr 15 22:40:39 2017

@author: m.irvine
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import imfilter
from scipy.ndimage.filters import convolve
from scipy.stats import expon
from numpy.matlib import repmat

expl = 1E-5
calculate_corrs = True

from bandingsim import shortenKernel, genKernel, calcBD, generateModelData

#determine means for weak priors (either through prior knowledge or correlated
#with size of lattice)
def rpriorfinder(A,i):
    betar = 10
    betak = 10
    width = 1
    betatheta = 0
    return betar,betak,width,betatheta
    
#BiB,BiD,DeB,DeD,N00,N01,N11  = BDcorrboundary( A,bbox );
def extend2Axis(a,n):
    a = np.repeat(a[:, np.newaxis], n, axis=1)
    a = np.repeat(a[:, :,np.newaxis], n, axis=2)
    return a
    
def BDcorrboundary(A,bbox):
    BiB,BiD,DeB,DeD,N00,N01,N11 = BDcorr(A)
    N,M = A.shape
    mid = int(N/2)
    m2 = int(bbox/2)
    BiB = BiB[mid-m2:mid+m2,mid-m2:mid+m2,:bbox]
    BiD = BiD[mid-m2:mid+m2,mid-m2:mid+m2,:bbox]
    DeB = DeB[mid-m2:mid+m2,mid-m2:mid+m2,:bbox]
    DeD = DeD[mid-m2:mid+m2,mid-m2:mid+m2,:bbox]
    N00 = N00[:bbox,:bbox,:bbox]
    N01 = N01[:bbox,:bbox,:bbox]
    N11 = N11[:bbox,:bbox,:bbox]
    return BiB,BiD,DeB,DeD,N00,N01,N11


def BDcorr( A ):

    # m parameter defines sizes of box over which to take correlations. Used
    # for real data where need to seriously worry about boundary conditions
    n,m=A.shape;
    rho1 = np.sum(A)
                  
    rho0 = np.sum(1-A)
    k = np.linspace(-n,n,n)
    yy,xx = np.meshgrid(k,k)
    r = np.sqrt(xx**2+yy**2)
    l = 2
    #declare correlations memory
    BiB = np.zeros((n,n,n)) #correlations in births locations to births
    DeB = np.zeros((n,n,n)) #correlations in death loacations to births
    BiD = np.zeros((n,n,n)) #correlations in births locations to deaths
    DeD = np.zeros((n,n,n)) #correlations in deaths locations to dearth
    D = np.arange(n)
    #normalisation constants for w_xx(d)
    
    N00=np.zeros(n)
    N01=np.zeros(n)
    N11=np.zeros(n)
    
    
    for d in D:
        K = ((d<=r)*(r<=(d+l))).astype(float)
        convd =  convolve(A,K,mode='wrap')
        N = np.sum(convolve(np.ones((n,n)),K,mode='wrap'))
        bib = convd*A
        deb = convd*(1-A);
        N11[d] = N*(rho1+rho1)
        N01[d] = N#%*(rho0+rho1);%sum(deb(:));
        BiB[:,:,d] = bib
        DeB[:,:,d] = deb
        
        
        convd =  convolve(1-A,K,mode='wrap')
        bid = convd*A
        ded = convd*(1-A)
        N00[d] = N*(rho0+rho0)
        BiD[:,:,d] = bid;
        DeD[:,:,d] = ded;
        
        
    N00 = np.transpose(extend2Axis(N00,n),(1, 2, 0))
    N01 = np.transpose(extend2Axis(N01,n),(1, 2, 0))
    N11 = np.transpose(extend2Axis(N11,n),(1, 2, 0))
    return BiB,BiD,DeB,DeD,N00,N01,N11
    

    
    
    

#import data
#A = np.random.rand(101,101)<0.5
#A = A.astype(float)
print('generating data...')
l1 = 1.0
l2 = 2.0
k = 1.0
r = 10.
theta = 0#np.pi/4
A = generateModelData(l1,l2,r,theta,k,N=100,T=500)

N,M = A.shape
N=min(N,M)
M=N

kx = np.arange(-(N-1)/2,(N-1)/2+1)#((-(N-1)/2):((N-1)/2)); 
ky = np.arange(-(M-1)/2,(M-1)/2+1)#((-(M-1)/2):((M-1)/2));
[xx,yy] = np.meshgrid(ky,kx)
N=max(N,M)

# define MCMC parameters
burnt,inft = 0,int(10000)

#declare memory for parameters and EDD records
rl2 = np.zeros(inft+1)
rl1 = np.zeros(inft+1)
rk = np.zeros(inft+1)
rr = np.zeros(inft+1)
rtheta = np.zeros(inft+1)
reta = np.zeros(inft+1)
rlambda = np.zeros(inft+1)
ais = np.zeros(inft+1)

EDD = np.zeros(inft+1)
prior = np.zeros(inft+1)


bbox = 80 #size of interior box over which zeta is calculated
#TODO : creat these functions
if calculate_corrs:
    print('Calculating correlations...')
    BiB,BiD,DeB,DeD,N00,N01,N11  = BDcorrboundary( A,bbox );
    betar,betak,width,betatheta  = rpriorfinder( A,0 );

#initialise parameters and EDD
print('Initialising MCMC chain...')
#r = np.random.exponential(betar);
#theta = np.mod(np.random.normal(betatheta,np.pi/16),np.pi)
#l2 = np.random.exponential(width); #calculate from distance between bands??
#l1 = np.random.exponential(width); #these need to be widths, calculate one from width of band
#k = 1.#np.random.exponential(betak); #hardest one to form a prior on. Maybe relate it to the strength of the signal in the power spectrum. i.e. if it is weak then competition probably is as well.

eps = 0.01
B,D = calcBD(A,N,l1,l2,r,theta,k,xx,yy,eps=eps)

mid = int(N/2);
m2 = int(bbox/2);
B=B[mid-m2:mid+m2,mid-m2:mid+m2];
D = D[mid-m2:mid+m2,mid-m2:mid+m2];
BB = np.repeat(B[:, :, np.newaxis], bbox, axis=2)
DB = np.repeat(D[:, :, np.newaxis], bbox, axis=2)
nn= N01;
H1 = (-BB*DeD+DB*BiD); #P00
H2 = (-DB*BiB+BB*DeB); #P11
H3 = (DB*BiB-DB*BiD-BB*DeB+BB*DeD); #H3 = (2*D.*BiB-D.*BiD-B.*DeB+2*B.*DeD); %P01
#H = -B*(DeD)+D*(BiD)-D*(BiB)+B*(DeB)-D*(BiD)+B*(DeD);
H=np.squeeze(np.sum(np.sum(H1,0),1)**2)/N00**2+np.squeeze(np.sum(np.sum(H3,0),1)**2)/N11**2;
#pEDD = exp(-sum(H(2:end))/temp);
EDD[0] = -np.sum(H[2:])/expl**2

prior[0]=-(l1/width+l2/width+k/betak+r/betar)+ np.log(r<25)


ind=0

for t in np.arange(1,(burnt+inft)):
    #initialise parameters and EDD
    pr,ptheta,pl1,pl2,pk = r,theta,l1,l2,k
    
    if t%4 == 0: pr = np.random.normal(r,0.5);
    if t%4 == 0: ptheta = np.mod(np.random.normal(theta,np.pi/16),np.pi)
    if t%4 == 1: pl2 = np.random.normal(l2,0.01) #calculate from distance between bands??
    if t%4 == 2: pl1 = np.random.normal(l1,0.01) #these need to be widths, calculate one from width of band
    if t%4 == 3: pk = np.random.normal(k,0.01) #hardest one to form a prior on. Maybe relate it to the strength of the signal in the power spectrum. i.e. if it is weak then competition probably is as well.

    B,D = calcBD(A,N,pl1,pl2,pr,ptheta,pk,xx,yy,eps=eps)

    mid = int(N/2);
    m2 = int(bbox/2);
    B=B[mid-m2:mid+m2,mid-m2:mid+m2];
    D = D[mid-m2:mid+m2,mid-m2:mid+m2];
    BB = np.repeat(B[:, :, np.newaxis], bbox, axis=2)
    DB = np.repeat(D[:, :, np.newaxis], bbox, axis=2)
    H1 = (-BB*DeD+DB*BiD); #P00
    H2 = (-DB*BiB+BB*DeB); #P11
    H3 = (DB*BiB-DB*BiD-BB*DeB+BB*DeD); #H3 = (2*D.*BiB-D.*BiD-B.*DeB+2*B.*DeD); %P01
    #H = -B*(DeD)+D*(BiD)-D*(BiB)+B*(DeB)-D*(BiD)+B*(DeD);
    H=np.squeeze(np.sum(np.sum(H1,0),1)**2)/N00**2+np.squeeze(np.sum(np.sum(H3,0),1)**2)/N11**2;
    #H=np.sum(np.sum(H1/nn,0),1)**2+ np.sum(np.sum(H2/nn,0),1)**2+np.sum(np.sum(H3/nn,0),1)**2;
    H = np.squeeze(H);
    #pEDD = exp(-sum(H(2:end))/temp);
    pEDD = -np.sum(H[2:])/expl**2
    
    pprior=-(pl1/width+pl2/width+pk/betak+pr/betar)+ np.log(pr<25) + np.log(pr>0) + np.log(pl1>0) + np.log(pl2>0) + np.log(pk>0)
    a = pEDD+ pprior - EDD[ind]-prior[ind];
    ais[ind] = a;

    if a>= 0: #a<rand
        EDD[ind+1] = pEDD #update new EDD
        prior[ind+1] = pprior
        #update to proposal parameters
        l1 = pl1;
        l2 = pl2; 
        
        r = pr;
        theta = ptheta;
        
        k = pk;
    elif np.log(np.random.rand())<a :
        EDD[ind+1] = pEDD #update new EDD
        prior[ind+1] = pprior
        #update to proposal parameters
        l1 = pl1;
        l2 = pl2; 
        
        r = pr;
        theta = ptheta;
        
        k = pk;
    else:
        prior[ind+1] = prior[ind];
        EDD[ind+1] = EDD[ind];
        #don't update proposal parameters

    if t> 0: #if greater than burn time record parameters.
        rl2[ind] = l2;
        rl1[ind] = l1;
        rk[ind] = k;
        rr[ind] = r;
        rtheta[ind]=theta;



    
    if np.mod(t,10)==0:
        print('index: {}, EDD: {}'.format(ind,EDD[ind+1]))
        print('Current parameters: r {}, l1 {}, l2 {}, k {}, theta {}'.format(r,l1,l2,k,theta))
    ind = ind+1;
    
    