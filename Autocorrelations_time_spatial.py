import os,sys
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import PlanarSubsets


#A. Reading in t sets of planar data 

Autocorrelation=PlanarSubsets.planar(path='MMC_Data/Planar_Subsets_Data')

Nt=Autocorrelation.Nt
t=Autocorrelation.t
n,N,R,U=Autocorrelation.getSample('horzArray','U')

#%%
def Spatial_AC(n,N,R,U,t,r): 

#performs a spatial autocorrelation on a planar set of velocity data (for a, x-range r) at time t and computes the integral length scale of the eddies in the flow at time t
    
    uavg=0
    vavg=0
    wavg=0 #averaged velocities over the whole plane

    for j in range(0,n):
        for i in range(0,n):
            uavg+=U[j,i,0,t]/n**2
            vavg+=U[j,i,1,t]/n**2
            wavg+=U[j,i,2,t]/n**2 #get average planar velocity components

    Uprime=np.zeros((n,n,3))
    Uprime[:,:,0]=U[:,:,0,t]-uavg 
    Uprime[:,:,1]=U[:,:,1,t]-vavg
    Uprime[:,:,2]=U[:,:,2,t]-wavg #get fluctuating velocity terms


    upsqavg=0

    for j in range(0,n):
        for i in range(0,n):
            upsqavg+=(Uprime[j,i,0])**2/N #calculate planar average of squared velocity for whole domain
            


   # r=n/2 #amount of domain sampled in x

    uprsqavg=0 #planar average of u for domain sized r*n (r in x, n in y)

    for j in range(0,n):
        for i in range(0,r):
            uprsqavg+=(Uprime[j,i,0])**2/(r*n) #calculate planar average of squared velocity for left half of the domain


               
            #B. Calculating spatial autocorrelations for half of the domain

    f_r=np.zeros((n-r))

    for k in range(0,n-r):
        for j in range(0,n):
            for i in range(0,r):
                f_r[k]+=(Uprime[j,i,0]*Uprime[j,i+k,0])/(r*n)

    f_r=f_r/uprsqavg
    L11=sp.integrate.simps(f_r[0:(n-r)], x=R[0,0:(n-r)]-R[0,0])
    
    return f_r, L11

r=n/2
f_r_spatial=np.zeros((Nt,n-r))
L11_spatial=np.zeros((Nt))

for i,ti in enumerate(t):
    f_r_spatial[i,:], L11_spatial[i]=Spatial_AC(n,N,R,U,i,r)


plt.figure()
plt.plot(t,L11_spatial[t])




#plt.figure()
#plt.plot(R[0,0:n-r]-R[0,0],f_r)
#plt.plot(R[0,0:n-r]-R[0,0],np.zeros(n-r))
#plt.ylabel('f(r)')
#plt.xlabel('r[m]')
#
#
#
#
#print 'L_11 =', L11, 'm'
