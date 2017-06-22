import os,sys
import numpy as np
from scipy.ndimage.filters import uniform_filter
import matplotlib.pyplot as plt
import scipy as sp


#A. Reading in t sets of planar data (just one for now)


path1='C:/Users/tambrico/Desktop/NWTC Internship Files/Python Assignments/MMC_Data/Planar_Subsets_Data/30100/horzArray_U.mesh'
path2='C:/Users/tambrico/Desktop/NWTC Internship Files/Python Assignments/MMC_Data/Planar_Subsets_Data/30100/horzArray_U.000.U'
with open(path1, 'r') as f:
    
    for i in range(1, 9):
        f.readline()

    N=int(f.readline())
    n=int(np.sqrt(N))
    R=np.zeros((n,n)) #xy coordinate array
    U=np.zeros((n,n,3))#velocity array, stores u,v,w at each point in R
    
    for i in range(0, n): #reads in the first 99 mesh coordinates. Since its a square, the points are fed into both the first row and the first col and then the col is repeated over
        R[0,i]=float(f.readline())
        if i==0: continue
        R[i,:]=R[0,i]
    
    

with open(path2, 'r') as f:
    
    for i in range(0,4):
        f.readline()
        
    for k in range(0,3):    
        for j in range(0,n):
            for i in range(0,n):
                U[j,i,k]=float(f.readline()) #stores u, v, w in a planar array corresponding to R. Scans x for fixed y, then scans y, then scans each velocity component

uavg=0
vavg=0
wavg=0 #averaged velocities over the whole plane

for j in range(0,n):
    for i in range(0,n):
         uavg+=U[j,i,0]/n**2
         vavg+=U[j,i,1]/n**2
         wavg+=U[j,i,2]/n**2 #get average planar velocity components

Uprime=np.zeros((n,n,3))
Uprime[:,:,0]=U[:,:,0]-uavg 
Uprime[:,:,1]=U[:,:,1]-vavg
Uprime[:,:,2]=U[:,:,2]-wavg #get fluctuating velocity terms


upsqavg=0

for j in range(0,n):
    for i in range(0,n):
        upsqavg+=(Uprime[j,i,0])**2/N #calculate planar average of squared velocity for whole domain



r=n/2 #amount of domain sampled in x

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


plt.figure()
plt.plot(R[0,0:n-r]-R[0,0],f_r)
plt.plot(R[0,0:n-r]-R[0,0],np.zeros(n-r))
plt.ylabel('f(r)')
plt.xlabel('r[m]')


L11=sp.integrate.simps(f_r[0:(n-r)], x=R[0,0:(n-r)]-R[0,0])

print 'L_11 =', L11, 'm'

