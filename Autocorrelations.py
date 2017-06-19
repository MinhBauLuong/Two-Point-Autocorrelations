import os,sys
import numpy as np
from scipy.ndimage.filters import uniform_filter
import matplotlib.pyplot as plt

#A. Reading in t sets of planar data (just one for now)


path1='C:/Users/tambrico/Desktop/NWTC Internship Files/Python Assignments/MMC_Data/Planar_Subsets_Data/horzArray_U.mesh'
path2='C:/Users/tambrico/Desktop/NWTC Internship Files/Python Assignments/MMC_Data/Planar_Subsets_Data/horzArray_U.000.U'
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


uref=np.zeros((1,n)) #ref u velocity
vref=np.zeros((1,n)) #ref v velocity
wref=np.zeros((1,n)) #ref w velocity
f_rt=np.zeros((n,n))
#g_rt=np.zeros((n,n))
#h_rt=np.zeros((n,n)) #unaveraged autocorrelations for each velocity component

for j in range(0,n):
    
    uref[0,j]=U[j,0,0]
    vref[0,j]=U[j,0,1]
    wref[0,j]=U[j,0,2] #at the left wall of the simulation
    #
    for i in range(0,n):
        f_rt[j,i]=uref[0,j]*U[j,i,0]
#       g_rt[j,i]=vref*U[j,i,1]
#       h_rt[j,i]=wref*U[j,i,2] #w0*w0, w0*w1...@y=0 ; #w0*w0, w0*w1...@y=y+yi 
#       
f_rt_avg=np.zeros((n))
g_rt_avg=np.zeros((n))
h_rt_avg=np.zeros((n)) #averaged autocorrelations
uavg=0
vavg=0
wavg=0 #averaged reference velocities

for j in range(0,n):
    for i in range(0,n):
        f_rt_avg[j]+=f_rt[i,j]/n
#        g_rt_avg[j]+=g_rt[i,j]/n
#        h_rt_avg[j]+=h_rt[i,j]/n #sum the column indices and compress to a row vector
        uavg+=U[j,i,0]/n**2
        vavg+=U[j,i,1]/n**2
        wavg+=U[j,i,2]/n**2
#%%
for j in range(0,n):
    f_rt_avg[j]=f_rt_avg[j]/uavg**2
#    g_rt_avg[j]=g_rt_avg[j]/v_ref_avg**2
#    h_rt_avg[j]=h_rt_avg[j]/w_ref_avg**2 #averages and normalizes according to eq 6.45a in Pope 

#f_rt_avg is very large throughout the whole domain--is Eliot's algorithm wrong? Ask Ryan tommorow to go over the algorithm with you

plt.figure()
plt.plot(R[0,:],f_rt_avg)
plt.ylabel('f(r)')


               

            
    
    
    
    
    