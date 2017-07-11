import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


#A. Reading in t sets of planar data (just one for now)


path1='C:/Users/tambrico/Desktop/NWTC_Internship_Files/Python_Assignments/MMC_Data/planarSubsets/41000/horzArray_U.mesh'
path2='C:/Users/tambrico/Desktop/NWTC_Internship_Files/Python_Assignments/MMC_Data/planarSubsets/41000/horzArray_U.000.U'

with open(path1, 'r') as f: #get x coordinates
    
    for i in range(1, 9): f.readline()

    N=int(f.readline())
    x=[]
    y=[]#xy coordinate array
     
    
    for i in range(0,N): 
        x.append(float(f.readline()))
        if x[i]==x[0] and i!=0: break
    del x[-1]
    xcoord=np.array(x)
    n_x=len(xcoord) 

with open(path1, 'r') as f:
    for i in range(1, 9):f.readline()
    
    n_y=0
    for i in range(0,N):
        if float(f.readline())==x[0]: n_y+=1 #xcoordinates are repeated n_y times
    
    f.readline()
    for j in range(0,n_y):
        y.append(float(f.readline()))
        for i in range(0,n_x):
            f.readline()
            if len(y)==n_y: break
        if len(y)==n_y: break
        
    
ycoord=np.array(y)
U=np.zeros((n_y,n_x,3)) #velocity components array

   #%% 

with open(path2, 'r') as f:
    
    for i in range(0,4):
        f.readline()
        
    for k in range(0,3):    
        for j in range(0,n_y):
            for i in range(0,n_x):
                U[j,i,k]=float(f.readline()) #stores u, v, w in a planar array corresponding to R. Scans x for fixed y, then scans y, then scans each velocity component
#%%
uavg=0
vavg=0
wavg=0 #averaged velocities over the whole plane

for j in range(0,n_y):
    for i in range(0,n_x):
         uavg+=U[j,i,0]/N
         vavg+=U[j,i,1]/N
         wavg+=U[j,i,2]/N #get average planar velocity components

Uprime=np.zeros((n_y,n_x,3))
Uprime[:,:,0]=U[:,:,0]-uavg 
Uprime[:,:,1]=U[:,:,1]-vavg
Uprime[:,:,2]=U[:,:,2]-wavg #get fluctuating velocity terms


upsqavg=0

for j in range(0,n_y):
    for i in range(0,n_x):
        upsqavg+=(Uprime[j,i,0])**2/N #calculate planar average of squared velocity for whole domain



r=n_x/4 #amount of domain sampled in x

uprsqavg=0 #planar average of u for domain sized r*n (r in x, n in y)

for j in range(0,n_y):
    for i in range(0,r):
        uprsqavg+=(Uprime[j,i,0])**2/(2*r*n_y)

for j in range(0,n_y):
    for i in range(n_x,n_x-r):
        uprsqavg+=(Uprime[j,i,0])**2/(2*r*n_y)


#calculate planar average of squared velocity for an rth of the domain on each side
        

               
#B. Calculating spatial autocorrelations for 2 rths of the domain; captures f(r) up to f(n-r)

f_r=np.zeros((n_x-r))

for k in range(0,n_x-r):
    print 'k=',k
    for j in range(0,n_y):
        print 'j=',j
        for i in range(0,r):
            f_r[k]+=(Uprime[j,i,0]*Uprime[j,i+k,0])/(2*r*n_y)

for k in range(0,n_x-r):
    for j in range(0,n_y):
        for i in range(n_x,n_x-r):
            f_r[k]+=(Uprime[j,i,0]*Uprime[j,i-k,0])/(2*r*n_y)

f_r=f_r/uprsqavg


plt.figure()
plt.plot(xcoord[0:n_x-r]-xcoord[0],f_r)
plt.plot(xcoord[0:n_x-r]-xcoord[0],np.zeros(n_x-r))
plt.ylabel('f(r)')
plt.xlabel('r[m]')


L11=sp.integrate.simps(f_r[0:(n_x-r)], x=xcoord[0:(n_x-r)]-xcoord[0])

print 'L_11 =', L11, 'm'

