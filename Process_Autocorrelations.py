import numpy as np
import matplotlib.pyplot as plt

def getfnames(fname):
    fnames=[]
    f=file(fname)
    for line in f:
        fnames.append(line.rstrip('\n'))
        
    return fnames
            

def getDimensions(fname):
    with open(fname, 'r') as f:
        n_x=int(float(f.readline()))
        n_y=int(float(f.readline()))
        N=int(float(f.readline()))
        Nt=int(float(f.readline()))
        r=int(float(f.readline()))
    dim=[n_x,n_y,3,Nt,N,r]
    return dim

def getArray(fname, dim, Xrank=1): #pick which array to retrieve and process
    
    with open(fname, 'r') as f:
        
        if Xrank==1:
            n=dim
            X=np.zeros((n))
            for i in range(0,n):
                X[i]=float(f.readline().rstrip('\n'))
        
        if Xrank==2:
            n_x=dim[0]
            n_y=dim[1]
            X=np.zeros((n_y,n_x))
            for j in range(0,n_y):
                for i in range(0,n_x):
                    X[j,i]=float(f.readline().rstrip('\n'))
                    
        elif Xrank==3:
            n_x=dim[0]
            n_y=dim[1]
            X=np.zeros((n_y,n_x,3))
            for k in range(0,n_y):
                for j in range(0,n_x):
                    for i in range(0,3):
                        X[k,j,i]=float(f.readline().rstrip('\n'))
        
        elif Xrank==4:
            n_x=dim[0]
            n_y=dim[1]
            Nt=dim[3]
            n=0
            X=np.zeros((n_y,n_x,3,Nt))
            for ti in range(0,n_y):
                for k in range(0,n_x):
                    for j in range(0,3):
                        for i in range(0,Nt):
                            X[ti,k,j,i]=float(f.readline().rstrip('\n'))
                t=t+1
                print 'n=',n
    return X


fnames=getfnames('MMC_Data/planarSubsets/fnames.txt')
dim=getDimensions(fnames[4])
n_x=dim[0]
n_y=dim[1]
Nt=dim[3]
N=dim[4]
r=dim[5]


#%%  Plotting L11(t)

x=getArray(fnames[5],dim[0])
y=getArray(fnames[6],dim[1])
t=getArray(fnames[9],dim[3])

#%%
f_r_spatial_dim=[dim[0]-dim[-1],dim[3]]
f_r_spatial=getArray(fnames[2],f_r_spatial_dim,Xrank=2)

L11_spatial=getArray(fnames[3],dim[3])

plt.figure()
plt.plot(x[0:n_x-r]-x[0],f_r_spatial[-1,:])
plt.plot(x[0:n_x-r]-x[0],np.zeros(n_x-r))
plt.ylabel('f(r)')
plt.xlabel('r[m]') #plot f(r) for one time instance
plt.title('f(r)_spatial at t=41000s for r=n_x/4')
plt.savefig('C:/Users/tambrico/Desktop/NWTC_Internship_Files/Python_Assignments/MMC_Data/Planar_Subsets_Plots_2/' + 'fr_tend' +'.png', format="png")

print 'L_11[t=41000s] =', L11_spatial[-1], 'm'

plt.figure()
plt.plot(t[:],L11_spatial[:])
plt.xlabel('t[s]')
plt.ylabel('L11_S(t) [m]') #plot integral length scale for all time instances
plt.title('L11(t) vs t for r=n_x/4')
plt.savefig('C:/Users/tambrico/Desktop/NWTC_Internship_Files/Python_Assignments/MMC_Data/Planar_Subsets_Plots_2/' + 'L11s(t)' +'.png', format="png")
#%%

# Plot the U' field for all t and save in a directory

Uprime=getArray(fnames[1],dim[0:4],Xrank=4)
#%%
for i, ti in enumerate(t):
    fig, (ax1) = plt.subplots(nrows=1)
    im = ax1.pcolormesh(x, y, Uprime[:,:,0,i], cmap='coolwarm')
    plt.colorbar(im, orientation='vertical')
    ax1.set_aspect('equal')
    ax1.autoscale(tight=True)
    plt.savefig('C:/Users/tambrico/Desktop/NWTC_Internship_Files/Python_Assignments/MMC_Data/Planar_Subsets_Plots_2/Uprime_field/' +  str(ti) +'.png', format="png") #save all the field data in a directory
    print 't=',ti
#U=getArray(fnames[0],dim[0:4])
#f_r_st=getArray(fnames[7],dim[0]-dim[-1])
#L11_st=getArray(fnames[8],1)
    
                            
    
        
            
        