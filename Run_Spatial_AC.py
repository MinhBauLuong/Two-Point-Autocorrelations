import os
import numpy as np
import PlanarSubsets as ps
import matplotlib.pyplot as plt

dpath='/scratch/equon/inflowOutflow/neutral/8mps_shear0.2_TI10.0_5km/postProcessing/planarSubsets'
fpath='/home/tambrico/MMC/ILS_Data'
AC=ps.planar_processing(path=dpath)


#%% Retrieve Parameters from txtfiles

vfnames=AC.getfnames(fpath+'/'+AC.case+'/vfnames.txt')
dim=AC.getDimensions(vfnames[0])
n_x=dim[0]
n_y=dim[1]
N=dim[2]
Nt=dim[3]
r=dim[4] 


#%%  Get plotting parameters

x=AC.getArray(vfnames[1],n_x)
y=AC.getArray(vfnames[2],n_y)
t=AC.getArray(vfnames[3],Nt) 
U=AC.getArray(vfnames[4],dim[0:4],Xrank=4) 
Uprime=AC.getArray(vfnames[5],dim[0:4], Xrank=4)

#%% Initialize Parameters

f_r_spatial=np.zeros((Nt,n_x-r))
L11_spatial=np.zeros((Nt))
timescale=np.zeros((Nt))
trange=np.zeros((Nt)) 
Uavg=np.zeros((Nt,3))


#%% Run Spatial_AC function for all t instances
for i,ti in enumerate(t):
    f_r_spatial[i,:], L11_spatial[i], Uavg[i,:]=AC.Spatial_AC(n_x,n_y,r,x,i,U,Uprime) #get f_r and L11 for all timepoints using Spatial_AC method
    timescale[i]=int(L11_spatial[i]/Uavg[i,0]) #calculate timescales for temporal AC method
    trange[i]=timescale[i]*5 #time range for the temporal_AC calculations. Arbitrary choice of 5x
    print "t=",ti,'s',' iteration=', i

#%% Save the Spatial_AC outputs into txtfiles for later use if want to look back and retrieve data faster
Sfnames=[]
Sfnames.append(fpath+'/'+AC.case+'/f_r_spatial.txt')
Sfnames.append(fpath+'/'+AC.case+'/L11_Spatial.txt')
Sfnames.append(fpath+'/'+AC.case+'/Uavg.txt')
Sfnames.append(fpath+'/'+AC.case+'/trange.txt')

with open(fpath+'/'+AC.case+'/Sfnames.txt','w') as f:
    f.writelines(["%s\n" % item  for item in Sfnames])


np.savetxt(Sfnames[0],f_r_spatial)
np.savetxt(Sfnames[1],L11_spatial)
np.savetxt(Sfnames[2],Uavg)
np.savetxt(Sfnames[3],trange)

#%% Make plot of f_r[-1]
plt.figure()
plt.plot(x[0:n_x-r]-x[0],f_r_spatial[-1,:])
plt.plot(x[0:n_x-r]-x[0],np.zeros(n_x-r))
plt.ylabel('f(r)')
plt.xlabel('r[m]') #plot f(r) for one time instance
plt.title('f(r)_spatial at t=41000s for r=n_x/4')
plt.savefig(fpath+'/'+AC.case+'/ILS_Plots/'+'f_r_tend.png', format="png")

#%% Make Uprime plots with auto color bar scale
os.makedirs(fpath+'/'+AC.case+'Up_plots')

for i, ti in enumerate(t):
    fig, (ax1) = plt.subplots(nrows=1)
    im = ax1.pcolormesh(x, y, Uprime[:,:,0,i], vmin=-2.5, vmax=2.5)
    plt.colorbar(im, orientation='vertical')
    ax1.set_aspect('equal')
    plt.savefig(fpath+'/'+AC.case+'/Up_plots/Uprime_field' +  str(ti) +'.png', format="png") #save all the field data in a directory
    print 't=',ti
    
#%%fft(U(x)) plots

Uyavg=np.zeros((Nt,n_x))
Uypavg=np.zeros((Nt,n_x))
#%%
for ti in range(0,Nt):
    for j in range(0,n_x):
        for i in range(0,n_y):
            Uyavg[ti,j]+=U[i,j,0,ti]/n_y
            Uypavg[ti,j]+=Uprime[i,j,0,ti]/n_y
#%%
for ti in range(0,Nt,Nt/10):
    labelstr = 't= %.1f s' % (ti)
    #else: labelstr = ''
    r = float(ti)/(Nt-1) # plot red -> blue
    Utrans=np.fft.fft(Uyavg[ti,:])
    freq=np.fft.fftfreq(len(Utrans))
    if ti==0 or ti==Nt-Nt/10:
        plt.loglog(freq,abs(Utrans),color=(1-r,0,r),label=labelstr)
    else:
        plt.loglog(freq,abs(Utrans),color=(1-r,0,r),label=None)
    plt.title('FFT of y-averaged u(x) for every 100s')
    plt.legend()
    plt.savefig(fpath+'/'+AC.case+'/ILS_Plots'+'/fft_uyavg.png', format="png")
    plt.savefig(fpath+'/'+AC.case+'/ILS_Plots'+'.png', format="png")
#%%    
for ti in range(0,Nt,Nt/10):
    labelstr = 't= %.1f s' % (ti)
    #else: labelstr = ''
    r = float(ti)/(Nt-1) # plot red -> blue
    Uptrans=np.fft.fft(Uypavg[ti,:])
    freq=np.fft.fftfreq(len(Uptrans))
    if ti==0 or ti==Nt-Nt/10:
        plt.loglog(freq,abs(Utrans),color=(1-r,0,r),label=labelstr)
    else:
        plt.loglog(freq,abs(Utrans),color=(1-r,0,r),label=None)
    plt.ylim([10**-4,10**4])
    plt.title('FFT of y-averaged uprime(x) for every 100s')
    plt.legend()
    plt.savefig(fpath+'/'+AC.case+'/ILS_Plots'+'/fft_uypavg.png', format="png")







    
