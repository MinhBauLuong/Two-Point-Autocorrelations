import numpy as np
import matplotlib.pyplot as plt
import PlanarSubsets as ps

cpath='C:/Users/tambrico/Desktop/NWTC_Internship_Files/Python_Assignments/'
dpath='C:/Users/tambrico/Desktop/NWTC_Internship_Files/Python_Assignments/MMC_Data/planarSubsets/'
AC=ps.planar_processing(path=dpath)
#%%
fnames=AC.getfnames(dpath+'fnames.txt')
dim=AC.getDimensions(cpath+fnames[4])
n_x=dim[0]
n_y=dim[1]
N=dim[2]
Nt=dim[3]
r=dim[4] #get sampling dimensions


#%%  Plotting L11(t)

x=AC.getArray(cpath+fnames[5],n_x)
y=AC.getArray(cpath+fnames[6],n_y)
t=AC.getArray(cpath+fnames[9],Nt) #get plotting parameters

#%%
f_r_spatial_dim=[n_x-r,Nt]
f_r_spatial=AC.getArray(fnames[2],f_r_spatial_dim,Xrank=2)

L11_spatial=AC.getArray(fnames[3],Nt)
L11_st=AC.getArray(fnames[8],1)
L11_st_plot=np.zeros(Nt)
L11_st_plot[:]=L11_st

plt.figure()
plt.plot(x[0:n_x-r]-x[0],f_r_spatial[-1,:])
plt.plot(x[0:n_x-r]-x[0],np.zeros(n_x-r))
plt.ylabel('f(r)')
plt.xlabel('r[m]') #plot f(r) for one time instance
plt.title('f(r)_spatial at t=41000s for r=n_x/4')
plt.savefig('C:/Users/tambrico/Desktop/NWTC_Internship_Files/Python_Assignments/MMC_Data/Planar_Subsets_Plots_2/' + 'fr_tend' +'.png', format="png")

print 'L_11[t=41000s] =', L11_spatial[-1], 'm'

plt.figure()
plt.plot(t[:],L11_spatial[:], label='L11_S(t)')
plt.plot(t[:],L11_st_plot, label='L11_st')
plt.xlabel('t[s]')
plt.ylabel('L11_S(t) [m]') #plot integral length scale for all time instances
plt.title('L11 vs t for r=n_x/4')
plt.savefig('C:/Users/tambrico/Desktop/NWTC_Internship_Files/Python_Assignments/MMC_Data/Planar_Subsets_Plots_2/' + 'L11s(t)' +'.png', format="png")
#%%

# Plot the U' field for all t and save in the directory Planar_Subsets_Plots_2/Uprime_field

Uprime=AC.getArray(cpath+fnames[1],dim[0:4],Xrank=4)
#%%
U=AC.getArray(fnames[0],dim[0:4],Xrank=4)
#%%
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
    Utrans=np.fft.fft(Uyavg[ti,:])
    freq=np.fft.fftfreq(len(Utrans))
    plt.plot(freq,abs(Utrans))
    plt.loglog(freq,abs(Utrans))
    plt.title('FFT of y-averaged u(x) for every 100s')
    plt.savefig('C:/Users/tambrico/Desktop/NWTC_Internship_Files/Python_Assignments/MMC_Data/Planar_Subsets_Plots_2/' + 'uavg(twindow)' +'.png', format="png")
    
for ti in range(0,Nt,Nt/10):
    Uptrans=np.fft.fft(Uypavg[ti,:])
    freq=np.fft.fftfreq(len(Uptrans))
    plt.plot(freq,abs(Uptrans))
    plt.loglog(freq,abs(Uptrans))
    plt.ylim([10**-4,10**4])
    plt.title('FFT of y-averaged uprime(x) for every 100s')
    plt.savefig('C:/Users/tambrico/Desktop/NWTC_Internship_Files/Python_Assignments/MMC_Data/Planar_Subsets_Plots_2/' + 'upavg(twindow)' +'.png', format="png")

#%%
f_r_st=AC.getArray(fnames[7],dim[0]-dim[-1])
plt.plot(x[0:n_x-r]-x[0],f_r_st[:])
plt.xlabel('r(m)')
plt.ylabel('f_r_st(r)')
plt.title('Space/time averaged f(r) vs r')
plt.savefig('C:/Users/tambrico/Desktop/NWTC_Internship_Files/Python_Assignments/MMC_Data/Planar_Subsets_Plots_2/' + 'f_r_st' +'.png', format="png")    

#%%
#for j in range(0,950,Nt/20):
#    for i, ti in enumerate(t[j:j+50]):
#        fig, (ax1) = plt.subplots(nrows=1)
#        im = ax1.pcolormesh(x, y, Uprime[:,:,0,i], cmap='coolwarm',vmin=np.min(Uprime[:,:,0,j:j+50]),vmax=np.max(Uprime[:,:,0,j:j+50]))
#        plt.colorbar(im, orientation='vertical')
#        ax1.set_aspect('equal')
#        plt.xlabel('x [m]')
#        plt.ylabel('y [m]')
#        plt.savefig('C:/Users/tambrico/Desktop/NWTC_Internship_Files/Python_Assignments/MMC_Data/Planar_Subsets_Plots_2/Uprime_field/' +  str(ti) +'.png', format="png") #save all the field data in a directory
#        print 't=',ti
#        break
#%% 

for i, ti in enumerate(t):
    fig, (ax1) = plt.subplots(nrows=1)
    im = ax1.pcolormesh(x, y, Uprime[:,:,0,i], cmap='coolwarm')
    plt.colorbar(im, orientation='vertical')
    ax1.set_aspect('equal')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.savefig('C:/Users/tambrico/Desktop/NWTC_Internship_Files/Python_Assignments/MMC_Data/Planar_Subsets_Plots_2/Uprime_field/' +  str(ti) +'.png', format="png") #save all the field data in a directory
    print 't=',ti
    break




    
                            
    
        
            
        