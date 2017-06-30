import numpy as np
import matplotlib.pyplot as plt
import PlanarSubsets as ps
#import imageio

#%% 
#A. Reading in t sets of planar data 
dpath='MMC_Data/planarSubsets'
AC=ps.planar_processing(path=dpath)

#B. Initialize process variables for AC functions
Nt=AC.Nt
t=AC.t
n_x,n_y,N,x,y,U=AC.getpsSample('horzArray','U')
r=n_x/4 #rth of domain to be sampled by Spatial_AC
f_r_spatial=np.zeros((Nt,n_x-r))
L11_spatial=np.zeros((Nt))
timescale=np.zeros((Nt))
trange=np.zeros((Nt)) 
uprms=np.zeros((Nt))
Uavg=np.zeros((Nt,3))
Uprime=np.zeros((n_y,n_x,3,Nt))

#C. Call spatial_AC function
#%%
for i,ti in enumerate(t):
    f_r_spatial[i,:], L11_spatial[i], Uavg[i,:], Uprime[:,:,:,i]=AC.Spatial_AC(n_x,n_y,N,x,U,i,r) #get f_r and L11 for all timepoints using Spatial_AC method
    timescale[i]=int(L11_spatial[i]/Uavg[i,0]) #calculate timescales for temporal AC method
    trange[i]=timescale[i]*5 #time range for the temporal_AC calculations. Arbitrary choice of 5x
    print "t=",ti,'s',' iteration=', i

#%%
#D. Save All the data and process in a new script
dimensions=[n_x,n_y,N,Nt,r]
name1=dpath+'/U'+'_planarSubsets'+'.txt'
name2=dpath+'/Uprime'+'_planarSubsets'+'.txt'
name3=dpath+'/f_r_spatial''_planarSubsets'+'.txt'
name4=dpath+'/L11_spatial'+'_planarSubsets'+'.txt'
name5=dpath+'/dimensions'+'_planarSubsets'+'.txt'
name6=dpath+'/x'+'_planarSubsets'+'.txt'
name7=dpath+'/y'+'_planarSubsets'+'.txt'
name8=dpath+'/f_r_st'+'_planarSubsets'+'.txt'
name9=dpath+'/L11_st'+'_planarSubsets'+'.txt'
name10=dpath+'/t.txt'

fnames=[name1,name2,name3,name4,name5,name6,name7,name8,name9,name10]
#%%
np.savetxt(name1,U.flatten())
np.savetxt(name2,Uprime.flatten())
np.savetxt(name3,f_r_spatial.flatten())
np.savetxt(name4,L11_spatial.flatten())
dim=np.array(dimensions)
np.savetxt(name5,dim.flatten())
np.savetxt(name6,x.flatten())
np.savetxt(name7,y.flatten())
np.savetxt(name10,t.flatten())
#%%
with open(dpath+'/fnames.txt','w') as f:
    f.writelines(["%s\n" % item  for item in fnames])



#%%
#E. Call temporal AC function

  
f_r_temporal=[],[]
L11_temporal=[]
tplot=[]  

for ti in range(0, Nt):
    f_r_temporal,L11_temporal,ti=AC.temporal_AC(n_x,n_y,N,x,U,ti,trange[ti],t)
    tplot.append(t[ti])

plt.figure()
plt.plot(tplot[:],L11_temporal[:])
plt.xlabel('t[s]')
plt.ylabel('L11_t(twindow) [m]')    

#%%
#F. Call space-time averaging AC function


f_r_st=np.zeros((n_x))
L11_st=0

f_r_st, L11_st=AC.Space_time_AC(n_x,n_y,N,Nt,x,U,t,r)
np.savetxt(name8,f_r_st)
with open(name9, 'w') as f:
    f.write(str(L11_st))


#%%
#with imageio.get_writer('C:/Users/tambrico/Desktop/NWTC Internship Files/Python Assignments/MMC_Data/Planar_Subsets_Plots/Uprime.gif', mode='I') as writer:
#    filenames='C:/Users/tambrico/Desktop/NWTC Internship Files/Python Assignments/MMC_Data/Planar_Subsets_Plots'
#    for filename in filenames:
#        image = imageio.imread(filename)
#        writer.append_data(image)






 


