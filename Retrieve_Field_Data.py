import os
import numpy as np
import PlanarSubsets as ps

#A. Reading in t sets of planar data 
dpath='/scratch/equon/inflowOutflow/neutral/8mps_shear0.2_TI10.0_5km/postProcessing/planarSubsets'
AC=ps.planar_processing(path=dpath)

#B. Initialize process variables for AC functions
Nt=AC.Nt
t=AC.t
n_x,n_y,N,x,y,U=AC.getpsSample('horzArray','U')
r=n_x/4 #rth of domain to be sampled by Spatial_AC
dimensions=[n_x,n_y,N,Nt,r]

Uprime=AC.getUprime(n_x,n_y,Nt,x,U)

fpath='/home/tambrico/MMC/ILS_Data'
os.makedirs(fpath+'/'+AC.case)
vfnames=[]
vfnames.append(fpath+'/'+AC.case+'/dimensions.txt')
vfnames.append(fpath+'/'+AC.case+'/x.txt')
vfnames.append(fpath+'/'+AC.case+'/y.txt')
vfnames.append(fpath+'/'+AC.case+'/t.txt')
vfnames.append(fpath+'/'+AC.case+'/U.txt')
vfnames.append(fpath+'/'+AC.case+'/Uprime.txt')


with open(fpath+'/vfieldfnames.txt','w') as f:
    f.writelines(["%s\n" % item  for item in vfnames])

np.savetxt(vfnames[0], dimensions)
np.savetxt(vfnames[1], x.flatten())
np.savetxt(vfnames[2], y.flatten())
np.savetxt(vfnames[3], t.flatten())
np.savetxt(vfnames[4], U.flatten())
np.savetxt(vfnames[5], Uprime.flatten())



