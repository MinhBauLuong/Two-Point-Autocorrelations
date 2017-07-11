import numpy as np
import PlanarSubsets as ps


dpath='/scratch/equon/inflowOutflow/neutral/8mps_shear0.2_TI10.0_5km/postProcessing/planarSubsets'
fpath='/home/tambrico/MMC/ILS_Data'
AC=ps.planar_processing(path=dpath)

#%% Retrieve Parameters from txtfiles

vfnames=AC.getfnames(fpath+'/'+AC.case+'/fnames.txt')
dim=AC.getDimensions(vfnames[0])
n_x=dim[0]
n_y=dim[1]
N=dim[2]
Nt=dim[3]
r=dim[4] 


#%%  Get plotting parameters

x=AC.getArray(vfnames[1],n_x)
t=AC.getArray(vfnames[3],Nt) 
U=AC.getArray(vfnames[4],dim[0:4],Xrank=4)
Uprime=AC.getArray(vfnames[5],dim[0:4], Xrank=4)

#%% Call space-time averaging AC function, have it average for every 100s

Ntwindow=int(Nt/10) #Average over Nt/10 time instances
f_r_st=np.zeros((Ntwindow/10,n_x-r))
L11_st=np.zeros((Ntwindow/10))
j=0

uprsqtavg_running=0
f_r_running=np.zeros((n_x-r))
L11_running=np.zeros((Nt))

for i in range(0,Nt,Nt/10):
    f_r_st[j,:], L11_st[j],f_r_running,L11_running,uprsqtavg_running=AC.Space_time_AC(n_x,n_y,r,Ntwindow,i,x,U,Uprime,uprsqtavg_running,f_r_running,L11_running,Nt)
    print 'i=',i
    j=j+1
    
#%%    
STfnames=[]
STfnames.append(fpath+'/'+AC.case+'/f_r_st.txt')
STfnames.append(fpath+'/'+AC.case+'/L11_st.txt')
STfnames.append(fpath+'/'+AC.case+'/f_r_st_all.txt')
STfnames.append(fpath+'/'+AC.case+'/L11_st_all.txt')

with open(fpath+'/'+AC.case+'/STfnames.txt','w') as f:
    f.writelines(["%s\n" % item  for item in STfnames])
    

np.savetxt(STfnames[0],f_r_st.flatten())
np.savetxt(STfnames[1],L11_st.flatten())
np.savetxt(STfnames[2],f_r_running.flatten())

with open(STfnames[3], 'w') as f:
    f.write(str(L11_st))



