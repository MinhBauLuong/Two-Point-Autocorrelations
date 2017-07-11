import numpy as np
import PlanarSubsets as ps
import matplotlib.pyplot as plt

dpath='/scratch/equon/inflowOutflow/neutral/8mps_shear0.2_TI10.0_5km/postProcessing/planarSubsets'
fpath='/home/tambrico/MMC/ILS_Data'
AC=ps.planar_processing(path=dpath)

#%% Retrieve Parameters from txtfiles

vfnames=AC.getfnames(fpath+'/'+AC.case+'/vfnames.txt')
Sfnames=AC.getfnames(fpath+'/'+AC.case+'/Sfnames.txt')
STfnames=AC.getfnames(fpath+'/'+AC.case+'/STfnames.txt')
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


f_r_st_dim=[n_x-r,Nt/100]
L11_spatial=AC.getArray(Sfnames[1],Nt)
f_r_st=AC.getArray(STfnames[0],f_r_st_dim,Xrank=2)
L11_st=AC.getArray(STfnames[1])
j=0

#%% Plot fr_st's
Ntwindow=Nt/10
j=0
for ti in range(0,Nt,Nt/10):
    plt.figure()
    plt.plot(x[0:n_x-r]-x[0],f_r_st[j,:])
    plt.xlabel('r[m]')
    plt.ylabel('f_r[x]')
    plt.title('ST Averaged f_r for 100s intervals, t=[%s, %s]' % (ti,ti+Ntwindow))
    plt.savefig(fpath+'/'+AC.case+'/ILS_Plots' + 'f_r(t=[%s, %s])' % (ti,ti+Ntwindow)+'.png', format="png")
    
    j+=1    

#%% Plot L11(t) and the ST Averaged L11s at 100s intervals
j=0
Ntwindow=Nt/10
Sfnames=AC.getfnames(fpath+'/'+AC.case+'/Sfnames.txt')
L11_spatial=AC.getArray(Sfnames[1],Nt)
f_r_st=AC.getArray(Sfnames[2],n_x-r)
with open(Sfnames[4], 'r') as f:
    L11_st_all=float(f.readline())
L11_st_all_plot=np.zeros((Nt))
L11_st_all_plot=L11_st_all
    
for ti in range(0,Nt,Nt/10):
    if ti>0:
       plt.plot(t[:],L11_spatial[:], label=None)
       plt.plot(t[ti],L11_st[j], 'ro',label=None)
       plt.plot(t[:],L11_st_all_plot[:],label=None)
    else:
        plt.plot(t[:],L11_spatial[:], label='L11(t)')
        plt.plot(t[ti],L11_st[j], 'ro',label='L11_st[t_window]')
        plt.plot(t[:],L11_st_all_plot[:],label='L11_st_all')
    plt.xlabel('time [s]')
    plt.ylabel('L11(t) [m]')
    plt.title("""Two-point autocorrelations of u'(x) as a function of time""" )
    plt.legend()
    plt.savefig(fpath+'/'+AC.case+'/ILS_Plots' + 'L11(t)' +'.png', format="png")
    j+=1
    
#%%Plot fr_st_all
plt.plot(x[0:n_x-r]-x[0],f_r_st[:])
plt.xlabel('r(m)')
plt.ylabel('f_r_st(r)')
plt.title('Space/time averaged f(r) vs r')
plt.savefig(fpath+'/'+AC.case+'/ILS_Plots' + 'f_r_st_all' +'.png', format="png")  

