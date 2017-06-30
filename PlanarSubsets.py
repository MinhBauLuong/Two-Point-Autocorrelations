#!/usr/bin/python
import sys,os
import numpy as np
import scipy as sp

class planar_processing:
    
    """ For collecting subsets of planar velocity data at t timepoints
    """

    sampleExt1 = 'mesh' #file that contains the coordinates of the sample
    sampleExt2 = 'U' #file that contains the velocity data at each time instance
    sampleExt3 = 'dat'

    def __init__(self,path='.'): #sets local directory as the path by default
        
        #self.ds=datastruct
        
        #if self.ds=='ps':
        """ Sets timeNames and sampleNames with the time directory and sample names, respectively (just the name of the file folder).
        Collects the time instances into an array t
        """
        self.path = path
        self.t = []
        self.timeNames = []
        self.sampleNames = []
            
        # get list of times from directory names
        try:
            dirs = next(os.walk(path))[1] #searches all branches of subdirectories within the specified path
        except StopIteration:
            print 'Path',path,'not found'
            return
        times = []
        timeNames = []
        for dname in dirs:
            try:
                tval = float(dname) #file names e.g. 40000 will be converted to a float and appended to the times array
                times.append( tval )
                timeNames.append( dname ) #stores the actual file name
            except ValueError: pass
        self.t = np.array(times)
        self.Nt = len(self.t)
        if self.Nt==0: 
            print 'No time directories found in', path
            return

# sort based on times
        order = np.argsort(self.t)
        self.t = self.t[order]
        self.timeNames = [ timeNames[i] for i in order ] #sorts the times and timenames from t0 to tend
        dt = np.diff(self.t)
        if np.max(dt)-np.min(dt) > 1e-14:
            print 'Warning: gaps detected in sampling times'
            self.dt = dt 
        else:
            self.dt = dt[0] #just makes sure the time step isn't too small

# gather sample names
        path0 = os.path.join(path,timeNames[0])
        flist = [ f for f in os.listdir(path0) if os.path.isfile(os.path.join(path0, f)) ] #puts all of the file names into a list
        for f in flist:
            fsplit = f.split('.')
            name = '.'.join( fsplit[:-1] )
            self.sampleNames.append( name )
            
            
        
    def __repr__(self):
        s=''
        s = 'Read times, t = %s\n' % (self.t) \
        + 'Sample names :'
        for name in sorted( self.sampleNames ):
            s += '\n  %s' % name
        return s #Displays which file folder is being sifted through
        
#        else: 
#            self.path=path
#            self.t=[]

    
    
    def getpsSample(self,name,field,verbose=True):
        
        #if self.ds!='ps': print("Warning: planar subset getSample function used")
        
        
        """ 1) Returns coordinates of the planes being analyzed. Returns a vector of shape (n_x,n_y,Nt), for an n_x*n_y plane at time t for Nt time points
        2) Returns a 4d array of velocity data of size (n_x*n_y*3*Nt) where n_x*n_y is a plane, 0-->3 are the velocity components, and Nt is the sampled field at each time instance t
        """
    #Retrieve the U.mesh filename
        
        found = False #Example file name: horzArray_U.mesh: name='horzArray', field='U'
        suffix1 = '_' + field
        for f in self.sampleNames:
            if f.startswith(name) and f.endswith(suffix1): #searches the sample names for the specified file
                found = True
                break
        if not found:
            print 'Sample',name,'with field',field,'not found'
            return

        fname1 = f + '.' + self.sampleExt1 #returns the full file name, 'horzArray_U.mesh'
        

    #Retrieve the U.000.U filename       
        
        found = False #Example file name: horzArray_U.000.U: name='horzArray', field='U'
        suffix2 = '_' + field + '.000'
        for f in self.sampleNames:
            if f.startswith(name) and f.endswith(suffix2):
                found = True
                break
        if not found:
            print 'Sample',name,'with field',field,'not found'
            return
        
        fname2=f+'.'+self.sampleExt2
    
    #initialize R and U files to write into and read from for reduced processing time
#    
#        Rfile = self.path + os.sep + f + '.R.npy'
#        Ufile = self.path + os.sep + f + '.'+field + '.npy'
#        Nfile = self.path + os.sep + f + '.N.npy'
#        
#        try:
#            # read from pre-processed numpy data file
#            #x,U = np.load( savefile )
#            R = np.load( Rfile )
#            U = np.load( Ufile )
#            N = np.load( Nfile )
#            n=int(np.sqrt(N))
#           
#            print 'Data read from',Ufile
#
#        except IOError:

            
        testfile1 = os.path.join(self.path, self.timeNames[0], fname1)

    #Retrieve the coordinates of the mesh from the first .mesh file

        with open(testfile1, 'r') as f:
                
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
        
        with open(testfile1, 'r') as f:
    
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
        U=np.zeros((n_y,n_x,3,self.Nt))
             
            # process U.000.U files in all time directories to fill the U array
            
        for it,tdir in enumerate(self.timeNames):
                sys.stdout.write('\r  reading t= %f' % self.t[it])
                with open(os.path.join(self.path, tdir, fname2), 'r') as f:
                     for i in range(0,4):
                         f.readline()
        
                     for k in range(0,3):    
                         for j in range(0,n_y):
                             for i in range(0,n_x):
                                 U[j,i,k,it]=float(f.readline()) #stores u, v, w in a planar array corresponding to R. Scans x for fixed y, then scans y, then scans each velocity component
                                    

#        print '  saving',Ufile
#        try:
#            np.save( Rfile, R )
#            np.save( Ufile, U )
#        except IOError as err:
#            print '  warning, unable to write out npy file:',err

        return n_x,n_y,N,xcoord,ycoord,U
    
#    def gettsSample(self,name,verbose=True):
#        
#        found=False
#        f=self.path+name+'.'+self.sampleExt3
#        
#        with open(f, 'r') as f:
#            
#            
            
        
    
    
    def Spatial_AC(self,n_x,n_y,N,xcoord,U,t,r): 

#performs a spatial autocorrelation on a planar set of velocity data (for a, x-range r) at time t and computes the integral length scale of the eddies in the flow at time t

        Uavg=np.zeros(3)

        for j in range(0,n_y):
            for i in range(0,n_x):
                Uavg[0]+=U[j,i,0,t]/N
                Uavg[1]+=U[j,i,1,t]/N
                Uavg[2]+=U[j,i,2,t]/N #get average planar velocity components

        Uprime=np.zeros((n_y,n_x,3))
        Uprime[:,:,0]=U[:,:,0,t]-Uavg[0] 
        Uprime[:,:,1]=U[:,:,1,t]-Uavg[1]
        Uprime[:,:,2]=U[:,:,2,t]-Uavg[2] #get fluctuating velocity terms


        upsqavg=0

        for j in range(0,n_y):
            for i in range(0,n_x):
                upsqavg+=(Uprime[j,i,0])**2/N #calculate planar average of squared velocity for whole domain
            


   # r=n/2 #amount of domain sampled in x

        uprsqavg=0 #planar average of u for domain sized r*n (r in x, n in y)
        uprms=0 #rms fluctuating velocity of an rth of the field. Use to do divide lengthscale by to get timescale
    
        for j in range(0,n_y):
            for i in range(0,r):
                uprsqavg+=(Uprime[j,i,0])**2/(2*r*n_y) #calculate planar average of squared velocity for rth of the domain being sampled
    
        for j in range(0,n_y):
            for i in range(n_x,n_x-r):
                uprsqavg+=(Uprime[j,i,0])**2/(2*r*n_y)
    
        uprms=np.sqrt(uprsqavg)    

               
            #B. Calculating spatial autocorrelations for an rth of the domain

        f_r=np.zeros((n_x-r))

        for k in range(0,n_x-r):
            #print 'k=',k
            for j in range(0,n_y):
                #print 'j=',j
                for i in range(0,r):
                    f_r[k]+=(Uprime[j,i,0]*Uprime[j,i+k,0])/(2*r*n_y)
                    
                
        for k in range(0,n_x-r):
            for j in range(0,n_x):
                for i in range(n_x,n_x-r):
                    f_r[k]+=(Uprime[j,i,0]*Uprime[j,i-k,0])/(2*r*n_y)
    

        f_r=f_r/uprsqavg
        L11=sp.integrate.simps(f_r[0:(n_x-r)], x=xcoord[0:(n_x-r)]-xcoord[0]) #get lengthscale at time t
    
        return f_r, L11, Uavg, Uprime
    
    def temporal_AC(self,n_x,n_y,N,xcoord,U,t0,t_window,t): 
      
        uavg=np.zeros(t_window)
        vavg=np.zeros(t_window)
        wavg=np.zeros(t_window) #averaged velocities over the whole plane for each t
        
        for ti in range(t0, t0+t_window): #t0 and t_window have to be array counters corresponding to actual times
            for j in range(0,n_y):
                for i in range(0,n_x):
                    uavg[ti]+=U[j,i,0,ti]/N
                    vavg[ti]+=U[j,i,1,ti]/N
                    wavg[ti]+=U[j,i,2,ti]/N #get average planar velocity components

        Uprime=np.zeros((n_y,n_x,3,t_window))
    
        for ti in range(t0, t0+t_window):
            Uprime[:,:,0,ti]=U[:,:,0,ti]-uavg[ti] 
            Uprime[:,:,1,ti]=U[:,:,1,ti]-vavg[ti]
            Uprime[:,:,2,ti]=U[:,:,2,ti]-wavg[ti] #get fluctuating velocity terms
        
        upt0sqavg=0
    
        for j in range(0,n_y):
            for i in range(0,n_x):
                upt0sqavg+=Uprime[j,i,0,0]**2/N #get planar average of fluctuating velocities at t0
    
        f_t=np.zeros(t_window)
    
        for k in range(t0,t0+t_window):
            for j in range(0,n_y):
                for i in range(0,n_x):
                    f_t[k]+=(Uprime[j,i,0,k]*Uprime[j,i,0,k+i])/N
                
        f_t=f_t/upt0sqavg
        L11_t=sp.integrate.simps(f_t[:],t[t0:t0+t_window])
    
        return f_t, L11_t, t0+t_window
    
    def Space_time_AC(self,n_x,n_y,N,Nt,xcoord,U,t,r): #calculates average lengthscale using spatial average approach, but also averages over all timesteps
     
        Uavg=np.zeros((3,Nt))
        for ti in range(0,Nt):
            for j in range(0,n_y):
                for i in range(0,n_x):
                    Uavg[0,ti]+=U[j,i,0,ti]/N
                    Uavg[1,ti]+=U[j,i,1,ti]/N
                    Uavg[2,ti]+=U[j,i,2,ti]/N #get average planar velocity components at each ti

        Uprime=np.zeros((n_y,n_x,3,Nt))
    
        for ti in range(0,Nt):
            Uprime[:,:,0,ti]=U[:,:,0,ti]-Uavg[0,ti] 
            Uprime[:,:,1,ti]=U[:,:,1,ti]-Uavg[1,ti]
            Uprime[:,:,2,ti]=U[:,:,2,ti]-Uavg[2,ti] #get fluctuating velocity terms


        upsqavg=0
    
        for ti in range(0,Nt):
            for j in range(0,n_y):
                for i in range(0,n_x):
                    upsqavg+=(Uprime[j,i,0,ti])**2/(N*Nt) #calculate planar/time average of squared velocity for whole domain
    
        uprsqtavg=0 #planar average of u for domain sized r*n (r in x, n in y)
    

        for ti in range(0,Nt):
            for j in range(0,n_y):
                for i in range(0,r):
                    uprsqtavg+=(Uprime[j,i,0,ti])**2/(2*r*n_y*Nt) #calculate planar average of squared velocity for left half of the domain
    
  
        for ti in range(0,Nt):
            for j in range(0,n_y):
                for i in range(n_x,n_x-r):
                    uprsqtavg+=(Uprime[j,i,0,ti])**2/(2*r*n_y*Nt)
                    
        f_r=np.zeros((n_x-r))

        for ti in range(0,Nt):
            for k in range(0,n_x-r):
                for j in range(0,n_y):
                    for i in range(0,r):
                        f_r[k]+=(Uprime[j,i,0,ti]*Uprime[j,i+k,0,ti])/(2*r*n_y*Nt)
                
        for ti in range(0,Nt):
            for k in range(0,n_x-r):
                for j in range(0,n_y):
                    for i in range(n_x,n_x-r):
                        f_r[k]+=(Uprime[j,i,0,ti]*Uprime[j,i-k,0,ti])/(2*r*n_y*Nt)
    

        f_r=f_r/uprsqtavg
        L11=sp.integrate.simps(f_r[0:(n_x-r)], x=xcoord[0:(n_x-r)]-xcoord[0])
    
        return f_r, L11




#===============================================================================
#===============================================================================
#===============================================================================
# test run
#if __name__=='__main__':
#
#    line0 = planar(path='linesTransverse')
#    print line0
#
#    x,U = line0.getSample('line_08km','U')