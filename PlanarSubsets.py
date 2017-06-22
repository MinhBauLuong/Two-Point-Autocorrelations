#!/usr/bin/python
import sys,os
import numpy as np

class planar:
    """ For post-processing subsets of planar velocity data at t timepoints
    """

    sampleExt1 = 'mesh' #file that contains the coordinates of the sample
    sampleExt2 = 'U' #file that contains the velocity data at each time instance

    def __init__(self,path='.'): #sets local directory as the path by default
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
        s = 'Read times, t = %s\n' % (self.t) \
          + 'Sample names :'
        for name in sorted( self.sampleNames ):
            s += '\n  %s' % name
        return s #Displays which file folder is being sifted through

    def getSample(self,name,field,verbose=True):
        
        """ 1) Returns coordinates of the planes being analyzed. Returns a vector of shape (n,n,Nt), for an nxn plane at time t for Nt time points
        2) Returns a 4d array of velocity data of size (nxnx3xNt) where nxn is a plane, 0-->3 are the velocity components, and Nt is the sampled field at each time instance t
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
    
                for i in range(1, 9):
                    f.readline()
                
                N=int(f.readline())
                n=int(np.sqrt(N))
                R=np.zeros((n,n)) #xy coordinate array
                U=np.zeros((n,n,3,self.Nt)) #velocity components array
    
                for i in range(0, n): #reads in the first 99 mesh coordinates. Since its a square, the points are fed into both the first row and the first col and then the col is repeated over
                    R[0,i]=float(f.readline())
                    if i==0: continue
                    R[i,:]=R[0,i]
             
            # process U.000.U files in all time directories to fill the U array
            
        for it,tdir in enumerate(self.timeNames):
                sys.stdout.write('\r  reading t= %f' % self.t[it])
                with open(os.path.join(self.path, tdir, fname2), 'r') as f:
                     for i in range(0,4):
                         f.readline()
        
                     for k in range(0,3):    
                         for j in range(0,n):
                             for i in range(0,n):
                                 U[j,i,k,it]=float(f.readline()) #stores u, v, w in a planar array corresponding to R. Scans x for fixed y, then scans y, then scans each velocity component
                                    

#        print '  saving',Ufile
#        try:
#            np.save( Rfile, R )
#            np.save( Ufile, U )
#        except IOError as err:
#            print '  warning, unable to write out npy file:',err

        return n,N,R,U 


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