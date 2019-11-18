import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import h5py 


plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.family']='STIXGeneral'


t0=120
tf=140

istep=6
c_omp=5

mysig=2
mybangle=50

sigstr = str(mysig)
banglestr=str(mybangle)

savedir = "bangle"+banglestr+"_sig"+sigstr

for t in range(t0,tf):
    tstr = "%03d" % t
    mybase = "/home/u21/davidrball/david/tristan_KSI/sig"+sigstr+"_bangle"+banglestr+"/output/spect." + tstr
    myf = h5py.File(mybase,'r')
    print(tstr)
    tspec=myf['spece']
    print('shape of spece ',np.shape(tspec))
    
    
    gam =myf['gamma']
    
    gamnum,spacenum = np.shape(tspec)

    x0_ind = int(spacenum/4.)
    xscan=3
    x0low=x0_ind-xscan
    x0up=x0_ind+xscan

    x1_ind = 3*x0_ind

    x1low=x1_ind-xscan
    x1up=x1_ind+xscan

    specx1 = np.zeros(gamnum)
    specx0 = np.zeros(gamnum)

    cold_low=0
    cold_up=4
    hot_ind=spacenum/2
    hot_low = int(hot_ind-2)
    hot_up = int(hot_ind+2)
    
    testspec_cold = np.zeros(gamnum)
    testspec_hot = np.zeros(gamnum)

    for i in range(gamnum):
        specx0[i]+=np.sum(tspec[i,x0low:x0up])
        specx1[i]+=np.sum(tspec[i,x1low:x1up])
        
        #testspec_cold[i]+=np.sum(tspec[i,cold_low:cold_up])
        #testspec_hot[i]+=np.sum(tspec[i,hot_low:hot_up])
        


    plt.plot(gam, gam*specx0,label='Left Interface')
    plt.plot(gam,gam*specx1,label='Right Interface')
    #plt.plot(gam,testspec_cold*gam,color="Blue",label='Cold')
    #plt.plot(gam,testspec_hot*gam,color="Red",label='Hot')
    
    plt.legend(frameon=False)
    plt.xlim(.1,70)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\gamma-1$')
    plt.ylabel('$(\gamma-1)dN/d\gamma$')
    plt.savefig(savedir+"/spec"+tstr+".png",dpi=300,bbox_inches='tight')
    plt.close()

    myf.close()
