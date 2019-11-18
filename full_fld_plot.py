import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from growth_rate_func import dens_analysis, return_dens_pert, fourier_analysis
import numpy as np
import h5py 

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 16})

mybangle=50
mysig=2

bstr=str(mybangle)
sigstr=str(mysig)
t0=34
tf=101

dstripe=200
istep=6
c_omp=5


#mybase = "../../relativistic_KSI/sig1/output/flds.tot.001"
#mybase = "/home/u21/davidrball/david/tristan_KSI/sig2_bangle0/output/flds.tot."
mybase = "/home/u21/davidrball/david/tristan_KSI/"
mybase += "sig"+sigstr +"_"+"bangle"+bstr+"/output/flds.tot."
for t in range(t0,tf):
    tstr = "%03d"%t
    print(tstr)
    newbase = mybase + tstr
    myf = h5py.File(newbase,'r')
    dens = myf['dens'][0,:,:] #can cut dens here to isolate one region and speed up calculation

    ylen,xlen=np.shape(dens)

    yext = ylen*istep/c_omp
    xext = xlen*istep/c_omp

    extlist=[-xext/2.,xext/2.,-yext/2.,yext/2.]

    x0=int(xlen/4.)
    x1=x0*3
    


    fig, ax1 = plt.subplots(1)
    ax1.imshow(dens,origin='lower',vmin=4,vmax=20,cmap='viridis',extent=extlist)

    
    savestr = "bangle"+bstr+"_sig"+sigstr+"/densplot"+tstr+".png"
    plt.savefig(savestr,bbox_inches='tight',dpi=300)
    plt.close()
    
