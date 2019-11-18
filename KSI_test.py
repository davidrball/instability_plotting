import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import h5py 
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.family'] = 'STIXGeneral'


t0=23
tf=150

istep=6
c_omp=5


#normalize to initial values
#mybase = "/home/u21/davidrball/david/tristan_KSI/nograv_nob/output/flds.tot."

mybase = "/home/u21/davidrball/david/tristan_nograv_inst/btest/output/flds.tot."

myf_1 = h5py.File(mybase+"001",'r')
bx=myf_1['bx'][0,:,:]
by=myf_1['by'][0,:,:]
bz=myf_1['bz'][0,:,:]
b2=bx**2+by**2+bz**2
bmin=np.min(b2)
bmax=np.max(b2)

ken=myf_1['ken'][0,:,:]
tmin=np.min(ken)
tmax=np.max(ken)

dens=myf_1['dens'][0,:,:]
densmin=np.min(dens)
densmax=np.max(dens)

myf_1.close()

for t in range(t0,tf):
    tstr = "%03d" % t
    mybase_tmp = mybase + tstr
    myf = h5py.File(mybase_tmp,'r')
    print(tstr)
    dens = myf['dens'][0,:,:]
    bz = myf['bz'][0,:,:]
    bx = myf['bx'][0,:,:]
    by = myf['by'][0,:,:]
    b2=bz**2 + by**2 + bx**2
    ken = myf['ken'][0,:,:]
    

    fig, (ax0,ax1,ax2) = plt.subplots(3,1,sharex=True)
    

    ylen, xlen = np.shape(dens)
    xext=xlen*istep/c_omp
    yext=ylen*istep/c_omp
    extlist=[-xext/2.,xext/2.,-yext/2.,yext/2.]


    x0 = -xext/4.
    x1 = xext/4.
    
    x0arr = x0*np.ones(100)
    x1arr = x1*np.ones(100)
    yarr = np.linspace(-yext/2.,yext/2.,100)
    im0=ax0.imshow(dens,cmap='viridis',origin='lower',vmin=densmin,vmax=densmax,extent=extlist)
    im1=ax1.imshow(b2, cmap = 'Blues',origin='lower',vmin=bmin,vmax=bmax,extent=extlist)
    im2=ax2.imshow(ken, cmap='hot',origin='lower',vmin=tmin,vmax=tmax,extent=extlist)

    ax1.plot(x0arr,yarr,color='red',linestyle='dashed')
    ax1.plot(x1arr,yarr,color='red',linestyle='dashed')

    ax0.set_title('Density',size=12)
    ax1.set_title('$B^{2}$',size=12)
    ax2.set_title('Total Energy',size=12)
    divider0 = make_axes_locatable(ax0)
    cax = divider0.append_axes("right",size="5%",pad=0.05)
    divider1 = make_axes_locatable(ax1)
    cax1=divider1.append_axes("right",size="5%",pad=0.05)
    divider2 = make_axes_locatable(ax2)
    cax2=divider2.append_axes("right",size="5%",pad=0.05)
    
    cb0 = fig.colorbar(im0,ax=ax0,cax=cax)
    cb1 = fig.colorbar(im1,ax=ax1,cax=cax1)
    cb2 = fig.colorbar(im2,ax=ax2,cax=cax2)
    
    cb0.ax.tick_params(labelsize=8)
    cb1.ax.tick_params(labelsize=8)
    cb2.ax.tick_params(labelsize=8)

    ax2.set_xlabel('$x \; (c/\omega_p)$')
    ax0.set_ylabel('$y \; (c/\omega_p)$')
    ax1.set_ylabel('$y \; (c/\omega_p)$')
    ax2.set_ylabel('$y \; (c/\omega_p)$')
    plt.savefig('nograv_b/testing_KSI'+tstr+'.png',bbox_inches='tight',dpi=300)


    plt.close()
