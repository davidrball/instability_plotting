import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import h5py 
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rcParams['figure.figsize']=25,5

t0=1
tf=150

istep=6
c_omp=5


#normalize to initial values
#mybase = "/home/u21/davidrball/david/tristan_KSI/nograv_nob/output/flds.tot."

mybase1 = "/home/u21/davidrball/david/tristan_nograv_inst/btest/output/flds.tot."
mybase2 = "/home/u21/davidrball/david/tristan_KSI/ple_runs/sig2_bangle0_nograv/flds.tot."
mybase3 = "/home/u21/davidrball/david/tristan_KSI/ple_runs/sig2_bangle0/flds.tot."


refbase = mybase1+"001"
myf_ref = h5py.File(refbase,'r')
densref = myf_ref['dens']
bref = np.array(myf_ref['bx'])**2+np.array(myf_ref['by'])**2+np.array(myf_ref['bz'])**2
kenref = myf_ref['ken']
densmin,densmax = np.min(densref), np.max(densref)
kenmin,kenmax = np.min(kenref),np.max(kenref)
bmin,bmax = np.min(bref),np.max(bref)
bmin=0


def get_flds(mybase,t):
    tstr = "%03d" % t
    mybase_tmp = mybase + tstr
    myf = h5py.File(mybase_tmp,'r')
    #print(tstr)
    
    dens = myf['dens'][0,:,:]
    bz = myf['bz'][0,:,:]
    bx = myf['bx'][0,:,:]
    by = myf['by'][0,:,:]
    b2=bz**2 + by**2 + bx**2
    ken = myf['ken'][0,:,:]

    

    
    ylen, xlen = np.shape(dens)
    myf.close()
    return dens, b2, ken,ylen,xlen
 


dens1min,dens1max,dens2min,dens2max,dens3min,dens3max,ken1min,ken1max,ken2min,ken2max,ken3min,ken3max,b1min,b1max,b2min,b2max,b3min,b3max=0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0



for t in range(t0,tf):
    fig, ((ax0,ax1,ax2),(ax3,ax4,ax5),(ax6,ax7,ax8)) = plt.subplots(3,3,sharex=True,sharey=True)
    tstr = "%03d"%t
    print(tstr)

    dens1,b1,ken1,ylen1,xlen1 = get_flds(mybase1,t)
    dens2,b2,ken2,ylen2,xlen2 = get_flds(mybase2,t)
    dens3,b3,ken3,ylen3,xlen3 = get_flds(mybase3,t)

    if t==1:
        dens1min,dens1max = np.min(dens1),np.max(dens1)
        b1min,b1max = np.min(b1),np.max(b1)
        ken1min,ken1max = np.min(ken1),np.max(ken1)

        dens2min,dens2max = np.min(dens2),np.max(dens2)
        b2min,b2max = np.min(b2),np.max(b2)
        ken2min,ken2max = np.min(ken2),np.max(ken2)

        dens3min,dens3max = np.min(dens3),np.max(dens3)
        b3min,b3max = np.min(b3),np.max(b3)
        ken3min,ken3max = np.min(ken3),np.max(ken3)

        b1min = 0
        b2min = 0
        b3min=0
    xext=xlen1*istep/c_omp
    yext=ylen1*istep/c_omp
    extlist=[-xext/2.,xext/2.,-yext/2.,yext/2.]

    ax0.imshow(dens1,origin='lower',cmap='viridis',extent=extlist,vmin=dens1min,vmax=dens1max)
    ax3.imshow(b1,origin='lower',cmap='Blues',extent=extlist,vmin=b1min,vmax=b1max)
    ax6.imshow(ken1,origin='lower',cmap='hot',extent=extlist,vmin=ken1min,vmax=ken1max)
    

    ax1.imshow(dens2,origin='lower',cmap='viridis',extent=extlist,vmin=dens2min,vmax=dens2max)
    ax4.imshow(b2,origin='lower',cmap='Blues',extent=extlist,vmin=b2min,vmax=b2max)
    ax7.imshow(ken2,origin='lower',cmap='hot',extent=extlist,vmin=ken2min,vmax=ken2max)

    ax2.imshow(dens3,origin='lower',cmap='viridis',extent=extlist,vmin=dens3min,vmax=dens3max)
    ax5.imshow(b3,origin='lower',cmap='Blues',extent=extlist,vmin=b3min,vmax=b3max)
    ax8.imshow(ken3,origin='lower',cmap='hot',extent=extlist,vmin=ken3min,vmax=ken3max)

    plt.subplots_adjust(hspace=0)
    plt.subplots_adjust(wspace=0.05)

    #divider0 = make_axes_locatable(ax0)
    #cax = divider0.append_axes("right",size="5%",pad=0.05)
    #divider1 = make_axes_locatable(ax1)
    #cax1=divider1.append_axes("right",size="5%",pad=0.05)
    #divider2 = make_axes_locatable(ax2)
    #cax2=divider2.append_axes("right",size="5%",pad=0.05)
    
    #cb0 = fig.colorbar(im0,ax=ax0,cax=cax)
    #cb1 = fig.colorbar(im1,ax=ax1,cax=cax1)
    #cb2 = fig.colorbar(im2,ax=ax2,cax=cax2)
    
    #cb0.ax.tick_params(labelsize=8)
    #cb1.ax.tick_params(labelsize=8)
    #cb2.ax.tick_params(labelsize=8)

    ax6.set_xlabel('$x \; (c/\omega_p)$')
    ax7.set_xlabel('$x \; (c/\omega_p)$')
    ax8.set_xlabel('$x \; (c/\omega_p)$')
    #ax0.set_ylabel('$y \; (c/\omega_p)$')
    #ax3.set_ylabel('$y \; (c/\omega_p)$')
    #ax6.set_ylabel('$y \; (c/\omega_p)$')
    ax0.set_ylabel('Density')
    ax3.set_ylabel('$B^2$')
    ax6.set_ylabel('Temperature')

    ax0.set_title('No Gravity, Constant B')
    ax1.set_title('No Gravity, Magnetically Confined')
    ax2.set_title('Gravity')
    plt.savefig('compare_flds/compare_test'+tstr+'.png',bbox_inches='tight',dpi=300)


    plt.close()
