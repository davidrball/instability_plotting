import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import h5py 
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 9})

t0=1
tf=120



for t in range(t0,tf):
    tstr = "%03d" % t
    mybase = "/home/u21/davidrball/david/tristan_KSI/ple_runs/sig2_bangle0_nograv2/flds.tot."+tstr

    #mybase = "/rigel/home/ls3326/home_david/tristan_KSI/smooth_offset_current_test/output/flds.tot."+tstr

    print(t)
    myf = h5py.File(mybase,'r')
    
    dens = myf['dens'][0,:,:]

    ey = myf['ey'][0,:,:]

    eyplus = np.roll(ey,1,axis=0)
    eyminus = np.roll(ey,-1,axis=0)

    eydiff = eyplus-eyminus

    #ez = myf['ez'][0,:,:]

    jy = myf['jy'][0,:,:]

    my0, mx0 = np.shape(jy)
    xscan=250
    x0 = int(mx0/4.)
    xlow=x0-xscan
    xup=x0+xscan
    


    #now just make cuts of only the quantities we actually want to plot
    jy = jy[:,xlow:xup]


    if t==1:
        jy0 = jy


    dens = dens[:,xlow:xup]

    ey = ey[:,xlow:xup]


    fig, (ax1, ax2, ax3) = plt.subplots(1,3,sharey=True)
    
    im1=ax1.imshow(eydiff,origin='lower',cmap='RdBu')
    im2=ax2.imshow(ey,origin='lower',cmap='RdBu',vmin=-.0075,vmax=.0075)
    im3=ax3.imshow(dens,origin='lower',cmap='viridis',vmin=4,vmax=24)

    plt.rcParams.update({'font.size': 4})

    cbar_ax1 = fig.add_axes([.125,.97,.22,.02])
    cb1 = fig.colorbar(im1,cax=cbar_ax1,orientation='horizontal')

    cbar_ax2 = fig.add_axes([.4,.97,.22,.02])
    cb2 = fig.colorbar(im2,cax=cbar_ax2,orientation='horizontal')
    
    cbar_ax3 = fig.add_axes([.675,.97,.22,.02])
    cb3 = fig.colorbar(im3,cax=cbar_ax3, orientation='horizontal')

    plt.rcParams.update({'font.size': 9})


    ax1.set_title('$10000j_y$')
    ax2.set_title('$E_y$')
    ax3.set_title('Density')

    plt.savefig("ple_bangle0_sig2_nograv/eydiff"+tstr+".png",dpi=300,bbox_inches='tight')
    plt.close()
