import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import h5py 

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})


#plt.register_cmap(name='viridis',cmap=cmaps.viridis)

istep=6
c_omp=5


t0=100
tf=155

#myfstart = h5py.File("~/home_david/tristan_KSI/smooth_offset_current_test/int1_out/flds.tot.001")


for t in range(t0,tf):
    tstr = "%03d" % t
    #mybase = "../../tristan_KSI/offset_test/int10_out/flds.tot." + tstr
    #mybase = "/rigel/home/ls3326/tigress/frb_reconnection/shear/kr/test_grav/output/flds.tot."+tstr
    print(t)
    #mybase = "../../lorenzo_setup/reg_grav/dstripe25_sigext.1/output/flds.tot."+tstr
    mybase = "/home/u21/davidrball/david/tristan_KSI/sig2_bangle50/output/flds.tot."+tstr
    myf = h5py.File(mybase,'r')
    
    dens = myf['dens'][0,:,:]
    densi = myf['densi'][0,:,:]
    dense = dens - densi
    ken = myf['ken'][0,:,:]
    bx = myf['bx'][0,:,:]
    by = myf['by'][0,:,:]
    bz = myf['bz'][0,:,:]
    ex = myf['ex'][0,:,:] 
    ey = myf['ey'][0,:,:]    
    vx = myf['v3x'][0,:,:]
    vy = myf['v3y'][0,:,:]
    vz = myf['v3z'][0,:,:]

    #ex_ext = myf['exext'][0,:,:]

    jy = myf['jy'][0,:,:]
    

    
    v2 = vx**2 + vy**2 + vz**2
    
    ylen, xlen = np.shape(dens)
    x0 = int(xlen/4.)
    #xscan = int(ylen/2.)
    xscan=125
    xlow = x0-xscan
    xup = x0+xscan
    
    yext = int(xscan*istep/c_omp)
    xext = int(ylen*istep/c_omp/2)
    extlist = [-xext,xext,-yext,yext]

    myf.close()

    fig, (ax0,ax1,ax2) = plt.subplots(1,3)

    dens = dens[:,xlow:xup]
    ken = ken[:,xlow:xup]
    b2 =by[:,xlow:xup]**2 +bx[:,xlow:xup]**2 + bz[:,xlow:xup]**2
    
    if t==t0:
        densmin=np.min(dens)
        kenmin=np.min(ken)
        bmin=np.min(b2)
        densmax=np.max(dens)
        kenmax=np.max(ken)
        bmax=np.max(b2)

    
    im0=ax0.imshow(np.rot90(dens,3),origin='lower',cmap='Greys',extent=extlist,vmin=densmin,vmax=densmax)
    im1=ax1.imshow(np.rot90(ken,3),origin='lower',cmap='hot',extent=extlist,vmin=kenmin,vmax=kenmax)
    im2=ax2.imshow(np.rot90(b2,3),origin='lower',cmap='Blues',extent=extlist,vmin=bmin,vmax=bmax)

    #im0=ax0.imshow(dens,origin='lower',cmap='summer',vmin=3.5,vmax=4.5)
    #im1=ax1.imshow(bz,origin='lower',cmap='RdBu',vmin=-.035,vmax=.035)
    #im2=ax2.imshow(ex,origin='lower',cmap='RdBu',vmin=-.0001,vmax=.0001)
    
    ax0.set_title('Density')
    ax1.set_title('Temperature')
    ax2.set_title('$B^{2}$')
    '''
    divider0 = make_axes_locatable(ax0)
    cax = divider0.append_axes('right',size='5%',pad=.05)
    fig.colorbar(im0,cax=cax,orientation='vertical')

    divider1 = make_axes_locatable(ax1)
    cax = divider1.append_axes('right',size='5%',pad=.05)
    fig.colorbar(im1,cax=cax,orientation='vertical')

    divider2 = make_axes_locatable(ax2)
    cax = divider2.append_axes('right',size='5%',pad=.05)
    fig.colorbar(im2,cax=cax,orientation='vertical')

    '''

    
    ax1.set_yticks([])
    ax2.set_yticks([])
    
    ax0.set_xlabel('$x \; (c/\omega_{p})$')
    ax0.set_ylabel('$z \; (c/\omega_{p})$')
    ax1.set_xlabel('$x \; (c/\omega_{p})$')
    ax2.set_xlabel('$x \; (c/\omega_{p})$')
    #plt.savefig('dstripe25_sigext.1/cut_flds_'+tstr+'.png')
    #plt.savefig('bangle50_sig2/cutflds_'+tstr+'.png',dpi=300,bbox_inches='tight')
    plt.savefig('hubble_app_fig'+tstr+'.png',dpi=300,bbox_inches='tight')
    plt.close()

