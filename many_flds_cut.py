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


t0=1
tf=63

mybase = "/home/u21/davidrball/david/tristan_KSI/"
plebase = "/home/u21/davidrball/david/tristan_KSI/ple_runs/"

mypost = "/output/flds.tot."
plepost = "/flds.tot."
#bangle0_base = mybase + 'sig2_test/output/flds.tot.'

bangle0_base = plebase + "sig2_bangle0"+plepost
bangle1_base = plebase + "sig2_bangle25"+plepost
bangle2_base = plebase + "sig2_bangle10"+plepost

#just for temporary comparison purposes
#bangle0_base = bangle200_base


def cutrot(mybase):
    myf = h5py.File(mybase,'r')
    dens = myf['dens'][0,:,:]
    ken = myf['ken'][0,:,:]
    bx = myf['bx'][0,:,:]
    by = myf['by'][0,:,:]
    bz = myf['bz'][0,:,:]

    ylen, xlen = np.shape(dens)
    x0 = int(xlen/4.)
    #xscan = int(ylen/2.)
    xscan=60
    xlow = x0-xscan
    xup = x0+xscan
    yext = int(xscan*istep/c_omp)
    xext = int(ylen*istep/c_omp/2)
    extlist = [-xext,xext,-yext,yext]
    myf.close()
    dens = dens[:,xlow:xup]
    ken = ken[:,xlow:xup]
    b2 =by[:,xlow:xup]**2 +bx[:,xlow:xup]**2 + bz[:,xlow:xup]**2
    return np.rot90(dens,3), np.rot90(ken,3), np.rot90(b2,3),extlist

for t in range(t0, tf):
    print(t)
    tstr = "%03d"%t
    base1=bangle0_base+tstr
    base2=bangle1_base+tstr
    base3=bangle2_base+tstr

    dens1, ken1, b1,extlist1 = cutrot(base1)
    dens2, ken2, b2,extlist2 = cutrot(base2)
    dens3, ken3, b3,extlist3 = cutrot(base3)

    logscale=0

    if t==1:
        bmax1 = np.max(b1)
        bmax2 = np.max(b2)
        bmax3 = np.max(b3)
        bmin1=np.min(b1)
        bmin2=np.min(b2)
        bmin3=np.min(b3)
        if logscale==1:
            bmax1=np.log10(bmax1)
            bmax2=np.log10(bmax2)
            bmax3=np.log10(bmax3)
            bmin1=np.log10(bmin1)
            bmin2=np.log10(bmin2)
            bmin3=np.log10(bmin3)

    if logscale==1:
        b1=np.log10(b1)
        b2=np.log10(b2)
        b3=np.log10(b3)
    fig, (ax0, ax1, ax2) = plt.subplots(1,3)
    
    plt.set_cmap('viridis')
    im0=ax0.imshow(b1,origin='lower',extent=extlist1,vmin=bmin1,vmax=bmax1/1.5)
    im1=ax1.imshow(b2,origin='lower',extent=extlist1,vmin=bmin2,vmax=bmax2/1.5)
    im2=ax2.imshow(b3,origin='lower',extent=extlist1,vmin=bmin3,vmax=bmax3/1.5)
    ax0.set_title('$\phi=0$')
    ax1.set_title('$\phi=\pi/25$')
    ax2.set_title('$\phi=\pi/10$ ($\phi>\phi_c$)')
    ax1.set_yticks([])
    ax2.set_yticks([])    
    ax0.set_xlabel('$x \; (c/\omega_{p})$')
    ax0.set_ylabel('$z \; (c/\omega_{p})$')
    ax1.set_xlabel('$x \; (c/\omega_{p})$')
    ax2.set_xlabel('$x \; (c/\omega_{p})$')

    plt.savefig('compare/phicrit_B'+tstr+'.png',dpi=300,bbox_inches='tight')
    plt.close()
    '''
    fig, (ax0,ax1,ax2)=plt.subplots(1,3)
    im0=ax0.imshow(dens1,origin='lower',cmap='Greys',extent=extlist1,vmin=5,vmax=20)
    im1=ax1.imshow(dens2,origin='lower',cmap='Greys',extent=extlist2,vmin=5,vmax=20)
    im2=ax2.imshow(dens3,origin='lower',cmap='Greys',extent=extlist3,vmin=5,vmax=20)

    plt.savefig('compare/ple_dens'+tstr+'.png',dpi=300,bbox_inches='tight')
    '''
