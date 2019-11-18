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

ppc0=4

c_omp=5
istep=6


for t in range(t0,tf):
    tstr = "%03d" % t
    mybase = "/home/u21/davidrball/david/tristan_KSI/ple_runs/sig2_delgam2_nograv_constb/flds.tot."+tstr

    #mybase = "/home/u21/davidrball/david/tristan_nograv_inst/btest/output/flds.tot."+tstr

   # mybase="/home/u21/davidrball/david/tristan_KSI/nograv_nob/output/flds.tot."+tstr
    print(t)
    myf = h5py.File(mybase,'r')
    
    dens = myf['dens'][0,:,:]
 
    ken = myf['ken'][0,:,:]

    ex = myf['ex'][0,:,:] 
    #exext = myf['exext'][0,:,:]

    ey = myf['ey'][0,:,:]
    ez = myf['ez'][0,:,:]

    bx = myf['bx'][0,:,:]
    by = myf['by'][0,:,:]
    bz = myf['bz'][0,:,:]

    vx = myf['v3x'][0,:,:]
 
 

    jx = myf['jx'][0,:,:]
    jy = myf['jy'][0,:,:]
    jz = myf['jz'][0,:,:]

    vy = myf['v3y'][0,:,:]



    vz = myf['v3z'][0,:,:]




    my0, mx0 = np.shape(bz)
    
    xlen_units = mx0*istep/c_omp #now in units of skin depths

    xarr = np.linspace(0,xlen_units,mx0)

    print(my0, mx0)


    #let's take some average quantities instead of just along one stripe
    densslice = np.mean(dens,axis=0)

    kenslice = np.mean(ken,axis=0)
    jxslice = np.mean(jx,axis=0)
    jyslice = np.mean(jy,axis=0)
    jzslice = np.mean(jz,axis=0)
    vxslice = np.mean(vx,axis=0)


    vyslice = np.mean(vy,axis=0)


    vzslice = np.mean(vz,axis=0)

    bxslice = np.mean(bx,axis=0)
    byslice = np.mean(by,axis=0)
    bzslice = np.mean(bz,axis=0)
    exslice = np.mean(ex,axis=0)
    eyslice = np.mean(ey,axis=0)
    ezslice = np.mean(ez,axis=0)
    #exextslice = np.mean(exext,axis=0)
    fig, (ax0,ax1,ax2,ax3,ax4,ax5) = plt.subplots(6,1,sharex=True)
    plt.subplots_adjust(hspace=.4)
    xcol="Blue"
    ycol="Red"
    zcol="Green"

    ax0.set_title('$B_{x \; y \; z}$')
    ax0.plot(xarr,bxslice,color=xcol,label='x')
    ax0.plot(xarr,byslice,color=ycol,label='y')
    ax0.plot(xarr,bzslice,color=zcol,label='z')
    #ax0.set_ylim(-.001,.003)
    ax0.legend(frameon=False,prop={'size':12},bbox_to_anchor=(1,1))

    ax1.set_title('$E_{x \; y \; z}$')
    ax1.plot(xarr,exslice,color=xcol)
    ax1.plot(xarr,eyslice,color=ycol)
    ax1.plot(xarr,ezslice,color=zcol)
    #ax1.plot(exextslice,color="Black",label="$g_{x}$")
    #ax1.set_ylim(-.0015,.0015)


    ax1.legend(frameon=False,prop={'size':10},loc='lower right')
    
    ax2.set_title('$v_{x \; y \; z}$')
    #ax2.plot(vxislice,color=xcol)
    #ax2.plot(vxeslice,color=xcol,linestyle='dashed')
    #ax2.plot(vyislice,color=ycol)
    ##ax2.plot(vyeslice,color=ycol,linestyle='dashed')
    #ax2.plot(vzislice,color=zcol)
    #ax2.plot(vzeslice,color=zcol,linestyle='dashed')
    #ax2.plot(vxslice,color=xcol)
    ax2.plot(xarr,vyslice,color=ycol)
    ax2.plot(xarr,vzslice,color=zcol)
    ax2.plot(xarr,vxslice,color=xcol)
    #ax2.set_ylim(-.1,.1)


    #ax3.set_title('$J_{x \; y \; z}$')
    #ax3.plot(jxslice,color=xcol)
    #ax3.plot(-jyslice,color=ycol)
    #ax3.plot(jzslice,color=zcol)
    #ax3.set_ylim(-.00015,.00015)
    ax3.plot(xarr,kenslice,color="Red")
    ax3.set_title('Total Energy (bulk kinetic + thermal)')
    #ax3.set_ylim(6,13)
    ax4.set_title('Density (ppc/ppc0)')
    #ax4.plot(densislice,color="Black")
    #ax4.plot(denseslice,color="Black",linestyle='--')
    ax4.plot(xarr,densslice/ppc0,color="Black")
    ax4.set_ylim(.5,6)
    #ax4.set_xlabel('x (output units)')
    


    ax5.set_title('Currents')
    ax5.plot(xarr,-jxslice,color=xcol)
    ax5.plot(xarr,-jyslice,color=ycol)
    ax5.plot(xarr,-jzslice,color=zcol)
    #ax5.set_ylim(-.0001,.0001)
    ax5.set_xlabel('$x \;  (c/\omega_{p}$)')

    ax0.set_xlim(0,xlen_units)
    ax1.set_xlim(0,xlen_units)
    ax2.set_xlim(0,xlen_units)
    ax3.set_xlim(0,xlen_units)
    ax4.set_xlim(0,xlen_units)
    ax5.set_xlim(0,xlen_units)

    plt.savefig('ple_sig2_delgam2_nograv_constb/int500_diag'+tstr+'.png',bbox_inches='tight',dpi=200)
    
