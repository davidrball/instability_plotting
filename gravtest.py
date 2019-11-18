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


t_list = []
vavg_list=[]

for t in range(t0,tf):
    tstr = "%03d" % t
    mybase = "/home/u21/davidrball/david/tristan_KSI/gravtest_hot_noB/output/flds.tot."+tstr

    #mybase = "/rigel/home/ls3326/home_david/tristan_KSI/smooth_offset_current_test/output/flds.tot."+tstr

    print(t)
    myf = h5py.File(mybase,'r')
    
    dens = myf['dens'][0,:,:]
    densi = myf['densi'][0,:,:]
    dense = dens - densi



    ex = myf['ex'][0,:,:] 
    ey = myf['ey'][0,:,:]

    vx = myf['v3x'][0,:,:]
    vxi = myf['v3xi'][0,:,:]
    vxe = 2*vx - vxi

    jx = myf['jx'][0,:,:]
    jy = myf['jy'][0,:,:]
    

    vy = myf['v3y'][0,:,:]
    vyi = myf['v3yi'][0,:,:]
    vye = 2*vy - vyi

    my0, mx0 = np.shape(vy)
    
    print(my0, mx0)


    #let's take some average quantities instead of just along one stripe
    densslice = np.mean(dens,axis=0)
    densislice = np.mean(densi,axis=0)
    denseslice = np.mean(dense,axis=0)

    jxslice = np.mean(jx,axis=0)
    jyslice = np.mean(jy,axis=0)
    

    vxslice = np.mean(vx,axis=0)
    vxislice = np.mean(vxi,axis=0)
    vxeslice = np.mean(vxe,axis=0)
    
    vyslice = np.mean(vy,axis=0)
    vyislice = np.mean(vyi,axis=0)
    vyeslice = np.mean(vye,axis=0)

    

    #bxslice = np.mean(bx,axis=0)
    


    exslice = np.mean(ex,axis=0)
    eyslice = np.mean(ey,axis=0)
    #ezslice = np.mean(ez,axis=0)
    #exextslice = np.mean(exext,axis=0)

    fig, (ax1,ax2,ax5) = plt.subplots(3,1,sharex=True)
    plt.subplots_adjust(hspace=.4)

    xcol="Blue"
    ycol="Red"
    zcol="Green"

    ax1.set_title('$E_{x \; y \; z}$')
    ax1.plot(exslice,color=xcol)
    ax1.plot(eyslice,color=ycol)
    #ax1.plot(ezslice,color=zcol)
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
    ax2.plot(vxslice,color=xcol)
    ax2.plot(vyslice,color=ycol)
   
    
    ax5.set_title('Currents')
    ax5.plot(-jxslice,color=xcol)
    ax5.plot(-jyslice,color=ycol)
    
    #ax5.set_ylim(-.0001,.0001)
    ax5.set_xlabel('X (4 $c/\omega_{p}$)')
    plt.savefig('gravtest_hot/int500_diag'+tstr+'.png',bbox_inches='tight',dpi=200)
    
