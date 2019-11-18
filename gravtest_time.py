import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import h5py 
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 16})

t0=2
tf=120


t_list = []
vavg_list=[]
vavg2_list=[]
v_expect_list = []
v_expect_list2 = []
Sigma=0.1
dstripe=200
c=.45
interval=500

delgam_hot=10
delgam_cold=1e-4

g=Sigma*2*np.pi*c**2/dstripe

for t in range(t0,tf):
    tstr = "%03d" % t
    mybase = "/home/u21/davidrball/david/tristan_KSI/gravtest_hot_noB/output/flds.tot."+tstr

    mybase2 = "/home/u21/davidrball/david/tristan_KSI/gravtest_cold_B/output/flds.tot."+tstr

    #mybase = "/rigel/home/ls3326/home_david/tristan_KSI/smooth_offset_current_test/output/flds.tot."+tstr


    print(t)
    myf = h5py.File(mybase,'r')    
    vx = myf['v3x'][0,:,:]
    myf.close()

    vx_expect = interval*(t-1)*g/(1+4*delgam_hot)
    v_expect_list.append(vx_expect)
    meanvx = np.mean(vx)
    vavg_list.append(meanvx)
    t_list.append(t)

    myf2 = h5py.File(mybase2,'r')
    vx2=myf2['v3x'][0,:,:]
    myf2.close()
    meanvx2 = np.mean(vx2)
    vavg2_list.append(meanvx2)
    vx_expect2 = 2.5*vx_expect
    v_expect_list2.append(vx_expect2)
    


    plt.plot(t_list,vavg_list,label='$\\theta=10 \; \sigma=0$',color="Red")
    plt.plot(t_list,vavg2_list,label='$\\theta=10^{-4} \; \sigma=10$',color="Blue")


    plt.xlabel('Time')
    plt.ylabel('$\\bar{v_x}$')
    plt.legend(frameon=False)
    #plt.plot(t_list,v_expect_list,color="Red",linestyle='--')
    #plt.plot(t_list,v_expect_list2,color="Blue",linestyle='--')
    plt.savefig('gravtest_both/vx_time'+tstr+'.png',bbox_inches='tight',dpi=200)

    plt.close()
    
