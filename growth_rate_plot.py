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


#plt.register_cmap(name='viridis',cmap=cmaps.viridis)



t0=1
tf=100

#myfstart = h5py.File("~/home_david/tristan_KSI/smooth_offset_current_test/int1_out/flds.tot.001")


for t in range(t0,tf):
    tstr = "%03d" % t
    #mybase = "../../tristan_KSI/offset_test/int10_out/flds.tot." + tstr
    #mybase = "/rigel/home/ls3326/tigress/frb_reconnection/shear/kr/test_grav/output/flds.tot."+tstr
    print(t)
    mybase = "../../lorenzo_setup/reg_grav/output/flds.tot."+tstr
    myf = h5py.File(mybase,'r')
    
    
    dens = myf['dens'][0,:,:]
    dens_pert = return_dens_pert(dens)

    fourier_amp_m1 = fourier_analysis(dens_pert,1)
    fourier_amp_m2 = fourier_analysis(dens_pert,2)
    fourier_amp_m3 = fourier_analysis(dens_pert,3)


    fig, (ax0,ax1) = plt.subplots(2,1)
    ax1.plot(fourier_amp_m1,label='$m=1$')
    ax1.plot(fourier_amp_m2,label='$m=2$')
    ax1.plot(fourier_amp_m3,label='$m=3$')
    plt.legend(frameon=False)
    ax1.set_xlabel('x')
    ax1.set_ylabel('Fourier Amplitude')
    ax1.set_ylim(0,1)
    ax0.imshow(dens,origin='lower',cmap='Greys')
    fig.savefig('lorenzo_setup_reg_grav/fourier'+tstr+'.png',bbox_inches='tight',dpi=300)
    
    '''
    plt.imshow(dens_pert,origin='lower',cmap='Greys')
    plt.colorbar()
    plt.savefig('dens_pert/out_test'+tstr+'.png',bbox_inches='tight',dpi=300)
    plt.close()
    '''


    '''
    dens_pert,dens_cut = dens_analysis(dens)
    print(dens_pert)
    fig, (ax0,ax1) = plt.subplots(2,1)
    ax1.plot(dens_pert)
    ax1.set_xlabel('x')
    ax1.set_ylabel('Average overdensity')#$\rho/\bar{\rho}$')
    ax1.set_ylim(0,.7)
    ax0.imshow(dens_cut,vmin=0,vmax=25,origin='lower',cmap='Greys')
    fig.savefig('dens_pert/'+tstr+'.png',bbox_inches='tight',dpi=300)
    plt.close()
    '''
