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
tf=50
c=.45
interval=500
my0=1200.
dstripe = 200.
istep = 6.
dstripe_down = dstripe/istep
my0_down = my0/istep
dstripe_m =  int(my0_down/dstripe_down)

print('m associated with dstripe : ', dstripe_m)
t_list = []
#m_list = [1,2,dstripe_m]#,6,10]
m_list = [2,6,12]
col_list = ["Blue","Red","Green","Orange","Magenta","Black","Cyan"]
mnum = len(m_list)
fm_list = []
#mybase = "../../relativistic_KSI/sig1/output/flds.tot.001"
mybase = "/home/u21/davidrball/david/tristan_KSI/"
base1 = mybase + "bangle8_flipcur/output/flds.tot."
base2 = mybase + "sig2_bangle200/output/flds.tot."
base3 = mybase + "sig2_bangle500/output/flds.tot."


growth_m = 2.
growth_k = 2*np.pi*growth_m/(ylen*istep)
growth_m2 = 12
growth_k2 = 2*np.pi*growth_m2/(ylen*istep)
for m in m_list:
    fm_list.append([])
for t in range(t0,tf):
    t_list.append(t)    
    tstr = "%03d" % t
    print(t)
    newbase1 = mybase1 + tstr
    newbase2 = mybase2 + tstr
    newbase3 = mybase3 + tstr
    myf1 = h5py.File(newbase1,'r')
    myf2 = h5py.File(newbase2,'r')
    myf3 = h5py.File(newbase3,'r')
    
    
    dens = myf['dens'][0,:,:] #can cut dens here to isolate one region and speed up calculation
    dens_pert = return_dens_pert(dens)
    ylen,xlen = np.shape(dens)

    x0=int(xlen/4.)
    x1=int(3*xlen/4.)
    xscan=int(dstripe/istep/2)
    xscan/=2
    print('my xscan is : ',xscan)
    xlow0 = int(x0-xscan+2)
    xup0 = int(x0+xscan+2)
    xlow1 = int(x1-xscan-2)
    xup1 = int(x1+xscan-2)
    print('xlow0, xup0',xlow0,xup0)

    for i in range(mnum):
        m = m_list[i]
        fourier_amp = fourier_analysis(dens_pert,m)
        fm_0 = np.mean(fourier_amp[xlow0:xup0])
        fm_1 = np.mean(fourier_amp[xlow1:xup1])

        fm_avg = fm_0
        #fm_avg = .5*(fm_0 + fm_1)
        fm_list[i].append(fm_avg)
        


    fig, (ax0,ax1) = plt.subplots(2,1)


    for i in range(mnum):
        #print(i)
        #print(fm_list)
        #print(col_list)
        m = m_list[i]
        ax1.plot(t_list,np.log10(fm_list[i]),label='$m='+str(m)+'$',color=col_list[i])

    ax1.plot(t_list,np.log10(growth_list),color="Blue",label="Linear Theory (m=2)",linestyle='dashed')
    ax1.plot(t_list, np.log10(growth_list2),color="Green",label="Linear Theory (m=12)",linestyle='dashed')

    ax1.set_xlim(t0,tf)
    ax1.set_ylim(-2.5,0)
    plt.legend(frameon=False,loc='lower right',prop={'size':12})
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Log($\\bar{f}_m$)')

    ax0.imshow(dens,origin='lower',cmap='Greys')
    #fig.savefig('dens_pert/fourier_time'+tstr+'.png',bbox_inches='tight',dpi=300)
    

    ax0.scatter(xlow0,100,color="Red",marker='x')
    ax0.scatter(xup0,100,color="Red",marker='x')

    ax0.scatter(xlow1,100,color="Blue",marker='x')
    ax0.scatter(xup1,100,color="Blue",marker='x')
    fig.savefig('compare/fourier_time'+tstr+'.png',bbox_inches='tight',dpi=300)
    plt.close()
    
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
