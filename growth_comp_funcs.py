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
tf=52
c=.45
interval=500
my0=1200.
dstripe = 200.
istep = 6.


m_list=[1,2,6]
col_list = ["Blue","Red","Green","Orange","Magenta","Black","Cyan"]
mnum = len(m_list)

#mybase = "../../relativistic_KSI/sig1/output/flds.tot.001"
mybase = "/home/u21/davidrball/david/tristan_KSI/"
plebase = "/home/u21/davidrball/david/tristan_KSI/ple_runs/"

mypost = "/output/flds.tot."
plepost = "/flds.tot."

base1 = mybase + "sig2_bangle0"+mypost
base2 = plebase + "sig2_bangle25" + plepost
base3 = plebase + "sig2_bangle10" + plepost

def return_growth_list(mybase, t0, tf,m_list):
    t_list=[]
    fm_list=[]
    for m in m_list:
        fm_list.append([])
    for t in range(t0,tf):
        t_list.append(t)    
        tstr = "%03d" % t
        print(t)
        newbase = mybase + tstr
        myf = h5py.File(newbase,'r')
        dens = myf['dens'][0,:,:] #can cut dens here to isolate one region and speed up calculation
        dens_pert = return_dens_pert(dens)
        ylen,xlen = np.shape(dens)
        x0=int(xlen/4.)
        x1=int(3*xlen/4.)
        xscan=int(dstripe/istep/2)
        xscan/=2
        print('my xscan is : ',xscan)
        xlow0 = int(x0-xscan+4)
        xup0 = int(x0+xscan+4)
        print('xlow0, xup0',xlow0,xup0)
        for i in range(mnum):
            m = m_list[i]
            fourier_amp = fourier_analysis(dens_pert,m)
            fm_0 = np.mean(fourier_amp[xlow0:xup0])
            #just calculate at one boundary
            fm_list[i].append(fm_0)
    return fm_list,t_list
            

growth1,t_list=return_growth_list(base1,t0,tf,m_list)
growth2,t_list = return_growth_list(base2,t0,tf,m_list)
growth3,t_list= return_growth_list(base3,t0,tf,m_list)


for i in range(len(m_list)):
    plt.plot(t_list,growth1[i],color=col_list[i],label='$\phi=0$ m='+str(m_list[i]))
    plt.plot(t_list,growth2[i],color=col_list[i],label='$\phi=\pi/25$ m='+str(m_list[i]),linestyle='dashed')
    plt.plot(t_list,growth3[i],color=col_list[i],label='$\phi=\pi/10$ m='+str(m_list[i]),linestyle='dotted')
    plt.legend(frameon=False,prop={'size':13})

#plt.yscale('log')
plt.xlabel('Time')
plt.ylabel('$\\bar{f}_{m}$')
plt.savefig('compare/growth_ple_comp.png',bbox_inches='tight',dpi=300)
