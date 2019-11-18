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


 



my0=1200.
dstripe = 25.
istep = 6.
dstripe_down = dstripe/istep
my0_down = my0/istep
dstripe_m =  my0_down/dstripe_down

print('m associated with dstripe : ', dstripe_m)



#myfstart = h5py.File("~/home_david/tristan_KSI/smooth_offset_current_test/int1_out/flds.tot.001")

t_list = [5,10,15,24]

tlabel_list = []
for t in t_list:
    tlabel_list.append('t='+str(t))


#m_list = [1,2,dstripe_m]#,6,10]
m_list = np.arange(1,dstripe_m)

col_list = ["Blue","Red","Green","Orange","Magenta","Black","Cyan"]
mnum = len(m_list)
fm_list = []


c=.45
dstripe = 25.
istep = 6
interval=500
dstripe_k = 2*np.pi/dstripe

mybase = "../../lorenzo_setup/reg_grav/dstripe25_sigext01/output/flds.tot.001"
myf = h5py.File(mybase,'r')
dens = myf['dens'][0,:,:]
ylen,xlen=np.shape(dens)
myf.close()

count=1
growth_m = 1.
growth_k = 2*np.pi*growth_m/(ylen*istep)
btheta=2. #not sure why this is factored in
sigma_ext = 0.1*btheta
g = sigma_ext*dstripe_k*c**2

#analytic_growth_rate = np.sqrt(g/2.)*growth_k**.75 * dstripe**.25
analytic_growth_rate = np.sqrt(growth_k*g / 3.)#thick sheet limit
analytic_growth_rate /= 5

print('analytic growth rate : ',analytic_growth_rate)
N=.0025
growth_list = []


tcount = 0

for t in t_list:
    

    fm_list = []
    tstr = "%03d"%t
    print(t)
    mybase = "../../lorenzo_setup/reg_grav/dstripe25_sigext01/output/flds.tot."+tstr
    myf = h5py.File(mybase,'r')
    
    
    dens = myf['dens'][0,:,:] #can cut dens here to isolate one region and speed up calculation
    dens_pert = return_dens_pert(dens)
    ylen,xlen = np.shape(dens)

    x0=int(xlen/4.)
    x1=3*x0
    xscan=int(dstripe/istep)/2
    print('my xscan is : ',xscan)
    xlow0 = x0-xscan
    xup0 = x0+xscan
    xlow1 = x1-xscan
    xup1 = x1+xscan
    print('xlow0, xup0',xlow0,xup0)

    for i in range(mnum):
        m = m_list[i]
        print(m)
        fourier_amp = fourier_analysis(dens_pert,m)
        fm_0 = np.mean(fourier_amp[xlow0:xup0])
        fm_1 = np.mean(fourier_amp[xlow1:xup1])

        
        fm_avg = .5*(fm_0 + fm_1)
        
        
        fm_list.append(fm_avg)
        


    #fig, (ax0,ax1) = plt.subplots(2,1)

    plt.plot(m_list, fm_list,label=tlabel_list[tcount],color=col_list[tcount])
    plt.xlabel('m')
    plt.ylabel('Average Fourier Amplitude')

    #plt.savefig('testing_growth_spec.png',dpi=300,bbox_inches='tight')
    

    tcount+=1
#plt.yscale('log')
plt.legend(frameon=False,loc='upper right')
plt.savefig('dstripe25_sigext01/growth_spec.png',dpi=300,bbox_inches='tight')
plt.close()

