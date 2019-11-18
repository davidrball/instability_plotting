import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import h5py 
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.family'] = 'STIXGeneral'
t0=1
tf=130

b_energy_list=[]
therm_energy_list=[]
tlist=[]
for t in range(t0,tf):
    tstr = "%03d" % t
    mybase = "/home/u21/davidrball/david/tristan_KSI/sig2_bangle50/output/flds.tot." + tstr
    myf = h5py.File(mybase,'r')
    print(tstr)
    dens = myf['dens'][0,:,:]
    bz = myf['bz'][0,:,:]
    bx = myf['bx'][0,:,:]
    by = myf['by'][0,:,:]
    b2=bz**2 + by**2 + bx**2
    b2/=(4*np.pi)
    
    btot = np.sum(b2)
    ken = myf['ken'][0,:,:]
    tlist.append(t)
    thermtot = np.sum(ken)
    b_energy_list.append(btot)
    therm_energy_list.append(thermtot)
#plt.plot(tlist,b_energy_list,color="Blue",label='Magnetic Energy')
plt.plot(tlist,therm_energy_list,color="Red",label="Particle Energy (bulk+thermal)")

plt.savefig('therm_energy.png',bbox_inches='tight')

    
