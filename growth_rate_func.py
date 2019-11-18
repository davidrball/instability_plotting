import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import h5py 
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 9})

 



def dens_analysis(dens):
    ylen,xlen = np.shape(dens) 
    #we really only want to look at a region within a width around the transition zone
    x0 = int(xlen/4.)
    xlow = x0 - xscan
    xup = x0 + xscan
    cut_dens = dens[:,xlow:xup]
    #ok now that we've focused in on the region we want, let's start to calculate the relevant values
    ylen, xlen = np.shape(cut_dens)
    ratio_list = [] #the quantity we're extracting that's a function of 
    for i in range(xlen):
        #loop in the direction of gravity
        dens_strip = cut_dens[:,i]
        avg_dens = np.mean(dens_strip)
        tmp_ratio_list = []
        for j in range(ylen):#loop perp to grav to measuredeviations from avg density
            mydens = dens_strip[j]
            ratio=mydens/avg_dens
            ratio_minus_one = np.abs(ratio-1)
            tmp_ratio_list.append(ratio_minus_one)
        mymean = np.mean(tmp_ratio_list) 
        
        ratio_list.append(mymean)

    return ratio_list,cut_dens

        #now we can average together the deviations from avg density at a given x, along y

    
def return_dens_pert(dens):
    ylen,xlen = np.shape(dens)
    dens_pert=np.zeros(np.shape(dens))

    for i in range(xlen):
        dens_strip = dens[:,i]
        avg_dens = np.mean(dens_strip)
        tmp_ratio = []
        for j in range(ylen):
            mydens = dens_strip[j]
            ratio = mydens/avg_dens - 1
            tmp_ratio.append(ratio)
        dens_pert[:,i]=tmp_ratio

    return dens_pert


def fourier_analysis(dens_pert,m):
    #give it the perturbed density
    ylen, xlen = np.shape(dens_pert)
    myk = 2.*np.pi*m/ylen

    yarr = np.arange(ylen)
    #calculate integral along y at each x
    
    fourier_amplitude = []

    for i in range(xlen):
        dens_strip = dens_pert[:,i]
        myexp = np.exp(-1j*myk*yarr)
        mysum = np.sum(dens_strip*myexp)
        myabs = np.absolute(mysum)/ylen
        fourier_amplitude.append(myabs)
    return fourier_amplitude #returns the fourier amplitude of mode m as a function of x (the direction of gravity)
