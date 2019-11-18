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
mybangle=50
mysig=2

bstr=str(mybangle)
sigstr=str(mysig)


t0=1
tf=101


c=.45
interval=500
my0=1200.
dstripe = 200.
istep = 6.
dstripe_down = dstripe/istep
my0_down = my0/istep
dstripe_m =  int(my0_down/dstripe_down)

print('m associated with dstripe : ', dstripe_m)



#myfstart = h5py.File("~/home_david/tristan_KSI/smooth_offset_current_test/int1_out/flds.tot.1")

t_list = []

#m_list = [1,2,dstripe_m]#,6,10]
m_list = [2,6,12]
col_list = ["Blue","Red","Green","Orange","Magenta","Black","Cyan"]
mnum = len(m_list)
fm_list = []



#mybase = "../../relativistic_KSI/sig1/output/flds.tot.001"
#mybase = "/home/u21/davidrball/david/tristan_KSI/sig2_bangle0/output/flds.tot."
mybase = "/home/u21/davidrball/david/tristan_KSI/"
#mybase += "sig"+sigstr +"_"+"bangle"+bstr+"/output/flds.tot."
mybase += "ple_runs/sig2_bangle0_nograv2/flds.tot."


myf = h5py.File(mybase+"001",'r')
dens = myf['dens'][0,:,:]
ylen,xlen=np.shape(dens)
print(ylen, xlen)
myf.close()


count=1

growth_m = 6.
growth_k = 2*np.pi*growth_m/(ylen*istep)


growth_m2 = 6
growth_k2 = 2*np.pi*growth_m2/(ylen*istep)

sigma_ext = 0.1
sigma=2.



dstripe_k = 2*np.pi/dstripe #k associated with sheet width
g = sigma_ext*dstripe_k*c**2

ppc0=4.
nstripe = 4 #actual value of nstripe
overdens = (nstripe+1)*ppc0


tburn = 9
#what k to use for growth rate?  let's just choose k of sheet thickness
atwood = (overdens-ppc0)/(overdens+ppc0)
RT_growth_rate = np.sqrt(atwood*g*dstripe_k)

#accounting for transition region
RT_growth_rate = np.sqrt(atwood*g*dstripe_k / (1+dstripe*dstripe_k))

#short wavelength limit of growth rate
RT_growth_rate = np.sqrt(atwood*g/dstripe)


#accounting for relativistic B fields / thermal pressure in enthalpy
atwood_enthalpy = ((nstripe+1)-1+1*sigma)/((nstripe+1)+1+3*sigma)


#accounting for fact that the peak value of g is only at one location, and the average of g is actually lower.  Lower by factor of pi if you average over delta/2 to delta/2

g/=np.pi

new_growth = np.sqrt(g*(nstripe+1)/dstripe)


#accounting for magnetic tension suppressing growth
phi=np.pi/50
sinphi=np.sin(phi)

tension_term = 2*c**2*sigma*growth_k**2*sinphi**2 /(nstripe+2)
RT_growth_rate = np.sqrt(atwood_enthalpy*g*growth_k/(1+dstripe*growth_k)-tension_term)


RT_growth_rate2 = np.sqrt(atwood_enthalpy*g*growth_k2/(1+dstripe*growth_k2))

print('new growth rate : ', new_growth)
print('RT growth : ', RT_growth_rate2)


N=.0025
growth_list = []
growth_list2=[]
analytic_growth_rate = RT_growth_rate
analytic_growth_rate2 = RT_growth_rate2

for m in m_list:
    fm_list.append([])
    #col_list.append('C'+str(count))
    #count+=1
for t in range(t0,tf):
    t_list.append(t)

    analytic_growth = N*np.exp(analytic_growth_rate*interval*(t-tburn))
    growth_list.append(analytic_growth)

    analytic_growth2 = N*np.exp(analytic_growth_rate2*interval*(t-tburn))
    growth_list2.append(analytic_growth2)
    tstr = "%03d" % t
    #mybase = "../../tristan_KSI/offset_test/int10_out/flds.tot." + tstr
    #mybase = "/rigel/home/ls3326/tigress/frb_reconnection/shear/kr/test_grav/output/flds.tot."+tstr
    print(t)
    #mybase = "../../relativistic_KSI/sig1/output/flds.tot."+tstr
    newbase = mybase + tstr
    myf = h5py.File(newbase,'r')
    
    
    dens = myf['dens'][0,:,:] #can cut dens here to isolate one region and speed up calculation
    dens_pert = return_dens_pert(dens)
    ylen,xlen = np.shape(dens)

    x0=int(xlen/4.)
    x1=int(3*xlen/4.)
    xscan=int(dstripe/istep/2)
    #xscan/=2
    xscan=1.5*dstripe/istep
    print('my xscan is : ',xscan)
    xlow0 = int(x0-xscan+6)
    xup0 = int(x0+xscan+6)
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

    #ax1.plot(t_list,np.log10(growth_list),color="Black",label="Linear Theory (m=2, with tension)",linestyle='dashed')
    #ax1.plot(t_list, np.log10(growth_list2),color="Green",label="Linear Theory (m=2, no tension)",linestyle='dashed')

    ax1.set_xlim(t0,tf)
    ax1.set_ylim(-3,-1)
    plt.legend(frameon=False,loc='upper left',prop={'size':12})
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Log($\\bar{f}_m$)')

    ax0.imshow(dens[:,xlow0-150:xup0+150],origin='lower',cmap='Blues',extent=[xlow0-150,xup0+150,0,200])
    #fig.savefig('dens_pert/fourier_time'+tstr+'.png',bbox_inches='tight',dpi=300)
    
    yarr = np.linspace(0,200,100)
    xlowarr = xlow0*np.ones(100)
    xuparr = xup0*np.ones(100)

    ax0.plot(xlowarr,yarr,color="Red",linestyle='--',linewidth=.5)
    ax0.plot(xuparr,yarr,color="Red",linestyle='--',linewidth=.5)

    ax0.set_xticks([])
    ax0.set_yticks([])
    #ax0.scatter(xlow1,100,color="Blue",marker='x')
    #ax0.scatter(xup1,100,color="Blue",marker='x')
    
    #savestr = "bangle"+bstr+"_sig"+sigstr+"/fourier_time"+tstr+".png"
    savestr="ple_bangle0_sig2_nograv/"+'fourier_time'+tstr+'.png'
    fig.savefig(savestr,bbox_inches='tight',dpi=300)
    

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
