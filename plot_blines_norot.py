import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 8})
plt.set_cmap("Blues")




#define the simulation matrix, right now it's 3x3
#first index picks out guide field strength, 2nd is temp




#guide field 0, delgam=.00005, .0005, .005
#myname = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_run2/output/flds.tot."

myname = "/home/u21/davidrball/david/tristan_KSI/sig2_bangle50/output/flds.tot."

fig, ax0  = plt.subplots(1,1)
t = 120
tstr = "%03d" % t

myname += tstr

print(myname)

#xscan = 50
istep = 6
c_omp = 5



           
mydens = h5py.File(myname,'r')['dens'][0,:,:]

bx = h5py.File(myname,'r')['bx'][0,:,:]
by = h5py.File(myname,'r')['by'][0,:,:]
bz = h5py.File(myname,'r')['bz'][0,:,:]
myshape = np.shape(mydens)

xlen = myshape[1]
ylen = myshape[0]

xscan=ylen/2

x0 = 3*int(xlen/4)
xlow=int(x0-xscan)
xup=int(x0+xscan)
mydens = mydens[:,xlow:xup]


bx = bx[:,xlow:xup]
by = by[:,xlow:xup]
bz = bz[:,xlow:xup]

b2 = bx**2 + by**2 + bz**2


xext = xscan*istep/c_omp
yext = ylen*istep/c_omp

print(xext,yext)

bmax = np.max(b2)

myim = ax0.imshow(b2,extent = [-yext/2, yext/2, -xext, xext],origin='lower',vmax=bmax/4.)
myshape = np.shape(bx)
ylen = myshape[0]
xlen = myshape[1]

#working version
yarr = np.linspace(-xext,xext,ylen)
xarr = np.linspace(-yext/2,yext/2,xlen)

xarr = np.linspace(-xext,xext,ylen)
yarr = np.linspace(-yext/2,yext/2,xlen)


xhlf=xlen/2
#seed_points = np.array([[0,0,0],[-50,0,50]])
#seed_points=np.array([[xarr[xhlf/4],xarr[xhlf],xarr[3*xhlf/4]],[yarr[0],yarr[2],yarr[4]]])

seed_points = np.array([[-100,-10,0,10,100],[0, 0, 0,0,0]])

ax0.streamplot(xarr, yarr,bx,by,color="Red",linewidth=0.5,arrowsize=.5,density=[.5,2.5])#,start_points=seed_points.T)#, density=[2,2])

#ax0.set_xlim(-xscan,xscan)
#ax0.set_ylim(-yext/2,yext/2)


#ax0.quiver(xarr, yarr, by, bx)

#ax0.set_adjustable('box-forced')

#cbar_ax = fig.add_axes([.905,.145,.015,.7])
#cb= fig.colorbar(myim,cax=cbar_ax,orientation="vertical")
#cbar_ax.set_ylabel("Density ($N_{ppc}$ / $N_{ppc0}$)",fontsize=10)#, labelpad=-10)

#cbytick_obj = plt.getp(cb.ax.axes,'xticklabels')
#plt.setp(cbytick_obj,fontsize='x-small')

#cbytick_obj = plt.getp(cb.ax.axes,'xticklabels')
#plt.setp(cbytick_obj,fontsize='x-small')


xstr = "$x \; (c/\omega_{p})$"
ystr = "$y \; (c/\omega_{p})$"
ax0.set_ylabel(ystr)                                        
ax0.set_xlabel(xstr)
plt.savefig('testing_blines.png',dpi=300,bbox_inches='tight')
