import scipy.io as sio
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
plt.rcParams.update({'font.size': 22})

plt.interactive(True)


dx=0.05
dz=0.025

ni=34
nj=49
nk=34


viscos=1./5200.

# loop over nfiles 
nfiles=2
#initialize fields
u3d_nfiles=np.zeros((ni,nj,nk,nfiles+1))
v3d_nfiles=np.zeros((ni,nj,nk,nfiles+1))
w3d_nfiles=np.zeros((ni,nj,nk,nfiles+1))
w3d_nfiles=np.zeros((ni,nj,nk,nfiles+1))
te3d_nfiles=np.zeros((ni,nj,nk,nfiles+1))
eps3d_nfiles=np.zeros((ni,nj,nk,nfiles+1))

for n in range(1,nfiles+1):
   print('file no: n=',n)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  read v_1 & transform v_1 to a 3D array (file 1)
   uvw = sio.loadmat('u'+str(n)+'_IDD_PANS.mat')
   ttu=uvw['u'+str(n)+'_IDD_PANS']
   u3d_nfiles[:,:,:,n]= np.reshape(ttu,(nk,nj,ni))  #v_1 velocity
# N.B.- We don't have to swich axex since python and fortran stores an array in the same way

   uvw = sio.loadmat('v'+str(n)+'_IDD_PANS.mat')
   tt=uvw['v'+str(n)+'_IDD_PANS']
   v3d_nfiles[:,:,:,n]= np.reshape(tt,(nk,nj,ni))  #v_2 velocity

   uvw = sio.loadmat('w'+str(n)+'_IDD_PANS.mat')
   tt=uvw['w'+str(n)+'_IDD_PANS']
   w3d_nfiles[:,:,:,n]= np.reshape(tt,(nk,nj,ni))  #v_3 velocity

   uvw = sio.loadmat('te'+str(n)+'_IDD_PANS.mat')
   tt=uvw['te'+str(n)+'_IDD_PANS']
   te3d_nfiles[:,:,:,n]= np.reshape(tt,(nk,nj,ni))  #modeled turbulent kinetic energy

   uvw = sio.loadmat('eps'+str(n)+'_IDD_PANS.mat')
   tt=uvw['eps'+str(n)+'_IDD_PANS']
   eps3d_nfiles[:,:,:,n]= np.reshape(tt,(nk,nj,ni))  #modeled turbulent dissipation


# merge nfiles. This means that new ni = nfiles*ni
u3d=u3d_nfiles[:,:,:,1]
v3d=v3d_nfiles[:,:,:,1]
w3d=w3d_nfiles[:,:,:,1]
te3d=te3d_nfiles[:,:,:,1]
eps3d=eps3d_nfiles[:,:,:,1]
for n in range(1,nfiles+1):
   u3d=np.concatenate((u3d, u3d_nfiles[:,:,:,n]), axis=0)
   v3d=np.concatenate((v3d, v3d_nfiles[:,:,:,n]), axis=0)
   w3d=np.concatenate((w3d, w3d_nfiles[:,:,:,n]), axis=0)
   te3d=np.concatenate((te3d, te3d_nfiles[:,:,:,n]), axis=0)
   eps3d=np.concatenate((eps3d, eps3d_nfiles[:,:,:,n]), axis=0)



# x coordinate direction = index 0, first index
# y coordinate direction = index 1, second index
# z coordinate direction = index 2, third index



ni=len(u3d)

x=dx*ni
z=dz*nk


umean=np.mean(u3d, axis=(0,2))
vmean=np.mean(v3d, axis=(0,2))
wmean=np.mean(w3d, axis=(0,2))
temean=np.mean(te3d, axis=(0,2))
epsmean=np.mean(eps3d, axis=(0,2))

# face coordinates
yc = np.loadtxt("yc.dat")

# cell cener coordinates
y= np.zeros(nj)
dy=np.zeros(nj)
for j in range (1,nj-1):
# dy = cell width
   dy[j]=yc[j]-yc[j-1]
   y[j]=0.5*(yc[j]+yc[j-1])

y[nj-1]=yc[nj-1]
tauw=viscos*umean[1]/y[1]
ustar=tauw**0.5
yplus=y*ustar/viscos

DNS=np.genfromtxt("LM_Channel_5200_mean_prof.dat", dtype=None,comments="%")
y_DNS=DNS[:,0]
yplus_DNS=DNS[:,1]
u_DNS=DNS[:,2]

DNS=np.genfromtxt("LM_Channel_5200_vel_fluc_prof.dat", dtype=None,comments="%")

u2_DNS=DNS[:,2]
v2_DNS=DNS[:,3]
w2_DNS=DNS[:,4]
uv_DNS=DNS[:,5]

k_DNS=0.5*(u2_DNS+v2_DNS+w2_DNS)

# find equi.distant DNS cells in log-scale
xx=0.
jDNS=[1]*40
for i in range (0,40):
   i1 = (np.abs(10.**xx-yplus_DNS)).argmin()
   jDNS[i]=int(i1)
   xx=xx+0.2

############################### U log
fig1,ax1 = plt.subplots()
plt.subplots_adjust(left=0.20,bottom=0.20)

plt.semilogx(yplus,umean/ustar,'b-')
plt.semilogx(yplus_DNS[jDNS],u_DNS[jDNS],'bo')
plt.axis([1, 8000, 0, 31])
plt.ylabel("$U^+$")
plt.xlabel("$y^+$")

############################### U linear plus zoom
fig1,ax1 = plt.subplots()
plt.subplots_adjust(left=0.20,bottom=0.20)

plt.plot(yplus,umean/ustar,'b-')
plt.plot(yplus_DNS[jDNS],u_DNS[jDNS],'bo')
plt.axis([1, 5200, 0, 31])
plt.ylabel("$U^+$")
plt.xlabel("$y^+$")


# make a zoom
axins1 = inset_axes(ax1, width="50%", height="50%", loc='lower right', borderpad=1) # borderpad = space to the main figure
# reduce fotnsize 
axins1.tick_params(axis = 'both', which = 'major', labelsize = 10)
#axins1.yaxis.set_label_position("left")
axins1.xaxis.set_label_position("top")  # put x labels at the top
axins1.xaxis.tick_top()                 # put x ticks at the top
#axins1.yaxis.tick_left()
#axins1.xaxis.tick_bottom()
plt.plot(yplus,umean/ustar,'b-')
plt.plot(yplus_DNS[jDNS],u_DNS[jDNS],'bo')
plt.axis([0, 30, 0, 14])

# Turn ticklabels of insets off
#axins1.tick_params(labelleft=False, labelbottom=False)
#axins1.tick_params(labelleft=False)


plt.savefig('u_linear-zoom_python.eps')

