import numpy as np
import matplotlib.pyplot as plt
import sys
#from scipy.signal import welch, hanning
from scipy.signal import welch, hann
plt.rcParams.update({'font.size': 22})

# ***** read u
# /chalmers/users/lada/python-DES-code/channel-5200-IDD-PANS-inlet-synt-MTF271/u-time-history-i70.dat

plt.interactive(True)

n1=7500
n2=15000

data = np.genfromtxt("u-time-history-i70.dat", dtype=None)
t=data[n1:n2,0] #time
u1=data[n1:n2,1]   #v_1 at point 1
u2=data[n1:n2,2]   #v_1 at point 2
u3=data[n1:n2,3]   #v_1 at point 3
u4=data[n1:n2,4]   #v_1 at point 4
u5=data[n1:n2,5]   #v_1 at point 5

data = np.genfromtxt("w-time-history-i70.dat", dtype=None)
w1=data[n1:n2,1]   #v_3 at point 1
w2=data[n1:n2,2]   #v_3 at point 2
w3=data[n1:n2,3]   #v_3 at point 3
w4=data[n1:n2,4]   #v_3 at point 4
w5=data[n1:n2,5]   #v_3 at point 5

dx=0.1
dz=0.05
dt=t[1]-t[0]

# plot time history 
fig1,ax1 = plt.subplots()
plt.subplots_adjust(left=0.20,bottom=0.20)

plt.plot(t,u2,'b-')
# plot every 30th time step
plt.plot(t[::30],u2[::30],'bo')
plt.plot(t,u4,'r-')
plt.xlabel('$t$')
plt.axis([10, 11, 9,25])
plt.ylabel('$u$')
plt.savefig('u-time.eps')

#autocoreelation
u1_fluct=u1-np.mean(u1)
two=np.correlate(u1_fluct,u1_fluct,'full')
# find index of max
nmax=np.argmax(two)
# and its value
two_max=np.max(two)
# Pick the right half and normalize
two_sym_norm=two[nmax:]/two_max

u2_fluct=u2-np.mean(u2)
two2=np.correlate(u2_fluct,u2_fluct,'full')
# find max
nmax2=np.argmax(two2)
# and its value
two_max2=np.max(two2)
# Pick the right half and normalize
two_sym_norm2=two2[nmax2:]/two_max2

u3_fluct=u1-np.mean(u3)
two3=np.correlate(u3_fluct,u3_fluct,'full')
# find max
nmax3=np.argmax(two3)
# and its value
two_max3=np.max(two3)
# Pick the right half and normalize
two_sym_norm3=two3[nmax3:]/two_max3

u4_fluct=u4-np.mean(u4)
two4=np.correlate(u4_fluct,u4_fluct,'full')
# find max
nmax4=np.argmax(two4)
# and its value
two_max4=np.max(two4)
# Pick the right half and normalize
two_sym_norm4=two4[nmax4:]/two_max4


u5_fluct=u5-np.mean(u5)
two5=np.correlate(u5_fluct,u5_fluct,'full')
# find max
nmax5=np.argmax(two5)
# and its value
two_max5=np.max(two5)
# Pick the right half and normalize
two_sym_norm5=two5[nmax5:]/two_max5




fig1 = plt.figure("Figure 1")
imax=500;
plt.plot(t[0:imax],two_sym_norm[0:imax],'b-')
plt.xlabel('t')
plt.ylabel('$B_{uu}$')

fig2 = plt.figure("Figure 2")
imax=500;
plt.plot(t[0:imax],two_sym_norm2[0:imax],'b-')
plt.xlabel('t')
plt.ylabel('$B_{uu}$')

fig3 = plt.figure("Figure 3")
imax=500;
plt.plot(t[0:imax],two_sym_norm3[0:imax],'b-')
plt.xlabel('t')
plt.ylabel('$B_{uu}$')

fig4 = plt.figure("Figure 4")
imax=500;
plt.plot(t[0:imax],two_sym_norm4[0:imax],'b-')
plt.xlabel('t')
plt.ylabel('$B_{uu}$')

fig5 = plt.figure("Figure 5")
imax=500;
plt.plot(t[0:imax],two_sym_norm5[0:imax],'b-')
plt.xlabel('t')
plt.ylabel('$B_{uu}$')

#integral timescale 
dt=t[1]-t[0]
int_T_1=np.trapz(two_sym_norm)*dt

dt2=t[1]-t[0]
int_T_2=np.trapz(two_sym_norm2)*dt

dt3=t[1]-t[0]
int_T_3=np.trapz(two_sym_norm3)*dt

dt4=t[1]-t[0]
int_T_4=np.trapz(two_sym_norm4)*dt

dt5=t[1]-t[0]
int_T_5=np.trapz(two_sym_norm5)*dt