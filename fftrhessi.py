from __future__ import division
import numpy as np
import healpy as hp
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from numpy import c_
from scipy import stats
import kuiper as kp
from scipy.optimize import curve_fit

#FT the data by hand and creating plots out of it

dat1=np.loadtxt('RHESSI_dat')

#dat[:,0] is arrival time
#dat[:,1] is energy in keV
#dat[:,2] is detector number

psat=4.07
pstar=7.54773
fsat=1/psat
fstar=1/pstar
fsat2=2*fsat
fstar2=2*fstar

#Look only at 100 s and 9 detectors
ind=np.where((dat1[:,0] > 32) & (dat1[:,0] < 132))
dat1=dat1[ind[0],:]
dat1=dat1[np.where((dat1[:,2] < 17.5) & (dat1[:,2] > 8.5)),:][0]

#The frequencies we're interested in, plus a background to compare to
freq=np.array([fsat,fsat2,fstar,fstar2,fsat+fstar,fsat2+fstar,fsat+fstar2,fsat2+fstar2,fsat-fstar,fsat2-fstar,-fsat+fstar2,fsat2-fstar2])
fbg=np.linspace(0,0.8,801)
fall=np.append(freq,fbg)
f=np.transpose(np.sort(fall))
w=2*np.pi*f

t=dat1[:,0]
detind=dat1[:,2]
alldets=np.linspace(9,17,9) #The detectors we're interested in

ft=[] #The Fourier transformed values

#Taking the fourier transform
for i in range(0,len(alldets)):
	tdet=t[np.where(detind == alldets[i])]
	fcs=np.outer(tdet,w)
	fcs2=np.exp(-1j*fcs)
	su=np.sum(fcs2,0)
	ft=np.append(ft,su)
	print('Detector ' + str(i+9) + ' done')

ftf=np.reshape(ft,(len(alldets),len(w)))
f2=abs(ftf)**2 #Getting the power spectrum

#Notes for all plots
# The lines appear on the frequencies we're interested in (only go to second harmonics)
# Lines are colour coded - green is a satellite-only harmonic, yellow is star-only harmonic, magenta is a combination
# None of these plots save automaticcally, they have to be done manually or have a savefig line added

#------------------------------------------------------Plotting the power spectrum for each detector

for i in range(0,len(alldets)):
	plt.plot(f,np.log10(f2[i,:]),label=str(i+9))

for i in range(0,2):
	plt.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(2,4):
	plt.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	plt.plot([freq[i],freq[i]],[0,10],'m--')

plt.legend(loc=3)
plt.xlim(0,0.8)
plt.ylim(3,10)
plt.show()

#------------------------------------------------------Plotting the power spectrum for the average and 8,9

avf2=f2[0:6,:]
avgf2=np.mean(avf2,0)

plt.plot(f,np.log10(avgf2),'k',label='Average')
plt.plot(f,np.log10(f2[7,:]),'b',label='Detector 16')
plt.plot(f,np.log10(f2[8,:]),'r',label='Detector 17')

for i in range(0,2):
	plt.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(2,4):
	plt.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	plt.plot([freq[i],freq[i]],[0,10],'m--')

plt.legend(loc=3)
plt.xlim(0,0.8)
plt.ylim(3,10)
plt.show()

#------------------------------------------------------Plotting the power spectrum for the average and 8,9 focussed on four sets of large peaks

asdf, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
#
ax1.plot(f,np.log10(avgf2),'k')
ax1.plot(f,np.log10(f2[7,:]),'b')
ax1.plot(f,np.log10(f2[8,:]),'r')

for i in range(0,2):
	ax1.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(2,4):
	ax1.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	ax1.plot([freq[i],freq[i]],[0,10],'m--')

ax1.set_xlim(0.1,0.14)
ax1.set_ylim(5,9)
#
ax2.plot(f,np.log10(avgf2),'k')
ax2.plot(f,np.log10(f2[7,:]),'b')
ax2.plot(f,np.log10(f2[8,:]),'r')

for i in range(0,2):
	ax2.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(2,4):
	ax2.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	ax2.plot([freq[i],freq[i]],[0,10],'m--')

ax2.set_xlim(0.22,0.28)
ax2.set_ylim(4,9)
#
ax3.plot(f,np.log10(avgf2),'k')
ax3.plot(f,np.log10(f2[7,:]),'b')
ax3.plot(f,np.log10(f2[8,:]),'r')

for i in range(0,2):
	ax3.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(2,4):
	ax3.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	ax3.plot([freq[i],freq[i]],[0,10],'m--')

ax3.set_xlim(0.37,0.41)
ax3.set_ylim(4,9)
#
ax4.plot(f,np.log10(avgf2),'k')
ax4.plot(f,np.log10(f2[7,:]),'b')
ax4.plot(f,np.log10(f2[8,:]),'r')

for i in range(0,2):
	ax4.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(2,4):
	ax4.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	ax4.plot([freq[i],freq[i]],[0,10],'m--')

ax4.set_xlim(0.48,0.54)
ax4.set_ylim(4,8)
#
plt.show()

#------------------------------------------------------ Plotting the power spectrum for the average and 8,9 focussed on each peak containing 2*fsat
asdf, ((ax1, ax2), (ax3, ax4),(ax5, ax6)) = plt.subplots(3, 2)
#
ax1.plot(f,np.log10(avgf2),'k')
ax1.plot(f,np.log10(f2[7,:]),'b')
ax1.plot(f,np.log10(f2[8,:]),'r')

for i in range(0,2):
	ax1.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(2,4):
	ax1.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	ax1.plot([freq[i],freq[i]],[0,10],'m--')

ax1.set_xlim(freq[11]-0.05,freq[11]+0.05)
ax1.set_ylim(4,9)
#
ax2.plot(f,np.log10(avgf2),'k')
ax2.plot(f,np.log10(f2[7,:]),'b')
ax2.plot(f,np.log10(f2[8,:]),'r')

for i in range(0,2):
	ax2.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(2,4):
	ax2.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	ax2.plot([freq[i],freq[i]],[0,10],'m--')

ax2.set_xlim(freq[9]-0.05,freq[9]+0.05)
ax2.set_ylim(3,9)
#
ax3.plot(f,np.log10(avgf2),'k')
ax3.plot(f,np.log10(f2[7,:]),'b')
ax3.plot(f,np.log10(f2[8,:]),'r')

for i in range(0,2):
	ax3.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(2,4):
	ax3.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	ax3.plot([freq[i],freq[i]],[0,10],'m--')

ax3.set_xlim(freq[1]-0.05,freq[1]+0.05)
ax3.set_ylim(4,8)
#
ax4.plot(f,np.log10(avgf2),'k')
ax4.plot(f,np.log10(f2[7,:]),'b')
ax4.plot(f,np.log10(f2[8,:]),'r')

for i in range(0,2):
	ax4.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(2,4):
	ax4.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	ax4.plot([freq[i],freq[i]],[0,10],'m--')

ax4.set_xlim(freq[5]-0.05,freq[5]+0.05)
ax4.set_ylim(4,8)
#
ax5.plot(f,np.log10(avgf2),'k')
ax5.plot(f,np.log10(f2[7,:]),'b')
ax5.plot(f,np.log10(f2[8,:]),'r')

for i in range(0,2):
	ax5.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(2,4):
	ax5.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	ax5.plot([freq[i],freq[i]],[0,10],'m--')

ax5.set_xlim(freq[7]-0.05,freq[7]+0.05)
ax5.set_ylim(4,8)
#
ax6.plot(f,np.log10(avgf2),'k')
ax6.plot(f,np.log10(f2[7,:]),'b')
ax6.plot(f,np.log10(f2[8,:]),'r')

for i in range(0,2):
	ax6.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(2,4):
	ax6.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	ax6.plot([freq[i],freq[i]],[0,10],'m--')

ax6.set_xlim(0,0.8)
ax6.set_ylim(3,9)
#
plt.show()
#------------------------------------------------------Plotting the real and imaginary parts of the FT for the average and 8,9

avftf=ftf[0:6,:]
avgftf=np.mean(avftf,0)

plt.plot(f,avgftf.real,'k-',label='Average')
plt.plot(f,ftf[7,:].real,'b-',label='Detector 16')
plt.plot(f,ftf[8,:].real,'r-',label='Detector 17')
plt.plot(f,avgftf.imag,'k--')
plt.plot(f,ftf[7,:].imag,'b--')
plt.plot(f,ftf[8,:].imag,'r--')

for i in range(0,2):
	plt.plot([freq[i],freq[i]],[-100000,100000],'g:')

for i in range(2,4):
	plt.plot([freq[i],freq[i]],[-100000,100000],'y:')

for i in range(4,len(freq)):
	plt.plot([freq[i],freq[i]],[-100000,100000],'m:')

plt.legend(loc=1)
plt.xlim(0,0.8)
plt.ylim(-20000,20000)
plt.show()

#------------------------------------------------------Plotting the real and imaginary parts of the FT for the average and 8,9 focussed on four sets of large peaks
 
asdf, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
#
ax1.plot(f,avgftf.real,'k-')
ax1.plot(f,ftf[7,:].real,'b-')
ax1.plot(f,ftf[8,:].real,'r-')
ax1.plot(f,avgftf.imag,'k--')
ax1.plot(f,ftf[7,:].imag,'b--')
ax1.plot(f,ftf[8,:].imag,'r--')

for i in range(2,4):
	ax1.plot([freq[i],freq[i]],[-1000000,1000000],'y--')

for i in range(4,len(freq)):
	ax1.plot([freq[i],freq[i]],[-1000000,1000000],'m--')

ax1.set_xlim(0.1,0.14)
ax1.set_ylim(-25000,25000)
#ax1.yscale("symlog")
#
ax2.plot(f,avgftf.real,'k-')
ax2.plot(f,ftf[7,:].real,'b-')
ax2.plot(f,ftf[8,:].real,'r-')
ax2.plot(f,avgftf.imag,'k--')
ax2.plot(f,ftf[7,:].imag,'b--')
ax2.plot(f,ftf[8,:].imag,'r--')

for i in range(0,2):
	ax2.plot([freq[i],freq[i]],[-1000000,1000000],'g--')

for i in range(2,4):
	ax2.plot([freq[i],freq[i]],[-1000000,1000000],'y--')

for i in range(4,len(freq)):
	ax2.plot([freq[i],freq[i]],[-1000000,1000000],'m--')

ax2.set_xlim(0.22,0.28)
ax2.set_ylim(-20000,20000)
#ax2.yscale("symlog")
#
ax3.plot(f,avgftf.real,'k-')
ax3.plot(f,ftf[7,:].real,'b-')
ax3.plot(f,ftf[8,:].real,'r-')
ax3.plot(f,avgftf.imag,'k--')
ax3.plot(f,ftf[7,:].imag,'b--')
ax3.plot(f,ftf[8,:].imag,'r--')

for i in range(4,len(freq)):
	ax3.plot([freq[i],freq[i]],[-1000000,1000000],'m--')

ax3.set_xlim(0.37,0.41)
ax3.set_ylim(4,9)
ax3.set_ylim(-15000,15000)
#ax3.yscale("symlog")
#
ax4.plot(f,avgftf.real,'k-')
ax4.plot(f,ftf[7,:].real,'b-')
ax4.plot(f,ftf[8,:].real,'r-')
ax4.plot(f,avgftf.imag,'k--')
ax4.plot(f,ftf[7,:].imag,'b--')
ax4.plot(f,ftf[8,:].imag,'r--')

for i in range(0,2):
	ax4.plot([freq[i],freq[i]],[-1000000,1000000],'g--')

for i in range(4,len(freq)):
	ax4.plot([freq[i],freq[i]],[-1000000,1000000],'m--')

ax4.set_xlim(0.48,0.54)
ax4.set_ylim(4,8)
ax4.set_ylim(-7500,7500)
#ax4.yscale("symlog")
#
plt.show()

#------------------------------------------------------Plotting the power spectrum for the average and 3,5

asdf, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
#
ax1.plot(f,np.log10(avgf2),'k')
ax1.plot(f,np.log10(f2[2,:]),'b')
ax1.plot(f,np.log10(f2[4,:]),'r')

for i in range(2,4):
	ax1.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	ax1.plot([freq[i],freq[i]],[0,10],'m--')

ax1.set_xlim(0.1,0.14)
ax1.set_ylim(5,9)
#
ax2.plot(f,np.log10(avgf2),'k')
ax2.plot(f,np.log10(f2[2,:]),'b')
ax2.plot(f,np.log10(f2[4,:]),'r')

for i in range(0,2):
	ax2.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(2,4):
	ax2.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	ax2.plot([freq[i],freq[i]],[0,10],'m--')

ax2.set_xlim(0.22,0.28)
ax2.set_ylim(4,9)
#
ax3.plot(f,np.log10(avgf2),'k')
ax3.plot(f,np.log10(f2[2,:]),'b')
ax3.plot(f,np.log10(f2[4,:]),'r')

for i in range(4,len(freq)):
	ax3.plot([freq[i],freq[i]],[0,10],'m--')

ax3.set_xlim(0.37,0.41)
ax3.set_ylim(4,9)
#
ax4.plot(f,np.log10(avgf2),'k')
ax4.plot(f,np.log10(f2[2,:]),'b')
ax4.plot(f,np.log10(f2[4,:]),'r')

for i in range(0,2):
	ax4.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(4,len(freq)):
	ax4.plot([freq[i],freq[i]],[0,10],'m--')

ax4.set_xlim(0.48,0.54)
ax4.set_ylim(4,8)
#
plt.show()

#------------------------------------------------------ Plotting the power spectrum for the average and 4,6

asdf, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
#
ax1.plot(f,np.log10(avgf2),'k')
ax1.plot(f,np.log10(f2[3,:]),'b')
ax1.plot(f,np.log10(f2[5,:]),'r')

for i in range(2,4):
	ax1.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	ax1.plot([freq[i],freq[i]],[0,10],'m--')

ax1.set_xlim(0.1,0.14)
ax1.set_ylim(5,9)
#
ax2.plot(f,np.log10(avgf2),'k')
ax2.plot(f,np.log10(f2[3,:]),'b')
ax2.plot(f,np.log10(f2[5,:]),'r')

for i in range(0,2):
	ax2.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(2,4):
	ax2.plot([freq[i],freq[i]],[0,10],'y--')

for i in range(4,len(freq)):
	ax2.plot([freq[i],freq[i]],[0,10],'m--')

ax2.set_xlim(0.22,0.28)
ax2.set_ylim(4,9)
#
ax3.plot(f,np.log10(avgf2),'k')
ax3.plot(f,np.log10(f2[3,:]),'b')
ax3.plot(f,np.log10(f2[5,:]),'r')

for i in range(4,len(freq)):
	ax3.plot([freq[i],freq[i]],[0,10],'m--')

ax3.set_xlim(0.37,0.41)
ax3.set_ylim(4,9)
#
ax4.plot(f,np.log10(avgf2),'k')
ax4.plot(f,np.log10(f2[3,:]),'b')
ax4.plot(f,np.log10(f2[5,:]),'r')

for i in range(0,2):
	ax4.plot([freq[i],freq[i]],[0,10],'g--')

for i in range(4,len(freq)):
	ax4.plot([freq[i],freq[i]],[0,10],'m--')

ax4.set_xlim(0.48,0.54)
ax4.set_ylim(4,8)
#
plt.show()

#------------------------------------------------------
