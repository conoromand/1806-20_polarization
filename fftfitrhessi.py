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
from scipy.optimize import fmin
import cmath

#FT the data by hand and solving for positions and differences of maxima in different detectors

datr=np.loadtxt('RHESSI_dat')

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
ind=np.where((datr[:,0] > 32) & (datr[:,0] < 132))
datc=datr[ind[0],:]
dat1=datc[np.where((datc[:,2] < 17.5) & (datc[:,2] > 8.5)),:][0]

#The frequencies we're interested in, plus a background to compare to
freq=np.array([fsat,fsat2,fstar,fstar2,fsat+fstar,fsat2+fstar,fsat+fstar2,fsat2+fstar2,fsat-fstar,fsat2-fstar,-fsat+fstar2,fsat2-fstar2])
fbg=np.linspace(0,0.8,801)
fall=np.append(freq,fbg)
f=np.transpose(np.sort(fall))
w=2*np.pi*f

t=dat1[:,0]
detind=dat1[:,2]
alldets=np.linspace(9,17,9) #The detectors we're interested in

ft=[]#The Fourier transformed values

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

#defining regions of interest
ind1=np.where((f >= 0.125) & (f <= 0.14))[0]
ind2=np.where((f >= 0.24) & (f <= 0.251))[0]
ind3=np.where((f >= 0.257) & (f <= 0.271))[0]
ind4=np.where((f >= 0.370) & (f <= 0.386))[0]
ind5=np.where((f >= 0.388) & (f <= 0.404))[0]
ind6=np.where((f >= 0.502) & (f <= 0.515))[0]
ind7=np.where((f >= 0.522) & (f <= 0.537))[0]

#fitting function, a 4th degree polynomial
def funcpara(x, a, b, c, d, e):
	return a*x**4+b*x**3+c*x**2+d*x+e

pars=[] #fit parameters

detinds=np.array([2,3,4,5,7,8]) #indices of the detectors we care about

for i in range(0,len(detinds)):
	popt, pcov = curve_fit(funcpara, f[ind1], np.log10(f2[detinds[i],ind1])) #curve fit
	mi=fmin(lambda x,a,b,c,d,e: -funcpara(x, a, b, c, d, e), 0.13,args=(popt[0], popt[1], popt[2],popt[3],popt[4]), xtol=0.0001) #find maximum
	pars=np.append(pars,detinds[i]) #write
	pars=np.append(pars,1) #write
	pars=np.append(pars,popt) #write
	pars=np.append(pars,mi) #write
	popt, pcov = curve_fit(funcpara, f[ind2], np.log10(f2[detinds[i],ind2])) #repeat for each region of interest
	mi=fmin(lambda x,a,b,c,d,e: -funcpara(x, a, b, c, d, e), 0.245,args=(popt[0], popt[1], popt[2],popt[3],popt[4]), xtol=0.0001)
	pars=np.append(pars,detinds[i])
	pars=np.append(pars,2)
	pars=np.append(pars,popt)
	pars=np.append(pars,mi)
	popt, pcov = curve_fit(funcpara, f[ind3], np.log10(f2[detinds[i],ind3]))
	mi=fmin(lambda x,a,b,c,d,e: -funcpara(x, a, b, c, d, e), 0.265,args=(popt[0], popt[1], popt[2],popt[3],popt[4]), xtol=0.0001)
	pars=np.append(pars,detinds[i])
	pars=np.append(pars,3)
	pars=np.append(pars,popt)
	pars=np.append(pars,mi)
	popt, pcov = curve_fit(funcpara, f[ind4], np.log10(f2[detinds[i],ind4]))
	mi=fmin(lambda x,a,b,c,d,e: -funcpara(x, a, b, c, d, e), 0.377,args=(popt[0], popt[1], popt[2],popt[3],popt[4]), xtol=0.0001)
	pars=np.append(pars,detinds[i])
	pars=np.append(pars,4)
	pars=np.append(pars,popt)
	pars=np.append(pars,mi)
	popt, pcov = curve_fit(funcpara, f[ind5], np.log10(f2[detinds[i],ind5]))
	mi=fmin(lambda x,a,b,c,d,e: -funcpara(x, a, b, c, d, e), 0.395,args=(popt[0], popt[1], popt[2],popt[3],popt[4]), xtol=0.0001)
	pars=np.append(pars,detinds[i])
	pars=np.append(pars,5)
	pars=np.append(pars,popt)
	pars=np.append(pars,mi)
	popt, pcov = curve_fit(funcpara, f[ind6], np.log10(f2[detinds[i],ind6]))
	mi=fmin(lambda x,a,b,c,d,e: -funcpara(x, a, b, c, d, e), 0.51,args=(popt[0], popt[1], popt[2],popt[3],popt[4]), xtol=0.0001)
	pars=np.append(pars,detinds[i])
	pars=np.append(pars,6)
	pars=np.append(pars,popt)
	pars=np.append(pars,mi)
	popt, pcov = curve_fit(funcpara, f[ind7], np.log10(f2[detinds[i],ind7]))
	mi=fmin(lambda x,a,b,c,d,e: -funcpara(x, a, b, c, d, e), 0.53,args=(popt[0], popt[1], popt[2],popt[3],popt[4]), xtol=0.0001)
	pars=np.append(pars,detinds[i])
	pars=np.append(pars,7)
	pars=np.append(pars,popt)
	pars=np.append(pars,mi)
	#plt.plot(f[ind1], np.log10(f2[detinds[i],ind1]))
	#plt.plot(f[ind1],funcpara(f[ind1], popt[0], popt[1], popt[2],popt[3],popt[4]))
	#plt.plot([mi,mi],[0,10])
	#plt.show()

parmat=np.reshape(pars,(42,8))

pfreqs=parmat[:,7]
wp=2*np.pi*pfreqs
pdets=parmat[:,0]

pft=[]

#finding the phases at each peak
for i in range(0,42):
	tdet=t[np.where(detind == int(pdets[i]+9))]
	fcs=tdet*wp[i]
	fcs2=np.exp(-1j*fcs)
	su=np.sum(fcs2,0)
	pft=np.append(pft,su)

phas=np.angle(pft)
parmat2=np.concatenate([parmat,phas],1)
parmat2=np.zeros((42,9))
parmat2[:,0:8]=parmat
parmat2[:,8]=phas
parmat2[:,[8,7]] = parmat2[:,[7,8]]

#By column: (detector #)-1, frequency region, fit parameters (5 columns), phase at maximum, maximum frequency
np.savetxt('fftfitpars.dat',parmat2,fmt='%01d %01d %+.5e %+.5e %+.5e %+.5e %+.5e %+1.6f %1.6f')

#finding difference in maxima (3 and 5, 4 and 6, 8 and 9)
p1a=parmat2[0:7,:]
p2a=parmat2[7:14,:]
p1b=parmat2[14:21,:]
p2b=parmat2[21:28,:]
p3a=parmat2[28:35,:]
p3b=parmat2[35:42,:]

p1id=p1a[:,0:2]
p2id=p2a[:,0:2]
p3id=p3a[:,0:2]

p1=p1a[:,7:9]-p1b[:,7:9]
p2=p2a[:,7:9]-p2b[:,7:9]
p3=p3a[:,7:9]-p3b[:,7:9]

pardif1=np.concatenate([p1,p2],0)
pardif2=np.concatenate([pardif1,p3],0)
pardif1id=np.concatenate([p1id,p2id],0)
pardif2id=np.concatenate([pardif1id,p3id],0)
pardif=np.concatenate([pardif2,pardif2id],1)
pardif[:,[0,1,2,3]] = pardif[:,[2,3,0,1]]
pardif[:,3]=1000*pardif[:,3]

#By column: (first detector #)-1, frequency region, difference in phase at maximum, difference in maximum frequency
np.savetxt('fftfitdifpars.dat',pardif,fmt='%01d %01d %+1.6f %+1.3f')
