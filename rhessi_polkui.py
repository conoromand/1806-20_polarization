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

#Getting the parameters for the data phased on the satellite rotation period

datr=np.loadtxt('RHESSI_dat')

#dat[:,0] is arrival time
#dat[:,1] is energy in keV
#dat[:,2] is detector number

rp=4.07 
bins=400
ind=np.where((datr[:,0] > 32) & (datr[:,0] < 132)) # We only look at these 100 seconds
datc=datr[ind[0],:]
dat=datc
dat[:,0]=np.remainder(datc[:,0], rp)
#t=np.linspace(rp/bins,rp,bins)

def func(x, a, b, c):
	return a*np.cos(2*np.pi/rp*x-b)+c

cums=[] #Will become cumulative distributions for eacch
pars=[] #will be fit paramteres for each
#plt.figure()
detnum=9 # number of detector we're interested in

#Fitting the data to a cos curve

for i in range(9,18): #Range of detector indices 
	dat9=dat[np.where(dat[:,2] == i),:2][0]
	his9=np.histogram(dat9, bins=bins, range=(0,rp))[0]
	t=np.histogram(dat9, bins=bins, range=(0,rp))[1]
	t2=(t[1:] + t[:-1])/2
	popt, pcov = curve_fit(func, t2, his9, bounds=(0., [400., 2*np.pi, 500.]))
	lovalues, lobase = np.histogram(dat9[:,0], bins=bins)
	locumulative = np.cumsum(lovalues)
	locum=locumulative/np.max(locumulative)
	cums=np.append(cums,locum)
	pars=np.append(pars,popt)
	if i < 20:
		plt.figure()
		plt.plot(t2-(popt[1]/(2*np.pi)*rp),his9, label=str(i))#+(popt[1]*(2*np.pi)/rp)
		plt.plot(t2,func(t2,popt[0],0,popt[2]), label='fit')
		plt.legend(loc='upper left')
#plt.text(0, 20, 'Parameters = ' + str(popt), fontsize=15)
#plt.plot(lobase[:-1], locum, label=str(i))

#plt.legend(loc='upper left')
#plt.plot([0,4.07],[0,1],'k--')
plt.show()

#reshape and save
cums=np.reshape(cums, (detnum, bins))
pars=np.reshape(pars, (detnum, 3))
np.savetxt('fitparsfsat.dat',pars,fmt='%-.6f')

pvals=[] #will become Kuiper p-values

#Do the Kuiper test - this is for all 18 detectors right now (code will fail)
for i in range(0,18):
	for j in range(0,18):
		st, p=kp.kuiper_two(cums[i,:],cums[j,:])
		pvals=np.append(pvals,p)	

#Save Kuiper values
pvals=np.reshape(pvals, (18, 18))
pvals[pvals != pvals]=1
np.savetxt('polkuivals.dat',pvals,fmt='%-.6f')
