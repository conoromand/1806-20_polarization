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

#Calculates parameters for fits for data phased on other time intervals

datr=np.loadtxt('RHESSI_dat')
#datp=np.loadtxt('fitparsfsat.dat')

#dat[:,0] is arrival time
#dat[:,1] is energy in keV
#dat[:,2] is detector number

psat=4.07
pstar=7.54773
bins=400
ind=np.where((datr[:,0] > 32) & (datr[:,0] < 132)) # We only look at these 100 s
datc=datr[ind[0],:]
dat=datc

#--------------------------------------------------------------------------------------2*fsat

rp=psat/2.
dat[:,0]=np.remainder(datc[:,0], rp)

#cums=[]
pars=[]

def func2fsat(x, a, b, c):
	return a*np.cos(2*np.pi/rp*x-b)+c

for i in range(0,18):
	dat9=dat[np.where(dat[:,2] == i),:2][0]
	his9=np.histogram(dat9, bins=bins, range=(0,rp))[0]
#	his=(his9-datp[i,2])/datp[i,0]
	t=np.histogram(dat9, bins=bins, range=(0,rp))[1]
	t2=(t[1:] + t[:-1])/2	
	popt, pcov = curve_fit(func2fsat, t2, his9, bounds=(0., [400., 2*np.pi, 500.]))	
	pars=np.append(pars,popt)
#popt, pcov = curve_fit(func, t2, his9, bounds=(0., [200., 2*np.pi, 500.]))
#	locum = np.cumsum(his)
#	cums=np.append(cums,locum)
#pars=np.append(pars,popt)
	#if i < 20:
		#plt.figure()
		#plt.plot(t2,his, label=str(i))
		#plt.plot(t2,func(t2,popt[0],popt[1],popt[2]), label='fit')
		#plt.legend(loc='upper left')
#plt.text(0, 20, 'Parameters = ' + str(popt), fontsize=15)
#plt.plot(lobase[:-1], locum, label=str(i))

#plt.legend(loc='upper left')
#plt.plot([0,4.07],[0,1],'k--')
#plt.show()
#cums=np.reshape(cums, (9, bins))
#plt.figure()

#for i in range(0,9):
#	plt.plot(t2,cums[i,:],label=str(i+9))

#plt.legend(loc='lower left')
#plt.show()

pars=np.reshape(pars, (18, 3))
np.savetxt('fitpars2fsat.dat',pars,fmt='%-.6f')

#--------------------------------------------------------------------------------------fsat+fstar

rp=1/(1/psat+1/pstar)
dat[:,0]=np.remainder(datc[:,0], rp)

pars=[]

def func(x, a, b, c):
	return a*np.cos(2*np.pi/rp*x-b)+c

for i in range(0,18):
	dat9=dat[np.where(dat[:,2] == i),:2][0]
	his9=np.histogram(dat9, bins=bins, range=(0,rp))[0]
	t=np.histogram(dat9, bins=bins, range=(0,rp))[1]
	t2=(t[1:] + t[:-1])/2.	
	popt, pcov = curve_fit(func, t2, his9, bounds=(0., [400., 2*np.pi, 500.]))	
	pars=np.append(pars,popt)

pars=np.reshape(pars, (18, 3))
np.savetxt('fitparsfsat+fstar.dat',pars,fmt='%-.6f')

#--------------------------------------------------------------------------------------2fsat+fstar

rp=1/(1/(2*psat)+1/pstar)
dat[:,0]=np.remainder(datc[:,0], rp)

pars=[]

def func(x, a, b, c):
	return a*np.cos(2*np.pi/rp*x-b)+c

for i in range(0,18):
	dat9=dat[np.where(dat[:,2] == i),:2][0]
	his9=np.histogram(dat9, bins=bins, range=(0,rp))[0]
	t=np.histogram(dat9, bins=bins, range=(0,rp))[1]
	t2=(t[1:] + t[:-1])/2.	
	popt, pcov = curve_fit(func, t2, his9, bounds=(0., [400., 2*np.pi, 500.]))	
	pars=np.append(pars,popt)

pars=np.reshape(pars, (18, 3))
np.savetxt('fitpars2fsat+fstar.dat',pars,fmt='%-.6f')

#--------------------------------------------------------------------------------------2fsat+2fstar

rp=1/(1/(2*psat)+1/(2*pstar))
dat[:,0]=np.remainder(datc[:,0], rp)

pars=[]

def func(x, a, b, c):
	return a*np.cos(2*np.pi/rp*x-b)+c

for i in range(0,18):
	dat9=dat[np.where(dat[:,2] == i),:2][0]
	his9=np.histogram(dat9, bins=bins, range=(0,rp))[0]
	t=np.histogram(dat9, bins=bins, range=(0,rp))[1]
	t2=(t[1:] + t[:-1])/2.	
	popt, pcov = curve_fit(func, t2, his9, bounds=(0., [400., 2*np.pi, 500.]))	
	pars=np.append(pars,popt)

pars=np.reshape(pars, (18, 3))
np.savetxt('fitpars2fsat+2fstar.dat',pars,fmt='%-.6f')

#--------------------------------------------------------------------------------------fstar

rp=pstar
dat[:,0]=np.remainder(datc[:,0], rp)

pars=[]

def func(x, a, b, c):
	return a*np.cos(2*np.pi/rp*x-b)+c

for i in range(0,18):
	dat9=dat[np.where(dat[:,2] == i),:2][0]
	his9=np.histogram(dat9, bins=bins, range=(0,rp))[0]
	t=np.histogram(dat9, bins=bins, range=(0,rp))[1]
	t2=(t[1:] + t[:-1])/2.	
	popt, pcov = curve_fit(func, t2, his9, bounds=(0., [400., 2*np.pi, 500.]))	
	pars=np.append(pars,popt)
#	plt.plot(t2, his9)
#	plt.show()

pars=np.reshape(pars, (18, 3))
np.savetxt('fitparsfstar.dat',pars,fmt='%-.6f')

#--------------------------------------------------------------------------------------2fstar

rp=pstar/2.
dat[:,0]=np.remainder(datc[:,0], rp)

pars=[]

def func(x, a, b, c):
	return a*np.cos(2*np.pi/rp*x-b)+c

for i in range(0,18):
	dat9=dat[np.where(dat[:,2] == i),:2][0]
	his9=np.histogram(dat9, bins=bins, range=(0,rp))[0]
	t=np.histogram(dat9, bins=bins, range=(0,rp))[1]
	t2=(t[1:] + t[:-1])/2.	
	popt, pcov = curve_fit(func, t2, his9, bounds=(0., [400., 2*np.pi, 500.]))	
	pars=np.append(pars,popt)

pars=np.reshape(pars, (18, 3))
np.savetxt('fitpars2fstar.dat',pars,fmt='%-.6f')
