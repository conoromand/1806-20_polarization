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

# Note: sample_rhessi.py for the prototype for this, it's not useful
# Generate sample data based on fits to the detectors

datr=np.loadtxt('RHESSI_dat')
parsat=np.loadtxt('fitparsfsat.dat')
par2sat=np.loadtxt('fitpars2fsat.dat')

psat=4.07
pstar=7.54773
fsat=1/psat
fstar=1/pstar
wsat=2*np.pi*fsat
wstar=2*np.pi*fstar
bins=100000

#Look only at 100 s and 9 detectors
ind=np.where((datr[:,0] > 32) & (datr[:,0] < 132))
datc=datr[ind[0],:]
dat1=datc[np.where((datc[:,2] < 17.5) & (datc[:,2] > 8.5)),:][0]
detind=dat1[:,2]

co=np.array([5,3.5,3,2.5,1.9,1.3,1.3]) #coefficients for the pulsar harmonics
so=np.sum(co) #normalizing factor

#satellite response * pulsar emission + noise
def satresp(t, a, b, c, d, phi, phi2):
	return a + (b + c + b*np.cos(wsat*t+phi) + c*np.cos(2*wsat*t+phi) + d*(np.cos(wsat*t+(wstar-phi2)))**2)*(co[0]+co[1]*np.cos(wstar*t)+co[2]*np.cos(2*wstar*t)+co[3]*np.cos(3*wstar*t)+co[4]*np.cos(4*wstar*t)+co[5]*np.cos(5*wstar*t)+co[5]*np.cos(6*wstar*t))/(so)

detid=[11,12,13,14,16,17] #detectors we're interested in
detlen=np.zeros(6)
for i in range(0,len(detid)):
	detlen[i]=len(detind[np.where(detind == detid[i])])

totnum=sum(detlen).astype(int)
seed=np.linspace(1,totnum,totnum)
np.random.seed(seed=seed.astype(int))
rand=np.random.rand(totnum) #seed for the random numbers

coords=[15.5+5j,-15.5+5j, 8+14j, -8+14j,-7.8-15.5j,7.8-15.5j] #coordinates for the detectors
phs=0.25 #possible additional phase difference

#parameters for each detector
a1=parsat[detid,2]-parsat[detid,0]
b1=parsat[detid,0]
c1=par2sat[detid,0]
d1=[0,0,0,0,b1[4]*0.2,b1[5]*0.2]
#phi1=parsat[detid,1]
phi1=np.angle(coords) #phases from geometry of detectors
phi2=[0,0,0,0,(phs)*np.pi,(phs+0.5)*np.pi] #extra phase due to scattering

tl=np.linspace(32,132,bins+1)
tphot=[] #photon arrival times
dets=[] #detector number

for i in range(0,len(detid)):

	#calculate normalized cumulative distributions and arrival functions
	funcr=satresp(tl,a1[i],b1[i],c1[i]*0,d1[i]*0,phi1[i],phi2[i])
	cfuncr = np.cumsum(funcr)	
	cfunc=(cfuncr-cfuncr[0])/(cfuncr[len(cfuncr)-1]-cfuncr[0])
	func=funcr/(cfuncr[len(cfuncr)-1]-cfuncr[0])

	#set up the iterative calculation
	to=((np.ones(detlen[i]))*tl[(bins/2)])
	dif=np.ones(detlen[i])*1000
	tol=0.001
	num=len(np.where(dif > tol)[0])
	count=0
	
	#iteratively calculate photon time by matching the cumulative distribtion
	#go for "count" iterations and until there are less than "num" within the tolerance
	while (num > 10) & (count < 20000):
		idx = np.round(((to-32)/100.*bins),decimals=0).astype(int)
		idx[np.where(idx > bins-1)]=bins-1
		idx[np.where(idx < 0)]=0
		t1=to-(cfunc[idx]-rand[0:detlen[i]])/(func[idx]*1000)
		dif=np.abs(t1-to)
		num=len(np.where(dif > tol)[0])
		to=t1
		count=count+1
		#print(to,num)

	#values, base = np.histogram(satresp(to,a1,b1,c1,d1,e1,f1,phi1), bins=bins)
	#cumulative = np.cumsum(satresp(to,a1[i],b1[i],c1[i],d1[i],phi1[i],phi2[i]))
	#cum=cumulative/np.max(cumulative)

	#plt.plot(tl,cfunc)
	#plt.plot(np.sort(to),cum)
	#plt.xlim(32,132)
	#plt.show()
	delind=np.linspace(0,detlen[i]-1,detlen[i]).astype(int)

	#set up the next detector
	rand=np.delete(rand,delind)
	tphot=np.append(tphot,to)
	newdet=np.ones(detlen[i].astype(int))*detid[i]
	dets=np.append(dets,newdet)

data=np.append(tphot,dets)
dataarr=np.reshape(data,(2,len(tphot)))
dataarr[[0,1],:]=dataarr[[1,0],:] #has arrival times and detector numbers

#np.savetxt('tsample2.dat',dataarr.T,fmt=('%-02d', '%-06.3f'))
#np.savetxt('tsamplenp.dat',dataarr.T,fmt=('%-02d', '%-06.3f'))

#name has the harmonics for the pulsa and the additional scattering phase
np.savetxt('tsamplep' + str(co) + str(phs) + '.dat',dataarr.T,fmt=('%-02d', '%-06.3f'))

