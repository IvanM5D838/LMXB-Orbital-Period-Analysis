
from cmath import phase
from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
def gen_peroidic_data(x,period=5.25,phase=0,amplitude=1,noise=0):
    y=amplitude*np.sin(2*np.pi*x/period-phase)
    dy=np.random.normal(0,np.sqrt(noise),size=len(y))
    return y+dy
#x=np.linspace(0,10,11)
#y=gen_peroidic_data(x,amplitude=2,period=1/0.7)
#x_signal=np.linspace(0,10,100)
#y_signal=gen_peroidic_data(x_signal,amplitude=2,period=1/0.7)
#fig, ax=plt.subplots(figsize=(8,4))
#ax.scatter(x,y)
#ax.plot(x_signal,y_signal)
 
#y_high=gen_peroidic_data(x,amplitude=2,period=1/2.7)
#y_signal_high=gen_peroidic_data(x_signal,amplitude=2,period=1/2.7)
#ax.scatter(x,y_high)
#ax.plot(x_signal,y_signal_high)    
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#fig.tight_layout()
#plt.show()
    

def chi2(theta,y,y_unc,x,f):
    a=theta[0]
    phi=theta[1]
    chi2=np.sum((y-a*np.sin(2*np.pi*x/f-phi))**2/y_unc**2)
    return chi2

def min_chi2(theta,y,y_unc,x,f):
    res=minimize(chi2,theta,args=(y,y_unc,x,f))
    return res.fun

def ls_periodogram(y,y_unc,x,f_grid):
    psd=np.empty_like(f_grid)
    chi2_0=np.sum(((y-np.mean(y))/y_unc)**2)
    for f_num,f in enumerate(f_grid):
        psd[f_num]=0.5*(chi2_0-min_chi2([0,0],y,y_unc,x,f))
    return psd

np.random.seed(185)
x=10*np.random.rand(100)
y=gen_peroidic_data(x,period=0.1,amplitude=7.4,noise=0.8)
y_unc=np.ones_like(x)*np.sqrt(0.8)
f_grid=np.linspace(1/np.ptp(x),10,10000)
psd_ls=ls_periodogram(y,y_unc,x,f_grid)
fig,ax=plt.subplots()
ax.plot(1/f_grid,psd_ls)    
ax.set_xlabel('Period')
ax.set_ylabel('P')
plt.tight_layout()
plt.show()









