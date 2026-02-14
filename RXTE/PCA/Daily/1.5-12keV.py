import matplotlib.pyplot as plt
import numpy as np
import statistics
from astropy.timeseries import LombScargle
import astropy.units as u
import lightkurve
import math
import csv

MAXI_PERIOD =  1 / 0.013842 #confirmed
RXTE_PERIOD = 1 / 0.005498 #confirmed

def div(countlist,timelist,errorlist,start):
    time = []
    count = []
    err = []
    for i in range(len(timelist)):
        if start <= timelist[i] - timelist[0] <= start + 1000:
            time.append(timelist[i])
            count.append(countlist[i])
            err.append(errorlist[i])
        else:
            pass
    return time,count,err

# ---------------------------------------------------------------------------

def data_red(timeListc,countRateListc,errorListc):
    countRateList = list(countRateListc)
    timeList = list(timeListc)
    errorList = list(errorListc)

    preDelList=[]

    for i in range(0,len(timeList)):
        # method 1
        if countRateList[i] < 0 or abs(countRateList[i] / errorList[i]) < 2:
            preDelList.append(i)
        
    preDelList = sorted(preDelList,reverse=True)
    #pop data in preDelList
    for i in preDelList:
        countRateList.pop(i)
        timeList.pop(i)
        errorList.pop(i)
    return timeList,countRateList,errorList

# ---------------------------------------------------------------------------

with open(f"D:/For Institute/NCU/高能天文實驗室/vscode/X2127+119/RXTE/ASM/Daily/daily_bin.tsv", "r") as fobj:
    time, count, error= np.loadtxt(fobj, skiprows=1, usecols=(0, 2, 3), unpack=True)

print(len(time))

TimeList_clear,CountRateList_clear,ErrorList_clear = data_red(time,count,error)

# plt.plot(TimeList_clear,CountRateList_clear)
# plt.legend(loc='upper right')
# plt.xlabel('MJDcenter')
# plt.ylabel('countrate')
# plt.show()

def LS_power(time=list, count=list):
    # frequency = cycles / year
    frequency = np.linspace(0.001, 1, 2000)
    power = np.zeros_like(frequency)
    N = len(time)
    h = count - np.mean(count)
    sigma_square = 1 / (N - 1) * np.sum((count - np.mean(count))**2)
    mean_freq = [] 

    for i, freq in enumerate(frequency):
        omega = 2 * np.pi * freq
        tau = np.arctan(np.sum(np.sin(2 * omega * np.array(time)))/np.sum(np.cos(2 * omega * np.array(time))))/ (2 * omega)
        t = time - tau
        sin = np.sum(h * np.sin(omega * t))
        cos = np.sum(h * np.cos(omega * t))
        cos2 = np.sum(np.cos(omega * t)**2)
        sin2 = np.sum(np.sin(omega * t)**2)
        power[i] = (cos**2 / cos2 + sin**2 / sin2) / (2 * sigma_square) 
    # for i in range(len(power)):
    #     if power[i] == np.max(power):
    #         print(frequency[i])
    return frequency * 365.25, power

# ----------------------------------------------------------------------- 

Frequency,Power = LS_power(TimeList_clear,CountRateList_clear)  

plt.plot(Frequency,Power,'b',label='1.5-12keV') 
plt.xlabel('frequency(cycle/year)')
plt.ylabel('power')
plt.title('LS_power spectrum(RXTE)')
plt.legend(loc='best')
plt.show()

# ---------------------------------------------------------------------------

def dynamic_ls_periodogram(t = list,ct = list,win = int,mov = int):
    index = []
    delta_t = t[-1] - t[0]
    n = math.floor((delta_t-win)/mov)
    tc = np.zeros(n)
    t2 = []
    ct2 = []
    dps = []
    #dps=dynamic power spectrum
    for i in range(0,n):
        tc[i] = np.array(int(t[0] + win * 0.5 + mov * i))
        tmax = np.array(t[0] + win + mov * i)
        tmin = np.array(t[0] + mov * i)
        index.append(list(np.where((t >= tmin) & (t < tmax))[0]))
        t2.append(list([t[index[i][k]] for k in range(0,len(index[i]))]))
        ct2.append(list([ct[index[i][k]] for k in range(0,len(index[i]))]))
        dps.append(LS_power(t2[i],ct2[i])[1])
    dps = np.rot90(dps)
    tc = tc - 50000
    return tc,dps

# TC,DPS = dynamic_ls_periodogram(TimeList_clear,CountRateList_clear,win = 1000,mov = 10)

# plt.imshow(DPS,aspect = 'auto',extent = [np.min(TC),np.max(TC),np.min(Frequency),np.max(Frequency)],cmap='jet')
# plt.xlabel(r'Time(MJD-50000)',fontsize = 12)
# plt.ylabel(r'Frequency(cycles/year)',fontsize = 12)
# plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(1))
# plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(1))
# plt.colorbar()
# plt.title('dynamic power spectrum(RXTE(1.5-12keV))')
# plt.show()

# -----------------------------------------------------------------------------

def fold_light_curve(countrate, time, error, period, num_bins):
    time = np.array(time)
    countrate = np.array(countrate)
    error = np.array(error)
    
    # 計算相位，將其範圍調整為0到2
    phase = np.mod(time - time[0], period) / period
    
    bin_width = 1 / num_bins  # binned_time點之間的間隔
    folded_countrate = np.zeros(num_bins)
    bin_weights = np.zeros(num_bins)  # binned_phase加權
    bin_errors = np.zeros(num_bins)  # binned_time誤差

    for i in range(len(phase)):
        bin_index = int(phase[i] // bin_width)
        weight = 1 / (error[i] ** 2)
        folded_countrate[bin_index] += countrate[i] * weight
        bin_weights[bin_index] += weight
        bin_errors[bin_index] += (error[i] ** -2)

    # 計算加權平均值
    folded_countrate /= bin_weights
    bin_errors = np.sqrt(1 / bin_errors)
    
    # 將相位0-1複製到1-2上
    folded_countrate_0_to_2 = np.concatenate((folded_countrate, folded_countrate))
    bin_errors_0_to_2 = np.concatenate((bin_errors, bin_errors))
    
    # 調整相位值
    phase_0_to_2 = np.linspace(bin_width / 2, 2 - bin_width / 2, 2 * num_bins)

    return folded_countrate_0_to_2, phase_0_to_2, bin_errors_0_to_2

# Weighted_Average,Phase,Error = fold_light_curve(CountRateList_clear,TimeList_clear,ErrorList_clear,RXTE_PERIOD,num_bins=32)

# plt.errorbar(Phase,Weighted_Average,yerr = Error,fmt = '^',label = 'RXTE_PERIOD = 182.083(day)',color = 'r')
# plt.title('fold light curve(RXTE)')
# plt.legend(loc = 'upper right')
# plt.xlabel('Phase')
# plt.ylabel('Weighted_Average_Countrate')
# # plt.ylim(2.3,3.6)
# plt.show()  

