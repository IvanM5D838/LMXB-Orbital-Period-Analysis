import matplotlib.pyplot as plt
import numpy as np
import statistics
from astropy.timeseries import LombScargle
import astropy.units as u
import lightkurve
import math

START = 0
DAYS = 6000
RXTE_PERIOD = 1 / 1.4191345672836417

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
    errorListStd = np.std(errorList,ddof=1)
    errorListAvg = np.mean(errorList)
    countRateListAvg = np.mean(countRateList)
    countRateListStd = np.std(countRateList,ddof=1)

    preDelList=[]

    for i in range(0,len(timeList)):
        # method 1
        if countRateList[i] < 0 or abs(countRateList[i] / errorList[i]) < 2:
            preDelList.append(i)
        # method 2
        # if abs(countRateList[i] - countRateListAvg) > 3 * countRateListStd or  abs(errorList[i] - errorListAvg) > 3 * errorListStd:
        #     preDelList.append(i)
        # method 3
        # if abs(errorList[i] - errorListAvg) > 3 * errorListStd:
        #     preDelList.append(i)
        
    preDelList = sorted(preDelList,reverse=True)
    #pop data in preDelList
    for i in preDelList:
        countRateList.pop(i)
        timeList.pop(i)
        errorList.pop(i)
    return timeList,countRateList,errorList

# ---------------------------------------------------------------------------

with open(f"D:/For Institute/NCU/高能天文實驗室/vscode/X2127+119/RXTE/ASM/Daily/daily_lc.tsv", "r") as fobj:
    time, count, error= np.loadtxt(fobj, skiprows=0, usecols=(0, 5, 6), unpack=True)

TimeList_clear,CountRateList_clear,ErrorList_clear = data_red(time,count,error)

print(len(TimeList_clear))

# plt.plot(TimeList_clear,CountRateList_clear)
# plt.legend(loc='upper right')
# plt.xlabel('MJDcenter')
# plt.ylabel('countrate')
# plt.show()

# plt.loglog(CountRateList_clear,ErrorList_clear,'.') #　可利用log圖觀察誤差極小的數據點並篩選出來
# plt.xlim(10**-3,10**2)
# plt.ylim(10**-3,10**2)
# plt.legend(loc='upper right')
# plt.xlabel('log(countrate)')
# plt.ylabel('log(error)')
# plt.title('log-log diagram(3-5keV,RXTE)')
# plt.show()

# ------------------------------------------------------------------------

# 利用histogram看誤差的分布並將誤差極大極小值列入篩選範圍內
# a=np.array(ErrorList_clear)
# b=np.array(CountRateList_clear)
# # x=[plt.hist(a,bins=len(ErrorList1))[0]]
# # print(x)
# plt.hist(b,bins=100)
# plt.xlim(-5,5)
# plt.title("histogram(3-5keV,RXTE)") 
# plt.show()

# ---------------------------------------------------------------------

def LS_power(time, count):
    # frequency = cycles / year
    frequency = np.linspace(1.35, 1.45, 2000)
    power = np.zeros_like(frequency)
    N = len(time)
    h = count - np.mean(count)
    sigma_square = 1 / (N - 1) * np.sum((count - np.mean(count))**2)

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

Frequency,Power = LS_power(time,count)  

plt.plot(Frequency,Power,'b',label='3-5keV') 
plt.xlabel('frequency(cycle/year)')
plt.ylabel('power')
# plt.annotate('1.003d', xy=(364.245, 30), xytext=(364.245, 34),
#              arrowprops=dict(facecolor='black', width=1, headwidth=4, headlength=3),
#              fontsize=10, ha='center')
plt.title('LS_power spectrum(RXTE ASM)')
# plt.ylim(0,40)
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
# plt.ylabel(r'Frequency(cycles/yr)',fontsize = 12)
# plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(1))
# plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(1))
# plt.colorbar()
# plt.title('dynamic power spectrum(RXTE(3-5keV))')
# plt.show()

# -----------------------------------------------------------------------------

def fold_light_curve(countrate, time, error, period, num_bins):
    
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

Weighted_Average,Phase,Error = fold_light_curve(count,time,error,RXTE_PERIOD,num_bins=32)

plt.errorbar(Phase,Weighted_Average,yerr = Error,fmt = '.',color = 'r',label = f'period = {format(RXTE_PERIOD,".3f")}d')
plt.title('fold light curve(RXTE ASM)')
plt.legend(loc = 'upper right')
plt.xlabel('Phase')
plt.ylabel('Countrate')
# plt.ylim(15,40)
plt.show()  
