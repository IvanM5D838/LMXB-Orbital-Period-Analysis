import matplotlib.pyplot as plt
import numpy as np
import statistics
from astropy.timeseries import LombScargle
import astropy.units as u
import lightkurve
import math

START = 0
DAYS = 6000
PERIOD = 0

def div(countlist,timelist,errorlist,start):
    time = []
    count = []
    err = []
    for i in range(len(timelist)):
        if start <= timelist[i] - timelist[0] <= start + DAYS:
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
        # if countRateList[i] < 0 or  abs(errorList[i] - errorListAvg) > 3 * errorListStd:
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

file = open('RXTE/ASM/Orbit_lc.txt','r')
timeList = []
countRateList = []
errorList = []
while True:
#一次讀一行
    str = file.readline()
#判斷讀進來的是不是空字串(如果是空字串代表讀到檔案結尾 就直接break)
    if not str:
        break
#把str(這個type是string) 用空白切割成一個一個的string list存到list裡
    Datalist = str.split()
#轉成float的形式丟到x,y陣列裡
    timeList.append(float(Datalist[0]))
    countRateList.append(float(Datalist[3]))
    errorList.append(float(Datalist[4]))

TimeList_clear,CountRateList_clear,ErrorList_clear = data_red(timeList,countRateList,errorList)

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
# plt.title('log-log diagram(1.5-3keV,RXTE)')
# plt.show()

# ------------------------------------------------------------------------

# 利用histogram看誤差的分布並將誤差極大極小值列入篩選範圍內
# a=np.array(ErrorList_clear)
# b=np.array(CountRateList_clear)
# # x=[plt.hist(a,bins=len(ErrorList1))[0]]
# # print(x)
# plt.hist(a,bins=100)
# plt.xlim(0,5)
# plt.title("histogram(1.5-3keV,RXTE)") 
# plt.show()

# ---------------------------------------------------------------------

def LS_power(time=list, count=list):
    frequency = np.linspace(0.001, 0.01, 2000)
    power = np.zeros_like(frequency)
    N = len(time)
    h = count - np.mean(count)
    sigma_square = 1 / (N - 1) * np.sum((count - np.mean(count))**2)
    mean_freq = []
    # print(len(h))

    for i, freq in enumerate(frequency):
        omega = 2 * np.pi * freq
        tau = np.arctan(np.sum(np.sin(2 * omega * np.array(time)))/np.sum(np.cos(2 * omega * np.array(time))))/ (2 * omega) 
        t = time - tau
        sin = np.sum(h * np.sin(omega * t))
        cos = np.sum(h * np.cos(omega * t))
        cos2 = np.sum(np.cos(t)**2)
        sin2 = np.sum(np.sin(t)**2)
        power[i] = (cos**2 / cos2 + sin**2 / sin2) / (2 * sigma_square) 
    # for idx,i in enumerate(power):
    #     if i > 200:
    #        mean_freq.append(frequency[idx])
    # print(np.mean(mean_freq))
    # print(len(power))
    return frequency * 365.25, power
# -----------------------------------------------------------------------

Frequency,Power = LS_power(TimeList_clear,CountRateList_clear)  

# plt.plot(Frequency,Power,'b',label='1.5-3keV') 
# plt.xlabel('frequency(1/day)')
# plt.ylabel('power')
# plt.vlines(x=0.210970464135021, ymin=0, ymax=12, colors='r', linestyles='dashed', label='period = 4.74(day)')
# plt.title('LS_power spectrum(RXTE)')
# plt.legend(loc='best')
# plt.show()

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
# plt.title('dynamic power spectrum(RXTE(1.5-3keV))')
# plt.show()

# -----------------------------------------------------------------------------

# def fold_light_curve(countrate, time, error, period, num_bins, time_start, days):
#     binned_time = []
#     Phase = []
#     for i in range(len(time)):
#         if 0 <= (time[i] - time[0] - time_start) <= days:
#             binned_time.append(time[i])
#         else:
#             continue

#     phase = [binned_time[i] / period for i in range(len(binned_time))]
            
#     for i in range(len(binned_time)):
#         if 0 < phase[i] < 2:
#             Phase.append(phase[i])
#         elif phase[i] >= 2:
#             phase[i] = phase[i] - 2 * int(phase[i] // 2) 
#             Phase.append(phase[i])
#             continue

#     # phase = np.mod(binned_time, period) / period  # 根據固定時間範圍內的點作圖
#     bin_width = 2 / num_bins  # binnend_time點之間的間隔
#     folded_countrate = np.zeros(num_bins)
#     bin_weights = np.zeros(num_bins)  # binned_phase加權
#     bin_errors = np.zeros(num_bins)  # binned_time誤差

#     for i in range(len(Phase)):
#         bin_index = int(Phase[i] // bin_width) 
#         weight = 1 / (error[i] ** 2)
#         folded_countrate[bin_index] += countrate[i] * weight
#         bin_weights[bin_index] += weight
#         bin_errors[bin_index] += (error[i] ** -2)

#     folded_countrate /= bin_weights
#     bin_errors = np.sqrt(1 / bin_errors)
#     bin_centers = np.linspace(bin_width / 2, 2 - bin_width / 2, num_bins)

#     return folded_countrate, bin_centers, bin_errors

# Weighted_Average,Phase,Error = fold_light_curve(CountRateList_clear,TimeList_clear,ErrorList_clear,PERIOD,num_bins=50,time_start=START,days=DAYS)

# plt.errorbar(Phase,Weighted_Average,yerr = Error,fmt = '^',label = 'X-ray dip',color = 'r')
# plt.title(f'fold light curve(swift)({START},{START+DAYS})')
# plt.legend(loc = 'upper right')
# plt.xlabel('Phase')
# plt.ylabel('Weighted_Average_Countrate')
# plt.show()  

import numpy as np
import matplotlib.pyplot as plt

def weighted_average_folded_light_curve(time, count, error, period, num_bins):
    # 將時間按週期折疊
    folded_time = time % period

    # 計算相位
    phase = (folded_time / period) % 1

    # 定義相位的數量和範圍
    bin_edges = np.linspace(0, 2, num_bins + 1)

    # 計算每個相位範圍的加權平均值和標準差
    weighted_averages = []
    errors = []
    for i in range(num_bins):
        mask = (phase >= bin_edges[i]) & (phase < bin_edges[i + 1])
        weighted_average = np.average(count[mask], weights=1 / error[mask]**2)
        weighted_averages.append(weighted_average)
        errors.append(np.sqrt(1 / np.sum(1 / error[mask]**2)))

    # 繪製加權平均的折疊光變曲線與錯誤條
    plt.figure(figsize=(10, 5))
    plt.errorbar(bin_edges[:-1], weighted_averages, yerr=errors, fmt='o-', color='blue')
    plt.title('Weighted Average Folded Light Curve with Error Bars')
    plt.xlabel('Phase')
    plt.ylabel('Weighted Average Flux')
    plt.show()

# 使用函數繪製加權平均的折疊光變曲線與錯誤條
weighted_average_folded_light_curve(TimeList_clear, CountRateList_clear, ErrorList_clear, period=123, num_bins=32)

