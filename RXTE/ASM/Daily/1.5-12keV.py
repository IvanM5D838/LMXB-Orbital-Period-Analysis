import matplotlib.pyplot as plt
import numpy as np
import statistics
from astropy.timeseries import LombScargle
import astropy.units as u
import lightkurve
import math

RXTE_PERIOD = 1 / 1.4024762381190594 #confirmed

def cut_data(time, count, error, start_time, end_time):
    start_idx = next(i for i, t in enumerate(time) if t >= start_time)
    end_idx = next(i for i, t in enumerate(time) if t > end_time)
    return time[start_idx:end_idx], count[start_idx:end_idx], error[start_idx:end_idx]

# ---------------------------------------------------------------------------

def data_red(timeListc,countRateListc,errorListc):
    countRateList = list(countRateListc)
    timeList = list(timeListc)
    errorList = list(errorListc)

    preDelList=[]

    for i in range(0,len(timeList)):
        # method 1
        if countRateList[i] == 0  or errorList[i] == 0: 
            preDelList.append(i)

    preDelList = sorted(preDelList,reverse=True)
    #pop data in preDelList
    for i in preDelList:
        countRateList.pop(i)
        timeList.pop(i)
        errorList.pop(i)
    return np.array(timeList),np.array(countRateList),np.array(errorList)

# ---------------------------------------------------------------------------

with open(f"D:/For Institute/NCU/高能天文實驗室/vscode/X2127+119/RXTE/ASM/Daily/daily_lc.tsv", "r") as fobj:
    time, count, error= np.loadtxt(fobj, skiprows=0, usecols=(0, 1, 2), unpack=True)

TimeList_clear,CountRateList_clear,ErrorList_clear = data_red(time,count,error)

cut_time, cut_count, cut_err = cut_data(TimeList_clear, CountRateList_clear, ErrorList_clear, 54500, 55500)
cut_time1, cut_count1, cut_err1 = cut_data(TimeList_clear, CountRateList_clear, ErrorList_clear, 50500, 50600)
cut_time2, cut_count2, cut_err2 = cut_data(TimeList_clear, CountRateList_clear, ErrorList_clear, 51200, 51250)
cut_time3, cut_count3, cut_err3 = cut_data(TimeList_clear, CountRateList_clear, ErrorList_clear, 51550, 51650)
cut_time4, cut_count4, cut_err4 = cut_data(TimeList_clear, CountRateList_clear, ErrorList_clear, 54225, 54325)
cut_time5, cut_count5, cut_err5 = cut_data(TimeList_clear, CountRateList_clear, ErrorList_clear, 54975, 55075)

def group_data(time, count_rate, interval=20):
    """
    将观测时间按照指定间隔分组，并计算每组的平均计数率。

    参数：
    - time: 观测时间的数组
    - count_rate: 计数率的数组，与观测时间对应
    - interval: 观测时间的间隔，默认为20（单位：天）

    返回值：
    - grouped_time: 每个分组的中间时间的数组
    - grouped_count_rate: 每个分组的平均计数率列表
    """
    grouped_time = []
    grouped_count_rate = []

    grouped_count = []
    start_index = 0
    for i, rate in enumerate(count_rate):
        grouped_count.append(rate)
        if (i + 1 - start_index) % interval == 0:  # 每interval个数据点为一组
            group_mid_time = (time[start_index] + time[i]) / 2  # 计算中间时间
            grouped_time.append(group_mid_time)
            grouped_count_rate.append(np.mean(grouped_count))
            grouped_count = []
            start_index = i + 1

    # 处理剩余的数据
    if start_index < len(time):
        group_mid_time = (time[start_index] + time[-1]) / 2
        grouped_time.append(group_mid_time)
        grouped_count_rate.append(np.mean(grouped_count))

    return grouped_time, grouped_count_rate

# 使用函数进行数据分组
grouped_time, grouped_count_rate = group_data(cut_time, cut_count, interval=20)

# # 绘制演变图
# plt.figure(figsize=(10, 6))
# plt.scatter(grouped_time, grouped_count_rate, color='blue', s=15, label='20day bin', alpha=0.6)
# plt.scatter(cut_time, cut_count, color='r', s=6, label='1day bin', alpha=0.3)
# plt.legend(loc='upper right')
# plt.xlabel('Time')
# plt.ylabel('Count Rate')
# plt.title('Count Rate Evolution with Time Binned (20-day intervals)(RXTE ASM)')
# plt.show()

# plt.scatter(grouped_time, grouped_count_rate, color='blue', s=15, label='20day bin', alpha=0.6)
# plt.scatter(cut_time, cut_count, color='black', s=6, label='1day bin')
# plt.scatter(cut_time1, cut_count1, color='r', s=6)
# plt.scatter(cut_time2, cut_count2, color='r', s=6)
# plt.scatter(cut_time3, cut_count3, color='r', s=6)
# plt.scatter(cut_time4, cut_count4, color='r', s=6)
# plt.scatter(cut_time5, cut_count5, color='r', s=6)
# plt.legend(loc='upper right')
# plt.xlabel('Time')
# plt.ylabel('Count Rate')
# plt.title('flare time of RXTE')
# plt.show()

# 将这些时间段内的数据从原始数据中删除
for cut_time, cut_count, cut_err in [(cut_time1, cut_count1, cut_err1),
                                     (cut_time2, cut_count2, cut_err2),
                                     (cut_time3, cut_count3, cut_err3),
                                     (cut_time4, cut_count4, cut_err4),
                                     (cut_time5, cut_count5, cut_err5)]:
    # 找到要删除的时间段在原始数据中的索引范围
    indices = np.where(np.logical_and(TimeList_clear >= cut_time[0], TimeList_clear <= cut_time[-1]))[0]
    
    # 删除这些索引对应的数据
    TimeList_clear = np.delete(TimeList_clear, indices)
    CountRateList_clear = np.delete(CountRateList_clear, indices)
    ErrorList_clear = np.delete(ErrorList_clear, indices)

# plt.scatter(TimeList_clear, CountRateList_clear, s=6)
# plt.xlabel('Time')
# plt.ylabel('Count')
# plt.title('Light Curve(RXTE)')
# plt.legend()
# plt.show()

# ---------------------------------------------------------------------------

def LS_power(time=list, count=list):
    # frequency = cycles / year
    frequency = np.linspace(1.38, 1.41, 2000)
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

Frequency,Power = LS_power(TimeList_clear,CountRateList_clear)  

# plt.plot(Frequency,Power,'b',label='1.5-12keV') 
# plt.xlabel('frequency(cycle/year)')
# plt.ylabel('power')
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

Weighted_Average,Phase,Error = fold_light_curve(CountRateList_clear, TimeList_clear, ErrorList_clear,RXTE_PERIOD,num_bins=32)

plt.errorbar(Phase,Weighted_Average,yerr = Error,fmt = '.',label = f'RXTE_PERIOD = {RXTE_PERIOD}(day)',color = 'r')
plt.title('fold light curve(RXTE)')
plt.legend(loc = 'upper right')
plt.xlabel('Phase')
plt.ylabel('Weighted_Average_Countrate')
# plt.ylim(2.3,3.6)
plt.show()  

