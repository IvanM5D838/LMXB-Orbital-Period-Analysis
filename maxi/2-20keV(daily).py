# XMM時間尺度短(s-ms)
# Swift有daily(一天)跟orbit(90min)
# RXTE有PCA(1ms)、ASM(約15min)、HEXTE(8ms)
# MAXI有daily(一天)跟orbit(96min)

import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle
import astropy.units as u
import math
from astropy.io import fits
from matplotlib.animation import FuncAnimation

START=6000
MAXI_PERIOD = 1 / 1.4024762381190594

def cut_data(time, count, error, start_time, end_time):
    start_idx = next(i for i, t in enumerate(time) if t >= start_time)
    end_idx = next(i for i, t in enumerate(time) if t > end_time)
    return time[start_idx:end_idx], count[start_idx:end_idx], error[start_idx:end_idx]


def data_red(timeListc,countRateListc,errorListc):
    countRateList = list(countRateListc)
    timeList = list(timeListc)
    errorList = list(errorListc)

    preDelList=[]

    for i in range(0,len(timeList)):
        # method 1
        if countRateList[i] <= 0 or errorList[i] <= 0:
            preDelList.append(i)

    preDelList = sorted(preDelList,reverse=True)
    #pop data in preDelList
    for i in preDelList:
        countRateList.pop(i)
        timeList.pop(i)
        errorList.pop(i)
    return timeList,countRateList,errorList

# with fits.open('maxi/6h.fits') as hdul:
#     data = hdul[5].data # 假設資料在第一個擴展HDU中

#     # 提取時間、計數和誤差資訊
#     time = data['MJD']
#     count = data['RATE']
#     error = data['RERR']

with open(f"D:/For Institute/NCU/高能天文實驗室/vscode/X2127+119/maxi/daily_lc.tsv", "r") as fobj:
    time, count, error= np.loadtxt(fobj, skiprows=0, usecols=(0, 1, 2), unpack=True)
TimeList_clear,CountRateList_clear,ErrorList_clear = data_red(time,count,error)

cut_time, cut_count, cut_err = cut_data(TimeList_clear, CountRateList_clear, ErrorList_clear, 59000, 60000)
cut_time1, cut_count1, cut_err1 = cut_data(TimeList_clear, CountRateList_clear, ErrorList_clear, 55425, 55500)
cut_time2, cut_count2, cut_err2 = cut_data(TimeList_clear, CountRateList_clear, ErrorList_clear, 55650, 55775)
cut_time3, cut_count3, cut_err3 = cut_data(TimeList_clear, CountRateList_clear, ErrorList_clear, 56500, 56600)
cut_time4, cut_count4, cut_err4 = cut_data(TimeList_clear, CountRateList_clear, ErrorList_clear, 59800, 59900)

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

# 绘制演变图
plt.figure(figsize=(10, 6))
# plt.scatter(grouped_time, grouped_count_rate, color='blue', s=15, label='20day bin', alpha=0.6)
plt.scatter(TimeList_clear, CountRateList_clear, color='black', s=6, label='1day bin')
plt.scatter(cut_time1, cut_count1, color='r', s=6)
plt.scatter(cut_time2, cut_count2, color='r', s=6)
plt.scatter(cut_time3, cut_count3, color='r', s=6)
plt.scatter(cut_time4, cut_count4, color='r', s=6)
plt.legend(loc='upper right')
plt.xlabel('Time')
plt.ylabel('Count Rate')
plt.title('flare time of MAXI')
plt.show()

# -------------------------------------------------------------------------

def LS_power(time=list, count=list):
    # frequency = cycles / year
    frequency = np.linspace(1.39, 1.41, 2000)
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
        cos2 = np.sum(np.cos(t)**2)
        sin2 = np.sum(np.sin(t)**2)
        power[i] = (cos**2 / cos2 + sin**2 / sin2) / (2 * sigma_square) 
    for i in range(len(power)):
        if power[i] == np.max(power):
            print(frequency[i])
    return frequency * 365.25, power

# Frequency,Power = LS_power(TimeList_clear,CountRateList_clear)  

# plt.plot(Frequency,Power,label='2-20 keV') 
# plt.xlabel('frequency(cycle/year)')
# plt.ylabel('power')
# plt.legend(loc = 'best')
# plt.title('LS_power spectrum(maxi)')
# plt.show()  

# -----------------------------------------------------------------------

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

Weighted_Average,Phase,Error = fold_light_curve(count,time,error,MAXI_PERIOD,num_bins=64)

# plt.errorbar(Phase,Weighted_Average,yerr = Error,fmt = '.',color = 'r')
# plt.title('fold light curve(maxi)')
# plt.legend(loc = 'upper right')
# plt.xlabel('Phase')
# plt.ylabel('Countrate')
# # plt.ylim(0.5,1.6)
# plt.show()  

# -----------------------------------------------------------------------

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
# plt.title('dynamic power spectrum(maxi(2-20keV))')
# plt.show()

# ----------------------------------------------------------------

def multiple_sinusoidal_model(x, y, err, n):
    num_param = 2 * n + 1
    Z = np.zeros((num_param, len(x)))
    Z[0, :] = 1.0
    for k in range(n):
        for i in range(len(x)):
            Z[2 * k + 1, i] = np.cos(2.0 * np.pi * (k + 1) * x[i])
            Z[2 * k + 2, i] = np.sin(2.0 * np.pi * (k + 1) * x[i])

    alpha = np.zeros((num_param, num_param))
    beta = np.zeros(num_param)

    for k in range(num_param):
        beta[k] += np.sum(y * Z[k] / err**2)
        for j in range(num_param):
            alpha[k, j] += np.sum(Z[k] * Z[j] / err**2)

    param = np.linalg.inv(alpha) @ beta

    result = param[0] * Z[0]
    if n > 0:
        for k in range(n):
            result += param[2 * k + 1] * Z[2 * k + 1]
            result += param[2 * k + 2] * Z[2 * k + 2]
    
    return result

n = 10 # order

fitted_curve = multiple_sinusoidal_model(Phase, Weighted_Average, Error, n)

# 计算自由度
chi_square = np.sum((Weighted_Average - fitted_curve)**2 / Error**2)
num_data_points = len(Phase)/2
dof = num_data_points - (2 * n + 1)
chi_square_reduced = chi_square / dof
min_y_index = np.where((Phase >= -0.5) & (Phase <= 0.5))[0][np.argmin(fitted_curve[(Phase >= -0.5) & (Phase <= 0.5)])]
fiducial_point_phase = Phase[min_y_index]
# residual = 

# plt.errorbar(Phase, Weighted_Average, yerr=Error, fmt='.')
# plt.plot(Phase, fitted_curve, linestyle='-', color='orange')
# plt.xlabel('phase')
# plt.ylabel('count rate')
# plt.title('RXTE ASM')
# plt.legend()
# plt.show()

# print(f"Chi-Square: {chi_square:.2f}")
# print(f"DOF: {dof:.2f}")
# print(f"Reduced Chi-Square: {chi_square_reduced:.2f}")
# print(f"fiducial point phase: {fiducial_point_phase:.5f}")
# print(f"f value: {f_value:.5f}")
# print(f"p value: {p_value:.5f}")