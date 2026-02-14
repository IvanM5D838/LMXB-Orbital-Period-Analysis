import matplotlib.pyplot as plt
import numpy as np
import statistics
from astropy.timeseries import LombScargle
import astropy.units as u
import lightkurve
import math
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import curve_fit

XMMNEWTON_PERIOD =  0.713014
LOCAL_PERIOD = 0.714272953030926

def cut_data(time, count, error, start_time, end_time):
    start_idx = next((i for i, t in enumerate(time) if t >= start_time), None)
    end_idx = next((i for i, t in enumerate(time) if t > end_time), None)
    
    if end_idx is None or end_idx >= len(time):
        end_idx = len(time) - 1
    
    return time[start_idx:end_idx], count[start_idx:end_idx], error[start_idx:end_idx]

# ---------------------------------------------------------------------------

def data_red(timeListc,countRateListc,errorListc):
    countRateList = list(countRateListc)
    timeList = list(timeListc)
    errorList = list(errorListc)

    preDelList=[]

    for i in range(0,len(timeList)):
        # method 1
        if countRateList[i] == 0  or errorList[i] == 0:# or abs(countRateList[i] / errorList[i]) < 2:
            preDelList.append(i)
        
    preDelList = sorted(preDelList,reverse=True)
    #pop data in preDelList
    for i in preDelList:
        countRateList.pop(i)
        timeList.pop(i)
        errorList.pop(i)
    return np.array(timeList),np.array(countRateList),np.array(errorList)

# ---------------------------------------------------------------------------

MJDREF = 50814

# with open(f"D:/For Institute/NCU/高能天文實驗室/vscode/X2127+119/xmmnewton/pn(1.5-12keV)/xmmnewton_merged_lc.tsv", "r") as fobj:
#     time, count, error= np.loadtxt(fobj, skiprows=1, usecols=(0, 1, 2), unpack=True)

# with fits.open('D:/For Institute/NCU/高能天文實驗室/vscode/X2127+119/xmmnewton/pn(1.5-12keV)/xmmnewton_merged_lc.fits') as hdul:
#     time = hdul[1].data.field('Time') # 假設資料在第一個擴展HDU中
#     count = hdul[1].data.field('Count')
#     error = hdul[1].data.field('Error')

with fits.open('D:/For Institute/NCU/高能天文實驗室/vscode/X2127+119/xmmnewton/pn(1.5-12keV)/0087350901_BC2.lc') as hdul:
    time = hdul[1].data.field('TIME') # 假設資料在第一個擴展HDU中
    count = hdul[1].data.field('RATE')
    error = hdul[1].data.field('ERROR')

time = MJDREF + time / 86400

TimeList, CountRateList, ErrorList = data_red(time,count,error)

plt.errorbar(TimeList - 53140, CountRateList, yerr=ErrorList, linestyle = 'none', fmt = '.')
plt.xlabel('MJD - 53140')
plt.ylabel('Count Rate(ct/s)')
plt.title('Light Curve(XMM-Newton/pn)')
plt.show()

# ----------------------------------------------------------------

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

# # 使用函数进行数据分组
# grouped_time, grouped_count_rate = group_data(time, count, interval=20)

# # 绘制演变图
# plt.figure(figsize=(10, 6))
# plt.scatter(grouped_time, grouped_count_rate, color='blue', s=15, label='20day bin', alpha=0.6)
# plt.scatter(time, count, color='r', s=6, label='1day bin', alpha=0.3)
# plt.legend(loc='upper right')
# plt.xlabel('Time')
# plt.ylabel('Count Rate')
# plt.title('Count Rate Evolution with Time Binned (20-day intervals)')
# plt.show()

# ----------------------------------------------------------------

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

# Frequency,Power = LS_power(TimeList,CountRateList)  

# plt.plot(Frequency,Power,'b') 
# plt.xlabel('frequency(cycle/year)')
# plt.ylabel('ls_power')
# plt.title('LS_power spectrum(xmmnewton)')
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

# TC,DPS = dynamic_ls_periodogram(time,count,win = 100,mov = 10) # we can see spin period when window is low

# plt.imshow(DPS,aspect = 'auto',extent = [np.min(TC),np.max(TC),np.min(Frequency),np.max(Frequency)],cmap='jet')
# plt.xlabel(r'Time(MJD-50000)',fontsize = 12)
# plt.ylabel(r'Frequency(cycles/year)',fontsize = 12)
# plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(1))
# plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(1))
# plt.colorbar()
# plt.title('dynamic power spectrum(xmmnewton)')
# plt.show()

# -----------------------------------------------------------------------------

def fold_light_curve(countrate, time, error, period, num_bins):
    time = np.array(time)
    countrate = np.array(countrate)
    error = np.array(error)
    
    # 計算相位，將其範圍調整為0到2
    phase = np.mod(time - 47790.463, period) / period

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
    phase_0_to_2 = np.linspace(-1 + bin_width / 2, 1 - bin_width / 2, 2 * num_bins)

    # 刪除folded_countrate_0_to_2 = nan 或 bin_errors_0_to_2 = inf的資料點
    valid_indices = np.where(~np.isnan(folded_countrate_0_to_2) & ~np.isinf(bin_errors_0_to_2))
    folded_countrate_0_to_2 = folded_countrate_0_to_2[valid_indices]
    bin_errors_0_to_2 = bin_errors_0_to_2[valid_indices]
    phase_0_to_2 = phase_0_to_2[valid_indices]

    return folded_countrate_0_to_2, phase_0_to_2, bin_errors_0_to_2

Weighted_Average,Phase,Error = fold_light_curve(CountRateList, TimeList, ErrorList, XMMNEWTON_PERIOD, num_bins=32)

# plt.errorbar(Phase,Weighted_Average,yerr = Error,fmt = '.', color = 'r')
# plt.axvspan(-0.5, -0.15, color='blue', alpha=0.15, label='Partial Eclipse Region')
# plt.axvspan(-0.05, 0.15, color='red', alpha=0.15, label='X-ray Dip Region')
# plt.title('Folded Light Curve(XMM-Newton/pn)')
# plt.legend(loc = 'upper right')
# plt.xlabel('Phase')
# plt.ylabel('Count Rate(ct/s)')
# plt.show()  

# ----------------------------------------------------------------

def multiple_sinusoidal_model(x, y, err, n):
    num_param = 2 * n + 1
    num_data_points = len(x) / 2
    Z = np.zeros((num_param, len(x)))
    Z[0, :] = 1.0
    for k in range(1, n+1):
        Z[2 * k - 1, :] = np.cos(2.0 * np.pi * k * x)
        Z[2 * k, :] = np.sin(2.0 * np.pi * k * x)

    alpha = np.zeros((num_param, num_param))
    beta = np.zeros(num_param)

    for k in range(num_param):
        beta[k] += np.sum(y * Z[k] / err**2)
        for j in range(num_param):
            alpha[k, j] += np.sum(Z[k] * Z[j] / err**2)

    param = np.linalg.inv(alpha) @ beta

    model_x = np.linspace(-1, 1, num=10000)
    model_Z = np.zeros((num_param, len(model_x)))
    model_Z[0, :] = 1.0
    for k in range(1, n+1):
        model_Z[2 * k - 1, :] = np.cos(2.0 * np.pi * k * model_x)
        model_Z[2 * k, :] = np.sin(2.0 * np.pi * k * model_x)

    model_y = param[0] * model_Z[0]
    for k in range(1, num_param):
        model_y += param[k] * model_Z[k]

    obs_Z = np.zeros((num_param, len(x)))
    obs_Z[0, :] = 1.0
    for k in range(1, n+1):
        obs_Z[2 * k - 1, :] = np.cos(2.0 * np.pi * k * x)
        obs_Z[2 * k, :] = np.sin(2.0 * np.pi * k * x)
    
    obs_y = param[0] * obs_Z[0]
    for k in range(1, num_param):
        obs_y += param[k] * obs_Z[k]

    return model_x, model_y, num_data_points, x, obs_y

n = 6 # order

model_phase, model_count, num_data_points, exp_phase, exp_count= multiple_sinusoidal_model(Phase, Weighted_Average, Error, n)

# 计算自由度
chi_square = np.sum((Weighted_Average - exp_count)**2 / Error**2)
dof = num_data_points - (2 * n + 1)
chi_square_reduced = chi_square / dof
min_y_index = np.where((model_phase >= -0.5) & (model_phase <= 0.5))[0][np.argmin(model_count[(model_phase >= -0.5) & (model_phase <= 0.5)])]
fiducial_point_phase = model_phase[min_y_index]

# plt.errorbar(Phase, Weighted_Average, yerr=Error, fmt='.', label = f'k = {n}')
# plt.plot(model_phase, model_count, linestyle='-', color='orange')
# plt.xlabel('phase')
# plt.ylabel('count rate(ct/s)')
# plt.title('multi-sinusoidal fitting(xmmnewton)')
# plt.legend()
# plt.show()

# print(f"Chi-Square: {chi_square:.2f}")
# print(f"DOF: {dof:.2f}")
# print(f"Reduced Chi-Square: {chi_square_reduced:.2f}")
# print(f"fiducial point phase: {fiducial_point_phase:.5f}")
# print(f"f value: {f_value:.5f}")
# print(f"p value: {p_value:.5f}")

# ----------------------------------------------------------------

def monte_carlo_simulation(count, phase, error):
    num_simulation = 10000
    simulated_phase = []

    for _ in range(num_simulation):
        # 生成带有高斯噪声的模拟计数
        count_sim = count + np.random.normal(loc=0, scale=error, size=len(count))
        
        # 筛选出 phase 在 [-0.5, 0.5] 范围内的索引
        valid_indices = np.where((phase >= -0.5) & (phase <= 0.5))[0]
        
        # 在有效索引中找到最小 count_sim 的索引
        min_y_index = valid_indices[np.argmin(count_sim[valid_indices])]
        
        # 添加对应的 phase 值到模拟相位列表
        simulated_phase.append(phase[min_y_index])

    # 计算模拟相位的标准差
    simulated_phase = np.array(simulated_phase)
    standard_deviation = np.sqrt(np.sum((simulated_phase - np.mean(simulated_phase))**2) / (num_simulation - 1))

    return simulated_phase, standard_deviation

def gaussian(phase, amp, mean, SDV): # 定义高斯函数
    return amp * np.exp(-(phase - mean)**2 / (2 * SDV**2))

simulated_phase_sin, standard_deviation = monte_carlo_simulation(exp_count, exp_phase, Error)
time_mid_point = np.sum(TimeList) / len(TimeList)
hist, bin_edges = np.histogram(simulated_phase_sin, bins=8) # 使用 bin 的中心作为 x 轴数据
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
max_index = np.argmax(hist) # 找到最大计数值对应的索引

x_values = np.linspace(-0.5, 0.5, 1000)
initial_guess = [np.max(hist), bin_centers[max_index], standard_deviation] # 初始猜测的振幅、均值和标准差
params, covariance = curve_fit(gaussian, bin_centers, hist, p0=initial_guess)
fitted_amp, fitted_mean, fitted_stddev = params

# plt.hist(simulated_phase_sin, bins=16, edgecolor='black')
# plt.plot(x_values, gaussian(x_values, *params), 'r-', lw=2, label='Fitted Gaussian Distribution')
# plt.xlabel('Phase')
# plt.ylabel('Counts/bin')
# plt.xlim(-0.5, 0.5)
# plt.title('Monte Carlo Simulation (XMM-Newton)')
# plt.legend()
# plt.show()

# print(f"minimum phase: {fitted_mean} (Monte Carlo simulation)")
# print(f"standard deviation: {fitted_stddev}")
# print(f"time mid point : {time_mid_point}")

# ----------------------------------------------------------------

# # 定義MJDREF
# MJDREF = 50814

# # 創建一個空的 DataFrame 來儲存合併後的資料
# merged_df = pd.DataFrame(columns=['Time', 'Count', 'Error'])

# # 處理每個檔案
# filenames = [
#     '0087350101_BC.lc',
#     '0087350601_BC.lc',
#     '0087350801_BC.lc',
#     '0087350901_BC.lc'
# ]

# for filename in filenames:
#     # 讀取 FITS 檔案
#     with fits.open(f'D:/For Institute/NCU/高能天文實驗室/vscode/X2127+119/xmmnewton/pn_timing_mode(1.5-12keV)/{filename}') as hdul:
#         data = hdul[1].data  # 假設資料在第一個 HDU 中
        
#         # 計算新的時間
#         time = data['TIME']
        
#         # 將處理後的資料加入到合併後的 DataFrame 中
#         df = pd.DataFrame({'Time': time, 'Count': data['RATE'], 'Error': data['ERROR']})
#         merged_df = pd.concat([merged_df, df])

# # 按時間列排序
# merged_df = merged_df.sort_values(by='Time')

# # 刪除Count和Error值為0的資料點
# merged_df = merged_df[(merged_df['Count'] != 0) & (merged_df['Error'] != 0)]

# # 將Time、Count和Error列舍入到小数点后一位
# merged_df['Time'] = merged_df['Time'].round(5)
# merged_df['Count'] = merged_df['Count'].round(5)
# merged_df['Error'] = merged_df['Error'].round(5)

# # 將 DataFrame 轉換為 Astropy Table
# table = Table.from_pandas(merged_df)

# # 將 Table 保存為 FITS 文件
# table.write('xmmnewton_merged_lc.fits', format='fits', overwrite=True)

# # 儲存合併後的資料
# merged_df.to_csv('xmmnewton_merged_lc.tsv', index=False, sep='\t')


