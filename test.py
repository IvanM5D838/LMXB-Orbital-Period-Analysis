
import os
import os.path
from os import path
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from matplotlib.colors import LogNorm
from scipy.stats import poisson
from IPython import display
import numpy as np
from scipy.stats import norm
from astropy.time import Time
from scipy.signal import convolve, correlate
from matplotlib.animation import FuncAnimation
from scipy.stats import t

#-------------------------------------------------------------------------

# 生成隨機數據（示例）
np.random.seed(0)
data = np.random.rand(100) * 10  # 生成0到10之間的隨機數據

# 定義移動框的大小（窗口大小）
window_size = 6

# 初始化平滑後的數據列表
smoothed_data = []

# 遍歷數據並應用移動框平滑
for i in range(len(data)):
    # 計算移動框的範圍
    start = max(0, i - window_size // 2)
    end = min(len(data), i + window_size // 2 + 1)
    
    # 取得移動框內的數據並計算平均值
    window_data = data[start:end]
    smoothed_value = np.mean(window_data)
    
    # 將平滑值添加到平滑後的數據列表中
    smoothed_data.append(smoothed_value)

# 繪製原始數據和平滑後的數據
# plt.figure(figsize=(10, 6))
# plt.plot(data, label='original data', marker='o')
# plt.plot(smoothed_data, label=f'moving box (window size={window_size})', linestyle='-', linewidth=2)
# plt.xlabel('time step')
# plt.ylabel('data')
# plt.legend()
# plt.title('moving box smoothing')
# plt.grid(True)
# plt.show()

#-------------------------------------------------------------------------

# 生成隨機數據（示例）
np.random.seed(0)
data = np.random.randn(100) * 2 + 5  # 生成均值為5、標準差為2的隨機數據

# 定義核函數（這裡使用標準正態分佈作為核函數）
kernel = norm(loc=0, scale=1) # loc = mean , scale = standdard deviation

# 初始化平滑後的數據列表
smoothed_data = []

# 遍歷數據並應用核平滑
for x in data:
    # 計算核函數在每個數據點的值
    kernel_values = kernel.pdf((data - x) / 1)  # 這裡的1是帶寬參數，控制平滑程度 pdf = probability density function
    
    # 使用核函數的加權平均來計算平滑值
    smoothed_value = np.sum(data * kernel_values) / np.sum(kernel_values)
    
    # 將平滑值添加到平滑後的數據列表中
    smoothed_data.append(smoothed_value)

# 繪製原始數據和平滑後的數據
# plt.figure(figsize=(10, 6))
# plt.plot(data, label='original data', marker='o')
# plt.plot(smoothed_data, label='kernel smooth', linestyle='-', linewidth=2)
# plt.xlabel('data point')
# plt.ylabel('data')
# plt.legend()
# plt.title('kernel smooth')
# plt.grid(True)
# plt.show()

#--------------------------------------------------------------------------

import random

def monte_carlo_pi(num_samples):
    inside_circle = 0
    
    # x,y範圍為(0,1)時在第一象限(1/4圓)
    for _ in range(num_samples):
        x = random.uniform(0, 1)
        y = random.uniform(0, 1)
        distance = x**2 + y**2
        
        if distance <= 1:
            inside_circle += 1
    
    pi_estimate = (inside_circle / num_samples) * 4
    return pi_estimate

# if __name__ == "__main__":
#     num_samples = 10000000  # 試驗次數，可以根據需要調整
#     estimated_pi = monte_carlo_pi(num_samples)
#     print(f"估計的π值為：{estimated_pi}")

#-------------------------------------------------------------------------

# class Student:
#     def __init__(self, name, score):
#         self.name = name
#         self.score = score

# # 自定义比较函数，按分数降序排列
# def compare_students(student1, student2):
#     return student2.score - student1.score

# students = [Student("Alice", 85), Student("Bob", 92), Student("Charlie", 78)]

# # 使用自定义比较函数进行排序
# sorted_students = sorted(students, key=lambda student: student.score, reverse=True)

# for student in sorted_students:
#     print(student.name, student.score)

#--------------------------------------------------------------------------

# # 生成随机数的数量
# num_samples = 10000

# # 生成均匀分布的随机数 U1 和 U2
# u1 = np.random.rand(num_samples)
# u2 = np.random.rand(num_samples)

# # 使用Box-Muller方法生成标准正态分布的随机数 V1 和 V2
# v1 = np.sqrt(-2 * np.log(u1)) * np.cos(2 * np.pi * u2)
# v2 = np.sqrt(-2 * np.log(u1)) * np.sin(2 * np.pi * u2)

# # 绘制生成的标准正态分布随机数的直方图
# plt.figure(figsize=(10, 5))

# plt.subplot(1, 2, 1)
# plt.hist(v1, bins=50, density=True, color='blue', alpha=0.7)
# plt.title('Histogram of V1')

# plt.subplot(1, 2, 2)
# plt.hist(v2, bins=50, density=True, color='red', alpha=0.7)
# plt.title('Histogram of V2')

# plt.tight_layout()
# plt.show()

#----------------------------------------------------------------------------

# # 设置泊松分布的参数
# lambda_parameter = 50.0  # 平均事件发生率

# # 生成泊松分布的随机数据
# data = np.random.poisson(lambda_parameter, size=10000)

# # 绘制直方图
# plt.hist(data, bins=50, density=True, color='blue', alpha=0.7)
# plt.title('Poisson Distribution Example')
# plt.xlabel('Number of Events')
# plt.ylabel('Probability Density')
# plt.show()

#----------------------------------------------------------------------------

# from scipy.stats import gamma

# # 设置伽玛分布的参数
# shape_parameter = 3.0  # 形状参数(k)
# scale_parameter = 1.0  # 尺度参数(θ)

# # 生成随机变量
# data = gamma.rvs(a=shape_parameter, scale=scale_parameter, size=1000)

# # 绘制直方图
# plt.hist(data, bins=30, density=True, alpha=0.7, color='blue', label='Histogram')

# # 绘制概率密度函数
# x = np.linspace(0, 10, 100)
# pdf = gamma.pdf(x, a=shape_parameter, scale=scale_parameter)
# plt.plot(x, pdf, 'r-', lw=2, label='PDF')

# # 添加标签和图例
# plt.title('Gamma Distribution Example')
# plt.xlabel('Value')
# plt.ylabel('Probability Density')
# plt.legend()

# # 显示图形
# plt.show()

#--------------------------------------------------------------------------

# # 采样率（样本频率）
# sampling_rate = 1000  # Hz
# # 信号的持续时间
# duration = 1.0  # 秒
# # 信号的时间数组
# t = np.linspace(0, duration, int(sampling_rate * duration), endpoint=False)

# # 生成信号
# signal = np.sin(2 * np.pi * 5 * t) + np.sin(2 * np.pi * 10 * t)

# # 进行DFT
# dft = np.fft.fft(signal)

# # 绘制信号
# plt.subplot(2, 1, 1)
# plt.plot(t, signal)
# plt.title('Signal')

# # 绘制频谱
# plt.subplot(2, 1, 2)
# frequencies = np.fft.fftfreq(len(t), 1.0 / sampling_rate)
# plt.plot(frequencies, np.abs(dft))
# plt.title('Spectrum')

# plt.tight_layout()
# plt.show()

#--------------------------------------------------------------------------

# student distribution 
# 也称为t分布，是一种概率分布，常用于统计学和假设检验中。
# 它的形状类似于正态分布（钟形曲线），但具有更宽的尾部。
# 以下是一个关于学生t分布的示例情景：

# 假设你是一家制药公司的质量控制工程师，负责检验一种新药的生产批次。
# 你要确定这批药品的平均药物含量是否符合规定，规定的平均含量为100毫克。

# 你从这批药品中随机抽取了10个样本，并测量了每个样本的药物含量。
# 现在，你想使用学生t分布来进行统计推断，以确定这批药品的平均药物含量是否与规定相符。

# 以下是具体的步骤：

# 1. **采样和数据收集：** 从批次中随机抽取10个样本，并测量每个样本的药物含量。
# 得到以下数据（单位：毫克）：

# 98.5, 99.0, 100.5, 101.0, 98.8, 100.2, 99.5, 100.0, 98.6, 99.8

# 2. **计算样本均值和样本标准差：** 计算这10个样本的均值和标准差。

#    样本均值 \( \bar{x} = 99.69 \) 毫克
#    样本标准差 \( s = 1.06 \) 毫克

# 3. **假设检验：** 假设你想进行双侧假设检验，检验是否存在统计显著的差异。
# 你使用学生t分布来计算t统计量。

#    假设零假设 \( H_0 \)：批次的平均药物含量等于100毫克。
#    对立假设 \( H_1 \)：批次的平均药物含量不等于100毫克。

# 4. **计算t统计量：** 使用以下公式计算t统计量：

#    \[ t = \frac{\bar{x} - \mu}{\frac{s}{\sqrt{n}}} \]

#    其中，\( \bar{x} \) 是样本均值，\( \mu \) 是假设的总体均值（100毫克），\( s \) 是样本标准差，\( n \) 是样本大小（10）。

#    计算得到 \( t \approx -0.628 \)

# 5. **查找临界值：** 使用学生t分布表或计算机软件，查找自由度为9（样本大小减1）的t分布在双侧显著性水平为0.05时的临界值。假设使用95%的置信水平。

#    对于双侧检验，临界值大约为 \(\pm 2.262\)。

# 6. **做出决策：** 比较计算得到的t统计量和临界值。
#    - 如果 \( |t| > 2.262 \)，则拒绝零假设，说明平均药物含量与规定不相符。
#    - 如果 \( |t| \leq 2.262 \)，则无法拒绝零假设，说明平均药物含量与规定相符。

# # 在这个示例中，计算得到的t统计量的绝对值小于2.262，因此无法拒绝零假设，表明这批药品的平均药物含量与规定相符。这是一个使用学生t分布进行假设检验的示例。
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy import stats

# # 样本均值、样本标准差和样本大小
# sample_mean = 99.69
# sample_std = 1.06
# sample_size = 10

# # 假设的总体均值
# population_mean = 100

# # 计算 t 统计量
# t_statistic = (sample_mean - population_mean) / (sample_std / np.sqrt(sample_size))

# # 创建 t 分布的概率密度函数（PDF）图表
# x = np.linspace(-5, 5, 400)
# pdf = stats.t.pdf(x, df = sample_size - 1)

# plt.plot(x, pdf, label='t distribution PDF', color='blue')
# plt.axvline(x=t_statistic, color='red', linestyle='--', label='t_statistic')

# # 标记临界值
# critical_value = 2.262  # 95% 置信水平的双侧临界值
# plt.axvline(x=critical_value, color='green', linestyle='--', label='critical value')
# plt.axvline(x=-critical_value, color='green', linestyle='--')

# # 添加标签和图例
# plt.xlabel('t_statistic ')
# plt.ylabel('probability density')
# plt.title('PDF and t_statistic')
# plt.legend()

# plt.show()

#---------------------------------------------------------

# sigma_f有兩種方法可以找

# 1.sigma_f = 1 / (2 * T) (T = time duration of the light curve)  
# 
# 2. Half Width Half Maximum : 
# (p_i - (1/2) * p_max) * (p_i + (1/2) * p_max) <=0 ( 判斷頻率的power是否在一半的上下 )

# ------------------------------------------------------------------------

# # 从文本数据创建一个NumPy数组
# column1_data = np.array(timeList, dtype=float)
# column2_data = np.array(countRateList[2], dtype=float)
# column3_data = np.array(errorList[2], dtype=float)
# column4_data = np.array(countRateList[3], dtype=float)
# column5_data = np.array(errorList[3], dtype=float)

# col1 = fits.Column(name='TIME', format='E', array=column1_data)
# # col2 = fits.Column(name='RATE', format='E', array=column2_data)
# # col3 = fits.Column(name='ERROR', format='E', array=column3_data)
# col4 = fits.Column(name='RATE', format='E', array=column4_data)
# col5 = fits.Column(name='ERROR', format='E', array=column5_data)

# columns = fits.ColDefs([col1, col4, col5])
# hdu = fits.BinTableHDU.from_columns(columns)
# hdulist = fits.HDUList([fits.PrimaryHDU(), hdu])

# # 设置TSTART和TSTOP的时间信息
# start_time = Time('2009-08-11T09:45:00')
# stop_time = Time('2023-08-11T11:14:57')
# hdulist[0].header['TSTART'] = start_time.isot
# hdulist[0].header['TSTOP'] = stop_time.isot

# hdulist.writeto('OrbitMAXIlc_Band3.fits', overwrite=True)

# -----------------------------------------------------------------------

# hdu1 = fits.open('asm_x0114650_lc.fits',mode='update')

# hdu2 = fits.open('4U0114+650/maxi/MAXI_Band2.fits',mode='update')

# 获取数据表（recarray）
# Field = hdu1[1].data.field
# count = Field('RATE')
# time = Field('TIME')
# error = Field('ERROR')

# def Std_fix(timeListc,countRateListc,errorListc,mode='default',times=int):
#     countRateList = list(countRateListc)
#     timeList = list(timeListc)
#     errorList = list(errorListc)
#     countRateListStd = np.std(countRateList,ddof=1)
#     errorListStd = np.std(errorList,ddof=1)
#     countRateListAvg = sum(countRateList) / len(countRateList)
#     errorListAvg = sum(errorList) / len(errorList)

#     preDelList=[]

#     for i in range(0,len(timeList)):
#         if mode == 'all' and abs(countRateList[i] - countRateListAvg) > times * countRateListStd or (errorList[i] - errorListAvg) > times * errorListStd or np.log(errorList[i]) < 10^(-3) :
#             preDelList.append(i)

#     #sort&reverse list
#     preDelList = sorted(preDelList,reverse=True)
#     #pop data in preDelList
#     for i in preDelList:
#         countRateList.pop(i)
#         timeList.pop(i)
#         errorList.pop(i)
#     return timeList,countRateList,errorList

# fix_time,fix_count,fix_error = Std_fix(time,count,error,'all',2)

# def LS_power(time=list, count=list):
#     frequency = np.linspace(0.01, 0.1, 1000)
#     power = np.zeros_like(frequency)
#     N = len(frequency)
#     mean_freq = []

#     for i, freq in enumerate(frequency):
#         omega = 2 * np.pi * freq
#         tau = np.arctan(np.sum(np.sin(2 * omega * np.array(time)))/np.sum(np.cos(2 * omega * np.array(time))))/ (2 * omega) 
#         cos_term = np.cos(omega * np.array(time - tau))
#         sin_term = np.sin(omega * np.array(time - tau))
#         num_sin = np.sum((np.array(count) - np.mean(np.array(count))) * sin_term)
#         num_cos = np.sum((np.array(count) - np.mean(np.array(count))) * cos_term)
#         den_cos = np.sum(cos_term**2)
#         den_sin = np.sum(sin_term**2)
#         power[i] = (num_cos**2 / den_cos + num_sin**2 / den_sin) / 2
#     normal_power = power/np.mean(power)
#     # for idx,i in enumerate(power):
#     #     if 0.025 < frequency[idx] < 0.03  and i > 1 * 1e-6:
#     #         mean_freq.append(frequency[idx])
#     # print(np.mean(mean_freq))

#     return frequency, normal_power

# Frequency,Power = LS_power(fix_time,fix_count)

# plt.plot(Frequency, Power, color = 'red')
# plt.show()

# Header = hdu1[1].header

# for i,idx in enumerate(Field):
#     if i == 0:
#         print(idx)
#     else:
#         pass

# print(Field)
# print(Header)

# 保存更改
# hdu1.flush()
# hdu1.close()

# print(hdu1[1].data.field('MJD'))
# print(hdu2[1].data.field('MJD'))

# -----------------------------------------------------------------------------

# rms_amplitude 

# mean,rms,percentage

# ---------------------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt

# # 生成示例數據（時間和測量值）
# t = np.linspace(0, 10, 1000)  # 時間
# y = 3.0 * np.sin(2 * np.pi * 1.0 * t) + 2.0 * np.sin(2 * np.pi * 2.5 * t) + np.random.normal(0, 0.5, t.shape)

# # 頻率範圍
# frequencies = np.linspace(0.1, 5.0, 1000)

# # 初始化 Lomb-Scargle 功率數組
# power = np.zeros_like(frequencies)

# # 計算 Lomb-Scargle 功率
# for i, freq in enumerate(frequencies):
#     phi = 2 * np.pi * freq * t
#     cos_term = np.cos(phi)
#     sin_term = np.sin(phi)
    
#     num = np.sum((y - np.mean(y)) * cos_term)
#     den_cos = np.sum(cos_term**2)
#     den_sin = np.sum(sin_term**2)
    
#     power[i] = 0.5 * ((num / den_cos)**2 + (num / den_sin)**2)

# # 繪製結果
# plt.plot(frequencies, power)
# plt.xlabel('Frequency')
# plt.ylabel('Lomb-Scargle Power')
# plt.title('Lomb-Scargle Periodogram')
# plt.show()

# -----------------------------------------------------------------

# import matplotlib.pyplot as plt
# import matplotlib.patches as patches

# # 创建一个新的图
# fig, ax = plt.subplots()

# # 创建一个圆形对象
# circle = patches.Circle((0.5, 0.5), 0.4, fill=False, color='blue')

# # 将圆形添加到图中
# ax.add_patch(circle)

# # 设置图的范围
# ax.set_xlim(0, 1)
# ax.set_ylim(0, 1)

# # 显示图
# plt.gca().set_aspect('equal', adjustable='box')
# plt.show()

# --------------------------------------------------------------------

# import matplotlib.pyplot as plt
# import numpy as np

# # 常数
# h = 6.626e-34  # 普朗克常数 (Joule秒)
# c = 3.0e8      # 光速 (米/秒)
# k = 1.38e-23   # 玻尔兹曼常数 (Joule/开尔文)
# T = 10000      # 温度（开尔文）

# # 波长范围（从0.1μm到100μm，可以根据需要调整范围）
# wavelengths = np.logspace(-7, -3, 10**6) # 波长范围

# # 计算黑体辐射能量密度
# def blackbody_radiation(wavelength, temperature):
#     numerator = 2 * h * c**2
#     denominator = wavelength**4 * (np.exp((h * c) / (wavelength * k * temperature)) - 1)
#     return numerator / denominator

# # 计算能量密度
# energy_density = blackbody_radiation(wavelengths, T)

# # 绘制SED曲线
# plt.figure(figsize=(10, 6))
# plt.plot(wavelengths * 1e6, energy_density, color='blue')  # 将波长转换为μm
# plt.xscale('log')
# plt.xlabel('wavelength(μm)')
# plt.xlim(0.1,1000)
# plt.yscale('log')
# plt.ylabel('energy density(W/m$^{2}$)')
# plt.title(f'{int(T)}K energy distribution')
# plt.show()

# ------------------------------------------------------------------

# def lambda_function(wavelength,Rv,Av):
#     x = 1 / wavelength
#     if 0.3 <= x <= 1.1:
#         a = 0.574 * x**1.61
#         b = -0.527 * x**1.61
#     elif 1.1 <= x <= 3.3:
#         y = (x - 1.82)
#         a = 1 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 + 0.72085 * y**4 + 0.01979 * y**5 - 0.77530 * y**6 + 0.32999 * y**7
#         b = 1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 5.38434 * y**4 - 0.62251 * y**5 + 5.3026 * y**6 - 2.09002 * y**7
#     else:
#         pass 
#     A_lambda = (a + b / Rv)*Av
#     return round(x,4),round(a,4),round(b,4),round(A_lambda,4)

# wavelength = [0.4357 , 0.5386 , 0.6430 , 0.9140 , 1.2350 , 1.6220 , 2.1590]
# A_lambda = []
# for idx,i in enumerate(wavelength):
#     A = lambda_function(i,3.1,0.421)
#     A_lambda.append(A)
#     print(A) 

# ----------------------------------------------------------------

# # 生成两个信号
# x = np.linspace(0, 10, 100)
# signal1 = np.sin(x)
# signal2 = np.exp(-x)

# # 计算信号的卷积
# convolution_result = convolve(signal1, signal2, mode='full')

# # 创建图形和轴
# fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 12))

# # 初始化动态显示的线
# line1, = ax1.plot([], [], lw=2, label='Signal 1')
# line2, = ax2.plot([], [], lw=2, label='Signal 2')
# line3, = ax3.plot([], [], lw=2, label='Convolution Result')
# ax1.set_title('Signal 1 Animation')
# ax2.set_title('Signal 2 Animation')
# ax3.set_title('Convolution Animation')

# # 设置图形的初始状态
# def init():
#     ax1.set_xlim(0, len(signal1))
#     ax1.set_ylim(min(signal1), max(signal1))
#     ax1.legend(loc='upper right')
#     ax2.set_xlim(0, len(signal2))
#     ax2.set_ylim(min(signal2), max(signal2))
#     ax2.legend(loc='upper right')
#     ax3.set_xlim(0, len(convolution_result))
#     ax3.set_ylim(min(convolution_result), max(convolution_result))
#     ax3.legend(loc='upper right')
#     return line1, line2, line3

# # 动态显示信号1的变化
# def animate_signal1(i):
#     y = np.zeros_like(signal1)
#     y[:i+1] = signal1[:i+1]
#     line1.set_data(range(len(y)), y)
#     return line1,

# # 动态显示信号2的变化
# def animate_signal2(i):
#     y = np.zeros_like(signal2)
#     y[:i+1] = signal2[:i+1]
#     line2.set_data(range(len(y)), y)
#     return line2,

# # 动态显示卷积结果的变化
# def animate_convolution(i):
#     if i >= len(signal1) + len(signal2) - 1:
#         return
#     y = np.zeros_like(convolution_result)
#     y[:i+1] = convolution_result[:i+1]
#     line3.set_data(range(len(y)), y)
#     return line3,

# # 设置初始状态
# def init():
#     max_len = max(len(signal1), len(signal2), len(convolution_result))
#     ax1.set_xlim(0, max_len)
#     ax1.set_ylim(min(signal1), max(signal1))
#     ax1.legend(loc='upper right')
#     ax2.set_xlim(0, max_len)
#     ax2.set_ylim(min(signal2), max(signal2))
#     ax2.legend(loc='upper right')
#     ax3.set_xlim(0, max_len)
#     ax3.set_ylim(min(convolution_result), max(convolution_result))
#     ax3.legend(loc='upper right')
#     return line1, line2, line3

# # 设置动画
# ani1 = FuncAnimation(fig, animate_signal1, frames=len(signal1), interval=100, init_func=init)
# ani2 = FuncAnimation(fig, animate_signal2, frames=len(signal2), interval=100, init_func=init)
# ani3 = FuncAnimation(fig, animate_convolution, frames=len(signal1) + len(signal2), interval=100, init_func=init)

# plt.show()

# ----------------------------------------------------------------

def cubic_spline_interpolation(x, y):
    n = len(x) - 1
    h = [x[i+1] - x[i] for i in range(n)]
    alpha = [3 * (y[i+1] - y[i]) / h[i] - 3 * (y[i] - y[i-1]) / h[i-1] for i in range(1, n)]

    l = [1] + [0] * (n - 1) + [1]
    mu = [0] * n
    z = [0] * (n + 1)

    for i in range(1, n):
        l[i] = 2 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i-1] - h[i-1] * z[i-1]) / l[i]

    b, c, d = [0] * n, [0] * (n + 1), [0] * n
    for j in range(n-1, -1, -1):
        c[j] = z[j] - mu[j] * c[j+1]
        b[j] = (y[j+1] - y[j]) / h[j] - h[j] * (c[j+1] + 2 * c[j]) / 3
        d[j] = (c[j+1] - c[j]) / (3 * h[j])

    a = y[:n]
    return a, b, c, d

def evaluate_spline(a, b, c, d, x, x_eval):
    n = len(a)
    y_eval = []
    for x_val in x_eval:
        for i in range(n):
            if x_val >= x[i] and x_val <= x[i+1]:
                dx = x_val - x[i]
                y_val = a[i] + b[i] * dx + c[i] * dx**2 + d[i] * dx**3
                y_eval.append(y_val)
                break
    return y_eval

# 示例数据点
x = [0, 1, 2, 3, 4, 5]
y = [0, 0.8, 0.9, 0.1, -0.8, -1]

# 计算三次样条插值的系数
a, b, c, d = cubic_spline_interpolation(x, y)

# # 生成插值点
# x_eval = np.linspace(min(x), max(x), 1000)
# y_eval = evaluate_spline(a, b, c, d, x, x_eval)

# # 绘制原始数据点和插值曲线
# plt.plot(x, y, 'o', label='Data points')
# plt.plot(x_eval, y_eval, '-', label='Cubic spline interpolation')
# plt.legend()
# plt.show()

# ----------------------------------------------------------------

def cubic_spline_coefficients(x, y):
    n = len(x)
    h = np.diff(x)
    alpha = np.zeros(n)
    for i in range(1, n - 1):
        alpha[i] = 3 * (y[i + 1] - y[i]) / h[i] - 3 * (y[i] - y[i - 1]) / h[i - 1]
    
    l = np.zeros(n)
    mu = np.zeros(n)
    z = np.zeros(n)
    l[0] = 1
    mu[0] = 0
    z[0] = 0
    
    for i in range(1, n - 1):
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]
    
    l[n - 1] = 1
    z[n - 1] = 0
    c = np.zeros(n)
    b = np.zeros(n)
    d = np.zeros(n)
    
    for j in range(n - 2, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]
        b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * h[j])
    
    return y[:-1], b, c[:-1], d

def evaluate_cubic_spline(x, a, b, c, d, x_eval):
    n = len(x) - 1
    y_eval = []
    
    for x_val in x_eval:
        for i in range(n):
            if x_val >= x[i] and x_val <= x[i + 1]:
                dx = x_val - x[i]
                y_val = a[i] + b[i] * dx + c[i] * dx**2 + d[i] * dx**3
                y_eval.append(y_val)
                break
    
    return np.array(y_eval)

def monte_carlo_simulation(x, y, error_std, num_simulations, x_fit):
    simulated_curves = []
    for _ in range(num_simulations):
        y_simulated = y + np.random.normal(0, error_std, len(y))
        a, b, c, d = cubic_spline_coefficients(x, y_simulated)
        y_fit_simulated = evaluate_cubic_spline(x, a, b, c, d, x_fit)
        simulated_curves.append(y_fit_simulated)
    return simulated_curves

# # 示例数据
# x = np.linspace(0, 2*np.pi, 10)
# y = np.sin(x)

# # 计算三次样条插值的系数
# a, b, c, d = cubic_spline_coefficients(x, y)

# # 设定插值点
# x_fit = np.linspace(0, 2*np.pi, 200)

# # 计算插值结果
# y_fit = evaluate_cubic_spline(x, a, b, c, d, x_fit)

# # Monte Carlo 模拟参数
# error_std = 0.1  # 噪声标准差
# num_simulations = 100  # 模拟次数

# # 进行 Monte Carlo 模拟
# simulated_curves = monte_carlo_simulation(x, y, error_std, num_simulations, x_fit)

# # 绘图
# plt.figure(figsize=(10, 6))
# plt.plot(x_fit, y_fit, label='Cubic Spline Interpolation', color='blue')
# plt.scatter(x, y, color='red', label='Data Points')

# # 绘制 Monte Carlo 模拟曲线
# for curve in simulated_curves:
#     plt.plot(x_fit, curve, color='gray', alpha=0.2)

# plt.title('Cubic Spline Interpolation with Monte Carlo Simulation')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.legend()
# plt.grid(True)
# plt.show()

# ----------------------------------------------------------------

def z_score_filter(time, count_rate, error, threshold=3):
    # 計算計數率數據的均值和標準差
    mean = np.mean(count_rate)
    std = np.std(count_rate)
    
    # 計算 Z-score
    z_scores = (count_rate - mean) / std
    
    # 篩選數據
    mask = np.abs(z_scores) < threshold
    filtered_time = time[mask]
    filtered_count_rate = count_rate[mask]
    filtered_error = error[mask]
    
    return filtered_time, filtered_count_rate, filtered_error

# ----------------------------------------------------------------

def iqr_filter(time, count_rate, error, factor=1.5):
    # 計算第一四分位數（Q1）和第三四分位數（Q3）
    Q1 = np.percentile(count_rate, 25)
    Q3 = np.percentile(count_rate, 75)
    IQR = Q3 - Q1

    # 定義上下限
    lower_bound = Q1 - factor * IQR
    upper_bound = Q3 + factor * IQR
    
    # 篩選數據
    mask = (count_rate >= lower_bound) & (count_rate <= upper_bound)
    filtered_time = time[mask]
    filtered_count_rate = count_rate[mask]
    filtered_error = error[mask]
    
    return filtered_time, filtered_count_rate, filtered_error

# ----------------------------------------------------------------

def moving_average_filter(time, count_rate, error, window_size):
    # 確保窗口大小是奇數，以便對稱
    if window_size % 2 == 0:
        window_size += 1
    
    # 計算移動平均值
    half_window = window_size // 2
    smoothed_count_rate = np.convolve(count_rate, np.ones(window_size)/window_size, mode='valid')
    
    # 對應更新時間和誤差數據
    smoothed_time = time[half_window: -half_window]
    smoothed_error = error[half_window: -half_window]
    
    return smoothed_time, smoothed_count_rate, smoothed_error

# ----------------------------------------------------------------

# # 样本数据
# data = np.array([2.1, 2.5, 3.3, 2.8, 3.0, 2.6, 2.9])
# # 已知的总体均值
# population_mean = 2.5

# # 样本均值
# sample_mean = np.mean(data)
# # 样本标准差
# sample_std = np.std(data, ddof=0)
# # 样本大小
# n = len(data)
# # 自由度
# df = n - 1

# # 计算 t 值
# t_value = (sample_mean - population_mean) / (sample_std / np.sqrt(n))

# # 计算 p 值（双尾检验）
# p_value = 2 * (1 - t.cdf(np.abs(t_value), df))

# print(f"t-value: {t_value:.3f}")
# print(f"p-value: {p_value:.3f}")

# # 绘制 t 分布并标出 t 值位置
# x = np.linspace(-4, 4, 1000)
# pdf = t.pdf(x, df)

# plt.figure(figsize=(10, 6))
# plt.plot(x, pdf, label=f't-distribution (df={df})')
# plt.axvline(t_value, color='r', linestyle='--', label=f't-value = {t_value:.3f}')
# plt.axvline(-t_value, color='r', linestyle='--')
# plt.fill_between(x, 0, pdf, where=(x >= t_value) | (x <= -t_value), color='red', alpha=0.3)
# plt.title('t-distribution with t-value')
# plt.xlabel('t')
# plt.ylabel('Density')
# plt.legend()
# plt.show()

# ----------------------------------------------------------------

# 设置自由度
df = 10

# 定义t值的范围
x = np.linspace(-4, 4, 1000)

# 计算PDF和CDF
pdf = t.pdf(x, df)
cdf = t.cdf(x, df)

# # 绘制PDF和CDF
# plt.figure(figsize=(12, 6))

# plt.subplot(1, 2, 1)
# plt.plot(x, pdf, label=f't-distribution (df = {df})')
# plt.title('Probability Density Function')
# plt.xlabel('t')
# plt.ylabel('Density')
# plt.legend()

# plt.subplot(1, 2, 2)
# plt.plot(x, cdf, label=f't-distribution (df = {df})')
# plt.title('Cumulative Distribution Function')
# plt.xlabel('t')
# plt.ylabel('Cumulative Probability')
# plt.legend()

# plt.tight_layout()
# plt.show()

# import numpy as np

# # 生成10的5次方個隨機點
# num_points = 10**5
# points = np.random.rand(num_points, 2)

# # 計算點到原點的距離
# distances = np.sqrt(points[:,0]**2 + points[:,1]**2)

# # 計算落在單位圓內的點的數量
# points_inside_circle = np.sum(distances <= 1)

# # 蒙地卡羅方法計算圓周率
# pi_estimate = 4 * points_inside_circle / num_points

# print(pi_estimate)

# import numpy as np
# import matplotlib.pyplot as plt

# # 參數設定
# e = 0.5       # 偏心率
# a = 1               # 半長軸 (可設定為 1 單位)
# theta = np.linspace(0, 2*np.pi, 1000)  # 角度範圍

# # 極坐標方程：r = a(1-e^2)/(1+e*cos(theta))
# r = a * (1 - e**2) / (1 + e * np.cos(theta))

# # 轉換成直角坐標
# x = r * np.cos(theta)
# y = r * np.sin(theta)

# # 繪圖
# plt.figure(figsize=(6,6))
# plt.plot(x, y, label="Orbit (e = 0.5)", color='blue')
# plt.scatter(0, 0, color='red', label="Focus (Central Star)")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.title("Elliptical Orbit with Eccentricity 0.5")
# plt.axis('equal')
# plt.grid(True)
# plt.legend()
# plt.show()

# ----------------------------------------------------------------

# import pyneb as pn
# import numpy as np
# import matplotlib.pyplot as plt

# # 建立 S II 原子對象
# S2 = pn.Atom('S', 2)

# # 固定電子溫度 Te
# Te = 10000  # 單位 K

# # 設定電子密度範圍（log scale）
# ne_vals = np.logspace(0, 5, 500)  # 從10^1 到 10^5 cm^-3

# # 計算每個 ne 對應的線比 I(6716) / I(6731)
# ratio_vals = [S2.getEmissivity(tem=Te, den=ne, wave=6716) / 
#               S2.getEmissivity(tem=Te, den=ne, wave=6731) for ne in ne_vals]

# # 計算對應的電子密度
# ne = S2.getTemDen(int_ratio=0.63, tem=Te, wave1=6716, wave2=6731)

# print(f"Estimated electron density: {ne:.1f} cm^-3")

# # 繪圖
# plt.figure(figsize=(8, 6))
# plt.plot(ne_vals, ratio_vals, label='[S II] 6716 / 6731')
# plt.axhline(y=0.63, color='red', linestyle='--', label='Observed Ratio = 0.63')
# plt.xscale('log')
# plt.xlabel('Electron Density $n_e$ (cm$^{-3}$)', fontsize=12)
# plt.ylabel('Line Ratio I(6716) / I(6731)', fontsize=12)
# plt.title('[S II] Density Diagnostic at $T_e$ = 10,000 K', fontsize=14)
# plt.grid(True, which='both', ls='--', alpha=0.3)
# plt.legend()
# plt.tight_layout()
# plt.show()

# ----------------------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.optimize import fsolve
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
# from skimage import measure

# # 星體質量參數
# M1 = 1.4
# M2 = 0.75
# mu = M2 / (M1 + M2)
# x1 = -mu
# x2 = 1 - mu

# # 拉格朗日點
# def L1_eq(x): return (1 - mu)*(x + mu)/abs(x + mu)**3 + mu*(x - 1 + mu)/abs(x - 1 + mu)**3 - x
# L1_x = fsolve(L1_eq, 0.7)[0]
# L4_x, L4_y = 0.5 - mu, np.sqrt(3)/2
# L5_x, L5_y = 0.5 - mu, -np.sqrt(3)/2

# # 建立 3D 網格
# grid_pts = 100
# x = np.linspace(-2, 2, grid_pts)
# y = np.linspace(-2, 2, grid_pts)
# z = np.linspace(-1, 1, grid_pts)
# X, Y, Z = np.meshgrid(x, y, z)

# # Roche 潛勢定義
# def roche_potential_3D(x, y, z):
#     r1 = np.sqrt((x - x1)**2 + y**2 + z**2)
#     r2 = np.sqrt((x - x2)**2 + y**2 + z**2)
#     r1 = np.maximum(r1, 1e-4)
#     r2 = np.maximum(r2, 1e-4)
#     omega2 = 1
#     return -(1 - mu)/r1 - mu/r2 - 0.5 * omega2**2 * (x**2 + y**2)

# # 計算潛勢
# Phi = roche_potential_3D(X, Y, Z)

# # 提取等位面
# iso_level = -2.5
# verts, faces, _, _ = measure.marching_cubes(Phi, level=iso_level, spacing=(x[1]-x[0], y[1]-y[0], z[1]-z[0]))

# # 畫圖
# fig = plt.figure(figsize=(10, 9))
# ax = fig.add_subplot(111, projection='3d')

# # 畫等位面
# mesh = Poly3DCollection(verts[faces], alpha=0.6)
# mesh.set_facecolor('skyblue')
# ax.add_collection3d(mesh)

# # 星體位置
# ax.scatter(x1, 0, 0, color='purple', s=100, label='Neutron Star')
# ax.scatter(x2, 0, 0, color='gray', s=100, label='White Dwarf')

# # 拉格朗日點
# for label, (lx, ly, lz) in zip(['L1', 'L4', 'L5'],
#                                [(L1_x, 0, 0), (L4_x, L4_y, 0), (L5_x, L5_y, 0)]):
#     ax.scatter(lx, ly, lz, color='red', s=50)
#     ax.text(lx, ly, lz + 0.05, label, color='red')

# # 軸設定
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# ax.set_title(f'3D Equipotential Surface (Top View)')
# ax.legend()

# # 設定俯視角度
# ax.view_init(elev=90, azim=0)  # 俯視圖：從 Z 軸往下看
# ax.set_box_aspect([1, 1, 0.4])  # 避免 z 軸太扁

# plt.tight_layout()
# plt.show()

# -----------------------------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation

# # 基本設定
# G = 6.67430e-8  # 重力常數 [cm^3 g^-1 s^-2]
# M_sun = 1.989e33  # 太陽質量 [g]
# R_sun = 6.96e10  # 太陽半徑 [cm]
# L_sun = 3.828e33  # 太陽光度 [erg/s]

# # 初始條件
# m_total = 0.8  # 紅巨星總質量 (Msun)
# m_core0 = 0.15  # 初始核心質量 (Msun)
# m_core_max = 0.45  # 結束核心質量
# m_NS = 1.4  # 中子星質量 (Msun)
# a0 = 3.0  # 初始軌道半徑 (R_sun)
# timesteps = 200
# dt = 1e6  # 單位為年

# # Webbink 公式中的參數 (近似)
# def stellar_radius(mc):
#     y = np.log(4 * mc)
#     return np.exp(0.5 * y) * R_sun  # 半徑以 cm 為單位

# # 計算軌道變化 (假設角動量守恆 & 質量轉移)
# def update_orbit(a, m1, m2):
#     return a * (m1 / m2)**2  # 非保守傳質簡化模型

# # 初始化數據
# m_core_list = np.linspace(m_core0, m_core_max, timesteps)
# r_star_list = [stellar_radius(mc) for mc in m_core_list]
# a_list = [a0 * R_sun]
# m_env_list = []

# for i in range(1, timesteps):
#     m_env = m_total - m_core_list[i]
#     m_env_list.append(m_env)
#     a_new = update_orbit(a_list[-1], m_total, m_NS)
#     a_list.append(a_new)

# # 轉為 AU 和 R_sun 單位方便繪圖
# r_star_list = np.array(r_star_list) / R_sun
# a_list = np.array(a_list) / R_sun
# m_core_list = np.array(m_core_list)

# # 建立動畫
# fig, ax = plt.subplots()
# star, = plt.plot([], [], 'ro', label='紅巨星', markersize=10) # 初始標記大小
# ns, = plt.plot([], [], 'bo', label='中子星', markersize=8)
# orb, = plt.plot([], [], 'k--', linewidth=0.5)

# ax.set_xlim(-5, 5)
# ax.set_ylim(-5, 5)
# ax.set_aspect('equal')
# ax.set_title("紅巨星–中子星 雙星系統演化")
# ax.set_xlabel("軌道半徑 (太陽半徑)")
# ax.set_ylabel("軌道半徑 (太陽半徑)")
# ax.legend()

# def init():
#     star.set_data([], [])
#     ns.set_data([], [])
#     orb.set_data([], [])
#     return star, ns, orb

# def animate(i):
#     r = a_list[i]
#     rs = r_star_list[i]
#     # 紅巨星在左，中子星在右 (近似質心位置)
#     star_x = -r * m_NS / (m_total + m_NS) # 更精確的質心位置近似
#     ns_x = r * m_total / (m_total + m_NS)  # 更精確的質心位置近似
#     star.set_data(star_x, 0)
#     star.set_markersize(rs * 10) # 根據紅巨星半徑調整大小
#     ns.set_data(ns_x, 0)

#     # 繪製軌道圓圈 (以原點為中心，半徑為 a_list[i])
#     theta = np.linspace(0, 2 * np.pi, 100)
#     orb_x = a_list[i] * np.cos(theta)
#     orb_y = a_list[i] * np.sin(theta)
#     orb.set_data(orb_x, orb_y)
#     return star, ns, orb

# ani = FuncAnimation(fig, animate, frames=timesteps, init_func=init, blit=True, interval=1000)
# plt.show()

# ----------------------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt

# # 定義 alpha 範圍
# alpha = np.linspace(0, 1, 1000)

# # 定義係數 k
# k_values = [1.91e-10,9.55e-10, 1.91e-9, 9.55e-9, 3.48e-8]
# c_values = [5.206663e-8 / k for k in k_values]
# labels = [f'k = {k:.2e}' for k in k_values]
# colors = ['yellow','blue', 'red', 'green', 'purple']

# # 計算 1 - beta
# def one_minus_beta(alpha, c):
#     return (c - 0.99) / (c + 0.076667 - alpha)

# # 繪圖
# plt.figure(figsize=(10, 8))
# for c, label, color in zip(c_values, labels, colors):
#     y = one_minus_beta(alpha, c)
#     plt.plot(alpha, y, label=label, color=color)

# plt.xlabel(r'$\alpha$')
# plt.ylabel(r'$1 - \beta$')
# plt.title(r'Plot of $1 - \beta$ vs $\alpha$ for Different Coefficients')
# plt.grid(True)
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.legend()
# plt.show()
# plt.savefig('plot_different_k.png')

# -----------------------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt

# # 設定一系列的光學深度值
# tau_values = np.linspace(0, 5, 100)  # 從 0 到 5 產生 100 個均勻分布的值

# # 設定五種不同的雲溫度 (K)
# temperatures = [100, 200, 300, 400, 500]
# colors = ['r-', 'g-', 'b-', 'c-', 'm-']  # 為每條曲線設定不同的顏色和線型

# # 繪製不同溫度的曲線
# for i, T in enumerate(temperatures):
#     T_obs_values = T * (1 - np.exp(-tau_values))
#     plt.plot(tau_values, T_obs_values, colors[i], label=f'T = {T} K')

# # 添加標籤和標題
# plt.xlabel('Optical Depth ($\\tau$)')
# plt.ylabel(r'$T_{obs}$ [K]')
# plt.title(r'$T_{obs}$ vs. $\tau$')
# plt.xlim(0, 5)
# plt.ylim(0, 500)  # 根據最高溫度調整 y 軸範圍
# plt.grid(True)
# plt.legend()
# plt.show()

# ------------------------------------------

with open(f"D:/For Institute/NCU/高能天文實驗室/vscode/4U1820-30/swift/daily_lc.tsv", "r") as fobj:
    time, count, error= np.loadtxt(fobj, skiprows=2, usecols=(0, 1, 2), unpack=True)
# with open(f"D:/For Institute/NCU/高能天文實驗室/vscode/X2127+119/maxi/Orbit_lc.tsv", "r") as fobj:
#     time, count, error= np.loadtxt(fobj, skiprows=1, usecols=(0, 1, 2), unpack=True)

import math

def compute_variance(data):
    """計算一組數據的方差（無偏估計）"""
    n = len(data)
    if n <= 1:
        return 0.0
    mean = sum(data) / n
    return sum((x - mean)**2 for x in data) / (n - 1)

def pdm(time, flux, min_period, max_period, step, bins=10):
    """
    PDM方法，找出最佳週期（不用任何外部套件）
    :param time: list, 觀測時間
    :param flux: list, 觀測值（如亮度）
    :param min_period: float, 最小週期
    :param max_period: float, 最大週期
    :param step: float, 週期搜尋步長
    :param bins: int, 相位分組數
    :return: tuple (best_period, theta_values, periods)
    """
    # 產生所有要測試的週期
    periods = []
    P = min_period
    while P < max_period:
        periods.append(P)
        P += step

    theta_values = []
    total_variance = compute_variance(flux)

    for P in periods:
        # 計算每個點的相位（0~1）
        phases = [(t % P) / P for t in time]
        # 根據相位排序
        sorted_pairs = sorted(zip(phases, flux), key=lambda x: x[0])
        sorted_phases, sorted_flux = zip(*sorted_pairs)
        # 分組
        bin_edges = [i / bins for i in range(bins + 1)]
        bin_flux = [[] for _ in range(bins)]
        for phase, val in zip(sorted_phases, sorted_flux):
            for b in range(bins):
                if bin_edges[b] <= phase < bin_edges[b+1]:
                    bin_flux[b].append(val)
                    break
        # 計算各組方差加權和
        s2 = 0.0
        n_total = 0
        for group in bin_flux:
            n = len(group)
            if n > 1:
                var = compute_variance(group)
                s2 += (n - 1) * var
                n_total += (n - 1)
        # 計算theta
        if n_total > 0 and total_variance > 0:
            theta = s2 / n_total / total_variance
        else:
            theta = 1.0  # 無方差縮減
        theta_values.append(theta)

    # 找出最佳週期
    best_index = min(range(len(theta_values)), key=lambda i: theta_values[i])
    best_period = periods[best_index]

    return best_period, theta_values, periods

# --- 呼叫 PDM 函數 ---
best_period, theta_values, periods = pdm(time, count, 150, 190, 0.01)

# --- 輸出結果 ---
print("最佳週期:", best_period)

plt.figure(figsize=(10, 5))
plt.plot(periods, theta_values, label='PDM theta')
plt.xlabel('Period')
plt.ylabel('Theta')
# plt.title('PDM for 4U2127+119')
plt.title('PDM for 4U1820-30')
# plt.axvline(best_period, color='red', linestyle='--', label='Best period')
plt.legend()
plt.show()

