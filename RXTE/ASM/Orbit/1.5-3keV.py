import matplotlib.pyplot as plt
import numpy as np
import statistics
from astropy.timeseries import LombScargle
import astropy.units as u
import lightkurve
import math

START = 0
DAYS = 6000
RXTE_PERIOD = 0.713014

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

with open(f"D:/For Institute/NCU/高能天文實驗室/vscode/X2127+119/RXTE/ASM/Orbit/DWD_lc.tsv", "r") as fobj:
    time, count, error= np.loadtxt(fobj, skiprows=1, usecols=(0, 3, 4), unpack=True)

TimeList_clear,CountRateList_clear,ErrorList_clear = data_red(time,count,error)

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

def LS_power(time, count):
    # frequency = cycles / year
    frequency = np.linspace(0.001, 0.1, 2000)
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

Frequency,Power = LS_power(TimeList_clear,CountRateList_clear)  

# plt.plot(Frequency,Power,'b',label='1.5-3keV') 
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
# plt.ylabel(r'Frequency(cycles/yr)',fontsize = 12)
# plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(1))
# plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(1))
# plt.colorbar()
# plt.title('dynamic power spectrum(RXTE(1.5-3keV))')
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

    return folded_countrate_0_to_2, phase_0_to_2, bin_errors_0_to_2

Weighted_Average,Phase,Error = fold_light_curve(CountRateList_clear, TimeList_clear, ErrorList_clear, RXTE_PERIOD, num_bins=32)

plt.errorbar(Phase,Weighted_Average,yerr = Error,fmt = '^',color = 'r')
plt.title('fold light curve(RXTE)')
plt.legend(loc = 'upper right')
plt.xlabel('Phase')
plt.ylabel('Countrate')
# plt.ylim(2.3,3.6)
plt.show()  

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

    model_x = np.linspace(-1, 1, num=100000)
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

n = 5 # order

model_phase, model_count, num_data_points, exp_phase, exp_count= multiple_sinusoidal_model(Phase, Weighted_Average, Error, n)

# 计算自由度
chi_square = np.sum((Weighted_Average - exp_count)**2 / Error**2)
dof = num_data_points - (2 * n + 1)
chi_square_reduced = chi_square / dof
min_y_index = np.where((model_phase >= -0.5) & (model_phase <= 0.5))[0][np.argmin(model_count[(model_phase >= -0.5) & (model_phase <= 0.5)])]
fiducial_point_phase = model_phase[min_y_index]

# plt.errorbar(Phase, Weighted_Average, yerr=Error, fmt='.', label = f'N = {n}')
# plt.plot(model_phase, model_count, linestyle='-', color='orange')
# plt.xlabel('phase')
# plt.ylabel('count rate')
# plt.title(f'multi-sinusoidal fitting(RXTE ASM)')
# plt.legend()
# plt.show()

# print(f"Chi-Square: {chi_square:.2f}")
# print(f"DOF: {dof:.2f}")
# print(f"Reduced Chi-Square: {chi_square_reduced:.2f}")
# print(f"fiducial point phase: {fiducial_point_phase:.6f}")
# print(f"f value: {f_value:.5f}")
# print(f"p value: {p_value:.5f}")

# ----------------------------------------------------------------

def monte_carlo_simulation(count, phase, error):
    num_simulation = 10000
    simulated_phase = []

    for _ in range(num_simulation):
        count_sim = count + np.random.normal(loc=0, scale=error, size=len(count))
        min_y_index = np.where((phase >= -0.5) & (phase <= 0.5))[0][np.argmin(count_sim[(phase >= -0.5) & (phase <= 0.5)])]
        simulated_phase.append(phase[min_y_index])

    standard_deviation = np.sqrt((1 / (num_simulation - 1)) * np.sum((simulated_phase - np.mean(simulated_phase))**2))
    return simulated_phase, standard_deviation

simulated_phase, standard_deviation = monte_carlo_simulation(exp_count,Phase,Error)

def factorial(n):
    if n==0:
        return 1
    else:
        return n*factorial(n-1)

def poisson_distribution(x,mean):
    return np.exp(-mean)*mean**x/factorial(x)

# plt.hist(simulated_phase, bins=400, edgecolor='black')
# plt.xlabel('phase')
# plt.ylabel('counts/bin')
# plt.title('Monte Carlo simulation(RXTE ASM)')
# plt.legend()
# plt.show()