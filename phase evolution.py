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
from scipy.special import gamma
from scipy.stats import f
from scipy.optimize import curve_fit

# 數據點
fiducial_point_ASCA_GIS = [-0.00629]
fiducial_point_ASCA_SIS = [-0.00071]
fiducial_point_bepposax = [0.03201]
fiducial_point_Chandra = [0.09450]
fiducial_point_exosat = [-0.04147]
fiducial_point_ginga = [-0.00472]
fiducial_point_heao_HED = [-0.04402]
fiducial_point_heao_MED = [-0.05248]
fiducial_point_MAXI_2_10keV_3point = [0.11548, 0.15633, 0.17467]
fiducial_point_RXTE_ASM_3point = [0.04496, 0.06130, 0.10807]
fiducial_point_RXTE_PCA = [0.05025]
fiducial_point_XMM = [0.07617]



error_ASCA_GIS = [0.01479]
error_ASCA_SIS = [0.01131]
error_bepposax = [0.01397]
error_Chandra = [0.01331]
error_exosat = [0.01371]
error_ginga = [0.00131]
error_heao_HED = [0.01413]
error_heao_MED = [0.01453]
error_MAXI_2_10keV_3point = [0.02406, 0.03379, 0.02893]
error_RXTE_ASM_3point = [0.02335, 0.02446, 0.03046]
error_RXTE_PCA = [0.00448]
error_XMM = [0.00482]



time_mid_ASCA_GIS = [49853.6]
time_mid_ASCA_SIS = [49853.7]
time_mid_bepposax = [51498.5]
time_mid_Chandra = [54047.3]
time_mid_exosat = [46021.3]
time_mid_ginga = [47456.7]
time_mid_heao_HED = [43611.4]
time_mid_heao_MED = [43634.2]
time_mid_MAXI_2_10keV_3point = [55876.9, 57603.3, 59490.3]
time_mid_RXTE_ASM_3point = [50885.8, 52705.8, 54596.9]
time_mid_RXTE_PCA = [51063.8]
time_mid_XMM = [52708.7]


# # 繪製 heao HED 數據
# plt.errorbar(time_mid_heao_HED, fiducial_point_heao_HED, yerr=error_heao_HED, fmt='s', capsize=5, capthick=1, elinewidth=1, label='HEAO-1 HED(2.5-60keV)', color='orange')

# # 繪製 heao MED 數據
# plt.errorbar(time_mid_heao_MED, fiducial_point_heao_MED, yerr=error_heao_MED, fmt='P', capsize=5, capthick=1, elinewidth=1, label='HEAO-1 MED(1.5-20keV)', color='purple')

# # 繪製 exosat 數據
# plt.errorbar(time_mid_exosat, fiducial_point_exosat, yerr=error_exosat, fmt='h', capsize=5, capthick=1, elinewidth=1, label='EXOSAT(0.05-2keV)', color='brown')

# # 繪製 ginga 數據
# plt.errorbar(time_mid_ginga, fiducial_point_ginga, yerr=error_ginga, fmt='X', capsize=5, capthick=1, elinewidth=1, label='Ginga(1.3-37keV)', color='black')

# # 繪製 ASCA SIS 數據
# plt.errorbar(time_mid_ASCA_SIS, fiducial_point_ASCA_SIS, yerr=error_ASCA_SIS, fmt='*', capsize=5, capthick=1, elinewidth=1, label='ASCA_SIS(0.4-10keV)')

# # 繪製 ASCA GIS 數據
# plt.errorbar(time_mid_ASCA_GIS, fiducial_point_ASCA_GIS, yerr=error_ASCA_GIS, fmt='*', capsize=5, capthick=1, elinewidth=1, label='ASCA_GIS(0.7-10keV)')

# # 繪製 RXTE PCA 數據
# plt.errorbar(time_mid_RXTE_PCA, fiducial_point_RXTE_PCA, yerr=error_RXTE_PCA, fmt='^', capsize=5, capthick=1, elinewidth=1, label='RXTE PCA(2-9keV)')

# # 繪製 bepposax 數據
# plt.errorbar(time_mid_bepposax, fiducial_point_bepposax, yerr=error_bepposax, fmt='.', capsize=5, capthick=1, elinewidth=1, label='BeppoSAX(1.3-10keV)')

# # 繪製 XMM 數據
# plt.errorbar(time_mid_XMM, fiducial_point_XMM, yerr=error_XMM, fmt='v', capsize=5, capthick=1, elinewidth=1, label='XMM-Newton(1.5-12keV)')

# # 繪製 Chandra 數據
# plt.errorbar(time_mid_Chandra, fiducial_point_Chandra, yerr=error_Chandra, fmt='d', capsize=5, capthick=1, elinewidth=1, label='Chandra(0.3-10keV)')

# # 繪製 RXTE ASM 數據
# plt.errorbar(time_mid_RXTE_ASM_3point, fiducial_point_RXTE_ASM_3point, yerr=error_RXTE_ASM_3point, fmt='o', capsize=5, capthick=1, elinewidth=1, label='RXTE ASM(1.5-12keV)')

# # 繪製 MAXI 數據
# plt.errorbar(time_mid_MAXI_2_10keV_3point, fiducial_point_MAXI_2_10keV_3point, yerr=error_MAXI_2_10keV_3point, fmt='x', capsize=5, capthick=1, elinewidth=1, label='MAXI(2-10keV)', color='green')

# 合併數據
fiducial_points = fiducial_point_RXTE_PCA + fiducial_point_XMM + fiducial_point_Chandra + fiducial_point_ASCA_GIS + fiducial_point_ASCA_SIS + fiducial_point_RXTE_ASM_3point + fiducial_point_MAXI_2_10keV_3point + fiducial_point_bepposax + fiducial_point_heao_HED + fiducial_point_heao_MED + fiducial_point_ginga + fiducial_point_exosat
errors = error_RXTE_PCA + error_XMM + error_Chandra + error_ASCA_GIS + error_ASCA_SIS + error_RXTE_ASM_3point + error_MAXI_2_10keV_3point + error_bepposax + error_heao_HED + error_heao_MED + error_ginga + error_exosat
time_mids = time_mid_RXTE_PCA + time_mid_XMM + time_mid_Chandra + time_mid_ASCA_GIS + time_mid_ASCA_SIS + time_mid_RXTE_ASM_3point + time_mid_MAXI_2_10keV_3point + time_mid_bepposax + time_mid_heao_HED + time_mid_heao_MED + time_mid_ginga + time_mid_exosat

# 转换为 numpy 数组
fiducial_points = np.array(fiducial_points)
errors = np.array(errors)
time_mids = np.array(time_mids)

# import pandas as pd

# # 创建 DataFrame 并保留五位小数
# data = {
#     'phase': np.round(fiducial_points, 5),
#     'error': np.round(errors, 5),
#     'time mid point': np.round(time_mids, 5)
# }
# df = pd.DataFrame(data)

# # 按 'phase' 列从小到大排序，其他列的值也会跟着排序
# df = df.sort_values(by='phase').reset_index(drop=True)

# # 保存为 .tsv 文件
# df.to_csv('phase_evolution.tsv', sep='\t', index=False)

# ----------------------------------------------------------------

P_fold = 0.713014
T_fold = 47790.463

# fitting a stright line

# phi(t) = a_0 + a_1 * (t - T_fold)

def GLS_model(x, y, err, n):
    if np.any(err == 0):
        raise ValueError("Error array `err` contains zero values, which will cause division by zero.")
    
    num_param = n + 1
    num_data_points = len(x)
    Z = np.zeros((num_param, len(x)))
    
    # Construct the design matrix
    for k in range(0, num_param):
        Z[k, :] = x ** k

    alpha = np.zeros((num_param, num_param))
    beta = np.zeros(num_param)

    # Build the alpha and beta matrices
    for k in range(num_param):
        beta[k] = np.sum(y * Z[k] / err**2)
        for j in range(num_param):
            alpha[k, j] = np.sum(Z[k] * Z[j] / err**2)

    # Solve for the parameters using a more stable approach
    param = np.linalg.solve(alpha, beta)

    covariance_matrix = np.linalg.inv(alpha)

    # Extract diagonal elements and calculate square roots
    diag_elements = np.diag(covariance_matrix)
    sqrt_diag_elements = np.sqrt(diag_elements)
    
    sigma_0, sigma_1 = sqrt_diag_elements[:2]
    sigma_01 = covariance_matrix[0][1]

    if n == 2:
        sigma_2 = sqrt_diag_elements[2]
        sigma_02 = covariance_matrix[0][2]
        sigma_12 = covariance_matrix[1][2]

    # Generate the fitted curve
    model_x = np.linspace(np.min(x), np.max(x), num=1000)
    model_Z = np.zeros((num_param, len(model_x)))
    for k in range(0, num_param):
        model_Z[k, :] = model_x ** k

    model_y = np.dot(param, model_Z)
    
    # Calculate the fitted values for the original data
    obs_Z = np.zeros((num_param, len(x)))
    for k in range(0, num_param):
        obs_Z[k, :] = x ** k

    if n == 1:
        a0, a1 = param
        return model_x, model_y, num_data_points, a0, a1, sigma_0, sigma_1, sigma_01
    elif n == 2:
        a0, a1, a2 = param
        return model_x, model_y, num_data_points, a0, a1, a2, sigma_0, sigma_1, sigma_2, sigma_01, sigma_02, sigma_12
    else:
        raise ValueError("Currently only supports models of degree 1 or 2.")

t_shifted = time_mids - T_fold

model_time_mids, model_fiducial_point, num_data_points, a0, a1, sigma_0, sigma_1, sigma_01= GLS_model(t_shifted, fiducial_points, errors, 1)

T0 = T_fold + a0 * P_fold
sigma_T0 = sigma_0 * P_fold
P0 = P_fold / (1 - a1 * P_fold)
sigma_P0 = sigma_1 * P0**2 
sigma_T0P0 = P_fold * P0**2 * sigma_01

# 生成擬合曲線
x_fit = np.linspace(np.min(t_shifted), np.max(t_shifted), 1000)
y_fit = a0 + a1 * x_fit
chi_square1 = np.sum(((fiducial_points - a1 * t_shifted - a0)**2) / errors**2)
N1 = len(t_shifted)
p1 = 2
DOF1 = (N1 - p1)
reduced_chi_square1 = chi_square1 / (N1 - p1) 

# print(f"intercept: a = {a0} ± {sigma_0}")
# print(f"slope: b = {a1} ± {sigma_1}")
# print(f"T0 = {T0} ± {sigma_T0}") 
# print(f"P0 = {P0} ± {sigma_P0}") 
# print(f"chi square = {chi_square1}")
# print(f"DOF1 = {N - p}")
# print(f"reduced chi square = {reduced_chi_square1}")
# print(f"cov(0,1) = {sigma_01}")
# print(f"r_01 = {cor_01}")

# plt.plot(model_time_mids + T_fold, model_fiducial_point, color='blue')
# plt.title('Phase evolution(linear model)')
# plt.xlabel('time mid-point(MJD)')
# plt.ylabel(r'Phase')
# plt.legend(fontsize=8)
# plt.show()

# ----------------------------------------------------------------

# fitting a quadratic model

# phi(t) = a_0 + a_1 * (t - T_fold) + a_2 * (t - T_fold)**2

model_time_mids, model_fiducial_point, num_data_points, a0, a1, a2, sigma_0, sigma_1, sigma_2, sigma_01, sigma_02, sigma_12= GLS_model(t_shifted, fiducial_points, errors, 2)

T0 = T_fold + a0 * P_fold
P0 = P_fold / (1 - a1 * P_fold)
Pdot = 2 * a2 * P0**2
sigma_T0 = sigma_0 * P_fold
sigma_P0 = sigma_1 * P0**2
sigma_Pdot = np.sqrt(4 * P0**4 * (sigma_2**2 + 2 * a2 * P0 * sigma_12 + 4 * a2**2 * P0**2 * sigma_1**2))
sigma_P0_Pdot = 2 * P_fold * P0**3 * (a2 * P0 * sigma_1**2 + sigma_12)
sigma_T0_P0 = P_fold * P0**2 * sigma_01 
sigma_T0_Pdot = 2 * P_fold**2 * P0 * (2 * a2 * P0 * sigma_01 + sigma_02)
sigma_1_2_P0_Pdot = np.sqrt(1 / 4 * (Pdot**2 * sigma_P0**2 + 2 * P0 * Pdot * sigma_P0_Pdot + P0**2 * sigma_Pdot**2))
sigma_Pdot_over_P0 = np.sqrt((1 / P0)**2 * sigma_Pdot**2 + (-Pdot / P0**2)**2 * sigma_P0**2 + 2 * (1 / P0) * (-Pdot / P0**2) * sigma_P0_Pdot)

# print(Pdot, sigma_Pdot)

# 生成擬合曲線
x_fit = np.linspace(np.min(t_shifted), np.max(t_shifted), 1000)
y_fit = a0 + a1 * x_fit + a2 * x_fit**2
chi_square2 = np.sum(((fiducial_points - a0 - a1 * t_shifted - a2 * t_shifted**2)**2) / errors**2)
N2 = len(t_shifted)
p2 = 3
DOF2 = (N2 - p2)
reduced_chi_square2 = chi_square2 / DOF2
F = ((chi_square1 - chi_square2) / (DOF1 - DOF2)) / (chi_square2 / DOF2)

# print(f"constant term: a0 = {a0} ± {sigma_0}")
# print(f"first order term: a1 = {a1} ± {sigma_1}")
# print(f"second order term: a2 = {a2} ± {sigma_2}")
# print(f"covariance: cov(0,1) = {sigma_01}, cov(0,2) = {sigma_02}, cov(1,2) = {sigma_12}")
# print(f"T_0 = {T0} ± {sigma_T0}")
# print(f"P_0 = {P0} ± {sigma_P0}")
# print(f"P_dot = {Pdot} ± {sigma_Pdot}")
# print(f"1/2 P0_Pdot = {1 / 2 * P0 * Pdot} ± {sigma_1_2_P0_Pdot}")
# print(f"Pdot_over_P0 = {Pdot / P0} ± {sigma_Pdot_over_P0}")
# print(f"chi square 1 = {chi_square1}")
# print(f"chi square 2 = {chi_square2}")
# print(f"DOF2 = {DOF2}")
# print(f"reduced chi square = {reduced_chi_square2}")
# print(f"f value = {F}")

error_p_dot_over_p = np.sqrt(Pdot**2 / P0**4 * sigma_P0**2 + 1 / P0**2 * sigma_Pdot**2)

# # 可视化结果
# plt.plot(x_fit + T_fold, y_fit, color='blue')
# plt.title('Phase evolution(quadratic model)')
# plt.xlabel('time mid-point(MJD)')
# plt.ylabel(r'Phase')
# plt.legend(fontsize=8)
# plt.show()

# ----------------------------------------------------------------

# F-distribution

# 定义概率密度函数
def probability_of_f(f_value, nu_1, nu_2):
    return (gamma((nu_1 + nu_2) / 2) / (gamma(nu_1 / 2) * gamma(nu_2 / 2)) *
            (nu_1 / nu_2)**(nu_1 / 2) * f_value**(nu_1 / 2 - 1) *
            (1 + (nu_1 / nu_2) * f_value)**(-(nu_1 + nu_2) / 2))

def probability_of_f_1(f_value, mu_2):
    return gamma((mu_2 + 1) / 2) / (gamma(1 / 2) * gamma(mu_2 / 2)) * mu_2**(-1 / 2) * f_value**(-1 / 2) * (1 + f_value / mu_2)**(-(mu_2 + 1) / 2)

# 使用梯形积分法计算累积分布函数
def cumulative_density_function(f_value, mu_2):
    num_points=1000
    x = np.linspace(0, f_value, num_points)
    y = probability_of_f_1(x, mu_2)
    
    # 實現梯形積分法
    def trapezoidal_rule(y, a, b, num_points):
        h = (b - a) / (num_points - 1)
        integral = h * (0.5 * y[0] + 0.5 * y[-1] + np.sum(y[1:-1]))
        return integral
    
    # 实现辛普森积分法
    def simpsons_rule(y, a, b, num_points):
        if num_points % 2 == 1:
            num_points += 1  # 辛普森积分法要求偶数个子区间
        h = (b - a) / (num_points - 1)
        integral = (h / 3) * (y[0] + y[-1] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-2:2]))
        return integral
    
    cdf_trapz = trapezoidal_rule(y, 0, f_value, num_points)
    cdf_simps = simpsons_rule(y, 0, f_value, num_points)
    return cdf_trapz, cdf_simps

# 计算累积分布函数值
cdf_trapz, cdf_simps = cumulative_density_function(F, DOF2)
p_value_trapz = 1 - cdf_trapz
p_value_simps = 1 - cdf_simps
p_value = 1 - f.cdf(F, DOF1 - DOF2, DOF2)
print(f'p_value = {p_value}')

# F 分布的 x 轴范围
x = np.linspace(0, 10, 1000)
pdf = probability_of_f(x, DOF1, DOF2)
pdf1 = probability_of_f_1(x, DOF2)

# # 绘制PDF和CDF
# plt.figure(figsize=(12, 6))

# plt.subplot(1, 2, 1)
# plt.plot(x, pdf1, label=f'F-distribution (df1={DOF1}, df2={DOF2})')
# plt.axvline(F, color='r', linestyle='--', label=f'F-value = {F:.3f}')
# plt.fill_between(x, 0, pdf1, where=(x >= F), color='red', alpha=0.3) 
# plt.title('F-distribution with F-value')
# plt.xlabel('F')
# plt.ylabel('Density')
# plt.legend()

plt.plot(x, pdf)
plt.axvline(F, color='r', linestyle='--', label=f'F-value = {F:.3f}')
plt.fill_between(x, 0, pdf, where=(x >= F), color='red', alpha=0.3) 
plt.title('F-distribution')
plt.xlabel('F')
plt.ylabel('Density')
plt.legend()
plt.show()

# ----------------------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.stats import f  # 從 SciPy 匯入 F-分佈

# # --- 請在此處設定您的參數 ---
# # 這是您觀察到的 F 統計量值
# F_statistic = F
# # 自由度 (Numerator and Denominator)
# DOF1 = (N1 - p1)
# DOF2 = (N2 - p2)
# # --------------------------

# # 1. 計算 p-value
# # 使用 f.sf() (Survival Function) 來計算 p-value，它相當於 1 - CDF
# # 這在數值上比 1 - f.cdf() 更精確
# p_value = f.sf(F_statistic, DOF1, DOF2)
# print(f"F-statistic = {F_statistic:.4f}")
# print(f"Degrees of Freedom = ({DOF1}, {DOF2})")
# print(f"p-value = {p_value:.4f}")

# # 2. 準備繪圖數據
# # 創建一個 x 軸的範圍，從 0 到一個適當的值
# x = np.linspace(0, 10, 1000)

# # *** 使用 SciPy 的 f.pdf() 函數計算機率密度 ***
# # 這取代了您自定義的 probability_of_f 函數
# pdf = f.pdf(x, DOF1, DOF2)

# # 3. 繪製圖形
# plt.figure(figsize=(10, 6)) # 創建一個圖形物件

# # 繪製 F-分佈的 PDF 曲線
# plt.plot(x, pdf, 'b-', label=f'F-distribution (df1={DOF1}, df2={DOF2})')

# # 標示出您的 F 統計量值
# plt.axvline(F_statistic, color='r', linestyle='--', label=f'F-value = {F_statistic:.3f}')

# # 填充 p-value 對應的右尾區域
# # 我們需要定義一個 x 的子集，從 F_statistic 開始
# x_fill = np.linspace(F_statistic, 10, 500)
# y_fill = f.pdf(x_fill, DOF1, DOF2)
# plt.fill_between(x_fill, y_fill, color='red', alpha=0.3, label=f'p-value area = {p_value:.3f}')

# # 設定圖形標題和軸標籤
# plt.title('F-distribution PDF')
# plt.xlabel('F-value')
# plt.ylabel('Probability Density')
# plt.legend() # 顯示圖例
# plt.grid(True, linestyle='--', alpha=0.6) # 添加網格線
# plt.show()