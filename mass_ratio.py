# import numpy as np
# import matplotlib.pyplot as plt

# # 定義 alpha 範圍
# alpha = np.linspace(0, 1, 1000)

# # 定義係數 k
# k_values = [9.55e-9, 1.91e-8, 2.87e-8, 3.48e-8]
# c_values = [5.206663e-8 / k for k in k_values]
# labels = [r'$L_x=5 \times 10^{37}ergs/s$', r'$L_x=10 \times 10^{37}ergs/s$', r'$L_x=15 \times 10^{37}ergs/s$', r'$L_x=18.2 \times 10^{37}ergs/s(Eddington)$']
# colors = ['blue', 'red', 'green', 'purple']

# # 計算 1 - beta
# def one_minus_beta(alpha, c):
#     return (c - 0.99) / (c + 0.076667 - alpha)

# # 繪圖
# plt.figure(figsize=(12, 8))  # 稍微加寬圖形以容納參數
# for c, label, color in zip(c_values, labels, colors):
#     y = one_minus_beta(alpha, c)
#     plt.plot(alpha, y, label=label, color=color)

# # 添加基本參數，使用 M_\odot 替代 solar mass
# parameters = (
#     r'$M_1 = 1.4 \, M_\odot$'
#     '\n' + r'$R_1 = 10 \, \text{km}$'
#     '\n' + r'$q = \frac{0.14}{1.4}$'
#     '\n' + r'$\dot{P}_{\text{orb}} / P_{\text{orb}} = 1.42 \times 10^{-7} \, \text{yr}^{-1}$'
# )
# plt.text(0.7, 0.1, parameters, transform=plt.gca().transAxes, fontsize=12, 
#          verticalalignment='center', horizontalalignment='left',
#          bbox=dict(facecolor='white', edgecolor='k', alpha=1.0))

# plt.xlabel(r'$\alpha$')
# plt.ylabel(r'$1 - \beta$')
# plt.title(r'When the companion star $M_2 = 0.14M_\odot$(This diagram is wrong)')
# plt.grid(True)
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.legend(loc='lower left')
# plt.show()

# -----------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt

# # 常數
# L_x_values = [5e30, 1e31, 1.5e31, 1.82e31]
# L_x_labels = [r'$L_x=5 \times 10^{37}ergs/s$', r'$L_x=10 \times 10^{37}ergs/s$', r'$L_x=15 \times 10^{37}ergs/s$', r'$L_x=18.2 \times 10^{37}ergs/s(Eddington)$']
# constant_term = 4.466e31
# beta_values = np.linspace(0.01, 0.99, 500)
# colors = ['blue', 'red', 'green', 'purple']

# # 繪圖
# plt.figure(figsize=(10, 8))
# for L_x, label, color in zip(L_x_values, L_x_labels, colors):
#     alpha_values = 1.5714 * (1 - 0.5714 * beta_values - constant_term * beta_values / L_x) / (1 - beta_values) - 0.1905
#     one_minus_beta_values = 1 - beta_values
#     plt.plot(alpha_values, one_minus_beta_values, label=label, color=color)

# # 添加基本參數，使用 M_\odot 替代 solar mass
# parameters = (
#     r'$M_1 = 1.4 \, M_\odot$'
#     '\n' + r'$R_1 = 10 \, \text{km}$'
#     '\n' + r'$q = \frac{0.8}{1.4}$'
#     '\n' + r'$\dot{P}_{\text{orb}} / P_{\text{orb}} = 1.42 \times 10^{-7} \, \text{yr}^{-1}$'
# )
# plt.text(0.65, 0.15, parameters, transform=plt.gca().transAxes, fontsize=12, 
#          verticalalignment='center', horizontalalignment='left',
#          bbox=dict(facecolor='white', edgecolor='k', alpha=1.0))

# # 設定軸標籤和標題
# plt.xlabel(r'$\alpha$')
# plt.ylabel(r'$1 - \beta$')
# plt.title(r'When the companion star $M_2 = 0.8M_\odot$')
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.grid(True)
# plt.legend()
# plt.show()

# ----------------------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt

# 常數
# L_x_values = [5e30, 1e31, 1.5e31, 1.82e31]
# L_x_labels = [r'$L_x=5 \times 10^{37}ergs/s$', r'$L_x=10 \times 10^{37}ergs/s$', r'$L_x=15 \times 10^{37}ergs/s$', r'$L_x=18.2 \times 10^{37}ergs/s(Eddington)$']
# constant_term = 2.769e31
# beta_values = np.linspace(0.01, 0.99, 500)
# colors = ['blue', 'red', 'green', 'purple']

# # 繪圖
# plt.figure(figsize=(10, 8))
# for L_x, label, color in zip(L_x_values, L_x_labels, colors):
#     alpha_values = 1.1 * (1 - 0.1 * beta_values - constant_term * beta_values / L_x) / (1 - beta_values) - 0.03333
#     one_minus_beta_values = 1 - beta_values
#     plt.plot(alpha_values, one_minus_beta_values, label=label, color=color)

# # 添加基本參數，使用 M_\odot 替代 solar mass
# parameters = (
#     r'$M_1 = 1.4 \, M_\odot$'
#     '\n' + r'$R_1 = 10 \, \text{km}$'
#     '\n' + r'$q = \frac{0.14}{1.4}$'
#     '\n' + r'$\dot{P}_{\text{orb}} / P_{\text{orb}} = 5 \times 10^{-7} \, \text{yr}^{-1}$'
# )
# plt.text(0.7, 0.15, parameters, transform=plt.gca().transAxes, fontsize=12, 
#          verticalalignment='center', horizontalalignment='left',
#          bbox=dict(facecolor='white', edgecolor='k', alpha=1.0))

# # 設定軸標籤和標題
# plt.xlabel(r'$\alpha$')
# plt.ylabel(r'$1 - \beta$')
# plt.title(r'When the companion star $M_2 = 0.14M_\odot$')
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.grid(True)
# plt.legend()
# plt.show()

# --------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt

# # 常數
# L_x_values = [5e30, 1e31, 1.5e31, 1.82e31]
# L_x_labels = [r'$L_x=5 \times 10^{37}ergs/s$', r'$L_x=10 \times 10^{37}ergs/s$', r'$L_x=15 \times 10^{37}ergs/s$', r'$L_x=18.2 \times 10^{37}ergs/s(Eddington)$']
# constant_term = 1.572e32
# beta_values = np.linspace(0.01, 0.99, 500)
# colors = ['blue', 'red', 'green', 'purple']

# # 繪圖
# plt.figure(figsize=(10, 8))
# for L_x, label, color in zip(L_x_values, L_x_labels, colors):
#     alpha_values = 1.5714 * (1 - 0.5714 * beta_values - constant_term * beta_values / L_x) / (1 - beta_values) - 0.1905
#     one_minus_beta_values = 1 - beta_values
#     plt.plot(alpha_values, one_minus_beta_values, label=label, color=color)

# # 添加基本參數，使用 M_\odot 替代 solar mass
# parameters = (
#     r'$M_1 = 1.4 \, M_\odot$'
#     '\n' + r'$R_1 = 10 \, \text{km}$'
#     '\n' + r'$q = \frac{0.8}{1.4}$'
#     '\n' + r'$\dot{P}_{\text{orb}} / P_{\text{orb}} = 5 \times 10^{-7} \, \text{yr}^{-1}$'
# )
# plt.text(0.7, 0.15, parameters, transform=plt.gca().transAxes, fontsize=12, 
#          verticalalignment='center', horizontalalignment='left',
#          bbox=dict(facecolor='white', edgecolor='k', alpha=1.0))

# # 設定軸標籤和標題
# plt.xlabel(r'$\alpha$')
# plt.ylabel(r'$1 - \beta$')
# plt.title(r'When the companion star $M_2 = 0.8M_\odot$')
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.grid(True)
# plt.legend()
# plt.show()

# ----------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt
# from astropy import constants as const
# from astropy import units as u
# from matplotlib.patches import FancyArrowPatch

# # 參數 (使用年作為時間單位)
# M1_solar = 1.4
# P_yr = 61605 / (365.25 * 24 * 3600)
# Lint = np.array([0.5, 1, 1.5, 1.82]) * 1e38  # 單位: erg/s，換算為真實數值
# Lint_labels = [r'$L_{int} = 0.5 \times 10^{38} \, \text{erg/s}$',
#                r'$L_{int} = 1.0 \times 10^{38} \, \text{erg/s}$',
#                r'$L_{int} = 1.5 \times 10^{38} \, \text{erg/s}$',
#                r'$L_{int} = 1.82 \times 10^{38} \, \text{erg/s}(L_{Edd})$'] # 使用近似愛丁頓光度標籤
# colors = ['blue', 'green', 'purple', 'red']
# sigma_P_dot_P = 1.01e-7

# # 轉換單位 (只需要質量到 kg)
# M_sun_kg = const.M_sun.value
# M1_kg = M1_solar * M_sun_kg
# R_1_m = 10000  # 轉換為米

# # 物理常數 (需要使用 SI 單位)
# G = const.G.value
# c = const.c.value

# # q 的範圍
# q_values = np.linspace(0.01, 1, 500)

# M2_dot_solar_per_yr_values = - (Lint * 1e-7 * R_1_m / (G * M1_kg * M_sun_kg) * 86400 * 365).astype(np.float64)

# # 繪圖
# plt.figure(figsize=(10, 6))

# for i, M2_dot_solar_per_yr in enumerate(M2_dot_solar_per_yr_values):
#     P_dot_over_P_values_per_yr = []
#     for q in q_values:
#         M2_kg = q * M1_kg
#         M_total_kg = M1_kg + M2_kg
#         P_sec = P_yr * (365.25 * 24 * 3600) # 將週期轉換回秒以代入 GR 公式
#         term_gr_over_P_per_sec = - (96 * (4 * np.pi**2)**(4/3) * G**(5/3) * q * M1_kg**(4/3) * (1 + q)**(-1/3)) / (5 * c**5 * P_sec**(8/3))
#         term_gr_over_P_per_yr = term_gr_over_P_per_sec * (365.25 * 24 * 3600) # 轉換為 per year

#         M2_dot_over_M1_per_yr = M2_dot_solar_per_yr / M1_solar
#         term_mt_over_P_per_yr = - 3 * (M2_dot_over_M1_per_yr / q) * (1 - q)

#         P_dot_over_P_per_yr = term_gr_over_P_per_yr + term_mt_over_P_per_yr
#         P_dot_over_P_values_per_yr.append(P_dot_over_P_per_yr)

#     # 判斷是否是艾丁頓光度
#     if i == 3:  # 第四條，也就是 Eddington 曲線
#         linestyle = '-'   # 實線
#     else:
#         linestyle = '--'  # 虛線

#     plt.plot(q_values, P_dot_over_P_values_per_yr, linestyle=linestyle, label = Lint_labels[i])

# # 黑色水平線標註
# plt.axhline(y=1.42e-7, color='black', linestyle='-', alpha=0.7, label=r'$(\dot{P} / P)_{obs}$')
# plt.text(0.8, 1.5e-7, r'$1.42 \times 10^{-7}$', color='black', fontsize=14)

# plt.axhline(y=1.42e-7 + sigma_P_dot_P, color='black', linestyle='-.', alpha=0.5)
# plt.text(0.8, 2.5e-7, r'$2.43 \times 10^{-7}$', color='black', fontsize=14)

# plt.axhline(y=1.42e-7 - sigma_P_dot_P, color='black', linestyle='-.', alpha=0.5)
# plt.text(0.8, 5e-8, r'$4.1 \times 10^{-8}$', color='black', fontsize=14)

# # 假設畫一個雙向箭頭：從 (x1, y1) 到 (x2, y2)
# x1, y1 = 0.7, 1.42e-7
# x2, y2 = 0.7, 1.42e-7 + sigma_P_dot_P
# x3, y3 = 0.7, 1.42e-7 - sigma_P_dot_P

# arrow1 = FancyArrowPatch(
#     (x1, y1), (x2, y2),
#     arrowstyle='<->',           # 雙向箭頭
#     mutation_scale=20,          # 箭頭大小，可調整
#     color='black',                # 顏色
#     linewidth=1
# )
# arrow2 = FancyArrowPatch(
#     (x1, y1), (x3, y3),
#     arrowstyle='<->',           # 雙向箭頭
#     mutation_scale=20,          # 箭頭大小，可調整
#     color='black',                # 顏色
#     linewidth=1
# )
# plt.gca().add_patch(arrow1)
# plt.text(0.72, (y1+y2)/2, r'$1\sigma$', color='black', ha='center', va='top') # 調整 y 座標和對齊方式
# plt.gca().add_patch(arrow2)  # 加上這一行
# plt.text(0.72, (y1+y3)/2, r'$1\sigma$', color='black', ha='center', va='top') # 調整 y 座標和對齊方式

# # # plt.plot(0.14/1.4, 2.96e-7, 'ko', markersize=6, label=r'$Mass \, ratio \, \frac{M_2}{M_1} = \frac{0.14}{1.4}$') # 'ko' 表示黑色圓圈
# # # plt.plot(0.8/1.4, 2.47e-8, 'ms', markersize=6, label=r'$Mass \, ratio \, \frac{M_2}{M_1} = \frac{0.8}{1.4}$') # 'ms' 表示洋紅色方塊
# # # plt.plot(0.14/1.4, 2.96e-7, '^', markersize=6, label=r'$Mass \, ratio \, \frac{M_2}{M_1} = \frac{0.14}{1.4} when \beta \approx 1$')
# # # plt.plot(0.8/1.4, 2.5e-8, '*', markersize=6, label=r'$Mass \, ratio \, \frac{M_2}{M_1} = \frac{0.8}{1.4} when \beta \approx 1$')
# plt.xlim(0, 1)
# plt.ylim(0, 5e-7)
# plt.xlabel(r'$q = M_2 / M_1$')
# plt.ylabel(r'$\dot{P} / P$ (yr$^{-1}$)')
# plt.title('Mass-conserved case')
# plt.grid(True)
# plt.yticks([1e-7, 2e-7, 3e-7, 4e-7]) # 設定 Y 軸的刻度值
# plt.legend(loc='upper right')
# plt.show()

# ----------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt

# # 定義 alpha 範圍
# alpha = np.linspace(0, 1, 1000)

# # 定義不同的 P_dot/P 值 (yr^-1)
# p_dot_p_values = [2.96e-7, 5e-7, 1e-6, 1e-5]
# colors = ['blue', 'red', 'green', 'purple']
# labels = [r'$\dot{P}_{orb} / P_{orb} = 2.96 \times 10^{-7} \, yr^{-1}$',
#           r'$\dot{P}_{orb} / P_{orb} = 5 \times 10^{-7} \, yr^{-1}$',
#           r'$\dot{P}_{orb} / P_{orb} = 10^{-6} \, yr^{-1}$',
#           r'$\dot{P}_{orb} / P_{orb} = 10^{-5} \, yr^{-1}$']

# # 計算 1 - beta
# def one_minus_beta_derived(alpha, p_dot_p):
#     C = 3.3466e6 * p_dot_p
#     return (0.99 - C) / (alpha - 0.076667 - C)

# # 繪圖
# plt.figure(figsize=(12, 8))
# for p_dot_p, color, label in zip(p_dot_p_values, colors, labels):
#     y = one_minus_beta_derived(alpha, p_dot_p)
#     plt.plot(alpha, y, label=label, color=color)

# # 添加基本參數
# parameters = (
#     r'$L_x = L_{Edd} \approx 1.82 \times 10^{38} \, ergs/s$'
#     '\n' + r'$M_1 = 1.4 \, M_\odot$'
#     '\n' + r'$R_1 = 10 \, \text{km}$'
#     '\n' + r'$q = 0.1$'
# )
# plt.text(0.7, 0.1, parameters, transform=plt.gca().transAxes, fontsize=12,
#          verticalalignment='center', horizontalalignment='left',
#          bbox=dict(facecolor='white', edgecolor='k', alpha=1.0))

# plt.xlabel(r'$\alpha$')
# plt.ylabel(r'$1 - \beta$')
# plt.title(r'Mass-nonconserved case')
# plt.grid(True)
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.legend(loc='lower left')
# plt.show()

# ----------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt

# # 定義 alpha 範圍
# alpha = np.linspace(0, 1, 1000)

# # 固定艾丁頓光度 (不需要直接使用，因為已經在公式推導中)

# # 定義不同的 P_dot/P 值 (yr^-1)
# p_dot_p_values = [2.5e-8, 1.42e-7, 1e-6, 1e-5]
# colors = ['blue', 'red', 'green', 'purple']
# labels = [r'$\dot{P}_{orb} / P_{orb} = 2.5 \times 10^{-8} \, yr^{-1}$',
#           r'$\dot{P}_{orb} / P_{orb} = 1.42 \times 10^{-7} \, yr^{-1}$',
#           r'$\dot{P}_{orb} / P_{orb} = 10^{-6} \, yr^{-1}$',
#           r'$\dot{P}_{orb} / P_{orb} = 10^{-5} \, yr^{-1}$']

# # 計算 1 - beta
# def one_minus_beta_derived(alpha, p_dot_p):
#     C = 2.714e7 * p_dot_p
#     return (0.6735 - C) / (alpha - 0.7073 - C)

# # 繪圖
# plt.figure(figsize=(12, 8))
# for p_dot_p, color, label in zip(p_dot_p_values, colors, labels):
#     y = one_minus_beta_derived(alpha, p_dot_p)
#     plt.plot(alpha, y, label=label, color=color)

# # 添加基本參數
# parameters = (
#     r'$L_x = L_{Edd} \approx 1.82 \times 10^{38} \, ergs/s$'
#     '\n' + r'$M_1 = 1.4 \, M_\odot$'
#     '\n' + r'$R_1 = 10 \, \text{km}$'
#     '\n' + r'$q = \frac{0.8}{1.4}$'
# )
# plt.text(0.7, 0.1, parameters, transform=plt.gca().transAxes, fontsize=12,
#          verticalalignment='center', horizontalalignment='left',
#          bbox=dict(facecolor='white', edgecolor='k', alpha=1.0))

# plt.xlabel(r'$\alpha$')
# plt.ylabel(r'$1 - \beta$')
# plt.title(r'Mass-nonconserved case')
# plt.grid(True)
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.legend(loc='lower left')
# plt.show()

# ----------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt

# # 資料
# L_int = [5, 10, 15, 18.2]  # 單位: 10^38 ergs/s
# q =   [0.06, 0.113, 0.161, 0.188]

# # 誤差棒上下限「位置」
# yerr_lower_pos = [0.036, 0.069, 0.1, 0.119]
# yerr_upper_pos = [0.181, 0.306, 0.398, 0.446]

# # 轉換為「相對中心值的誤差長度」
# yerr_lower = [q[i] - yerr_lower_pos[i] for i in range(len(q))]
# yerr_upper = [yerr_upper_pos[i] - q[i] for i in range(len(q))]
# yerr = [yerr_lower, yerr_upper]

# for i in range(len(L_int)):
#     if i == len(L_int) - 1:  # 判斷是否為最後一個點
#         plt.errorbar(L_int[i], q[i],
#                      yerr=[[yerr_lower[i]], [yerr_upper[i]]],
#                      fmt='*', capsize=5, color='red')  # 紅色星點
#     else:
#         plt.errorbar(L_int[i], q[i],
#                      yerr=[[yerr_lower[i]], [yerr_upper[i]]],
#                      fmt='o', capsize=5)  # 其他點為藍色圓點

# plt.text(3, 0.07, r'$L_{37} = 5$', color='black')
# plt.text(7.5, 0.12, r'$L_{37} = 10$', color='black')
# plt.text(12.5, 0.17, r'$L_{37} = 15$', color='black')
# plt.text(15.5, 0.2, r'$L_{37} = 18.2$', color='black')

# plt.axhline(y=0.8/1.4, color='black', linestyle='-.', alpha=0.5)
# plt.text(0.02, 0.74/1.4, r'$q=0.57$', color='black')
# plt.axhline(y=0.14/1.4, color='black', linestyle='-.', alpha=0.5)
# plt.text(0.02, 0.18/1.4, r'$q=0.1$', color='black')
# plt.xlabel(r'$L_{int} \, (10^{37} \, \mathrm{erg/s})$')
# plt.ylabel('mass ratio  $q$')
# plt.xlim(0, 20)
# plt.ylim(0, 0.7)
# plt.legend(loc='upper right')
# plt.title(r'Mass-conserved case')
# plt.grid(True)
# plt.tight_layout()
# plt.show()

# -----------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt

# # 定義 alpha 範圍
# alpha = np.linspace(0, 1, 1000)

# # 固定 P_dot/P 值 (yr^-1)
# p_dot_p_fixed = 1.42e-7

# # 定義不同的光度值 (單位 10^37 ergs/s) 和對應的 c 值
# lx_values = [5, 10, 15, 18.2]
# c_values = [5.54, 2.72, 1.82, 1.496]
# labels = [r'$L_x = 5 \times 10^{37} \, ergs/s$',
#           r'$L_x = 10 \times 10^{37} \, ergs/s$',
#           r'$L_x = 15 \times 10^{37} \, ergs/s$',
#           r'$L_x = 18.2 \times 10^{37} \, ergs/s (L_{Edd})$']
# colors = ['blue', 'red', 'green', 'purple']

# # 計算 1 - beta
# def one_minus_beta(alpha, c):
#     return (c - 0.99) / (c + 0.076667 - alpha)

# # 繪圖
# plt.figure(figsize=(12, 8))
# for c, label, color in zip(c_values, labels, colors):
#     y = one_minus_beta(alpha, c)
#     plt.plot(alpha, y, label=label, color=color)

# # 添加基本參數
# parameters = (
#     r'$\dot{P}_{\text{orb}} / P_{\text{orb}} = 1.42 \times 10^{-7} \, yr^{-1}$'
#     '\n' + r'$M_1 = 1.4 \, M_\odot$'
#     '\n' + r'$R_1 = 10 \, \text{km}$'
#     '\n' + r'$q = \frac{0.14}{1.4}$'
# )
# plt.text(0.7, 0.1, parameters, transform=plt.gca().transAxes, fontsize=12,
#          verticalalignment='center', horizontalalignment='left',
#          bbox=dict(facecolor='white', edgecolor='k', alpha=1.0))

# plt.xlabel(r'$\alpha$')
# plt.ylabel(r'$1 - \beta$')
# plt.title(r'When the companion star $M_2 = 0.14M_\odot$ (at fixed $\dot{P}_{orb}/P_{orb}$)')
# plt.grid(True)
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.legend(loc='lower left')
# plt.show()

# -----------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt

# P_dot_over_P = 2.43e-7
# # P_dot_over_P = 1.42e-7
# # P_dot_over_P = 4.1e-8
# q = 0.14/1.4
# # q = 0.8/1.4
# # 常數
# L_x_values = [5e30, 1e31, 1.5e31, 1.82e31]
# L_x_labels = [r'$L_x=5 \times 10^{37}ergs/s$', r'$L_x=10 \times 10^{37}ergs/s$', r'$L_x=15 \times 10^{37}ergs/s$', r'$L_x=18.2 \times 10^{37}ergs/s(Eddington)$']
# constant_term = ((1/3) * P_dot_over_P + 1.26e-11 * q * (1 + q)**(-1/3)) * 1.66e39 * q
# beta_values = np.linspace(0.0001, 0.9999, 500)
# colors = ['blue', 'green', 'purple', 'red']

# # 繪圖
# plt.figure(figsize=(10, 8))
# for L_x, label, color in zip(L_x_values, L_x_labels, colors):
#     alpha_values = (1 + q) * (1 - q * beta_values - constant_term * beta_values / L_x) / (1 - beta_values) - q/3
#     one_minus_beta_values = 1 - beta_values
#     plt.plot(alpha_values, one_minus_beta_values, label=label, color=color)

# # 添加基本參數，使用 M_\odot 替代 solar mass
# parameters = (
#     r'$M_1 = 1.4 \, M_\odot$'
#     '\n' + r'$R_1 = 10 \, \text{km}$'
#     '\n' + r'$q = \frac{0.14}{1.4}$'
#     '\n' + r'$\dot{P}_{\text{orb}} / P_{\text{orb}} = 2.43 \times 10^{-7} \, \text{yr}^{-1}$')

# plt.text(0.7, 0.15, parameters, transform=plt.gca().transAxes, fontsize=12, 
#          verticalalignment='center', horizontalalignment='left',
#          bbox=dict(facecolor='white', edgecolor='k', alpha=1.0))

# # 設定軸標籤和標題
# plt.xlabel(r'$\alpha$')
# plt.ylabel(r'$1 - \beta$')
# plt.title('Mass non-conserved case')
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.grid(True)
# plt.legend()
# plt.show()

# --------------------------------------------------------

# import numpy as np
# from scipy.interpolate import interp1d
# from astropy import constants as const
# from astropy import units as u

# # 參數 (與原程式碼相同)
# M1_solar = 1.4
# P_yr = 61605 / (365.25 * 24 * 3600)
# Lacc = np.arange(0.5, 1.82, 0.1) * 1e38  # erg/s
# M_sun_kg = const.M_sun.value
# M1_kg = M1_solar * M_sun_kg
# R_1_m = 10000
# G = const.G.value
# c = const.c.value
# q_values = np.linspace(0.01, 1, 500)

# # 計算質量傳輸率
# M2_dot_solar_per_yr_values = - (Lacc * 1e-7 * R_1_m / (G * M1_kg * M_sun_kg) * 86400 * 365).astype(np.float64)

# # 目標週期變化率
# target_P_dot_P = np.array([0.41e-7, 1.42e-7, 2.43e-7])

# # 儲存結果
# results = {L: {t: None for t in target_P_dot_P} for L in Lacc}

# # 計算每條曲線的 q 值
# for i, (L_acc, M2_dot_solar_per_yr) in enumerate(zip(Lacc, M2_dot_solar_per_yr_values)):
#     P_dot_over_P_values_per_yr = []
#     for q in q_values:
#         M2_kg = q * M1_kg
#         M_total_kg = M1_kg + M2_kg
#         P_sec = P_yr * (365.25 * 24 * 3600)
#         term_gr_over_P_per_sec = - (96 * (4 * np.pi**2)**(4/3) * G**(5/3) * q * M1_kg**(4/3) * (1 + q)**(-1/3)) / (5 * c**5 * P_sec**(8/3))
#         term_gr_over_P_per_yr = term_gr_over_P_per_sec * (365.25 * 24 * 3600)
#         M2_dot_over_M1_per_yr = M2_dot_solar_per_yr / M1_solar
#         term_mt_over_P_per_yr = - 3 * (M2_dot_over_M1_per_yr / q) * (1 - q)
#         P_dot_over_P_per_yr = term_gr_over_P_per_yr + term_mt_over_P_per_yr
#         P_dot_over_P_values_per_yr.append(P_dot_over_P_per_yr)

#     # 插值找到目標 P_dot/P 對應的 q
#     P_dot_over_P_values_per_yr = np.array(P_dot_over_P_values_per_yr)
#     interp_func = interp1d(P_dot_over_P_values_per_yr, q_values, bounds_error=False, fill_value=np.nan)
    
#     for target in target_P_dot_P:
#         q_interpolated = interp_func(target)
#         if not np.isnan(q_interpolated) and 0.01 <= q_interpolated <= 1:
#             results[L_acc][target] = q_interpolated

# # 輸出結果
# print("不同 L_acc 在目標 P_dot/P 下的 q 值：")
# for L_acc in Lacc:
#     print(f"\nL_acc = {L_acc / 1e38:.1f} x 10^38 erg/s:")
#     for target in target_P_dot_P:
#         q_val = results[L_acc][target]
#         if q_val is not None:
#             print(f"  P_dot/P = {target:.2e} yr^-1: q = {q_val:.3f}")
#         else:
#             print(f"  P_dot/P = {target:.2e} yr^-1: 無解 (q 超出範圍或無交點)")

# --------------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt

# # L_x
# L_x = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.82]

# # 週期變化率中間的 q (對應 \dot{P}/P = 1.42)
# q_middle = [0.061, 0.072, 0.083, 0.093, 0.104, 0.114, 0.124, 0.134, 0.143, 0.153, 0.162, 0.171, 0.180, 0.188, 0.190]

# # 週期變化率上限的 q (對應 \dot{P}/P = 0.41)
# q_upper = [0.182, 0.211, 0.238, 0.263, 0.286, 0.309, 0.329, 0.349, 0.367, 0.384, 0.401, 0.417, 0.431, 0.445, 0.448]

# # 週期變化率下限的 q (對應 \dot{P}/P = 2.43)
# q_lower = [0.036, 0.043, 0.050, 0.057, 0.063, 0.070, 0.076, 0.083, 0.089, 0.095, 0.101, 0.107, 0.113, 0.119, 0.120]

# # 轉換為「相對中心值的誤差長度」
# yerr_lower = [q_middle[i] - q_lower[i] for i in range(len(q_middle))]
# yerr_upper = [q_upper[i] - q_middle[i] for i in range(len(q_middle))]
# yerr = [yerr_lower, yerr_upper]

# for i in range(len(L_x)):
#     if i == len(L_x) - 1:  # 判斷是否為最後一個點
#         plt.errorbar(L_x[i], q_middle[i],
#                      yerr=[[yerr_lower[i]], [yerr_upper[i]]],
#                      fmt='*', capsize=5, color='red')  # 紅色星點
#     else:
#         plt.errorbar(L_x[i], q_middle[i],
#                      yerr=[[yerr_lower[i]], [yerr_upper[i]]],
#                      fmt='o', capsize=5)  # 其他點為藍色圓點
        
# plt.axhline(y=0.8/1.4, color='black', linestyle='-.', alpha=0.5)
# plt.text(0.05, 0.74/1.4, r'$\frac{0.8}{1.4}$', color='black')
# plt.axhline(y=0.14/1.4, color='black', linestyle='-.', alpha=0.5)
# plt.text(0.05, 0.18/1.4, r'$\frac{0.14}{1.4}$', color='black')
# plt.xlabel(r'$L_x \, (10^{38} \, \mathrm{erg/s})$')
# plt.ylabel('mass ratio  $q$')
# plt.xlim(0, 2)
# plt.ylim(0, 0.7)
# plt.legend(loc='upper right')
# plt.title(r'Mass-conserved case in Different $L_{acc}$')
# plt.grid(True)
# plt.tight_layout()
# plt.show()

# --------------------------------------------------------

# # 改變週期變化率固定質量比
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation

# # 常數和初始設置
# q = 0.14 / 1.4
# L_int_values = [5e30, 1e31, 1.5e31, 1.82e31]
# L_int_labels = [r'$L_{int}=5 \times 10^{37} \, \text{erg/s}$', 
#               r'$L_{int}=10 \times 10^{37} \, \text{erg/s}$', 
#               r'$L_{int}=15 \times 10^{37} \, \text{erg/s}$', 
#               r'$L_{int}=18.2 \times 10^{37} \, \text{erg/s (Eddington)}$']
# beta_values = np.linspace(0.0001, 0.9999, 500)
# colors = ['blue', 'green', 'purple', 'red']

# # 設置圖表
# fig, ax = plt.subplots(figsize=(10, 8))

# # 初始化四條曲線
# lines = [ax.plot([], [], label=label, color=color)[0] for label, color in zip(L_int_labels, colors)]

# # 設置軸範圍和標籤
# ax.set_xlim(0, 1)
# ax.set_ylim(0, 1)
# ax.set_xlabel(r'$\alpha$')
# ax.set_ylabel(r'$1 - \beta$')
# ax.set_title('Mass non-conserved case')
# ax.grid(True)
# ax.legend()

# # 初始化參數框的文本
# parameters_text = ax.text(0.7, 0.15, '', transform=ax.transAxes, fontsize=12, 
#                          verticalalignment='center', horizontalalignment='left',
#                          bbox=dict(facecolor='white', edgecolor='k', alpha=1.0))

# # 定義 P_dot_over_P 的範圍
# P_dot_over_P_values = np.linspace(1e-8, 3e-7, 100)  # 從 1e-8 到 3e-7，100 幀

# # 更新函數（每一幀調用）
# def update(frame):
#     P_dot_over_P = P_dot_over_P_values[frame]
#     constant_term = ((1/3) * P_dot_over_P + 1.26e-11 * q * (1 + q)**(-1/3)) * 1.66e39 * q
    
#     # 更新每條曲線
#     for line, L_int in zip(lines, L_int_values):
#         alpha_values = (1 + q) * (1 - q * beta_values - constant_term * beta_values / L_int) / (1 - beta_values) - q/3
#         one_minus_beta_values = 1 - beta_values
#         line.set_data(alpha_values, one_minus_beta_values)
    
#     # 更新參數框
#     parameters = (
#         r'$M_1 = 1.4 \, M_\odot$' + '\n' +
#         r'$R_1 = 10 \, \text{km}$' + '\n' +
#         r'$q = \frac{0.14}{1.4}$')
#     # + '\n' +
#     #    r'$\frac{\dot{P}_{\text{orb}}}{P_{\text{orb}}} \, range \, from \, 10^{-8} \, to \, 3 \times 10^{-7} \text{yr}^{-1}$'
#     parameters_text.set_text(parameters)
    
#     return lines + [parameters_text]

# # 創建動畫
# ani = FuncAnimation(fig, update, frames=len(P_dot_over_P_values), interval=50, blit=True)

# # 保存動畫（需要安裝 ffmpeg）
# ani.save('mass_non_conserved_animation(0.14).mp4', writer='ffmpeg', fps=20)

# # 顯示動畫（可選）
# plt.show()