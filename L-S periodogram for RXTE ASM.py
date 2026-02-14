import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

with open("C:/Users/jerry/Desktop/junleiwu d/4U1820-30/data analyze/data/RXTE ASM/filter/RXTEASM 1.5-12keV_light_curve filter.tsv", "r") as fobj:
    Time, Rate = np.loadtxt(fobj, usecols=(0, 1), unpack=True)

N = len(Time)
frequencies = np.linspace(0, 0.02, 1000)
hbar = np.average(Rate)
sig2 = (1 / (N - 1)) * np.sum((Rate - hbar)**2)
h = (Rate - hbar)

P = []

for freq in frequencies:
    omg = 2 * np.pi * freq
    tou = (1 / (2 * omg)) * np.arctan(np.sum(np.sin(2 * omg * Time)) / np.sum(np.cos(2 * omg * Time)))
    t = (Time - tou)
    cos = np.sum(h * np.cos(omg * t))
    sin = np.sum(h * np.sin(omg * t))
    cos2 = np.sum(np.cos(t)**2)
    sin2 = np.sum(np.sin(t)**2)
    power = (1 / (2 * sig2)) * ((cos**2 / cos2) + (sin**2 / sin2))
    P.append(power)

max_index = np.argmax(P)

# 获取最高点的x和y值
max_x = frequencies[max_index]
max_y = P[max_index]

# 繪製結果
p = 1 / max_x

print(max_x)
print(max_y)

print(p)


plt.figure(figsize=(10, 6))
fig, ax = plt.subplots()
plt.plot(frequencies, P, color='blue')
plt.xlabel('Frequency (1/d)')
plt.ylabel('Power')
plt.title('RXTE ASM 1.5-12 keV \n Lomb-Scargle Periodogram ')
plt.xticks(np.arange(0, 0.02, 0.002))
plt.annotate(f'Max Point\nFrequency: {maxx:.7f}\nPower: {maxy:.3f}\nPeriod: {p:.3f} days',
             xy=(maxx, maxy), xycoords='data',
             xytext=(50, -50), textcoords='offset points',
             arrowprops=dict(arrowstyle="->"))
plt.ylim(0, 8000)
formatter = ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-1, 1))
ax.xaxis.set_major_formatter(formatter)
plt.show()
# fig.savefig('RXTE ASM 1.5-12 keV Lomb-Scargle Periodogram with little peaks.eps')
# fig.savefig('RXTE ASM 1.5-12 keV Lomb-Scargle Periodogram with little peaks.pdf')


