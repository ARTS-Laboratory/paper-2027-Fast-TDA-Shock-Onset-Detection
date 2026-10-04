#%% import necessary modules
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

#%% setting up font for figures
plt.rcParams.update({'image.cmap': 'viridis'})
cc = plt.rcParams['axes.prop_cycle'].by_key()['color']
plt.rcParams.update({'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif',
                                    'Bitstream Vera Serif', 'Computer Modern Roman', 'New Century Schoolbook',
                                    'Century Schoolbook L',  'Utopia', 'ITC Bookman', 'Bookman',
                                    'Nimbus Roman No9 L', 'Palatino', 'Charter', 'serif']})
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'font.size': 10})
plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.rm': 'serif'})
plt.rcParams.update({'mathtext.it': 'serif:italic'})
plt.rcParams.update({'mathtext.bf': 'serif:bold'})
plt.close('all')

#%% uploatd dataset
xlsx_path = r"signal_parameters.xlsx"

df = pd.read_excel(xlsx_path)

signal=[]
'''

Each number corresponds to the respective signal
1- AMP_3-VIB_0.1-SH_0.5
2- AMP_8-VIB_0.1_SH_1
3- AMP_10-VIB_0.1-SH_0.5
4- AMP_10-VIB_0.4-SH_0.4

'''
a=[]
b=[]
c=[]
d=[]
e=[]
for i in range(1,5):
    signal.append(int(df.iloc[i,0]))
    a.append(df.iloc[i,1])
    b.append(df.iloc[i,2])
    c.append(df.iloc[i,3])
    d.append(df.iloc[i,4])
    e.append(df.iloc[i,5])

#%% plot barcharts

# Bar width and x locations
w = .15
x = np.arange(len(signal))

fig, ax = plt.subplots()
ax.bar(x-2*w, a, width=w, label='a',color=(0,0.45,0.72)) # a parameter
ax.bar(x-w, b, width=w, label='b',color=(0.85,0.33,0.10)) # b parameter
ax.bar(x,c,width=w,label="c",color=(0.93,0.69,0.13)) # c parameter
ax.bar(x+w,d,width=w,label="d",color=(0.49,0.18,0.56)) # d parameters
ax.bar(x+2*w,e,width=w,label="e",color=(0.47,0.67,0.19)) # e parameter

ax.set_xticks(x)
ax.set_xticklabels(signal)
ax.set_ylabel('latency (s)') 
ax.legend()

fig.set_figwidth(6.5)   # width in inches
fig.set_figheight(2.5)  # height in inches
plt.ylim([-.75,.125])   # set y-limits
plt.axhline(0,color='black')

plt.tight_layout()
plt.savefig('fast_tda_latency_bar_chart.png') # saves .png file
plt.savefig('fast_tda_latency_bar_chart.svg') # saves .svg file
