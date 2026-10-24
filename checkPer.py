#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 17:17:41 2026

@author: fedora
"""

import matplotlib.pyplot as plt
import numpy as np


def sigma_clip(y, lvl=3, n=1):
    med = np.nanmedian(y)
    std = np.nanstd(y)
    u = np.zeros(len(y), dtype=bool)
    for i in range(n):
        u |= abs(y - med) > lvl * std
    y[u] = np.nan
    return y
def fold(t, per, fi_0=0):
    t_f = (t%per / per + fi_0)%1
    n = np.int32(t//per) + 1
    return t_f, n


dir_SLC = '/home/fedora/astronomy/STEP/main/results/lightcurves/V7/F5/SLC/SLC_F_st0_sc31/'
#dir_SLC = '/home/fedora/astronomy/STEP/main/results/lightcurves/V7/F4/SLC/SLC_JDnorm_sc51/'

g_inx = 12347
SLC = readLC(dir_SLC, g_inx, field, ver_LC)

p_win = 351
f = SLC['f']
t = SLC['t']
t_f, n_f = fold(t, p_win)
plt.figure()
plt.plot(t_f, f, '.')


plt.figure()
for n in np.unique(n_f):
    u = n_f==n
    plt.plot(t_f[u], f[u], '.')


plt.figure()
plt.plot(SLC.JD, f+0.15, '.')




for p_win in np.arange(278.5, 279.5, 0.1):
    t_f, n_f = fold(t, p_win)
    plt.figure()
    plt.plot(t_f, f, '.')
    plt.title(p_win)

f = sigma_clip(f, 3, 1)
fi0 = 0.136
dfi = 0.01
u = (t_f > (fi0 - dfi)) & (t_f < (fi0 + dfi))



plt.figure()
plt.plot(t_f, f, '.')
plt.plot(t_f[u], f[u], '.')


plt.figure()
plt.plot(t, f, '.')
plt.plot(t[u], f[u], '.')



plt.plot(t, f+0.1, '.')


x0 = 11055
x_min = t.min()
x_max = t.max()
y_min = f.min()
y_max = f.max()

per = p_win/2
c = 'g'
N = int((x_max - x0) / per) + 1
for i in range(N):
    x = x0 + i * per
    plt.plot([x, x], [y_min, y_max], color=c)

N = int((x0 - x_min) / per) + 1
for i in range(N):
    x = x0 - i * per
    plt.plot([x, x], [y_min, y_max], color=c)



x1 = 4996.22
per = x1 - x0


g_inx = 3520
SLC = readLC(dir_CLC, g_inx, field, ver_LC)
f = SLC['f']
t = SLC['t']



dir_CLC = '/home/fedora/astronomy/STEP/main/results/lightcurves/V7/F5/CLC/CLCF/'

g_inx = 12957

file = dir_CLC + 'CLC_' +  str(g_inx) + '.npy'
LC = np.load(file).T
JD = LC[6]

plt.figure()
for i in range(6):
    plt.plot(JD, 1 - LC[i] + i*1, '.')







