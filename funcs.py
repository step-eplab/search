#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 10:53:43 2026

@author: fedora
"""

import os
import json
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from astropy.timeseries import BoxLeastSquares

import matplotlib.pyplot as plt

import matplotlib
matplotlib.use('Agg')  
print('QT5 off')
#matplotlib.use('qt5agg') 

def readConfig(config_name, keys=[]):
    print(f'Read: {config_name}')
    with open(config_name, 'r') as file:
        configs = json.load(file)

    if keys==[]:
        keys = configs.keys()
    print(f'Keys: {keys}')
    V = []
    for k in keys:
        V.append(configs[k])
    return V

def getGInxLC(dir_LC, field, ver_LC):
    if field=='F1':
        if ver_LC=='V6':
            files = pd.Series(os.listdir(dir_LC))
            G_inx = files.str.split('_').str[1].str.split('.').str[0]
            G_inx = G_inx.astype(int).values
    return G_inx

#############################################

def pGrid(P_ranges):
    P = []
    for r in P_ranges:
        P = np.concatenate((P, np.arange(r[0], r[1], r[2])))
    return P

def readLC(dir_LC, g_inx, field, ver_LC):
    if field=='F1':
        if ver_LC=='V6':
            file = dir_LC + 'SLC_' + str(g_inx) + '.csv'
            
            if os.path.isfile(file):
                LC = pd.read_csv(file)
                if 'med_S' in LC.columns:
                    LC = LC[['JD', 'med_S']]
                else:
                    cols = LC.columns[1:]
                    col = cols[LC[cols].isna().sum(axis=0).argmin()]
                    LC = LC[['JD', col]]
                LC['t'] = (LC.JD - LC.JD.min()) * 24
                LC.columns = ['JD', 'f', 't']
                LC = LC.dropna(subset=['f'])
                return LC

def BLS(t, f, Per, Dur, stat=0):
    model = BoxLeastSquares(t, f)
    res = model.power(Per, Dur)
    pwr = res.power
    if stat:
        inx = np.argmax(pwr)
        S = {'per': Per[inx], 
             'dep': res.depth[inx],
             'dur': res.duration[inx],
             't0': res.transit_time[inx],
             'pwr': res.power}
        return S
    else:
        return pwr

#############################################

def readRes(g_inx, dir_BLS):
    file = dir_BLS + 'BLS_' + str(g_inx) + '.npy'
    if os.path.isfile(file):
        res = np.load(file, allow_pickle=True)
        return res
    else:
        print('Not found: ' + str(g_inx), end=', ')

def fit_bin(x, y, n, sig_ret=0):
    X = np.linspace(min(x), max(x), n)
    X_av = X[:-1] + (X[1]-X[0])/2
    Y_av = np.zeros(len(X_av)) * np.nan

    for i in range(n-1):
        u = (x > X[i]) & (x < X[i+1])
        if np.sum(u)>0:
            if np.sum(~np.isnan(y[u]))==0:
                Y_av[i] = np.nan
            else:
                Y_av[i] = np.nanmedian(y[u])     
    return X_av, Y_av   

def clearP(Per, p_0, dp, i0=1):
    i1 = int(max(Per) / p_0) + 2
    for i in range(i0, i1):
        u = (Per > (p_0 * i - dp)) & (Per < (p_0 * i + dp))
        Per = Per[~u]
    return Per

def dropna(x, y):
    mask = ~np.isnan(y)
    return x[mask], y[mask]

def get_peaks(y, n=1):
    n_max, peaks = find_peaks(y, height=0)
    ns = n_max[np.argsort(peaks['peak_heights'])[::-1]][:n]
    return ns


#############################################
def plotFLC(SLC, pers, cols):
    t = (SLC.JD - min(SLC.JD))*24
    Nc = len(cols)          
    fig, Ax = plt.subplots(Nc, 1, figsize=(16, 10))
    if Nc==1:
        Ax = [Ax]
    for c, ax, p_win in zip(cols, Ax, pers):
        f = SLC[c]
        f = trimSC(f, 3, 1)
        t_f, n_f = fold(t, p_win)
        ax.plot(t_f, f, '.', alpha=0.5)
        for n, clr, lnw in zip([64, 164, 512], ['k', 'r', 'lime'], [3, 2, 1]):
            x_av, y_av = fit_bin(t_f, f, n)
            ax.plot(x_av, y_av, color=clr, linewidth=lnw)
            ax.grid()
        ax.set_title(c + ', per = ' + str(round(p_win, 4)) + 'h')


