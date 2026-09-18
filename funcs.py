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
from scipy.stats import binned_statistic

from astropy.timeseries import BoxLeastSquares

import matplotlib.pyplot as plt



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
    X = np.linspace(np.min(x), np.max(x), n)

    Y_av, _, _ = binned_statistic(
        x,
        y,
        statistic='median',
        bins=X
    )

    X_av = (X[:-1] + X[1:]) / 2

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


def trimSC(m0, k, n=1):
    m = m0.copy()
    for i in range(n):
        med = np.nanmedian(m)
        std = np.nanstd(m)
        u = abs(m - med) > (std * k)
        m[u] = np.nan
    return m

    
def fold(t, per, fi_0=0):
    t_f = (t%per / per + fi_0)%1
    n = np.int32(t//per) + 1
    return t_f, n


#############################################

def saveFLC(dir_save, g_inx, n_plot):
    plt.suptitle(str(g_inx))
    plt.tight_layout()
    plt.savefig(dir_save + '/' + str(g_inx) + '_' + str(n_plot) + '.png')
    plt.close()
    
def plotFLC(dir_save, g_inx, SLC, pers, N_plot=3):
    t = SLC['t']
    N_per = len(pers)
    N_fig = (N_per - 1) // N_plot + 1
    
    N_plot2 = N_per % N_plot 
    N_plots = np.full(N_fig, N_plot)
    if N_plot2>0:
        N_plots[-1] = N_plot2
    
    f = SLC['f']
    std = np.nanstd(f)
    ylim = [(1 - 5*std), (1 + std)]
    i_per = 0
    for i, n_plot in enumerate(N_plots):
        fig, Ax = plt.subplots(n_plot, 1, figsize=(16, 10))
        if n_plot==1:
            Ax = [Ax]
        
        for ax in Ax:
            p_win = pers[i_per]
            #f = trimSC(f, 3, 1)
            t_f, n_f = fold(t, p_win)

            ax.plot(t_f, f, '.', alpha=0.5)
            print(f'plot p = {p_win}')
            for n, clr, lnw in zip([64, 164, 512], ['k', 'r', 'lime'], [3, 2, 1]):
                x_av, y_av = fit_bin(t_f, f, n)
                ax.plot(x_av, y_av, color=clr, linewidth=lnw)
                ax.grid()
                ax.set_ylim(ylim)

            ax.set_title('per = ' + str(round(p_win, 4)) + 'h')
            
            i_per += 1
        saveFLC(dir_save, g_inx, i)






