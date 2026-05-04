#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 13:22:30 2024

@author: fedora
"""

import os
import time
#import lzma

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import Process, Manager

from scipy.signal import find_peaks
###############################################################################

def calcRo(g_inx, Cat_g, Cat_ro):
    d_ra = Cat_ro.ra - Cat_g.loc[g_inx, 'ra']
    d_dec = Cat_ro.dec - Cat_g.loc[g_inx, 'dec']
    ro = (d_ra**2 + d_dec**2)**0.5
    if g_inx in ro.index:
        ro.loc[g_inx] = np.nanmax(ro)
    return ro

def trimSC(m0, k, n=1):
    m = m0.copy()
    for i in range(n):
        med = np.nanmedian(m)
        std = np.nanstd(m)
        u = abs(m - med) > (std * k)
        m[u] = np.nan
    return m


def run_multi(func, list_arg, args, dicts=[]):
    zero = time.time()
    n = os.cpu_count()
    N = len(list_arg)//n + 1
    print('len(Procs) = ' + str(N), '\n')
    
    if len(dicts)>0:
        manager = Manager()
        D_use = []
        for d in dicts:
            D = manager.dict()
            D_use.append(D)
 
    for i in range(N):
        list_arg_i = list_arg[i * n : (i + 1) * n]
        print(i, end=',')
        processes = []        
        for arg_i in list_arg_i:
            y = args[::-1]
            y.append(arg_i)
            if len(dicts)>0:
                for D in D_use[::-1]:                
                    y.append(D)
            args_i = y[::-1]
            proc = Process(target=func, args=(args_i))
            processes.append(proc)
            proc.start()
        for proc in processes:
            proc.join()
    if len(dicts)>0:
        for i in range(len(dicts)): 
            dicts[i].update(D_use[i])
        return dicts
    print(time.time() - zero)
    
#######################################

def readLC(dir_LC, pref, g_inx, date=''):
    post = ''
    if len(date)==8:
        post += '_' + date
    file = dir_LC + pref + '_' + str(g_inx) + post + '.csv'

    if os.path.isfile(file):
        LC = pd.read_csv(file)
        return LC
    else:
        '''
        file = dir_LC + pref + '_' + str(g_inx) + post + '.csv.xz'
        if os.path.isfile(file):
            file_open = lzma.open(file)
            LC = pd.read_csv(file_open)
            return LC
        else:
        '''
        print(pref + '_' + str(g_inx) + ' doesnt exist')
#######################################

def plotRLC(RLC, ax=plt, parts=[], new_fig=True, drop_med=False, dm=0):
    if (ax==plt) & new_fig:
        plt.figure()  
    
    if len(parts)==0:
        parts = RLC.columns[1:]
    for p in parts:
        m = RLC[p].copy()
        if drop_med:
            m = m - np.nanmedian(m) + dm
        ax.plot(RLC.JD, m, 'o', alpha=0.5)

def plotCLC(CLC, x_JD=True, ax=plt, new_fig=True, dm=0, dm_0=0, drop_med=False):
    if x_JD:
        x = CLC.JD
    else:
        x = (CLC.JD - min(CLC.JD))*24  
    if (ax==plt) & new_fig:
        plt.figure()  
    for i, c in enumerate(CLC.columns[1:]):
        m = CLC[c]
        if drop_med:
            m = m - np.nanmedian(m)
        ax.plot(x, m + dm_0 + dm * i, '.')

        
#############################################################
    
def clearP(Per, p_0, dp, i0=1):
    i1 = int(max(Per) / p_0) + 2
    for i in range(i0, i1):
        u = (Per > (p_0 * i - dp)) & (Per < (p_0 * i + dp))
        Per = Per[~u]
    return Per
        
def get_peaks(y, n=1):
    n_max, peaks = find_peaks(y, height=0)
    ns = n_max[np.argsort(peaks['peak_heights'])[::-1]][:n]
    return ns


def dropna(x, y):
    u_nan = np.isnan(y)
    x = x[~u_nan]
    y = y[~u_nan]
    return x, y

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
    
def fold(t, per, fi_0=0):
    t_f = (t%per / per + fi_0)%1
    n = np.int32(t//per) + 1
    return t_f, n

def prep(JD, m, JD0=0, jdm=0):
    u = ~np.isnan(m)
    JD1 = JD[u]
    if JD0==0:
        JD0 = np.nanmin(JD1)

    y = 10**(-m[u]/2.5)
    t = (JD1 - JD0) * 24
    inx = np.argsort(t)
    if jdm:
        return t[inx], y[inx], JD[u][inx], m[u][inx]
    else:
        return t[inx], y[inx]
    
def clearJD(JD, JD_del):
    if len(JD_del) > 0:
        for i in range(len(JD_del)):
            jd = JD_del[i]
            if type(jd)!=list:
                JD_del[i] = [jd, jd+1]
                
        u_del = np.zeros(len(JD), dtype=bool)
        for jd in JD_del:
            u_del |= (JD > jd[0]) & (JD < jd[1])
        u_cl = ~u_del
        return u_cl
    else:
        return True

def normJD(JD, m0):
    m = m0.copy()
    JDf = np.floor(JD)
    JDfu = np.unique(JDf[~np.isnan(m)])
    for jd in JDfu:
        u = jd==JDf
        if sum(u)>0:
            m[u] = m[u] - np.nanmedian(m[u])
    return m

def detrendJD(JD, t, f, deg=1, min_points=30):
    JDf = np.floor(JD)
    JDfu = np.unique(JDf)
    Y = np.zeros(len(f)) * np.nan
    for jd in JDfu:
        u = jd==JDf
        y = f[u]
        if sum(u)>min_points:
            f0 = trimSC(f[u], 1, 3)
            u1 = ~np.isnan(f0)
            
            if sum(u1) > (min_points/2):
                p = np.polyfit(t[u][u1], f0[u1], deg)               
                fit = np.polyval(p, t[u])
                y_dtr = f[u] - fit
                if np.nanstd(y_dtr) < np.nanstd(f[u]):
                    y = y_dtr                    
        Y[u] = y - np.nanmedian(y) + 1
    return Y


###### Search

def splitG(d):
    files = os.listdir(d)
    G_inx = pd.Series(files).str.split('_').str[1].str.split('.').str[0]
    G_inx = np.unique(np.int32(G_inx))
    return G_inx

def readRes(g_inx, dir_BLS):
    file = dir_BLS + 'BLS_' + str(g_inx) + '.npy'
     
    if os.path.isfile(file):
        res = np.load(file, allow_pickle=True).tolist()
        return res
    else:
        '''
        file = file + '.xz'
        if os.path.isfile(file):
            file_open = lzma.open(file)
            try:
                res = np.load(file_open, allow_pickle=True).tolist()
                return res
            except EOFError:
                print('EOFError ' + str(g_inx), end=', ')
        else:
        '''
        print('de ' + str(g_inx), end=', ')


