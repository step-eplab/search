#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 12:40:04 2024

@author: fedora
"""

import os
import json

import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle 
from astropy.timeseries import BoxLeastSquares 

from funcs import readLC, run_multi

###############################################################################
def p_grid(P_ranges):
    P = []
    for r in P_ranges:
        P = np.concatenate((P, np.arange(r[0], r[1], r[2])))
    return P

#############################################################

def LS(t, f, Per):
    model = LombScargle(t, f)
    pwr = model.power(1/Per)
    return pwr

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
    
#############################################################

def get_cols_main(SLC, JD_ret=True):
    cols_mag = pd.Series(SLC.columns[1:])
    cols_S = cols_mag[cols_mag.str[-1]=='S'].values
    cols_F = cols_mag[cols_mag.str[-1]=='F'].values
    col_med = cols_mag[cols_mag.str[-3:]=='med'].values
    if JD_ret:
        cols = np.concatenate((['JD'], cols_S, cols_F, col_med))
    else:
        cols = np.concatenate((cols_S, cols_F, col_med))
    return cols

def calcPDG(LC, Pers_LS, Pers_BLS, Durs_BLS, min_dt, min_N, reg='both', cols_reg='short'):
    # reg = BLS / LS / both
    LS_res_g = {}
    BLS_res_g = {}

    cols_mag = LC.columns[1:]
    cols_main_0 = get_cols_main(LC, JD_ret=False)
    cols_ind = cols_mag[~np.isin(cols_mag, cols_main_0)]

    if cols_reg=='long':
        if len(cols_main_0) < 5:
            '''
            n = np.min([len(cols_ind), 6-len(cols_main_0)])
            cols_add = np.random.choice(cols_oth, n, replace=False)
            '''
            n = 5 - len(cols_main_0)
            L = len(LC) - LC[cols_ind].isna().sum(axis=0)
            inx = np.argsort(L)[::-1]
            cols_add = cols_ind[inx[:n]]
            cols_main = np.append(cols_main_0, cols_add)

    if cols_reg=='short':
        cols_main = np.array(['med', 'med_S', 'med_F'])
        cols_main = cols_main[np.isin(cols_main, cols_mag)]

        if len(cols_main) < 3:
            n = 3 - len(cols_main)                                                                                                                                                                                                            
            L = len(LC) - LC[cols_ind].isna().sum(axis=0)                                                                                                                                                                                     
            inx = np.argsort(L)[::-1]                                                                                                                                                                                                         
            cols_add = cols_ind[inx[:n]]                                                                                                                                                                                                      
            cols_main = np.concatenate((cols_main, cols_add))

    cols_main = cols_main[np.isin(cols_main, cols_mag)]
    LC['t'] = (LC.JD - LC.JD.min()) * 24
    for c in cols_main:
        TF = LC.loc[~LC[c].isna(), ['t', c]].values
        t, f = TF[:, 0], TF[:, 1]

        dt = max(t) - min(t)
        dur = len(t)
        if (dt > min_dt) & (dur > min_N):
            if reg!='BLS':
                pwr_ls = LS(t, f, Pers_LS)
                LS_res_g[c] = pwr_ls

            if reg!='LS':
                pwr_bls = BLS(t, f, Pers_BLS, Durs_BLS, stat=0)
                BLS_res_g[c] = pwr_bls

    return LS_res_g, BLS_res_g

#############################################################

def createPDG(g_inx, Pers_LS, Pers_BLS, Durs_BLS,
                dir_SLC, dir_search, min_dt, min_N, reg):
    LC = readLC(dir_SLC, 'SLC', g_inx)
    dir_BLS = dir_search + 'BLS/'
    dir_LS = dir_search + 'LS/'
    if not LC is None:
        LS_res, BLS_res = calcPDG(LC, Pers_LS, Pers_BLS, Durs_BLS, min_dt, min_N, reg)
        
        if len(LS_res)>0:
            np.save(dir_LS + 'LS_' + str(g_inx) + '.npy', LS_res)
        if len(BLS_res)>0:
            np.save(dir_BLS + 'BLS_' + str(g_inx) + '.npy', BLS_res)

###############################################################################
with open('configs.json', 'r') as file:
    configs = json.load(file)

run_proc = configs['run_proc']
multi_proc = configs['multi_proc']

dir_main = configs['dir_main']
cat_path = configs['cat_path']
ver_slc = configs['ver_slc']
ver_search = configs['ver_search']

min_dt = configs['min_dt']
min_N = configs['min_N']
mag_max = configs['mag_max']

N_g_inx = configs['N_g_inx']

Pers_LS = p_grid(configs['Pers_LS'])
Pers_BLS = p_grid(configs['Pers_BLS'])
Durs_BLS = p_grid([configs['Durs_BLS']])

reg = configs['reg']

###############################################################################

dir_SLC0 = dir_main + 'SLC/'
dir_SLC = dir_SLC0 + 'SLC' + ver_slc + '/'

dir_search0 = dir_main + 'search/'
dir_search = dir_search0 + 'search' + ver_search + '/'
dir_BLS = dir_search + 'BLS/'
dir_LS = dir_search + 'LS/'
dir_conf = dir_search + 'configs/'

for d in [dir_search0, dir_search, dir_BLS, dir_LS,  dir_conf]:
    os.makedirs(d, exist_ok=True)

for D, name in zip([Pers_LS, Pers_BLS, Durs_BLS], ['P_LS', 'P_BLS', 'D_BLS']):
    file = dir_conf + name + ver_search + '.npy'
    if not os.path.isfile(file):
        np.save(file, D)
        print('save ' + file)
    else:
        print(file + ' exists')


G_del = configs['G_del']
G_select = configs['G_select']
###############################################################################
# LC file format: '*_{g_inx}.*'
G_SLC = pd.Series(os.listdir(dir_SLC)).str.split('_').str[1].str.split('.').str[0].astype(int).values
if cat_path:
    F = pd.read_csv(cat_path, sep=' ', names=['ra', 'dec', 'mag', 'source_id'])
    F_inx = F.loc[F.mag<mag_max].index

if len(G_select)==0:
    G_select = G_SLC

G_run = G_select[~np.isin(G_select, G_del) & np.isin(G_select, F_inx)]

for G, txt in zip([G_del, G_select, G_run], ['del', 'select', 'run']):
    print(f'len G_inx {txt} : ' + str(len(G)))

if N_g_inx>0:
    G_run = G_run[:N_g_inx]

if run_proc:
    if multi_proc:
        run_multi(createPDG, G_run, Pers_LS, Pers_BLS, Durs_BLS, dir_SLC, dir_search, min_dt, min_N, reg)
    else:
        print('createPDG: ', ', ')
        for g_inx in G_run:
            print(g_inx, end=', ')
            createPDG(g_inx, Pers_LS, Pers_BLS, Durs_BLS, dir_SLC, dir_search, min_dt, min_N, reg)

#######################################
L1 = len(G_run)
L2 = [len(os.listdir(dir_LS)), len(os.listdir(dir_BLS))]
X = ['len of objects (input)', 'len of objects (output: LS, BLS)',
    'path to LC for search', 'directory for result',
    'minum lenght of LC for analysis (hours)']
Y = [L1, L2, dir_SLC, dir_search, min_dt, min_N]
S = ''
for x, y in zip(X, Y):
    s = x + ' : ' + str(y)
    print(s)
    S += s + '\n'

with open(dir_conf + 'search' + ver_search + '_configs', 'w') as x:
    x.write(S)


###############################################################################

'''
import matplotlib.pyplot as plt

res1 = np.load('/home/fedora/astronomy/STEP/main/copy/F1/v6/search/search_4/BLS/BLS_18986.npy', allow_pickle=True).tolist()
plt.figure()
for c in res1.keys():
    pwr = res1[c]
    plt.plot(Pers_BLS, pwr,'.-', alpha=0.2)
'''


'''
G_var = np.array([0,  2331,  2875,  3737,  3883,  5466,  7199, 9285, 9651, 10780,
       11122, 11372, 14694])

for g_inx in G_var:
    LCS[g_inx] = CLCS
    if not CLCS is None:
        plotCLC(CLCS)
        plt.title(g_inx)


for g_inx in LS_res.keys():
    LS_rg = LS_res[g_inx]
    plt.figure()
    for c in LS_rg.keys():
        pwr = LS_rg[c]
        plt.plot(Pers_LS, pwr, '.-')
    plt.title(g_inx)

'''


'''
LCS = {}
for g_inx in G_var:
    SLC = readLC(dir_SLC, 'SLC', g_inx)
    if not SLC is None:
        LCS[g_inx] = SLC
        
for g_inx in LS_res.keys():
    LS_rg = LS_res[g_inx]
    ps = []
    for c in LS_rg.keys():
        pwr = LS_rg[c]
        p_win = Pers_LS[np.argmax(pwr)]
        ps.append(p_win)
    pv = pd.Series(ps).value_counts()
    p_win = pv.index[pv.argmax()]

    LC = LCS[g_inx]
    JD = LC.JD
    f = LC['med']
    t = (JD - min(JD)) * 24
    t_f, n_f = fold(t, p_win)
    plt.figure()
    plt.plot(t_f, f, '.')


'''

'''
ns = get_peaks(pwr, 1)
for i_max in ns:
    per = P[i_max]
    t_f, n_f = fold(t, per * 2)
    plt.figure()
    plt.plot(t_f, f, '.')
    plt.title(per)

'''