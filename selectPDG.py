#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:33:30 2026

@author: fedora
"""

import os
import sys
import numpy as np
import pandas as pd

from funcs import readConfig, getGInxLC, get_peaks

if len(sys.argv)>1:
    config_name = sys.argv[1]
else:
    config_name = '../configs/configs_search_V7F5/selectPDG.json'
#############################################################
# Preparation
#############################################################
#### read configs
(dir_res, field, ver_LC, ver_search,
 ver_cor,
 mode, 
 d_rel, N_peak, N_per, mean,
 v_max, n_peaks, 
run) = readConfig(config_name)
if N_per==0:
    N_per = N_peak
    print('N_per = N_per')
#### create folders
dir_search = dir_res + 'search/search' + ver_search + '/'
dir_conf = dir_search + 'configs/'
dir_BLS = dir_BLS = dir_search + 'BLS/' 
dir_BLS_cor = dir_BLS[:-1] + ver_cor + '/'
dir_res = dir_search + 'results/'
os.makedirs(dir_res, exist_ok=True)

#### read Period grid, G_inx from BLS_cor folder = G_run
Per = np.load(dir_conf + 'P_BLS' + ver_search + ver_cor + '.npy')
G_run, _ = getGInxLC(dir_BLS_cor, field, ver_LC)

print(f'len G_run: {len(G_run)}')
print(f'mean: {mean}')
#############################################################
# Run
#############################################################
### !!!! pd -> np in cycle
if run:
    print('\n SELECT PDG')
    S = pd.DataFrame(columns=['g_inx', 'per', 'val', 'i_st'])
    S.loc[-1] = [-1, -1, -1, -1]
    for ij, g_inx in enumerate(G_run):
        file = dir_BLS_cor + 'BLS_' + str(g_inx) + '.npy'
        PWR = np.load(file, allow_pickle=True)
        
        if mean:
            PWR = [np.median(PWR, axis=0)]

        for i_st in range(len(PWR)):
            pwr = PWR[i_st] #['pwr']
            pwr = pwr / np.max(pwr)
            n_max = get_peaks(pwr, N_peak)
            p_wins = Per[n_max] # np.round(Per[n_max], 3)
            val_wins = pwr[n_max] # np.round(pwr[n_max], 1)
            #p_wins = np.array([13.45, 93.90, 80.45])
            if mode=='Drel':
                if d_rel > 0:
                    for i in range(N_per):
                        rel1 = p_wins / p_wins[i]
                        rel2 = 1 / rel1
                        u1 = (rel1 >= 2) & ((rel1%1) < d_rel)
                        u2 = (rel2 >= 2) & ((rel2%1) < d_rel)
                        u3 = (sum(val_wins >= v_max)  <= n_peaks ) & (val_wins[i] >= v_max)
                        u = (u1 | u2) & u3
                        if sum(u)>0:
                            S.loc[len(S)] = [g_inx, p_wins[i], val_wins[i], i_st]
                       #     for j in np.argwhere(u==True).flatten():
                       #         S.loc[len(S)] = [g_inx, p_wins[j], val_wins[j], i_st]
                else:
                    S = pd.concat((S, pd.DataFrame({'g_inx':g_inx, 
                                                    'per':p_wins[:N_per], 
                                                    'val':val_wins[:N_per],
                                                    'i_st':i_st
                                                    })
                                   ))
            elif mode=='Peaks':
                u = val_wins >= v_max
                if sum(u) <= n_peaks:
                    for j in np.argwhere(u==True).flatten():
                        S.loc[len(S)] = [g_inx, p_wins[j], val_wins[j], i_st]           
        if ij%100==0:
            print(ij, end=', ')

    ver_sel = f'_{d_rel}-{N_peak}-{N_per}_{v_max}-{n_peaks}'

    
    S = S.drop(-1, axis=0)
    S['g_inx'] = S['g_inx'].astype(int)
    S = S.drop_duplicates()
    
    name = 'Periods_BLS' + ver_search + ver_cor + ver_sel + '.csv'
    S.to_csv(dir_res + name)
    print(f'\n save {name} table to {dir_res}')
    print(S)
else:
    print('\n\nrun: false\n\n')










































