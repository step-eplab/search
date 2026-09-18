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

from funcs import readConfig, getGInxLC, readRes, get_peaks

if len(sys.argv)>1:
    config_name = sys.argv[1]
else:
    config_name = 'configs/selectPDG.json'
#############################################################
# Preparation
#############################################################
#### read configs
(dir_res, field, ver_LC, ver_search,
 ver_cor, 
 d_rel, N_peak, N_per) = readConfig(config_name)

#### create folders
dir_search = dir_res + 'search/search' + ver_search + '/'
dir_conf = dir_search + 'configs/'
dir_BLS = dir_BLS = dir_search + 'BLS/' 
dir_BLS_cor = dir_BLS[:-1] + ver_cor + '/'
dir_res = dir_search + 'results/'
os.makedirs(dir_res, exist_ok=True)

#### read Period grid, G_inx from BLS_cor folder = G_run
Per = np.load(dir_conf + 'P_BLS' + ver_search + ver_cor + '.npy')
G_run = getGInxLC(dir_BLS_cor, field, ver_LC)

#############################################################
# Run
#############################################################
S = pd.DataFrame(columns=['g_inx', 'per', 'val'])
for ij, g_inx in enumerate(G_run):
    pwr = readRes(g_inx, dir_BLS_cor)

    pwr_s = pwr / np.nanstd(pwr)
    n_max = get_peaks(pwr_s, N_peak)
    
    p_wins = Per[n_max]
    val_wins = pwr_s[n_max]
    for i in range(N_per):
        rel1 = p_wins / p_wins[i]
        rel2 = 1 / rel1
        u1 = (rel1 >= 2) & (rel1%1 < d_rel)
        u2 = (rel2 >= 2) & (rel2%1 < d_rel)
        r = sum(u1) + sum(u2)
        if r > 0:
            S.loc[len(S)] = [g_inx, round(p_wins[i], 3), round(val_wins[i], 1)]
    if ij%100==0:
        print(ij, end=', ')

S.to_csv(dir_res + 'res_BLS' + ver_search + ver_cor + '.csv')
print(len(S))
