#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 19:02:02 2024

@author: fedora
"""

import os
import sys
import numpy as np

from funcs import readConfig, getGInxLC, readRes, fit_bin, clearP, dropna

def procPDG(g_inx, dir_BLS, dir_BLS_cor, u_per):
    try:
        res = readRes(g_inx, dir_BLS)
    except:
        res = None
    if not res is None:
        pwr = res[u_per]
        x_av, y_av = fit_bin(Per, pwr, 24)
        x_av, y_av = dropna(x_av, y_av)
        p = np.polyfit(x_av, y_av, 1)
        y = np.polyval(p, Per)
        pwr2 = pwr/y - 1

        np.save(dir_BLS_cor + 'BLS_' + str(g_inx) + '.npy', pwr2)

if len(sys.argv)>1:
    config_name = sys.argv[1]
else:
    config_name = 'configs/procPDG.json'
#############################################################
# Preparation
#############################################################
#### read configs
(dir_res, field, ver_LC, ver_search,
 ver_cor, 
 P_min, P_max, P_del, P_del_range) = readConfig(config_name)

#### create folders
dir_search = dir_res + 'search/search' + ver_search + '/'
dir_conf = dir_search + 'configs/'
dir_BLS = dir_BLS = dir_search + 'BLS/' 
dir_BLS_cor = dir_BLS[:-1] + ver_cor + '/'
os.makedirs(dir_BLS_cor, exist_ok=True)

### read Period grid, G_inx from BLS folder = G_run
P_b = np.load(dir_conf + 'P_BLS' + ver_search + '.npy') 
G_run = getGInxLC(dir_BLS, field, ver_LC)

### Period grid correction 
### P_grid=(P_min<P<P_max)&(P!=(p in n*[p_del-p_del_range, p_del+p_del_range]))
Per = P_b.copy()
Per = Per[(Per > P_min) & (Per < P_max)]
for p_del, dp in zip(P_del, P_del_range):
    Per = clearP(Per, p_del, dp)
    print('del Per = ' + str(p_del) + '*n +- ' + str(dp))
u_per = np.isin(P_b, Per)

#### save configs
np.save(dir_conf + 'P_BLS' + ver_search + ver_cor + '.npy', Per)
print('len init. Per = ' + str(len(P_b)))
print('len corr. Per = ' + str(len(Per)))

#############################################################
# Run
#############################################################
print('\n PROC PDG')
for g_inx in G_run:
    print(g_inx, end=', ')
    procPDG(g_inx, dir_BLS, dir_BLS_cor, u_per)

print('\n\n\nlen init. PDG : ' + str(len(G_run)))
print('len corr. PDG : ' + str(len((os.listdir(dir_BLS_cor)))))
