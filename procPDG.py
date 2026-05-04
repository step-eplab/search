#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 19:02:02 2024

@author: fedora
"""

import os
import sys
import json
import numpy as np

from funcs import splitG, readRes, fit_bin, clearP, dropna, run_multi

######################################
def procPDG(g_inx, dir_BLS, dir_BLS_cor, u_per):
    try:
        res = readRes(g_inx, dir_BLS)
    except:
        res = None
    if not res is None:
        res2 = {}
        for c in res.keys():
            pwr = res[c][u_per]
            x_av, y_av = fit_bin(Per, pwr, 24)
            x_av, y_av = dropna(x_av, y_av)
            p = np.polyfit(x_av, y_av, 1)
            y = np.polyval(p, Per)
            pwr2 = pwr/y - 1
            res2[c] = pwr2
            
        np.save(dir_BLS_cor + 'BLS_' + str(g_inx) + '.npy', res2)


###############################################################################
if len(sys.argv)>1:
    config_name = sys.argv[1]
else:
    config_name = 'configs_default.json'
    

with open(config_name, 'r') as file:
    configs = json.load(file)

run_proc = configs['run_proc']
multi_proc = configs['multi_proc']

dir_res = configs['dir_res']

ver_search = configs['ver_search']
ver_cor = configs['ver_cor']

P_min = configs['P_min']
P_max = configs['P_max']
P_del = configs['P_del']
P_del_range = configs['P_del_range']


dir_search = dir_res + 'search/search' + ver_search + '/'
dir_conf = dir_search + 'configs/'
dir_BLS = dir_BLS = dir_search + 'BLS/' 
dir_BLS_cor = dir_BLS[:-1] + ver_cor + '/'
os.makedirs(dir_BLS_cor, exist_ok=True)

P_b = np.load(dir_conf + 'P_BLS' + ver_search + '.npy') 
G_inx = splitG(dir_BLS)

Per = P_b.copy()
Per = Per[(Per > P_min) & (Per < P_max)]
for p_del, dp in zip(P_del, P_del_range):
    Per = clearP(Per, p_del, dp)
    print('del Per = ' + str(p_del) + '*n +- ' + str(dp))

u_per = np.isin(P_b, Per)

np.save(dir_conf + 'P_BLS' + ver_search + ver_cor + '.npy', Per)
print('len init. Per = ' + str(len(P_b)))
print('len corr. Per = ' + str(len(Per)))
######################################
print('correct PDG')
G_run = G_inx

if run_proc:
    if multi_proc:
        run_multi(procPDG, G_run, [dir_BLS, dir_BLS_cor, u_per], dicts=[])
    else:
        print('procPDG: ', end='')
        for g_inx in G_run:
            print(g_inx, end=', ')
            procPDG(g_inx, dir_BLS, dir_BLS_cor, u_per)

print('\n\n\nlen init. PDG : ' + str(len(G_inx)))
print('len corr. PDG : ' + str(len((os.listdir(dir_BLS_cor)))))

###############################################################################

'''
import matplotlib.pyplot as plt
for g_inx in G_inx[:10]:
    res = readRes(g_inx, dir_BLS)
    res2 = readRes(g_inx, dir_BLS_cor)
    
    plt.figure()
    plt.plot(P_b, res['med'], alpha=0.5)
    plt.plot(Per, res2['med'], alpha=0.5)
'''
