#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 14:50:08 2026

@author: fedora
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from funcs import readConfig, getGInxLC, savePNG

if len(sys.argv)>1:
    config_name = sys.argv[1]
else:
    config_name = 'configs/plotPDG.json'
#############################################################
# Preparation
#############################################################
#### read configs
dir_res, ver_search, ver_cor, G_del, G_select = readConfig(config_name)

#### create folder
dir_search0 = dir_res + 'search/'
dir_search = dir_search0 + 'search' + ver_search + '/'
dir_conf = dir_search + 'configs/'
dir_BLS = dir_search + 'BLS' + ver_cor + '/'
dir_BLSpng = dir_BLS[:-1] + 'png/'  
os.makedirs(dir_BLSpng, exist_ok=True)

files = os.listdir(dir_BLS)

G_inx, files = getGInxLC(dir_BLS)

G_run = G_inx[np.isin(G_inx, G_select)]

pers = np.load(dir_conf + 'P_BLS' + ver_search + ver_cor + '.npy')

### G_inx filter -> G_run = (G_select & G_inx & F_inx) - G_del
if len(G_select)==0:
    G_select = G_inx
else:
    G_select = np.array(G_select)
G_run = G_select[~np.isin(G_select, G_del)]

#############################################################
# Run
#############################################################
print('\n PLOT PDG: ')
for g_inx, file in zip(G_run, files):
    print(g_inx, end=', ')
    name = str(g_inx)
    plt.figure(figsize=(18,12))
    pwr = np.load(dir_BLS + file)
    plt.plot(pers, pwr, '.-')
    plt.grid()
    plt.title(name)
    plt.tight_layout()
    
    savePNG(dir_BLSpng, name)
























