#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 22:54:54 2026

@author: fedora
"""

import os
import sys
import numpy as np
import pandas as pd

from funcs import readConfig, getGInxLC, pGrid, readLC, BLS
    
def calcPDG(LC, Pers_BLS, Durs_BLS):
    TF = LC[['t', 'f']].values
    t, f = TF[:, 0], TF[:, 1]

    pwr_bls = BLS(t, f, Pers_BLS, Durs_BLS, stat=0)
    return pwr_bls


def createPDG(g_inx, dir_LC, field, ver_LC, dir_search,
              Pers_BLS, Durs_BLS, 
              min_dt, min_N):
    LC = readLC(dir_LC, g_inx, field, ver_LC)
    if not LC is None:
        t = LC['t']
        if ((max(t) - min(t)) > min_dt) & (len(t) > min_N):
            dir_BLS = dir_search + 'BLS/'
            BLS_res = calcPDG(LC, Pers_BLS, Durs_BLS) #, min_dt, min_N)
            if not BLS_res is None:
                np.save(dir_BLS + 'BLS_' + str(g_inx) + '.npy', BLS_res)

if len(sys.argv)>1:
    config_name = sys.argv[1]
else:
    config_name = 'configs/createPDG.json'
#############################################################
# Preparation
#############################################################
#### read configs
(dir_LC, dir_res, cat_path, field, ver_LC, ver_search, 
 ver_slc, 
 rnd_ch, G_del, G_select, N_g_inx,  mag_min, mag_max, 
 pers, durs, min_dt, min_N) = readConfig(config_name)

#### create folders
dir_search0 = dir_res + 'search/'
dir_search = dir_search0 + 'search' + ver_search + '/'
dir_BLS = dir_search + 'BLS/'
dir_conf = dir_search + 'configs/'
   
dir_SLC0 = dir_LC + 'SLC/'
dir_SLC = dir_SLC0 + 'SLC' + ver_slc + '/'

for d in [dir_search0, dir_search, dir_BLS, dir_conf]:
    os.makedirs(d, exist_ok=True)
    
#### read reference catalog, star indexes G_inx from folder with LC
F = pd.read_csv(cat_path, sep=' ', names=['ra', 'dec', 'mag', 'source_id'])
F_inx = F.loc[(F.mag>mag_min) & (F.mag<mag_max)].index

G_inx = getGInxLC(dir_SLC, field, ver_LC)    
if rnd_ch:
    G_inx = np.random.choice(G_inx, len(G_inx))
    
### create grids for search (period, duration)
Pers_BLS = pGrid(pers)
Durs_BLS = pGrid([durs])

### G_inx filter -> G_run = (G_select & G_inx & F_inx) - G_del
if len(G_select)==0:
    G_select = G_inx
else:
    G_select = np.array(G_select)
G_run = G_select[~np.isin(G_select, G_del) & np.isin(G_select, F_inx)]

if N_g_inx>0:
    G_run = G_run[:N_g_inx]

#### save configs
for G, txt in zip([G_del, G_select, G_run], ['del', 'select', 'run']):
    print(f'len G_inx {txt} : ' + str(len(G)))

for D, name in zip([Pers_BLS, Durs_BLS], ['P_BLS', 'D_BLS']):
    file = dir_conf + name + ver_search + '.npy'
    if not os.path.isfile(file):
        np.save(file, D)
        print('save ' + file)
    else:
        print(file + ' exists')

#############################################################
# Run
#############################################################
print('\n CREATE PDG: ')
for g_inx in G_run:
    print(g_inx, end=', ')
    createPDG(g_inx, dir_SLC, field, ver_LC, dir_search, 
              Pers_BLS, Durs_BLS, min_dt, min_N
              )

#############################################################
# Save configs
#############################################################
L1 = len(G_run)
L2 = len(os.listdir(dir_BLS))
X = ['len of objects (input)', 'len of objects (output: BLS)',
    '\n\npath to LC for search', 'directory for result',
    'minum lenght of LC for analysis (hours)',
    'minum lenght of LC for analysis (points)',]
Y = [L1, L2, dir_SLC, dir_search, min_dt, min_N]
S = ''
for x, y in zip(X, Y):
    s = x + ' : ' + str(y)
    print(s)
    S += s + '\n'

with open(dir_conf + 'search' + ver_search + '_configs', 'w') as x:
    x.write(S)



