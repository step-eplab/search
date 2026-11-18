#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 19:26:52 2024

@author: fedora
"""

import os
import sys
import numpy as np
import pandas as pd

from funcs import readConfig, readLC, plotFLC

if len(sys.argv)>1:
    config_name = sys.argv[1]
else:
    config_name = '../configs/configs_search_V7F5/searchOI.json'
#############################################################
# Preparation
#############################################################
#### read configs
(dir_LC, dir_res, field, ver_LC, ver_search, 
 ver_slc, ver_cor, ver_sel,
 reg_search, 
 max_num,
 v_min, v_max,  p_top, G_select, G_del, N, sgp_min,
 N_max, 
 run, ask, check) = readConfig(config_name)

#### create folders
dir_search = dir_res + 'search/search' + ver_search + '/'
dir_conf = dir_search + 'configs/'
dir_res = dir_search + 'results/'
dir_SLC = dir_LC + 'SLC/SLC' + ver_slc + '/'
dir_FLC = dir_res + 'FLC_BLS' + ver_search + ver_cor + ver_sel+ '/'
dir_FLC_conf = dir_FLC + 'configs/'


for d in [dir_FLC, dir_FLC_conf]:
    os.makedirs(d, exist_ok=True)

#### read selectPDG table
S = pd.read_csv(dir_res + 'Periods_BLS' + ver_search + ver_cor + ver_sel + '.csv', index_col=0)

#############################################################
# Run
#############################################################

import matplotlib
if False:
    matplotlib.use('qt5agg') 
    print('QT5 on')
else:
    matplotlib.use('Agg')
    print('QT5 off')


if len(G_select)>0:
    S = S[S.g_inx.isin(G_select)]

print(len(S))
sv = S.per.value_counts()
print(sv)
if max_num>0: 
    S = S.loc[S.per.isin(sv[sv <= max_num].index)]
S = S.reset_index(drop=True)


if reg_search=='vote':          
    SGP = S.value_counts(['g_inx', 'per'])
    
    inx = np.array(SGP[SGP>=sgp_min].index.tolist())
    if len(inx)>0:
        G_win = np.int32(inx[:, 0])
        P = inx[:, 1]
        S = pd.DataFrame({'g_inx':np.int32(G_win), 'per':P})
    else:
        S = pd.DataFrame(columns=S.columns)
    ver_flc = f'{reg_search}_{len(S)}'
elif reg_search=='Vmin':
    u = (S['val'] >= v_min) & (S['val'] < v_max)
    S = S[u]
    ver_flc = f'{reg_search}_{v_min}_{v_max}_{len(S)}'
elif reg_search=='top':
    n = int(len(S)*(p_top/100))
    inx = S['val'].sort_values(ascending=False).index[:n]
    S = S.loc[inx]
    ver_flc = f'{reg_search}_{p_top}_{len(S)}'
elif reg_search=='max':
    S = S.drop_duplicates('g_inx')
    ver_flc = f'{reg_search}_{len(S)}'
elif reg_search=='N':
    S = S.sort_values('val').loc[S.index[:N]]
    ver_flc = f'{reg_search}_{N}_{len(S)}'
elif reg_search=='all':
    ver_flc = f'{reg_search}_{len(S)}'

dir_save = dir_FLC + ver_flc + f'_maxrep{max_num}/'
os.makedirs(dir_save, exist_ok=True)

G_win = S['g_inx'].unique()
G_win = G_win[~np.isin(G_win, G_del)]
Vars = [ 386,  2300,  2694,  3355,  4041,  5543,  5789, 6469, 7601,
        9588,  9638,  9829, 10142, 10939, 13635, 125]
EA = [ 23,  6651, 5789]
# sin: 386, 2300, 
# EA 23,  6651, 5789


if len(G_win)<100:
    print(G_win)
print(f'\n len(G_WIN):{len(G_win)}')
print(f' len(Vars):{len(G_win[np.isin(G_win,Vars)])}')
print(f' len(EA):{len(G_win[np.isin(G_win,EA)])}')
G_win = G_win[~np.isin(G_win,Vars)]
print(f' len(G_run):{len(G_win)}\n')

if check:
    check_path = dir_conf + 'checked.npy'
    if os.path.isfile(check_path):   
        Check = np.load(check_path, allow_pickle=True).tolist()
    else:
        Check = {}


if ask:
    ans = input('press q for quit\n')
    if ans=='q':
        run = False
S.i_st = S.i_st.astype(int) + 1
i = 0  
if run:
    for g_inx in G_win:
        print(g_inx, end=', ')
        SLC = readLC(dir_SLC, g_inx, field, ver_LC)
        pers, st_is = S.loc[S.g_inx==g_inx, ['per', 'i_st']].T.values
        Fs = []
        for st_i in st_is:
            Fs.append(SLC[f'f{int(st_i)}'].values)
        pers_u = np.unique(pers)
        t = SLC['t'].values
        if check:
            if g_inx in Check:
                pers_u = pers_u[~np.isin(pers_u, Check[g_inx])]
                if len(pers_u)>0:
                    Check[g_inx.item()] += pers_u.tolist()
                else:
                    print('checked')
            else:
                Check[g_inx.item()] = pers_u.tolist()
        if len(pers_u)>0:
            plotFLC(t, Fs, dir_save, g_inx, SLC, pers_u)
            i+=1
else:
    print('\n\nrun: false\n\n')

print('\n len saved: ', i)
if check:
    print(' len checked: ', len(Check))
    np.save(check_path, Check)

'''
S = pd.read_csv(dir_res + 'Periods_BLS' + ver_search + ver_cor + ver_sel + '.csv', index_col=0)
ss = S.groupby('g_inx')['per'].value_counts()
G_win = np.int32(np.array(ss[ss==3].index.tolist())[:, 0])
P = np.array(ss[ss==3].index.tolist())[:, 1]
S = pd.DataFrame({'g_inx':np.int32(G_win), 'per':P})
print(len(S))
sv = S.per.value_counts()
print(sv)
if max_num>0: 
    S = S.loc[S.per.isin(sv[sv <= 10].index)]
S = S.reset_index(drop=True)
'''


'''
#GD = {}
if start:
    for g_inx in G_win:
        pers = S[S.g_inx==g_inx].per.values
        if g_inx in GD:
            pers = pers[~np.isin(pers, GD[g_inx])]
        if len(pers)>0:
            if g_inx in GD:
                GD[g_inx] = np.concatenate((GD[g_inx], pers))
            print(g_inx, end=', ')
            SLC = readLC(dir_SLC, g_inx, field, ver_LC)
            plotFLC(dir_save, g_inx, SLC, pers)
        else:
            print(g_inx)

'''
'''
import pandas as pd
from funcs import getGInxLC

cat_path = '/home/fedora/astronomy/STEP/main/results/catalogs/GAIA/F5.list'

F = pd.read_csv(cat_path, sep=' ', names=['ra', 'dec', 'mag', 'source_id'])
G_win = F.index #[F.mag<15].index
G_here = getGInxLC(dir_SLC)[0]
G_win = G_win[np.isin(G_win, G_here)]


from multiprocessing import Pool, cpu_count

def funpl(g_inx):
    if g_inx%1000==0:
        print(f'{g_inx}  ')
    SLC = readLC(dir_SLC, g_inx, 'F5', 'V7')
    pers = S[S.g_inx==g_inx].per.values[:3]
    plotFLC(dir_save, g_inx, SLC, pers)


args = [(g, ) for g in G_win]
    
with Pool(processes=cpu_count() - 1) as pool:
    pool.starmap(funpl, args)
'''



'''

g_inx = 4608 
SLC = readLC(dir_SLC, g_inx, field, ver_LC)
t, f = SLC['t'], SLC['f']
p_win = 82
t_f, n_f = fold(t, p_win)
for n in np.unique(n_f):
    u = n_f==n
    plt.plot(t_f[u], f[u], '.')

f0, f1 = 0.725, 0.75

u_f = (t_f > f0) & (t_f < f1)

plt.figure()
plt.plot(t, f, '.')
plt.plot(t[u_f], f[u_f], 'x')


 files = pd.Series(os.listdir(dir_LC))
 G_inx = files.str.split('_').str[0]
 G_inx = G_inx.astype(int).values



'''




