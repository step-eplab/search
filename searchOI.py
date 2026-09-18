#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 19:26:52 2024

@author: fedora
"""

import os
import sys
import pandas as pd

from funcs import readConfig, readLC, plotFLC

if len(sys.argv)>1:
    config_name = sys.argv[1]
else:
    config_name = 'configs/searchOI.json'
#############################################################
# Preparation
#############################################################
#### read configs
(dir_LC, dir_res, field, ver_LC, ver_search, 
 ver_slc, ver_cor, ver_sel,
 reg_search, v_min,  p_top, G_select, N, N_max,show) = readConfig(config_name)

#### create folders
dir_search = dir_res + 'search/search' + ver_search + '/'
dir_conf = dir_search + 'configs/'
dir_res = dir_search + 'results/'
dir_SLC = dir_LC + 'SLC/SLC' + ver_slc + '/'
dir_FLC = dir_res + 'FLC_BLS' + ver_search + ver_cor + '/'
dir_FLC_conf = dir_FLC + 'configs/'

for d in [dir_FLC, dir_FLC_conf]:
    os.makedirs(d, exist_ok=True)

#### read selectPDG table
S = pd.read_csv(dir_res + 'Periods_BLS' + ver_search + ver_cor + ver_sel + '.csv', index_col=0)


#############################################################
# Run
#############################################################

import matplotlib
if show:
    matplotlib.use('qt5agg') 
    print('QT5 on')
else:
    matplotlib.use('Agg')
    print('QT5 off')


if len(G_select)>0:
    S = S[S.g_inx.isin(G_select)]

if reg_search=='Vmin':
    u = S.val > v_min
    S = S[u]
    ver_flc = f'{reg_search}_{v_min}_{len(S)}'
elif reg_search=='top':
    S = S.loc[S.val.nlargest(int(len(S)*(p_top/100))).index]
    ver_flc = f'{reg_search}_{p_top}_{len(S)}'
elif reg_search=='max':
    S = S.drop_duplicates('g_inx')
    ver_flc = f'{reg_search}_{len(S)}'
elif reg_search=='N':
    S = S.sort_values('val').loc[S.index[:N]]
    ver_flc = f'{reg_search}_{N}_{len(S)}'

dir_save = dir_FLC + ver_flc
os.makedirs(dir_save, exist_ok=True)

G_win = S['g_inx'].unique()
for g_inx in G_win:
    print(f'\ng_inx: {g_inx}')
    SLC = readLC(dir_SLC, g_inx, field, ver_LC)
    pers = S[S.g_inx==g_inx].per.values
    plotFLC(dir_save, g_inx, SLC, pers)





'''
if reg_search!='S':
    N = S[['g_inx','per']].groupby(['g_inx','per']).value_counts()
    
    pv = pd.Series(N.values).value_counts()
    inx = pv[(pv < N_max)==True].index
    
    if len(inx)==0:
        inx = [max(N)]
    
    print(pv)
    print(inx)
    print('n = ' + str(n))
    
    dir_save = dir_FLC + 'N_' + str(n)
    os.makedirs(dir_save, exist_ok=True)
    
    u_n = N==n
    s = N[u_n]
    GP = N.index[u_n]
    G = []
    Sgs = []
    for gp in GP:
        g_inx, per = gp
        sg = S[(S.g_inx==g_inx) & (S.per==per)]  
        G.append(g_inx)
        Sgs.append(sg)
    
    
    print('len obj : ' + str(len(GP)))
    print('save FLC to ' + dir_save)
    ij = 0
    for g_inx, sg in zip(G, Sgs):
        print(g_inx, end=',')
        SLC = readLC(dir_SLC, 'SLC', g_inx)
        
        if not SLC is None:
            cols, pers = sg['col'].values, sg['per'].values
            
            plotFLC(SLC, pers, cols)
            plt.suptitle(str(g_inx))
            plt.tight_layout()
            plt.savefig(dir_save + '/' + str(g_inx) + '_' + str(round(pers[0], 2)) + '.png')
            plt.close()
    
            if ij%10==0:
                print(ij, end=', ')
            ij+=1
    
    print('len obj G_N : ' + str(len(G)))
    print('save FLC to ' + dir_save)
    
    np.save(dir_FLC_conf + 'G_N_' + str(n) + '.npy', G)

######################################

if reg_search!='G':
    pv = pd.Series(np.round(S.val.values)).value_counts()
    print(pv)
    
    pc = pv[::-1].cumsum()
    
    
    u_s = (S.val > v_min)
    print(sum(u_s))
    print(len(np.unique(S.g_inx[u_s])))
    
    dir_save = dir_FLC + 'S_' + str(v_min)
    if not os.path.isdir(dir_save):
        os.mkdir(dir_save)
        
    G_s = np.unique(S.loc[u_s, 'g_inx'])
    print('len obj G_S : ' + str(len(G_s)))
    print('save FLC to ' + dir_save)
    
    
    ij = 0
    for g_inx in G_s:
        SLC = readLC(dir_SLC, 'SLC', g_inx)
        if not SLC is None:
            sg = S[(S.g_inx==g_inx) & u_s]
            
            cols, pers = sg['col'].values, sg['per'].values
    
            plotFLC(SLC, pers, cols)
            plt.suptitle(str(g_inx))
            plt.tight_layout()
            plt.savefig(dir_save + '/' + str(g_inx) + '_' + str(round(pers[0], 2)) + '.png')
            plt.close()  
                
        if ij%10==0:
            print(ij, end=', ')
    
    np.save(dir_FLC_conf + 'G_S_' + str(v_min) + '.npy', G_s)
'''
