#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 19:26:52 2024

@author: fedora
"""

import os
import json
import numpy as np
import pandas as pd


import matplotlib.pyplot as plt

import matplotlib
matplotlib.use('Agg')  
print('QT5 off')
#matplotlib.use('qt5agg') 

from funcs import trimSC, fold, fit_bin, readLC

def plotFLC(SLC, pers, cols):
    t = (SLC.JD - min(SLC.JD))*24
    Nc = len(cols)          
    fig, Ax = plt.subplots(Nc, 1, figsize=(16, 10))
    if Nc==1:
        Ax = [Ax]
    for c, ax, p_win in zip(cols, Ax, pers):
        f = SLC[c]
        f = trimSC(f, 3, 1)
        t_f, n_f = fold(t, p_win)
        ax.plot(t_f, f, '.', alpha=0.5)
        for n, clr, lnw in zip([64, 164, 512], ['k', 'r', 'lime'], [3, 2, 1]):
            x_av, y_av = fit_bin(t_f, f, n)
            ax.plot(x_av, y_av, color=clr, linewidth=lnw)
            ax.grid()
        ax.set_title(c + ', per = ' + str(round(p_win, 4)) + 'h')

###############################################################################
with open('configs.json', 'r') as file:
    configs = json.load(file)

dir_LC = configs['dir_LC']
dir_res = configs['dir_res']

ver_search = configs['ver_search']
ver_slc = configs['ver_slc']
ver_cor = configs['ver_cor']

N_max = configs['N_max']

reg_search = configs['reg_search']

n = configs['n'] #inx[0]
v_min = configs['v_min'] #int(pc[pc <= N_max].index[::-1][0])

dir_search = dir_res + 'search/search' + ver_search + '/'
dir_conf = dir_search + 'configs/'
dir_res = dir_search + 'results/'
dir_SLC = dir_LC + 'SLC/SLC' + ver_slc + '/'
dir_FLC = dir_res + 'FLC_BLS' + ver_search + ver_cor + '/'
dir_FLC_conf = dir_FLC + 'configs/'
for d in [dir_FLC, dir_FLC_conf]:
    os.makedirs(d, exist_ok=True)

S = pd.read_csv(dir_res + 'res_BLS' + ver_search + ver_cor + '.csv')
##################################

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

###############################################################################
    