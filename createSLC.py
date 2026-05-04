#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 16:37:01 2024

@author: fedora
"""

import os
import sys
import numpy as np
import pandas as pd

from _funcs import detrendJD, trimSC, prep, readLC, run_multi, clearJD

###############################################################################

def createSLCR(g_inx, regs, min_points, JD_del, CLCF=0, CLCS=0):
    # regs = 'F' / 'S'
    # CLCF!=0 | CLCS!=0
    SLC_list = []
    for reg in regs:
        if reg=='F':
            CLC = CLCF
        else:
            CLC = CLCS
        if len(CLC) > 0:
            JD0 = CLC.loc[0, 'JD']

            # individual
            cols_mag = CLC.columns[1:]
            SLC = CLC.copy()
            SLC[cols_mag] = np.nan
            for c in cols_mag:
                mJD = CLC.loc[~CLC[c].isna(), ['JD', c]].values
                if len(mJD) > min_points:
                    JD = mJD[:, 0]
                    m = mJD[:, 1]  
                    u_cl = clearJD(JD, JD_del)
                    JD = JD[u_cl]
                    m = m[u_cl]

                    # (!!!)
                    m = trimSC(m, 10)
                    t, f, JD, m = prep(JD, m, JD0=JD0, jdm=1)
                    f_dtr = detrendJD(JD, t, f, deg=1, min_points=min_points)  
                    f_dtr = trimSC(f_dtr, 5, 3)
                    # (!!!)
                    SLC.loc[SLC.JD.isin(JD), c] = f_dtr
                else:
                    SLC = SLC.drop(c, axis=1)

            # average by part
            cols_mag_s = SLC.columns[1:]
            parts = cols_mag_s.str[0]
            parts_u = np.unique(parts)
            N = []
            cols_mag_p = []
   
            for p in parts_u:
                if len(parts_u)==1:
                    col_p = 'med_' + reg
                else:
                    col_p = p + '_' + reg
                u = parts==p
                med = SLC[cols_mag_s[u]].median(axis=1)
                SLC[col_p] = med
                N.append(np.sum(~np.isnan(med)))
                cols_mag_p.append(col_p)
            
            # merge parts
            if len(parts_u) > 1:
                col_med = 'med_' + reg
            
                cols_mag_p = pd.Series(cols_mag_p)
                parts = cols_mag_p.str[0]
                parts_u = np.unique(parts)
                inx = np.argsort(N)[::-1]
                
                SLC[col_med] = SLC[cols_mag_p[inx[0]]] 
                for c in cols_mag_p[inx[1:]]:
                    u_nan = SLC[col_med].isna()
                    SLC.loc[u_nan, col_med] = SLC.loc[u_nan, c]
        
            SLC[SLC.columns[1:]] = SLC[SLC.columns[1:]].astype(np.float32)
            SLC_list.append(SLC)

    return SLC_list

def get_regs(CLCF, CLCS):
    u1 = CLCF is not None
    u2 = CLCS is not None

    if u1 & u2 == False:
        if u1==True:
            regs = ['F']
        elif u2==True:
            regs = ['S']
        else:
            regs = []
    else:
        regs = ['F', 'S']
    return regs

def get_cols_main(SLC):
    cols_mag = pd.Series(SLC.columns[1:])
    cols_S = cols_mag[cols_mag.str[-1]=='S'].values
    cols_F = cols_mag[cols_mag.str[-1]=='F'].values
    col_med = cols_mag[cols_mag.str[-3:]=='med'].values
    cols = np.concatenate((['JD'], cols_S, cols_F, col_med))
    return cols


def createSLC(g_inx, dir_CLCF, dir_CLCS, dir_SLC, min_points, JD_del,
              pref_CLCF='CLCF', pref_CLCS='CLCS'):
    CLCF = readLC(dir_CLCF, pref_CLCF, g_inx)
    CLCS = readLC(dir_CLCS, pref_CLCS, g_inx)
    regs = get_regs(CLCF, CLCS)
    
    if len(regs)>0:
        SLC_list = createSLCR(g_inx, regs, min_points, JD_del, CLCF, CLCS)
        
        if len(SLC_list)==2:
            SLC = pd.merge(SLC_list[0], SLC_list[1], on='JD', how='outer')
            SLC['med'] = SLC[['med_F', 'med_S']].median(axis=1)
        else:
            SLC = SLC_list[0]
                    
        if len(SLC.columns)>33:
            cols = get_cols_main(SLC)
            SLC = SLC[cols]
        
        SLC.to_csv(dir_SLC + 'SLC_' + str(g_inx) + '.csv', index=None)


###############################################################################
if len(sys.argv)>1:
    name_ver = sys.argv[1]
    dir_res = sys.argv[2]
else:
    name_ver = '_1'
    dir_res = '/home/fedora/astronomy/STEP/copy/F2/v6/'

dir_CLCS = dir_res + 'CLCS/'
dir_CLCF = dir_res + 'CLCF/'
dir_SLC0 = dir_res + 'SLC/'
dir_SLC = dir_SLC0 + 'SLC' + name_ver + '/'

for d in [dir_SLC0, dir_SLC]:
    if not os.path.isdir(d):
         os.mkdir(d)
         print('create ' + str(d))


min_points = 50
print('__________________CONFIGS__________________')
print('min points for LC : ' + str(min_points))

'''
if field==2:
    JD_del = pd.read_csv('JD_del.csv')
elif field==3:
    JD_del = [[2460417.2485, 2460417.3068],
    [2460424.498, 2460425],
    [2460436.4442, 2460436.4584],
    [2460446.4118, 2460446.4602],
    [2460545.4459, 2460545.559],
    [2460551.5451, 2460551.5619], 
    [2460552.4354, 2460552.4656],
    [2460558.3525, 2460558.4793], 
    [2460559.3729, 2460559.4503], 
    [2460580.2258, 2460580.3139]]
'''
JD_del = []
###############################################################################

files = pd.Series(os.listdir(dir_CLCF))
files2 = pd.Series(os.listdir(dir_CLCS))

G_inx_F = np.int32(files.str[5:-4])
G_inx_S = np.int32(files2.str[5:-4])
G_inx = np.unique(np.append(G_inx_S, G_inx_F))

print(str(len(G_inx)) + ' SLC')

#######################################

print('\n\n\n' + 'create SLC: CLCS, CLCF ---> SLC')

pref_CLCF='CLCF'
pref_CLCS='CLCS'

G_var = [579,   769,  1109,  1993,  2682,  3448,  4062,  4086,  4304,  4615,
    6741,  6889,  9472, 10982, 11669, 11712, 12290, 12890, 12967, 14700,
   14711, 15044, 15747, 15782, 16106, 18856, 22642]
G_run = G_inx #np.random.choice(G_inx, 30, replace=False)

import time
ch = time.time()
run_multi(createSLC, G_run, [dir_CLCF, dir_CLCS, dir_SLC, min_points, JD_del,
                                  pref_CLCF, pref_CLCS], dicts=[])
print(time.time() - ch)
'''
###############################################################################

SLC = pd.read_csv('/home/fedora/astronomy/STEP/copy/F2/v6/SLC/SLC_1/SLC_0.csv')

plotCLC(SLC, dm=0.1)

G_here = np.int32(pd.Series(os.listdir('/home/fedora/astronomy/STEP/copy/F2/v6/SLC/SLC_1/')).str[4:-4])
G = G_inx[~np.isin(G_inx, G_here)]
for g_inx in G:
    createSLC(g_inx, dir_CLCF, dir_CLCS, dir_SLC, min_points, JD_del,
                                      pref_CLCF, pref_CLCS)
    
'''
    
'''
def readLC(dir_LC, pref, g_inx, date=''):
    post = ''
    if len(date)==8:
        post += '_' + date
    file = dir_LC + pref + '_' + str(g_inx) + post + '.csv'
    if os.path.isfile(file):
        LC = pd.read_csv(file)
        if len(LC) > 0:
            return LC
        else:
            print(file + ' len = 0')
    else:
        print(file + ' doesnt exist !!!')
'''   



