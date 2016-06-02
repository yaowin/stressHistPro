# -*- coding: utf-8 -*-
"""
Created on Mar 15 12:22:04 2016


@author: Yao
"""

import os
import numpy as np
import time
import xlsxwriter
#%%

def filename_screen(fn):
    if fn.find('Loop') < 0:
        cat = 'Background'
    else:
        cat = 'Loop'
    
    n = fn[fn.find('Case')+4 : fn.find('.')]
    
    return [cat, n]
    
#%%
#-- Obtain the sim file list and print the total number of sim files to be processed
#---- this'll walk your directories, subfolders recursively and return all absolute pathnames to matching .sim files.
#

#curpath = os.getcwd()		# use the current folder as where this code is
curpath = ["Y:\\Odd Job VIV Fatigue_Lankhorst_w_str_dmg\\Background_Current_VIV_3500ft\\",
           "Y:\\Odd Job VIV Fatigue_Lankhorst_w_str_dmg\\Loop_Current_VIV_3500ft\\"]# use your root folder where all sim files are

arclength = [2.96, 6.25, 6.71, 10.67, 35.0, 1561.20, 1595.70, 1650.70, 1966.42]
# arc-length of interest (m), consistent with model in (m)

fl1 = [os.path.join(dirpath, f)
    for dirpath, dirnames, files in os.walk(curpath[0])
    for f in files if f.endswith('.str')]

fl2 = [os.path.join(dirpath, f)
    for dirpath, dirnames, files in os.walk(curpath[1])
    for f in files if f.endswith('.str')]

fl = fl1 + fl2

fn = [f[f.rfind('\\')+1:]
    for f in fl]

md1 = [os.path.join(dirpath, f)
    for dirpath, dirnames, files in os.walk(curpath[0])
    for f in files if f.endswith('.mds')]

md2 = [os.path.join(dirpath, f)
    for dirpath, dirnames, files in os.walk(curpath[1])
    for f in files if f.endswith('.mds')]

md = md1 + md2

mds = [f[f.rfind('\\')+1:]
    for f in md]
        
segment = np.genfromtxt('RiserArcLength.dat')
x = np.array([])
for sta, end, nstep in segment:
    x = np.append(x , np.arange(sta + (end-sta)/nstep/2, end, (end-sta)/nstep))
x = np.append(x, end)

x_L = x/np.max(x)
arc_L = arclength/np.max(x)

exp_hour = np.genfromtxt('VIV_Exposure_Hours.dat', usecols =(1))

print 'Total VIV case number: {0:3d} Background Current + {1:3d} Loop Current'.format(len(fl1),len(fl2))
#%%
#------ Main body of stage 1 process - create stress history files *.CSV
#
f_ind = 0
x_id = 0
print 'data extraction begins...'	


stress_results = np.array([]).reshape(0, len(arc_L) + 4)
workbook = xlsxwriter.Workbook('VIV stress results.xlsx')
worksheet = workbook.add_worksheet('results')
header = ['currentType', 'currentCase#', 'VIV_Mode#','Mode Freq','Mode Prob','Number of Cycles']
worksheet.write_row(1, 0, header)
worksheet.write_row(1, 6, arclength)

for filename in (fl):
    f_ind +=1
    start = time.clock()
    caseCat, caseN = filename_screen(filename)
    md_fn = filename.replace('S7', 'Common')
    md_fn = md_fn.replace('.str', '.mds')
    
    md_range = np.genfromtxt(filename,max_rows = 1)
        # read the lowest and highest mode numbers for which Sr,rms is listed
    md_prob = np.genfromtxt(filename, skip_header = 1, max_rows = 1)
        # read the listing of the time sharing probabilities for each of the above modes
    stress_op = np.genfromtxt(filename, skip_header = 2, max_rows = len(x_L))
    stress_op = stress_op / 1E6 # Pa -> MPa
        # read the listing Sr,rms array for the whole line
    md_freq = np.genfromtxt(md_fn, skip_header = 1, usecols =(0,1), max_rows = 200 )
        # read the full 200 modes frequency
    md_ind = np.arange(md_range[0], md_range[1]+1)
    if caseCat == 'Loop':
        T_hr = exp_hour[int(caseN)-1 + 8]
    else:
        T_hr = exp_hour[int(caseN)-1]

    m_id = 0
    for mi in md_ind:
        freq = md_freq[mi-1, 1]
        N = T_hr * 3600 * freq / 2 / np.pi * md_prob[m_id]
        stress_arc = np.interp(arc_L, x_L, stress_op[:, m_id + 1])
        stress_arc = np.append([mi, freq, md_prob[m_id], N,], stress_arc)
        worksheet.write_row(x_id + 2 , 2, stress_arc)
        worksheet.write_string(x_id + 2 , 1, caseN)
        worksheet.write_string(x_id + 2 , 0, caseCat)        
        stress_results = np.vstack((stress_results, stress_arc))
        m_id +=1
        x_id +=1
    print 'Extreme case processed: {0:3d} out of {1:3d}  <==> Time elapsed: {2:2d} min; Time to go: {3:4d} min'.format(f_ind, len(fl), int((time.clock()-start)/60), int((time.clock()-start)*(len(fl)-f_ind)/60))
np.savetxt('VIV_stress_results.csv', stress_results, delimiter=",")
workbook.close()

#
#----[FINAL]------- Obtain xx arcLen X 2 inner/outer = xx total subfolders under '/Stress Histograms/,
#------------------ each subfolders will contain xxx (wave bins) csv files
#------------------ each csv file will contain 1 time trace, and 8 (angular positions) columns of stress history data.
#
