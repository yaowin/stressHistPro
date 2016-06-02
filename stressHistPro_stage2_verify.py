# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 17:18:04 2015

@author: Yao

TO-DO: add in log file tracebility, add in exposure hour to file name tracebility (screen and log)
"""

import os
import time
import numpy as np
import OrcFxAPI as orc
import pandas as pd
import xlsxwriter

#%%
#-- Obtain the stress history csv files folder list
#---- this'll walk your directories, subfolders recursively and return all folders that contains files.
#

#curpath = os.path.join(os.getcwd(), 'Stress Histograms')		# use the current folder as where this code is
curpath = "Y:\\OJ Wave Fatigue_1311mWD\\Stress Histograms"			# use your root folder where all sim files are

dirlist = [os.path.join(curpath, 'arcLen ' + str(arcnum) +'m_' + iopos_str)
    for arcnum in [22,]
	for iopos_str in ['Outer',]]

#
#-- Obtain the wave bins information with annually exposure hours per wave bin
#---- use pandas to read in the excel files that defines each wave bin and its associated exposure hours
#
waveData = pd.read_excel(os.path.join(os.getcwd(), 'Fatigue Bins + ProbOccurrence - WindWave_OFI.xlsx'), 
                         'Sheet1', header=1, index_col=0, skip_footer=1, 
                         parse_cols=18, na_values=np.nan)
waveData = waveData.dropna(axis = 1, how='all')
#print waveData.dtypes

Wavebin_hrs = waveData['Annual Exposure Time - Directional.1']
#--- check the column field name to ensure the correct exposure hours

#%%
#-------- Define the necessary inputs for histogram calculations
#-------- Angular positions must be the same as stage 1 angle definition!
#
angle = [0, 90, 180, 270]               

rBins = np.concatenate((np.arange(0,10,0.5),np.arange(10,30,1),
                  np.arange(30,120,5),np.arange(120,250,10))) # Note the unit is MPa

unit_c = 1    # stress history is in 1 MPa = 1 MPa
simhr = 3 # each simulation file has 3 hours
#
#---------------- Specify the output file name for Stress Histogram, Binned/Original rainflow half count distribution
#
OutputDir = 'Y:\\Python Drivers\\Verification'

#%%
#---------------- 1. Extract zz stress time history by column (angular position) for each wave bin (load case).
#		  2. Run rain flow half counts at the time history.
#		  3. Get histogram of the sorted Half Cycle Counts and bin the results.
#		  4. Output the binned results and statistics into Excel file, each angular position is an added worksheet.
#		  5. Output the annual occurence hours weighted counts

print '*** histogram generation begins...'	

cdir = dirlist[0]
fl = [os.path.join(dirpath, f)
    for dirpath, dirnames, files in os.walk(cdir)
    for f in files if f.endswith('.csv')]
    
for col in range(len(angle)):    #scan through all 8 columns in the array
    start = time.clock()
    cycleCounts = np.array([]).reshape(0,len(rBins)-1)
    outfile = os.path.join(OutputDir, 'Stage2 angle_{0:1d}_part3.xlsx'.format(angle[col]))
    workbook   = xlsxwriter.Workbook(outfile)
    cr = 301
    for filename in (fl):
        waveid = int(filename[filename.rfind('_')+1:-4])        
        worksheet = workbook.add_worksheet('waveid_{0:1d}'.format(waveid))
        bin_hrs = Wavebin_hrs[waveid]
        history_data = np.genfromtxt(filename, delimiter=',')
        # read the stress history data from each cvs file
        stressHistory = history_data[:, col + 1] * unit_c    # read and convert stress history into ksi, skip the time stamp column
        halfcyc = orc.RainflowHalfCycles(stressHistory)   # perform rain flow half count
        counts, bins = np.histogram(halfcyc, bins=rBins)  # histogram the stress half cycles
        probCounts = counts / simhr * (bin_hrs)  # exposure hour weighted histogram for each wave bin and stack up
        cycleCounts = np.vstack((cycleCounts, probCounts))
        worksheet.write_row(0, 0, ['waveid','bin hours','stress history', 'rainflow halfcyc','bins','histogram counts','weighted counts'])
        worksheet.write(1, 0, waveid)
        worksheet.write(1, 1, bin_hrs)
        worksheet.write_column(1, 2, stressHistory)
        worksheet.write_column(1, 3, halfcyc)
        worksheet.write_column(1, 4, bins)
        worksheet.write_column(1, 5, counts)
        worksheet.write_column(1, 6, probCounts)
        cr +=1
        print 'Waveid processed: {0:1d} out of {1:1d}  <==> Time elapsed: {2:3d} min'.format(cr, len(fl), int((time.clock()-start)/60))  
    workbook.close()
    sumCounts = np.sum(cycleCounts,axis=0)
print('*****Complete sweeping ' + cdir + '*****')
