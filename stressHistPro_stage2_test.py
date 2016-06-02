# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 17:18:04 2015

@author: Yao

TO-DO: add in log file tracebility, add in exposure hour to file name tracebility (screen and log)
"""

import os
import numpy as np
import OrcFxAPI as orc
import pandas as pd
import xlsxwriter
import time

#%%
#-- Obtain the stress history csv files folder list
#---- this'll walk your directories, subfolders recursively and return all folders that contains files.
#
curpath = os.path.join(os.getcwd(), 'Stress Histograms')
#curpath = os.path.join(os.getcwd(), 'Stress Histograms')		# use the current folder as where this code is
curpath = "Y:\\Stress Histograms"			# use your root folder where all sim files are

dirlist = [os.path.join(curpath, 'arcLen ' + str(arcnum) +'m_' + iopos_str)
    for arcnum in [3, 6, 7,]
	for iopos_str in ['Inner', 'Midwall',]]

#
#-- Obtain the wave bins information with annually exposure hours per wave bin
#---- use pandas to read in the excel files that defines each wave bin and its associated exposure hours
#
waveData = pd.read_excel('Fatigue Bins + ProbOccurrence - WindWave_OFI.xlsx', 
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
angle = [0, 45, 90, 135, 180, 225, 270, 315]               

rBins = np.concatenate((np.arange(0,10,0.5),np.arange(10,30,1),
                  np.arange(30,120,5),np.arange(120,250,10))) # Note the unit is MPa

unit_c = 1    # stress history is in 1 MPa = 1 MPa
simhr = 3 # each simulation file has 3 hours
#
#---------------- Specify the output file name for Stress Histogram, Binned/Original rainflow half count distribution
#
OutputStr = 'Histogram Binned Annual Probability Weighed test.xlsx'

#%%
#---------------- 1. Extract zz stress time history by column (angular position) for each wave bin (load case).
#		  2. Run rain flow half counts at the time history.
#		  3. Get histogram of the sorted Half Cycle Counts and bin the results.
#		  4. Output the binned results and statistics into Excel file, each angular position is an added worksheet.
#		  5. Output the annual occurence hours weighted counts

print '*** histogram generation begins...'	
workbook = xlsxwriter.Workbook(os.path.join(curpath,OutputStr))
bold = workbook.add_format({'bold': True})
boldred = workbook.add_format({'bold': True, 'font_color': 'red'})
scifm = workbook.add_format({'num_format': 0x0b})
intfm = workbook.add_format({'num_format': '#,##0'})
digfm = workbook.add_format({'num_format': '0.00'})

worksheet = workbook.add_worksheet('Cyclic Stress Histogram')
worksheet.set_column('A:D', 15)
worksheet.write('D1', 'Arclength (ft from hang-off)', bold)
worksheet.write('D2', 'Angular Position', bold)
worksheet.write('A3', 'Bin number:', bold)
worksheet.write('B3', 'Min Bin Stress (MPa)', bold)
worksheet.write('C3', 'Max Bin Stress (MPa)', bold)
worksheet.write('D3', 'Bin Size (MPa)', bold)

i_col = 4
for cdir in dirlist:
    fl = [os.path.join(dirpath, f)
        for dirpath, dirnames, files in os.walk(cdir)
        for f in files if f.endswith('.csv')]
    
    worksheet.write(0, i_col, cdir[cdir.rfind('\\')+1:], bold)

    for col in range(len(angle)):    #scan through all 8 columns in the array
        start = time.clock()
        worksheet.write_number(1, i_col, angle[col], intfm)
        worksheet.write(2, i_col, 'Annual Counts', bold)        
        cycleCounts = np.array([]).reshape(0,len(rBins)-1)
        for filename in (fl):
            waveid = int(filename[filename.rfind('_')+1:-4])
            bin_hrs = Wavebin_hrs[waveid]
            history_data = np.genfromtxt(filename, delimiter=',')
            
            # read the stress history data from each cvs file
            stressHistory = history_data[:, col + 1] * unit_c    # read and convert stress history into ksi, skip the time stamp column
            halfcyc = orc.RainflowHalfCycles(stressHistory)   # perform rain flow half count
            counts, bins = np.histogram(halfcyc, bins=rBins)  # histogram the stress half cycles
            probCounts = counts / simhr * (bin_hrs)  # exposure hour weighted histogram for each wave bin and stack up
            cycleCounts = np.vstack((cycleCounts, probCounts))
            #print 'wave bin ' + str(waveid) + ': done'
        sumCounts = np.sum(cycleCounts,axis=0)
            #------------------ print data column headers and statistical results
            #------------------ print binned histogram results
        i_row = 3
        for bind in range(len(sumCounts)):
            bmin = rBins[bind]
            bmax = rBins[bind+1]
            bmean = np.mean([bmin, bmax])
            bSize = bmax - bmin
            bN = sumCounts[bind]
            if i_col < 5:
                worksheet.write_number(i_row, 0,     bind+1,  intfm )
                worksheet.write_number(i_row, 1,     bmin,   digfm )
                worksheet.write_number(i_row, 2,     bmax,   digfm )
                worksheet.write_number(i_row, 3,     bSize,   digfm )
            worksheet.write_number(i_row, i_col,     bN,    scifm )

            i_row += 1
        i_col +=1

        print 'Angular Position processed: {0:1d} out of {1:1d}  <==> Time elapsed: {2:3d} min'.format(col+1, len(angle), int((time.clock()-start)/60))
    print('*****Complete sweeping ' + cdir + '*****')

workbook.close()