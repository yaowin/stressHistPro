# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 12:22:04 2015

@author: Yao
"""

import os
import numpy as np
import OrcFxAPI as orc
import time
#%%
#---------------- Define arc length and angle for time history extraction 
#
arclength = [2.96, 6.25, 6.71, 10.67, 35.0, 1561.20, 1595.70, 1650.70, 1966.42]
# arc-length of interest (m), but model is built on (m), pay attention to Unit!
angle = [0, 45, 90, 135, 180, 225, 270, 315]               # angular position (degrees)
rploc = [0, 1]                          # location for stress history; 0 is inner surface and 1 is outer surface and 0.5 is midwall
rpstr =['_Inner', '_Midwall', '_Outer']				# the string text go with the position

#%%
#-- Obtain the sim file list and print the total number of sim file to be processed
#---- this'll walk your directories, subfolders recursively and return all absolute pathnames to matching .sim files.
#
curpath = os.getcwd()
curpath = "Y:\\"
fl = [os.path.join(dirpath, f)
    for dirpath, dirnames, files in os.walk(curpath)
    for f in files if f.endswith('.sim')]

fn = [f[f.rfind('\\')+1:]
    for f in fl]

print 'Total sim file number in ' + curpath + ' is ' + str(len(fl))

#
#------ Define period of simulation and variable for data to be extracted, and bin number in the following histogram calculation
#
simPeriod = orc.Period(1)   # use the stage 1 period of simulation
            # pnSpecifiedPeriod, pnLatestWave, pnWholeSimulation, pnStaticState
VariableName = "ZZ Stress"          # set the ZZ stress component extraction target
LineName = "8.625x1.13in_EN3"       # set the riser/flowline name in the model

#%%
#------ Main body of stage 1 process - create stress history files
#
ind = 0
print 'data extraction begins...'	
for filename in (fl[:1]):
    ind +=1
    start = time.clock()
    waveid = filename[filename.rfind('_')+1:-4]
    model = orc.Model(filename, threadCount=8)
    riserLine = model[LineName]
    timeHist = model.general.TimeHistory("Time", period=simPeriod)
    #vmstressRG = riserLine.RangeGraph("Max von Mises Stress", period=simPeriod)
    for arcl in (arclength):
        for rp in (rploc):
            Histories = np.array(timeHist)
            for ang in (angle):
                stressHist = riserLine.TimeHistory(VariableName, period=simPeriod, 
                objectExtra=orc.oeLine(
                ArcLength= arcl,
                Theta=ang,
                RadialPos=rp)
                )
                Histories = np.vstack((Histories, stressHist))   
                # for each eight angular position, extract the stress history, append to the time history column, end up with 9 columns in total
            Histories = np.transpose(Histories)         
            # get the stress history array in shape
            subfoldpath = curpath +'\\Stress Histograms' + '\\arcLen ' + str(int(arcl)) +'ft' + rpstr[rp]
            if not os.path.exists(subfoldpath):
                os.makedirs(subfoldpath)
            textFstr = subfoldpath + '\\StressHistory at arcLen=' + str(int(arcl)) +'ft' + rpstr[rp] + ' wavebin_' + waveid + '.csv'
            np.savetxt(textFstr, Histories, delimiter=",")
    print 'Simulation file number processed: ' + str(ind)
    print '\nCompletion time: ', int(time.clock()-start), "s"	
#
#----[FINAL]------- Obtain xx arcLen X 2 inner/outer = xx total subfolders under '/Stress Histograms/,
#------------------ each subfolders will contain xxx (wave bins) csv files
#------------------ each csv file will contain 1 time trace, and 8 (angular positions) columns of stress history data.
#
