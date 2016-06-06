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
arclength = [0.00, 2.96, 5.79, 7.39, 7.70, 9.52, 11.66, 42.06, 2344.91, 2663.82, 2781.70, 2917.40]
# arc-length of interest (m), but model is built on (ft), pay attention to Unit!
seed = [9, 9, 9, 9, 9, 9, 9, 9, 1, 2, 2, 2]
angle = [0, 45, 90, 135, 180, 225, 270, 315]               # angular position (degrees)
rploc = [0, 1]                          # location for stress history; 0 is inner surface and 1 is outer surface and 0.5 is midwall
rpstr =['_Inner', '_Outer']				# the string text go with the position

#%%
#-- Obtain the sim file list and print the total number of sim file to be processed
#---- this'll walk your directories, subfolders recursively and return all absolute pathnames to matching .sim files.
#

curpath = "Z:\\SouthSantaCruz\\Wave Fatigue\\Short Term"

fl = [os.path.join(dirpath, f)
    for dirpath, dirnames, files in os.walk(curpath)
    for f in files if f.endswith('.sim')]

print 'Total sim file number is ' + str(len(fl))

#
#------ Define period of simulation and variable for data to be extracted, and bin number in the following histogram calculation
#
simPeriod = orc.Period(1)   # use the stage 1 period of simulation
            # pnSpecifiedPeriod, pnLatestWave, pnWholeSimulation, pnStaticState
VariableName = "ZZ Stress"          # set the ZZ stress component extraction target
LineName = "SSC_Riser"       # set the riser/flowline name in the model

#%%
#------ Main body of stage 1 process - create stress history files
#
ind = 0
print 'data extraction begins...'	
for arcl in (arclength):
    start = time.clock()
    curfl = [f 
            for f in fl if f.endswith(str(seed[ind]) + '.sim')] 
    for filename in (curfl):
        wind = filename[filename.find('Wind_'):filename.find('Wind_')+8]
        binid = filename[filename.find('Bin'):filename.find('Bin')+6]
        model = orc.Model(filename, threadCount=8)
        riserLine = model[LineName]
        timeHist = model.general.TimeHistory("Time", period=simPeriod)
		#vmstressRG = riserLine.RangeGraph("Max von Mises Stress", period=simPeriod)
        for rp in (rploc):
            Histories = np.array(timeHist)
            for ang in (angle):
                stressHist = riserLine.TimeHistory(VariableName, period=simPeriod, 
                objectExtra=orc.oeLine(
                ArcLength= arcl*3.281,
                Theta=ang,
                RadialPos=rp)
                )
                Histories = np.vstack((Histories, stressHist*0.04788))   
                # for each eight angular position, extract the stress history,
                # append to the time history column, end up with 9 columns in total
                # unit conversion from meter to ft in arc length
                # unit conversion from ksf to MPa
            Histories = np.transpose(Histories)         
            # get the stress history array in shape
            subfoldpath = curpath +'\\Stress Histograms' + '\\arcLen ' + str(int(arcl)) +'m' + rpstr[rp]
            if not os.path.exists(subfoldpath):
                os.makedirs(subfoldpath)
            textFstr = os.path.join(subfoldpath, "_".join((wind,binid)) + '.csv')
            np.savetxt(textFstr, Histories, delimiter=",")
    ind +=1 
    print 'Simulation file number processed: ' + str(ind)
    print 'Completion time in this arc length step: ', int(time.clock()-start), "s\n"	
#
#----[FINAL]------- Obtain xx arcLen X 2 inner/outer = xx total subfolders under '/Stress Histograms/,
#------------------ each subfolders will contain xxx (wave bins) csv files
#------------------ each csv file will contain 1 time trace, and 8 (angular positions) columns of stress history data.
#
