# -*- coding: utf-8 -*-
"""
Created on Feb 13 12:22:04 2016
Python code for riser fatigue stress histogram of ECA analysis
Stress histogram is done by 2 stages.
Stage 1: obtain the ZZ stress history at all locations, stored history to CSV files
Stage 2: run rainflow half counts on the stress time history, run histogram analysis, 
combine wave bins based on exposure hours of each wave bin.
This code is the Stage 1 operation.

@author: Yao
"""

import os
import numpy as np
import OrcFxAPI as orc
import time

#%%
#---------------- Define arc length and angle for stress time history extraction 
#
arclength = [2.96, 6.25, 6.71, 10.67, 35.0, 1561.20, 1595.70, 1650.70, 1966.42]
# arc-length of interest (m), consistent with model in (m)
Ey = [114E3, 114E3, 207E3, 207E3, 207E3, 207E3, 207E3, 207E3, 207E3]
# Young's modulus in (MPa)
angle = [0, 45, 90, 135, 180, 225, 270, 315]
# angular position (degrees)
rploc = [0, 0.5, 1.0]
# location for stress history; 0 is inner surface and 1 is outer surface and 0.5 is midwall
rpstr =['_Inner', '_Midwall', '_Outer']
# the string text go with the position
od = np.array(np.tile(0.219075,len(arclength)))
# the riser pipe OD in (m)
od[0]=0.262631    # the TSJ taper section OD calculated in (m)
d = 0.161671    # the riser pipe ID in (m)
t = (od-d)/2    # the riser wall thickness calculated
Ac = np.pi/4 * (od**2 - d**2)    # the riser cross section area
Ixy = np.pi/64 * (od**4 - d**4)   # the riser area moment of inertia

#%%
#-- Obtain the sim file list and print the total number of sim files to be processed
#---- this'll walk your directories, subfolders recursively and return all absolute pathnames to matching .sim files.
#

#curpath = os.getcwd()		# use the current folder as where this code is
curpath = "Y:\\"			# use your root folder where all sim files are
fl = [os.path.join(dirpath, f)
    for dirpath, dirnames, files in os.walk(curpath)
    for f in files if f.endswith('.sim')]

fn = [f[f.rfind('\\')+1:]
    for f in fl]

print 'Total Wave Bin (sim files) number in ' + curpath + ' is ' + str(len(fl))

#
#------ Define period of simulation and variable for data to be extracted
#
simPeriod = orc.Period(1)   # use the stage 1 period of simulation
            # pnSpecifiedPeriod, pnLatestWave, pnWholeSimulation, pnStaticState

LineName = "8.625x1.13in_EN3"       # use the riser/flowline name in the model

#%%
#------ Main body of stage 1 process - create stress history files *.CSV
#
ind = 0
print 'data extraction begins...'	
for filename in (fl):
    ind +=1
    start = time.clock()
    waveid = filename[filename.rfind('_')+1:-4]
    model = orc.Model(filename, threadCount=8)
    riserLine = model[LineName]
    timeHist = model.general.TimeHistory("Time", period=simPeriod)
    timeHist = np.array([timeHist])
    
    for arcl in (arclength):
        Eps = riserLine.TimeHistory("Direct Tensile Strain", period=simPeriod, 
                                    objectExtra=orc.oeArcLength(arcl)
                                    )
        Cx = riserLine.TimeHistory("x-Curvature", period=simPeriod, 
                                    objectExtra=orc.oeArcLength(arcl)
                                    )
        Cy = riserLine.TimeHistory("y-Curvature", period=simPeriod, 
                                    objectExtra=orc.oeArcLength(arcl)
                                    )
		
        Cr = np.outer(Cx, np.sin(np.radians(angle))) - np.outer(Cy, np.cos(np.radians(angle)))
        # Resultant curverature based on theta and two x- and y- component, -->
		# reference OrcaFlex Help: Tension and Curvature Stress Factors (stress factors only) <--
        for rp in (rploc):
            TenS = np.dot(Ey[arclength.index(arcl)], Eps/100) 
			#Direct Tensile stress = tensile strain times modulus, the tensile strain was in percentage
            BenS = np.dot(Cr, Ey[arclength.index(arcl)]*(d/2 + rp*t[arclength.index(arcl)]))
			#Bending stress = resultant curvarture times modulus times distance from neutrual axis to stress fibre
            stressHist = np.transpose(np.tile(TenS, (np.size(angle), 1))) + BenS
			#ZZ stress = Tensile Stress + Bending Stress
            Histories = np.hstack((timeHist.T, stressHist))
            # Append to the time history column, end up with 9 columns data matrix
			
            subfoldpath = curpath +'\\Stress Histograms' + '\\arcLen ' + str(int(np.rint(arcl))) +'m' + rpstr[rploc.index(rp)]
            if not os.path.exists(subfoldpath):
                os.makedirs(subfoldpath)
            textFstr = subfoldpath + '\\StressHistory at arcLen=' + str(int(np.rint(arcl))) +'m' + rpstr[rploc.index(rp)] + ' wavebin_' + waveid + '.csv'
            np.savetxt(textFstr, Histories, delimiter=",")
	print 'Wavebin processed: {0:3d} out of {1:3d}  <==> Time elapsed: {2:2d} min; Time to go: {3:4d} min'.format(ind, len(fl), int((time.clock()-start)/60), int((time.clock()-start)*(len(fl)-ind)/60))

#
#----[FINAL]------- Obtain xx arcLen X 2 inner/outer = xx total subfolders under '/Stress Histograms/,
#------------------ each subfolders will contain xxx (wave bins) csv files
#------------------ each csv file will contain 1 time trace, and 8 (angular positions) columns of stress history data.
#
