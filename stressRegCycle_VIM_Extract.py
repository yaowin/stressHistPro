# -*- coding: utf-8 -*-
"""
Created on Mar 15 12:22:04 2016


@author: Yao
"""

import os
import numpy as np
import OrcFxAPI as orc
import time
import xlsxwriter
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
curpath = ["Y:\\Odd Job VIM FAtigue_1311mWD\\Background_Current_VIM\\",
           "Y:\\Odd Job VIM FAtigue_1311mWD\\Loop_Current_VIM\\"]# use your root folder where all sim files are
fl1 = [os.path.join(dirpath, f)
    for dirpath, dirnames, files in os.walk(curpath[0])
    for f in files if f.endswith('.sim')]

fl2 = [os.path.join(dirpath, f)
    for dirpath, dirnames, files in os.walk(curpath[1])
    for f in files if f.endswith('.sim')]
    
fl = fl1 + fl2

fl = np.sort(fl)

fn = [f[f.rfind('\\')+1:]
    for f in fl]

print 'Total VIM case number: ' + str(len(fl))

#
#------ Define period of simulation and variable for data to be extracted
#
simPeriod = orc.pnLatestWave   # use the latest wave period of simulation
            # Period(1), pnSpecifiedPeriod, pnLatestWave, pnWholeSimulation, pnStaticState

LineName = "8.625x1.13in_EN3"       # use the riser/flowline name in the model

#%%
#------ Main body of stage 1 process - create stress history files *.CSV
#
ind = 0
print 'data extraction begins...'	
ext_results_max = np.array([]).reshape(0,8*3+1)
ext_results_min = np.array([]).reshape(0,8*3+1)

for filename in (fl):
    ind +=1
    start = time.clock()
    waveid = filename[filename.rfind('\\')+1:filename.rfind('_VIM')]
    model = orc.Model(filename, threadCount=8)
    riserLine = model[LineName]
    
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

        TenS = np.dot(Ey[arclength.index(arcl)], Eps/100) 

        # Resultant curverature based on theta and two x- and y- component, -->
		# reference OrcaFlex Help: Tension and Curvature Stress Factors (stress factors only) <--
        maxSS = np.array([waveid,])
        minSS = np.array([waveid,])
        
        for rp in (rploc):

			#Direct Tensile stress = tensile strain times modulus, the tensile strain was in percentage
            BenS = np.dot(Cr, Ey[arclength.index(arcl)]*(d/2 + rp*t[arclength.index(arcl)]))
			#Bending stress = resultant curvarture times modulus times distance from neutrual axis to stress fibre
            stressHist = np.transpose(np.tile(TenS, (np.size(angle), 1))) + BenS
			#ZZ stress = Tensile Stress + Bending Stress
            
            # Append to the time history column, end up with 9 columns data matrix
            maxSS = np.append(maxSS, np.amax(stressHist, axis=0))
            minSS = np.append(minSS, np.amin(stressHist, axis=0))
        ext_results_max = np.vstack((ext_results_max, maxSS))
        ext_results_min = np.vstack((ext_results_min, minSS))
    print 'VIM case processed: {0:3d} out of {1:3d}  <==> Time elapsed: {2:2d} min; Time to go: {3:4d} min'.format(ind, len(fl), int((time.clock()-start)/60), int((time.clock()-start)*(len(fl)-ind)/60))
np.savetxt('VIM_results_max.csv', ext_results_max, fmt='%15s '*25, delimiter=",")
np.savetxt('VIM_results_min.csv', ext_results_min, fmt='%15s '*25, delimiter=",")

#%%
#------ Main body of stage 2 process - create stress history files *.CSV
#

workbook = xlsxwriter.Workbook('VIM results.xlsx')

for ind in range(len(arclength)):
    worksheet = workbook.add_worksheet('Arclength ' + str(int(np.rint(arclength[ind]))))
    irow = 0
    maxSS = ext_results_max[ind::9,:]
    minSS = ext_results_min[ind::9,:]
    for row in maxSS:
        worksheet.write_row(irow, 1, row)
        irow+=1
    irow+=1
    for row in minSS:
        worksheet.write_row(irow, 1, row)
        irow+=1
workbook.close()
#
#----[FINAL]------- Obtain xx arcLen X 2 inner/outer = xx total subfolders under '/Stress Histograms/,
#------------------ each subfolders will contain xxx (wave bins) csv files
#------------------ each csv file will contain 1 time trace, and 8 (angular positions) columns of stress history data.
#
