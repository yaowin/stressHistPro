# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 15:29:31 2015
Range Graph Statistics Post-processing
@author: yyao
"""
#%
import os
import sys
import glob
import time
import numpy as np
import pandas as pd

def Usage():
    errorString = """
The expected command line for this script has the form:

    >> python SpaceToTabs.py <input folder> [number of parallel processes]

The parallel process count is optional.
An output folder will be created next to the input folder.
"""
    print errorString

#%%
def get_frame(f, arcs, colnames):
    narc = arcs.shape[0]
    arcname = arcs.index
    df = pd.read_csv(f, sep = '\s', skiprows=1, names = colnames, engine = 'python')
    ds_f = pd.DataFrame(np.zeros([narc,7]), columns = colnames, index =arcname)
    for ind in range(narc):    
        ds = df[(df['Arc_Length'] >= arcs['from'][ind]) & 
            (df['Arc_Length'] < arcs['to'][ind])]
        ds_ext = np.amax(ds)
        ds_ext['Min_Eff_Ten'] = np.amin(ds)['Min_Eff_Ten']
        ds_ext = ds_ext.reshape(1,7)
        ds_f.iloc[ind] = ds_ext # append the arclength ranges max/min results to a frame
    return ds_f

#%%
#-- Obtain the stress history folder list
#---- this'll walk your directories, subfolders recursively and return all folders that contains files.
#
def ProcessCases(curpath):

	fl = [os.path.join(dirpath, f)
		for dirpath, dirnames, files in os.walk(curpath)
		for f in files if f.endswith('.rangegraph.out')]
	fl = np.sort(fl)

	fn = [os.path.split(f)[1].split('.')[0]
		for f in fl]
	#
	#-- Obtain the wave bins information with annually exposure hours per wave bin
	#---- use pandas to read in the excel files that defines each wave bin and its associated exposure hours
	#
	colnames = ['Arc_Length', 'Max_Eff_Ten', 'Min_Eff_Ten', 'Max_Bend_Strs', 'Max_vMis_Strs'
				, 'Max_2RD_Strs', 'Max_2RD_Util']

	arcs = pd.read_csv('test_arcs_input.csv', sep=',', index_col = 0)
	arcname = arcs.index
	print arcs
	frames = [ get_frame(f, arcs, colnames) for f in fl ]
	result = pd.concat(frames, keys = fn)
	result.drop('Arc_Length', axis=1, inplace=True)

	# Create a Pandas Excel writer using XlsxWriter as the engine.
	output = pd.ExcelWriter(os.path.join(curpath, 'Riser RangeGraph Stats.xlsx'), engine='xlsxwriter')

	for an in arcname:
		All = slice(None)
		result_out = result.loc[(All, an), All]
		result_out.to_excel(output, sheet_name= an)
	output.save()
	
	
if __name__ == '__main__':
    # We've come in from a command line, running this script as a program.
    try:
        inputFolder = sys.argv[1]
        goal_dir = os.path.join("Z:\SouthSantaCruz\Strength\Iter#10", inputFolder)
        print goal_dir
    except:
        Usage()
        exit()
    start = time.clock()
    ProcessCases(goal_dir)
    print '\n\nComplete:', int(time.clock()-start), "s"	
