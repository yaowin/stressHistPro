# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 15:29:31 2015
Script to Create a folder 'Associated Stats' to have all csv files
containing the linked stats from the list of 
['walTenRoot', 'shearFRoot', 'bendMomRoot'
                , 'walTenTip', 'shearFTip', 'bendMomTip']
in each of the *.tsj.out file and major maximum value are in that order!
@author: yyao
"""

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
#-- Obtain the stress history folder list
#---- this'll walk your directories, subfolders recursively and return all folders that contains files.
#
def ProcessCases(curpath):

	fl = [os.path.join(dirpath, f)
		for dirpath, dirnames, files in os.walk(curpath)
		for f in files if f.endswith('.tdp.out')]
	fl = np.sort(fl)

	fn = [os.path.split(f)[1].split('.')[0]
		for f in fl]
	#
	#-- Obtain the wave bins information with annually exposure hours per wave bin
	#---- use pandas to read in the excel files that defines each wave bin and its associated exposure hours
	#
	colnames = ['time', 'arclen', 'x', 'y', 'z', 'effTen', 'walTen']

	dout = np.array([]).reshape(0,7)
	for f in fl:
		fname = os.path.split(f)[1].split('.')[0]
		df = pd.read_csv(f, names=colnames, sep = '\s')
		ds = df[df.effTen==np.amax(df.effTen)]
		dout = np.vstack([dout, ds])
		print fname + ' is completed'
	dout = pd.DataFrame(dout, columns = colnames, index=fn)
	fout = 'TDP Timehistory MaxTen Summary.csv'
	dout.to_csv(os.path.join(curpath,fout))

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
