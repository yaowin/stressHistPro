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
		for f in files if f.endswith('.tsj.out')]
	fl = np.sort(fl)

	fn = [os.path.split(f)[1].split('.')[0]
		for f in fl]
	#
	#-- Obtain the wave bins information with annually exposure hours per wave bin
	#---- use pandas to read in the excel files that defines each wave bin and its associated exposure hours
	#
	colnames = ['time', 'effTenRoot', 'walTenRoot', 'shearFRoot', 'bendMomRoot'
				, 'effTenTip', 'walTenTip', 'shearFTip', 'bendMomTip']
	coltogo = ['max wall tension', 'associated shear', 'associated bending moment',
			   'associated wall tension', 'max shear', 'associated bending moment',
			   'associated wall tension', 'associated shear', 'max bending moment']
	dout = np.array([]).reshape(0,9)
	for f in fl:
		fname = os.path.split(f)[1].split('.')[0]
		df = pd.read_csv(f, names=colnames, sep = '\s')
		ds = np.hstack([df[df.walTenRoot==np.amax(df.walTenRoot)].loc[:,['walTenRoot', 'shearFRoot', 'bendMomRoot']],
						 df[df.shearFRoot==np.amax(df.shearFRoot)].loc[:,['walTenRoot', 'shearFRoot', 'bendMomRoot']],
						 df[df.bendMomRoot==np.amax(df.bendMomRoot)].loc[:,['walTenRoot', 'shearFRoot', 'bendMomRoot']]])
		dout = np.vstack([dout, ds])
		print fname + ' is completed'
	dout = pd.DataFrame(dout, columns = coltogo, index=fn)
	fout = 'TSJ Root force linked stats.csv'
	dout.to_csv(os.path.join(curpath,fout))

	dout = np.array([]).reshape(0,9)
	for f in fl:
		fname = os.path.split(f)[1].split('.')[0]
		df = pd.read_csv(f, names=colnames, sep = '\s')
		ds = np.hstack([df[df.walTenTip==np.amax(df.walTenTip)].loc[:,['walTenTip', 'shearFTip', 'bendMomTip']],
						 df[df.shearFTip==np.amax(df.shearFTip)].loc[:,['walTenTip', 'shearFTip', 'bendMomTip']],
						 df[df.bendMomTip==np.amax(df.bendMomTip)].loc[:,['walTenTip', 'shearFTip', 'bendMomTip']]])
		dout = np.vstack([dout, ds])
		print fname + ' is completed'
	dout = pd.DataFrame(dout, columns = coltogo, index=fn)
	fout = 'TSJ Tip force linked stats.csv'
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
