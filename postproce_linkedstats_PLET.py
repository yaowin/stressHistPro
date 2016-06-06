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

    >> python postproce_***.py <input folder> [number of parallel processes]

The parallel process count is optional.
An output folder will be created next to the input folder.
"""
    print errorString

#%%
#-- Subroutine to get the absolute max value of tension
#
def process_file(f):
    fname = os.path.split(f)[1].split('.')[0]
    colnames = ['time', 'Shear_X', 'Shear_Y', 'Tension', 'bendMom_X', 'bendMom_Y', 'Torsion']
    df = pd.read_csv(f, names=colnames, sep = '\s', engine = 'python')
    ten = np.absolute(df.Tension)
    ds = df[ten==np.amax(ten)]
    print fname + ' is completed for ' + str(ds.shape[0])
    if ds.shape[0] >1:
        ds = ds[0:1]
        print 'this repitition is truncated'
    return ds
    

#%%
#-- Obtain the stress history folder list
#---- this'll walk your directories, subfolders recursively and return all folders that contains files.
#
def ProcessCases(curpath):

	fl = [os.path.join(dirpath, f)
		for dirpath, dirnames, files in os.walk(curpath)
		for f in files if f.endswith('.plet.out')]
	fl = np.sort(fl)

	fn = [os.path.split(f)[1].split('.')[0]
		for f in fl]
	#
	#-- Obtain the wave bins information with annually exposure hours per wave bin
	#---- use pandas to read in the excel files that defines each wave bin and its associated exposure hours
	#

	coltogo = ['associated Shear_X', 'associated Shear_Y', 'Max Tension', 'associated bendMom_X'
				, 'associated bendMom_Y', 'associated Torsion']
	frames = [process_file(f) for f in fl]
	dout = pd.concat(frames)
	dout = dout.drop('time',axis =1)
	dout.columns = coltogo
	dout.index = fn
	fout = 'PLET_Anchor Timehistory MaxTen Summary.csv'
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
