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
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, figure, subplots, scatter, show, tight_layout
from matplotlib.pyplot import imshow, colorbar, grid

#%% Calculate the Titanium Fatigue
def calc_TiFatigue(cyc, t, scf):
    
    sig_f = 0.00008356*t**2 - 0.1337*t + 182.59
    b =(-1)*(np.log10(sig_f)-1.4971)/7.301
    n = 1/b
    A = (0.5**n)*(1/(sig_f*(2**b)))**n
    R_wavg = np.sum(cyc['R_xcyc'])/np.sum(cyc['cyc'])
    ## calculate cycles to failure based on individual R ratio    
    N = A*(np.sqrt(2/(1-cyc['R']))*scf*cyc['cyc'])**n
    ## calculate cycles to failure based on weighted average R ratio
    #N = A*(np.sqrt(2/(1-R_wavg))*scf*cyc['cyc'])**n

    halfcount_dam = 1/N
    dam = halfcount_dam / 2    

    return np.sum(dam)
#%% Rainflow Half Counts
def get_Rainflow(s, res_s=0, thresh_s=10):
    # perform rain flow half count using associated values
    #-----THIS MUST BASE ON LASTEST ORCAFLEX API
    # default threshhold for the stress cycle range is 10 MPa
    # the delta sigma prime will only be factored by 1/sqrt(1-R) if delta sigma
    # is above the threshhold.
    
    halfcyc = orc.RainflowAssociatedValues(s)   # perform rain flow half count
    cycrange = halfcyc.Range
    meantab = halfcyc.AssociatedValue + np.tile(res_s, len(cycrange))
    mintab = meantab - cycrange / 2
    maxtab = meantab + cycrange / 2
    badindex = np.absolute(maxtab) < np.absolute(mintab)
    goodindex = np.logical_not(badindex)

    mintab_adj = np.append(mintab[goodindex], maxtab[badindex])
    maxtab_adj = np.append(maxtab[goodindex], mintab[badindex])
    R = mintab_adj / maxtab_adj
    cycrange_adj = np.absolute(maxtab_adj - mintab_adj)
    R_d = (np.tile(np.mean(s),len(cycrange)) - cycrange_adj/2) / (np.tile(np.mean(s),len(cycrange)) + cycrange_adj/2)
    
    cycrange_prime = np.tile(0.0, len(cycrange))
    threshindex = cycrange_adj >=thresh_s
    belowindex = np.logical_not(threshindex)
    cycrange_prime[threshindex] = cycrange_adj[threshindex] / np.sqrt(1.0 - R[threshindex])
    cycrange_prime[belowindex] = cycrange_adj[belowindex]

    A = pd.DataFrame({ 'cyc' : cycrange_adj,
                       'cyc_p' : cycrange_prime,
                       'minstress' :mintab_adj,
                       'maxstress' :maxtab_adj,
                         'R' : R,
                         'R_d': R_d,
                         'R_xcyc': np.multiply(cycrange_adj, R)})

    A_sort = A.sort_values(['cyc'])

    return A_sort
#%% Calculate the Titanium Fatigue
def output_Details(A, trace_name = 'test', outlev = 0):
    from matplotlib.pyplot import plot, figure, subplots, scatter, show, tight_layout
    from matplotlib.pyplot import imshow, colorbar, grid

    mintab = A['minstress']
    maxtab = A['maxstress']
    cycrange = A['cyc']
    R_d = A['R_d']
    R = A['R']
    
    if outlev < 1:
        return
    if outlev < 2:
        fig, (ax0, ax1) = subplots(ncols=2, figsize=(16, 8))
        ax0.scatter(mintab, maxtab)
        nlim = 100.0
        t1 = np.arange(-nlim, nlim, 0.1)
        ax0.plot(t1, t1, 'r-')
        ax0.plot(t1, -t1, 'r-')
        ax0.set_xlim(-nlim, nlim)
        ax0.set_ylim(-nlim, nlim)
        ax0.grid('on')
        ax0.set_xlabel('Min value in half cycle (MPa)')
        ax0.set_ylabel('Max value in half cycle (MPa)')
        ax0.set_title('Rainflow Half Cycle Amplitude scatter plot', fontsize=16)
        ax1.scatter(cycrange, R)
        ax1.plot(cycrange, R_d, c='red')
        ax1.grid('on')
        ax1.set_xlabel('Cycle Range (MPa)')
        ax1.set_ylabel('R value')
        ax1.set_title('Cycle Range vs R scatter plot', fontsize=16)
        tight_layout()
        fig.savefig("-".join([trace_name, "Adjusted HalfCycle DataScatter.png"]))
    
        fig, (ax0, ax1) = subplots(ncols=2, figsize=(16, 8))
        # Create a histogram by providing the bin edges (unequally spaced).
        ax0.hist(cycrange, 20, normed=0, histtype='stepfilled', facecolor='g', alpha=0.75)
        ax0.set_title('Cycle Range Historgram stepfilled', fontsize=16)
        ax1.hist(R, 20, normed=0, histtype='stepfilled', facecolor='g', alpha=0.75)
        ax1.set_title('R value Historgram stepfilled', fontsize=16)
        tight_layout()
        fig.savefig("-".join([trace_name, "Histogram.png"]))

        X, Y = np.meshgrid(rBins, sBins)
        im = imshow(hist,interpolation = 'nearest',cmap='cool')
        colorbar(im, orientation='horizontal')
        grid(True)
        fig.savefig("-".join([trace_name, "R-cyc 2Dhistogram colormap.png"]))
        
    else:
        A.to_csv("-".join([trace_name, "cycrange Rvalue.csv"]), index=False)
        hist.to_csv("-".join([trace_name, "2dhistogram.csv"]), index=False)

    return

#%% running average routine
def runningMean(x, N):
    y = np.zeros((len(x)-1,))
    for ctr in range(len(x)-1):
         y[ctr] = np.sum(x[ctr:(ctr+N)])
    return y/N
#%%
#-- Obtain the stress history csv files folder list
#---- this'll walk your directories, subfolders recursively and return all folders that contains files.
#

#curpath = os.path.join(os.getcwd(), 'Stress Histograms')		# use the current folder as where this code is
curpath = "Y:\\OJ Wave Fatigue_1311mWD\\Stress History"			# use your root folder where all csv files are
outpath = "Y:\\OJ Wave Fatigue_1311mWD\\Stress Histogram"			# spec your root folder where all output files are

dirlist = [os.path.join(curpath, 'arcLen ' + str(arcnum) +'m_' + iopos_str)
    for arcnum in [0, 3, 5]
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
angle = [0, 45, 90, 135, 180, 225, 270, 315]               

s1 = np.append([0,1.0], np.linspace(5, 140, 28))
s2 = np.append(np.linspace(155, 275,9), np.linspace(310,1010,21)) # Note the unit is MPa
s3 = np.append(s1, s2)
sBins = np.append(s3,  np.linspace(1080, 1500, 7))
rBins = np.linspace(-1.0, 1.0, 11)

unit_c = 1    # stress history is in 1 MPa = 1 MPa
simhr = 3 # each simulation file has 3 hours
#
#---------------- Specify the output file name for Stress Histogram, Binned/Original rainflow half count distribution
#
OutputStr = 'Histogram Binned Annual Probability Weighed.xlsx'

#%%
#---------------- 1. Extract zz stress time history by column (angular position) for each wave bin (load case).
#		  2. Run rain flow half counts at the time history.
#		  3. Get histogram of the sorted Half Cycle Counts and bin the results.
#		  4. Output the binned results and statistics into Excel file, each angular position is an added worksheet.
#		  5. Output the annual occurence hours weighted counts

print '*** histogram generation begins...'	

for cdir in dirlist:
    fl = [os.path.join(dirpath, f)
        for dirpath, dirnames, files in os.walk(cdir)
        for f in files if f.endswith('.csv')]
   
    for col in range(len(angle)):    #scan through all 8 columns in the array
        start = time.clock()
        rawR = np.array([]).reshape(0,0)
        rawCyc = np.array([]).reshape(0,0)
        rawHist = np.zeros((len(rBins)-1, len(sBins)-1))
        outdir = os.path.join(outpath, cdir[cdir.rfind('\\')+1:], 'Angular_' + str(angle[col]))
        if not os.path.exists(outdir):
            os.makedirs(outdir)    
        totDamage = 0.0
        
        ind_f = 1
        for filename in (fl):
            waveid = int(filename[filename.rfind('_')+1:-4])
            bin_hrs = Wavebin_hrs[waveid]
            if np.count_nonzero(bin_hrs) > 0:
                history_data = np.genfromtxt(filename, delimiter=',')
                # read the stress history data from each cvs file
                stressHistory = history_data[:, col + 1] * unit_c    # read and convert stress history into ksi, skip the time stamp column
                cycrange = get_Rainflow(stressHistory, 68.948)
    
                #damage = calc_TiFatigue(cycrange, 138, 1.0)
                #probDamage = damage / simhr * (bin_hrs)
                #totDamage =  totDamage + probDamage
                rawR = np.append(rawR, cycrange['R'])
                rawCyc = np.append(rawCyc, cycrange['cyc'])
                # run the 2D historgram on R ratio and cyc range
                hist, xedges, yedges = np.histogram2d(cycrange['R'], cycrange['cyc'], (rBins, sBins))
                # exposure hour weighted histogram for each wave bin
                rawHist = rawHist + hist / simhr * (bin_hrs)

                rx = runningMean(xedges, 2)
                ry = runningMean(yedges, 2)
                hist_out = pd.DataFrame(hist.T, index = ry, columns = rx)
                hist_out_str = "_".join(('Histogram waveid', str(waveid)))
                hist_out.to_csv(os.path.join(outdir, hist_out_str+'.csv'))
                cyc_out_str = "_".join(('CycRange_R waveid', str(waveid)))
                cycrange.to_csv(os.path.join(outdir, cyc_out_str+'.csv'))
            ind_f += 1
            # call the output routine for detail level control
            # ---TODO--- get full position info into the output
            #output_Details(cycrange, waveid, 1)
            print ind_f, filename[filename.rfind('wavebin'):]
        hist2d, xedges, yedges = np.histogram2d(rawR, rawCyc, (rBins, sBins))
        rx = runningMean(xedges, 2)
        ry = runningMean(yedges, 2)
        hist_out = pd.DataFrame(hist2d.T, index = ry, columns = rx)
        hist_out_str = 'Histogram Allwavecases'
        hist_out.to_csv(os.path.join(outdir, hist_out_str+'.csv'))
        
        print 'plotting starts...'
        fig1, (ax0, ax1) = subplots(ncols=2, figsize=(16, 8))
        ax0.scatter(rawCyc[rawCyc>10], rawR[rawCyc>10])
        ax0.set_ylim(-1.0, 1.0)        
        ax0.grid('on')
        ax0.set_xlabel('Cycle Range (MPa)')
        ax0.set_ylabel('R ratio')
        ax0.set_title('Cycle Range vs R scatter plot', fontsize=16)
        
        ax1.hexbin(rawCyc[rawCyc>10], rawR[rawCyc>10], bins='log', cmap='cool')
        ax1.set_ylim(-1.0, 1.0)
        ax1.grid('on')
        ax1.set_title('Cycle Range vs R hexbin plot', fontsize=16)
        ax1.set_xlabel('Cycle Range (MPa)')
        tight_layout()
        trace_name = cdir[cdir.rfind('\\')+1:] + str(angle[col]) + 'deg'
        fig1.savefig(os.path.join(outdir, "-".join([trace_name, "DataScatter.png"])))
       
        ra = np.tile(rx,len(sBins)-1).reshape((len(sBins)-1,len(rBins)-1)).T
        # sum the number counts of stress histogram
        counts = np.sum(rawHist,axis=0)                    
        # calculate the count weighted average of R ratio for all delta_s
        hist_R = np.sum(np.multiply(ra, rawHist),axis=0) / np.sum(rawHist,axis=0)
        # for zero count, zero stress range, R = 1.0 by definition 
        #hist_R[np.isnan(hist_R)] = 0.0
        
        fig2, (ax0, ax1) = subplots(ncols=2, figsize=(16, 8))
        ax0_waf = ax0.twinx()
        ax0.hist(rawCyc[rawCyc>10], 50, normed=0, histtype='stepfilled', facecolor='g', alpha=0.75)
        ax0_waf.plot(ry,hist_R, 'r', lw=1.0)
        ax0_waf.set_ylim(-1.0, 1.0)
        ax0_waf.yaxis.grid(True)
        ax0.set_title('Cycle Range Histogram with weighted mean R', fontsize=16)
        ax0.set_xlabel('Cycle Range (MPa)')
        all_index = [np.all([rawCyc>yedges[y], rawCyc<yedges[y+1]],axis=0) for y in range(len(ry))]
        all_Rdata = [rawR[index] for index in all_index]
        ax1.boxplot(all_Rdata,
                         vert=True,   # vertical box aligmnent
                         patch_artist= True, # fill with color
                         showfliers = False)   # don't show fliers
        plt.setp(ax1, xticks=[y+10 for y in np.arange(0, len(all_Rdata), 10)],
         xticklabels=['100', '200', '300', '400', '500', '600'])
        ax1.plot(hist_R,'r', lw=1.0)
        ax1.set_ylim(-1.0, 1.0)
        ax1.grid('on')
        ax1.set_title('Cycle Range vs R distribution box plot', fontsize=16)
        ax1.set_xlabel('Cycle Range (MPa)')
        fig2.savefig(os.path.join(outdir, "-".join([trace_name, "2DHistogram Hexbin.png"])))
        
    print '*** Complete Sweeping Arclength <==> Time elapsed in this sweep: {0:3d} min'.format(int((time.clock()-start)/60))
