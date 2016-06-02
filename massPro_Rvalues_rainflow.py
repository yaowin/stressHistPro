# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 10:49:27 2016

@author: yyao
"""

import os
import time
import numpy as np
import pandas as pd
import OrcFxAPI as orc
from matplotlib.pyplot import plot, figure, subplots, scatter, show, tight_layout

def dataProcess(s, trace_name = 'test'):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """

    halfcyc = orc.RainflowAssociatedValues(s)   # perform rain flow half count
    cycrange = halfcyc.Range
    meantab = halfcyc.AssociatedValue
    mintab = meantab - cycrange / 2
    maxtab = meantab + cycrange / 2
    badindex = np.absolute(maxtab) < np.absolute(mintab)
    goodindex = np.logical_not(badindex)

    mintab_adj = np.append(mintab[goodindex], maxtab[badindex])
    maxtab_adj = np.append(maxtab[goodindex], mintab[badindex])
    R = mintab_adj / maxtab_adj
    cycrange_adj = np.absolute(maxtab_adj - mintab_adj)
    R_d = (np.tile(np.mean(s),len(cycrange)) - cycrange_adj/2) / (np.tile(np.mean(s),len(cycrange)) + cycrange_adj/2)
    
    A = pd.DataFrame({ 'cycrange' : cycrange_adj,
                         'R' : R})
    B = pd.DataFrame({ 'cycrange' : cycrange_adj,
                         'R_d' : R_d})
    
    A_sort = A.sort(columns = 'cycrange')
    B_sort = B.sort(columns = 'cycrange')
    
    fig, (ax0, ax1) = subplots(ncols=2, figsize=(16, 8))
    ax0.scatter(mintab_adj, maxtab_adj)
    nlim = 150.0    
    t1 = np.arange(-nlim, nlim, 0.1)
    ax0.plot(t1, t1, 'r-')
    ax0.plot(t1, -t1, 'r-')
    ax0.set_xlim(-nlim, nlim)
    ax0.set_ylim(-nlim, nlim)
    ax0.grid('on')
    ax0.set_xlabel('Min value in half cycle (ksi)')
    ax0.set_ylabel('Max value in half cycle (ksi)')
    ax0.set_title('Rainflow Half Cycle Amplitude scatter plot', fontsize=16)
    ax1.scatter(cycrange_adj, R)
    #ax1.plot(B_sort['cycrange'], B_sort['R_d'], c='red')
    ax1.grid('on')
    ax1.set_xlabel('Cycle Range (ksi)')
    ax1.set_ylabel('R value')
    ax1.set_title('Cycle Range vs R scatter plot', fontsize=16)
    tight_layout()
    show()
    fig.savefig("-".join([trace_name, "Adjusted HalfCycle DataScatter.png"]))

    fig, (ax0, ax1) = subplots(ncols=2, figsize=(16, 8))
    # Create a histogram by providing the bin edges (unequally spaced).
    ax0.hist(cycrange, 20, normed=0, histtype='stepfilled', facecolor='g', alpha=0.75)
    ax0.set_title('Cycle Range Historgram stepfilled', fontsize=16)
    ax1.hist(R, 20, normed=0, histtype='stepfilled', facecolor='g', alpha=0.75)
    ax1.set_title('R value Historgram stepfilled', fontsize=16)
    tight_layout()
    show()    
    fig.savefig("-".join([trace_name, "Adjusted Ravlue Histogram.png"]))
    
    A_sort.to_csv("-".join([trace_name, "cycrange Rvalue.csv"]), index=False)

    return

if __name__=="__main__":
    trace = ["Bin001Rep001",]
             
    curpath = "Z:\\SouthSantaCruz\\Wave Fatigue\\Long Term\\TS_SD\\Rep001"			# use your root folder where all sim files are
    
    fl = [os.path.join(dirpath, f)
        for dirpath, dirnames, files in os.walk(curpath)
        for f in files if f[:-4] in trace and f.endswith('.sim')]

    for fn in fl:
        tc = fn[fn.rfind('\\')+1:-4]
        print(": ".join(["opening file", tc]))
        start = time.clock()  
        model = orc.Model(fn, threadCount=8)
        print(": ".join(["file is opened", tc]))
        print 'File reading time: ', int(time.clock()-start), "s\n"	
        LineName = "SSC_Riser"
        simPeriod = 1
        VariableName = "ZZ Stress"
        arcl = 0.0
        ang = 315
        rp = 1
        riserLine = model[LineName]
        timeHist = model.general.TimeHistory("Time", period=simPeriod)
        stressHist = riserLine.TimeHistory(VariableName, period=simPeriod, 
                    objectExtra=orc.oeLine(
                    ArcLength= arcl,
                    Theta=ang,
                    RadialPos=rp)
                    )
        dataProcess(stressHist*0.04788, tc)
        print(" ".join([tc, "is completed"]))