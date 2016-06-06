# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:53:39 2016

@author: Y.Yao
All rights reserved, any usage in public domain is prohibited.
"""

import pandas as pd
import xlsxwriter as xw
import xlrd as xr
import os
import time
import OrcFxAPI as orc
import numpy as np

#%%
def get_ParameterInput(input_file = 'test_ECA StressHistPro Input Driver_new.xlsx'):
    """
    Returns all input paramters in one data structure
    
    function para = get_ParameterInput(input_file)
    % GET_PARAMETERINPUT Fetch and group all input parameters from input xls
    %        PARA = GET_PARAMETERINPUT(INPUT_FILE) input the parameters
    %        defined in the spreadsheet of a pre-set format and tabulated
    %        data. PARA consists of a dictionary data structure with multiple
    %        fields that are useful to the stress history post-processing.
    %      
    %        All stress units in raw data extraction from OrcaFlex is meant to
    %        be transferred to MPa.
    %        
    
    % Yao Yao, 5.16.16 (copyrighted OFI).
    % This function is restricted to the public; Any use requires permit.
    
    """
    book = xr.open_workbook(input_file)
    sheet = book.sheet_by_index(0)
    
    data = {}
    parasetting = {}
    
    with pd.ExcelFile(input_file) as xls:
        data['InputSheet'] = pd.read_excel(xls, 'InputSheet', header = 13, na_values=['NA'])
        parasetting['WaveBins'] = pd.read_excel(xls, 'WaveBins', header = 7, parse_cols = "B:D", index_col= 0)
        parasetting['HistBins'] = pd.read_excel(xls, 'WaveBins', parse_cols = "F", index_col= None)
    
    mat_Elist = {'Titanium': 114E3, 'Steel': 207E3} # Unit - MPa!
    rpLoc_list = {'Inner': 0.0, 'Midwall': 0.5, 'Outer': 1.0}
    simp_list = {'Stage 1': 1, 'Stage 2': 2, 'Latest Wave': 'pnLatestWave', 
                 'Whole Simulation': 'pnWholeSimulation', 'Static State': 'pnStaticState'}
    unitcon_list = {'kPa -> Mpa': 1.0, 'psf -> Mpa': 3.281 , 'psi -> Mpa': 39.37, 
                 'none': 1.0} 
                 # Note this unit conversion is at curvature data, so only the
                 # length unit is convert to meter, for MPa results
    outdet_list = {'Minimal': 0, 'Moderate': 1, 'Full Log Detail': 2}
    
    parasetting['path'] = sheet.cell_value(rowx=5, colx=3)
    parasetting['outputfn'] = sheet.cell_value(rowx=6, colx=3)
    parasetting['linename'] = sheet.cell_value(rowx=7, colx=3)
    parasetting['temperature'] = sheet.cell_value(rowx=8, colx=3)
    parasetting['tothrs'] = sheet.cell_value(rowx=10, colx=5)    
    parasetting['unitstr'] = sheet.cell_value(rowx=10, colx=6)
    parasetting['arclength'] = data['InputSheet'].arclen.dropna()
    parasetting['outerd'] = data['InputSheet'].od.dropna()
    parasetting['innerd'] = data['InputSheet'].id.dropna()
    parasetting['thick'] = (parasetting['outerd'] - parasetting['innerd']) / 2.0
    parasetting['matstr'] = data['InputSheet'].mat.dropna()
    parasetting['matval'] = [mat_Elist[matstr] for matstr in parasetting['matstr']]
    parasetting['snstr'] = data['InputSheet'].sn.dropna()
    parasetting['residual'] = data['InputSheet'].res.dropna()    
    parasetting['SCF'] = data['InputSheet'].scf.dropna()
    parasetting['angular'] = data['InputSheet'].ang.dropna()
    parasetting['radstr'] = data['InputSheet'].rad.dropna()
    parasetting['radval'] = [rpLoc_list[rdstr] for rdstr in parasetting['radstr']]
    simperiod = data['InputSheet'].simp.iloc[0]
    unitcon = data['InputSheet'].unitcon.iloc[0]
    outdstr = data['InputSheet'].outdetail.iloc[0]
    
    parasetting['simPeriod'] = orc.Period(simp_list[simperiod])   # use the stage 1 period of simulation
                # pnSpecifiedPeriod, pnLatestWave, pnWholeSimulation, pnStaticState
    parasetting['unit'] = unitcon_list[unitcon]    # convert the model ouput unit to Mega Pa
    parasetting['outlev'] = outdet_list[outdstr]    # control the output details level

    return parasetting

#%%
def save_StressHistories(parasetting, processCount):
#------ Main body of stage 1 process - create stress history files *.CSV
#       the stress history data will be in MPa unit

#- 1.0 Fetch the file list
#
    fl = [os.path.join(dirpath, f)
        for dirpath, dirnames, files in os.walk(parasetting['path'])
        for f in files if f.endswith('.sim')]
    
    #fn = [f[f.rfind('\\')+1:]
    #    for f in fl]
    
    print 'Total sim files number in ' + parasetting['path'] + ' is ' + str(len(fl))
    
    ind = 0
    maxhr = 0.0
    max_s = 0.0
    min_s = 0.0

#- 2.0 Go through the sim files and extract curvature data
#
    print 'data extraction begins...'	
    for filename in (fl):
        ind +=1
        start = time.clock()
        # ---ATTENTION--- the waveid convention in full file name
        waveid = filename[filename.rfind('\\')+1:-4]
        model = orc.Model(filename, threadCount=processCount)
        riserLine = model[parasetting['linename']]
        timeHist = model.general.TimeHistory("Time", period=parasetting['simPeriod'])
        timeHist = np.array([timeHist])
        maxhr = max(maxhr, np.max(timeHist)/3600.0)
        
        for arcl in (parasetting['arclength']):
            # the arc length is converted to the model unit, used for data extraction only
            arc = arcl * parasetting['unit']
            # the index of current processing arc length
            arcind = np.nonzero(parasetting['arclength'] == arcl)[0][0]
            # the tensile strain extraction (unitless)
            Eps = riserLine.TimeHistory("Direct Tensile Strain", period=parasetting['simPeriod'], 
                                        objectExtra=orc.oeArcLength(arc)
                                        )
            # the x-curvature (convert unit from model to 1/m)
            Cx = riserLine.TimeHistory("x-Curvature", period=parasetting['simPeriod'], 
                                        objectExtra=orc.oeArcLength(arc)
                                        ) * parasetting['unit']
            # the y-curvature (convert unit from model to 1/m)
            Cy = riserLine.TimeHistory("y-Curvature", period=parasetting['simPeriod'], 
                                        objectExtra=orc.oeArcLength(arc)
                                        ) * parasetting['unit']
            # the resultant-curvature (computed for all angular positions)
            Cr = np.outer(Cx, np.sin(np.radians(parasetting['angular']))) - np.outer(Cy, np.cos(np.radians(parasetting['angular'])))
            # Resultant curverature based on theta and two x- and y- component, -->
    		# reference OrcaFlex Help: Tension and Curvature Stress Factors (stress factors only) <--
            for rp in (parasetting['radval']):
                rpind = parasetting['radval'].index(rp)
                TenS = np.dot(parasetting['matval'][arcind], Eps/100) 
    			#Direct Tensile stress = tensile strain times modulus, the tensile strain was in percentage
                BenS = np.dot(Cr, parasetting['matval'][arcind]*(parasetting['innerd'][arcind+1]/2 + rp*parasetting['thick'][arcind+1]))
    			#Bending stress = resultant curvarture times modulus times distance from neutrual axis to stress fibre
                stressHist = np.transpose(np.tile(TenS, (np.size(parasetting['angular']), 1))) + BenS
    			#ZZ stress = Tensile Stress + Bending Stress
                Histories = np.hstack((timeHist.T, stressHist))
                # Append to the time history column, end up with 9 columns data matrix
                maxs = np.amax(stressHist)
                mins = np.amin(stressHist)
                max_s = max(max_s, maxs)
                min_s = min(min_s, mins)
                subfoldpath = parasetting['path'] +'\\Stress Histograms' + '\\arcLen ' + str(int(np.rint(arcl))) + 'm' + parasetting['radstr'][rpind+1]
                if not os.path.exists(subfoldpath):
                    os.makedirs(subfoldpath)
                textFstr = os.path.join(subfoldpath, waveid + '.csv')

#- 3.0 Save the stress histories data in text file *.csv
#                
                np.savetxt(textFstr, Histories, delimiter=",")
    	print 'Wavebin processed: {0:3d} out of {1:3d}  <==> Time elapsed: {2:2d} min; Time to go: {3:4d} min'.format(ind, len(fl), int((time.clock()-start)/60), int((time.clock()-start)*(len(fl)-ind)/60))

    return [maxhr, max_s, min_s]
#
#----[FINAL]------- Obtain xx arcLen X inner/midwall/outer = xxx total subfolders under '/Stress Histograms/,
#------------------ each subfolders will contain xxx (wave bins) csv files
#------------------ each csv file will contain 1 time trace, and 8 (angular positions) columns of stress history data.
#

#%% Calculate the Titanium Fatigue
def calc_TiFatigue(cyc, snstr, t, scf, rt = False):
    
    unit_c = 0.14504   # convert MPa to ksi
    
    try:    
        if snstr == 'RTI R-ratio':
            sig_f = 0.00008356*t**2 - 0.1337*t + 182.59
            b =(-1)*(np.log10(sig_f)-1.4971)/7.301
            n = 1/b
            A = (0.5**n)*(1/(sig_f*(2**b)))**n
            R_wavg = np.sum(cyc['R_xcyc'])/np.sum(cyc['cyc'])
            ## calculate cycles to failure based on individual R ratio
            if rt:
                N = A*(np.sqrt(2/(1-cyc['R']))*scf*cyc['cyc']*unit_c)**n
            else:
                N = A*(scf*cyc['cyc']*unit_c)**n
            ## calculate cycles to failure based on weighted average R ratio
            #N = A*(np.sqrt(2/(1-R_wavg))*scf*cyc['cyc'])**n
        if snstr == 'RTI 1G':
            N = 6.79986*(10**19) *  (scf*cyc['cyc'])**(-6)
    except snstr not in ['RTI R-ratio', 'RTI 1G']:
        print("Unexpected S-N curve name error")

    halfcount_dam = 1 / N
    dam = halfcount_dam / 2    

    return np.sum(dam)
    
#%% Calculate the Titanium Fatigue
def calc_SteelFatigue(cyc, snstr, t, scf):
    
    try:    
        if snstr == 'BS 7608C':    # Stress Joint Forging
            ## calculate cycles to failure based on individual cyc range
            t_b = 0.016 # ref thickness is 16mm
            S = cyc['cyc']*scf*(np.max([t, t_b])/t_b)**0.25
            N = 4.22669*(10**13) * S**(-3.5)
        if snstr == 'DNV D':    #Pipe Girth Weld[DnV D, CP Protected in Water]
            t_b = 0.032 # ref thickness is 32mm
            S = cyc['cyc']*scf*(np.max([t, t_b])/t_b)**0.1
            N = 5.80764*(10**11) * S**(-3.0)
    except snstr not in ['BS 7608C', 'DNV D']:
        print("Unexpected S-N curve name error")

    halfcount_dam = 1 / N
    dam = halfcount_dam / 2    

    return np.sum(dam)
    
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
        nlim = 100.00
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
        fig.savefig("-".join([trace_name, "Adjusted Ravlue Histogram.png"]))

        X, Y = np.meshgrid(rBins, sBins)
        im = imshow(hist,interpolation = 'nearest',cmap='cool')
        colorbar(im, orientation='horizontal')
        grid(True)
        fig.savefig("-".join([trace_name, "R-cyc 2Dhistogram colormap.png"]))
        
    else:
        A.to_csv("-".join([trace_name, "cycrange Rvalue.csv"]), index=False)
        hist.to_csv("-".join([trace_name, "2dhistogram.csv"]), index=False)

    return

#%% Rainflow Half Counts
def get_Rainflow(s, res_s=0, thresh_s=68.948):
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
#%% running average routine
def runningMean(x, N):
    y = np.zeros((len(x)-1,))
    for ctr in range(len(x)-1):
         y[ctr] = np.sum(x[ctr:(ctr+N)])
    return y/N
#%% running average routine
def Usage():
    errorString = """
The expected command line for this script has the form:

    >> python StressHistPro.py <input excel file> [number of stage processes]

The stage process number is optional, default = 3

Stage 0: fetch input parameters from input excel
Stage 1: Stress History Extraction and Storage to *.CSV file
Stage 2: Run rainflow halfcount, stress histogram, and fatigue damage calc.
[number of stage processes] is the sum of your desired stage numbers to run. e.g. 3 is running all stages

All input parameters are controled in the input excel file.
All output folders will be created inside the working folder.

"""
    print errorString
	
#%% Main Program
if __name__=="__main__":
    import sys
    import matplotlib.pyplot as plt
    # We've come in from a command line, running this script as a program.
    try:
        inputFile = sys.argv[1]
        try:
            processStage = int(sys.argv[2])
        except IndexError:
            processStage = 3
    except:
        Usage()
        exit()
    
    
    print '*** Stage 0: Input Parameters reading...'	
    para = get_ParameterInput(inputFile)
    print 'Completed.'
    
    if processStage in [1, 3]:
        print '*** Stage 1: Stress History Calculation and Storage begins...'	
        [simhr, maxs, mins] = save_StressHistories(para, 8)
        print 'All Stress Time History Data: hours in 1 simulation = {0:.1f}, max stress = {1:.1f} MPa, min stress = {2:.1f} MPa'.format(simhr, maxs, mins)
        if processStage == 1:
            exit()
    
    print '*** Stage 2: Stress Histogram Generation begins...'	

    debug_dir ="\\".join((para['path'], 'Stress Histograms', 'debug log'))
    if not os.path.exists(debug_dir):
        os.makedirs(debug_dir)    
    sBins = para['HistBins']['Histogram Bins (Mpa)']
    #sBins = np.linspace(0, np.ceil(maxs), 21)
    rBins = np.linspace(-1.0, 1.0, 11)
    
    workbook = xw.Workbook(os.path.join(para['path'],para['outputfn']))
    bold = workbook.add_format({'bold': True})
    boldred = workbook.add_format({'bold': True, 'font_color': 'red'})
    scifm = workbook.add_format({'num_format': 0x0b})
    redscifm = workbook.add_format({'num_format': 0x0b, 'font_color': 'red'})
    intfm = workbook.add_format({'num_format': '#,##0'})
    digfm = workbook.add_format({'num_format': '0.00'})
    
    for arcl in (para['arclength']):    #step through all acrlength position
        arcind = np.nonzero(para['arclength'] == arcl)[0][0]
        start = time.clock()
        i_col = 4
        worksheet = workbook.add_worksheet(" ".join(('ArcLen=',str(int(np.rint(arcl))), 'm')))
        worksheet.set_column('A:D', 15)
        worksheet.write('D1', 'from hang-off', bold)
        worksheet.write('D2', 'Angular Position', bold)
        worksheet.write('A3', 'Bin number:', bold)
        worksheet.write('B3', 'Min Bin Stress (MPa)', bold)
        worksheet.write('C3', 'Max Bin Stress (MPa)', bold)
        worksheet.write('D3', 'Bin Size (MPa)', bold)                

        for rp in (para['radval']):    #step through radial position
            rpind = para['radval'].index(rp)
            cdir = "\\".join((para['path'], 'Stress Histograms', 'arcLen ' + str(int(np.rint(arcl))) + 'm' + para['radstr'][rpind+1]))
            
            fl = [os.path.join(dirpath, f)
                for dirpath, dirnames, files in os.walk(cdir)
                for f in files if f.endswith('.csv')]

            worksheet.write(0, i_col, cdir[cdir.rfind('\\')+1:], bold)
            for ang in (para['angular']): #scan through all 8 columns in the array
                col = np.nonzero(para['angular'] == ang)[0][0]
                worksheet.write_number(1, i_col, ang, intfm)
                worksheet.write(2, i_col, 'Annual Counts', bold)        
                worksheet.write(2, i_col+1, 'R ratio', bold) 
                cycleCounts = np.array([]).reshape(0,len(sBins)-1)
                cycleR = np.array([]).reshape(0,len(sBins)-1)
                totDamage = 0.0
                for filename in (fl):
                    #----- ATTENTION ----                    
                    waveid = filename[filename.rfind('\\')+1:-4]
                    #-----check the wave load case name
                    bin_hrs = para['WaveBins']['Annual Exposure Hours'].loc[waveid]
                    r_toggle = para['WaveBins']['R-ratio Toggle'].loc[waveid]
                    # read the stress history data from each cvs file
                    history_data = np.genfromtxt(filename, delimiter=',')
                    # read and convert stress history into MPa, skip the time stamp column
                    #stressHistory = history_data[:, col + 1] * para['unit']
                    stressHistory = history_data[:, col + 1]
                    # load stress time trace into rainflow half counts 
                    cycrange = get_Rainflow(stressHistory, para['residual'][arcind+1])
                    #cycrange = get_Rainflow(stressHistory)

                    # calculate the fatigue damage
                    if para['matstr'][arcind+1] == 'Titanium':
                        damage = calc_TiFatigue(cycrange, para['snstr'][arcind+1], para['temperature'], para['SCF'][arcind+1], r_toggle)
                        probDamage = damage / simhr * (bin_hrs)
                        totDamage =  totDamage + probDamage
                    if para['matstr'][arcind+1] == 'Steel':
                        thick = (para['outerd'][arcind+1] - para['innerd'][arcind+1])/2
                        damage = calc_SteelFatigue(cycrange, para['snstr'][arcind+1], thick, para['SCF'][arcind+1])
                        probDamage = damage / simhr * (bin_hrs)
                        totDamage =  totDamage + probDamage                    
                    # histogram the stress half cycles
                    hist, xedges, yedges = np.histogram2d(cycrange['R'], cycrange['cyc'], (rBins, sBins))
                    rx = runningMean(xedges, 2)
                    ry = runningMean(yedges, 2)
                    ra = np.tile(rx,len(sBins)-1).reshape((len(sBins)-1,len(rBins)-1)).T
                    if para['outlev'] > 1:
                        hist_out = pd.DataFrame(hist.T, index = ry, columns = rx)
                        hist_out_str = "_".join(('Histgram arcl='+str(int(arcl)), 'ang='+str(int(ang)), 'rad='+para['radstr'][rpind+1], waveid))
                        hist_out.to_csv(os.path.join(debug_dir, hist_out_str+'.csv'))
                        cyc_out_str = "_".join(('CycRange_R arcl='+str(int(arcl)), 'ang='+str(int(ang)), 'rad='+para['radstr'][rpind+1], waveid))  
                        cycrange.to_csv(os.path.join(debug_dir, cyc_out_str+'.csv'))
                    
					#-- get the array of R ratio bin middle value in shape

                    # sum the number counts of stress histogram
                    counts = np.sum(hist,axis=0)                    
                    # calculate the count weighted average of R ratio for all delta_s
                    hist_R = np.sum(np.multiply(ra, hist),axis=0) / np.sum(hist,axis=0)
                    # for zero count, zero stress range, R = 1.0 by definition 
                    hist_R[np.isnan(hist_R)] = 1.0   
                    # exposure hour weighted histogram for each wave bin
                    probCounts = counts / simhr * (bin_hrs)
                    wR = hist_R * bin_hrs / para['tothrs']
                    # stack up across all wave load cases
                    cycleCounts = np.vstack((cycleCounts, probCounts))
                    cycleR = np.vstack((cycleR, wR))
                    # call the output routine for detail level control
                    # ---TODO--- get full position info into the output
                    #output_Details(cycrange, waveid, para['outlev'])
                    
                sumCounts = np.sum(cycleCounts,axis=0)
                wavgR = np.sum(cycleR, axis=0)
                    #------------------ print data column headers and statistical results
                    #------------------ print binned histogram results
                i_row = 3
                for bind in range(len(sumCounts)):
                    bmin = sBins[bind]
                    bmax = sBins[bind+1]
                    bmean = np.mean([bmin, bmax])
                    bSize = bmax - bmin
                    bN = sumCounts[bind]
                    bR = wavgR[bind]

                    worksheet.write_number(i_row, 0,     bind+1,  intfm )
                    worksheet.write_number(i_row, 1,     bmin,   digfm )
                    worksheet.write_number(i_row, 2,     bmax,   digfm )
                    worksheet.write_number(i_row, 3,     bSize,   digfm )
                    worksheet.write_number(i_row, i_col,     bN,    scifm )
                    worksheet.write_number(i_row, i_col + 1,     bR,    digfm )
                    i_row += 1
                worksheet.write(i_row, 3,     'Fatigue Damage',    bold )
                worksheet.write_number(i_row, i_col,     totDamage,    redscifm )
                i_col +=2
                print 'Angular Position: {0:1d} -- Radial Position: {1:.1f} -- Arclength: {2:.2f}'.format(int(ang), rp, arcl)
        print '*** Complete Sweeping Arclength {0:.2f} <==> Time elapsed in this sweep: {1:3d} min'.format(arcl, int((time.clock()-start)/60))
    workbook.close()
