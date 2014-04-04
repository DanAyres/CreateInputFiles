#!/usr/bin/env python

from Python_Modules.WriteFiles import *
from Python_Modules.ReadWimsCrossSections import ReadDataFile as WimsData
from Python_Modules.ReadNJOYCrossSections import ReadDataFile as NjoyData
from Python_Modules.ReadNJOYCrossSections import WIMSMacrosopicData, FormatNJOYData

import os

#from Python_Modules.WriteSampyFiles import WriteParameterFiles, Sampy_control, EVENT_template, GemFile, GemEVENT_scripts
#from Python_Modules.Plotting import PlotCorrelationMat, PlotCorrelationData, WriteEigenSpectra
from numpy import set_printoptions, zeros
#import os


if __name__ == '__main__':
 
    set_printoptions(edgeitems=3,infstr='inf', linewidth=150, nanstr='nan', precision=8, suppress=False, threshold=1000) 
    
    wims_file='InputData/FU.radmats'
#    wims_file='InputData/Godiva6.radmats'
    njoy_file='InputData/3group_lib.dat'
    name = 'FLATTOP_U'


    
    # Read the wims interface    
    wims10=True
    WimsXS=WimsData(open(wims_file,'r'), wims10)
    


    NJOY_Extract = ['U235']
    MT_EXTRACT=[18]
#    MT_EXTRACT=[18,452]

    norm=True    # Assume the covariance matrices are realtive - multiply by the wims cross sections
    PCR=False    # Perfectly Correlated Reactions


    # Read the raw NJOY data from file for selected nuclides
    NJOYXS = NjoyData(open(njoy_file,'r'), NJOY_Extract, MT_EXTRACT)
    
    
    # Combine covariances to a form compatible with WIMS
    NJOYXS = FormatNJOYData(NJOYXS, NJOY_Extract)


    # Normalise the covariance data with wims data
    WimsXS=WIMSMacrosopicData(WimsXS, NJOYXS)


    #
    #   Output Files
    #
    if not os.path.exists(name):
        os.makedirs(name)    
    
    GemFile(name, WimsXS)
    
    WriteParameterFiles(name, NJOY_Extract, WimsXS, False)

    EVENT_template(name, WimsXS)    
    
    GemEVENT_scripts(name)

    Sampy_control(name)
#    
#    WriteParameterFiles(WimsXS.NG, filename, NJOY_Extract, WimsXS, MACRO)
#    
#    
#    """
#        Output any extra data: 
#                                    
#            - Correlation matrices
#            - Cumulative Eigenvalues
#    """
#        
#    PlotCorrelationData(WimsXS, NJOY_Extract, folder)
#    WriteEigenSpectra(WimsXS, NJOY_Extract, folder)
##    quit()
#    
#    
