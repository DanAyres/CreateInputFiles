#!/usr/bin/env python

from Python_Modules.ReadWimsCrossSections import ReadDataFile as WimsData
from Python_Modules.ReadNJOYCrossSections import ReadDataFile as NjoyData
from Python_Modules.ReadNJOYCrossSections import WIMSMacrosopicData,NormWithWimsData,FindMissing, FormatNJOYData
from Python_Modules.WriteSampyFiles import WriteParameterFiles, Sampy_control, EVENT_template, GemFile, GemEVENT_scripts
from Python_Modules.Plotting import PlotCorrelationMat, PlotCorrelationData, WriteEigenSpectra
from numpy import set_printoptions, zeros
import os

from Python_Modules.WimsTestData import Sood3Data, GemTestMacro, Sood2

if __name__ == '__main__':
 
    set_printoptions(edgeitems=3,infstr='inf', linewidth=150, nanstr='nan', precision=8, suppress=False, threshold=1000) 
    
#    wims_file='InputData/rowlands_3grp_MOX.radmats'
    wims_file='InputData/FU.radmats'
    njoy_file='InputData/3group_lib.dat'
#    njoy_file='InputData/test.dat'
    name = 'Rowlands_MOX416'


    MACRO= False
        
    WimsXS=WimsData(open(wims_file,'r'), MACRO)

#    WimsXS=Sood2(MACRO, WimsXS)
    
##
#    WimsXS=Sood3Data(MACRO)
##
    
    norm=True    # Assume the covariance matrices are realtive - multiply by the wims cross sections
    PCR=False    # Perfectly Correlated Reactions

    
    NJOY_Extract = ['U235']
#    NJOY_Extract = ['U235', 'U238', 'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242', 'Am241']
#    NJOY_Extract = ['U235', 'U238', 'Pu239', 'Pu240']
#    NJOY_Extract = WimsXS.nuclides
    

    # Read the raw NJOY data from file for selected nuclides
    NJOYXS = NjoyData(open(njoy_file,'r'), NJOY_Extract)
    
    
    # Combine covariances to a form compatible with WIMS
    NJOYXS = FormatNJOYData(NJOYXS, NJOY_Extract)
    

    WimsXS=WIMSMacrosopicData(WimsXS, NJOYXS, norm, MACRO, PCR)
    
    
    """
        Produce files for SAMPY and for the associated deterministic code.
        The following codes are "supported":
            - EVENT
    """
    folder=name+'_'+str(WimsXS.NG)+'grp'
    if MACRO:folder+='Macro'
    else: folder+='Micro'
    if PCR: folder+='PCR'
    print 'Writing Files to " ',folder,' "'
    
    filename = folder+'/'+name

    if not os.path.exists(folder):
        os.makedirs(folder)
    
    
    GemEVENT_scripts(filename)
    
    GemFile(filename, WimsXS, MACRO)
#    GemTestMacro(filename, WimsXS, MACRO)
    
    
    EVENT_template(filename, WimsXS)    
    
    Sampy_control(filename)
    
    WriteParameterFiles(WimsXS.NG, filename, NJOY_Extract, WimsXS, MACRO)
    
    
    """
        Output any extra data: 
                                    
            - Correlation matrices
            - Cumulative Eigenvalues
    """
        
    PlotCorrelationData(WimsXS, NJOY_Extract, folder)
    WriteEigenSpectra(WimsXS, NJOY_Extract, folder)
#    quit()
    
    
