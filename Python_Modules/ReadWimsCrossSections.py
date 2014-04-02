"""
    This module reads a WIMS cross section file and returns
    a list containing either the macroscopic or microscopic XS
    
    
    Author:
            Dan Ayres
            
    
"""

from numpy import zeros


class Data:

    """
        A class that contains all of the necessary cross section
        data from the WIMS file.
    """

    def __init__(self, NG, SPCTRM, MicroList, MacroList, nuclides, matspec, MACRO):
        self.NG=NG
        self.SPECTRUM=SPCTRM
        self.Micro=MicroList
        self.Macro=MacroList
        self.nuclides=nuclides
        self.MatSpecs=matspec
        
        if MACRO:        
            self.CrossSections=MacroList
        else:
            self.CrossSections=MicroList
        

class XS:
    """
        This is a class for the cross section data.
    """
    
    def __init__(self, name, nuclide,num_den,mats, SIGABS, SIGTRAN, scatter, SIGFISS,  FISSNU, FISSCHI):
        self.name=name              # A reference name of your choosing
        self.nuclide=nuclide        # Nuclide (s)
        self.ND=num_den             # Number density (s)
        self.mats=mats
        self.SIGABS=SIGABS          # Absorption
        self.SIGTRAN=SIGTRAN        # Transport
        self.scatter=scatter        # g'-> g scatter matrix - 0 and 1 (if given) moments 
        self.SIGFISS=SIGFISS        # Fission
        self.FISSNU=FISSNU          # Average number of fission neutrons
        self.FISSCHI=FISSCHI        # Fission spectrum
        
        # Calculate the total XS
        self.SIGTOT = [SIGABS[i]+sum(scatter[0][i,:]) for i in range(len(SIGABS))]
        
        # Calculate the capture XS
        if SIGFISS:
            self.SIGCAP=[a-f for a,f in zip(SIGABS,SIGFISS)]
        else:
            self.SIGCAP=SIGABS

################################################################################
################################################################################
        
def AdjustFISSCHI(FISSCHI, NG):
    """
        Wims produces a fission spectrum as a function of incident neutron
        energy. This routine removes all but the first NG values. It then 
        remove any trailing values which make the normalisation integral 
        greater than one.
            
        Input:
                NG          - number of groups.
                FISSCHI     - Fission spectrum
                
        Return: 
                FISSCHI     - Adjusted fission spectrum
            
    """
    if len(FISSCHI) > NG:
        FISSCHI = [FISSCHI[i] for i in range(NG) ]    
    elif len(FISSCHI) < NG:
        for i in range(NG-len(FISSCHI)):
            FISSCHI.append(0.0) 
    integral = 0.0          
    for i in range(len(FISSCHI)):
        integral+= FISSCHI[i]
        if integral > 1.0:
            FISSCHI[i] = 0.0
    return FISSCHI


################################################################################
################################################################################

def MaterialLists(MATERIALS, matnum):
    """
        Produces a list of nuclides and number densities for given a material 
        number from the list MATERIALS calculated in ReadDataFile.
    
        Input:
                MATERIALS       - list of the specifications produced by WIMS
                matnum          - the material number
                
        Return:
                
    """
    ND=[]; NU=[]
    for row in MATERIALS:
        if row[0]==matnum:
            NU.append(row[1])
            ND.append(row[3])
    return ND, NU

################################################################################
################################################################################

def ReadDataFile(wims):
    """
        Read the wims file and return the micro and macro cross sections.
        
        Input:
                wims      - The open wims file
        Return
                WIMSData  - type Data
    """
    
    
    wims10=True
        
    MATERIALS=[]
    XS_List=[]
    XS_List_Macro=[]
    XS_List_Micro=[]
    NUCLIDES=[]
    
    line = wims.readline()
    while line:

        if(line.split()[0].lower() == 'energy_bounds'):
            SPCTRM=[]
            while True:
                split = wims.readline().split()
                if split[0].lower() in ['nuclide_names']: 
                    wims.seek(wims.tell()-len(line))
                    break
                for e in split:
                    SPCTRM.append(float(e))
            NUM_GRPS=len(SPCTRM)-1
            
        elif(line.split()[0].lower() == 'nuclide_names'):
            line = wims.readline()
            NUCLIDES=line.split()
            
        elif(line.split()[0].lower() == 'material'):
            mat=int(line.split()[1])
            while True:
                line = wims.readline()
                split=line.split()
                if(split[0].lower() == 'micro'):
                    wims.seek(wims.tell()-len(line))
                    break            
                if (split[0].lower() == 'material'): 
                    mat=int(line.split()[1])
                    continue
                
                if wims10:
                    for j in range(0,len(split),4):
                        MATERIALS.append([mat, split[j], int(split[j+1]), float(split[j+2])])
                else:
                    for j in range(0,len(split),3):
                        MATERIALS.append([mat, split[j], int(split[j+1]), float(split[j+2])])

        elif(line.split()[0].lower() == 'micro'):
#            if not MICRO:
#                line = wims.readline()
#                continue
            iso_lbl = int(line.split()[1])
            SIGABS=[]; SIGTRAN=[]; scatter=[]; SIGFISS=[]; FISSNU=[]; FISSCHI=[]
            while True:        
                line = wims.readline()
                if not line: break
                if line.split()[0].lower() in ['micro', 'macro'] :
                    wims.seek(wims.tell()-len(line))
                    break
                elif(line.split()[0].lower() == 'sigabs'):
                    while True:
                        line = wims.readline()
                        for val in line.split():
                            SIGABS.append(float(val))
                        if len(SIGABS) >= NUM_GRPS: break
                elif(line.split()[0].lower() == 'sigtran'):
                    while True:
                        line = wims.readline()
                        for val in line.split():
                            SIGTRAN.append(float(val))
                        if len(SIGTRAN) >= NUM_GRPS: break
                elif(line.split()[0].lower() == 'sigfiss'):
                    while True:
                        line = wims.readline()
                        for val in line.split():
                            SIGFISS.append(float(val))
                        if len(SIGFISS) >= NUM_GRPS: break
                elif(line.split()[0].lower() == 'fissnu'):
                    while True:
                        line = wims.readline()
                        for val in line.split():
                            FISSNU.append(float(val))
                        if len(FISSNU) >= NUM_GRPS: break
                elif(line.split()[0].lower() == 'fisschi'):
                    while True:
                        line = wims.readline()
                        for val in line.split():
                            FISSCHI.append(float(val))
                        if len(FISSCHI) >= NUM_GRPS: break
                elif(line.split()[0].lower() == 'scatter'):
                    scat_mat=zeros([NUM_GRPS,NUM_GRPS],float)
                    for grp in range(NUM_GRPS):
                        line = wims.readline()
                        tmp=line.split()
                        current = int(tmp[1]); start=int(tmp[3]); finish=int(tmp[5])
                        j=start-1
                        while True:
                            line = wims.readline()
                            for val in line.split():
                                scat_mat[grp,j] = float(val)            # grp to current
                                j+=1
                            if j==finish :break                    
                    scatter.append(scat_mat)
#           Make sure we associate the correct data with each isotope according to MATERIALS
            for row in MATERIALS :
                if row[2] == iso_lbl:
                    FISSCHI = AdjustFISSCHI(FISSCHI,NUM_GRPS)
                    XS_List_Micro.append( XS('M'+str(row[0])+row[1], row[1], row[3],'M'+str(row[0])+row[1], SIGABS, SIGTRAN, scatter, SIGFISS, FISSNU, FISSCHI)  )
                    break


        elif(line.split()[0].lower() == 'macro'):
            iso_lbl = int(line.split()[1])
            SIGABS=[]; SIGTRAN=[]; scatter=[]; SIGFISS=[]; FISSNU=[]; FISSCHI=[]        
            while True:  
                line = wims.readline()
                if not line: break

                if line.split()[0].lower() in ['micro', 'macro','fluxes'] :
                    wims.seek(wims.tell()-len(line))
                    break
                if(line.split()[0].lower() == 'sigabs'):
                    while True:
                        line = wims.readline()
                        for val in line.split():
                            SIGABS.append(float(val))
                        if len(SIGABS) >= NUM_GRPS: break
                elif(line.split()[0].lower() == 'sigtran'):
                    while True:
                        line = wims.readline()
                        for val in line.split():
                            SIGTRAN.append(float(val))
                        if len(SIGTRAN) >= NUM_GRPS: break
                elif(line.split()[0].lower() == 'sigfiss'):
                   while True:
                        line = wims.readline()
                        for val in line.split():
                            SIGFISS.append(float(val))
                        if len(SIGFISS) >= NUM_GRPS: break
                elif(line.split()[0].lower() == 'fissnu'):
                    while True:
                        line = wims.readline()
                        for val in line.split():
                            FISSNU.append(float(val))
                        if len(FISSNU) >= NUM_GRPS: break
                elif(line.split()[0].lower() == 'fisschi'):
                    while True:
                        line = wims.readline()
                        if len(FISSCHI) >= NUM_GRPS: break
                        if line.split()[0].lower() in ['micro', 'macro','fluxes'] :
                            wims.seek(wims.tell()-len(line))
                            break
                        for val in line.split():
                            FISSCHI.append(float(val))
                elif(line.split()[0].lower() == 'scatter'):
                    scat_mat=zeros([NUM_GRPS,NUM_GRPS],float)
                    for grp in range(NUM_GRPS):
                        line = wims.readline()
                        tmp=line.split()
                        current = int(tmp[1]); start=int(tmp[3]); finish=int(tmp[5])
                        j=start-1
                        while True:
                            line = wims.readline()
                            for val in line.split():
                                scat_mat[grp,j] = float(val)            # grp to current
                                j+=1
                            if j==finish : break
                    scatter.append(scat_mat)
            name='ma'+str(iso_lbl)
            ND,NU=MaterialLists(MATERIALS, iso_lbl)
#            print ND, NU, iso_lbl
            FISSCHI = AdjustFISSCHI(FISSCHI,NUM_GRPS)
            XS_List_Macro.append(  XS(name, NU, ND, ['M'+str(iso_lbl)+n for n in NU], SIGABS, SIGTRAN, scatter, SIGFISS, FISSNU, FISSCHI)  )

        elif(line.split()[0].lower() == 'fluxes'): break
           

        line = wims.readline()        

#       -------------
#       End of file reading

    Material_Specs=[]
    for i in range(max([row[0] for row in MATERIALS])):
        ND,NU=MaterialLists(MATERIALS, i+1)
        for j in range(len(NU)):
            NU[j] = 'M'+str(i+1)+NU[j]
        Material_Specs.append([ND,NU])

    WIMSData=Data(NUM_GRPS, SPCTRM, XS_List_Micro, XS_List_Macro, NUCLIDES, Material_Specs, True)
    return WIMSData


