"""
    This module reads an NJOY data file and extracts the mean values,
    standard deviations and covarainces for any specified materials and 
    reactions.

    Author:
            Dan Ayres

"""


from itertools import product as tensor_prod
from operator import mul
from math import ceil as ceiling
from numpy import zeros, sqrt
from numpy import power, set_printoptions

from Python_Modules.Plotting import PlotCorrelationMat, minmax

"""
    __combinations defines how to combine all of the sub-reactions into the total reactions such as capture.    
"""

__combinations=[]
__combinations.append([[102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,16,17],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-2]])    # Capture first
__combinations.append([[2,4,16,17],[1,1,2,3]])    # Scatter
__combinations.append([[18],[1]])                 # Fission
__combinations.append([[452],[1]])                # nu-bar



class Data():

    def __init__(self,cov,mean, reactions):
        self.cov=cov
        self.mean=mean
        self.reactions=reactions

################################################################################
################################################################################    

class XS():
    """
        This is a class which conatins the cross section data
        for each MAT number
    """

    def __init__(self,name):
        self.name=name          # Name of the nuclide e.g U235
        self.covariances={}     # Dictionary of covariances indexed by reaction (MT number)
                                # The key is as follows - covariance beteen elastic and (n,g): key=2_102
        reactions=[]            # Contains the reactions present for this nuclide
        self.means={}           # Dictionary of means indexed as above
        self.stds={}            # Dictionary of standard deviations indexed as above
        self.NG=None


################################################################################
################################################################################


def NormWithWimsData(Wims, NJOY, MACRO):

    for nuc in Wims.CrossSections:
        if nuc.SIGFISS==[]:
            Wmeans=zeros([2*Wims.NG],float)
            Wmeans[0:Wims.NG] = nuc.SIGCAP
            Wmeans[Wims.NG:2*Wims.NG] = [nuc.scatter[0][i,i] for i in range(Wims.NG)]
        else:
            Wmeans=zeros([4*Wims.NG],float)
            Wmeans[0:Wims.NG] = nuc.SIGCAP
            Wmeans[Wims.NG:2*Wims.NG] = [nuc.scatter[0][i,i] for i in range(Wims.NG)]
            Wmeans[2*Wims.NG:3*Wims.NG] = nuc.SIGFISS
            Wmeans[3*Wims.NG:4*Wims.NG] = nuc.FISSNU

        if MACRO:
            nuc.CovMat=NormCovMat(len(NJOY.cov[nuc.name]),NJOY.cov[nuc.name],Wmeans)
        else:
            nuc.CovMat=NormCovMat(len(NJOY.cov[nuc.nuclide]),NJOY.cov[nuc.nuclide],Wmeans)


    return Wims      
     


def WIMSMacrosopicData(Wims,data):
    """
        Takes the microscopic data (exception of Zr) produced by NJOY and combine
        them using the number densities provided to produce macroscopic quantities.
        
        The relative uncertainty information (cov,std) are then normalised
        appropriately
        
        Input:
                data        - NJOY miroscopic infomation stored as Data()
                Wims        - wims cross section information
                norm        - if True assume cov and std are relative to mean
        
        Return:
                data        - Wims maroscopic cross sections with covariance data
    """
    
    
    microcovs={};matcov={}
    for nuc in Wims.Micro:
        
        if nuc.nuclide in data.cov:
            WimsMeans=[]
            for r, sig in zip(['102','2','18','452'],[list(nuc.SIGCAP), [nuc.scatter[0][i,i] for i in range(Wims.NG)],list(nuc.SIGFISS), list(nuc.FISSNU)]): 
                if r in data.reactions[nuc.nuclide]: 
                    WimsMeans.extend(sig)
            nuc.CovMat=NormCovMat(len(data.cov[nuc.nuclide]),data.cov[nuc.nuclide],WimsMeans)    
            
    return Wims


################################################################################
################################################################################

def Correlated_Reactions(C, N, NG):
    std=[sqrt(C[i,i]) for i in range(N)]
    ct=zeros([N,N],float)
    for i in range(N):
        for j in range(N):
            ct[i,j] = std[i]*std[j]
    for k in range(4):
        ct[k*NG:(k+1)*NG,k*NG:(k+1)*NG] = C[k*NG:(k+1)*NG,k*NG:(k+1)*NG]
    return ct
    


################################################################################
################################################################################


def NormCovMat(NG,C,mean):
    """
    This is a routine which multiplies a relative covariance matirx C to
    form an absolute covariance matrix.
    """
    CC=zeros([NG,NG],float)
    for i in range(NG):
        for j in range(NG):
            CC[i,j] = mean[i]*mean[j]*C[i,j]
    return CC

################################################################################
################################################################################

def ZeroExpand(Cin, N, Ns):
    """
        Puts zeros into a matrix inplace of fission and nubar which only has 
        capture and scatter.
    """
    C = zeros([N,N],float)
    C[0:Ns,0:Ns] = Cin
    return C

################################################################################
################################################################################

def IsotopeToMAT(wims_names):
    """
        Converts the nuclide name as defined by WIMS into the ENDF format
        MAT number. e.g. U238=9237
        
        Input:
        
        Return:
        
    """

    ENDF_MAT={'H':125, 'O':825,
              'Zr90':4025, 'Zr91':4028, 'Zr92':4031, 'Zr94':4037, 'Zr96':4043,
              'U235':9228, 'U238':9237,
              'Pu238':9434, 'Pu239':9437, 'Pu240':9440, 'Pu241':9443, 'Pu242':9446,
              'Am241':9543}
    Zirc={'Zr90':4025, 'Zr91':4028, 'Zr92':4031, 'Zr94':4037, 'Zr96':4043}

    MAT_NUMS=[]
    NAMES=[]
    for name in wims_names:
        if name in ENDF_MAT:    
            MAT_NUMS.append(ENDF_MAT[name])
            NAMES.append(name)
        if name == 'Zr':
            for key in Zirc:
                MAT_NUMS.append(Zirc[key])
                NAMES.append(key)
      
    return MAT_NUMS, NAMES      


################################################################################
################################################################################

def chunk_float(s, n):
    """
        Split the sting s into n character chunks and 
        return as a list of floats.
    """
    vals=[]
    for start in range(0, len(s)-1, n):
        vals.append(float(s[start:start+n].strip()) )
    return vals

################################################################################
################################################################################

def chunk_int(s, n):
    """
        Split the sting s into n character chunks and 
        return as a list of floats.
    """
    vals=[]
    for start in range(0, len(s)-1, n):
        vals.append(int(s[start:start+n].strip()) )
    return vals    

################################################################################
################################################################################
    
def vecboxer(NG,values,boxes):
    """
        Convert a vector stored in boxer format into
        dense format.
        
        Input:
                NG        - Length of the dense vector
                values    - Value array
                boxes     - Location Information
                
        Return:
                C         - Dense vector  (numpy)
    
    """
    C=zeros([NG],float)
    row=0
    val_pos=0
    for b in boxes:
        if b<0:
            for i in range(abs(b)):
                C[row]=values[val_pos];row+=1
                if row==NG: return C
            val_pos+=1
        else:
            for i in range(b):
                C[row]=0.0;row+=1
                if row==NG: return C
    return C    
    
################################################################################
################################################################################    
    
def boxer(NG,values,boxes,sym_flag):
    """
        Convert a matrix stored in boxer format into
        dense format.
        
        Input:
                NG        - Length of the dense vector
                values    - Value array
                boxes     - Location Information
                
        Return:
                C         - Dense matrix  (numpy)
    
    """    
    C=zeros([NG,NG],float) #C[row,col]
    if sym_flag==0:
        col=0
        row=0
        row_above=C[:,row]
        val_pos=0; box_pos=0
        
        for b in boxes:
            if b<0:
                for i in range(abs(b)):
                    C[row,col] = values[val_pos]
                    C[col,row] = C[row,col]
                    if col==NG-1:
                        if row==NG-1: return C
                        row+=1
                        row_above=C[:,row-1] #
                        col=row
                    else:
                        col+=1
                val_pos+=1
            else:
                if row>0: row_above=C[:,row-1]
                for i in range(b):
                    C[row,col] = row_above[col]
                    C[col,row] = C[row,col]
                    if col==NG-1:
                        if row==NG-1: return C
                        row+=1
                        row_above=C[:,row-1] #
                        col=row
                    else:
                        col+=1                    
    else:
        col=0
        row=0
        row_above=C[:,row]
        val_pos=0; box_pos=0
        
        for b in boxes:
            if b<0:
                for i in range(abs(b)):
                    C[row,col] = values[val_pos]
                    if col==NG-1:
                        if row==NG-1: return C
                        row+=1
                        row_above=C[:,row-1] #
                        col=0
                    else:
                        col+=1
                val_pos+=1
            else:
                if row>0: row_above=C[:,row-1]
                for i in range(b):
                    C[row,col] = row_above[col]
                    if col==NG-1:
                        if row==NG-1: return C
                        row+=1
                        row_above=C[:,row-1] #
                        col=0
                    else:
                        col+=1  
    
################################################################################
################################################################################
    
    
def CalcVectorValues(values, NG):
    """
        Combine component reactions to form the following
            - total capture
            - scatter
            - fission
            - nu-bar
               
        Input:
                values      - dictionary of reactions (mean or standard deviations)
                NG          - number of values in each reactions (energy groups)
                  
        Return:
                vals        - vector of values in the following order
                                  [capture,scatter,fission,nubar]
        
    """
#    if str(18) not in values:
#        maxi = 2
#        vals=zeros([maxi*NG],float)
#    else:
#        maxi=4
#        vals=zeros([maxi*NG],float)

    reacs=['102','2','18','452']

    maxi=0
    comb=[]
    present=[]
    for r,c in zip(reacs,__combinations):
        if r in values: 
            maxi+=1
            comb.append(c)
            present.append(r)

    vals=zeros(maxi*NG,float)
    if maxi==0:
        print 'No reactions'
        
    for i in range(maxi):
        ct=zeros([NG],float)
            
        for r,w in zip(comb[i][0],comb[i][1]):
            if str(r) not in values: continue
            ct+=w*values[str(r)]
                
        vals[ i*NG:(i+1)*NG]=ct    
       
    return vals, present
    

    
################################################################################
################################################################################


def weights_and_keys(wts,kys,wts1,kys1):
    """
    return a list of the square of wts and kys
    """
    weight=[]
    for i in list(tensor_prod(wts,wts1)) :
        weight.append(reduce(mul,i))  
    
    keys=[]
    flags=[]
    for i in list(tensor_prod(kys,kys1)) :
        if i[0]>i[1]:
            flags.append(True)
        else:
            flags.append(False)
        i=sorted(i)
        keys.append(str(i[0])+'_'+str(i[1]))
        
    return weight,keys,flags


################################################################################
################################################################################


def ReverseGroupOrder(NG, C):
    """Change the group ordering of a covariance matrix"""
    
    temp=zeros([NG,NG],float)
    for i in range(NG):
        for j in range(NG):
            temp[i,j] = C[NG-i-1,NG-j-1]
    return temp
 

################################################################################
################################################################################


def ReverseGroupOrderVec(NG, C):
    """Change the group ordering of a covariance matrix"""
    
    temp=zeros([NG],float)
    for i in range(NG):
        temp[i] = C[NG-i-1]
    return temp 
    
################################################################################
################################################################################      



def RelativeCov(C,Mi,Mj,stdi,stdj):
    """
        Create an absolute covariance matrix from a relative one with
        the mean values Mi and Mj
    """  
        
    bodge=True
        
    N=len(C)
    CC=zeros([N,N],float)
    if N!=len(Mi) or N!=len(Mj):
        print 'Mean vector and Cov Mat do not match'
        raise SystemExit(0)

    if bodge:
        for i in range(N):
            for j in range(N):
                if stdj[j]==0.0 or stdi[i]==0.0 : continue
                else:
                    if C[i,j]/(stdi[i]*stdj[j]) > 1.1:
                        C[i,j] = 0.9*stdi[i]*stdj[j]
                    if C[i,j]/(stdi[i]*stdj[j]) < -1.1:
                        C[i,j] = -0.9*stdi[i]*stdj[j]

            
    for i in range(N):
        for j in range(N):
            CC[i,j] = Mi[i]*Mj[j]*C[i,j]
            
    return CC

            
################################################################################
################################################################################
    
def ReadDataFile(njoy,wims_names, reactions):
    """
        Read in all of the specified data from the njoy file.
        
        Input:
                njoy         - open NJOY data file
                wims_names   - names of the nuclides to extract
    """    

    mean=True
    std=True
    cov=True

    materials,names=IsotopeToMAT(wims_names)

#    reactions=[2,18,102,452]
#    reactions=[2,4,16,17,18,102,452]
#    reactions=[2,4,16,18,102,452]
#    reactions=[2,4,16,17,18,102,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,452]


    nuclides={} #Dictionary of all the isotopes
    for mat,name in zip(materials,names):
        nuclides[str(mat)] = XS(name)

    line=njoy.readline()
    while line:

        split=line.split()
        if int(split[0])==1: #Read the mean values
                  
            NG=int(split[2]);MAT=int(split[6]);MT=int(split[7]);MAT1=int(split[8]);MT1=int(split[9])
            val_lines=int(ceiling(float(split[10])*10/80.0))
            box_lines=int(ceiling(float(split[12])*float(split[13])/80.0))
            sym_flag=int(split[16])            
            
            
            
            if not mean:
                for i in range(val_lines): njoy.readline()
                for i in range(box_lines): njoy.readline()
            
            elif MAT in materials and MT in reactions and MT1 in reactions:
                
                values=[];boxes=[]
                
                for i in range(val_lines):        
                    line=njoy.readline()
                    for val in chunk_float(line,10): values.append(val)
                for i in range(box_lines):
                    line=njoy.readline()
                    for box in chunk_int(line,3): boxes.append(box)
                    
                nuclides[str(MAT)].means[str(MT)]=vecboxer(NG,values,boxes)
                if nuclides[str(MAT)].NG==None: nuclides[str(MAT)].NG=NG

            else: 
                for i in range(val_lines): njoy.readline()
                for i in range(box_lines): njoy.readline()         


        elif int(split[0])==2: #Read the standard deviations (may be relative)
            NG=int(split[2]);MAT=int(split[6]);MT=int(split[7]);MAT1=int(split[8]);MT1=int(split[9])
            val_lines=int(ceiling(float(split[10])*10/80.0))
            box_lines=int(ceiling(float(split[12])*float(split[13])/80.0))
            sym_flag=int(split[16])     
            
            if not std:
                for i in range(val_lines): njoy.readline()
                for i in range(box_lines): njoy.readline()            
                   
            elif MAT in materials and MT in reactions and MT1 in reactions:
                
                values=[];boxes=[]
                
                for i in range(val_lines):        
                    line=njoy.readline()
                    for val in chunk_float(line,10): values.append(val)
                for i in range(box_lines):
                    line=njoy.readline()
                    for box in chunk_int(line,3): boxes.append(box)
                    
                nuclides[str(MAT)].stds[str(MT)]=vecboxer(NG,values,boxes)
                if nuclides[str(MAT)].NG==None: nuclides[str(MAT)].NG=NG

            else: 
                for i in range(val_lines): njoy.readline()
                for i in range(box_lines): njoy.readline()      

        elif int(split[0])==3: #Read the covariances
            val_lines=int(ceiling(float(split[10])*10/80.0))
            box_lines=int(ceiling(float(split[12])*float(split[13])/80.0))
            sym_flag=int(split[16])
            NG=int(split[2]);MAT=int(split[6]);MT=int(split[7]);MAT1=int(split[8]);MT1=int(split[9])
    

            if not cov:
                for i in range(val_lines): njoy.readline()
                for i in range(box_lines): njoy.readline()
            
            if MAT in materials and MT in reactions and MT1 in reactions:
                if MAT != MAT1:
                    print 'correlated materials!!!'
                    continue
                values=[];boxes=[]
                
                for i in range(val_lines):        
                    line=njoy.readline()
                    for val in chunk_float(line,10): values.append(val)
                for i in range(box_lines):
                    line=njoy.readline()
                    for box in chunk_int(line,4): boxes.append(box)
                    
                nuclides[str(MAT)].covariances[str(MT)+'_'+str(MT1)]=boxer(NG,values,boxes,sym_flag)
#                if MT not in nuclides[str(MAT)].reactions:
#                    nuclides[str(MAT)].reactions.append(MT)
                if nuclides[str(MAT)].NG==None: nuclides[str(MAT)].NG=NG
                
            else: 
                for i in range(val_lines): njoy.readline()
                for i in range(box_lines): njoy.readline()         
        else: 
            MAT=int(split[6])
            val_lines=int(ceiling(float(split[10])*10/80.0))
            box_lines=int(ceiling(float(split[12])*float(split[13])/80.0))
            for i in range(val_lines): njoy.readline()
            for i in range(box_lines): njoy.readline()
            
        line=njoy.readline()

#       -------------------------
#       End of file reading

    return nuclides
    
    
    
def FormatNJOYData(nuclides, wims_names):

    """
        Zirconium is used in wims in natural form i.e a composition of
        all of the naturally occuring isotopes. The data processed by NJOY
        is isotopic so we must combine the isotopes according to their
        relative abundance.
    """

    materials,names=IsotopeToMAT(wims_names)

    Zirc={4025:0.5144999435258152, \
          4028:0.112200005753799 , \
          4031:0.17150002742767534,\
          4037:0.17380002070257586,\
          4043:0.0280000025901346}


    cross_reactions=True#False
    meanxs={};stdxs={};covxs={}; reacs={}


    missed_mats=[]
    reactions=['2_2','18_18','102_102','452_452']
    missed_reactions={}
    ignore_list=[key for key in Zirc]
    ignore_list.append(125) #H
    ignore_list.append(825) #O
    

    for mat in materials:
    
        if nuclides[str(mat)].means=={}: 
            missed_reactions[nuclides[str(mat)].name]=reactions
            continue
    
    
        # First, reverse the order of the group stucture to match WIMS
        for key in nuclides[str(mat)].means:
            nuclides[str(mat)].means[key]=ReverseGroupOrderVec(nuclides[str(mat)].NG, nuclides[str(mat)].means[key])
        for key in nuclides[str(mat)].stds:
            nuclides[str(mat)].stds[key]=ReverseGroupOrderVec(nuclides[str(mat)].NG, nuclides[str(mat)].stds[key])
        for key in nuclides[str(mat)].covariances:
            nuclides[str(mat)].covariances[key]=ReverseGroupOrder(nuclides[str(mat)].NG, nuclides[str(mat)].covariances[key])

        temp=[key for key in nuclides[str(mat)].covariances]
        temp2=[r for r in reactions if r not in temp]
        if temp2 != [] and mat not in ignore_list:
            missed_reactions[nuclides[str(mat)].name]=temp2


        # Calculate the absolute covariances
        for key in nuclides[str(mat)].covariances:
            ki=key.split('_')[0]; kj=key.split('_')[1]
            Mi=nuclides[str(mat)].means[ki]; Mj=nuclides[str(mat)].means[kj]
            stdi=nuclides[str(mat)].stds[ki]; stdj=nuclides[str(mat)].stds[kj]
            nuclides[str(mat)].covariances[key]=RelativeCov(nuclides[str(mat)].covariances[key],Mi,Mj,stdi,stdj)


        # Calculate the 'Big' covariance matrix and mean vector
        meanvec, present_r = CalcVectorValues(nuclides[str(mat)].means, nuclides[str(mat)].NG)
        covmat=CalcCovMat(nuclides[str(mat)].covariances, nuclides[str(mat)].means, nuclides[str(mat)].NG, cross_reactions)
        reacs[nuclides[str(mat)].name]=present_r

        # Combine all of the zirconium isotopes
        if mat in Zirc:      
            if 'Zr' in meanxs:
                meanxs['Zr'] += Zirc[mat]*meanvec
            else:
                meanxs['Zr'] = Zirc[mat]*meanvec
                
            if 'Zr' in covxs:
                covxs['Zr'] += power(Zirc[mat],2)*covmat
            else:
                covxs['Zr'] = power(Zirc[mat],2)*covmat
                
        else:
            meanxs[nuclides[str(mat)].name] = meanvec
            covxs[nuclides[str(mat)].name] = covmat
        
    # Create relative covariances using NJOY means
    for mat in covxs:
        covxs[mat]=RelCov(covxs[mat],meanxs[mat])
          
        
    NJOYData = Data(covxs,meanxs, reacs)         

    return NJOYData


################################################################################
################################################################################      


def TestCorrelation(C,Mi,Mj):
    from numpy import ravel, amin,amax, isnan
    """
        Form the correlation matrix from the covariance matrix and check that
        it is bounded.
    """  
        
    N=len(C)
    CC=zeros([N,N],float)
    rsd=[sqrt(abs(C[i,i])) if abs(C[i,i])>0.0 else 0.0 for i in range(N)]
    diag = [C[i,i] for i in range(N)]




    for i in range(N):
        for j in range(N):
            if Mj[j]==0.0 or Mi[i]==0.0 :
                CC[i,j]==0.0
            else:
                CC[i,j] = C[i,j]/(Mi[i]*Mj[j])

#            if rsd[j]==0.0 or rsd[i]==0.0 :
#                CC[i,j]==0.0
#            else:
#                CC[i,j] = C[i,j]/(rsd[i]*rsd[j]M


            if CC[i,j] > 1.01: print CC[i,j], C[i,j], Mi[i],Mj[j], i+1,j+1 
            if CC[i,j] < -1.01: print CC[i,j], C[i,j], Mi[i],Mj[j], i+1,j+1 
                            
    print Mi
    print
    print Mj
    print

#    for i,j,d,r in zip(Mi,Mj,diag,rsd):
#        print "%15.8E%15.8E%15.8E"%(i*j, d,d/(i*j))

#    print '\n\n'
#    print C

                
    dd=ravel(CC)
#    print CC
    return amin(dd), amax(dd)

################################################################################
################################################################################

def FindMissing(Wims_nuclides, NJOY_Data):
    """
        Find the nuclides and reactions that are missing from the NJOY
        covariance data.
    """
    
    
    NJOY_nuclides=[cov for cov in NJOY_Data.cov]
    for cov in Wims_nuclides:
        print cov, cov in NJOY_nuclides
        
    print NJOY_Data.cov['U235']

    
################################################################################
################################################################################

def CalcCovMat(covariances, means, NG, cross):
    """
        Combine component reactions to form the following covarinces and 
        the cross covariances between if specified
            - total capture
            - scatter
            - fission
            - nu-bar
                
        Input:
                values      - dictionary of reactions (mean or standard deviations)
                NG          - number of values in each reactions (energy groups)
                cross       - if True evaluate the cross covariances between reactions
                    
        Return:
                C           - matrix of values in the following order
                              [ [c-c,c-s,c-f,c-n],
                                [s-c,s-s,s-f,s-n],
                                [f-c,f-s,f-f,f-n],
                                [n-c,n-s,n-f,n-n] ]
       
    """    

    set_printoptions(edgeitems=5,infstr='inf', linewidth=200, nanstr='nan', precision=3, suppress=False, threshold=1000) 


    reacs=['102_102','2_2','18_18','452_452']

    maxi=0
    comb=[]
    for r,c in zip(reacs,__combinations):
        if r in covariances: 
            maxi+=1
            comb.append(c)
            
    C=zeros([maxi*NG,maxi*NG],float)
    if maxi==0:
        print 'No reactions'
        
    for i in range(maxi):
        for j in range(i+1):

            """ Ignore cross correlations """        
            if not cross: 
                if i!=j: continue
                
            ct=zeros([NG,NG],float)
            weights,keys,flags=weights_and_keys(comb[i][1],comb[i][0],comb[j][1],comb[j][0])
            
            for ctr in range(len(keys)):
                
                k=keys[ctr]         # Key for the covariances e.g. 2_18
                w=weights[ctr]      # Weight for the combination
                flag=flags[ctr]     # True: take tranpose

                if k in covariances: 
#                    print k, i,j,flag
                    ct=zeros([NG,NG],float)
#                    ki=k.split('_')[0]; kj=k.split('_')[1]
                    if flag: ct=w*covariances[k].T
                    else: ct=w*covariances[k]
                    
                    C[i*NG:(i+1)*NG, j*NG:(j+1)*NG]+=ct
                    if i!=j:C[j*NG:(j+1)*NG, i*NG:(i+1)*NG]+=ct.T
                    
    return C


def ZeroCovs(C,NG):
    CC=zeros([NG,NG],float)
    for i in range(NG):
        CC[i,i]=C[i,i]
    return CC

    
def RelCov(C,mean):
    N=len(C)
    CC=zeros([N,N],float)
    for i in range(N):
        for j in range(N):
            if mean[i]==0.0 or mean[j]==0.0:
                CC[i,j] == 0.0
            else:
                CC[i,j] = C[i,j]/(mean[i]*mean[j])

    return CC
    
