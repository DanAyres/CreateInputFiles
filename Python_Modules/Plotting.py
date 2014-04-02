from matplotlib.pyplot import figure, show, axes, sci, subplot
from matplotlib import cm, colors
from matplotlib.font_manager import FontProperties
from numpy import amin, amax, ravel, sqrt, zeros, arange, set_printoptions
from numpy import linalg
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker

class ImageFollower:
    'update image in response to changes in clim or cmap on another image'
    def __init__(self, follower):
        self.follower = follower
    def __call__(self, leader):
        self.follower.set_cmap(leader.get_cmap())
        self.follower.set_clim(leader.get_clim())

################################################################################
################################################################################

def minmax(C):
    N=len(C)
    
    Cor=zeros([N,N],float)
    rsd=[]
    
    for i in range(N):
        if C[i,i] >0.0:
            rsd.append(sqrt(C[i,i]))
        else:
            rsd.append(None)
    for i in range(N):
        if rsd[i]==None:
            Cor[i,:] = [0.0 for j in range(N)]   
        else:
            for j in range(N):
                if rsd[j] ==None:
                    Cor[i,j] = 0.0 
                else:
                    Cor[i,j] = C[i,j]/(rsd[i]*rsd[j])
    
    dd=ravel(Cor)
    return amin(dd),amax(dd)

################################################################################
################################################################################

def WriteEigenSpectra(WimsData, nuclides, folder):
    """
        Write the Cumulative eigenvalues of each covariance matrix to file
    """
    
    filename = folder+'/'+'EigenSpectra.eigs'
    
    e=[]
    tag=[]
    for iso in WimsData.CrossSections:
#        if iso.nuclide=='Am241': continue
        
        
        if iso.CovMat!=None:
            U, eigs, V = linalg.svd(iso.CovMat)
            eigs=eigs/sum(eigs)
            new=[]
            for i in range(len(eigs)):
                new.append(sum(eigs[:i+1]))
            e.append(new)
            tag.append(iso.name)

    maxnum=max([len(row) for row in e])

    for i in range(len(e)):
        if len(e[i]) < maxnum: 
            for j in range( maxnum-len(e[i])):
                e[i].append(1.0)
   
    f=open(filename,'w')
    f.write('#   ')
    for name in tag: f.write('%15s'%name)
    f.write('\n')
    for j in range(maxnum):
        f.write('%4i'%(j+1))
        for i in range(len(e)):
            f.write("%15.4E" %e[i][j])
        f.write('\n')
    
    f.close()

################################################################################
##################1##############################################################

def PlotCorrelationData(Wims,nuclides,folder):

    for nuc in Wims.CrossSections:   
        if nuc.CovMat != None:
            name=folder+'/'+nuc.name+'.pdf'
            PlotCorrelationMat(name, Wims.NG,  nuc.CovMat)

################################################################################
################################################################################
    
def PlotCorrelationMat(name, NG,  C):
    from math import isnan

#    set_printoptions(edgeitems=3,infstr='inf', linewidth=200, nanstr='nan', precision=3, suppress=False, threshold=1000) 
    
    N=len(C)
    
    Cor=zeros([N,N],float)
    
    rsd=[sqrt(C[i,i]) for i in range(N)]
    
    for i in range(N):
        for j in range(N):
            if rsd[i]==0.0 or rsd[j]==0.0:
                Cor[i,j] ==0.0
            elif isnan(rsd[j]) or isnan(rsd[i]):
                Cor[i,j] ==0.0
            else:
                Cor[i,j] = C[i,j]/(rsd[i]*rsd[j])
                
#                if Cor[i,j] > 1.0001: print Cor[i,j],i,j
    dd=ravel(Cor)
    out=name.split('/')[1].split('.')[0]
#    if out=='M1Pu242': 
#        print Cor[0:26,0:26]
#        print 
#        print [Cor[i,i] for i in range(N)]
    print out.ljust(10),'Min Correlation:', '%10.5E'%amin(dd), '     Max Correlation: ','%10.5E'%amax(dd)
    
    fig = figure()
    cmap = cm.seismic    
    
    ax=fig.add_subplot(111)
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    ax.set_xticklabels([])
    ax.set_yticklabels([])    
    cax=ax.imshow(Cor,  interpolation='nearest', cmap=cmap, vmin=-1.0,vmax=1.0)
    bar=fig.colorbar(cax)
    bar.set_ticks([-1.0,0.0,1.0])#,2.0,3.0])
    

    pp = PdfPages(name+'_Cov.pdf')
    pp.savefig()
    pp.close()

