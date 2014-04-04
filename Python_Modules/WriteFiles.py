
import numpy as np

#from os import chmod
#from stat import S_IXUSR

def GemEVENT_scripts(name):
    run=open(name+'/'+name+'.sh','w')
    run.write("#!/bin/bash \n")
    run.write("\n")
    run.write("event %s.event" %name)
    run.close()
    
    setup=open(name+'/'+name+'_Setup.sh','w')
    setup.write("#!/bin/bash \n")
    setup.write("\n")
    setup.write("gem %s.gem" %name)
    setup.close()
    

def Sampy_control(name):

    filename=name+'/'+name+'.control'
    outfile=open(filename,'w')

    outfile.write("# \n")
    outfile.write("# General parameters \n")
    outfile.write("# \n")
    outfile.write("Run Ref:            %s \n" %name)
    outfile.write("Param Infile:       %s.paramfile \n" %name)
    outfile.write("Template file:      %s.template \n" %name)
    outfile.write("input file:         %s.event/data/matxs/%s \n" %(name,name))
    outfile.write("output file:        %s.event/out/totals/%s \n" %(name,name))
    outfile.write("Run script:         %s.sh \n" %name)
    outfile.write("setup script:       %s_Setup.sh \n" %name)
    outfile.write("setup files:        %s.gem \n" %name)
    outfile.write("Num procs:          2 \n")
    outfile.write("# \n")
    outfile.write("# Monte Carlo parameters \n")
    outfile.write("# \n")
    outfile.write("Use MC:         yes \n")
    outfile.write("MC Method:      latinhyper \n")
    outfile.write("MC Runs:        10 \n")
    outfile.write("MC dump:        10         # ouput results every <MC Dump> realisations \n")
    outfile.write("MC PDF:         yes \n")
    outfile.write("MC PDF bins:    100 \n")
    outfile.write("# \n")
    outfile.write("#   SCM Parameters \n")
    outfile.write("# \n")
    outfile.write("Use SCM:            no \n")
    outfile.write("SCM SG level:       4 \n")
    outfile.write("SCM Quad Adapt:     yes  \n")
    outfile.write("SCM Quad Tol:       1.0E-5  \n")    
    outfile.write("SCM PCE:            no \n")
    outfile.write("SCM Poly Order:     2 \n")
    outfile.write("SCM PDF:            no \n")
    outfile.write("SCM PDF samp:       100000 \n")
    outfile.write("SCM PDF bins:       100 \n")    
    outfile.write("# \n")
    outfile.write("#   HDMR Parameters \n")
    outfile.write("# \n")
    outfile.write("Use HDMR:               no \n")
    outfile.write("HDMR SG level:          4,4 \n")
    outfile.write("HDMR Quad Adapt:        yes  \n")
    outfile.write("HDMR Quad Tol:          1.0E-4  \n")
    outfile.write("HDMR PCE:               yes \n")
    outfile.write("HDMR Poly Order:        3       # Maximum polynomial order \n")
    outfile.write("HDMR PDF:               no \n")
    outfile.write("HDMR PDF samp:          100000 \n")
    outfile.write("HDMR PDF bins:          100 \n")
    outfile.write("HDMR Adapt:             yes \n")
    outfile.write("HDMR metric:            1E-3 \n")
    
    outfile.close()




####################################################################################################
####################################################################################################

def WriteConstant(openfile, const, NG, name):

    openfile.write("CONSTANT   %i   %s" %(NG, name) )
    openfile.write("   [")
    for i in range(len(const)-1):
        openfile.write("%11.5E," %const[i])
    openfile.write("%11.5E" %const[-1])
    openfile.write("] ")
    openfile.write("\n")

####################################################################################################
####################################################################################################

def WriteParameterFiles(outfile, nuclides, WimsData, MACRO):
    """
        Creates a parameter file for each of the nuclides in the problem.
        
        Input:
                 outfile     -
                 nuclides    -
        
        Return:
    
    """

    writefile=outfile+'/'+outfile+'.paramfile'
    f=open(writefile,'w')
    f.write("# \n")
    f.write("#   Name:       A tag describing the group. (e.g. Sig_f_Fuel) \n")
    f.write("#   Groups:     The number of parameters associated with that name \n")
    f.write("#   Means:      The nominal values of each parameter. \n")
    f.write("#   sigma:      The standard deviation of each parameter or a single value \n")
    f.write("#               followed by '%' describing R=sigma/mean. \n")
    f.write("#   depend:     If 'yes' then all parameters in a group will have a coefficient  \n")
    f.write("#               of correlation, rho = 1. If 'no' then rho = 0. If 'corr' then the  \n")
    f.write("#               correlation and the uncertainties are described by a covariance matrix. \n")
    f.write("#   dist:       The type of probability distibution. Normal and LogNormal are currently  \n")
    f.write("#               supported. \n")
    f.write("#   Cov_File:   Name of the covarince file. If NULL then the covariance matrix is \n")
    f.write("#               generated using the exponential kernel and the values of sigma. \n")
    f.write("#   Cov_len:    Correlation length parameter for the exponential kernel. \n")
    f.write("#   decomp:     chol = Cholesky fatorisation. KL = Karhunen Loeve. \n")
    f.write("# \n")
    f.write("#Name   Groups  Means             sigma     depend  dist        Cov_File    Cov_len     decomp \n")    
    
    
    for n in nuclides:
        f.write('FILE   Data.%s\n' %n)
    
    
    f.write("# \n")
    f.write("#   Constraints. \n")
    f.write("#   Must relate to those parameters defined in the template file \n")
    f.write("#   Must be in the form X = something \n")
    f.write("#   Must define constants before constraints (if any needed) \n")
    f.write("# \n")    
    
    flag=False
    
    
    for nuc in WimsData.Micro:
        if nuc.CovMat==None: continue
        const = [sum(nuc.scatter[0][i,:]) - nuc.scatter[0][i,i]  for i in range(WimsData.NG)]
        
        if len(nuc.reactions)==3:
        
            if 's' not in nuc.reactions:
                for i in range(WimsData.NG):
                    const[i] += nuc.scatter[0][i,i]
            elif 'f' not in nuc.reactions:
                for i in range(WimsData.NG):
                    const[i] += nuc.SIGFISS[i]
                constname='C'+nuc.name+'_f'
                WriteConstant(f, nuc.SIGFISS, WimsData.NG, constname)
            elif 'c' not in nuc.reactions:
                for i in range(WimsData.NG):
                    const[i] += nuc.SIGCAP[i]
                constname='C'+nuc.name+'_c'
                WriteConstant(f, nuc.SIGCAP, WimsData.NG, constname)
            elif 'n' not in nuc.reactions:
                constname='C'+nuc.name+'_n'
                WriteConstant(f, nuc.FISSNU, WimsData.NG, constname)    
        
        elif len(nuc.reactions)==2:
        
            if ('s' not in nuc.reactions) and ('f' not in nuc.reactions):
                for i in range(WimsData.NG):
                    const[i] += (nuc.scatter[0][i,i] + nuc.SIGFISS[i])
                constname='C'+nuc.name+'_f'
                WriteConstant(f, nuc.SIGFISS, WimsData.NG, constname) 
            elif ('s' not in nuc.reactions) and ('n' not in nuc.reactions):
                for i in range(WimsData.NG):
                    const[i] += nuc.scatter[0][i,i]
                constname='C'+nuc.name+'_n'
                WriteConstant(f, nuc.FISSNU, WimsData.NG, constname)     
            elif ('s' not in nuc.reactions) and ('c' not in nuc.reactions):   
                for i in range(WimsData.NG):
                    const[i] += (nuc.scatter[0][i,i] + nuc.SIGCAP[i])
                constname='C'+nuc.name+'_c'
                WriteConstant(f, nuc.SIGCAP, WimsData.NG, constname)
            elif ('f' not in nuc.reactions) and ('n' not in nuc.reactions):
                for i in range(WimsData.NG):
                    const[i] += nuc.SIGFISS[i]
                constname='C'+nuc.name+'_f'
                WriteConstant(f, nuc.SIGFISS, WimsData.NG, constname) 
            elif ('f' not in nuc.reactions) and ('c' not in nuc.reactions):
                for i in range(WimsData.NG):
                    const[i] += (nuc.SIGCAP[i] + nuc.SIGFISS[i])
                constname='C'+nuc.name+'_f'
                WriteConstant(f, nuc.SIGFISS, WimsData.NG, constname)
            elif ('n' not in nuc.reactions) and ('c' not in nuc.reactions):
                for i in range(WimsData.NG):
                    const[i] += nuc.SIGCAP[i]
                constname='C'+nuc.name+'_c'
                WriteConstant(f, nuc.SIGCAP, WimsData.NG, constname)
                constname='C'+nuc.name+'_n'
                WriteConstant(f, nuc.FISSNU, WimsData.NG, constname)
                
        elif len(nuc.reactions)==1:
            if 's' in nuc.reactions:
                for i in range(WimsData.NG):
                    const[i] += (nuc.SIGCAP[i] + nuc.SIGFISS[i])
            elif 'c' in nuc.reactions:
                for i in range(WimsData.NG):
                    const[i] += (nuc.scatter[0][i,i] + nuc.SIGFISS[i])
                constname='C'+nuc.name+'_f'
                WriteConstant(f, nuc.SIGFISS, WimsData.NG, constname)
            elif 'f' in nuc.reactions:
                for i in range(WimsData.NG):
                    const[i] += (nuc.scatter[0][i,i] + nuc.SIGCAP[i])
                constname='C'+nuc.name+'_n'
                WriteConstant(f, nuc.FISSNU, WimsData.NG, constname)
                constname='C'+nuc.name+'_c'
                WriteConstant(f, nuc.SIGCAP, WimsData.NG, constname)
            elif 'n' in nuc.reactions:
                constname='C'+nuc.name+'_f'
                WriteConstant(f, nuc.SIGFISS, WimsData.NG, constname)
                
        if len(nuc.reactions)==1 and 'n' in nuc.reactions:
            continue
        else:
            constname='C'+nuc.name
            WriteConstant(f, const, WimsData.NG, constname)
        
        if not nuc.SIGFISS:
            flag=True

    if(flag):
        f.write("CONSTANT   %i   zero  [" %WimsData.NG)
        for i in range(WimsData.NG-1):
            f.write("0.0,")
        f.write("0.0]\n")
    
    #
    # Write constraints
    #
    for nuc in WimsData.Micro:
        if nuc.CovMat==None: continue
        if len(nuc.reactions)==4:    
            f.write("CONSTRAINT  %i   %s_s = %s_s \n" %(WimsData.NG, nuc.name, nuc.name) )
            f.write("CONSTRAINT  %i   %s_t = %s_c + %s_s + %s_f + C%s \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name, nuc.name, nuc.name) )
            f.write("CONSTRAINT  %i   %s_a = %s_c + %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
            f.write("CONSTRAINT  %i   %s_nf = %s_n * %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )

        elif len(nuc.reactions)==3:
            if 's' not in nuc.reactions:
                f.write("CONSTRAINT  %i   %s_t = %s_c + %s_f + C%s \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_a = %s_c + %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_nf = %s_n * %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
            elif 'f' not in nuc.reactions:
                f.write("CONSTRAINT  %i   %s_t = %s_c + %s_s + C%s \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_a = %s_c + C%s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_nf = %s_n * C%s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )            
            elif 'c' not in nuc.reactions:
                f.write("CONSTRAINT  %i   %s_t = %s_s + %s_f + C%s \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_a = C%s_c + %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_nf = %s_n * %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
            elif 'n' not in nuc.reactions:
                f.write("CONSTRAINT  %i   %s_t = %s_c + %s_s + %s_f + C%s \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_a = %s_c + %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_nf = C%s_n * %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                
        elif len(nuc.reactions)==2:
            if ('s' not in nuc.reactions) and ('f' not in nuc.reactions):
                f.write("CONSTRAINT  %i   %s_t = %s_c + C%s \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_a = %s_c + C%s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_nf = %s_n * C%s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )            
            elif ('s' not in nuc.reactions) and ('n' not in nuc.reactions):
                f.write("CONSTRAINT  %i   %s_t = %s_c + %s_f + C%s \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_a = %s_c + %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_nf = C%s_n * %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
            elif ('s' not in nuc.reactions) and ('c' not in nuc.reactions):
                f.write("CONSTRAINT  %i   %s_t = %s_f + C%s \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_a = C%s_c + %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_nf = %s_n * %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )      
            elif ('f' not in nuc.reactions) and ('n' not in nuc.reactions):
                f.write("CONSTRAINT  %i   %s_t = %s_c + %s_s + C%s \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_a = %s_c + C%s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
            elif ('f' not in nuc.reactions) and ('c' not in nuc.reactions):
                f.write("CONSTRAINT  %i   %s_t = %s_s + C%s \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_nf = %s_n * C%s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
            elif ('n' not in nuc.reactions) and ('c' not in nuc.reactions):
                f.write("CONSTRAINT  %i   %s_t = %s_s + %s_f + C%s \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_a = C%s_c + %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_nf = C%s_n * %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )

        elif len(nuc.reactions)==1:
            if 's' in nuc.reactions:
                f.write("CONSTRAINT  %i   %s_t = %s_s + C%s \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
            elif 'c' in nuc.reactions:
                f.write("CONSTRAINT  %i   %s_t = %s_c + C%s \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_a = %s_c + C%s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
            elif 'f' in nuc.reactions:
                f.write("CONSTRAINT  %i   %s_t = %s_f + C%s \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_a = C%s_c + %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
                f.write("CONSTRAINT  %i   %s_nf = C%s_n * %s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
            elif 'n' in nuc.reactions:
                f.write("CONSTRAINT  %i   %s_nf = %s_n * C%s_f \n" %(WimsData.NG, nuc.name, nuc.name, nuc.name) )
    f.close()    
    
    exts=['_c','_s','_f','_n']

    for n in nuclides:
        numsets=0
        itemsperset=0
        for iso in WimsData.CrossSections:
            if iso.nuclide == n:
                numsets+=1
                itemsperset=len(iso.CovMat)
        
        writefile=outfile+'/Data.'+n
        f=open(writefile,'w') 
        f.write('Number of Groups: %i \n' %WimsData.NG)       
        f.write('Distribution: normal \n')
        f.write('Dependency: corr \n')
        f.write('Decomposition: KL \n')
        f.write('KL theta: 1.0 \n')
        f.write('Cov Mat: yes \n')
        f.write('Cov Len: 0.1 \n')
        f.write('Number of sets: %i \n' %numsets)
        f.write('Items per set: %i \n' %itemsperset)

        counter=1
            
        for iso in WimsData.CrossSections:
            if iso.nuclide == n:
                    
                f.write('Set %5i   %s \n' %(counter,iso.name))
                counter+=1
                    
                for val, lab in  zip(iso.MP,iso.labels):
                    f.write("%s"%lab)
                    f.write("%20.7E"%val)
                    f.write("    1%")
                    f.write("\n")                    

                length=len(iso.CovMat)

                f.write('CovMat  %i  %s\n' %(length,iso.name))
                for i in range(length):
                    for j in range(length):
                        f.write("%13.5E" %iso.CovMat[i,j])
                    f.write("\n")
        f.close()


 
def GemFile(filename, WimsData):

    format_str1 = "%16.7E%16.7E%16.7E ^ \n"
    outfile=open(filename + '/' + filename+'.gem','w')


    string1='Title: FLATTOP Uranium \n \n\
problem radiation\nmode direct\ncase eigenvalue\ngeometry spherical\nangle 3\nscatter 0\ngroups %i\nupscatter yes\nmonitor eigenvalue\n\n\
pnt p1 0.0      0.0\n\
pnt p2 6.1156   0.0\n\
pnt p3 24.1242  0.0\n\n\
line l1 p1 p2\n\
line l2 p2 p3\n\n\
divide 10 l1\n\
divide 10 l2\n\n\
region r1 l1\n\
region r2 l2\n'% WimsData.NG

    
#    string1='Title: Rowlands Pincell (case 1) \n \nproblem radiation\nmode direct\ncase eigenvalue\ngeometry xy\nangle 3\nscatter 0\ngroups %i\nupscatter yes\nmonitor eigenvalue\n\npnt p1  0.0           0.0\npnt p2  0.400205702   0.0\npnt p3  0.450231403   0.0\npnt p4  0.60          0.0\npnt p5  0.28298817    0.28298817\npnt p6  0.31836168    0.31836168\npnt p7  0.60          0.60\n\nline l1  p1 p2\nline l2  p2 p3\nline l3  p3 p4\nline l4  p1 p5\nline l5  p5 p6\nline l6  p6 p7\nline l7  p2 p5 p1\nline l8  p3 p6 p1\nline l9  p4 p7\n\ndivide 10 l1 l2 l3\ndivide 10 l4 l5 l6\ndivide 10 l7 l8 l9\n\nregion r1 l1 l7 l4\nregion r2 l2 l8 l5 l7\nregion r3 l3 l9 l6 l8\n' % WimsData.NG
    
    string2='\nproperties r1 ma1\nproperties r2 ma2\n\nboundary vacuum l2 \n\nfill\n\ndata\n\nstop'

    useTransport=True#False
    
    outfile.write(string1)
   
    for nuc in WimsData.Micro:
    
        if useTransport:
            sigmaT = nuc.SIGTRAN
            for i in range(WimsData.NG):
                nuc.scatter[0][i,i] = nuc.SIGTRAN[i] - nuc.SIGABS[i] - sum(nuc.scatter[0][i,:]) + nuc.scatter[0][i,i]
        else:
            sigmaT = nuc.SIGTOT
    
#        print [sigmaT[i] - sum(nuc.scatter[0][i,:]) - nuc.SIGABS[i] for i in range(WimsData.NG)]
    
        outfile.write("\n")
        outfile.write("material %s" %nuc.name)
        if nuc.SIGFISS:
            outfile.write(format_str1  %(sigmaT[0],nuc.SIGABS[0],nuc.SIGFISS[0]*nuc.FISSNU[0]) )
        else:
            outfile.write(format_str1  %(sigmaT[0],nuc.SIGABS[0],0.0 ) )
        GemScatter(nuc.scatter[0][0,:], outfile)
        if(WimsData.NG>1):
            for grp in range(1,WimsData.NG):
                outfile.write("             ")
                if nuc.SIGFISS:
                    outfile.write(format_str1  %(sigmaT[grp], nuc.SIGABS[grp], nuc.SIGFISS[grp]*nuc.FISSNU[grp]) )
                else:
                    outfile.write(format_str1  %(sigmaT[grp], nuc.SIGABS[grp], 0.0 ))
                GemScatter(nuc.scatter[0][grp,:], outfile)            


    outfile.write("\n")
    c=1
    for mat in WimsData.MatSpecs:
        outfile.write("mix ma%s \n" %c)
        for nd, iso in zip(mat[0],mat[1]):
            outfile.write("mix ma%s  %s  %13.6e \n" %(c,iso,nd))   
        c+=1                
                

    outfile.write("\n")
    flag=False
    for nuc in WimsData.CrossSections:
        if nuc.SIGFISS:
            GemSpectrum(nuc.FISSCHI,outfile)
            flag=False
            break
        else:
            flag=True
    if(flag):
        print 'No Fission spectrum'   
   
    outfile.write(string2)
            
    outfile.close()
 
 

def GemScatter(scat,outfile):
    numperline=3
    totalnum=len(scat)
      
    delta=[]
    for i in range( int(totalnum)/int(numperline) ):
        delta.append(numperline)    
    last=int(np.math.fmod(totalnum,numperline))
    if(last != 0):
        delta.append(last)    
    
    start=0
    finish=0
    for wvec in delta:
        finish+=wvec
        outfile.write("             ")
        for i in range(start,finish):
            outfile.write("%16.7E" %scat[i])
        outfile.write(" ^ \n")
        start+=wvec

 
def GemSpectrum(sptrm,outfile):
    numperline=4
    totalnum=len(sptrm)
      
    delta=[]
    for i in range( int(totalnum)/int(numperline) ):
        delta.append(numperline)    
    last=int(np.math.fmod(totalnum,numperline))
    if(last != 0):
        delta.append(last)    
    
    start=0
    finish=0
    outfile.write("spectrum   ")

    finish+=delta[0]
    for i in range(start,finish):
        outfile.write("%13.6e" %sptrm[i])
    outfile.write(" ^ \n")
    start+=delta[0] 

    for j in range(1,len(delta)):
        outfile.write("           ")
        finish+=delta[j]
        for i in range(start,finish):
            outfile.write("%13.6e" %sptrm[i])
        outfile.write(" ^ \n")
        start+=delta[j]         
 
 
 
def EVENT_template(filename, WimsData):

    outfile=open(filename+'/'+filename+'.template','w')
    useTransport = False
    
    for nuc in WimsData.CrossSections:

        if useTransport:
            sigmaT = nuc.SIGTRAN
            for i in range(WimsData.NG):
                nuc.scatter[0][i,i] = nuc.SIGTRAN[i] - nuc.SIGABS[i] - sum(nuc.scatter[0][i,:]) + nuc.scatter[0][i,i]
        else:
            sigmaT = nuc.SIGTOT
    
        # If not sampling output mean values
        if nuc.CovMat==None: 
        
            for i in range(WimsData.NG):
                if nuc.SIGFISS:
                    outfile.write("'%s' %14.7E %14.7E %14.7E \n"  %(nuc.name.ljust(8), sigmaT[i],nuc.SIGABS[i],nuc.SIGFISS[i]*nuc.FISSNU[i]) )
                else:
                    outfile.write("'%s' %14.7E %14.7E %14.7E \n"  %(nuc.name.ljust(8), sigmaT[i],nuc.SIGABS[i],0.0 ) )            
                WriteScatter2(nuc.name, nuc.scatter[0][i,:], i, outfile)
        
        # Print sampling tags
        else:
        
            if nuc.reactions != []:
                # Total cross section will always vary

                if len(nuc.reactions)==4:
                    for i in range(WimsData.NG):
                        outfile.write("'%s' $%s_t%i $%s_a%i $%s_nf%i\n" %(nuc.name.ljust(8),nuc.name,i+1,nuc.name,i+1,nuc.name,i+1  ))
                        WriteScatter(nuc.name, nuc.scatter[0][i,:], i, outfile)
                elif len(nuc.reactions)==3:
                    if 's' not in nuc.reactions:
                        for i in range(WimsData.NG):
                            outfile.write("'%s' $%s_t%i $%s_a%i $%s_nf%i\n" %(nuc.name.ljust(8),nuc.name,i+1,nuc.name,i+1,nuc.name,i+1  ))
                            WriteScatter2(nuc.name, nuc.scatter[0][i,:], i, outfile)
                    else:
                        for i in range(WimsData.NG):
                            outfile.write("'%s' $%s_t%i $%s_a%i $%s_nf%i\n" %(nuc.name.ljust(8),nuc.name,i+1,nuc.name,i+1,nuc.name,i+1  ))
                            WriteScatter(nuc.name, nuc.scatter[0][i,:], i, outfile)                    
                elif len(nuc.reactions)==2:
                    if 's' in nuc.reactions:
                        for i in range(WimsData.NG):
                            if 'c' in nuc.reactions:
                                outfile.write("'%s' $%s_t%i $%s_a%i %14.7E\n" %(nuc.name.ljust(8),nuc.name,i+1,nuc.name,i+1,nuc.SIGFISS[i]*nuc.FISSNU[i] ))
                            elif 'n' in nuc.reactions:
                                outfile.write("'%s' $%s_t%i %14.7E $%s_nf%i\n" %(nuc.name.ljust(8),nuc.name,i+1,nuc.SIGABS[i],nuc.name,i+1))
                            else:
                                outfile.write("'%s' $%s_t%i $%s_a%i $%s_nf%i\n" %(nuc.name.ljust(8),nuc.name,i+1,nuc.name,i+1,nuc.name,i+1  ))
                            WriteScatter(nuc.name, nuc.scatter[0][i,:], i, outfile)                    
                    else:
                        for i in range(WimsData.NG):
                            outfile.write("'%s' $%s_t%i $%s_a%i $%s_nf%i\n" %(nuc.name.ljust(8),nuc.name,i+1,nuc.name,i+1,nuc.name,i+1  ))
                            WriteScatter2(nuc.name, nuc.scatter[0][i,:], i, outfile)
                elif len(nuc.reactions)==1:
                    if 's' in nuc.reactions:
                        for i in range(WimsData.NG):
                            outfile.write("'%s' $%s_t%i %14.7E %14.7E\n" %(nuc.name.ljust(8),nuc.name,i+1,nuc.SIGABS[i],nuc.SIGFISS[i]*nuc.FISSNU[i] ))
                            WriteScatter(nuc.name, nuc.scatter[0][i,:], i, outfile)                    
                    elif 'n' in nuc.reactions:
                        for i in range(WimsData.NG):
                            outfile.write("'%s' %14.7E %14.7E $%s_nf%i\n" %(nuc.name.ljust(8),sigmaT[i],nuc.SIGABS[i],nuc.name,i+1))
                            WriteScatter2(nuc.name, nuc.scatter[0][i,:], i, outfile)
                    elif 'c' in nuc.reactions:
                        for i in range(WimsData.NG):
                            outfile.write("'%s' $%s_t%i $%s_a%i %14.7E\n" %(nuc.name.ljust(8),nuc.name,i+1,nuc.name,i+1,nuc.SIGFISS[i]*nuc.FISSNU[i] ))
                            WriteScatter2(nuc.name, nuc.scatter[0][i,:], i, outfile)
                    elif 'f' in nuc.reactions:
                        for i in range(WimsData.NG):
                            outfile.write("'%s' $%s_t%i $%s_a%i $%s_nf%i\n" %(nuc.name.ljust(8),nuc.name,i+1,nuc.name,i+1,nuc.name,i+1  ))
                            WriteScatter2(nuc.name, nuc.scatter[0][i,:], i, outfile)
                            
            else:
                print 'No reactions present in ', nuc.name
            
        
    
    
    flag=False
    for nuc in WimsData.CrossSections:
#        if nuc.CovMat==None: continue
        if nuc.SIGFISS:
            WriteSpectrum(nuc.FISSCHI,outfile)
            flag=False
            break
        else:
            flag=True
    if(flag):
        print 'No Fission spectrum - Fixed source problem??'

    outfile.close()

def WriteSpectrum(sptrm,outfile):
    numperline=4
    totalnum=len(sptrm)
      
    delta=[]
    for i in range( int(totalnum)/int(numperline) ):
        delta.append(numperline)    
    last=int(np.math.fmod(totalnum,numperline))
    if(last != 0):
        delta.append(last)  

    start=0
    finish=0
    outfile.write("'SPCTRM'   ")

    finish+=delta[0]
    for i in range(start,finish):
        outfile.write(" %13.7E " %sptrm[i])
    outfile.write("\n")
    start+=delta[0] 

    for j in range(1,len(delta)):
        outfile.write("'SPCTRM'   ")
        finish+=delta[j]
        for i in range(start,finish):
            outfile.write(" %13.7E " %sptrm[i])
        outfile.write("\n")
        start+=delta[j]        
        
def WriteScatter(name, scat, grp, outfile):
    numperline=4
    totalnum=len(scat)
      
    delta=[]
    for i in range( int(totalnum)/int(numperline) ):
        delta.append(numperline)    
    last=int(np.math.fmod(totalnum,numperline))
    if(last != 0):
        delta.append(last)    
    
    start=0
    finish=0
    for wvec in delta:
        finish+=wvec
        outfile.write("'%s'" %name.ljust(8))
        for i in range(start,finish):
            if i==grp:
                outfile.write(" $%s_s%i" %(name,grp+1))
            else:
                outfile.write(" %14.7E" %scat[i])
        outfile.write("\n")
        start+=wvec


def WriteScatter2(name, scat, grp, outfile):
    numperline=4
    totalnum=len(scat)
      
    delta=[]
    for i in range( int(totalnum)/int(numperline) ):
        delta.append(numperline)    
    last=int(np.math.fmod(totalnum,numperline))
    if(last != 0):
        delta.append(last)    
    
    start=0
    finish=0
    for wvec in delta:
        finish+=wvec
        outfile.write("'%s'" %name.ljust(8))
        for i in range(start,finish):
            outfile.write(" %14.7E" %scat[i])
        outfile.write("\n")
        start+=wvec            
