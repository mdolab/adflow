#!/usr/bin/env python 
#
# Collection of results from the CBD output and preparing polara sum-up files for plots
#
from numpy import *
from myFuncs import parLocator,stringLocator,setPolaraType,setVelDir,readParameter,setContribution
import os
import sys 
from optparse import OptionParser

#   import numarray   # numarray was found only for python 2.4, while cluster 5 has python 2.3

def main():

    parser = OptionParser()
#
# control parameters: define command line parser parameters
#
    usage = "usage: %prog [options] (enter %prog -h for a list of valid otions) "
    version = '2.04'    
    parser = OptionParser(usage=usage,version=version)
    parser.set_defaults(verbose=False)
    parser.add_option("-v","--verbose", action="store_true", dest="verbose",
                  help="Print debug messages (def=false)" )
    (options, args) = parser.parse_args()
    verbose=options.verbose 
#
# control parameters: read CBd dat file
#
    fcbd=open('.cbdDat','r')
    fcbdR=fcbd.readlines()
    ctrlFile =fcbdR[0].strip(' ');ctrlFile=ctrlFile.strip('\n');ctrlFile=ctrlFile.strip(' ');
    fcbd.close()
    if verbose:
        print 'CBD data read' 
    
#--------------- now read the control file and parse it
    
    fc=open(ctrlFile,'r')
    ctrl=fc.readlines()
    nc=size(ctrl)
    fc.close()

#------------- determine polara type

    PA,polarSweepType,nAalpha,nBeta,nPhi,alpha,beta,phi,polarVar=setPolaraType(ctrl,nc,verbose)
    nPolara = max(nAalpha,nPhi)

    print '-----PolarSweepType = '+str(polarSweepType)+' Polar sweep in '+polarVar+' using '+str(nPolara)+' angles ----'
  
 
#------------- base SUnb input file ----------------------- 
    keyWordInputbaseFile='input base file'
    iBaseInputF = parLocator(keyWordInputbaseFile,ctrl,nc,-1,verbose)
    bIFLine=ctrl[iBaseInputF]
    icol=bIFLine.index(':')
    sBIF=bIFLine[icol+1:]
    inputbaseFile=sBIF.strip(' ')
    inputbaseFile=inputbaseFile.strip('\n')
    if verbose:
        print '--------------   Using base MB input file: '+inputbaseFile+' ------------'

# Identify the base surface for cxBase
#---------------------------------------- 
    keyWordBaseSur='base surface'
    iBaseSur= parLocator(keyWordBaseSur,ctrl,nc,-1,verbose)
    if iBaseSur == -1:
        baseSurName='notDefinedInInput'
    else:
        icol=ctrl[iBaseSur].index(':')
        baseSurName=ctrl[iBaseSur][icol+1:-1].strip(' ').lower()
    if verbose:
        print 'base Surface Name ->'+baseSurName+'<-'

# Seconday loads section
#------------------------------
    keyWordCPR='components polara required'
    iCPR = parLocator(keyWordCPR,ctrl,nc,-1,verbose)
    if iCPR == -1:
        print 'No component was designated for secondary polara'
        nCPR=0
    else:
        icol=ctrl[iCPR].index(':')
        CPR_line=ctrl[iCPR][icol+1:-1]
        CPRdataRaw=CPR_line.split(',')
        CPRdata=[x.strip(' ').lower() for x in CPRdataRaw]
        nCPR=size(CPRdata)
        print '---'+str(nCPR)+' components were designated for secondary polar sweep '
        if verbose:
            for k in range(0,nCPR):
                print 'CPR '+str(k)+': ->'+CPRdata[k]+'<-'
    
#--------- end of reading from control file
#------- some administrative issues like file name and keyWord
    
    tt=inputbaseFile.split('.')
    baseName=tt[0]

#--------- read the base input file
    
    f=open(inputbaseFile,'r')
    b=f.readlines()
    n=size(b)
    f.close()

# open the destination polara file : coeff file
#------------------------------------------------
    coeffFileName=baseName+'PolaraCoeff.dat'
    foc=open(coeffFileName,'w')
    foc.write('%  \n%         Main coefficients of a polara run  \n%  \n')
#
# locate reference values from the base input file
#
# geometrical data
#
    keyWordAr='reference surface'
    Ar,ipr  =  readParameter(b,n,keyWordAr,-1,verbose)
    keyWordLr='reference length'
    Lr,ipr  =  readParameter(b,n,keyWordLr,-1,verbose)
    keyWordXr='moment reference point x'
    Xr,ipr  =  readParameter(b,n,keyWordXr,-1,verbose)
    keyWordYr='moment reference point y'
    Yr,ipr  =  readParameter(b,n,keyWordYr,-1,verbose)
    keyWordZr='moment reference point z'
    Zr,ipr  =  readParameter(b,n,keyWordZr,-1,verbose)

    sRefVal=[Ar,Lr]
    refVal=map(float,sRefVal)
    sRefP=[Xr,Yr,Zr]
    refP=map(float,sRefP)
#
# physical data
#
# Non-dimensional numbers
#
    keyWord='mach for coefficients'
    iparMcoeff = parLocator(keyWord,b,n,-1,verbose)
    keyWord='mach'
    MachNum,ipr=readParameter(b,n,keyWord,iparMcoeff,verbose)  # look for Mach, but avoid Mach for coefficients
    keyWord='reynolds length (in meter)'
    ReNumRefLength,iprDRe=readParameter(b,n,keyWord,-1,verbose)
    keyWord='reynolds'
    ReNum,ipr=readParameter(b,n,keyWord,iprDRe,verbose)

    sNonDimNum=[MachNum,ReNum,ReNumRefLength]
    nonDimNum=map(float,sNonDimNum)
#
# Dimensional flow conditions
#
    keyWord='Reference pressure (in Pa)'
    pRef,ipr=readParameter(b,n,keyWord,-1,verbose)
    keyWord='Reference density (in kg/m^3)'
    rhoRef,ipr=readParameter(b,n,keyWord,-1,verbose)
    keyWord='Reference temperature (in K)'
    TRef,ipr=readParameter(b,n,keyWord,-1,verbose)

    sDimCond=[pRef,rhoRef,TRef]
    dimCond=map(float,sDimCond)

#
# Thermodynamic properties
#
    keyWord='Constant specific heat ratio'
    gamma,ipr=readParameter(b,n,keyWord,-1,verbose)
    keyWord='Gas constant (J/(kg K))'
    rGas,ipr=readParameter(b,n,keyWord,-1,verbose)
    keyWord='Free stream temperature (in K)'
    TFreeS,ipr=readParameter(b,n,keyWord,-1,verbose)

    sThermoPar=[gamma,rGas,TFreeS]
    thermoPar=map(float,sThermoPar)
#
# mesh file
#
    keyWord='grid file'
    gridFile,ipr=readParameter(b,n,keyWord,-1,verbose)
    gridFile=gridFile.strip(' ')
    gridFile=gridFile.strip('\n')

#
# Non contributing components
#
    keyWord='Contribute to forces'
    ipr=1
    ipro=0
    numOfNonContributingParts=0
    listOfNonContributingParts=[]
    bNext=b;
    secMarkedForRemoval=[]  # prepare a list of component indices identified for removal
    while ipr > -1:
        nNext=size(bNext)
        partName,removeContribution,ipr=setContribution(bNext,nNext,keyWord,-1,verbose)
        ipro=ipr+ipro
        bNext=b[ipro:]
        if ipr > -1:
            if removeContribution == True:
                numOfNonContributingParts=numOfNonContributingParts+1
                listOfNonContributingParts.append(partName)
                
    if numOfNonContributingParts > 0:

#  now chexk if components that are assigned for contribution removal are part of component print out.
#  add if required
        cprAdded = False
        if nCPR == 0:
            CPRdata=[x for x in listOfNonContributingParts]
            nCPR=numOfNonContributingParts
            cprAdded = True
            secMarkedForRemoval=range(nCPR)
        else:
            for i in range(numOfNonContributingParts):
                try:
                    secMarkedForRemoval.append(CPRdata.index(listOfNonContributingParts[i].lower()))
#   do nothing, since already the contribution of this component is assigned to be listed
                except ValueError:
#   oops: the contribution of this component is NOT assigned to be listed
                    CPRdata.append(listOfNonContributingParts[i].lower())
                    nCPR=nCPR+1
                    secMarkedForRemoval.append(nCPR-1)
                    cprAdded = True
        if cprAdded:
            print 'List of componenents for separated contribution increased to '+str(nCPR)+' components (non-contribution check)'

                

# write down sweep variable name and single angle name

    if polarSweepType == 1:
        sweepAngleName='alpha'; singleAngleName='phi'; singleAngle=phi[0]
    elif polarSweepType == 2:
        sweepAngleName='alpha'; singleAngleName='beta' ; singleAngle = beta[0]
    else:
        sweepAngleName='phi' ; singleAngleName='alpha' ; singleAngle = alpha[0]

    polaraHeader = '  polarSweepType : %i  sweepAngle : %s  singleAngle : %s = %.1f deg \n'%(polarSweepType,sweepAngleName,singleAngleName,singleAngle)

    foc.write('% '+polaraHeader)
#
# write down the reference values
#
    refLine1=    '  Reference point for moments   : [ %10.5f , %10.5f, %10.5f ] '%(refP[0],refP[1],refP[2])+' % [m] \n'
    refLine2=    '  Reference area and length: Ar : %14.7f , Lr : %10.5f    '%(refVal[0],refVal[1])+' % [m] \n'
    noDLine1=    '  Mach                          : %7.2f ,    '%(nonDimNum[0])
    noDLine2=    '  Reynolds Number               : %10.3e , '%(nonDimNum[1])
    noDLine3=    '  Reynolds length :  %12.5f '%(nonDimNum[2])+' % [m] \n'
    noDLineAll = noDLine1+noDLine2+noDLine3
    dimCondLine1='  Reference pressure (in Pa)    : %10.2f , '%(dimCond[0])
    dimCondLine2='  Reference density (in kg/m^3) : %10.5f ,  '%(dimCond[1])
    dimCondLine3='  Reference temperature (in K) : %9.2f  \n'%(dimCond[2])
    dimCondLineAll=dimCondLine1+dimCondLine2+dimCondLine3
    thermoLine1 ='  Constant specific heat ratio  : %6.2f,  '%(thermoPar[0])
    thermoLine2 ='      Gas constant (J/(kg K))       :  %10.2f,'%(thermoPar[1])
    thermoLine3 ='  Free stream temperature (in K) : %9.2f  \n'%(thermoPar[2])
    thermoLineAll=thermoLine1+thermoLine2+thermoLine3
    gridLine    ='  Grid file                     :  %s \n'%(gridFile)

    foc.write('% \n')
    foc.write('% ========================= Reference parameters ===============================================================\n')
    foc.write('% \n')
    foc.write('% '+refLine1)
    foc.write('% '+refLine2)
    foc.write('% '+noDLineAll)
    foc.write('% '+dimCondLineAll)
    foc.write('% '+thermoLineAll)
    foc.write('% '+gridLine)
    foc.write('% \n')
    foc.write('% ========================= Header lines from input file ======================================================\n')
    foc.write('% \n')

# now write in the first 5 lines of the base input file, assuming they include some header
    for j in range(0,8):
        foc.write('% '+b[j])
        
#  now write in the header line for the polara


    tableHeader='%  '+sweepAngleName+'     cFx          cFy            cFz         cMx          cMy          cMz         Cxbase '
    tableHeader=tableHeader+'       CQ '
    foc.write(tableHeader+'\n')
    foc.write('% -------------------------___------------------------------------------------------------------------------\n')

#   Secondary components polar sweep files
#-----------------------------------------------
    if nCPR > 0 :
        scpf=[]
        for k in range(0,nCPR):
            scpFileNme=baseName+'_secPolarCoeff_'+CPRdata[k]+'.dat'
            fncpf='fncpf'+str(k)
            scpf.append(fncpf)
            scpf[k]=open(scpFileNme,'w')
#  now write in the header line for the secondary polara
            scpfHeader = '     Secondary polar sweep coefficients for %s  \n' % (CPRdata[k])
            scpfHeader='% \n% '+scpfHeader+'% \n'
            scpf[k].write(scpfHeader)
            scpf[k].write('% '+polaraHeader)
            scpf[k].write('% \n')
            scpf[k].write('% ========================= Reference parameters ==+============================================================\n')
            scpf[k].write('% \n')
            scpf[k].write('% '+refLine1)
            scpf[k].write('% '+refLine2)
            scpf[k].write('% '+noDLineAll)
            scpf[k].write('% '+dimCondLineAll)
            scpf[k].write('% '+thermoLineAll)
            scpf[k].write('% '+gridLine)
            scpf[k].write('% \n')
            scpf[k].write('% ========================= Header lines from input file ==============================================================\n')
            scpf[k].write('% \n')
# now write in the first 5 lines of the base input file, assuming they include some header
            for j in range(0,8):
                scpf[k].write('% '+b[j])
            tableHeader='%  '+sweepAngleName+'     cFx          cFy            cFz         cMx          cMy          cMz  '
            scpf[k].write(tableHeader+'\n')
            scpf[k].write('% -------------------------------------------------------------------------------------------\n')


#--- A for loop over polara cases
#    Prepare a run input file and a cbd file for each case
    
    suffix=[]

    nPolaraActual = 0            # the actual number of cases whose output exists
    mainCoeffData = []           # this set carries all the main polar sweep data lines
    removeCoeffData = []        # this set carries the contribution of components that should be removed

# start main loop over the polar sweep angles

    for j in range(0,nPolara):
        if j > 0 :
            angle_jm1=angle

        if polarSweepType == 1 :
            alpha_case=alpha[j] ; phi_case=phi[0]
            angle=alpha_case; keyw='a'
        elif polarSweepType == 2 :
            alpha_case=alpha[j] ; beta_case=beta[0] 
            angle=alpha_case; keyw='a'
        else:
            alpha_case=alpha[0] ; phi_case=phi[j]
            angle= phi_case  ; keyw='phi'
            
#
#------------- now find out if the given angle is integer or has a fractin residue 
#
        residue=abs(angle)-floor(abs(angle))
        if residue == 0 :
            sSufifix       = '%s%.0f' % (keyw,angle)
        else:
            if angle > 0 :
                angleName=floor(angle)
            else:
                angleName=ceil(angle)
            sSufifix       = '%s%.0fp%.0f' % (keyw,angleName,10.0*residue)

        suffix.append(sSufifix)
        cbdOutput=baseName+'_'+suffix[j]+'_CBDOUT.dat'

        try:
            fd=open(cbdOutput,'r');
            d=fd.readlines()
            nd=size(d)
            fd.close()

            CBDexists = 1
            nPolaraActual = nPolaraActual + 1

            if verbose:
                if polarSweepType == 1 :
                    st1 = ' Polara input file No %i; alpha = %.1f deg;  phi =  %.1f deg \n' %(j+1,alpha_case,phi_case)
                elif polarSweepType == 2 :
                    st1 = ' Polara input file No %i; alpha = %.1f deg;  beta =  %.1f deg \n' %(j+1,alpha_case,beta_case)
                else:
                    st1 = ' Polara input file No %i; alpha = %.1f deg;  phi =  %.1f deg \n' %(j+1,alpha_case,phi_case)
                print st1

        except IOError:
            st1 = ' Polara input file No %i; alpha = %.1f deg;  phi =  %.1f deg is missing. loop goes on \n' %(j+1,alpha_case,phi_case)
            print st1
            CBDexists = 0
#            break

#------ determine convergence quality

        if CBDexists == 1:
            keyWord='convergenceQuality'
            sConvergenceQuality,icq=readParameter(d,nd,keyWord,-1,verbose)
            if icq == -1:
                convergenceQuality=-1
            else:
                convergenceQuality=map(int,[sConvergenceQuality])[0]

            if verbose:
                ttpr= ' file No %i; angle = %.1f deg;  convergenceQuality = %i\n' %(j+1,angle,convergenceQuality)
                print ttpr
            
 #----  now search for the total coefficients line
            keyTotalCoefLine='#   forces/moments sum (total)'
            iTotC = stringLocator(keyTotalCoefLine,d,nd,verbose)
            if iTotC == -1:
                print 'ERROR: Something wrong with the CBDOUT file of angle = '+str(angle)+' Failed to find the total coefficients line '
                raise SystemExit('You can remove file '+cbdOutput+' and rerun polaraPlotPrep without it')
            
            iTotC_Data=iTotC+4
# locate results of base surface for cxbase
            iBaseRes=stringLocator(baseSurName,d,nd,verbose)
            if iBaseRes == -1:
                cxBase=0.0   # such surface results not found
            else:
                sCxBase=d[iBaseRes+4].split()[0]
                cxBase=float(sCxBase)

            dataLine='  %.1f   %s     %.5f      %i \n' % (angle,d[iTotC_Data][:-1],cxBase,convergenceQuality)
            mainCoeffData.append(dataLine)
#
# -------- now look at seconday coefficients sweeps
#
           
            if nCPR > 0 :
                for k in range(0,nCPR):
                    iCPR_C=stringLocator(CPRdata[k],d,nd,verbose) 
                    iCPR_C_data=iCPR_C+4    # go the the data line of component k in the CBDOUT file
                    dataLine='  %.1f   %s  ' % (angle,d[iCPR_C_data])
                    scpf[k].write(dataLine)

            removeLine=[0.0,0.0,0.0,0.0,0.0,0.0]
            if  numOfNonContributingParts > 0:
                for i in range(numOfNonContributingParts):
                    k=secMarkedForRemoval[i]
                    iCPR_C=stringLocator(CPRdata[k],d,nd,verbose)
                    iCPR_C_data=iCPR_C+4    # go the the data line of component k in the CBDOUT file
#                                remove blanks, strip leading and trailing blanks and then replace single blanks with commas
                    modLine=d[iCPR_C_data].replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').strip().replace(' ',', ').split(',')
                    smodLineK=map(float,modLine)
                    removeLine=[x+y for x,y in zip(removeLine,smodLineK)]
                            
                removeCoeffData.append(removeLine)
                            
# end of j loop on lines

    fd.close()
    if verbose: 
        print 'file '+cbdOutput+' closed '

#
# Now remove the contribution of thoes components that were marked to be non-contributing
#
    if numOfNonContributingParts > 0:
        for j in range(0,nPolaraActual):
            mainLine= mainCoeffData[j].replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').strip().replace(' ',',').split(',')

            try:
                mainLineJ=map(float,mainLine)
            except ValueError:
                print 'Large blanks in primary data are creating problems in manipulating. I try to fix that '
                mainLine= mainCoeffData[j].replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').\
                replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').replace('  ',' ').strip().replace(' ',',').split(',')
                mainLineJ=map(float,mainLine)

            mainLineJ[1:7]=[x - y for x, y in zip(mainLineJ[1:7],removeCoeffData[j])]
            mainCoeffData[j]= '  %.1f     %.5f       %.5f        %.5f       %.5f      %.5f      %.5f      %.5f     %i \n' % (mainLineJ[0],mainLineJ[1],mainLineJ[2],
                                                                                               mainLineJ[3],mainLineJ[4],mainLineJ[5],
                                                                                               mainLineJ[6],mainLineJ[7],mainLineJ[8])
                
    for j in range(0,nPolaraActual):
        foc.write(mainCoeffData[j])
    foc.close()
    print 'Polara file '+coeffFileName+' completed and closed with '+str(nPolaraActual)+' data lines'

    if nCPR > 0 :
        for k in range(0,nCPR):
            scpf[k].close()
            scpFileNme=baseName+'_secPolarCoeff_'+CPRdata[k]+'.dat'
            print 'Secondary Polara file '+scpFileNme+' completed and closed '
            

if __name__ == "__main__":
    main()

    

