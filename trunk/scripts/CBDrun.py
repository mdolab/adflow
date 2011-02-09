#!/usr/bin/env python 

from numpy import *
from myFuncs import parLocator,setPolaraType,setVelDir 
import os
import sys 
from optparse import OptionParser

#   import numarray   # numarray was found only for python 2.4, while cluster 5 has python 2.3

def main():

    parser = OptionParser()
#
# control parameters: read CBd dat file
#
    if os.path.exists('.cbdDat'):
        fcbd=open('.cbdDat','r')
        fcbdR=fcbd.readlines()
        ctrlFile =fcbdR[0].strip(' ');ctrlFile=ctrlFile.strip('\n');ctrlFile=ctrlFile.strip(' ');
        nprocs=int(fcbdR[1]); nprocsCBD=min(nprocs,6) ;
        machinefileName=fcbdR[2].strip(' '); machinefileName=machinefileName.strip('\n') 
        machinefileName=machinefileName.strip(' ')
        envpar=fcbdR[3].strip('\n')
        if 'local' in envpar:
            envpar=envpar.strip(' ')

        fcbd.close()

    else:
        raise SystemExit('ERROR: CBD run requires .cbdDat file (created by polaraRun). Missing')

    print 'CBD data read' 
#
# control parameters: define command line parser parameters
#
    usage = "usage: %prog [options] (enter %prog -h for a list of valid otions) "
    version = '1.09'
    parser = OptionParser(usage=usage,version=version)
    parser.set_defaults(filename=ctrlFile)
    parser.set_defaults(verbose=False)
    parser.set_defaults(nprocsCBDR=8)
    parser.set_defaults(machinefilename=machinefileName)
    parser.set_defaults(env=envpar)
    
    parser.add_option("-f", "--file",action="store", type="string", dest="filename",
                  help="reads polara control parameters from FILE (default: runPolara control file) ", 
                  metavar="FILE")
    parser.add_option("-v","--verbose", action="store_true", dest="verbose",
                  help="Print debug messages (def=false)" )
    parser.add_option("-n",'--nprocsCBD', type="int", help='Number of MPI procs for CBD (def=8)' ,
                  dest="nprocsCBDR")
    parser.add_option("-m", "--machinefile",action="store", type="string", dest="machinefilename",
                  help="machinefile name (default: machinefile of polaraRun) ", 
                  metavar="MACHINEFILE")
    parser.add_option("-e","--env",action="store", type="string", dest="env",
                  help="MPI parameters, such as device (def: env of polaraRun) ", 
                  metavar="ENV")

    (options, args) = parser.parse_args()

    ctrlFile=options.filename 
    verbose=options.verbose 
    nprocsCBDR=options.nprocsCBDR
    nprocsCBD=min(nprocs,nprocsCBDR)
    machinefileName=options.machinefilename
    envpar=options.env

#--- check if machinefile includes directory name. If not, add local directory to it
    try:
        islash=machinefileName.index('/')    
        pass   # do nothing. directory is included
        
    except ValueError:
        wd=os.getcwd()
        machinefileName=wd+'/'+machinefileName
        

#--------------- now read the parameters control file and parse it
    
    fc=open(ctrlFile,'r')
    ctrl=fc.readlines()
    nc=size(ctrl)
    fc.close()

#------------- determine polara type

    PA,polarSweepType,nAalpha,nBeta,nPhi,alpha,beta,phi,polarVar=setPolaraType(ctrl,nc,verbose)
    nPolara = max(nAalpha,nPhi)

    print '-----PolarSweepType = '+str(polarSweepType)+' Polar sweep in '+polarVar+' using '+str(nPolara)+' angles ----'
    

#------------- base SUnb input file ----------------------- 
    inputbaseFileString='input base file'
    keyWordInputbaseFile=inputbaseFileString.lower()
    iBaseInputF = parLocator(keyWordInputbaseFile,ctrl,nc,-1,verbose)
    bIFLine=ctrl[iBaseInputF]
    icol=bIFLine.index(':')
    sBIF=bIFLine[icol+1:]
    inputbaseFile=sBIF.strip(' ')
    inputbaseFile=inputbaseFile.strip('\n')
    print '--------------   Using base MB input file: '+inputbaseFile+' ------------'

#--------- end of reading from control file
#------- some administrative issues like file name and keyWord
    
    tt=inputbaseFile.split('.')
    baseName=tt[0]

#-------- set flow direction vector
    
    d2r=pi/180
 #   alphar=numarray.multiply(alpha,d2r) # numarray was found only for python 2.4, while cluster 5 has python 2.3
    alphar  = [d2r*x for x in alpha]
    phir    = [d2r*x for x in phi]
    betar   = [d2r*x for x in beta]

    dv1,dv2,dv3=setVelDir(polarSweepType,PA,alphar,phir,betar)

#--------- read the base input file
    
    f=open(inputbaseFile,'r')
    b=f.readlines()
    n=size(b)
    f.close()

#----------- Locate keywords for modifications in input file
#            In case a line is not present, the relevant line will be added
#            at the end of the input file

    keyWordAoA='free stream velocity direction'
    iFocus = parLocator(keyWordAoA,b,n,-1,verbose)
    if iFocus == -1:
        b.append(' '); iFocus=n ;  n=n+1;
    keyWordSur='surface solution file'
    iSurSol = parLocator(keyWordSur,b,n,-1,verbose)
    if  iSurSol== -1:
        b.append(' '); iSurSol=n ;  n=n+1;
    keyWordSol='solution file'
    iSol = parLocator(keyWordSol,b,n,iSurSol,verbose)
    if iSol == -1:
        b.append(' '); iSol=n ;  n=n+1;
    keyWordRestartFile='restart file'
    iRestartFile = parLocator(keyWordRestartFile,b,n,-1,verbose)
    if iRestartFile == -1:
        b.append(' '); iRestartFile=n ;  n=n+1;
    keyWordRestart='restart'
    iRestart = parLocator(keyWordRestart,b,n,iRestartFile,verbose)
    if iRestart == -1:
        b.append(' '); iRestart=n ;  n=n+1;
    keyWordCBD = 'components break down'
    iCBD =  parLocator(keyWordCBD,b,n,iRestartFile,verbose)
    if iCBD == -1:
        b.append(' '); iCBD=n ;  n=n+1;
    
#--- A for loop over polara cases
#    Prepare a run input file and a cbd file for each case
    
    bCBD=[x for x in b]
    suffix=[]
    for j in range(0,nPolara):
        if j > 0 :
            angle_jm1=angle

        if polarSweepType == 1 :
            alpha_case=alpha[j] ; phi_case=phi[0]
            angle=alpha_case; keyw='a'
            st1 = ' Polara input file No %i; alpha = %.1f deg;  phi =  %.1f deg \n' %(j+1,alpha_case,phi_case)
            parFocLine ='     %s : %.2f  %.5f %.5f  # alpha = %.1f deg phi = %.1f deg \n'% (keyWordAoA,dv1,dv2[j],dv3[j],alpha_case,phi_case)
        elif polarSweepType == 2 :
            alpha_case=alpha[j] ; beta_case=beta[0] 
            angle=alpha_case; keyw='a'
            st1 = ' Polara input file No %i; alpha = %.1f deg;  beta =  %.1f deg \n' %(j+1,alpha_case,beta_case)
            parFocLine ='     %s : %.2f  %.5f %.5f  # alpha = %.1f deg beta = %.1f deg \n'% (keyWordAoA,dv1,dv2[j],dv3[j],alpha_case,beta_case)
        else:
            alpha_case=alpha[0] ; phi_case=phi[j]
            angle= phi_case  ; keyw='phi'
            st1 = ' Polara input file No %i; alpha = %.1f deg;  phi =  %.1f deg \n' %(j+1,alpha_case,phi_case)
            parFocLine ='     %s : %.2f  %.5f %.5f  # alpha = %.1f deg phi = %.1f deg \n'% (keyWordAoA,dv1,dv2[j],dv3[j],alpha_case,phi_case)
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
        
        outCbdFileName    = baseName+'_'+suffix[j]+'.cbd' 
        
        
        fd=open(outCbdFileName,'w');
        
        fd.write(' \n')
        fd.write('=================================================== \n ')
        fd.write(' \n')
        fd.write(st1)
        fd.write(' \n')
        fd.write('=================================================== \n ')
        fd.write(' \n')
        
# insert the required modificatins to input file
        
        bCBD[iFocus]=parFocLine
        solLine=      '                     Solution file : Sol.'+suffix[j]+'.cgns  \n '
        surLine=      '             Surface solution file : Sol_surface.'+suffix[j]+'.cgns  \n '
        bCBD[iSol] = solLine
        bCBD[iSurSol] = surLine
        
        RestartLineCBD  = '                     Restart file  : Sol.'+suffix[j]+'.cgns  \n '
        bCBD[iRestartFile] = RestartLineCBD
        bCBD[iRestart] ='                           Restart: yes \n'
        bCBD[iCBD]     ='                 components break down : yes \n '

# --------- now write down the input file

        for i in range(0, n):
            fd.write(str(bCBD[i]))
        
        fd.close()
        if verbose: 
            print 'file '+outCbdFileName+' closed '
 #
# find out which version of mb we intend to use
#options: adl version : if env var USE_ADL=Y
#         Previous version (BDF, ver 180) else
    useADL=os.getenv('USE_ADL')

    if useADL == 'Y':
        print 'running ADL version'
        codeName='/usr/local/mb/mbadl'
    else:
        print 'running BDF (180) version'
        codeName='/usr/local/mb/mb'
           
#-------now run the polara (CBD mode)

    for j in range(0,nPolara):
        outCbdFileName    = baseName+'_'+suffix[j]+'.cbd'
        cbdOutput=baseName+'_'+suffix[j]+'_CBDOUT.dat'

        if envpar =='local':
            cbdCommanLine = 'mpiexec -machinefile %s -n %i %s %s >  %s \n ' % (machinefileName,nprocsCBD,codeName,outCbdFileName ,cbdOutput)

        else:
            cbdCommanLine = 'mpiexec -machinefile %s -n %i -env %s %s %s >  %s \n ' % (machinefileName,nprocsCBD,envpar,codeName,outCbdFileName ,cbdOutput)

        print cbdCommanLine
        os.system(cbdCommanLine)
        
# now  collect the results onto a polara output file
    os.system('polaraPlotPrep.py')

if __name__ == "__main__":
    main()

    

