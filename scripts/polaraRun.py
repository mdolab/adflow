#!/usr/bin/env python 

from numpy import *
from myFuncs import parLocator, readList,setPolaraType,setVelDir,processAddAngle, updatedControlFile 
import os
import sys 
from optparse import OptionParser
#  import subprocess  # missing in python 2.3 in cluster5 
#   import numarray   # numarray was found only for python 2.4, while cluster 5 has python 2.3

def main():

    debug = True
    parser = OptionParser()

#
# control parameters: define command line parser parameters
#
    usage = "usage: %prog [options] (enter %prog -h for a list of valid otions) "
    version = '2.03'
    parser = OptionParser(usage=usage,version=version)
    parser.set_defaults(filename='polaraCtrl.dat')
    parser.set_defaults(verbose=False)
    parser.set_defaults(nprocs=32)
    parser.set_defaults(machinefilename='machinefile')
    parser.set_defaults(env='I_MPI_DEVICE rdma')
    parser.set_defaults(cbd=False)
    parser.set_defaults(oldVersion=False)
    parser.set_defaults(addRunStr='')

    parser.add_option("-f", "--file",action="store", type="string", dest="filename",
                  help="reads polara control parameters from FILE (default:polaraCtrl.dat) ", 
                  metavar="FILE")
    parser.add_option("-v","--verbose", action="store_true", dest="verbose",
                  help="Print debug messages (def=false)" )
    parser.add_option("-n",'--nprocs', type="int", help='Number of MPI procs (def=32)' ,
                  dest="nprocs")
    parser.add_option("-m", "--machinefile",action="store", type="string", dest="machinefilename",
                  help="machinefile name (default: machinefile in current directory) ", 
                  metavar="MACHINEFILE")
    parser.add_option("-e","--env",action="store", type="string", dest="env",
                  help="MPI parameters, such as device (def: rdma, local : ignored) ", 
                  metavar="ENV")
    parser.add_option("-c","--cbd", action="store_true", dest="cbd",
                  help="Perform Component Break Down (def=false)" )
    parser.add_option("-o","--old", action="store_true", dest="oldVersion",
                  help="Use older (previous) version of solver (def=false)" )
    parser.add_option("-a","--add", action="store",type="string",  dest="addRunStr",
                  help="List of additional/Rerun angles to be added/to improve convergence of already completed list (a1,a2,...)" )

    (options, args) = parser.parse_args()

    ctrlFile=options.filename 
    verbose=options.verbose 
    nprocs=options.nprocs
    nprocsCBD=nprocs
    machinefileName=options.machinefilename
    envpar=options.env
    cbd=options.cbd
    oldVersion=options.oldVersion
    addRunStr=options.addRunStr

    if addRunStr == '':
        print 'no Rerun '
        reRunMode=False
    else:
        reRunMode=True
        
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

#------------- base SUmb input file ----------------------- 
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
    keyWordGenCBD = 'generate cbd-output file';
    iGenCBDOUT =  parLocator(keyWordGenCBD,b,n,iRestartFile,verbose);
    if iGenCBDOUT == -1:
        b.append(' ');  iGenCBDOUT=n; n=n+1;
        
#
# find out which version of mb we intend to use
#options: adl version : if env var USE_ADL=Y
#         Previous version (BDF, ver 180) else
    useADL=os.getenv('USE_ADL')

    if useADL == 'Y':
       
        if oldVersion:
            codeName='/usr/local/mb/mbadl_previousVersion'
            print 'running ADL Previous (old) version'
        else:
            codeName='/usr/local/mb/mbadl'
            print 'running ADL version'
    else:
        print 'running BDF (180) version'
        codeName='/usr/local/mb/mb'

#
# ------ add/rerun angles 
#
    nPolaraOrig=nPolara
    if reRunMode:
        print 'running in rerun mode '
        if polarSweepType < 3:
            nPolara,alpha,computeCase,rerunCase=processAddAngle(addRunStr,nPolara,alpha)
        else:
            nPolara,phi,computeCase,rerunCase=processAddAngle(addRunStr,nPolara,phi)

#-- Since number of anges hae been modified, recompute flow directions

        alphar  = [d2r*x for x in alpha]
        phir    = [d2r*x for x in phi]
        betar   = [d2r*x for x in beta]

        dv1,dv2,dv3=setVelDir(polarSweepType,PA,alphar,phir,betar)

#--- Prepare an updated control file
        if nPolara > nPolaraOrig:
            if polarSweepType < 3:
                updatedControlFile(ctrl,nc,alpha,ctrlFile,verbose)
            else:
                updatedControlFile(ctrl,nc,phi,ctrlFile,verbose)
    
    else:
#--------  compute all cases in control file
        print 'running in regular mode '
        computeCase=[True for j in range(0,nPolara)]
        rerunCase=[False for j in range(0,nPolara)]


#--- A for loop over polara cases
#    Prepare a run input file and a cbd file for each case
    suffix=[]
    for j in range(0,nPolara):
        if j > 0 :
            angle_jm1=angle
        
        if polarSweepType == 1 :
            alpha_case=alpha[j] ; phi_case=phi[0]
            angle=alpha_case; keyw='a'
            st1 = ' Polara input file No %i; alpha = %.1f deg;  phi =  %.1f deg \n' %(j+1,alpha_case,phi_case)
        elif polarSweepType == 2 :
            alpha_case=alpha[j] ; beta_case=beta[0] 
            angle=alpha_case; keyw='a'
            st1 = ' Polara input file No %i; alpha = %.1f deg;  beta =  %.1f deg \n' %(j+1,alpha_case,beta_case)
        else:
            alpha_case=alpha[0] ; phi_case=phi[j]
            angle= phi_case  ; keyw='phi'
            st1 = ' Polara input file No %i; alpha = %.1f deg;  phi =  %.1f deg \n' %(j+1,alpha_case,phi_case)
#
#------------- now find out if the given angle is integer or has a fraction residue 
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
        outFileName    = baseName+'_'+suffix[j]+'.inp' 
        if polarSweepType == 2 :
            parFocLine ='     %s : %.2f  %.5f %.5f  # alpha = %.1f deg beta = %.1f deg \n'% (keyWordAoA,dv1,dv2[j],dv3[j],alpha_case,beta_case)
        else:
            parFocLine ='     %s : %.2f  %.5f %.5f  # alpha = %.1f deg phi = %.1f deg \n'% (keyWordAoA,dv1,dv2[j],dv3[j],alpha_case,phi_case)


        if computeCase[j]:

            fo=open(outFileName,'w'); 
       
            fo.write(' \n')
            fo.write('=================================================== \n ')
            fo.write(' \n')
            fo.write(st1)
            fo.write(' \n')
            fo.write('=================================================== \n ')
            fo.write(' \n')


# insert the required modificatins to input file
        
            b[iFocus]=parFocLine; 
            solLine=      '                     Solution file : Sol.'+suffix[j]+'.cgns  \n '
            surLine=      '             Surface solution file : Sol_surface.'+suffix[j]+'.cgns  \n '
            b[iSol] = solLine  ; 
            b[iSurSol] = surLine;

            if rerunCase[j]:
#----  this is a rerun of previous cases for convergence improvement

                RestartFileName = 'Sol.'+suffix[j]+'.cgns'
# -----check that the solution file of present angle exists
                if os.path.exists(RestartFileName):
                    print 'Rerun case: '+suffix[j]
                else:
                    raise SystemExit('ERROR in rerun case '+suffix[j]+': file '+RestartFileName+' is missing')
            
                RestartLine = '                     Restart file  :  '+RestartFileName+' \n'
                b[iRestart] ='                           Restart: yes \n'
                b[iRestartFile] = RestartLine

            else:
        
                if j > 0:
                    b[iRestart] ='                           Restart: yes \n' 
#
                    RestartFileName = 'Sol.'+suffix[j-1]+'.cgns'
                    RestartLine = '                     Restart file  :  '+RestartFileName+' \n'
                    b[iRestartFile] = RestartLine

            b[iGenCBDOUT] = '   Generate CBD-output file : yes \n' ;
        
# --------- now write down the input file

            for i in range(0, n):
                fo.write(str(b[i]))
        
            fo.close() ; 
            if verbose: 
                print 'file '+outFileName+' closed '
        
#-------now run the polara (computation mode)
# start by cleaning remanants
    
    if os.path.exists('typescript'):
        os.system('rm typescript')

    for j in range(0,nPolara):
        if computeCase[j]:
            inputfileFile=baseName+'_'+suffix[j]+'.inp'

            if envpar =='local':
                commanLine='"mpiexec -machinefile '+machinefileName+' -n '+str(nprocs)+' '+codeName+' '+inputfileFile+' "'
            else:
                commanLine='"mpiexec -machinefile '+machinefileName+' -n '+str(nprocs)+' -env '+envpar+' '+codeName+' '+inputfileFile+' "'
            
            print commanLine
            os.system('script -c '+commanLine)
#
#-----  try to verify that no error occured during present run
#
            if debug == False :
                abort = False
                wc_cmd=os.popen('wc typescript')
                wc_str=wc_cmd.read() ; wc_set=wc_str.split() 
                if verbose:
                    print 'check typescript length: '+wc_str+' first term: '+str(wc_set[0])
                if wc_set[0].isdigit():
                    wci=int(wc_set[0])
                if wci < 40 :
                    abort=True
                else:
                    abort=False
                if abort :
                    raise SystemExit('ERROR in case no '+str(j)+': typescript too short or missing. Abort run ')
        
#  all is well (Hopefully... )
            typescriptName='typescript_'+suffix[j] 
            os.system('mv typescript '+typescriptName)
        
#------ prepare a data control file for CBD run

            fcbd=open('.cbdDat','w')
            fcbd.write(ctrlFile+' \n') ; 
            fcbd.write(str(nprocs)+' \n')
            fcbd.write(machinefileName+' \n')
            fcbd.write(envpar+' \n')
            fcbd.close()

#-------now run the polara (CBD mode)

# now  collect the resultsfrom the CBDOUT files  onto a polara output file
    os.system('polaraPlotPrep.py')

# direct CBD
    if cbd:
        os.system('CBDrun.py')
        pass # for dPolara
    else:
        print 'for direct Component Break Down sweep, enter CBDrun.py (no parameter needed)'

if __name__ == "__main__":
    main()

    

