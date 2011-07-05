# list of utility functions 

def parLocator(keyWord,b,n,iDoNot,verbose):

#---- -- locate the relevant line in base input file
# --- do not select line iDoNot (unless it is -1)
#
    keyWord=keyWord.lower()
    iFocus=-1
    icol = -1
    for i in range(1, n):

        lineString=str(b[i]).lower()
# check if : exist in line
        try:
            icol=lineString.index(':')
        except ValueError:
             pass   # do nothing
        if icol > -1 :
# verify that this line was not commented out
            try:
                ii=lineString[:icol-1].index('#')
                pass   # do nothing
            except ValueError:
# This line wasn't commented out
                try:
                    ii=lineString.index(keyWord)
                    if i != iDoNot:
# string.index and not string.find is used here, since index raises
# exception when search is failed
                        if verbose: 
                            print 'parLocator: '+str(i)+' found:  '+str(b[i])
                        iFocus=i
                        break
                    else:
                        iFocus=-1
                except ValueError:
                    pass   # do nothing
            
    if iFocus == -1:
        if verbose: 
            print 'parLocator: Keyword ->'+str(keyWord)+'<-  not found'
    return iFocus
    
def stringLocator(keyWord,b,n,verbose):

#---- -- locate the relevant line in a file
#
    keyWord=keyWord.lower()
    iFocus=-1
    for i in range(1, n):
        lineString=str(b[i]).lower()
        try:
            ii=lineString.index(keyWord)
            if verbose: 
                print 'parLocator: '+str(i)+' found:  '+str(b[i])        
            iFocus=i
            break        
        except ValueError:
            pass   # do nothing
            
    if iFocus == -1:
        if verbose: 
            print 'parLocator: Keyword ->'+str(keyWord)+'<-  not found'
    
    return iFocus
    

def readList(dataFile,iLine,verbose):

    from numpy import size
#
#----read list from file to a local float list
#
    listDataLine=dataFile[iLine]
    icol=listDataLine.index(':')
    Data=listDataLine[icol+1:]
    lData=Data.split(',')
    nData=size(lData)     
 
    if verbose:
        print 'readList nData = '+str(nData)
    fData=map(float,lData)
    return fData,nData

def readParameter(dataFile,nLines,keyWord,iDoNot,verbose):

    from numpy import size
#
#----read a parameter from a file-list
#
    keyWord=keyWord.lower()
    ipar = parLocator(keyWord,dataFile,nLines,iDoNot,verbose)
    if ipar == -1:
        if verbose:
            print ' failed to locate '+keyWord+' in base input file; Set value to 1'
        paVal = 1
    else:
        paLine=dataFile[ipar]
        icol=paLine.index(':')
        try:
            iComment=paLine.index('#')
            paVal=paLine[icol+1:iComment-1].lower()
        except ValueError:
            paVal=paLine[icol+1:].lower()

    if verbose:
        if ipar != -1:
            print keyWord+' = '+paVal

    return paVal,ipar

def setContribution(dataFile,nLines,keyWord,iDoNot,verbose):

    from numpy import size
    from myFuncs import parLocator
    import string
#
# default values
#
    nameText=''
    removeContribution = False
#
#----Determine if a given amily contribute to force
#
# Start by locating lines setting contribution
#
    keyWord=keyWord.lower()
    ipar = parLocator(keyWord,dataFile,nLines,iDoNot,verbose)
    if ipar == -1:
        if verbose:
            print ' failed to locate '+keyWord+' in base input file; Set value to 1'
        paVal = 1
    else:
        paLine=dataFile[ipar]
        icol=paLine.index(':')

# Now identify the first part of this line
        firstPart = paLine[0:icol]

# now find out where the standard text ends
        iBF=firstPart.lower().index('family')+6
        nameText=string.join( firstPart[iBF:].split(), "")      

# component name located. Now check about its contribution
        try:
            iComment=paLine.index('#')
            secondPart=paLine[icol+1:iComment-1].lower()
        except ValueError:
            secondPart=paLine[icol+1:].lower()

# find the second colon of this line
        icol2=secondPart.index(':')
        try:
            iComment=secondPart.index('#')
            yesNoText=secondPart[icol2+1:iComment-1].lower()
        except ValueError:
            yesNoText=secondPart[icol2+1:].lower()

        try:
            noFound=yesNoText.lower().index('no')
            removeContribution = True
        except ValueError:
             removeContribution = False

    if verbose:
        if ipar != -1:
            print ' part: '+nameText+' remove contribution: '+str(removeContribution)

    return nameText,removeContribution,ipar


def setPolaraType(ctrl,nc,verbose):

# scan the control file and determine polara type and angles
    from myFuncs import parLocator, readList
# Determine pitch direction from control file
# ---------------------------------------------------

    keyWordPitchAxis='pitch axis'
    iPA = parLocator(keyWordPitchAxis,ctrl,nc,-1,verbose)

    if iPA == -1:
        PA='z' # This is the default
    else:
        paLine=ctrl[iPA]
        icol=paLine.index(':')
        paVal=paLine[icol+1:].lower()
        zFound = 'z' in paVal
        if zFound:
            PA='z'
        else:
            yFound = 'y' in paVal
            if yFound:
                PA='y'
            else:
                raise SystemExit('ERROR in control file: only Y or Z can be given for control keyWord  ->'+keyWordPitchAxis+'<-')
            
    if verbose:
        print 'Pitch axis is '+PA.upper()
#
# angles definitions:
#   alpha ... angle of attack
#   beta  ... side-slip angle
#   phi   ... roll angle
#
# Note: Actually alpha here is the angle of rotation about the above-defined pitch axis
#       Thus, by replacing the pitch-axis, all that is said here about alpha is actually for beta
#
# Several combinations of angles are possible:
#------------------------------------------------
# 1. Polar-sweep in alpha per given phi                        ...... polarVar   = aoa
# 2. Polar-sweep in alpha per given beta (side slip angle)     ...... polarVar   = aoa
# 3. Polar-sweep in phi per given alpha                        ...... polarVar   = phi
#
# Note: Seting both phi and beta is impossible
#
# Now let us find out which angles are specified in the control file, to figure out polarSweepType and polarVar
#
    keyWordListAoA='angles of attack'
    iListAOA = parLocator(keyWordListAoA,ctrl,nc,-1,verbose)
    keyWordListPhi='roll angles'
    iListPhi = parLocator(keyWordListPhi,ctrl,nc,-1,verbose)
    keyWordListBeta='side slip angle'
    iListBeta = parLocator(keyWordListBeta,ctrl,nc,-1,verbose)

    if iListPhi == -1 :
        if  iListBeta == -1 :
            polarSweepType=1 ; 
            polarVar='aoa' ;  # phi/beta not found. Polar sweep in alpha for phi=beta=0
            phi = [0.0] ; nPhi = 1; 
            nBeta=0 ; beta=[ ];
        else:
#         beta was found in control file, phi is not there; check how about alpha
            nPhi = 0 ; phi=[ ];
            polarSweepType=2 ; 
            polarVar='aoa' ; 
            beta,nBeta=readList(ctrl,iListBeta,verbose)
            if nBeta > 1 :
                raise SystemExit('ERROR in control file: nBeta > 1. For polar sweep in beta exchange pitch-axis and use aoa')
        
        if  iListAOA == -1 :
            raise SystemExit('ERROR in control file: phi and alpha are missing. Polar sweep not defined')

        alpha,nAalpha=readList(ctrl,iListAOA,verbose)
        
    else:   
#      phi was found in control file, so beta must not be there
        if iListBeta > -1:
            raise SystemExit('ERROR in control file: both phi and beta specified.  Polar sweep not defined ')

        nBeta=0 ; beta=[ ];
#     Check now if alpha appears
        if  iListAOA == -1 :
#     phi found in control file, but alpha is missing, so it is a polar-sweep in phi with alpha=0
            polarSweepType=3 ; 
            polarVar='phi' ; alpha =[0.0]  ; nAalpha=1  

        else:
#
#     Both alpha and phi found in control file. Find out which one is a list
#
            alpha,nAalpha=readList(ctrl,iListAOA,verbose)

        phi,nPhi=readList(ctrl,iListPhi,verbose)

        if nAalpha == 1:
            if nPhi > 1 :
                polarSweepType=3  
                polarVar='phi' 
            else:
                polarSweepType=1 
                polarVar='aoa'

            nBeta=0; beta=[ ];
        else:
#-----that is nAlpha > 1
            if nPhi > 1 :
                 raise SystemExit('ERROR in control file: read lists in both alpha and phi.  Polar sweep not defined ')

            polarSweepType=1 ; 
            polarVar='aoa'; nBeta=0; beta=[ ];

#-------------------------------------------------------------------------------------------
    if verbose:
        if polarSweepType == 1 :
            print 'Sweep type: '+str(polarSweepType)+' in alpha. nAalpha = '+str(nAalpha)+' phi = '+str(phi)
        elif polarSweepType == 2 :
            print 'Sweep type: '+str(polarSweepType)+' in alpha. nAalpha = '+str(nAalpha)+' beta = '+str(beta)
        else:
            print 'Sweep type: '+str(polarSweepType)+' in phi. nPhi = '+str(nPhi)+' alpha = ',str(alpha)

    nPolara = max(nAalpha,nPhi)
    return PA,polarSweepType,nAalpha,nBeta,nPhi,alpha,beta,phi,polarVar
   
def setVelDir(polarSweepType,PA,alphar,phir,betar):

# set the velocity direction
    from numpy import sin,cos,tan
    from myFuncs import parLocator, readList

    if polarSweepType == 2 :

        if PA == 'z':
            dv1 = cos(betar)
            dv2 = tan(alphar)*cos(betar)  
            dv3 = [sin(betar) for x in alphar]
        else:
            dv1 = cos(betar)
            dv2 = [sin(betar)  for x in alphar]
            dv3 = tan(alphar)*cos(betar)
    else:
        dv1=1.0
        if PA == 'z':
            dv2=tan(alphar)*cos(phir)  
            dv3=tan(alphar)*sin(phir) 
        else:
            dv2=tan(alphar)*sin(phir)
            dv3=tan(alphar)*cos(phir) 
    
    return dv1,dv2,dv3

def processAddAngle(addRunStr,nPolara,parAngle):
#
#---------------------------------------------------------------------
# Process the list of interactively added angles
#--------------------------------------------------------------------

    from numpy import size, sort 

#------- criterion for identical angles
    angleEqualCriterion=0.1

#------------ create a list out of entered angles

    addRunList=addRunStr.split(',')
    fAddRunListRaw=map(float,addRunList)
    fAddRunList=sort(fAddRunListRaw)
    nAddRun=size(fAddRunList)

#------- By default, do not compute the cases in input file, since they were computed already

    computeCase=[False for j in range(0,nPolara+nAddRun)]
    rerunCase = [False for j in range(0,nPolara+nAddRun)]


#----- Now, for each  new angle verify if it is a rerun or inserted new value

    closestAngle=[0 for i in range(0,nAddRun)]
    for i in range(0,nAddRun):

        diff=[abs(fAddRunList[i]-x) for x in parAngle]
        iClose=diff.index(min(diff))
        closestAngle[i]=parAngle[iClose]
        if min(diff) < angleEqualCriterion:
            computeCase[iClose]=True
            rerunCase[iClose]=True
            
        else:
                
            if fAddRunList[i] > closestAngle[i]:
                ii=iClose+1
            else:
                ii=iClose
                    
            tmpAng1=parAngle[:ii]
            tmpAng1.append(fAddRunList[i])
            tmpAng2=parAngle[ii:]
            tmpAng1.extend(tmpAng2)

            parAngle=tmpAng1
            nPolara=nPolara+1
            computeCase[ii]=True


    return  nPolara,parAngle,computeCase,rerunCase   
     
def updatedControlFile(ctrl,nc,parAngle,ctrlFile,verbose):

# generate a modified control file for case with addRun options

    from myFuncs import parLocator, readList
    import os

#
#-- get a proper list of updated parameter-angle
    st1=str(parAngle)
    updaedAngleList=st1[1:-1]
# Now let us find out which angles are specified in the control file, to figure out polarSweepType and polarVar
#
    keyWordListAoA='angles of attack'
    iListAOA = parLocator(keyWordListAoA,ctrl,nc,-1,verbose)
    keyWordListPhi='roll angles'
    iListPhi = parLocator(keyWordListPhi,ctrl,nc,-1,verbose)
    keyWordListBeta='side slip angle'
    iListBeta = parLocator(keyWordListBeta,ctrl,nc,-1,verbose)

    if iListPhi == -1 :
        if  iListBeta == -1 :
            polarSweepType=1 ; 
            polarVar='aoa' ;  # phi/beta not found. Polar sweep in alpha for phi=beta=0
            phi = [0.0] ; nPhi = 1; 
            nBeta=0 ; beta=[ ];
        else:
#         beta was found in control file, phi is not there; check how about alpha
            nPhi = 0 ; phi=[ ];
            polarSweepType=2 ; 
            polarVar='aoa' ; 
            beta,nBeta=readList(ctrl,iListBeta,verbose)
            if nBeta > 1 :
                raise SystemExit('ERROR in control file: nBeta > 1. For polar sweep in beta exchange pitch-axis and use aoa')
        
        if  iListAOA == -1 :
            raise SystemExit('ERROR in control file: phi and alpha are missing. Polar sweep not defined')

        alpha,nAalpha=readList(ctrl,iListAOA,verbose)
        ctrl[iListAOA]=' angles of attack :  '+updaedAngleList+'\n'
        
    else:   
#      phi was found in control file, so beta must not be there
        if iListBeta > -1:
            raise SystemExit('ERROR in control file: both phi and beta specified.  Polar sweep not defined ')

        nBeta=0 ; beta=[ ];
#     Check now if alpha appears
        if  iListAOA == -1 :
#     phi found in control file, but alpha is missing, so it is a polar-sweep in phi with alpha=0
            polarSweepType=3 ; 
            polarVar='phi' ; alpha =[0.0]  ; nAalpha=1 

        else:
#
#     Both alpha and phi found in control file. Find out which one is a list
#
            alpha,nAalpha=readList(ctrl,iListAOA,verbose)
            if nAalpha > 1 :
                ctrl[iListAOA]=' angles of attack :  '+updaedAngleList+'\n'

        phi,nPhi=readList(ctrl,iListPhi,verbose)
        if nPhi > 1:
            ctrl[iListPhi]=' roll angles :  '+updaedAngleList+'\n'
       
# Prepare a backup of control file

    st1='cp '+ctrlFile+' '+ctrlFile+'.bck'
    os.system(st1)
#
# --- Write down the updated file
    fc=open(ctrlFile,'w')
    fc.writelines(ctrl)
    fc.close()

    print 'More cases were added. Original ctrl file saved at '+ctrlFile+'.bck File '+ctrlFile+' updated'


    return

   



