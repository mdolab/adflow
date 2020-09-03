"""these are a set of default options for the regression tests"""



defaultFuncList = ['lift', 'drag', 'cl', 'cd', 'fx', 'fy', 'fz', 'cfx', 'cfy', 'cfz',
                   'mx', 'my', 'mz', 'cmx', 'cmy', 'cmz', 'sepsensor', 'sepsensoravgx',
                   'sepsensoravgy', 'sepsensoravgz']

defaultAeroDVs = ['alpha', 'beta', 'mach', 'P', 'T', 'xRef', 'yRef', 'zRef']

# Note that the option keys here are all consistently LOWERCASE.
adflowDefOpts = {
    # Common Paramters
    'gridfile':'default.cgns',
    'restartfile':None,
    'outputdirectory':'../output_files',

    # Output Parameters
    'storerindlayer':True,
    'writesurfacesolution':True,
    'writevolumesolution':True,
    'nsavevolume':1,
    'nsavesurface':1,
    'solutionprecision':'double',
    'gridprecision':'double',
    'isosurface':{},
    'isovariables':[],
    'viscoussurfacevelocities':True,

    # Physics Paramters
    'discretization':'central plus scalar dissipation',
    'coarsediscretization':'central plus scalar dissipation',
    'limiter':'vanalbeda',
    'smoother':'runge kutta',
    'equationtype':'euler',
    'equationmode':'steady',
    'flowtype':'external',
    'turbulencemodel':'sa',
    'turbulenceorder':'first order',
    'turbresscale':10000.0,
    'usewallfunctions':False,
    'useapproxwalldistance':True,
    'eulerwalltreatment':'linear pressure extrapolation',
    'dissipationscalingexponent':0.67,
    'vis4':0.0156,
    'vis2':0.25,
    'vis2coarse':0.5,
    'restrictionrelaxation':.80,
    'liftindex':2,
    'lowspeedpreconditioner':False,

    # Common Paramters
    'ncycles':500,
    'ncyclescoarse':500,
    'nsubiterturb':1,
    'nsubiter':1,
    'cfl':1.7,
    'cflcoarse':1.0,
    'mgcycle':'3w',
    'mgstartlevel':-1,
    'resaveraging':'alternateresaveraging',
    'smoothparameter':1.5,
    'cfllimit':1.5,

    # Unsteady Paramters
    'timeintegrationscheme':'bdf',
    'timeaccuracy':2,
    'ntimestepscoarse':48,
    'ntimestepsfine':400,
    'deltat':.010,

    # Grid motion parameters
    'useale':True,
    'usegridmotion':False,

    # Time Spectral Paramters
    'timeintervals':1,
    'alphamode':False,
    'betamode':False,
    'machmode':False,
    'pmode':False,
    'qmode':False,
    'rmode':False,
    'altitudemode':False,
    'windaxis':False,
    'alphafollowing':True,
    'tsstability':False,

    # Convergence Paramters
    'l2convergence':1e-6,
    'l2convergencerel':1e-16,
    'l2convergencecoarse':1e-2,
    'maxl2deviationfactor':1.0,

    # Newton-Krylov Paramters
    'usenksolver':False,
    'nkswitchtol':2.5e-4,
    'nksubspacesize':60,
    'nklinearsolvetol':0.3,
    'nkuseew':True,
    'nkadpc':False,
    'nkviscpc':False,
    'nkasmoverlap':1,
    'nkpcilufill':2,
    'nkjacobianlag':20,
    'rkreset':False,
    'nrkreset':5,
    'applypcsubspacesize':10,
    'nkinnerpreconits':1,
    'nkouterpreconits':1,
    'nkls':'cubic',
    'nkcfl0':1000000000000.0,

    # Approximate Newton-Krylov Parameters
    'useanksolver':False,
    'ankswitchtol':1e-2,
    'anksubspacesize':5,
    'anklinearsolvetol':0.5,
    'ankasmoverlap':1,
    'ankpcilufill':1,
    'ankjacobianlag':20,
    'ankinnerpreconits':1,

    # LoadBalance/partitioning Parameters
    'blocksplitting':False,
    'loadimbalance':0.1,
    'loadbalanceiter':10,
    'partitiononly':False,

    # Misc Paramters
    'autosolveretry':False,
    'autoadjointretry':False,
    'numbersolutions':True,
    'printiterations':True,
    'printtiming':True,
    'setmonitor':True,
    'printwarnings':True,
    'monitorvariables':['cpu', 'resrho', 'resturb', 'cl', 'cd', 'cmz','totalr'],
    'surfacevariables':['cp', 'vx', 'vy', 'vz', 'mach'],
    'volumevariables':[],

    # Multidisciplinary Coupling Parameters:
    'forcesastractions':True,

    # Adjoint Paramters
    'adjointl2convergence':1e-6,
    'adjointl2convergencerel':1e-16,
    'adjointl2convergenceabs':1e-16,
    'adjointdivtol':1e5,
    'approxpc':True,
    'adpc':False,
    'viscpc':False,
    'usediagtspc':True,
    'restartadjoint':True,
    'adjointsolver':'gmres',
    'adjointmaxiter':500,
    'adjointsubspacesize':100,
    'adjointmonitorstep':10,
    'dissipationlumpingparameter':6.0,
    'preconditionerside':'right',
    'matrixordering':'rcm',
    'globalpreconditioner':  'additive schwartz',
    'localpreconditioner':'ilu',
    'ilufill':2,
    'asmoverlap':1,
    'innerpreconits':1,
    'outerpreconits':3,
    'applyadjointpcsubspacesize':20,
    'frozenturbulence':True,
    'usematrixfreedrdw':True,

    # ADjoint debugger
    'firstrun':True,
    'verifystate':True,
    'verifyspatial':True,
    'verifyextra':True,

    # Function Parmeters
    'sepsensoroffset':0.0,
    'sepsensorsharpness':10.0,
}


pyWarpDefOpts = {}

IDWarpDefOpts = {
    'fileType':'cgns',
    'specifiedSurfaces':None,
    'symmetrySurfaces':None,
    'symmetryPlanes':None,
    'aExp': 3.0,
    'bExp': 5.0,
    'LdefFact':1.0,
    'alpha':0.25,
    'errTol':0.0005,
    'evalMode':'fast',
    'symmTol':1e-6,
    'useRotations':True,
    'zeroCornerRotations':True,
    'cornerAngle':30.0,
    'restartFile':None,
    'bucketSize':8,
}
