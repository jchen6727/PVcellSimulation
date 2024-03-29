"""
batch.py 

Batch simulation for M1 model using NetPyNE

Contributors: salvadordura@gmail.com
"""
import numpy as np
from netpyne import specs
from netpyne.batch import Batch

""" Not used here, but may be useful if we modify the model
# ----------------------------------------------------------------------------------------------
# Weight Normalization Exc
# ----------------------------------------------------------------------------------------------

#pops =  ['IT2', 'IT4', 'IT5A', 'IT5B', 'PT5B', 'IT6', 'CT6', 'PV2', 'SOM2'],
def weightNormE(pops = ['PV5B'], secs = None, locs = None,
                allSegs = True, rule = 'PV_reduced', weights =list(np.arange(0.01, 0.2, 0.01)/100.0)):

    # Add params
    from netParams import netParams
    from cfg import cfg

    excludeSegs = ['axon']
    if not secs:
        secs = []
        locs = []
        for secName, sec in netParams.cellParams[rule]['secs'].items():
            if secName not in excludeSegs:
                if allSegs:
                    nseg = sec['geom']['nseg']
                    for iseg in list(range(nseg)):
                        secs.append(secName) 
                        locs.append((iseg+1)*(1.0/(nseg+1)))
                else:
                    secs.append(secName) 
                    locs.append(0.5)

    params = specs.ODict()
    params[('NetStim1', 'pop')] = pops
    params[('NetStim1', 'loc')] = locs
    params[('NetStim1', 'sec')] = secs
    params[('NetStim1', 'weight')] = weights

    groupedParams = [('NetStim1', 'sec'), ('NetStim1', 'loc')]


    initCfg = {}
    initCfg['duration'] = 1.0*1e3
    initCfg[('analysis','plotTraces','timeRange')] = [0, 1000]
    initCfg['weightNorm'] = False
    initCfg['stimSubConn'] = False
    initCfg['addNetStim'] = True
    initCfg[('NetStim1', 'synMech')] = ['AMPA','NMDA']
    initCfg[('NetStim1','synMechWeightFactor')] = [0.5,0.5]
    initCfg[('NetStim1', 'start')] = 700
    initCfg[('NetStim1', 'interval')] = 1000
    initCfg[('NetStim1','ynorm')] = [0.0, 1.0]

    initCfg[('NetStim1', 'pop')] = ['PV5B']

    initCfg[('NetStim1', 'noise')] = 0
    initCfg[('NetStim1', 'number')] = 1
    initCfg[('NetStim1', 'delay')] = 1
    #initCfg[('GroupNetStimW1', 'pop')] = 'None'
    initCfg['addIClamp'] = 0

    b = Batch(params=params, netParamsFile='sim/netParams.py', cfgFile='sim/cfg.py', initCfg=initCfg, groupedParams=groupedParams)

    return b


# ----------------------------------------------------------------------------------------------
# EPSPs via NetStim
# ----------------------------------------------------------------------------------------------
def EPSPs():
    params = specs.ODict()

    params['groupWeight'] = [x*0.05 for x in np.arange(1, 8, 1)]
    params['ihGbar'] = [0.0, 1.0]
 
    
    # initial config
    initCfg = {}
    initCfg['duration'] = 0.5*1e3
    initCfg['addIClamp'] = False
    initCfg['addNetStim'] = True
    initCfg[('GroupNetStimW1', 'pop')] = 'PV5B'
    initCfg[('analysis','plotTraces','timeRange')] = [0, 500]
    initCfg['excTau2Factor'] = 2.0
    initCfg['weightNorm'] = True
    initCfg['stimSubConn'] = False
    initCfg['ihGbarZD'] = None

    groupedParams = [] 

    b = Batch(params=params, netParamsFile='sim/netParams.py', cfgFile='sim/cfg.py', initCfg=initCfg,
              groupedParams=groupedParams)

    return b
"""

# ----------------------------------------------------------------------------------------------
# f-I curve
# ----------------------------------------------------------------------------------------------
def fIcurve(): #TODO: Change values for the IClamp1 so it respects the experimental protocol
    params = specs.ODict()

    params[('IClamp1', 'pop')] = ['PV5B']
    params[('IClamp1', 'amp')] = list(np.arange(0.0, 10.0, 1)/10.0)
    #params['ihGbar'] = [0.0, 1.0, 2.0]
    # params['axonNa'] = [5, 6, 7, 8] 
    # params['gpas'] = [0.6, 0.65, 0.70, 0.75] 
    # params['epas'] = [1.0, 1.05] 
    # params['ihLkcBasal'] = [0.0, 0.01, 0.1, 0.5, 1.0] 

    # initial config
    initCfg = {}
    initCfg['duration'] = 1.0*1e3
    initCfg['addIClamp'] = True
    initCfg['addNetStim'] = False
    initCfg['weightNorm'] = True
    initCfg[('IClamp1','sec')] = 'soma'
    initCfg[('IClamp1','loc')] = 0.5
    initCfg[('IClamp1','start')] = 200
    initCfg[('IClamp1','dur')] = 600
    initCfg[('analysis','plotTraces','timeRange')] = [0, 1000]

    groupedParams = []

    b = Batch(params=params, netParamsFile='sim/netParams.py', cfgFile='sim/cfg.py', initCfg=initCfg, groupedParams=groupedParams)
    b.method = 'grid'

    return b

# ----------------------------------------------------------------------------------------------
# Run configurations
# ----------------------------------------------------------------------------------------------
def setRunCfg(b, type='mpi_bulletin', nodes=1, coresPerNode=8):
    if type=='mpi_bulletin':
        b.runCfg = {'type': 'mpi_bulletin',
            'script': 'sim/init.py',
            'skip': False,
            'skipCustom': '_data.json'}

    elif type=='mpi_direct':
        b.runCfg = {'type': 'mpi_direct',
            'cores': 1,
            'script': 'sim/init.py',
            'mpiCommand': 'mpiexec',
            'skip': False,
            'skipCustom': '_data.pkl'}

# ----------------------------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------

if __name__ == '__main__': 

    b = fIcurve()
    b.batchLabel = 'fIcurve' #'wscale2'
    b.saveFolder = 'data/'+b.batchLabel
    #b.method = 'grid'  # For weight normalization and fI curve
    setRunCfg(b, 'mpi_direct')
    b.run() # run batch 
