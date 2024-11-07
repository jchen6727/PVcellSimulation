from netpyne.batchtools import specs, comm
from netpyne import sim
import json
import numpy
#from netpyne.specs.dicts import Dict
from pandas import Series, DataFrame
comm.initialize()

"""
see targets.py
"""

AMPS  = [0.02, 0.04, 0.08, 0.1, 0.12  , 0.14  , 0.16  , 0.18 , 0.2   , 0.22  ]
RATES = [0.0 , 0.0 , 2.85, 7.5, 12.125, 16.425, 19.275, 22.85, 26.075, 28.575]

# ---- cfg creation & batch update ---- #

cfg = specs.SimConfig()

cfg.simLabel = 'fi'
cfg.saveFolder = '.'

cfg.x0 = 0.048
cfg.x1 = 0.0001

cfg.recordCells = []
cfg.duration = 800
cfg.dt = 0.050
cfg.recordTraces = {'Vsoma': {'sec': 'soma', 'loc': 0.5, 'var': 'v'}}
cfg.recordStep = 0.1
cfg.analysis['plotTraces'] = {'include': [0,1,2,3,4], 'saveFig': True, 'showFig': False}
cfg.update_cfg()

# ---- FS3 cell loading and customization ---- #
netParams = specs.NetParams()
with open('FS3.json', 'r') as f:
    cell = json.load(f)

#cell['conds']['cellType'] = 'FS3'
cell['secs']['soma']['mechs']['kapin']['gkabar'] = cfg.x0
cell['secs']['soma']['mechs']['kctin']['gkcbar'] = cfg.x1

netParams.cellParams['FS3'] = cell
# ---- create population of cells ---- #
targets = {}
netParams.synMechParams['exc'] = {'mod': 'Exp2Syn', 'tau1': 0.1, 'tau2': 1.0, 'e': 0}
for i, target in zip(AMPS, RATES): # current injection populations
    netParams.popParams['pop_{}'.format(i)] = {'cellType': 'FS3', 'numCells': 1} # create pop
    netParams.stimSourceParams['ic_{}'.format(i)] = {'type': 'IClamp', 'delay': 400, 'dur': 400, 'amp': i} # create current
    netParams.stimTargetParams['ic_{}->FS3'.format(i)] = {'source': 'ic_{}'.format(i), 'conds': {'pop': 'pop_{}'.format(i)}, 'sec': 'soma', 'loc': 0.5} # inject current

    targets['pop_{}'.format(i)] = target

sim.createSimulateAnalyze(netParams=netParams, simConfig=cfg)
print("rates")
rates = Series(dict(sim.analysis.popAvgRates(show=False))) #because this prints data...

targets = Series(targets)

data = DataFrame({'rate': rates, 'target': targets, 'error': abs(rates-targets)})
print(data)

mean_error = numpy.mean(data['error'])

print('Mean error: {}'.format(mean_error))

comm.send(json.dumps({'mean_error': mean_error}))

comm.close()
