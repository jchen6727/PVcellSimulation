from netpyne import specs
import json

netparams = specs.NetParams()

cell = netparams.importCellParams(label='FS3', cellName='FScell', fileName='FS3.hoc') # template name is FScell

#with open('FS3_err.json', 'w') as f:
#    json.dump(cell, f, indent=4)
