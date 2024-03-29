from netpyne.analysis.tools import loadData
from netpyne.plotting.plotRaster import plotRaster
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import numpy.ma as ma
matplotlib.use("QtAgg")

filename = "PT5B_full_cellParams.pkl"

filePKL = loadData(filename)

print([i for i in filePKL['conds']])
print([i for i in filePKL['secs']])
print([i for i in filePKL['secLists']])
print([i for i in filePKL['globals']])
