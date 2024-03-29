# Usage python analysis/fIcurve_calculateAndPlot.py
import json
from netpyne.analysis.tools import loadData
from netpyne.plotting.plotRaster import plotRaster
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys, glob
import time
import numpy.ma as ma
#matplotlib.use("QtAgg") #I need it in Windows to make the plot to appear

def CalculateFI(NetPyNE_dict, DeltaFreqStat=3):
    IStep = NetPyNE_dict['simConfig']['IClamp1']['amp']
    StepDur = NetPyNE_dict['simConfig']['IClamp1']['dur']
    StepStart = NetPyNE_dict['simConfig']['IClamp1']['start']
    SpikeTimes = NetPyNE_dict['simData']['spkt']
    FreqAverage = len(SpikeTimes)/(StepDur/1000)
    FreqInstantaneous = np.ndarray.tolist(1000/np.diff(SpikeTimes)) # To see adaptation of the spike train

    if len(SpikeTimes)>10: # This is important to avoid transitory spikes
        FreqAverageAux = len(SpikeTimes)/((SpikeTimes[-1]-SpikeTimes[0])/1000) # Since there is a big latency for the first spike, we calculate the average frequency here as the number of spikes between the first and the last one
        if abs(FreqInstantaneous[-1]-FreqAverageAux)<DeltaFreqStat:
            FreqStationary = np.mean(FreqInstantaneous[-5:-1])   # Mean of the last N spikes
        else:
            FreqStationary = 0
    else:
        FreqStationary = 0
    return IStep, FreqAverage, FreqInstantaneous, FreqStationary

# ----------------------------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------

if __name__ == '__main__':
    folderSave = 'data/fIcurve/'
    fileName = '*_data.json'

    IStep = []
    FreqAverage = []
    FreqInstantaneous = []
    FreqStationary = []

    for name in glob.glob(folderSave + fileName):
        with open(name) as f:
            fileJSON = json.load(f)
        IStepAux, FreqAverageAux, FreqInstantaneousAux, FreqStationaryAux = CalculateFI(fileJSON)
        IStep.append(IStepAux); FreqAverage.append(FreqAverageAux)
        FreqInstantaneous.append(FreqInstantaneousAux)
        FreqStationary.append(FreqStationaryAux)
    DictResults = {'IStep' : IStep, 'FreqAverage' : FreqAverage, 'FreqInstantaneous' : FreqInstantaneous,
                   'FreqStationary' : FreqStationary}
    with open('analysis/Results.json', 'w', encoding='utf-8') as f:
        json.dump(DictResults, f, ensure_ascii=False, indent=4)
    print(IStep, '\n', FreqAverage, '\n', FreqInstantaneous, '\n', FreqStationary)