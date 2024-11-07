import numpy

ExpAmps = numpy.around(numpy.arange(0.08, 0.54, 0.02), decimals=2) # amplitudes
ExpTargetRates = [1.14, 3., 4.85, 6.57, 7.71, 9.14, 10.43, 11.43, 12.43, 12.86, 13.86, 14.43, 15.15, 15.86, 16.29,
                  17.29, 17.57, 18.57, 19.14, 19.57, 19.43, 19.71, 20.2, 20.]
amps = ExpAmps.tolist()[:8]  # [::3]
amps.insert(0, 0.04)
amps.insert(0, 0.02)
times = list(numpy.arange(1000, 2000 * len(amps), 2000))  # start times
dur = 400  # ms
targetRates = ExpTargetRates[:8]  # [::3]
targetRates.insert(0, 0)
targetRates.insert(0, 0)
targetRates = [i * 1000 / 400 for i in targetRates]
# fitness function
fitnessFuncArgs = {}
fitnessFuncArgs['target'] = {'rates': targetRates}