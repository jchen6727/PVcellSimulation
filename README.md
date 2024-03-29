# PVcell
## Description
A model of mouse primary motor cortex (M1) parvalbumin-positive (PV) corticospinal cell

Developed using NetPyNE (https://www.netpyne.org/)

## Setup and execution

Requires NEURON with Python and MPI support. 

1. Type or `./compile` (make sure it is an executable) or the equivalent `nrnivmodl mod`. This should create a directory called either i686 or x86_64, depending on your computer's architecture.
2. To run single simulation type: `./runsim [num_proc]` (make sure it is an executable) or the equivalent `mpiexec -np [num_proc] nrniv -python -mpi sim/init.py`
3. To run batch script type: `python sim/batch.py`
4. To run analysis type:  `python analysis\fIcurve_calculateAndPlot.py`

Use the flag `-u` after python if you want real-time printing on console

## Overview of file structure:

* sim/init.py: Main executable; calls functions from other modules. Sets what parameter file to use.

* sim/netParams.py: Network parameters

* sim/cfg.py: Simulation configuration

* sim/batch.py: Batch file

For further information please contact: romanbaravalle@gmail.com

