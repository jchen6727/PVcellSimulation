from batchtk.runtk.trial import trial, LABEL_POINTER, PATH_POINTER
from netpyne.batchtools import dispatchers, submits
from cmaes import CMA
import numpy
import os
Dispatcher = dispatchers.INETDispatcher
Submit = submits.SHSubmitSOCK

cwd = os.getcwd()

# evaluation
def eval_rosenbrock(x0, x1, tid):
    cfg = {
        'x0': x0,
        'x1': x1,
    }
    submit = Submit()
    submit.update_templates(**{'command': 'python fi_curve.py',})
    label = 'fi_curve'
    return float(trial(cfg, label, tid, Dispatcher, cwd, '../cma', submit)['mean_error'])




# suggestor
optimizer = CMA(mean=numpy.zeros(2), bounds=numpy.array([[0,1], [0,1]]), sigma=1.0) # create an array of zeros with the number of parameters
for generation in range(3):
    solutions = []
    for cand in range(optimizer.population_size):
        x = optimizer.ask()
        value = eval_rosenbrock(x[0], x[1], "{}_{}".format(cand, generation))
        solutions.append((x, value))
        print(f"#{generation} {value} (x0={x[0]}, x1={x[1]})")
    optimizer.tell(solutions)
