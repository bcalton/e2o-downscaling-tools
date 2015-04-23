
# coding: utf-8

# In[1]:

import sys
sys.path.append("../../e2o_dstools/")
import e2o_calculateEvaporation


# In[2]:

# import message passing interface for python
from mpi4py import MPI


# In[1]:

# import message passing interface for python
from mpi4py import MPI

# import for memory use
from pympler import tracker
tr = tracker.SummaryTracker()
tr.print_diff() 


comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size

if rank == 0:
    # todo = [(row, col) for row in range(nrow) for col in range(ncol)]
    # Reorganize work a bit so we can scatter it
    keyfunc = lambda x:x[0] % size
    work = itertools.groupby(sorted(enumerate(list_pars), key=keyfunc), keyfunc)
    # Data is now in the format:
    # Expand the work so we get lists of row, col per node
    workpernode = [[x[1] for x in val] for (key, val) in work]
else:
    workpernode = None


# Now distribute the workload accross all processes
workpernode = comm.scatter(workpernode, root=0)
# workpernode = workpernode[0]
#logging.info("Got the following work in process rank {} : {}".format(rank, workpernode))

# Each node can now do it's own work. The main advantage is that we can do a gather at the end to collect all results.
# Keep track of all the runs per node in scores
scores = []

# before starting any runs, make sure that you know in which folder we run this MPI run routine. 
# Always return to this folder before the next run
curdir = os.getcwd()
for i, par_set in workpernode:
    logging.info("rank %02.f computing scores for parameter set nr %04.f" % (rank, i))

    runId = '%s_%04.f' % (R, i)
    keys=[]
    values=[]
    for key, value in par_set.items():
        keys.append ( key.replace('\r',''))
        values.append(value)
        par_set2 = dict(zip(keys, values))
    parmult = '%s' %(par_set2)
    argv = ['-C', caseFolder, '-R',runId, '-c', inifile, '-I', '-T', str(timeSteps), '-s 3600', '-P', parmult, '-l ERROR']
#    argv = ['-C', caseFolder, '-R',runId, '-c', inifileNew, '-T', str(timeSteps), '-s 86400']    
    # run model, perhaps return something to store in a list
    wf.main(argv)
    tr.print_diff()
    scores.append(i)
    # Now make sure that you return to the right directory, in case the model has changed the path
    os.chdir(curdir)

# Wait here so we can collect all runs
# Because we distributed the work evenly all processes should be here at approximately the same time
comm.Barrier()
# Great, we're all here. Now let's gather the scores...
# Collect values from all the processes in the main root
scores = comm.gather(scores, root=0)
# Only the root node should collect all the data
logging.debug("Rank {} has scores {}".format(rank, scores))




# In[ ]:



