"""
doloop.pyx

"""

import cython
import time

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np


#cdef extern from "gperftools/profiler.h":
#    int ProfilerStart( char* fname )
#    int ProfilerStop()
 
#def profiler_start(fname):
#    ProfilerStart(<char *>fname)
 
#def profiler_stop():
#    ProfilerStop()

cdef extern void c_map_f(char treestring[], double* A, double* B, int shape_i, int shape_j, double S0_start, double S1_start,double S0_end, double S1_end, double * forcing, int forcing_dim)

cdef extern void c_nextgen(int my_number,int qsubs,int runlen,int popsize, int compgridsize, double mutationrate, int tour)

cdef extern void c_single_tree(char treestr[], int ts_factor, double * S_init_array)


@cython.boundscheck(False)
@cython.wraparound(False)
def map_f(char * treestring, np.ndarray[double, ndim=2, mode="c"] A not None,np.ndarray[double, ndim=2, mode="c"] B not None ,np.ndarray[double, ndim=1, mode="c"] forcing not None,forcing_dim=1, S0_start=-2., S1_start=-2., S0_end=2., S1_end=2.):
  """
    Calculate d/dt of S_0 and S_1 
   
    A, B: 2 dim np array containing results for S_0 and S_1
 
    forcing: 1 dim np array. single value for each forcing. e.g. np.array([1,0.2]) for f_1 = 1, f_2 = 0.2

  """
 
  cdef int shape_i,shape_j

  shape_j, shape_i = A.shape[0], A.shape[1]

  c_map_f(treestring, &A[0,0], &B[0,0], shape_i, shape_j, S0_start, S1_start, S0_end, S1_end, &forcing[0],forcing_dim)

  return None


@cython.boundscheck(False)
@cython.wraparound(False)
def nextgen( my_number, qsubs, runlen,popsize, compgridsize = 10, mutationrate = 0.08, tour = 15, S_init = -100):
  
  start_time = time.time()  
#  profiler_start("yo.log")

  c_nextgen(my_number,qsubs,runlen,popsize, compgridsize, mutationrate, tour)

#  profiler_stop()
  print "Time: %d"%(time.time() - start_time)

  return None

@cython.boundscheck(False)
@cython.wraparound(False)
def single_tree(char * treestr, np.ndarray[double, ndim=1, mode="c"] S_init_array not None, ts_factor=1):
   
  c_single_tree(treestr, ts_factor, &S_init_array[0])

  return None
