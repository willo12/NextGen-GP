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

cdef extern void c_nextgen(double* ffs,double* result, double* old_score_vals, double* score_vals,double* obs, long* I, char treefile[], char treeoutfile[], int ffs0, int ffs1, int result0, int my_number,int qsubs,int obs0, int obs1,int runlen,int popsize,double* model,int model0, int ts_factor, int startscore_i, double S_init)

cdef extern void c_single_tree(double* ffs,double* result,double* obs, long* I, char treestr[], int ffs0, int ffs1, int result0,int obs0, int obs1,double* model,int model0, int ts_factor, int startscore_i, double S_init)


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
def nextgen(np.ndarray[double, ndim=2, mode="c"] ffs not None,np.ndarray[double, ndim=2, mode="c"] result not None,np.ndarray[double, ndim=1, mode="c"] old_score_vals not None,np.ndarray[double, ndim=1, mode="c"] score_vals not None,np.ndarray[double, ndim=2, mode="c"] obs not None,np.ndarray[long, ndim=1, mode="c"] I not None,char * treefile,char * treeoutfile, my_number, qsubs, runlen,popsize,np.ndarray[double, ndim=1, mode="c"] model not None, ts_factor=1, startscore_i=0, S_init=-1000):
 
  cdef int ffs0, ffs1, result0, score_vals0, obs0, I0

  ffs0, ffs1 = ffs.shape[0],ffs.shape[1]
  result0 = result.shape[0]  # result should have shape time length x SPACEDIM
  score_vals0=len(score_vals)
  obs0, obs1 = obs.shape[0],obs.shape[1]
  I0 = len(I)
  model0 = len(model)

  start_time = time.time()  
#  profiler_start("yo.log")

  c_nextgen(&ffs[0,0],&result[0,0],&old_score_vals[0],&score_vals[0],&obs[0,0],&I[0], treefile, treeoutfile, ffs0, ffs1,result0,my_number,qsubs,obs0,obs1,runlen,popsize,&model[0],model0, ts_factor, startscore_i, S_init)

#  profiler_stop()
  print "Time: %d"%(time.time() - start_time)

  return None

@cython.boundscheck(False)
@cython.wraparound(False)
def single_tree(np.ndarray[double, ndim=2, mode="c"] ffs not None,np.ndarray[double, ndim=2, mode="c"] result not None,np.ndarray[double, ndim=2, mode="c"] obs not None,np.ndarray[long, ndim=1, mode="c"] I not None,char * treestr,np.ndarray[double, ndim=1, mode="c"] model not None, ts_factor=1, startscore_i=0, S_init=-1000):
 
  cdef int ffs0, ffs1, result0, obs0, I0

  ffs0, ffs1 = ffs.shape[0],ffs.shape[1]
  result0 = result.shape[0]  # result should have shape time length x SPACEDIM

  obs0, obs1 = obs.shape[0],obs.shape[1]
  I0 = len(I)
  model0 = len(model)
  
  c_single_tree(&ffs[0,0],&result[0,0],&obs[0,0],&I[0], treestr, ffs0, ffs1,result0,obs0,obs1,&model[0],model0, ts_factor, startscore_i, S_init)

  return None
