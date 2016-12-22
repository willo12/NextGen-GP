import os.path
from math import sqrt
import numpy as np
import os
import h5py as hdf
import sys
import time

from copy import deepcopy
import subprocess
from configobj import ConfigObj, ConfigObjError, flatten_errors
from validate import Validator
from my_data import glacial

import nextgen



def load_config():
  # get directory component of pathname using dirname (would get the file with basename)
  dir_path = os.path.dirname(os.path.realpath(__file__)) # get the canonical path to current file (shortest absolute).
  config_name = 'config.ini'
  try:
    config = ConfigObj(config_name, configspec=os.path.join(dir_path,'configspec.ini'), file_error=True)
  except (ConfigObjError, IOError), e:
    raise Exception('Could not read "%s": %s' % (config_name, e))

  validator = Validator()
  results = config.validate(validator)
  
  if results != True:
    for (section_list, key, _) in flatten_errors(config, results):
      if key is not None:
        raise Exception('The "%s" key in the section "%s" failed validation' % (key, ', '.join(section_list)) )
      else:
        raise Exception('The following section was missing:%s ' % ', '.join(section_list) )

  return config

def getf(fname):
  f = open(fname,"r")
  content = f.read()
  f.close()
  return content

def putf(fname,content):
  f = open(fname,"w")
  f.write(content)
  f.close()





def launch(com='python test.py',pbs_name = None,num_proc=1,vmem="1gb",walltime=6,workdir="$PBS_O_WORKDIR",
  prologue="",epilogue=""):

  pbs_content = """#!/bin/bash
 
#PBS -N %s
#PBS -l nodes=1:ppn=%d
#PBS -l vmem=%s
#PBS -l walltime=%d:00:00
#PBS -j oe

%s
%s

cd %s

module load python/2.7.9
 
%s

"""

#  pbs_name = None
 
  print "launching %s"%com
  if pbs_name == 'None':
    if not os.path.isfile("proc.lck"): # external safeguard
      print 'pbs_name None, running job in shell, no pbs.'
      subprocess.Popen(com, shell=True)
  else:
#    if not os.path.isfile(pbs_name): # create pbs file
    putf(pbs_name,pbs_content%(pbs_name,num_proc,vmem,walltime,prologue,epilogue,workdir,com))
    if not os.path.isfile("sub.lck"): # external safeguard
      os.system("qsub %s"%pbs_name)
    else:
      print "Wrote %s, no submit due to sublock"%pbs_name

def landscape(max_x,max_y,sigma = 2):
  """
  must have max and min at equator
  """

  y1 = 1.2*max_y/4    
  x1 = 1.2*max_x/4

  y2 = 2.8*max_y/4    
  x2 = 2.8*max_x/4

  f = np.zeros((max_y,max_x))

  for y in xrange(0,max_y):
    for x in xrange(0,max_x):
      f[y,x] += np.exp( -(0.5/sigma**2)*( (x-x1)**2 + (y-y1)**2  ) )
      f[y,x] -= np.exp( -(0.5/sigma**2)*( (x-x2)**2 + (y-y2)**2  ) )

  return f

def landscape2(max_x,max_y,sigma = 2):
  """
  landscape of amplitudes
  """

  y1 = max_y/2    
  x1 = max_x/2

  f = np.zeros((max_y,max_x))

  for y in xrange(0,max_y):
    for x in xrange(0,max_x):
      f[y,x] += np.exp( -(0.5/sigma**2)*( (x-x1)**2 + (y-y1)**2  ) )

  return f

def landscape_val(x,y,max_x,max_y,sigma = 2):
  """
  must have max and min at equator
  """

  y1 = 1.2*max_y/4    
  x1 = 1.2*max_x/4

  y2 = 2.8*max_y/4    
  x2 = 2.8*max_x/4

  return np.exp( -(0.5/sigma**2)*( (x-x1)**2 + (y-y1)**2  ) ) - np.exp( -(0.5/sigma**2)*( (x-x2)**2 + (y-y2)**2  ) )


def landscape_val2(x,y,max_x,max_y,sigma = 2):
  """
  must have max and min at equator
  """

  y1 = max_y/2    
  x1 = max_x/2

  return np.exp( -(0.5/sigma**2)*( (x-x1)**2 + (y-y1)**2  ) ) 

def prep_landscape():

  pass

#my_number_y = (my_number/qcols)%qrows
#my_number_x = my_number%qcols

#qrows = int(sqrt(qsubs));
#qcols = qsubs/qrows;

#Ampl = 1
#sigma = 1
#t_A_cutoff = 75

#if (timestamp*runlen < t_A_cutoff):
#  A=Ampl*(1-  (1-float(timestamp*runlen)/t_A_cutoff)*landscape_val2(my_number_x,my_number_y,qcols,qrows,sigma) )
#else:
#  A=Ampl

#print 'ampl: %f timestamp: %f\n'%(Ampl*(1.0-float(timestamp)/t_A_cutoff) ,timestamp)
#print Ampl*(1.0-float(timestamp)/float(t_A_cutoff))

#lenresult = (0-t0)/dt+1
#mid=-(t0+1.e6)/dt
#result=np.zeros((0-t0)/dt+1 , dtype=np.float64)  # no initial field given

#result = 1 + (2*A*np.arctan(  (4*dt*1e-5)*(np.arange(lenresult) - mid)  )/np.pi).astype(np.float64) 

#if (abs(A) > 0.01) and (timestamp*runlen < t_A_cutoff):
#  result = 1 + (2*A*np.arctan(  (4*dt*1e-5)*(np.arange(lenresult) - mid)  )/np.pi).astype(np.float64) 
#else:
#  result=np.ones((0-t0)/dt+1 , dtype=np.float64)

#model = np.polyfit(otime,d18o, 2).astype(np.float64)

#result = np.polyval(model, np.arange(t0,dt,dt) ).astype(np.float64)


def smoothGaussian(array,strippedXs=False,degree=5):  

  if degree==0:
    return array

  window=degree*2-1  

  weight=np.array([1.0]*window)  

  weightGauss=[]  

  for i in range(window):  

    i=i-degree+1  

    frac=i/float(window)  

    gauss=1/(np.exp((4*(frac))**2))  

    weightGauss.append(gauss)  

  weight=np.array(weightGauss)*weight  

  smoothed=[0.0]*(len(array)-window)  

  for i in range(len(smoothed)):  

    smoothed[i]=sum(np.array(array[i:i+window])*weight)/sum(weight)  

  return np.array(smoothed)  

def prep_run(config, SPACEDIM=2):
  """
  Prepare to run tree from c

  
  """

  S_init = config['driver']['S_init']
  compgridsize = config['driver']['compgridsize']
  mutationrate = config['driver']['mutationrate'] 
  tour = config['driver']['tour'] 

  Iffs, obs,I, ts_factor, startscore_i = glacial(**(config['glacial'].dict()))

  # raise Exception('TEST: %s'%str(startscore_i))

  otime = obs[startscore_i:,0]

  obs = obs[:,1:].copy(order='C'); # removing time
  I = I.copy(order='C')

  aspect = -1

  result=np.ones((Iffs.shape[0] , SPACEDIM ), dtype=np.float64)

  model = np.zeros(3).astype(np.float64)
  model[0]=aspect

  return Iffs, obs,I, ts_factor, startscore_i, S_init, obs, I, aspect, result, model, otime, compgridsize, mutationrate, tour

def batch_run(old_score_vals,score_vals, treefile,treeoutfile,my_number,qsubs,runlen,popsize,config=None,SPACEDIM=2):
  """ Run the main GP loop.

      Prepares data, reads config and calls netgen in c module. 

  """
  if config is None:
    config = load_config()

  Iffs, obs,I, ts_factor, startscore_i, S_init, obs, I, aspect, result, model, otime, compgridsize, mutationrate, tour = prep_run(config,SPACEDIM)

  # Start c code. This will read treefile and save tree file with filename treeoutfile
  nextgen.nextgen(Iffs,result,old_score_vals,score_vals,obs,I, treefile,treeoutfile,my_number,qsubs,runlen,popsize,model,ts_factor,startscore_i, S_init, compgridsize, mutationrate, tour)

def stability_test(tree,config=None,SPACEDIM=2,ts_factor=[2,8]):
  """
  Test whether tree integration is stable with respect to time step

  config: run configuration to use
  ts_factor: list of time ts_factors (in c code) to test. time step is forcing time step divided by ts_factor: higher ts_factor is finer integration time step.


  """

  if isinstance(tree,str):
    tree = [tree,]

  ts_factor_orig = ts_factor

  #tree = tree[:10]
  
  tree_list = []

  for treestr in tree:
    for i in range(len(ts_factor)):
      tree_list.append(treestr)    

  Iffs,result,obs,I, tree_list,model,ts_factor,startscore_i, S_init, otime, scores = single_run(tree_list,ts_factor=[2,8])

  n_perturbations = len(ts_factor_orig) - 1
  stab_scores = []
  for j,treestr in enumerate(tree):
    for i in range( n_perturbations ):
      stab_scores.append( scores[j*n_perturbations+i+1]-scores[j*n_perturbations] )      
     

  return stab_scores

def single_run(tree,config=None,SPACEDIM=2,ts_factor=None):

  if isinstance(tree,str):
    tree = [tree,]

  if config is None:
    config = load_config()

  Iffs, obs,I, ts_factor_config, startscore_i, S_init, obs, I, aspect, result, model, otime = prep_run(config,SPACEDIM)

  if ts_factor is None:
    ts_factor = ts_factor_config

  # Start c code. This will read treefile and save tree file with filename treeoutfile
  if isinstance(ts_factor,int):
    ts_factor = [ts_factor,]

  if len(ts_factor) < len(tree):
    if (len(tree)%len(ts_factor) == 0):
      ts_factor = ts_factor*int(len(tree)/len(ts_factor)) # recurring pattern
    else:
      ts_factor = ts_factor[0]*int(len(tree)) # constant value assumed

  scores = []
  results = []

  for i, treestr in enumerate(tree):    
    nextgen.single_tree(Iffs,result,obs,I, treestr,model,ts_factor[i],startscore_i, S_init)
    results.append(deepcopy(result))
    scores.append(model[1])

  if (len(tree) == 1):
    results = results[0]

  return Iffs,results,obs,I, tree,model,ts_factor,startscore_i, S_init, otime, scores
  


