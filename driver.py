#!/usr/bin/env python

import os.path
from math import sqrt
import numpy as np

import os
from runtools import launch, getf, putf, load_config, batch_run

import h5py as hdf
import sys
import time

SPACEDIM = 2  # must correspond with c_nextgen.c value for SPACEDIM

timestamp_flag = False

if os.path.isfile("disp.lck"):
  print "Driver locked by disp.lck file, exiting..."
  sys.exit();

my_number = int(sys.argv[1])
popsize = int(sys.argv[2])
qsubs = int(sys.argv[3])
runlen = int(sys.argv[4])
qsub_name = sys.argv[5]

if timestamp_flag:
  timestampfile = 'ts%d'%my_number

  try:
    timestamp = int(getf(timestampfile))
  except:
    timestamp = 0

config = load_config()



tree_names = "popout%d"
score_names = "scores%d"

treefile ="pop%d"%my_number
treeoutfile =treefile

score_file=score_names%my_number

# init memory
old_score_vals=1e5*np.ones(popsize+1, dtype=np.float64)
score_vals=1e5*np.ones(popsize+1, dtype=np.float64)

if os.path.isfile(score_file):
  infile = hdf.File(score_file, "r")
  old_score_vals = infile["/scv"][:]        
  infile.close()
  print "Read scores of length %d"%len(old_score_vals)


#raise Exception("obs[0,0]: %g"%obs[0,0])
     
# Start c code. This will read treefile and save tree file with filename treeoutfile
#nextgen.nextgen(Iffs,result,old_score_vals,score_vals,obs,I, treefile,treeoutfile,my_number,qsubs,runlen,popsize,model,ts_factor,startscore_i, S_init)

batch_run(old_score_vals,score_vals, treefile,treeoutfile,my_number,qsubs,runlen,popsize,config,SPACEDIM)

# save score ndarray
outfile = hdf.File(score_file,"w")
outfile.create_dataset("scv",data=score_vals)
outfile.close()

if timestamp_flag:
  print timestamp
  putf(timestampfile,'%d'%(timestamp+1))

#print min(score_vals)

# bash command expects a symlink driver to driver.py, accessible via the PATH (e.g. ~/bin) 
command = "driver %d %d %d %d %s > pr%d"%(my_number,popsize, qsubs, runlen,qsub_name,my_number)

if os.path.isfile("disp.lck"):
  print "Driver locked by disp.lck file, no launch."
else:
  try:
    launch(com=command, pbs_name=qsub_name,vmem="2gb",walltime=10)
  except:
    time.sleep(30)
    launch(com=command, pbs_name=qsub_name,vmem="2gb",walltime=10)



