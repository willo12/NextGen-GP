#!/usr/bin/env python

import os.path
from math import sqrt
import numpy as np

import os
from runtools import launch, getf, putf, load_config, batch_run

import sys
import time

timestamp_flag = False

if os.path.isfile("disp.lck"):
  print "Driver locked by disp.lck file, exiting..."
  sys.exit();

my_number = int(sys.argv[1])

qsub_name = sys.argv[2]

if timestamp_flag:
  timestampfile = 'ts%d'%my_number

  try:
    timestamp = int(getf(timestampfile))
  except:
    timestamp = 0

config = load_config()

#raise Exception("obs[0,0]: %g"%obs[0,0])
     
# Start c code. 

batch_run(my_number,config)


if timestamp_flag:
  print timestamp
  putf(timestampfile,'%d'%(timestamp+1))

# bash command expects a symlink driver to driver.py, accessible via the PATH (e.g. ~/bin) 
command = "driver %d  %s > pr%d"%(my_number,qsub_name,my_number)

if os.path.isfile("disp.lck"):
  print "Driver locked by disp.lck file, no launch."
else:
  try:
    launch(com=command, pbs_name=qsub_name,vmem="2gb",walltime=10)
  except:
    time.sleep(30)
    launch(com=command, pbs_name=qsub_name,vmem="2gb",walltime=10)



