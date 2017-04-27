#!/usr/bin/env python
'''
Dispatcher, launches qsub jobs and sets main config.

'''
import os.path
from runtools import launch, load_config, prep_data
import sys

if os.path.isfile("disp.lck"):
  print "Dispatcher locked by disp.lck file, exiting..."
  sys.exit();

pbs_name = "j_"

#popsize = 40000
#qsubs = 256
#runlen = 6

#popsize = 80000
#qsubs = 289
#runlen = 4

#popsize = 240000
popsize=2000
#qsubs = 256
#qsubs = 2500
qsubs=1
runlen = 4

config = load_config()
qsubs = config['driver']['qsubs']
max_inter_qsubs = config['driver']['max_inter_qsubs']

pbs_name = config['driver']['pbs_name']

#print('disp with popsize %d, qsubs %d, runlen %d'%(popsize, qsubs, runlen))

if len(sys.argv) > 1:
    if sys.argv[1]=='cli':
        use_pbs = False
    else:
        use_pbs = True
else:
    use_pbs = True

prep_data(config) # create local data files

if not use_pbs:
  if qsubs > max_inter_qsubs:
    raise Exception("Refusing to launch %d interactive jobs: exceeds set max of %d"%(qsubs,max_inter_qsubs))

# submit the jobs
for i in xrange(qsubs):
    if use_pbs:
        my_pbs_name = "%s%d.pbs"%(pbs_name,i)
    else:
        my_pbs_name = 'None'
    command = "driver %d %s > pr%d"%(i,my_pbs_name,i)
    launch(com=command, pbs_name=my_pbs_name,vmem="2gb",walltime=10)


