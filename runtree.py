#!/usr/bin/env python

import os.path
from math import sqrt
import numpy as np
#import nextgen
import os
from treetools import find_initial_conditions
from runtools import launch, getf, putf, load_config, single_run
from my_data import synth, epica, glacial
import h5py as hdf
import sys
import time
import matplotlib.pyplot as plt
import pickle
from configobj import ConfigObj, ConfigObjError, flatten_errors
from validate import Validator

#from spacegrids import dlmread

MIS = [('19',-787.0),('17',-707.0),('15c',-624.4),('15a',-579.6), ('13',-499.0), ('11',-424.8),('9',-335.5),('7e',-243.8),('7c',-214.7),('5',-131.4),('1',-11.7)]

HOME = os.getenv("HOME")
fig_name = 'runtree_series'
fig_type = 'eps'

print_flag = False
show_flag = True
iteration = '30'
obs = 'SL'

if __name__ == "__main__":

#  print(sys.path)
  dir_path = os.path.dirname(os.path.realpath(__file__))
#  sys.path.remove(dir_path)

  result_file = 'result'

  print('run locally from %s'%dir_path)

  i = 0

  i += 1
  if len(sys.argv) > i:
    tree = sys.argv[i].strip()

    i += 1
    if len(sys.argv) > i:
      if (sys.argv[i] == 'print'):
        print_flag = True     

      i += 1
      if len(sys.argv) > i:
        iteration = sys.argv[i]

        i += 1
        if len(sys.argv) > i:
          obs = int(sys.argv[i])

          i += 1
          if len(sys.argv) > i:
            config['driver']['S_init'] = float(sys.argv[i])

  else:
    print("WARNING: NO TREE SUPPLIED. Using simple test tree.")
    tree='((0,p1)V,(0,p1)V)M'


  fig_name = '_'.join([fig_name,obs,'it%s'%iteration])
  fig_name = '.'.join([fig_name, fig_type])
  out_path = os.path.join(HOME,'Dropbox','paper_algo',fig_name) 

  popsize = 1

  config = load_config()

  (tree, S_init_array) = find_initial_conditions(tree)    

  Iffs,result,obs_other,I, tree,ts_factor,startscore_i, S_init, otime, scores = single_run(tree=tree,config = config, S_init_array = S_init_array)





  obs = np.loadtxt("obs")
  if len(obs.shape) == 1:
    obs = obs.reshape( (len(obs),1) )

#  print('Score: %g, Starting obs time: %d, I[startscore_i]: %d, t0: %d, S_init: %g, ts_factor: %s'%(scores[0],obs[0,0] , I[startscore_i] ,Iffs[0,0] , S_init, str(ts_factor)  ) )

#  print ("startscore_i: %d %d"%(startscore_i, result[I[startscore_i]]))
 
  #pickle.dump({'result':result, 'Iffs':Iffs,'I':I,'startscore_i':startscore_i,'otime':otime,'obs':obs},open("result.p","wb"))

#  F = lambda x: 2*(x - np.nanmin(x))/(np.nanmax(x)-np.nanmin(x))
#  A=np.loadtxt('/home/wim/PROJECTS/paper_glac/edc-dust2008.txt');A[:,2] = F(A[:,2])

  fig = plt.figure(1)

#  plt.plot(otime*1e-3,obs[startscore_i:,1],'g'); # note: obs has been re-indexed
#  plt.plot(-A[:,1],A[:,2],'g')
  c1, = plt.plot(otime*1e-3,obs[startscore_i:,0],'b'); # plot obs
  c2, = plt.plot(otime*1e-3,result[I[startscore_i:],0],'r'); 
  c3, = plt.plot(otime*1e-3,result[I[startscore_i:],1],'g'); 
#  plt.plot(otime*1e-3,result[I[startscore_i:],2],'y'); 


  for mis in MIS:
    plt.plot(mis[1],-1.8,'ko')

    plt.text(mis[1]-10,-2.1,mis[0])

  plt.xlabel('Time kyr BP');plt.ylabel('(scaled)');
  axes = plt.gca()
  axes.set_ylim([-2.2,2.5])  

  dt = 80
  plt.xticks(np.arange(-800,dt,dt))
  plt.legend([c1,c2,c3],['obs','$S_1$','$S_2$'],loc=2)

  plt.grid();

  if out_path and print_flag:
    print("Saving fig to %s"%out_path)
    fig.savefig(out_path)
  else:
    print("Fig not saved.")

  plt.show()

