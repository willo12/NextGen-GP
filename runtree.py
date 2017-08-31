#!/usr/bin/env python

import os.path
from math import sqrt
import numpy as np
#import nextgen
import os
from treetools import parse_params_tree_str
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

colors = ['b','r','darkgreen','y','c']

descr = ['obs','$S_1$','$S_2$','$S_3$','$S_4$']

MIS = [('19',-787.0),('17',-707.0),('15c',-624.4),('15a',-579.6), ('13',-499.0), ('11',-424.8),('9',-335.5),('7e',-243.8),('7c',-214.7),('5',-131.4),('1',-11.7)]

HOME = os.getenv("HOME")
fig_name = 'runtree_series'
fig_type = 'eps'

show_f = False
print_flag = False
show_flag = True
iteration = '30'
obs = 'SL'

def find_nearest(array,value):
    return (np.abs(array-value)).argmin()
#    return array[idx]


if __name__ == "__main__":

#  print(sys.path)
  dir_path = os.path.dirname(os.path.realpath(__file__))
#  sys.path.remove(dir_path)

  result_file = 'result'
  show_MIS = False

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
          obs = sys.argv[i]

          i += 1
          if len(sys.argv) > i:
            config['driver']['S_init'] = float(sys.argv[i])

  else:
    print("WARNING: NO TREE SUPPLIED. Using simple test tree.")
    tree='((0,p1)V,(0,p1)V)M'


  fig_name = '_'.join([fig_name,obs,'it%s'%iteration])
  fig_name = '.'.join([fig_name, fig_type])
  out_path = os.path.join(HOME,'Dropbox','paper_algo',fig_name) 

# run tree
  config = load_config()
  (tree, S_init_array, scalars) = parse_params_tree_str(tree)    
  Iffs,result,obs_other,I, tree,ts_factor,startscore_i, S_init, otime, scores = single_run(tree=tree,config = config, S_init_array = S_init_array, scalars = scalars)
# end run tree

  otime = otime*1e-3

  if len(result.shape)==1:
    result.shape = (result.shape[0],1)

  result_plot = result[I[startscore_i:],:]

  for i in range(result_plot.shape[1]):
    print "series %d: max: %g min: %g"%(i,max(result_plot[:,i]) , min(result_plot[:,i]) )

  obs = np.loadtxt("obs")
  if len(obs.shape) == 1:
    obs = obs.reshape( (len(obs),1) )

  fig = plt.figure(1)

  c, = plt.plot(otime,obs[startscore_i:,0],colors[0]); # plot obs

  handles = [c,]

  for i in range(result.shape[1]):
    c, = plt.plot(otime,result_plot[:,i],colors[i+1]); 
    handles.append(c)


  if show_f:
    plt.plot(Iffs[:,0]*1e-3 , Iffs[:,1] ,color='grey')


  if show_MIS:
    for i, mis in enumerate(MIS):
      i_mis = find_nearest(otime,mis[1])
      plt.plot(mis[1],result_plot[i_mis,1] ,'ko')

      plt.text(mis[1]-10,-2.1+0.12*(i%2),mis[0],fontsize=9)

  plt.xlabel('Time kyr BP');plt.ylabel('(scaled)');
  axes = plt.gca()
  axes.set_ylim([-2.2,2.5])  

  dt = (otime[-1]-otime[0])/10
  plt.xticks(np.arange(otime[0],dt,dt))
  plt.legend(handles,[e for e in descr[:result.shape[1]+1]],loc=2)

  plt.grid();

  if out_path and print_flag:
    print("Saving fig to %s"%out_path)
    fig.savefig(out_path)
  else:
    print("Fig not saved.")

  plt.show()

