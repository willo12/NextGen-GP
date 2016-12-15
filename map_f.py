#!/usr/bin/env python

import numpy as np
from nextgen import map_f
import matplotlib.pyplot as plt

from tree2latex import newick2Latex
import os, sys

# Plot the tree as a function of its inputs
# argument 1: tree, in string format
# argument 2: print_flag, boolean. prints fig to file if true
# argument 3: iteration: iteration to use in figure name
# argument 4: observations name to use in figure name

HOME = os.getenv("HOME")
fig_name = 'map_f_S_dot'
fig_type = 'eps'

print_flag = False
show_flag = True
iteration = '30'
obs = 'SL'

if __name__ == "__main__":

  i = 0

  i += 1
  if len(sys.argv) > i:
    tree = sys.argv[i]

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
  else:
    print("Please provide tree.")
    exit()


  fig_name = '_'.join([fig_name,obs,'it%s'%iteration])
  fig_name = '.'.join([fig_name, fig_type])
  out_path = os.path.join(HOME,'Dropbox','paper_algo',fig_name) 

  A=np.zeros((100,100))
  B=np.zeros((100,100))

  print(newick2Latex(tree, index_start=1))

  S0_range = np.arange(-2,2,0.04)
  S1_range = np.arange(-2,2,0.04)
  ff_range = np.arange(-2,2,0.04)

  AA = []
  BB = []

  # loop over forcing
  for ff in ff_range:

    map_f(tree,A,B,np.array([ff,]),forcing_dim=1, S0_start=S0_range[0], S0_end=S0_range[-1], S1_start=S1_range[0], S1_end=S1_range[-1])

    AA.append(A.copy())  # S0
    BB.append(B.copy())  # S1

  # order: indices ff, j, i  <=> ff, S1, S0
  AA = np.array(AA)  # first dim will now be forcing
  BB = np.array(BB)


  results = [AA,BB]

  names = [r'f',r'$S_2$',r'$S_1$']
  ranges = [ff_range,S1_range,S0_range]
  titles = [r'$\dot{S}_1$',r'$\dot{S}_2$']

  i_ff = 0
  i_S1 = 1
  i_S0 = 2

  i_S0_dot = 0
  i_S1_dot = 1

  scale = 1e4

  #levels = np.arange(-6,7,0.5)
  levels = None

  show_vars = {i_S0_dot: [i_ff, i_S0], i_S1_dot: [i_ff, i_S1] }
  #keep_fixed_vars = [i_ff, i_S0]

  # --- start figure ---

  fig = plt.figure(1)

  lbl = ord('a')

  pan = 0

  height = 2
  width = 2

  rows = 2
  cols = 2

  for y, i_dot_show in enumerate(show_vars):
    for x, i_fixed in enumerate(show_vars[i_dot_show]):
      result = results[i_dot_show]
      i_slice = 51
      slices = [slice(None)]*3
      slices[i_fixed] = i_slice

      labels = [e for i, e in enumerate(names) if i != i_fixed]
      show_ranges = [e for i, e in enumerate(ranges) if i != i_fixed]

      ax = plt.subplot2grid((height, width), (int(np.floor(pan/cols)), pan%cols) )

      plt.contourf(show_ranges[0],show_ranges[1],scale*result[slices].T, levels = levels , cmap=plt.get_cmap('bwr') )

      cs=plt.contour(show_ranges[0],show_ranges[1],scale*result[slices].T, levels = levels ,colors='k')
      plt.clabel(cs,fmt='%1.1f');

#plt.colorbar()

      plt.xlabel(labels[0])

      if y == len(show_vars)-1:
        plt.xlabel(labels[0])
      else:
        plt.tick_params(axis='x',labelbottom='off')
        plt.xlabel('')

      plt.ylabel(labels[1])
      if x == 0:
        plt.ylabel(labels[1])
      else:
        plt.tick_params(axis='y',labelleft='off')



      plt.title(' '.join( ['(%s)'%chr(lbl), titles[i_dot_show],'at %s = %s'%(names[i_fixed], str(S0_range[i_slice])    )  ] ) )

#plt.plot(ff_range,BB[:,50,0].T)

      lbl += 1
      pan += 1

  if out_path and print_flag:
    print("Saving fig to %s"%out_path)
    fig.savefig(out_path)
  else:
    print("Fig not saved.")

  if show_flag:
    plt.show()



