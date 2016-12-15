#!/usr/bin/env python

import os.path
from math import tanh
import numpy as np
import matplotlib.pyplot as plt

# Script to check whether dynamical system formula is correct by comparing results to original tree
# Enter your functions below as lambda functions, and compare output of this script to that of tree with mapf

S0_range = np.arange(-2,2,0.04)
S1_range = np.arange(-2,2,0.04)
ff_range = np.arange(-2,2,0.04)

c_0 = 7.20581e-05 
c_1 = 7.08437e-05 
c_2 = -3.14001e-02

# enter your tree functions here as lambda functions
F = lambda S_0, S_1, f: c_0*( f-c_2*S_0-S_1  )
G = lambda S_0, S_1, f : c_1 *( S_0*(1-tanh(S_1) -S_1 ) )


S0_dot = np.zeros((len(ff_range),len(S1_range),len(S0_range)))
S1_dot = np.zeros((len(ff_range),len(S1_range),len(S0_range)))

for k, ff in enumerate(ff_range):
  for j, S1 in enumerate(S1_range):
    for i, S0 in enumerate(S0_range):
      S0_dot[k,j,i] = F(S0,S1,ff)
      S1_dot[k,j,i] = G(S0,S1,ff)

results = [S0_dot, S1_dot]

names = [r'f',r'$S_1$',r'$S_0$']
titles = [r'$\dot{S}_0$',r'$\dot{S}_1$']

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

    ax = plt.subplot2grid((height, width), (int(np.floor(pan/cols)), pan%cols) )

    plt.contourf(ff_range,S1_range,scale*result[slices].T, levels = levels , cmap=plt.get_cmap('bwr') )

    cs=plt.contour(ff_range,S1_range,scale*result[slices].T, levels = levels ,colors='k')
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

plt.show()



