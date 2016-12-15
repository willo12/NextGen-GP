#!/usr/bin/env python

import os.path

import numpy as np
import matplotlib.pyplot as plt
import os
from runtools import stability_test
import sys

def get_elite(name,qsubs=400,i_end=-1):

  HOME = os.environ['HOME']

  path = os.path.join(HOME,'DATA',name,'elite%d')

  L=[]

  for i in xrange(qsubs):

    series=np.loadtxt(path%i,dtype={'names':('score','tree'),'formats':('f4','S10000')})
    try:
      L.append((series['score'][i_end],series['tree'][i_end],len(series)))
    except:
      print "Problem for %s"%series

  return min(L)

def get_elite_stab(name,qsubs=400,i_end=-1, substring_filter=''):

  """
  Get the lowest score numerically stable tree.

  i_end: iteration to look at, default -1
  substring_filter (str): omit trees where this substring_filter is contained in the tree string

  Returns:
    the tree string, stability score, the number of attempts to find a stable tree
     

  """

  HOME = os.environ['HOME']

  path = os.path.join(HOME,'DATA',name,'elite%d')

  L=[] # will be list of 3 tuples (score, tree string, series length) to be sorted according to score

  for i in xrange(qsubs):

    series=np.loadtxt(path%i,dtype={'names':('score','tree'),'formats':('f4','S10000')})
    try:
      L.append((series['score'][i_end],series['tree'][i_end],len(series)))
    except:
      print "Problem for %s"%series

  if substring_filter:
    L = [e for e in L if substring_filter not in e[1]]
    

  L.sort()

  for i, e in enumerate(L):
# test tree (e[1]) stability with respect to time step 
    stab_scores = stability_test(e[1],config=None,SPACEDIM=2,ts_factor=[2,8]) 
    if stab_scores[0] < 10:
      break # found stable enough tree

#  L = zip(scores, trees, stab_scores)
#  print stab_scores

  return e, stab_scores[0],i

if __name__ == "__main__":

  if len(sys.argv) > 1:
    name = sys.argv[1]
    if len(sys.argv) > 2:
      qsubs = int(sys.argv[2])
    else:
      qsubs = 256

    if len(sys.argv) > 3:
        i_end = int(sys.argv[3])
    else:
        i_end = -1

    print get_elite_stab(name, qsubs,i_end,'')
  

#    print get_elite(name, qsubs)
