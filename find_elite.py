#!/usr/bin/env python

import os.path

from treetools import *

import numpy as np
import matplotlib.pyplot as plt
import os
from runtools import stability_test
import sys

def get_elite(name,qsubs=400,i_end=-1, dtype={'names':('score','col','S_init','tree'),'formats':('f4','S10000','S10000','S10000')}  ):

  HOME = os.environ['HOME']
  path = os.path.join(HOME,'DATA',name,'elite%d')
  L=[]

  for i in xrange(qsubs):

    series=np.loadtxt(path%i,dtype= dtype )

    try:
      L.append( (series['score'][i_end],series['S_init'][i_end], series['tree'][i_end], len(series) )   )
    except:
      print "Problem for %s"%series

  return min(L)

def extract_elite(name,qsubs=400,i_end=-1, substring_filter='' , dtype={'names':('score','col','S_init','tree'),'formats':('f4','S10000','S10000','S10000')}, names = None , data_dir = 'DATA', verbose=False):

  """
  Compile score (f4), col (string), S_init (string) and tree (string) from collection of elite NextGen output files.

  Inputs:
	name: directory name where elite files are kept
	qsubs: number of elite files
	i_end: index at which to record data in extraction output
	substring_filter: word that should not occur in tree
	dtype: standard dtype argument to np.loadtxt
	names: names to override names key in dtype
	data_dir: dir where data is kept

  Often, col is not included in dtype argument, yielding just: ('score','S_init','tree').

  Change selection of columns by providing names argument, which corresponds to np.loadtxt dtype names.

  Returns list L of small lists, with elements L[i] being [score, col, S_init, tree] at end of elite file i.

  """

  HOME = os.environ['HOME']
  path = os.path.join(HOME,data_dir,name,'elite%d')

  if verbose:
    print "Extracting from %s"%path

  if names is None:
    names = dtype['names']

  L=[] # will be list of 3 tuples (score, tree string, series length) to be sorted according to score

  for i in xrange(qsubs):

    series=np.loadtxt(path%i,dtype = dtype )
    try:
      L.append( [series[name][i_end] for name in names ]  )
    except:
      print "Problem for %s"%series

  if substring_filter:
    L = [e for e in L if substring_filter not in e[2]]

  return L


def get_elite_stab(name,qsubs=400,i_end=-1, substring_filter='', tolerance=10, dtype={'names':('score','col','S_init','tree'),'formats':('f4','S10000','S10000','S10000')}, names = None , data_dir = os.path.join("DATA"), stab_test=False):

  """
  Get the lowest score numerically stable tree.

  i_end: iteration to look at, default -1
  substring_filter (str): omit trees where this substring_filter is contained in the tree string

  Returns:
    the tree string, stability score, the number of attempts to find a stable tree
  """

  L = extract_elite(name,qsubs,i_end, substring_filter, dtype=dtype, names = names, data_dir = data_dir)  

  L.sort()

  if stab_test:
    for i, e in enumerate(L):
# test tree (e[1]) stability with respect to time step

      S_init_array = np.array([float(ee) for ee in e[1][1:-2].split(',') ])
      print S_init_array

      stab_scores = stability_test(e[2],config=None,ts_factor=[2,8],S_init_array=S_init_array) 
      if stab_scores[0] < tolerance:
        break # found stable enough tree

#  L = zip(scores, trees, stab_scores)
#  print stab_scores

    return e, stab_scores[0],i

  else:

    return L[0][:],0,0



def elite_cluster(name,qsubs=400,i_end=-1, substring_filter='', tolerance=10, dtype={'names':('score','col','S_init','tree'),'formats':('f4','S10000','S10000','S10000')} , d_min=50):

  """
  Cluster elite according to distance
     

  """

  L = extract_elite(name,qsubs,i_end, substring_filter, dtype = dtype)  
  L = [(e[0], "%s %s"%(e[1],e[2]) , Node.from_newick(e[2])   ) for e in L  ]
  L.sort()
  items_min = [] 

  while (len(L)>0):

    items_min.append(L[0])
    L_new = []
    for e in L:

      distances = [e[2].distance(ee[2]) for ee in items_min]
      if (min(distances) > d_min):
        L_new.append(e)

    L = L_new    

  return items_min


def freq_array(x):

  y = np.bincount(x)
  ii = np.nonzero(y)[0]
  return zip(ii,y[ii])


def sanitize_elite(L):

  return [ [e[0], int(e[1].strip("()")), " ".join([e[2],e[3]]) ]  for e in L ]

def best_per_int_par(L, n_int_par=21):

  L = sanitize_elite(L)

  L.sort()

  result = []
  for i_int_par in range(n_int_par):

    B = [e for e in L if e[1] == i_int_par]
    if B:
      result.append([B[0][0], B[0][1] ] )

  return np.array(result)


def elite_scalar_pars(name,qsubs=400,i_end=-1, substring_filter='', tolerance=10, dtype={'names':('score','col','S_init','tree'),'formats':('f4','S10000','S10000','S10000')} , d_min=50, SPACEDIM=2):

  """
  Find frequencies of parameters for lat runs when force_params_node is called due to SCALARPARS>0 in my_ops.c in NextGen, and so all forcing terms are made equal to int_par in trees.

  Example: [(16, 3), (17, 1), (18, 4), (19, 2), (20, 3), (21, 87)]

  Returns frequency array of columns (index as in j_40north_trunc.txt*1e3^6 etc) used across all elite files at i_end (often i_end = -1).

  """

  # start index for each latitude file by latitude: used to produce correct column
  # e.g. for lat 40, start index is 6, corresponding to forcing_files = 'j_40north_trunc.txt*1e3^6-7-8-9-...'
  I_start = {0:14, 10:12, 15:11, 20:10, 25:9, 30:7 ,35:7, 40:6, 45:5, 50:3, 55:2, 60:1, 65:1, 70:1 }

  # improve later with regex and checks
  name_dir = name.split('/')[-1]
  lat = int(name_dir.replace("TEST",""))

  i_start = I_start[lat]

  # obtain long list L (of length qsubs) of small lists of format [ [score, col, S_init, tree], ... ]
  L = extract_elite(name,qsubs,i_end, substring_filter, dtype = dtype , data_dir = os.path.join("DATA","TEST_LAT"))  

  # extract recorded col (parameter index) entry and correct for first columns used up by space variables in output (which counts S0 and S1 for p0 and p1, and then assigns p2 to column i_start in forcing file)
  # list has length qsubs, and records col of tree for each elite file 
  preferred_cols = np.array([ int(e[1].strip("()")) + i_start -SPACEDIM for e in L  ])

  return freq_array(preferred_cols)

def L2struct_array(L,dtype={'names':('score','col','S_init','tree'),'formats':('f4','S10000','S10000','S10000')}):
  """
  Convert list output from extract_elite to structured Numpy array.

  Contracts initial conditions string with tree string. Converts column string to int (after removing brackets).

  Inputs:
	L: list of short lists, each containing score, col, S_init and tree.

  """

  return np.array([ (e[0], int(e[1].strip("()")), "%s %s"%(e[2], e[3]) ) for e in L ], dtype = [('score', 'f4'), ('col','i4'), ('ic_tree','S10000') ])


def L2array(L):
  """

  """

  return np.array([ [ e[0], float(e[1].strip("()"))  ] for e in L ])


if __name__ == "__main__":

  i_end = -1
  qsubs = 256
  dtype={'names':('score','S_init','tree'),'formats':('f4','S10000','S10000')}

  if len(sys.argv) > 1:
    name = sys.argv[1]
    if len(sys.argv) > 2:
      qsubs = int(sys.argv[2])

      if len(sys.argv) > 3:
        i_end = int(sys.argv[3])

        if len(sys.argv) > 4:
          if sys.argv[4] == 'scalar_pars':
            print "using scalar_pars format"
            dtype={'names':('score','col','S_init','tree'),'formats':('f4','S10000','S10000','S10000')}
        
    result = get_elite_stab(name, qsubs,i_end,'' , dtype=dtype, names = None )

    print "".join(["%s "%str(e) for e in result[0]])

   
#    print get_elite(name, qsubs)
