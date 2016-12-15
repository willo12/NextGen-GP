#!/usr/bin/env python

import os.path

import numpy as np
import matplotlib.pyplot as plt
import os

import sys

def get_report(name,qsubs=400):

  HOME = os.environ['HOME']

  path = os.path.join(HOME,'DATA',name,'report%d')

  L=[]

  not_founds = []
  actual_qsubs = 0
  for i in xrange(qsubs):
    try:
      item = np.loadtxt(path%i)
      L.append(item)
      actual_qsubs += 1
    except:
      not_founds.append(i)

  LL = [len(item) for item in L]
  M=max(LL)
  m=min(LL)

  shpe = L[0].shape

  data = np.nan*np.ones((M,qsubs,shpe[1]))
  if not_founds:
    print("warning, qsubs %s not found."%' '.join([str(nf) for nf in not_founds]))

  print("Compiling data from %d qsubs."%actual_qsubs)
  for i in xrange(actual_qsubs):
    data[:LL[i],i,:] = np.array(L[i])

  return data

def normalize(series):

  return series/series[0]




# disk dir names vs display names
names = {'SL_t800k_q484_p80_col1_12':'SL1','JUN_dust':'current'} 

show_min = True
if len(sys.argv)>1:
   if sys.argv[1] == "nomin":
     show_min = False

cutoff = 100

#exps = [ ('SL_t800k_q484_p80_col1_12',484)]
exps = [ ('JUN_dust',625)]

#exps = [ ('TEST',4),('TEST_05',4),('TEST_095',4)]

data = [(get_report(exp,size),exp,size) for exp,size in exps]



#i=0;
#k=0;
#k+=1;print "min data%d: %g"%(k,np.nanmin(data1[:,:,i]))
#k+=1;print "min data%d: %g"%(k,np.nanmin(data2[:,:,i]))

colors = ['k','r','g','b','y','c']

pan = 0

height = 2
width = 2

rows = 2
cols = 2

ax = plt.subplot2grid((height, width), (int(np.floor(pan/cols)), pan%cols) )

i_field=0
for i,(report, exp,size) in enumerate(data):
  color=colors[i%len(colors)]
  plt.plot(np.nanmean(report[:cutoff,:,i_field],1),'--',color=color);
  if exp in names:
    label = names[exp]
  else:
    label = exp
    print("Using raw exp name")

  plt.plot(np.mean(report[:,:,i_field],1),color=color, label=label);
  if show_min is True:
    plt.plot(np.nanmin(report[:cutoff,:,i_field],1),'.',color=color);

#plt.plot(np.nanmean(data2[:,:,i],1),'--',color='r', label="data2");
#plt.plot(np.mean(data2[:,:,i],1),color='r', label="data2");


#plt.plot(np.nanmin(data1[:,:,i],1),'.',color='k', label="data1");
#plt.plot(np.nanmin(data2[:,:,i],1),'.',color='r', label="data2");

#plt.plot(np.mean(data1[:,:,i],1),color='k', label="data1");
#plt.plot(np.mean(data2[:,:,i],1),color='r', label="data2");

plt.xlabel('Iteration')
plt.ylabel('Score')

plt.legend()
plt.grid();plt.show()

