#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

L=[]

for i in range(4):
  L.append(np.loadtxt('/home/wim/Dropbox/ALGOS/TEST_SIMPLE/report%d'%i))


[plt.plot(e[:,0]) for e in L];plt.grid();plt.show()
