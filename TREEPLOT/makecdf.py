import os.path
import spacegrids as sg
import numpy as np
import matplotlib.pyplot as plt
import os

import sys

def get_report(name,qsubs=400):

  HOME = os.environ['HOME']

  path = os.path.join(HOME,'DATA',name,'report%d')

  L=[]

  for i in xrange(qsubs):

    item = np.loadtxt(path%i)
    L.append(item)

  LL = [len(item) for item in L]
  M=max(LL)
  m=min(LL)

  shpe = L[0].shape

  data = np.nan*np.ones((M,qsubs,shpe[1]))

  for i in xrange(qsubs):
    data[:LL[i],i,:] = np.array(L[i])

  return data

data = get_report('M_NOSCAPE_b441',441)

W=data[:,:,0].reshape((data.shape[0],21,21)).astype(np.float64)

X=sg.Ax('X')
Y=sg.Ax('Y')
T=sg.Ax('T')

t=sg.Coord(name='t',value=np.arange(W.shape[0]).astype(np.float64),axis=T)
y=sg.Coord(name='y',value=np.arange(W.shape[1]).astype(np.float64),axis=Y)
x=sg.Coord(name='x',value=np.arange(W.shape[2]).astype(np.float64),axis=X)

F = sg.Field(name='score',value=W,grid=t*y*x)

mF = F/(F/(Y*X))
mF.write()

plt.contourf(mF.value[60,:,:])
plt.colorbar()

plt.show()

