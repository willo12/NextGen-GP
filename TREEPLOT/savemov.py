import matplotlib
matplotlib.use("Agg")

from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import os.path
import spacegrids as sg
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import animation
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

def print_there(x, y, text):
     sys.stdout.write("\x1b7\x1b[%d;%df%s\x1b8" % (x, y, text))
     sys.stdout.flush()


data = get_report('M_PS_b100',100)
maxtime = data.shape[0]


lens = []
for r in xrange(data.shape[1]):

  nans = np.argwhere(np.isnan(np.squeeze(data[:,r,0])))
  if len(nans)>0:    
    i_nan = nans[0][0]
    step = float(i_nan-1)/maxtime
    Icoarse = np.arange(0.,float(i_nan))
    Ifine = np.arange(0.,float(i_nan-1),step)

#    print Icoarse
#    print Ifine
 
    

    fI = interp1d(Icoarse,np.squeeze(data[:i_nan,r,0]))(Ifine)
    if len(fI) > data.shape[0]:
      data[:,r,0] = fI[:-1]
    else:
      data[:,r,0] = fI[:]

  else:

    pass


#sys.exit()

W=data[:,:,0].reshape((data.shape[0],10,10)).astype(np.float64)

X=sg.Ax('X')
Y=sg.Ax('Y')
T=sg.Ax('T')

t=sg.Coord(name='t',value=np.arange(W.shape[0]).astype(np.float64),axis=T)
y=sg.Coord(name='y',value=np.arange(W.shape[1]).astype(np.float64),axis=Y)
x=sg.Coord(name='x',value=np.arange(W.shape[2]).astype(np.float64),axis=X)

tfine=sg.Coord(name='t',value=np.arange(0,W.shape[0],0.3).astype(np.float64),axis=T)

F = sg.Field(name='score',value=W,grid=t*y*x)

mF = (F/(F/(Y*X))).regrid(tfine*y*x)

#print len(tfine.value)
#mF.write()

tlen=mF.shape[0]

Y, X = (y*x).meshgrid()


FFMpegWriter = animation.writers['mencoder']
metadata = dict(title='Movie Test', artist='Matplotlib',comment='Movie support!')
writer = FFMpegWriter(fps=30,metadata=metadata)

fig = plt.figure()
#ax = fig.gca(projection='3d')
ax = fig.gca()

#  cont = plt.pcolor(mF.value[0,:,:])

with writer.saving(fig,'writer_test.mp4',100):

  for i in range(tlen):
    print_there(0,0,"%d"%i)
    ax.cla()
    cont = plt.contourf(mF.value[i,:,:], cmap=cm.coolwarm)
#    surf = ax.plot_surface(X, Y, mF.value[i,:,:], rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    writer.grab_frame()


