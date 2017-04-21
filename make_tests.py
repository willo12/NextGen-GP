#!/usr/bin/env python

import numpy as np
from runtools import *

from random import random,randint,choice
from copy import deepcopy
from math import log

def hiddenfunction(x,y):
    return x**2+2*y+3*x+5

def buildhiddenset():
  rows=[]
  for i in range(200):
    x=randint(0,40)
    y=randint(0,40)
    rows.append([x,y,hiddenfunction(x,y)])

  return np.array(rows)

if __name__ == "__main__":

  rows = buildhiddenset()

  write_array('rows',rows, fmt='%d')


