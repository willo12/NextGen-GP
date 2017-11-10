#!/usr/bin/env python

import re
from math import sqrt
import numpy as np
import warnings

def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str


# ---- Latex tools -----

def bracket_block(text):
  """
  finds start and end of the first bracketed block
  """

  istart=text.find('(')

  cursor=istart

  if cursor != -1:

    counter=1
    while (cursor<len(text)) and counter>0:
      cursor+=1
      if text[cursor]=='(':
         counter += 1
      elif text[cursor]==')':
         counter -= 1  

  return istart, cursor


def nw_split(text):
  
  istart,iend=bracket_block(text)
  if istart == -1:
    return text.split(',')
  else:
    # include the opcode at the end
    result=[text[istart:iend+2],]

    if iend<len(text)-2:
      # following text may contain further brackets
      result = result + nw_split(text[iend+3:])

    if istart>0:
      # preceding text does not contain brackets by construction
      # so can be readily included without further special splitting
      result = text[:istart-1].split(',') + result

    return result



def node2Latex(node, index_start=1):


  pieces = {}
  output = '**** LATEX *******\n\n'
  
  for i, child in enumerate(node.children):
    output += '$ \\dot{S_%d} = %s $ \\newline \\newline \n'%(i+1,child.to_latex())



  return output


def newick2Latex(text, index_start=1, SPACEDIM=2):

  try:
    node=Node.from_newick(text)
    SPACEDIM = len(node.children)
  except:
    warnings.warn('Tree cannot be parsed using Node.from_newick. Leaving SPACEDIM = %d'%SPACEDIM)
    SPACEDIM = 2


  if ' ' in text:
    raise ValueError("No spaces allowed in tree string %s"%text)

  pieces = {}
  output = '**** LATEX *******\n\n'
  
  output += '$ \\dot{S} = %s $ \\newline \\newline \n'%newick2human(text,pieces, index_start=index_start, SPACEDIM=SPACEDIM)

  coef_list = pieces.items()
  coef_list.sort()

  # list coeficients (constants) in Latex
  for coef_tuple in coef_list:
    output += '$ %s = %s $ \\newline \n'%(coef_tuple[0],coef_tuple[1])


  output += '\n**** END LATEX *******\n\n'

  # list them again for copy paste into Python
  for coef_tuple in coef_list:
    output += '%s = %s\n'%(coef_tuple[0],coef_tuple[1])


  return output

def newick2human(text, pieces = None, SPACEDIM=2, index_start=1,dtype=float):

  strings = {'V':'[' + ',\n '.join(['%s']*SPACEDIM) + ']','A':'%s + %s','S':'%s - %s','M':'(%s) (%s)','D':'(%s) / (%s)','Q': '\sqrt{%s}','I':'%s $ for $ %s>0 $ and $ %s $ otherwise $ ', 'T':'tanh(%s)'}

  leaf_ops = {'Q':'\sqrt','E':'e','O':'1/','T':'tanh','L':'log'}

  if pieces == None:
    pieces = {}

  if '(' in text:
    opcode = text[-1]
    childrentext = text[1:-2]

    childtexts = nw_split(childrentext)

    if (opcode == 'Q') and ('(' not in childtexts[0]) and ('p' not in childtexts[0]):
      newvarname='c_%d'%len([el for el in pieces if 'c' in el])
      pieces[newvarname] = '%.2e'%sqrt(dtype(childtexts[0]))
      return newvarname
#    print strings[opcode]
#    print tuple([newick2human(ct) for ct in childtexts])

    if (opcode == 'I'):
#      print [el for el in pieces if 'a' in el], len([el for el in pieces if 'a' in el])
      newvarname='a_%d'%len([el for el in pieces if 'a' in el])
      pieces[newvarname] = 'placeholder'

      pieces[newvarname] = strings[opcode]%tuple([newick2human(text=ct,pieces=pieces,SPACEDIM=SPACEDIM, index_start=index_start,dtype=dtype) for ct in [childtexts[1], childtexts[0],childtexts[2] ]  ])

      return newvarname

    return strings[opcode]%tuple([newick2human(text=ct,pieces=pieces,SPACEDIM=SPACEDIM, index_start=index_start,dtype=dtype) for ct in childtexts])

  else:
 
    if text in leaf_ops:
      return leaf_ops[text]

    if 'p' in text:
      if text[1:] in [str(e) for e in range(SPACEDIM)]:
        return 'S_%d'%(int(text[1:])+index_start)
      else:
        return 'f_%d'%(int(text[1:])-SPACEDIM+index_start)
    else:
      newvarname='c_%d'%(len([el for el in pieces if 'c' in el])+index_start)
      pieces[newvarname] = '%.5e'%dtype(text)
      return newvarname

# -------- end Latex tools -----


def parse_params_tree_str(tree_str, dtype=float):

  S_init_array = None
  scalars = None

  tree_str = re.sub('\s+',' ',tree_str).strip()

  L = tree_str.split(" ")
 
  if L:
   tree = L.pop()
  else:
   raise Exception("Provide at least a tree string.")

  if L:
    S_init_array = np.array( [dtype(e) for e in  L.pop().strip("()").split(",") ] )

  if L:
    scalars = np.array( [dtype(e) for e in  L.pop().strip("()").split(",") ] )

  return (tree, S_init_array, scalars )
    

def protected_mem(i,numpar,spacedim=6):

  if i<spacedim:
    i=spacedim
  elif (i>numpar-7):
    i=numpar-7

  return i;


def nopFun(reg):
  return 0

def addFun(reg,x,y):

  return x + y

def subFun(reg,x,y):

  return x - y

def mulFun(reg,x,y):

  return x * y

def ifFun(reg,x,y,z):

  if (x>0):
    return y
  else:
    return z

def isgreaterFun(reg,x,y):

  if (x>y):
    return 1
  else:
    return 0


def iseqFun(reg, x,y):

  if (x==y):
    return 1
  else:
    return 0


def copyFun(reg,x,y):

  reg[protected_mem(x,len(reg))] = reg[y%len(reg)]  

  return x

def storeFun(reg,x,y):

  reg[protected_mem(y,len(reg))] = x

  return x

def ldaFun(reg,x):

  return reg[x%len(reg)] 


def lda_xFun(reg,x):

  return reg[(x+reg[len(reg)-6])%len(reg)] 

def ldxFun(reg, x,y):

  reg[len(reg)-6] = reg[y%len(reg)]    

  return x

def sta_xFun(reg,x,y):

  reg[protected_mem(y+reg[len(reg)-6],len(reg))] = x

  return x

def stxFun(reg,x,y):

  reg[protected_mem(y,len(reg))] = reg[len(reg)-6]  

  return x

def inxFun(reg,x):

  reg[len(reg)-6] += 1
#  if (reg[len(reg)-6] >= len(reg)):
#    reg[len(reg)-6] -= len(reg)

  reg[len(reg)-6] = reg[len(reg)-6]%len(reg)

  return x

def dexFun(reg,x):

  reg[len(reg)-6] -= 1
  if (reg[len(reg)-6] < len(reg)):
    reg[len(reg)-6] = 0

  return x


def reverse_sign(op_sign):
 
  if (op_sign == '+'):
   
    return '-'
  else:
    return '+'

funcs = {'A':addFun, 'S':subFun,'M':mulFun,'I':ifFun, 'G':isgreaterFun, 'E':iseqFun, 'N':copyFun,'R':storeFun,
         'L':ldaFun, 'l':lda_xFun, 'X':ldxFun, 's':sta_xFun, 'x':stxFun, 'i':inxFun, 'd':dexFun}


additive = {'A':'+','S':'-'}
multiplicative = {'M':'','S':'/'}

functions = {'T':'tanh(%s)','Q':'\sqrt(%s)'}


class Node(object):


#  strings = {'V':'[' + ', '.join(['%s']*SPACEDIM) + ']','A':'%s + %s','S':'%s - %s','M':'(%s) (%s)','D':'(%s) / (%s)','Q': '\sqrt{%s}','I':'%s $ for $ %s>0 $ and $ %s $ otherwise $ ', 'T':'tanh(%s)'}

#  leaf_ops = {'Q':'\sqrt','E':'e','O':'1/','T':'tanh','L':'log'}






  def __repr__(self):
    return self.op

  def __init__(self,op,children=[], op_table = {'A':'+'}): 

    self.op = op
    self.children = children


  def __call__(self, reg):
    return funcs[self.op](reg,*[c(reg) for c in self.children])

  @classmethod
  def from_newick(cls,text,dtype=float):

    if '(' in text:
      opcode = text[-1]
      childrentext = text[1:-2]

      childtexts = nw_split(childrentext)
    
      return cls(opcode, [cls.from_newick(ct,dtype=dtype) for ct in childtexts])
    else:
      if 'p' in text:
        return ParamNode(int(text[1:]))
      elif 'T' in text:
        return Node(text,[])
      else:
        return ConstNode(dtype(text))

  def walk(self,indent=0):

    print (' '*indent) + self.__repr__()
    for c in self.children:
      c.walk(indent+1)


  def mbrack(self):

    if self.op in functions:
      return self.to_latex()

    if self.__class__ is TermNode:
      return self.to_latex()

    if self.children:
      return '(%s)'%self.to_latex()
    else:
      return self.to_latex()


  def to_latex(self):

    
    if self.op in additive:
      op_sign = additive[self.op]
#      if isinstance(self.children[1],ConstNode):
#        if (self.children[1].value < 0):
#          op_sign = reverse_sign(op_sign)

      return ' '.join([self.children[0].to_latex(), op_sign ,  self.children[1].to_latex() ] ).replace('+ -','- ')
    elif self.op in multiplicative:
      return ' '.join([self.children[0].mbrack(), multiplicative[self.op] ,  self.children[1].mbrack() ] )
    elif self.op in functions:
      return functions[self.op]%self.children[0].to_latex()
    else:
      return '%s'%self


  def is_equal(self, other):

    if (self.op != other.op):
      return False

    elif (len(self.children) > 0):
      return  reduce( lambda x,y: x and y, [e.is_equal(other.children[i]) for i, e in enumerate(self.children) ] )
    else:

      return True


  def distance(self, other, d=0):

    if (other is None):

      d += 1 
      for c in self.children:
        d =  c.distance(None,d)
 
      return d


    if (self.op != other.op):
      d += 1

    if (len(other.children) > len(self.children) ):
      left = other
      right = self

    else:
      left = self
      right = other   

    m = min( len(left.children), len(right.children) )
 
    for i in range(m):

      d = left.children[i].distance(right.children[i],d)

    for i in range( len(left.children) -m   ):
      d = left.children[i].distance(None,d)

    return d


class TermNode(Node):
  """
  Terminal node
  """

  def distance(self, other, d=0):

    if (other is None):
      
      return d+1

    elif (not isinstance(other,self.__class__)):

      d += 1 # account for node

      for e in other.children:
        d = e.distance(None,d)

    elif not self.is_equal(other):
      d += 1

    return d

class ParamNode(TermNode):


  def __repr__(self):
    return "S_%d"%(self.index+1)

  def __init__(self,index ): 

    self.op = 'P'
    self.children = []
    self.index = index;

  def __call__(self, reg):

    try:
      return reg[self.index]
    except:
      raise Exception("error fetching p node: index=%d vs len(reg)=%d\n"%(self.index,len(reg)))


  def is_equal(self, other):

    return isinstance(other,self.__class__) and (self.index == other.index)



class ConstNode(TermNode):

  def __repr__(self, dtype=float):
    if dtype is float:
      return latex_float(self.value)
    else:
      return '%d'%self.value

  def __init__(self,value ): 

    self.op = 'C'
    self.children = []
    self.value = value;

  def __call__(self, reg):
    return self.value

  def is_equal(self, other):

    return isinstance(other,self.__class__) and (self.value == other.value)


