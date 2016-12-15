
from math import sqrt


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

def parse_newick(text):

  if '(' in text:
    opcode = text[-1]
    childrentext = text[1:-2]

    childtexts = nw_split(childrentext)
    
    return node(code2wrapper[opcode],[parse_newick(ct) for ct in childtexts])
  else:
    if 'p' in text:
      return paramnode(int(text[1:]))
    else:
      return constnode(float(text))


def newick2Latex(text, index_start=0):

  pieces = {}
  output = '**** LATEX *******\n\n'
  
  output += '$ \\dot{S} = %s $ \\newline \\newline \n'%newick2human(text,pieces, index_start=index_start)

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

def newick2human(text, pieces = None, SPACEDIM=2, index_start=0):

  strings = {'V':'[%s , %s]','A':'(%s + %s)','S':'(%s - %s)','M':'%s %s','D':'%s / %s','Q': '\sqrt{%s}','I':'%s $ for $ %s>0 $ and $ %s $ otherwise $ '}

  leaf_ops = {'Q':'\sqrt','E':'e','O':'1/','T':'tanh','L':'log'}

  if pieces == None:
    pieces = {}

  if '(' in text:
    opcode = text[-1]
    childrentext = text[1:-2]

    childtexts = nw_split(childrentext)

    if (opcode == 'Q') and ('(' not in childtexts[0]) and ('p' not in childtexts[0]):
      newvarname='c_%d'%len([el for el in pieces if 'c' in el])
      pieces[newvarname] = '%.2e'%sqrt(float(childtexts[0]))
      return newvarname
#    print strings[opcode]
#    print tuple([newick2human(ct) for ct in childtexts])

    if (opcode == 'I'):
#      print [el for el in pieces if 'a' in el], len([el for el in pieces if 'a' in el])
      newvarname='a_%d'%len([el for el in pieces if 'a' in el])
      pieces[newvarname] = 'placeholder'

      pieces[newvarname] = strings[opcode]%tuple([newick2human(ct,pieces, index_start=index_start) for ct in [childtexts[1], childtexts[0],childtexts[2] ]  ])

      return newvarname

    return strings[opcode]%tuple([newick2human(ct,pieces, index_start=index_start) for ct in childtexts])

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
      pieces[newvarname] = '%.5e'%float(text)
      return newvarname



