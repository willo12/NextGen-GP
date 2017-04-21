#!/usr/bin/env python

import sys
from random import random,randint,choice
from treetools import *
import curses
import time
from curses import wrapper
import copy

stdscr = curses.initscr()

def gridgame(p, aliases, random=False):
  direction_aliases = {0:"up",1:"down",2:"left",3:"right"}
 
  # Board size
  max=(3,3)
  
  # Remember the last move for each player
  lastmove=[-1,-1]
  
  # Remember the player's locations
  location=[[randint(0,max[0]),randint(0,max[1])]]
  
  # Put the second player a sufficient distance from the first
  location.append([(location[0][0]+2)%4,(location[0][1]+2)%4])

  update_board(location, location)     

  # Maximum of 50 moves before a tie
  for o in range(50):

    location_last = copy.deepcopy(location)
  
    # For each player
    for i in range(2):
      locs=location[i][:]+location[1-i][:]
      locs.append(lastmove[i])

      if random:
        locs.append(randint(0,1))

      move=p[i](locs)%4

      msg = "%s moved %s                        "%(aliases[i],direction_aliases[move])
      stdscr.addstr(7+i,10, "%s "%msg)
      stdscr.refresh()

      # You lose if you move the same direction twice in a row
      if lastmove[i]==move:
#        print "twice last move %d"%i 
        msg = "%s MOVED TWICE %s                      "%(aliases[i],direction_aliases[move])
        stdscr.addstr(7+i,10, "%s "%msg)
        stdscr.refresh()

        return 1-i
      lastmove[i]=move
      if move==0: 
        location[i][0]-=1
        # Board wraps
        if location[i][0]<0: location[i][0]=0
      if move==1: 
        location[i][0]+=1
        if location[i][0]>max[0]: location[i][0]=max[0]
      if move==2: 
        location[i][1]-=1
        if location[i][1]<0: location[i][1]=0
      if move==3: 
        location[i][1]+=1
        if location[i][1]>max[1]: location[i][1]=max[1]

      update_board(location_last, location)           
      curses.napms(150)
 
    # If you have captured the other player, you win
      if location[i]==location[1-i]: 
        msg = "%s CAPTURED %s moving %s                       "%(aliases[i],aliases[1-i],direction_aliases[move])
        stdscr.addstr(7+i,10, "%s "%msg)
        stdscr.refresh()

        return i

  return -1


def update_board(location_last, location):

  me_last=tuple(location_last[1])
  others_last=[tuple(location_last[0]),]    

  me=tuple(location[1])
  others=[tuple(location[0]),]    

#  stdscr.addstr(20,10, "%s %s"%(location, location_last))

  for i in range(4):
    for j in range(4):
      i_color_pair = 1
      if (i,j)==me:
        if (i,j) in others:
          square = "*"
        else:
          i_color_pair = 1
          square = "O"
      elif (i,j) in others:
        i_color_pair = 1
        square = "X"
      elif (i,j) == me_last:
        i_color_pair = 2
        square = "O"      
      elif (i,j) in others_last:
        i_color_pair = 2
        square = "X"

      else:
        square = "."

      stdscr.addstr(i+10,j+10, "%s "%square, curses.color_pair(i_color_pair))

  stdscr.refresh()

class humanplayer:
  def __call__(self,board):

    # Get my location and the location of other players
#    me=tuple(board[0:2])
#    others=[tuple(board[x:x+2]) for x in range(2,len(board)-1,2)]
    
    # Display the board

    x = 0

    while True:
 
#      update_board(me, others, stdscr)
#      stdscr.nodelay(1)
      while True:
        x = stdscr.getch()
 
        if x == curses.KEY_UP:
          return 0
        elif x == curses.KEY_DOWN:
          return 1
        elif x == curses.KEY_LEFT:
          return 2
        elif x == curses.KEY_RIGHT:
          return 3
        elif x == ord('q'):
          curses.endwin()
          exit(0)

      
        time.sleep(0.1)

def play_game(tree, random=False):

  aliases = {-1:'tie',0:'computer', 1:'you'}

  result = aliases[gridgame([tree, humanplayer()], aliases, random)  ] + " won."

#  stdscr.addstr(15,10, "%s "%result)  
#  stdscr.refresh()

  return result

def main(stdscr, argv):

  if len(argv) > 1:
    treestr = argv[1]
  else:
    raise Exception("Provide tree");

  result = "Welcome."

  tree = Node.from_newick(treestr)  

  stdscr.nodelay(1)

#  curses.start_color() # taken care of by wrapper
  curses.init_pair(1, curses.COLOR_WHITE, curses.COLOR_BLACK)
  curses.init_pair(2, curses.COLOR_RED, curses.COLOR_BLACK)

# turn off input echoing
  curses.noecho()

# respond to keys immediately (don't wait for enter)
  curses.cbreak()

# map arrow keys to special values
  stdscr.keypad(True)

  stdscr.clear()
  stdscr.border(0)

  stdscr.addstr(17,10, "p to play game ")  
  stdscr.addstr(18,10, "q to quit ")  
  stdscr.refresh()

  x = 0
  while x != ord('q'):
    x = stdscr.getch()

    if x == ord('p'):
      stdscr.addstr(15,10, "                      ")  
      result = play_game(tree, True);    
      stdscr.addstr(15,10, "%s "%result)  
      stdscr.refresh()


  curses.endwin()



if __name__ == "__main__":

  curses.wrapper(main, sys.argv)
 

