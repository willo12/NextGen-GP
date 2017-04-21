
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <compile_options.h>
#include <states.h>
#include <node_structures.h>
#include <fields.h>
#include <basic_ops.h>
#include <my_ops.h>
#include <tree_io.h>
#include <model.h>
#include <c_nextgen.h>

// local preprocessor directives

#define RANDOM
#define INTERP_FORCING  // interpolate forcing while solving

#define HILLTHRESH 0.02
#define MAXHILL 8 // max hill climbing steps

#define RNDGRAIN 0.005
#define MAXSEARCH 0

#define BLOWUP 4  // values for S_i beyond which the solution is considered unusable
#define MAXCHANGE 0.1

char* instrset = "ASMIGE"; /* used at runtime */
char* instrset_rare = "T"; /* used at runtime */
char* instrsetswap = "ASMGE"; /* instr that can swap */
char* instrsetzero = "\0"; /* function subleaf nodes */
/* char* instrsetzero = "QEOTL";  */

char* instrsetstem = "V";
char* instrsetsuperleaf = "V";
char* instrsetleaf = "CP";

char const_op_char = 'C';
char par_op_char = 'P';

int numpar;

// SPECIFIC TO ops_muldim!! CHANGE LATER
double buffers[STACKSIZE][SPACEDIM];

// local function definitions


int init_reg(Experiment Exp);
Population play_all_others(Population pop, int i_ns, Experiment Exp);
int gridgame(Node *players[], Experiment Exp);

// functions

int model_init(Experiment Exp)
{
  init_tables();
  init_reg(Exp);

  return 0;
}


int free_experiment(Experiment E)
{

//  free_int_field(E.rows);

  free_state(E.S);
  free_state(E.S_i);

  return 0;
}



int init_reg(Experiment Exp)
{
  int i;

//  int ffs1 = Exp.rows.dims.cols;

  for (i=0;i<MAXPAR;i++)
  {
    reg[i] = 0;
  }

#ifdef RANDOM
  numpar = 6;
#else
  numpar = 5;
#endif

  return check_numpar(numpar);

}

Population play_all_others(Population pop, int i_ns, Experiment Exp)
{
  int i, winner;

  Node *players[2];

  players[0] = pop.pop[i_ns].node;
  
  for (i=0;i<pop.popsize;i++) 
  {
//    pop.pop[i_ns].score += gridtraining(pop.pop[i_ns].node, Exp);

    if (i==i_ns)
    {
      continue;
    }

    players[1] = pop.pop[i].node;
    winner = gridgame(players, Exp);

    if (winner == 0)
    {
      pop.pop[i].score += 2;
    }
    else if (winner == 1)
    {
      pop.pop[i_ns].score += 2;
    } 
    else if (winner == -1)
    {
      pop.pop[i].score += 1;
      pop.pop[i_ns].score += 1;
    }
  }

  return pop;
}

Population score_pop(Population pop_old,Population pop, Experiment Exp, int compgridsize)
{
#ifndef INLINESCORING
  int i;
  int *current = make_int(0);

  for (i=0;i<pop.popsize;i++) /* fill pop_next. i is only incremented once new_ns.node is succesfully added, not discarded.  */
  {
    pop.pop[i].score = 0.0;
  }

  for (i=0;i<pop.popsize;i++) /* fill pop_next. i is only incremented once new_ns.node is succesfully added, not discarded.  */
  {
    pop = play_all_others(pop, i, Exp);
  }

  for (i=0;i<pop.popsize;i++) /* fill pop_next. i is only incremented once new_ns.node is succesfully added, not discarded.  */
  {
    pop.toterror += pop.pop[i].score;
    if (pop.pop[i].score < pop.minerror)
    {
      pop.minerror = pop.pop[i].score;
      pop.i_min = i; // keep track of best score
    }
  }

  free(current);

#endif

  return pop;
}


int gridgame(Node *players[], Experiment Exp)
{
  int i,o;
  int max_board[2] = {3,3};
  int last_move[2] = {-1,-1};

  int location[4];
  int move;

  State S_result = Exp.S;

  location[0] = rand()%(max_board[0]+1);
  location[1] = rand()%(max_board[1]+1);

  location[2] = (location[0]+2)%4;
  location[3] = (location[1]+2)%4;

  for (o=0;o<50;o++)
  {
    for (i=0;i<2;i++)
    {
      reg[0] = location[2*i];
      reg[1] = location[2*i+1];

      reg[2] = location[2*(1-i)];
      reg[3] = location[2*(1-i)+1];
      reg[4] = last_move[i];

      evaluate_tree(players[i], S_result);
      move = S_result.data[0]%4;

      if (last_move[i] == move)
      {
        return 1-i;
      }
      last_move[i] = move;    

      if (move==0)
      {
        location[2*i] -= 1;
        if (location[2*i] < 0)
        {
          location[2*i] = 0;
        }
      }

      else if (move==1)
      {
        location[2*i] += 1;
        if (location[2*i] > max_board[0])
        {
          location[2*i] = max_board[0];
        }
      }
      else if (move==2)
      {
        location[2*i+1] -= 1;
        if (location[2*i+1] < 0)
        {
          location[2*i+1] = 0;
        }
      }

      else if (move==3)
      {
        location[2*i+1] += 1;
        if (location[2*i+1] > max_board[1])
        {
          location[2*i+1] = max_board[1];
        }
      }

      if ( (location[2*i] == location[2*(1-i)] ) && (location[2*i+1] == location[2*(1-i)+1] ) )
      {
        return i;
      }

    }
  }

  return -1;
}


double get_score(Node *newtree, State S_i, Experiment Exp)
{

  return 0.0;
}



Experiment make_experiment(int flag)
{
/* Create Experiment struct and load data from file to insert into members.

   Read user config file here too. At the moment: array file called params

*/
//  char f1[] = "rows";

  Experiment E;

  // these members are fields. Fields are created and filled from file, and malloc is used.

  // E.rows = read_array_int(f1);


  // check if read was succesful, exit if not.
  //check_int_field(E.rows);

  E.S_i = init_state(SPACEDIM, NULL);  
  E.S = init_state(SPACEDIM, NULL);

  return E;
}


