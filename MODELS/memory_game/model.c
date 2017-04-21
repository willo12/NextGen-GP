
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <compile_options.h>
#include <my_states.h>
#include <states.h>
#include <node_structures.h>
#include <fields.h>
#include <basic_ops.h>
#include <my_ops.h>
#include <tree_io.h>
#include <model.h>
#include <c_nextgen.h>

// local preprocessor directives

#define INTERP_FORCING  // interpolate forcing while solving

#define HILLTHRESH 0.02
#define MAXHILL 8 // max hill climbing steps

#define RNDGRAIN 0.005
#define MAXSEARCH 0

#define BLOWUP 4  // values for S_i beyond which the solution is considered unusable
#define MAXCHANGE 0.1



char* instrset; /* used at runtime */
char* instrset_rare = "T"; /* used at runtime */
char* instrsetswap = "ASMGE"; /* instr that can swap */
//char* instrsetswap = "Z"; // turn off swap
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


char * add_op_code(const char * a, int op)
{
  char tmp_str[3];
  sprintf(tmp_str,"%c",op);
  return stradd(a,tmp_str);
}

int init_reg(Experiment Exp);
Population play_all_others(Population pop, int i_ns, Experiment Exp);
int detect_positions(int location[4], int i_player, int i_insert);
int gridgame(NodeScore players[], Experiment Exp);

// functions

int model_init(Experiment Exp)
{
  char tmp_str[3];

  instrset = "\0";

//  instrset = add_op_code(instrset,0x69);

  instrset = stradd(instrset,"A");
  instrset = stradd(instrset,"S");
  instrset = stradd(instrset,"M");
  instrset = stradd(instrset,"I");
  instrset = stradd(instrset,"G");
  instrset = stradd(instrset,"E");
  instrset = stradd(instrset,"N");
  instrset = stradd(instrset,"R");

  instrset = stradd(instrset,"L");
  instrset = stradd(instrset,"l");
  instrset = stradd(instrset,"x");
  instrset = stradd(instrset,"s");
  instrset = stradd(instrset,"x");
  instrset = stradd(instrset,"i");
  instrset = stradd(instrset,"d");


  printf("Initialized instrset %s\n",instrset);

  op_table['F'] = &funFun; 

  op_table['R'] = &storeFun; 
  op_table['N'] = &copyFun; 

  op_table['L'] = &ldaFun;  
  op_table['l'] = &lda_xFun;  
  op_table['X'] = &ldxFun;  
  op_table['s'] = &sta_xFun;  
  op_table['x'] = &stxFun;  
  op_table['i'] = &inxFun;  
  op_table['d'] = &dexFun;  

  arg_table['L'] = 1;
  arg_table['l'] = 1;
  arg_table['X'] = 2;
  arg_table['s'] = 2;
  arg_table['x'] = 2;
  arg_table['i'] = 1;
  arg_table['d'] = 1; 	

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
    reg[i].value = 0;
  }

  numpar = 20;



  return check_numpar(numpar);

}


IntPars makerandomint_pars(Experiment Exp)
{
  int i;
 
  IntPars ip = make_int_pars(INTPARS);

  // initialize
  for (i=0;i<INTPARS;i++)
  {
    *(ip.data+i) = rand()%MAXINT; 
  }

  return ip;
}


IntPars mut_int_pars(IntPars ip, Experiment Exp)
{

  int i;
  int i_forc, jump_range;
  int rnd;

  IntPars new_ip = make_int_pars(INTPARS);


  for (i=0;i<INTPARS;i++)
  {

    rnd = rand()%10;

    if (rnd < 3)
    {
      new_ip.data[i] = rand()%MAXINT; 
    }
    else
    {
      new_ip.data[i] = ip.data[i];
    }
  }

  return new_ip;
}



Population play_all_others(Population pop, int i_ns, Experiment Exp)
{
  int i, winner;

  NodeScore players[2];

  players[0] = pop.pop[i_ns];
  
  for (i=0;i<pop.popsize;i++) 
  {
//    pop.pop[i_ns].score += gridtraining(pop.pop[i_ns].node, Exp);

    if (i==i_ns)
    {
      continue;
    }

    players[1] = pop.pop[i];
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
  int i, leafcount;
  int *current = make_int(0);
  double error;

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

int protected_mem(int i)
{
  int max_i=numpar -7;

  if (i<SPACEDIM)
  {
    i=SPACEDIM;
  }
  else if (i>max_i)
  {
    i=max_i;
  }
  return i;

}


int detect_positions(int location[4], int i_player, int i_insert)
{

  reg[i_insert].value = location[2*i_player];
  reg[i_insert+1].value = location[2*i_player+1];

  reg[i_insert+2].value = location[2*(1-i_player)];
  reg[i_insert+3].value = location[2*(1-i_player)+1];

  return 0;
}

int update_regs(State S, int location[], int i_player)
{

  int i;
  for (i=0;i<S.size;i++)
  {
    reg[i] = S.data[i];
  }

  detect_positions(location, i_player, numpar-5);

  return 0;
}



int gridgame(NodeScore players[], Experiment Exp)
{
  int i,j,o;
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
      update_regs(players[i].S_i,location,i);      

//      reg[4].value = last_move[i];


      evaluate_tree(players[i].node, S_result);
      move = S_result.data[0].value%4;

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
  int i;

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

StateComponent storeFun(Node *tree)
{
  /* STA. store result of first child and return first child */


  StateComponent to_store = run_node(((Node *) tree->children[0])); // the accumulator
  StateComponent dest = run_node(((Node *) tree->children[1]));

  reg[protected_mem(dest.value%numpar)] = to_store;

  return to_store;

};

StateComponent sta_xFun(Node *tree)
{
  /* STA. store result of first child and return first child */

  StateComponent to_store = run_node(((Node *) tree->children[0])); // the accumulator
  StateComponent dest = run_node(((Node *) tree->children[1]));

  reg[protected_mem((dest.value+reg[numpar-6].value)%numpar)] = to_store;

  return to_store;

};

StateComponent stxFun(Node *tree)
{
  /* STA. store result of first child and return first child */


  StateComponent to_store = reg[numpar-6]; // the accumulator
  StateComponent dest = run_node(((Node *) tree->children[1]));

  reg[protected_mem(dest.value%numpar)] = to_store;

  return run_node(((Node *) tree->children[0]));

};

StateComponent ldxFun(Node *tree)
{
  StateComponent or = run_node(((Node *) tree->children[1])); 

  reg[numpar-6] = reg[or.value%numpar];

  return run_node(((Node *) tree->children[0]));
}

StateComponent ldaFun(Node *tree)
{
  StateComponent or = run_node(((Node *) tree->children[0])); 

  return reg[or.value%numpar];
}

StateComponent lda_xFun(Node *tree)
{
  StateComponent or = run_node(((Node *) tree->children[0])); 

  return reg[(or.value+reg[numpar-6].value)%numpar];
}

StateComponent inxFun(Node *tree)
{

  reg[numpar-6].value += 1;
  reg[numpar-6].value = reg[numpar-6].value%numpar;

  return run_node(((Node *) tree->children[0]));
}

StateComponent dexFun(Node *tree)
{

  reg[numpar-6].value -= 1;

  if (reg[numpar-6].value < 0)
  {
    reg[numpar-6].value = 0;
  }


  return run_node(((Node *) tree->children[0]));
}



StateComponent copyFun(Node *tree)
{
  /* move p_i, p_j. copy and return result of first child */

  StateComponent or = run_node(((Node *) tree->children[0]));
  StateComponent dest = run_node(((Node *) tree->children[1]));

  reg[protected_mem(dest.value)] = reg[or.value%numpar];

  return or;

};


