
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

char* instrset = "ASM"; /* used at runtime */
char* instrset_rare = "T"; /* used at runtime */
char* instrsetswap = "ASM"; /* instr that can swap */
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
void update_regs(int n, int i_rows, long int *rows);

// functions

int model_init(Experiment Exp)
{
  init_tables();
  init_reg(Exp);

  return 0;
}


int free_experiment(Experiment E)
{

  free_int_field(E.rows);

  free_state(E.S);
//  free_state(E.S_i);

  return 0;
}





void update_regs(int n, int i_rows, long int *rows)
{ /* deposit S values and forcing in registers reg[] */

  /* n: number of forcing terms plus time (so n=2 for one forcing term) */

  int i;

  /* remaining registers reserved for forcing */
  for (i=0;i<n-1;i++)
  {
    reg[i]= rows[i_rows*n+i];
  }
}






int init_reg(Experiment Exp)
{
  int i;

  int ffs1 = Exp.rows.dims.cols;

  for (i=0;i<MAXPAR;i++)
  {
    reg[i] = 0;
  }

  numpar = ffs1 - 1;

  return check_numpar(numpar);
}

Population score_pop(Population pop_old,Population pop, Experiment Exp, int compgridsize)
{
#ifndef INLINESCORING
  int i, i_tree_old, leafcount;
  int *current = make_int(0);
  double error;

  Point my_loc; // to be used in internal grid

  for (i=0;i<pop.popsize;i++) /* fill pop_next. i is only incremented once new_ns.node is succesfully added, not discarded.  */
  {
    my_loc = i2point(i, pop.max_loc.x);

    i_tree_old = does_tree_exist(pop.pop[i],pop_old, my_loc, compgridsize);
    if (i_tree_old > -1)
    {
      error = pop_old.pop[i_tree_old].score;
      pop.trees_reused++;            
      pop.pop[i].score = error;  // parsimony already included previously
    }
    else
    {
 #ifdef EVOLVEIC
      error = get_score(pop.pop[i].node, pop.pop[i].S_i, Exp);
 #else
      error = get_score(pop.pop[i].node, Exp.S_i, Exp);
 #endif
      leafcount = conparcount(pop.pop[i].node, current);
      *current = 0;

      pop.pop[i].score = error + ((double) PARSIMONY*((double) leafcount)); 
    }

/*  parsimony not included in error reporting to report and elite output files */

    pop.toterror += error;
    if (error < pop.minerror)
    {
      pop.minerror = error;
      pop.i_min = i; // keep track of best score
    }
  }

  free(current);

#endif

  return pop;
}



double get_score(Node *newtree, State S_i, Experiment Exp)
{
  double error=0;

  int i,j;

  State S_result = make_state(SPACEDIM);

  int m = Exp.rows.dims.rows;
  int n = Exp.rows.dims.cols;

  for (i=0;i<m;i++)
  {
    update_regs(n, i, Exp.rows.data);
    evaluate_tree(newtree, S_result);

    for (j=0;j<S_result.size;j++)
    {
      error += (S_result.data[j] - Exp.rows.data[i*n+j+2])*(S_result.data[j] - Exp.rows.data[i*n+j+2]);
    }
//    fprintf(stderr,"(%li, %g) ",(S_result.data[j] - Exp.rows.data[i*n+j+2])*(S_result.data[j] - Exp.rows.data[i*n+j+2]), error  );

  }

  return error;
}



Experiment make_experiment(int flag)
{
/* Create Experiment struct and load data from file to insert into members.

   Read user config file here too. At the moment: array file called params

*/
//  int i;

  char f1[] = "rows";

  Experiment E;

  // these members are fields. Fields are created and filled from file, and malloc is used.

  E.rows = read_array_int(f1);


  // check if read was succesful, exit if not.
  check_int_field(E.rows);
  
  E.S = init_state(SPACEDIM, NULL);

  return E;
}


