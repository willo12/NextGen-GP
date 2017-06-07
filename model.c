
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <states.h>
#include <compile_options.h>
#include <node_structures.h>
#include <fields.h>
#include <basic_ops.h>
#include <my_ops.h>
#include <tree_io.h>
#include <model.h>
#include <c_nextgen.h>

// local preprocessor directives

//#define INTERP_FORCING  // interpolate forcing while solving

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
char treebuffer[TREEBUFFER];

// SPECIFIC TO ops_muldim!! CHANGE LATER
//double buffers[STACKSIZE][SPACEDIM];

// local function definitions
void update_regs(State S, int n, int steps_forc, double *ffs);
Experiment RK4(Node *tree, State S_i, Experiment Exp);
int init_reg(Experiment Exp);
double score_fun_basic(Experiment Exp);


// functions

int model_init(Experiment Exp)
{
  init_tables();
  init_reg(Exp);

  return 0;
}


int free_experiment(Experiment E)
{

  free_field(E.Iffs);
  free_field(E.obs);
  free_int_field(E.I);
  free_field(E.result);  
  free_state(E.S);
  free_state(E.S_i);

  return 0;
}


IntPars makerandomint_pars(Experiment Exp)
{
  int i;
 
  IntPars ip = make_int_pars(INTPARS);

  // initialize
  for (i=0;i<INTPARS;i++)
  {
    *(ip.data+i) = 0; 
  }

  // enter specific detailed lines here
  *(ip.data) = SPACEDIM + rand()%(Exp.Iffs.dims.cols-1); 

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

    if (rnd < 5)
    {
      i_forc = ip.data[i] - SPACEDIM;

      if (rnd < 3)
      {
        jump_range = 3;
      }
      else
      {
        jump_range = rand()%30+1;
      }

      i_forc += (rand()%jump_range - ( (int) ( (double) jump_range)/2.0 )  );
      i_forc = i_forc%(Exp.Iffs.dims.cols-1);
      
      if (i_forc < 0)
      {
        i_forc = i_forc + (Exp.Iffs.dims.cols-1);
      }
        
      new_ip.data[i] = i_forc + SPACEDIM;

    }
    else
    {
      new_ip.data[i] = ip.data[i];
    }
  }

  return new_ip;
}


State handle_S_init(double S_init, Experiment Exp)
{ /* set initial value array S_init_array for S array, based on signal S_init 

     S_init_array must be initialized to values before calling this function (e.g. [0,0] )

     Set S_init in config.ini to following values:
     S_init < -10 indicates that initial value of obs must be used to initialize all state variables
     S_init > 100 that all state variables init must have a random component around initial obs value
     S_init > 10 that all state variables init init must be entirely random

*/

  int i;
  int obs1 = Exp.obs.dims.cols;

//  double * S_init_array;

  State S_i = init_state(SPACEDIM, NULL); // create and initialize data with 0.0


  if (obs1 > SPACEDIM)
  {
    obs1 = SPACEDIM;  
  }

  if (S_init < -10.0)
  { 
    /* remember time is removed from obs in python code! hence obs[0]*/
/*    S_init = obs[0];  */

//    fill_state(S_i, Exp.obs.data);
    (*S_i.data) = *Exp.obs.data;

  }
  else if (S_init > 100.0)
  {
    for (i=0;i<obs1;i++)
    {
      (*(S_i.data + i)) =  Exp.obs.data[i] + (S_init*1e-4)*(( (double) (rand() -0.5*RAND_MAX) )/RAND_MAX);

    }
  }
  else if (S_init > 10.0)
  {
    for (i=0;i<SPACEDIM;i++)
    {
      (*(S_i.data + i)) = 4*(( (double) (rand() -0.5*RAND_MAX) )/RAND_MAX);
    }
  }
  else
  {
    for (i=0;i<SPACEDIM;i++)  // THIS COULD BE READ FROM FILE IN FUTURE: INTRODUCE STATE RESTARTS
    {
      (*(S_i.data + i)) = S_init;
    }
  }

  return S_i;
}


//void update_regs(State S, int steps_forc, double *ffs)

void update_regs(State S,int n, int steps_forc, double *ffs)

// uncomment for use with interp_forcing:
//void update_regs(State S,int n, int steps_forc, double t_elapsed, double dt_forc, double *ffs)


{ /* deposit S values and forcing in registers reg[] */

  /* n: number of forcing terms plus time (so n=2 for one forcing term) */


  /* first registers reserved for state vector S */

  int i;

#if SPACEDIM == 1
    reg[0]=S.data[0];
#elif SPACEDIM == 2
    reg[0]=S.data[0];
    reg[1]=S.data[1];
#else

  for (i=0;i<SPACEDIM;i++)  // SPACEDIM for speed
  {
    reg[i]=S.data[i];
  }
#endif

//  forcing_t = ((int) floor(t_elapsed/dt_forc));  /* time index to forcing array */


#if FORCING_TERMS > 0 // optimization for known hard coded number of forcing terms
 #if FORCING_TERMS == 1
    reg[SPACEDIM] = ffs[steps_forc*(FORCING_TERMS+1)+1]; // cols in Iffs forcing file are FORCING_TERMS+1. We pick 2nd col.
 #else
  for (i=0;i<FORCING_TERMS;i++)
  {
    reg[SPACEDIM+i] = ffs[steps_forc*(FORCING_TERMS+1)+1+i];
  }
 #endif

#else // use forcing dimension data from function arg
  /* remaining registers reserved for forcing */
  for (i=0;i<n-1;i++)
  {

#if defined INTERP_FORCING && FORCING_TERMS==0 
    double remainder;
    remainder = t_elapsed - forcing_t*dt_forc;
    reg[SPACEDIM+i] = ffs[steps_forc*n+1+i] + remainder*(ffs[(steps_forc+1)*n+1+i] - ffs[steps_forc*n+1+i])/dt_forc;
#else
    reg[SPACEDIM+i] = ffs[steps_forc*n+1+i];
#endif
  }

#endif

}


Experiment RK4(Node *tree, State S_i, Experiment Exp)
{
/*
  Runge Kutta 4 solver of tree equation. 

  Trees are executed recursively by calling (*op_table[tree->op])(tree,0), where the 0 arguments tells functions to deposit results in vector buffers[0][:]. Calling the tree reads values from reg[:].

  reg[2] is forcing, i.e. the explicit time dependence of the function f.

  result is on the forcing grid, indexed with steps_forc

  m: shape[0] of forcing ffs, time dimension array length of forcing
  n: shape[1] of forcing ffs, the amount of forcing terms plus time (so n=2 for one forcing term)
*/

  int m = Exp.Iffs.dims.rows;
  int n = Exp.Iffs.dims.cols;

//  int t0=Exp.Iffs.data[0];
//  double tend=Exp.Iffs.data[(m-1)*n];
  double dt;
//  int steps=0;
//  double t_elapsed=0;
  int i;
  int i_rel;
  double dt_forc;
  int steps_forc=0;

/*  double S_init=0.3; 
  double S_init=0.345;
*/
//  State S_return = init_state(SPACEDIM, NULL);
/* Runge Kutta k's */

  State k1 = init_state(SPACEDIM, NULL);
  State k2 = init_state(SPACEDIM, NULL);
  State k3 = init_state(SPACEDIM, NULL);
  State k4 = init_state(SPACEDIM, NULL);

  State S = Exp.S; // shorthand
  State S_arg = init_state(SPACEDIM, NULL); // argument to update_regs in intermediate steps
  State S_dot_dt = init_state(SPACEDIM, NULL);

/*  S[0]=S_init; */

  copy_state(S_i, Exp.S); // copy initial conditions (from arg) to state

  dt_forc = Exp.Iffs.data[n]-Exp.Iffs.data[0]; /* time step for forcing */
  dt=dt_forc/Exp.ts_factor; /* time step for integration loop. Typically dt_forc = 2*dt */

/*  fprintf(stderr,"%d %d ", dt_forc, dt);
  fprintf(stderr,"%g ", tend - t0);
*/

  while (steps_forc<m) /* integration loop */
  {

//    fill result
    for (i=0;i<SPACEDIM;i++)
    {
      *(Exp.result.data+(steps_forc*SPACEDIM)+i ) = S.data[i]; // RESULT IS ON FORCING TIME STEP
    }

#if TS_FACTOR > 0
    for (i_rel=0;i_rel<TS_FACTOR;i_rel++)
#else
    for (i_rel=0;i_rel<Exp.ts_factor;i_rel++)
#endif

    {

//    t_elapsed += dt;

// to keep track of integration time, uncomment:
//    steps = step_forc*ts_factor+ts_rel
//    t_elapsed = steps*dt;


    /* Runge Kutta term k1 */

//    update_regs(S, n, steps_forc, t_elapsed, dt_forc, Exp.Iffs.data);
      update_regs(S, n, steps_forc, Exp.Iffs.data);
      evaluate_tree(tree, k1);


      /* Runge Kutta term k2 */
      for (i=0;i<SPACEDIM;i++)
      {
        S_arg.data[i] = S.data[i] + k1.data[i]*dt/2;
      }
      /* deposit S values and forcing in registers reg[] and run tree */
  //    update_regs(S_arg, n, steps_forc, t_elapsed + dt/2, dt_forc, Exp.Iffs.data);  
      update_regs(S_arg, n, steps_forc, Exp.Iffs.data);  
      evaluate_tree(tree, k2);

      /* Runge Kutta term k3 */
      for (i=0;i<SPACEDIM;i++)
      {
        S_arg.data[i] = S.data[i] + k2.data[i]*dt/2;
      }
  //    update_regs(S_arg, n, steps_forc, t_elapsed + dt/2, dt_forc, Exp.Iffs.data);
      update_regs(S_arg, n, steps_forc, Exp.Iffs.data);
      evaluate_tree(tree, k3);

      /* Runge Kutta term k4 */
      for (i=0;i<SPACEDIM;i++)
      {
        S_arg.data[i] = S.data[i] + k3.data[i]*dt;
      }
  //    update_regs(S_arg, n, steps_forc, t_elapsed + dt, dt_forc, Exp.Iffs.data); 
      update_regs(S_arg, n, steps_forc, Exp.Iffs.data); 
      evaluate_tree(tree, k4);

      /* increment S */
      for (i=0;i<SPACEDIM;i++)
      {
        S_dot_dt.data[i] = (k1.data[i] + 2.0*k2.data[i] + 2.0*k3.data[i] + k4.data[i])*dt/6.0;
        S.data[i] += S_dot_dt.data[i];
      }


  /* degeneracy checks */
      for (i=0;i<SPACEDIM;i++)
      {
        if (fabs(S.data[i])>BLOWUP || isnan(S.data[i]) || (fabs( S_dot_dt.data[i] ) > MAXCHANGE  ) )
        {

  /*      fprintf(stderr,"degeneracy at i=%d",i);  */
  //      t_elapsed = tend-t0;
          i = SPACEDIM;

          free_state(k1);
          free_state(k2);
          free_state(k3);
          free_state(k4);

          free_state(S_arg);
          free_state(S_dot_dt);
    
          *Exp.result.data=1e19;
          return Exp;
        }
      }

    }
    steps_forc++;
  };

  free_state(k1);
  free_state(k2);
  free_state(k3);
  free_state(k4);

  free_state(S_arg);
  free_state(S_dot_dt);
 
  return Exp;
}


Experiment euler(Node *tree, State S_i, Experiment Exp)
{
// RESULT IS ON FORCING TIME STEP

  int m = Exp.Iffs.dims.rows;
  int n = Exp.Iffs.dims.cols;

  double dt;

  int i;
  int i_rel;
  double dt_forc;
  int steps_forc=0;

  State k1 = init_state(SPACEDIM, NULL);
  State S = init_state(SPACEDIM, NULL);
  State S_dot_dt = init_state(SPACEDIM, NULL);

  copy_state(S_i, S); // copy initial conditions (from arg) to state

  dt_forc = Exp.Iffs.data[n]-Exp.Iffs.data[0]; /* time step for forcing */
  dt=dt_forc/TS_FACTOR; /* time step for integration loop. Typically dt_forc = 2*dt */

  while (steps_forc<m) 
  {
    for (i=0;i<SPACEDIM;i++)
    {
      *(Exp.result.data+(steps_forc*SPACEDIM)+i ) = S.data[i]; 
    }

    for (i_rel=0;i_rel<TS_FACTOR;i_rel++)
    {
      update_regs(S, n, steps_forc, Exp.Iffs.data);
      evaluate_tree(tree, k1);

      for (i=0;i<SPACEDIM;i++)
      {
        S.data[i] += k1.data[i]*dt;
      }

  /* degeneracy checks */
      for (i=0;i<SPACEDIM;i++)
      {
        if ( (fabs(S.data[i])>BLOWUP) || isnan(S.data[i])  )
        {
          i = SPACEDIM;
          free_state(k1);

          *Exp.result.data=1e19;
          return Exp;
        }
      }
    }
    steps_forc++;
  };

  free_state(k1);
  free_state(S_dot_dt);
 
  return Exp;
}




void c_map_f(char *treestring, double *A, double *B, int shape_i, int shape_j, double S0_start, double S1_start, double S0_end, double S1_end, double *forcing, int forcing_dim)
{
  /* For specific forcing value, create two 2D arrays of tree values for S vector values between -2 and 2. */

  Node *tree;

  int i,j;
  double dS0, dS1;
 
  char tmp_str[MAXTREESTR];

  State S_result = init_state(SPACEDIM, NULL);

  init_tables();

  reg[0] = S0_start;
  reg[1] = S1_start;
  for (i=0;i<forcing_dim;i++)
  {
    reg[2+i] = *(forcing+i);
  }

  dS0 = ((double)  (S0_end-S0_start)/shape_i);
  dS1 = ((double)  (S1_end-S1_start)/shape_j);

/*  fprintf(stderr,"%s\n",treestring);
*/
  strcpy(tmp_str,treestring);  

  if ((check_brackets(treestring,'(') == -1) || (check_brackets(treestring,'[') == -1) )
  {
    fprintf(stderr,"Error in tree string: bracket tree integrity check failed. \n");
    /* abort and trigger failure */
    return;    
  }

  tree = str2node(tmp_str,'(',',');

  if (check_node(tree) == 1)
  {
    fprintf(stderr,"Error: tree integrity check failed.");
        /* abort and trigger failure */
    return;  
  }

/*  printf("dS0: %g, dS1: %g \n",dS0,dS1); */

  for (j=0;j<shape_j;j++)
  {
    reg[0] = S0_start;
    for (i=0;i<shape_i;i++)
    {
      
/*      printf("%g %g %g \n",reg[0],reg[1], reg[2]); */
//      (*op_table[tree->op])(tree,0); /* run tree and deposit in buffers[:] */

      evaluate_tree(tree, S_result);

      A[j*shape_i+i] =  S_result.data[0];
      B[j*shape_i+i] =  S_result.data[1];

      reg[0] += dS0;
    }
    reg[1] += dS1;
  }
  free_node(tree);
}



int init_reg(Experiment Exp)
{
  int i;

  int ffs1 = Exp.Iffs.dims.cols;

  for (i=0;i<MAXPAR;i++)
  {
    reg[i] = 0.0;
  }

  numpar = ffs1 - 1 + SPACEDIM;

  return check_numpar(numpar);
}


double score_fun_basic(Experiment Exp)
{

  
  int j;
  double error, tmp_error0;
 

  error = 0;

#if OBSCOLS == 1
  for (j=Exp.startscore_i;j<Exp.I.dims.rows;j++)
  {      
 /*   fprintf(stderr,"%g # ",obs[j*Exp.obs.cols+1]);  */

# ifdef OBSERROR

    tmp_error0 = Exp.obs.data[j*Exp.obs.dims.cols] - Exp.result.data[Exp.I.data[j]*SPACEDIM]; // error relative to lower bound
    if (tmp_error0 > 0) // result lies below lower bound obs
    {
      tmp_error0 = Exp.obs.data[j*Exp.obs.dims.cols] - Exp.result.data[Exp.I.data[j]*SPACEDIM];
      error += tmp_error0*tmp_error0;
    }
    else 
    {
      tmp_error0 = Exp.obs.data[j*Exp.obs.dims.cols+1] - Exp.result.data[Exp.I.data[j]*SPACEDIM]; // error relative to upper bound
      if (tmp_error0 < 0) // result lies above upper bound obs
      {
        error += tmp_error0*tmp_error0;
      }
    }

# else
    tmp_error0 = Exp.obs.data[j*Exp.obs.dims.cols+k] - Exp.result.data[Exp.I.data[j]*SPACEDIM+k];
    error += tmp_error0*tmp_error0;
# endif

  }

#else
  int k; // k to cycle through obs columns.
  for (k=0;k<Exp.obs.dims.cols;k++)  /* number of obs columns must not exceed SPACEDIM! */
  {
    for (j=Exp.startscore_i;j<Exp.I.dims.rows;j++)
    {      
      tmp_error0 = Exp.obs.data[j*Exp.obs.dims.cols+k] - Exp.result.data[Exp.I.data[j]*SPACEDIM+k];
      error += tmp_error0*tmp_error0;
    }
  }

#endif

  return 1000.0*error/(Exp.I.dims.cols*Exp.I.dims.rows);
}






int make_itg(Node *tree, char trstr[])
{ // make c code integration function
  int i;
  Node * child;
  char tmp_str[MAXTREESTR];
  char tmp_str2[MAXTREESTR];

  strcpy(trstr,"");

  for (i=0;i<SPACEDIM;i++)
  {
    child = ((Node *) tree->children[i]);
    (*code2c_str_table[(int) child->op])(child,tmp_str);
    sprintf(tmp_str2,"S.data[%d]+=%s*dt;\n",i,tmp_str);
    strcat(trstr, tmp_str2); 
  }

//  fprintf(stderr,"made itg: %s\n",trstr);
  return 0;
}

void make_ic(NodeScore ns, char trstr[])
{

//  S.data[0] = -0.4614;
//  S.data[1] = 0.378;
  int i;
  char tmp_str[MAXTREESTR];

  strcpy(trstr,"");  

  for (i=0;i<SPACEDIM;i++)
  {
    sprintf(tmp_str,"S.data[%d]=%g;\n",i, ns.S_i.data[i]);  
    strcat(trstr, tmp_str); 
  }
}

int make_itg_fun(NodeScore ns, char template[], char total_str[])
{

  char filename[30];
  char tmp_str[MAXTREESTR];
  char calc_str[MAXTREESTR];
  char ic_str[MAXTREESTR];
 
  sprintf(filename,"templates/itg.c");
  read_text(filename, template);
 
  make_ic(ns,ic_str); // obtain initial conditions statements S[i] = something;
//  fprintf(stderr,"ic:\n%s\n",tmp_str); 

  make_itg(ns.node,calc_str); // obtain calculations inside integration loop
//  fprintf(stderr,"itg:\n%s\n",tmp_str); 

  sprintf(total_str,template,ic_str,calc_str);

  return 0;
}

Population pop2c(Population pop)
{

  int i;

  char filename[30];
  char file_out[30];

  char tmp_str[MAXTREESTR];

  char template[MAXTEMPLATESIZE];
  char tmp_str_long[MAXTEMPLATESIZE];

  treebuffer[0] = '\0';

  sprintf(filename,"templates/euler1.c");
  read_text(filename, template);
  sprintf(tmp_str_long,template,pop.popsize);

  strcat(treebuffer, tmp_str_long);
  
  sprintf(filename,"templates/itg.c");

  read_text(filename, template);
  for (i=0;i<pop.popsize;i++) /* fill pop_next. i is only incremented once new_ns.node is succesfully added, not discarded.  */
  {
    make_itg_fun(pop.pop[i], template, tmp_str_long);
    strcat(treebuffer, tmp_str_long);
  }
 
  sprintf(filename,"templates/euler2.c");
  read_text(filename, template);
  strcat(treebuffer, template);

  sprintf(file_out,"test_euler.c");
  store_data(file_out,treebuffer,"w");

}

Population score_pop(Population pop_old,Population pop, Experiment Exp, int compgridsize, int check_exist)
{
#ifndef INLINESCORING
  int i, i_tree_old, leafcount;
  double error;

  Point my_loc; // to be used in internal grid

# ifdef COMPILE
  char buff[512];
  pop2c(pop);
  int status = system("clang -o  test_euler  states.c fields.c test_euler.c -I include -lm -O0");

  FILE *in=popen("./test_euler","r");

  i=0;
  while((fgets(buff, sizeof(buff), in)!=NULL) && (i<pop.popsize)  )
  {
    error = atof(buff);

    pop.pop[i].score = error;
    i++;

    pop.toterror += error;
    if (error < pop.minerror)
    {
      pop.minerror = error;
      pop.i_min = i; // keep track of best score
    }
  }
  pclose(in);

# else
  int *current = make_int(0);
  for (i=0;i<pop.popsize;i++) /* fill pop_next. i is only incremented once new_ns.node is succesfully added, not discarded.  */
  {
    if (check_exist == 1)
    {
      my_loc = i2point(i, pop.max_loc.x);
      i_tree_old = does_tree_exist(pop.pop[i],pop_old, my_loc, compgridsize);
    }
    if ( (1 == 0) && (check_exist == 1) && (i_tree_old > -1) )
    {
      error = pop_old.pop[i_tree_old].score;
      pop.trees_reused++;            
      pop.pop[i].score = error;  // parsimony already included previously
    }
    else
    {
#  ifdef EVOLVEIC
      error = get_score(pop.pop[i].node, pop.pop[i].S_i, Exp);
#  else
      error = get_score(pop.pop[i].node, Exp.S_i, Exp);
#  endif
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
# endif

#endif
  return pop;
}


double get_score(Node *newtree, State S_i, Experiment Exp)
{
  double error;

#ifdef DEBUG
  if (check_node(newtree) == 1)
  {
    fprintf(stderr,"Debug Mode: Illegal node encountered in get_score. Aborting. \n");
    exit(-1);
  }
#endif

#ifdef CUSTOM_NODES
  Node *runtree;
  runtree = make_run_node(newtree);
  RK4(runtree,S_i, Exp);
  free_node(runtree);
#else

# ifdef RK4
  RK4(newtree,S_i, Exp);
# else
  euler(newtree,S_i, Exp);
# endif

#endif


  if (Exp.result.data[0] <1e18) 
  {
    error = score_fun_basic( Exp);
 
    return error;
  }
  else
  {
    return 1e19;
  }
  return 1e19;
}



Experiment make_experiment(int ts_factor)
{
/* Create Experiment struct and load data from file to insert into members.

   Read user config file here too. At the moment: array file called params

*/
  int i;

  double S_init;

  char f1[] = "Iffs";
  char f2[] = "obs";
  char f3[] = "I";
  char f4[] = "params";

  Field tmp;

  Experiment E;

  // these members are fields. Fields are created and filled from file, and malloc is used.
  E.Iffs = read_array(f1);
  E.obs = read_array(f2);
  E.I = read_array_int(f3);
  tmp = read_array(f4);

  // check if read was succesful, exit if not.
  check_field(E.Iffs);
  check_field(E.obs);
  check_int_field(E.I);
  check_field(tmp);

  if (E.I.dims.rows  > E.obs.dims.rows )
  {
    fprintf(stderr,"Error: len I %d > len obs %d. Must be smaller. \n", E.I.dims.rows  , E.obs.dims.rows);
    exit(-1);
  }
  
  E.S = init_state(SPACEDIM, NULL);

//  E.S_i = init_state(SPACEDIM, NULL);

  // these members are numbers
  E.startscore_i = *tmp.data; // this is cast to an integer
  S_init = *(tmp.data+1); // blueprint for initial values of S

  if (ts_factor <= -1)
  {
    E.ts_factor = *(tmp.data+2);
  }
  else
  {
    E.ts_factor = ts_factor;
  }

  // construct result from size of forcing Iffs. Results will be saved at same time steps as Iffs
  E.result = make_field(make_dims(SPACEDIM,E.Iffs.dims.rows));

  for (i=0;i<E.I.dims.rows;i++)
  {
    if (E.I.data[i] >= E.result.dims.rows)
    {
      fprintf(stderr,"Error: I contains index %ld greater than result array len %d.\n", E.I.data[i], E.result.dims.rows);
      exit(-1);
    }
  }


//  fprintf(stderr,"Result: %d   I: %d \n",E.result.dims.rows   , E.I.dims.rows );

  E.S_i = handle_S_init(S_init, E); // creates State

  free_field(tmp);

  return E;
}


void c_single_tree(char treestr[], int ts_factor, double * S_init_array)
{
  int i;
  char tmp_str[MAXTREESTR];

  char result_file[] = "result";
  char score_file[] = "score_val";

  Field score_field = make_field(make_dims(1,1));

  Node *newtree;
  
  double error;

  State S_i;

  srand( getpid()+time(NULL)  );
  Experiment Exp = make_experiment(ts_factor);

  init_tables();
  init_reg(Exp);

  if (abs(*S_init_array) > 5.0 )
  {
    S_i = handle_S_init(*S_init_array, Exp);
  }
  else
  {
    S_i = init_state(SPACEDIM, S_init_array);    
  }

/*  sprintf(tmp_str,"(-7.75295e-07,E)V"); */

  strcpy(tmp_str,treestr); /* leave original in tact: string will be altered on interpretation */

  if ((check_brackets(treestr,'(') == -1) || (check_brackets(treestr,'[') == -1) )
  {
    fprintf(stderr, "Bracket check failed, exiting. \n");
    exit(-1);
  }

  newtree = str2node(tmp_str,'(',',');
  
  for (i=0;i<SPACEDIM;i++)
  {
    printf("S%d: %g, ",i,S_i.data[i]);
  }

  printf("ts_factor: %d \n",ts_factor);

  error = get_score(newtree, S_i, Exp);

  printf("score: %g \n",error);

  *score_field.data = error;

  write_array(score_field, score_file);
  write_array(Exp.result, result_file);

  free_experiment(Exp);

/*  printf("tree score: %g with S_init: %g, numpar: %d \n", error,S_init, numpar);
*/
}

