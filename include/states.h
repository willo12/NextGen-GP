#ifndef STATESH
#define STATESH

#define WANDERSTEP 0.04


typedef struct state {

  int size;

#ifdef INTSTATES
  int* data;
#else
  double* data;
#endif

} State;

typedef struct int_pars {

  int* data;

} IntPars;


int free_state(State S);
State make_state(int size);

#ifdef INTSTATES
State fill_state(State S, int *A);
State init_state(int size, int *A);
#else
State fill_state(State S, double *A);
State init_state(int size, double *A);
#endif

State make_copy_state(State S);
State copy_state(State S_from, State S_to);
int compare_state(State S1, State S2);  
State makerandomstate(State S_i);
State crossover_states(State S1, State S2);
State mut_state(State S);
int state2bufferline(char treebuffer[], State S );
State str2state(char * text);

#endif


