
#ifndef MODELH
#define MODELH

#define INTCONSTS
//#define STEM_NODES
#define SPACEDIM 1   // spatial dimension of dynamical system
#define MAXINT 20

//#define INSTRSET_RARE

// exported vars 

extern int numpar;
extern int maxint;

extern double buffers[STACKSIZE][SPACEDIM];


typedef struct experiment {

  FieldInt rows;

  State S;
  State S_i;

} Experiment;


Population score_pop(Population pop_old,Population pop, Experiment Exp, int compgridsize);

int model_init(Experiment Exp);

int free_experiment(Experiment E);

Node *copy_node(Node *tree);

double get_score(Node *newtree, State S_i, Experiment Exp);

Experiment make_experiment(int flag);

#endif

