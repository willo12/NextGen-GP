
#ifndef MODELH
#define MODELH

#define INTCONSTS
//#define STEM_NODES
#define SPACEDIM 6   // spatial dimension of dynamical system
#define MAXINT 20

//#define INSTRSET_RARE

// exported vars 

extern int numpar;
extern int maxint;

extern double buffers[STACKSIZE][SPACEDIM];

StateComponent storeFun(Node *tree);
StateComponent copyFun(Node *tree);

StateComponent ldxFun(Node *tree);
StateComponent lda_xFun(Node *tree);
StateComponent sta_xFun(Node *tree);
StateComponent stxFun(Node *tree);
StateComponent ldaFun(Node *tree);
StateComponent inxFun(Node *tree);
StateComponent dexFun(Node *tree);


typedef struct experiment {

  State S;
  State S_i;

} Experiment;

int protected_mem(int i);
int model_init(Experiment Exp);
int free_experiment(Experiment E);
IntPars makerandomint_pars(Experiment Exp);
IntPars mut_int_pars(IntPars ip, Experiment Exp);
Node *copy_node(Node *tree);
Population score_pop(Population pop_old,Population pop, Experiment Exp, int compgridsize);
double get_score(Node *newtree, State S_i, Experiment Exp);
Experiment make_experiment(int flag);

#endif

