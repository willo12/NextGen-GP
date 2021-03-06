
#ifndef MODELH
#define MODELH

//#define INTCONSTS
#define STEM_NODES

//#define MAXINT 20

#define INSTRSET_RARE


// exported vars 

extern int numpar;
extern int maxint;

extern double buffers[STACKSIZE][SPACEDIM];


typedef struct experiment {

  int ts_factor;
  int startscore_i;

  Field Iffs;
  Field obs;
  FieldInt I;
  Field result;

  State S;
  State S_i;

} Experiment;


int model_init(Experiment Exp);
int make_itg(Node *tree, char trstr[]);
void make_ic(NodeScore ns, char trstr[]);
int make_itg_fun(NodeScore ns, char template[], char total_str[]);
int free_experiment(Experiment E);
ScalarPars makerandomscalar_pars(Experiment Exp);
ScalarPars mut_scalar_pars(ScalarPars ip, Experiment Exp);
State handle_S_init(double S_init, Experiment Exp);
void c_map_f(char *treestring, double *A, double *B, int shape_i, int shape_j, double S0_start, double S1_start, double S0_end, double S1_end, double *forcing, double *scalars, int forcing_dim);
Population score_pop(Population pop_old,Population pop, Experiment Exp, int compgridsize, int check_exist);
double get_score(NodeScore ns, State S_i, Experiment Exp);

Experiment make_experiment(int ts_factor);
void c_single_tree(char treestr[], int ts_factor, double * S_init_array, double * scalars);

#endif

